"""Replay the archived OpenSim static-optimization gait forces on the human MuJoCo model.

Default outputs are the engineering CSV and force plot beside this script.
--details DIR also runs the historical ablation and three separate gait motions.
"""
import argparse
from contextlib import contextmanager
import csv
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MODELS = ROOT / "Solid_Models/OpenSim/Gait2392_Robotbody"
HUMAN = MODELS / "mjc/gait2392_simbody"
SETUP = MODELS / "zz_bsolve_analyze.xml"
SIGNS = HERE / "hunt_training/runs/ik_audit/fk_signs.json"
os.environ.setdefault("CONDA_PREFIX", str(Path(sys.executable).parent))
sys.path[:0] = [str(HERE / "functional_gait2392"), str(HERE / "physics_fidelity")]

import mujoco
import numpy as np
import opensim
from scipy.interpolate import CubicSpline

import bsolve_ik
import runner
from mechanics import CoupledArms, JOINTS
from static_optimization_force import StaticOptimizationForce

MUSCLES = ("soleus_r", "med_gas_r", "lat_gas_r", "tib_ant_r", "vas_lat_r",
           "rect_fem_r", "semimem_r", "bifemsh_r", "psoas_r", "glut_max2_r",
           "tfl_r", "per_brev_r")


@contextmanager
def working_directory(directory):
    prior = Path.cwd()
    try:
        os.chdir(directory)
        yield
    finally:
        os.chdir(prior)


@dataclass
class Motion:
    time: np.ndarray
    coordinates: tuple
    position: np.ndarray
    speed: np.ndarray


def read_so(quantity):
    path = MODELS / "ResultsBSolve" / f"zz_bsolve_StaticOptimization_{quantity}.sto"
    t, names, values = bsolve_ik.read_mot(path)
    return t, values[:, [names.index(n) for n in MUSCLES]]


def archived_motion(times):
    """Ask AnalyzeTool for the same filtered, differentiated states used by SO."""
    with working_directory(MODELS):
        tool = opensim.AnalyzeTool(str(SETUP))
        model = tool.getModel()
        state = model.initSystem()
        tool.loadStatesFromFile(state)
        states = tool.getStatesStorage()
        labels = states.getColumnLabels()
        names = [labels.get(k) for k in range(labels.getSize())][1:]
        columns_q = [j for j, name in enumerate(names) if name.endswith("/value")]
        coordinates = tuple(names[j].split("/")[-2] for j in columns_q)
        columns_u = [names.index(names[j].rsplit("/", 1)[0] + "/speed") for j in columns_q]
        source_times = np.array([states.getStateVector(k).getTime() for k in range(states.getSize())])
        indices = np.searchsorted(source_times, times - 1e-8)
        if np.any(indices >= len(source_times)) or not np.allclose(source_times[indices], times, atol=1e-8, rtol=0):
            raise ValueError("An archived force time has no matching OpenSim analysis state")
        rows_q, rows_u = [], []
        for index in indices:
            data = states.getStateVector(int(index)).getData()
            rows_q.append([data.get(j) for j in columns_q])
            rows_u.append([data.get(j) for j in columns_u])
    return Motion(times, coordinates, np.array(rows_q), np.array(rows_u))


def independent_motion(path):
    t, names, values = bsolve_ik.read_mot(path)
    if not np.all(np.diff(t) > 0):
        raise ValueError(f"Motion times must strictly increase in {path}")
    for j, name in enumerate(names):
        if name not in bsolve_ik.TRANS_DRIVERS:
            values[:, j] = np.deg2rad(values[:, j])
    return Motion(t, tuple(names), values, CubicSpline(t, values, axis=0)(t, 1))


class SourceModel:
    def __init__(self, path):
        self.model = opensim.Model(str(path))
        self.state = self.model.initSystem()
        self.muscles = [opensim.Thelen2003Muscle.safeDownCast(self.model.getMuscles().get(n))
                        for n in MUSCLES]

    def replay(self, motion, activation, law=None):
        coords = [self.model.getCoordinateSet().get(n) for n in motion.coordinates]
        result = np.empty_like(activation)
        equation_error = 0.
        for k in range(len(motion.time)):
            for j, coordinate in enumerate(coords):
                coordinate.setValue(self.state, motion.position[k, j], False)
                coordinate.setSpeedValue(self.state, motion.speed[k, j])
            self.model.realizeVelocity(self.state)
            result[k] = activation[k] * np.array([
                muscle.calcInextensibleTendonActiveFiberForce(self.state, 1.) for muscle in self.muscles
            ])
            if law is not None:
                independent = law.from_lengths([m.getLength(self.state) for m in self.muscles],
                                               [m.getLengtheningSpeed(self.state) for m in self.muscles],
                                               activation[k])
                equation_error = max(equation_error, np.max(np.abs(independent - result[k])))
        return result, float(equation_error)


class ConvertedModel:
    def __init__(self, generic_source):
        if runner.MODEL.resolve() != (HUMAN / "gait2392_simbody_cvt3.xml").resolve():
            raise ValueError("This replay requires the human Gait2392 conversion")
        raw = mujoco.MjModel.from_xml_path(str(runner.MODEL))
        self.model = runner.apply_harness(raw, mujoco.MjData(raw))
        self.data = mujoco.MjData(self.model)
        self.couplers = CoupledArms(self.model)
        self.signs = {row["coordinate"]: row["sign"] for row in json.loads(SIGNS.read_text())}
        self.law = StaticOptimizationForce.from_models(self.model, generic_source.model, MUSCLES)
        self.ids = self.law.actuator_ids
        self.activation_slots = np.array([self.model.actuator(n).actadr[0] for n in MUSCLES])
        if np.any(self.activation_slots < 0):
            raise ValueError("Expected an activation state per selected muscle")

    def replay(self, motion, activation, signs=None, legacy_position=False):
        signs = self.signs if signs is None else signs
        result = {key: np.empty_like(activation) for key in
                  ("native_active", "source_law", "total", "passive", "stationary_total")}
        # Muscle routes depend on the anatomical joints, independently of a
        # common rigid translation or rotation of the entire pelvis.
        source = [motion.coordinates.index(n) if n in motion.coordinates else None for n in JOINTS]
        for k in range(len(motion.time)):
            mujoco.mj_resetDataKeyframe(self.model, self.data, 0)
            velocity = np.zeros(len(JOINTS))
            for j, index in enumerate(source):
                if index is not None:
                    self.data.qpos[self.couplers.addresses[j]] = signs[JOINTS[j]] * motion.position[k, index]
                    velocity[j] = signs[JOINTS[j]] * motion.speed[k, index]
            if legacy_position:
                historical_cubic_projection(self.model, self.data)
            else:
                bsolve_ik.apply_eq_followers(self.model, self.data)
            self.couplers.update(self.data)
            self.data.qvel[:] = self.couplers.E @ velocity
            self.data.act[:] = 0.
            mujoco.mj_forward(self.model, self.data)
            passive = -self.data.actuator_force[self.ids].copy()
            self.data.act[self.activation_slots] = activation[k]
            mujoco.mj_fwdActuation(self.model, self.data)
            total = -self.data.actuator_force[self.ids]
            result["native_active"][k] = total - passive
            result["passive"][k] = passive
            result["total"][k] = total
            result["source_law"][k] = self.law.evaluate(self.data, activation[k])
            self.data.qvel[:] = 0.
            mujoco.mj_forward(self.model, self.data)
            result["stationary_total"][k] = -self.data.actuator_force[self.ids]
        return result

    def historical(self, qpos, activation):
        predicted = np.empty_like(activation)
        for k, pose in enumerate(qpos):
            self.data.qpos[:] = pose
            self.data.qvel[:] = 0.
            self.data.act[:] = 0.
            self.data.act[self.activation_slots] = activation[k]
            historical_cubic_projection(self.model, self.data)
            mujoco.mj_forward(self.model, self.data)
            predicted[k] = np.abs(self.data.actuator_force[self.ids])
        return predicted


def historical_cubic_projection(model, data):
    """Reproduce the archived comparison only. Its legacy polynomial omits x⁴."""
    for row in np.flatnonzero(model.eq_type == mujoco.mjtEq.mjEQ_JOINT):
        dependent, parent = model.jnt_qposadr[[model.eq_obj1id[row], model.eq_obj2id[row]]]
        x = data.qpos[parent] - model.qpos0[parent]
        c0, c1, c2, c3 = model.eq_data[row, :4]
        data.qpos[dependent] = model.qpos0[dependent] + c0 + c1*x + c2*x*x + c3*x*x*x


def stats(reference, predicted):
    err = predicted - reference
    residual = np.sum(err**2)
    variance = np.sum((reference - reference.mean())**2)
    return {"R2": float(1 - residual / variance) if variance else None,
            "RMSE_N": float(np.sqrt(np.mean(err**2))),
            "bias_N": float(err.mean()),
            "mean_N": float(predicted.mean()),
            "reference_mean_N": float(reference.mean())}


def compare_table(reference, predictions):
    return {name: {tag: stats(reference[:, j], estimates[:, j]) for tag, estimates in predictions.items()}
            for j, name in enumerate(MUSCLES)}


def interpolate(t, rows, request):
    return np.column_stack([np.interp(request, t, column) for column in rows.T])


def write_csv(path, table):
    rows = []
    for name, methods in table.items():
        row = {"muscle": name, "OpenSim_mean_N": methods["native_active"]["reference_mean_N"]}
        for method, values in methods.items():
            row.update({method + "_" + key: value for key, value in values.items() if key != "reference_mean_N"})
        rows.append(row)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def force_plot(path, time, reference, alternatives):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "Arial", "font.size": 10})
    fig, axes = plt.subplots(3, 2, sharex=True, figsize=(7.2, 6.7))
    selected = (("rect_fem_r", "Rectus femoris"), ("bifemsh_r", "Biceps femoris short head"),
                ("soleus_r", "Soleus"), ("med_gas_r", "Medial gastrocnemius"),
                ("vas_lat_r", "Vastus lateralis"), ("tib_ant_r", "Tibialis anterior"))
    for ax, (name, title) in zip(axes.flat, selected):
        j = MUSCLES.index(name)
        ax.plot(time, reference[:, j], color="black", linewidth=1.7, label="OpenSim")
        ax.plot(time, alternatives["native_active"][:, j], "--", color="#4477AA", linewidth=1.5, label="MuJoCo native active")
        ax.plot(time, alternatives["source_law"][:, j], ":", color="#AA3377", linewidth=1.6, label="MuJoCo with OpenSim law")
        ax.set_title(title, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("Active force (N)", fontweight="bold")
    for ax in axes[-1]:
        ax.set_xlabel("Time (s)", fontweight="bold")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=True, title="Legend")
    fig.tight_layout(rect=(0, 0, 1, .91))
    path.parent.mkdir(exist_ok=True, parents=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".png", delete=False) as file:
        pending = Path(file.name)
    try:
        fig.savefig(pending, dpi=180)
        os.replace(pending, path)
    finally:
        pending.unlink(missing_ok=True)
        plt.close(fig)


def diagnostics(folder, motion, activation, reference, converted, primary, subject, generic):
    folder.mkdir(exist_ok=True, parents=True)
    with np.load(HERE / "bsolve_out.npz", allow_pickle=False) as archived:
        old_t, old_q, old_sign = archived["t"], archived["qpos"], json.loads(str(archived["signs"]))
    t = motion.time
    old_a = interpolate(t, activation, old_t)
    old_y = interpolate(t, reference, old_t)
    old = converted.historical(old_q, old_a)
    old_project = converted.replay(motion, activation, old_sign, legacy_position=True)
    old_signs = converted.replay(motion, activation, old_sign)
    ablation = {
        "old_signs_and_cubic_followers_stationary_total": old_project["stationary_total"],
        "old_signs_exact_followers_stationary_total": old_signs["stationary_total"],
        "FK_signs_exact_followers_stationary_total": primary["stationary_total"],
        "FK_signs_exact_followers_moving_total": primary["total"],
        "FK_signs_exact_followers_moving_active": primary["native_active"],
        "same_converted_geometry_source_active_law": primary["source_law"],
    }
    common = (old_t >= t[0]) & (old_t <= t[-1])
    so_common = (t >= old_t[0]) & (t <= old_t[-1])
    alternatives = {k: interpolate(t, primary[k], old_t) for k in ("native_active", "source_law")}
    original_ik = {"historic": old, **alternatives}
    audit = {"SO_time_s": [float(t[0]), float(t[-1])], "SO_samples": len(t),
             "IK_time_s": [float(old_t[0]), float(old_t[-1])], "IK_samples": len(old_t),
             "strict_shared_IK_samples": int(common.sum()),
             "strict_shared_SO_samples": int(so_common.sum()),
             "historic_121_same_endpoint_clamp_only_for_comparison": compare_table(old_y, original_ik),
             "strict_shared_120_IK": compare_table(old_y[common], {k: v[common] for k, v in original_ik.items()}),
             "strict_shared_125_SO": compare_table(reference[so_common], {k: primary[k][so_common] for k in alternatives}),
             "ablation_on_exact_SO_samples": compare_table(reference, ablation),
             "holdouts": {}}
    for filename in ("normal.mot", "crouch1.mot", "crouch4.mot"):
        path = MODELS / "Tutorial1" / filename
        holdout = independent_motion(path)
        unit = np.ones((len(holdout.time), len(MUSCLES)))
        # Fresh default states keep coordinates omitted in a tutorial file at their model defaults.
        source = SourceModel(HUMAN / "gait2392_simbody.osim")
        expected, error = source.replay(holdout, unit, converted.law)
        check = converted.replay(holdout, unit)
        audit["holdouts"][filename] = {
            "case": "Full-activation active capacity on independent motion, not measured force",
            "samples": len(unit), "motion_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "maximum_independent_law_vs_native_OpenSim_N": error,
            "muscles": compare_table(expected, {k: check[k] for k in alternatives})}
        if error > 1e-5:
            raise ValueError(f"Source equation check failed on {filename}: {error} N")
    generic_ref, error = generic.replay(motion, activation, converted.law)
    scaled_ref, _ = subject.replay(motion, activation)
    audit["native_scaled_source_vs_archived_max_N"] = dict(zip(MUSCLES, np.max(abs(scaled_ref-reference), axis=0)))
    audit["generic_opensim_vs_scaled_SO"] = compare_table(reference, {"generic_source": generic_ref})
    audit["independent_law_vs_native_generic_max_N"] = error
    np.savez_compressed(folder / "force_traces.npz", time=t, muscles=MUSCLES, activation=activation,
                        archived_OpenSim_N=reference, native_OpenSim_generic_N=generic_ref, **ablation,
                        old_time=old_t, old_interpolated_OpenSim_N=old_y, old_total_N=old)
    (folder / "force_audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE)
    parser.add_argument("--details", type=Path)
    args = parser.parse_args()
    opensim.Logger.setLevelString("error")
    t, reference = read_so("force")
    ta, activation = read_so("activation")
    np.testing.assert_array_equal(t, ta)
    motion = archived_motion(t)
    subject = SourceModel(MODELS / "subject01_simbody.osim")
    reconstructed, _ = subject.replay(motion, activation)
    if np.max(abs(reference - reconstructed)) > 1e-4:
        raise ValueError("OpenSim states, activation, and force files do not reproduce the archived SO output")
    generic = SourceModel(HUMAN / "gait2392_simbody.osim")
    converted = ConvertedModel(generic)
    traces = converted.replay(motion, activation)
    alternatives = {k: traces[k] for k in ("native_active", "source_law")}
    table = compare_table(reference, alternatives)
    args.output.mkdir(exist_ok=True, parents=True)
    write_csv(args.output / "muscle_force_compare.csv", table)
    force_plot(args.output / "figures/muscle_force_compare.png", t, reference, alternatives)
    if args.details:
        diagnostics(args.details, motion, activation, reference, converted, traces, subject, generic)
    for name in MUSCLES:
        print(f"{name:12}  native {table[name]['native_active']['R2']:.4f}  "
              f"source law {table[name]['source_law']['R2']:.4f}")


if __name__ == "__main__":
    main()
