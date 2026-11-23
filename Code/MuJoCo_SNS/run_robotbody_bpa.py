"""Run the mirrored Gait2392 robotbody with pressure-driven pneumatic muscles.

From the repository root, in the MuJoCo 2.3.7 environment:
    python Code/MuJoCo_SNS/run_robotbody_bpa.py --output <results-directory>

The default runs matched right knee-extensor, knee-flexor, and control trials
with the pelvis fixed and the feet clear. Use --pressures for a new time series,
--spinal for live sensory feedback through the repository's spinal circuit, or
--parameters to substitute actuator dimensions and measured pneumatic dynamics.
Use --engine opensim to run directly on the original mirrored OpenSim paths.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import xml.etree.ElementTree as ET

import numpy as np
from scipy.io import loadmat

os.environ.setdefault("CONDA_PREFIX", str(Path(sys.executable).parent))
import mujoco

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SPINAL = HERE / "spinal"
sys.path.insert(0, str(HERE / "robot_torso/bpa_neural_coordination"))
sys.path.insert(0, str(SPINAL / "functional_gait2392"))
from cad_interface import RoutedBPABank
from mechanics import CoupledArms, JOINTS
from numba import njit

SOURCE = ROOT / "Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_robot/gait2392_robot_cvt3.xml"
ROUTES = HERE / "robot_torso/robotbody_routes.json"
SPECIMENS = ROOT / "Testing_Data/kmax_20.mat"


@njit(cache=True)
def transmission_limit(length, rest, stroke_fraction, other_length, stiffness):
    """Maximum cable tension compatible with a sleeve at its shortest length."""
    return stiffness * np.maximum(0., length - other_length - rest * (1. - stroke_fraction))


def route_assignments():
    routes = {r["name"]: r["actuator_assignment"]
              for r in json.loads(ROUTES.read_text())["routes"] if r["actuator_assignment"]}
    ids = {side: {v["actuator_id_per_side"] for k, v in routes.items() if k.endswith("_" + side)}
           for side in ("l", "r")}
    if len(routes) != 54 or any(i != set(range(1, 28)) for i in ids.values()):
        raise ValueError("The robot requires one route on each side for each of its 27 design IDs.")
    return routes


def model_and_routes(timestep, fixture):
    """Keep the converted spatial paths and their equalities, replacing muscle actuators."""
    source = mujoco.MjModel.from_xml_path(str(SOURCE))
    pose = mujoco.MjData(source)
    mujoco.mj_resetDataKeyframe(source, pose, 0)
    for side in ("l", "r"):
        pose.qpos[source.joint("knee_angle_" + side).qposadr] = np.radians(-30.)
    CoupledArms(source).project(pose)
    mujoco.mj_forward(source, pose)

    root = ET.parse(SOURCE).getroot()
    routes = route_assignments()
    actuators = root.find("actuator")
    for entry in list(actuators):
        name = entry.attrib["name"]
        if name in routes:
            # A unit-gain transmission applies the signed BPA tension in newtons.
            # No fitted biological forces, activation dynamics, or ideal capacities remain.
            ET.SubElement(actuators, "motor", name=name, tendon=entry.attrib["tendon"])
        actuators.remove(entry)
    if {a.get("name") for a in actuators} != set(routes):
        raise ValueError("A registered robot route is absent from the converted model.")

    key = root.find("./keyframe/key")
    if fixture == "fixed-pelvis":
        pelvis = root.find("./worldbody/body[@name='pelvis']")
        joints = pelvis.findall("joint")
        if len(joints) != 6 or [int(source.joint(j.attrib["name"]).qposadr[0]) for j in joints] != list(range(6)):
            raise ValueError("Unexpected pelvis layout in the converted model.")
        for j in joints:
            pelvis.remove(j)
        body = source.body("pelvis").id
        pelvis.set("pos", " ".join(map(str, pose.xpos[body])))
        pelvis.set("quat", " ".join(map(str, pose.xquat[body])))
        key.set("qpos", " ".join(map(str, pose.qpos[6:])))
        root.remove(root.find("contact"))
    else:
        key.set("qpos", " ".join(map(str, pose.qpos)))
        contacts = root.find("contact")
        for pair in list(contacts):
            geoms = (pair.get("geom1", ""), pair.get("geom2", ""))
            if "ground-plane" not in geoms or not any(g.startswith(("calcn", "talus", "toes")) for g in geoms):
                contacts.remove(pair)
        if len(contacts) != 6:
            raise ValueError("Expected bilateral talus, calcaneus, and toe contact pairs.")
    key.attrib.pop("act", None)
    followers = {j.get("joint1") for j in root.findall("./equality/joint")}
    for j in root.findall(".//body/joint"):
        if j.get("name") in followers:
            # The polynomial equality prescribes its full path. The shallow sampled
            # slide bounds in the conversion conflict with that path during flexion.
            j.set("limited", "false")
        elif j.get("name") in JOINTS:
            j.set("limited", "true")
            j.set("solreflimit", "0.001 1")
            j.set("solimplimit", "0.999 0.9999 0.001")
    root.find("compiler").set("meshdir", str(SOURCE.parent))
    root.find("option").set("integrator", "implicitfast")
    root.find("option").set("timestep", str(timestep))
    xml = ET.tostring(root, encoding="unicode")
    model = mujoco.MjModel.from_xml_string(xml)
    names = [model.actuator(i).name for i in range(model.nu)]
    if model.nu != 54 or model.na or np.any(model.actuator_gaintype) or np.any(model.actuator_biastype):
        raise ValueError("Unexpected muscle or native force source in the pneumatic model.")
    if not np.all(model.actuator_gainprm[:, 0] == 1) or not np.all(model.actuator_gear[:, 0] == 1):
        raise ValueError("Robotbody pressure interface requires unit gains and transmission gears.")
    return model, names, routes, xml


def specimen_configuration(names, lengths, routes):
    """Measured sleeve and contraction pairs, modeled installations at the initial pose."""
    specimens = loadmat(SPECIMENS)["kmax20"]
    # The 20 mm pressure–contraction calibration used 300–509 mm sleeves.
    specimens = specimens[(specimens[:, 0] >= 30.) & (specimens[:, 0] <= 50.9)]
    configuration = {"command_mode": "direct_pressure", "actuators": {}}
    for i, name in enumerate(names):
        opposite = name[:-1] + ("l" if name.endswith("r") else "r")
        paired_lengths = lengths[[i, names.index(opposite)]]
        # Leave at least 30 mm of route for fittings and cable. Fit the same
        # sleeve and the same fixed-length offset to both sides of the model.
        choices = specimens[specimens[:, 0] / 100. <= paired_lengths.min() - .030]
        if len(choices):
            length_cm, free_cm, _ = choices[-1]
            rest, free, diameter = length_cm / 100., free_cm / 100., 20
            reference = "Testing_Data/kmax_20.mat (kmax20); Testing_Data/BPA_Test_Results-Moe.xlsx"
        else:
            # StraightForceTest, 11.2 cm worksheet, B1/B2: measured clamp-to-clamp in mm.
            rest, free, diameter = .112, .094, 10
            reference = "Testing_Data/StraightForceTest.xlsx, sheet '11.2 cm, no T', B1:B2"
        if paired_lengths.min() < rest + .030:
            raise ValueError(f"No reference specimen fits the route pair for {name}.")
        max_contraction = (rest - free) / rest
        fixed_length = paired_lengths.mean() - rest * (1. - .4 * max_contraction)
        configuration["actuators"][name] = dict(
            design_id=routes[name]["actuator_id_per_side"], diameter_mm=diameter,
            resting_length_m=rest, free_length_at_620kPa_m=free,
            maximum_contraction_fraction=max_contraction, nonmuscle_length_m=fixed_length,
            parallel_count=1, pressure_max_kpa=620., series_stiffness_N_m=2e5,
            tau_fill_s=.040, tau_vent_s=.060, delay_s=0.,
            provenance=reference + ". Specimen assigned to this simulated route. Cable/fittings set for 40% of free contraction at the starting pose. Count, stiffness, and pressure time constants are modeling choices.",
        )
    return configuration


def read_pressures(path, names, duration):
    """A complete set of chamber targets in kPa, held until the next sample."""
    with open(path, newline="", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file)
        if not reader.fieldnames or len(reader.fieldnames) != len(names) + 1 or set(reader.fieldnames) != {"time_s", *names}:
            raise ValueError("Expected time_s and exactly one column in kPa for each of the 54 routes.")
        rows = list(reader)
    data = np.array([[float(row[key]) for key in ["time_s", *names]] for row in rows])
    if not len(data) or not np.isfinite(data).all() or data[0, 0] != 0 or np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError("Pressure times must be finite, start at zero, and increase strictly.")
    if data[-1, 0] < duration:
        raise ValueError("Pressure input does not cover the requested duration.")
    return data


class SpinalPressure:
    """Direct motor-output–pressure mapping. Afferents use actual BPA length, speed, and force."""

    def __init__(self, names, initial_length, dt, sensory_feedback=True):
        sys.path.insert(0, str(SPINAL))
        import build_network
        import params
        self.dt = .002
        self.sensory_feedback = sensory_feedback
        self.stride = round(self.dt / dt)
        if self.stride < 1 or not np.isclose(self.stride * dt, self.dt, atol=1e-12, rtol=0):
            raise ValueError("The physics step must divide the 2 ms neural interval.")
        self.network = build_network.SpinalNetwork(
            muscles={n: build_network.classify(n) for n in names}, sides=("l", "r"), interleg=True)
        self.network.compile(dt=self.dt)
        self.currents = self.network.make_inputs()
        ports = {n: i for i, n in enumerate(self.network.inputs)}
        self.ports = {k: [ports[k + "_" + n] for n in names] for k in ("Ia", "II", "Ib")}
        self.currents[ports["DRIVE"]] = 2.
        self.currents[ports["POSTURE"]] = .15
        self.motor = [self.network.idx["MN_" + n] for n in names]
        self.initial_length = initial_length.copy()
        self.voltage_range = params.E_HI
        self.output = np.zeros(len(names))
        self.motor_voltage = np.zeros(len(names))

    def step(self, i, length, velocity, tension, force_scale):
        if i % self.stride == 0:
            # The encoder scales follow the CAD model's sensory interface. Tension
            # comes from the preceding BPA update; no MN currents are overwritten.
            if self.sensory_feedback:
                self.currents[self.ports["Ia"]] = np.maximum(velocity, 0.)
                self.currents[self.ports["II"]] = .3 * np.maximum(length - self.initial_length, 0.) / np.maximum(.1 * self.initial_length, .01)
                self.currents[self.ports["Ib"]] = .15 * np.maximum(tension, 0.) / force_scale
            voltage = self.network.step(self.currents)
            self.motor_voltage[:] = voltage[self.motor]
            self.output[:] = np.clip(self.motor_voltage / self.voltage_range, 0., 1.)
        return self.output


class PressureExperiment:
    """One command protocol, pneumatic model, and set of records for both engines."""

    def __init__(self, names, configuration, dt, steps, condition, lengths, pressures, sensory_feedback=True):
        if set(configuration["actuators"]) != set(names) or configuration.get("command_mode") != "direct_pressure":
            raise ValueError("Specify direct pressure for all and only the 54 robot routes.")
        self.bank = RoutedBPABank(names, configuration, dt)
        specs = [configuration["actuators"][n] for n in names]
        self.gains = np.array([r.get("transmission_gain", 1.) for r in specs])
        self.counts = np.array([r["parallel_count"] for r in specs])
        self.input_ref = np.array([r.get("bpa_side_reference_m", 0.) for r in specs])
        self.output_ref = np.array([r.get("output_reference_m", 0.) for r in specs])
        self.limits = self.bank.parameters[:, 8]
        self.dt = dt
        self.pressures = pressures
        self.pulse = names.index(condition) if condition in ("vas_med_r", "bifemsh_r") else None
        if pressures is not None and (np.any(pressures[:, 1:] < 0) or np.any(pressures[:, 1:] > self.limits)):
            raise ValueError("Input pressure lies outside a route's configured pressure range.")
        if condition in ("control", "vas_med_r", "bifemsh_r") and np.min(self.limits) < 600.:
            raise ValueError("The matched example requires a 600 kPa pressure capacity.")
        self.spinal = SpinalPressure(names, lengths, dt, sensory_feedback) if condition == "spinal" else None
        self.baseline = np.full(len(names), 0. if condition == "zero_pressure" else 350.)
        initial = np.zeros(len(names)) if self.spinal else (pressures[0, 1:] if pressures is not None else self.baseline)
        self.bank.pressure[:] = initial
        # Hold the specified initial pressure through any route-specific input delay.
        self.bank.history[:] = initial[:, None]
        self.bank.step_count = self.bank.history.shape[1]
        keys = ["command_kpa", "pressure_kpa", "tension_n", "path_length_m", "path_velocity_m_s",
                "bpa_side_path_m", "single_sleeve_tension_n", "sleeve_length_m"]
        if self.spinal:
            keys += ["motor_voltage_mV", "Ia_input_nA", "II_input_nA", "Ib_input_nA"]
        self.trace = {key: np.empty((steps, len(names))) for key in keys}
        self.trace.update(route_names=np.array(names), initial_pressure_kpa=initial.copy())
        if self.spinal:
            self.trace["spinal_cell_count"] = len(self.spinal.network.idx)

    def step(self, i, lengths, velocity):
        if self.spinal:
            command = self.spinal.step(i, lengths, velocity, self.bank.force/self.gains,
                                       self.bank.parameters[:, 3]/self.gains) * self.limits
            self.trace["motor_voltage_mV"][i] = self.spinal.motor_voltage
            for group in ("Ia", "II", "Ib"):
                self.trace[group + "_input_nA"][i] = self.spinal.currents[self.spinal.ports[group]]
        elif self.pressures is not None:
            row = np.searchsorted(self.pressures[:, 0], i * self.dt, side="right") - 1
            command = self.pressures[row, 1:]
        else:
            command = self.baseline.copy()
            if self.pulse is not None and i * self.dt >= .005:
                command[self.pulse] = 600.
        bpa_path = self.input_ref + (lengths - self.output_ref)/self.gains
        try:
            bundle = self.bank.step(command / self.limits, bpa_path)
        except ValueError as error:
            # A sleeve pulled past its measured rest length cannot sustain an
            # arbitrary cable load. Name the incompatible installations.
            failed = [n for n, f in zip(self.bank.names, self.bank.force) if not np.isfinite(f)]
            if failed:
                raise ValueError(f"At {i*self.dt:g} s the BPA calibration cannot support these route "
                                 f"lengths at the current pressure: {', '.join(failed)}") from error
            raise
        # At the minimum sleeve length a residual in the fitted surface is
        # reacted by the contraction stop. Only the available cable stretch
        # transmits tension to the joints.
        params = self.bank.parameters
        limit = transmission_limit(bpa_path, params[:,0], params[:,1], params[:,2], params[:,4])
        bundle[:] = np.minimum(bundle, limit)
        force = bundle/self.gains
        sleeve_length = bpa_path-self.bank.parameters[:,2]-bundle/self.bank.parameters[:,4]
        for key, value in (("command_kpa", command), ("pressure_kpa", self.bank.pressure),
                           ("tension_n", force), ("path_length_m", lengths), ("path_velocity_m_s", velocity),
                           ("bpa_side_path_m", bpa_path), ("single_sleeve_tension_n", bundle/self.counts),
                           ("sleeve_length_m", sleeve_length)):
            self.trace[key][i] = value
        return force


def step_count(duration, dt):
    steps = round(duration / dt)
    if steps < 1 or not np.isclose(steps * dt, duration, rtol=0, atol=1e-12):
        raise ValueError("Duration must be a positive, whole number of physics steps.")
    return steps


def simulate(model, experiment, steps):
    dt = model.opt.timestep
    d = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, d, 0)
    mujoco.mj_forward(model, d)
    joint_ids = [model.joint(n).id for n in JOINTS]
    addresses = model.jnt_qposadr[joint_ids]
    q = np.empty((steps + 1, len(JOINTS)))
    q[0] = d.qpos[addresses]
    generalized = np.empty((steps, model.nv))
    contacts = np.empty(steps, dtype=int)
    max_error = 0.
    max_parity = 0.
    for i in range(steps):
        mujoco.mj_step1(model, d)
        lengths, velocity = d.actuator_length, d.actuator_velocity
        force = experiment.step(i, lengths, velocity)
        d.ctrl[:] = -force
        expected = -d.actuator_moment.T @ force
        contacts[i] = d.ncon
        mujoco.mj_step2(model, d)
        if not np.isfinite(d.qpos).all() or not np.isclose(d.time, (i + 1) * dt, atol=1e-10, rtol=0):
            raise ValueError(f"Integration failed at step {i + 1}.")
        max_error = max(max_error, float(np.max(abs(d.actuator_force + force))))
        max_parity = max(max_parity, float(np.max(abs(d.qfrc_actuator - expected))))
        generalized[i] = d.qfrc_actuator
        q[i + 1] = d.qpos[addresses]
    if max_error > 1e-9 or max_parity > 1e-9 or np.any(d.qfrc_applied) or np.any(d.xfrc_applied):
        raise RuntimeError("Unexpected force source or mismatch in the native transmission.")
    trace = experiment.trace
    trace.update(time_s=np.arange(steps + 1) * dt, joint_angle_rad=q, joint_names=np.array(JOINTS),
                 generalized_force=generalized, contacts=contacts)
    return trace, {"max_tensile_force_error_N": max_error, "max_generalized_force_error": max_parity,
                   "max_ground_contacts": int(contacts.max())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--engine", choices=("mujoco", "opensim"), default="mujoco")
    parser.add_argument("--duration", type=float, default=.15)
    parser.add_argument("--dt", type=float, default=.00025)
    parser.add_argument("--fixture", choices=("fixed-pelvis", "ground"), default="fixed-pelvis")
    parser.add_argument("--parameters", type=Path, help="Load a complete saved or measured route configuration.")
    parser.add_argument("--design", choices=("baseline","selected","reference-specimens"), default="baseline",
                        help="Appendix C survey design. reference-specimens reproduces the earlier single-sleeve test.")
    parser.add_argument("--mount", choices=("transfer","literal","survey-datum"), default="transfer",
                        help="Use an explicit robotbody transmission datum, a direct fit, or the source survey pose.")
    parser.add_argument("--reference-contraction", type=float,
                        help="Use one starting fraction of the saved stroke for every transfer mount. "
                             "Default: 0.4, with 0.7 on tibialis anterior and peroneus tertius for lengthening reserve.")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--pressures", type=Path, help="CSV of time_s and 54 pressure-target columns in kPa.")
    source.add_argument("--spinal", action="store_true", help="Run the existing nonspiking spinal network with BPA sensory feedback.")
    parser.add_argument("--silence-afferents", action="store_true",
                        help="With --spinal, keep descending and postural drive and silence its Ia, II, and Ib inputs.")
    args = parser.parse_args()
    if not np.isfinite([args.dt, args.duration]).all() or min(args.dt, args.duration) <= 0:
        parser.error("Time step and duration must be finite and positive.")
    if args.reference_contraction is not None and not 0. <= args.reference_contraction <= 1.:
        parser.error("Reference contraction must be within the documented sleeve stroke.")
    if args.engine == "opensim" and args.fixture != "fixed-pelvis":
        parser.error("The native OpenSim model has no foot contact installation. Use its fixed-pelvis fixture.")
    if args.silence_afferents and not args.spinal:
        parser.error("--silence-afferents requires --spinal.")
    args.output.mkdir(exist_ok=True, parents=True)
    # An out-of-range input or integration failure must preserve a previous run.
    # Publish the manifest last, once all requested conditions have completed.
    with tempfile.TemporaryDirectory(dir=args.output, prefix=".robotbody-") as staging:
        output = Path(staging)
        result = run(args, output)
        for file in sorted(output.iterdir(), key=lambda entry: entry.name == "summary.json"):
            file.replace(args.output / file.name)
    print(result)


def run(args, output):
    steps = step_count(args.duration, args.dt)
    supplied = json.loads(args.parameters.read_text()) if args.parameters else None
    survey_sizing = (supplied.get("design") in ("baseline","selected")) if supplied else args.design != "reference-specimens"
    if args.engine == "opensim":
        import robotbody_opensim_bpa as native
        native.osim.Logger.setLevelString("error")
        routes = route_assignments()
        source_model = ROOT / "Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_robot.osim"
        model, names, lengths = native.prepare(source_model, output / "robotbody.osim", routes, JOINTS, survey_sizing)
        version = native.osim.GetVersion()
        coordinates = model.getCoordinateSet()
        joint_limits = np.array([[coordinates.get(n).getRangeMin(), coordinates.get(n).getRangeMax()] for n in JOINTS])
    else:
        model, names, routes, xml = model_and_routes(args.dt, args.fixture)
        d = mujoco.MjData(model)
        mujoco.mj_resetDataKeyframe(model, d, 0)
        mujoco.mj_forward(model, d)
        lengths = d.actuator_length.copy()
        (output / "robotbody.xml").write_text(xml, encoding="utf-8")
        source_model = SOURCE
        version = mujoco.__version__
        joint_limits = np.array([model.jnt_range[model.joint(n).id] for n in JOINTS])
    if supplied:
        configuration = supplied
    elif args.design == "reference-specimens":
        configuration = specimen_configuration(names, lengths, routes)
    else:
        import robotbody_survey
        common_fraction = .4 if args.reference_contraction is None else args.reference_contraction
        reserve = {"tib_ant": .7, "per_tert": .7} if args.reference_contraction is None else {}
        configuration = robotbody_survey.configuration(names, lengths, routes, args.design, args.mount,
                                                      common_fraction, reserve)
    pressures = read_pressures(args.pressures, names, args.duration) if args.pressures else None
    if args.spinal:
        conditions = ["zero_pressure", "spinal"]
    elif pressures is not None:
        conditions = ["supplied"]
    else:
        conditions = ["control", "vas_med_r", "bifemsh_r"]
    (output / "parameters.json").write_text(json.dumps(configuration, indent=2), encoding="utf-8")
    summary = {"engine": args.engine, "engine_version": version, "model_source": str(source_model.relative_to(ROOT)),
               "model_source_sha256": hashlib.sha256(source_model.read_bytes()).hexdigest(),
               "design": configuration.get("design","supplied" if args.parameters else "reference-specimens"),
               "mount": configuration.get("mount","specified route parameters"),
               "physical_sleeves":sum(s["parallel_count"] for s in configuration["actuators"].values()),
               "fixture": args.fixture, "dt_s": args.dt, "duration_s": args.duration, "conditions": {}}
    baseline = None
    for condition in conditions:
        experiment = PressureExperiment(names, configuration, args.dt, steps, condition, lengths, pressures,
                                        sensory_feedback=not args.silence_afferents)
        if args.engine == "opensim":
            trace, result = native.simulate(model, experiment, steps, args.dt, JOINTS)
        else:
            trace, result = simulate(model, experiment, steps)
        tension = trace["tension_n"]
        angles = trace["joint_angle_rad"]
        excursion = np.maximum(0., np.maximum(joint_limits[:, 0] - angles.min(axis=0),
                                              angles.max(axis=0) - joint_limits[:, 1]))
        result.update(routes_producing_tension=int(np.count_nonzero(tension.max(axis=0) > 1e-6)),
                      initial_positive_tensions=int(np.count_nonzero(tension[0] > 1e-6)),
                      peak_route_tension_N=float(tension.max()),
                      final_right_knee_deg=float(np.degrees(angles[-1, JOINTS.index("knee_angle_r")])),
                      maximum_joint_stop_excursion_deg=float(np.degrees(excursion.max())))
        if baseline is None:
            baseline = trace["joint_angle_rad"].copy()
        if len(conditions) > 1:
            change = np.degrees(trace["joint_angle_rad"][-1] - baseline[-1])
            for side, label in (("r", "right"), ("l", "left")):
                result[f"final_{label}_knee_change_from_{conditions[0]}_deg"] = float(change[JOINTS.index("knee_angle_" + side)])
        if condition == "spinal":
            result["spinal_cells"] = trace["spinal_cell_count"]
            result["afferent_inputs"] = "silenced" if args.silence_afferents else "active"
        np.savez_compressed(output / (condition + ".npz"), **trace)
        summary["conditions"][condition] = result
    result = json.dumps(summary, indent=2)
    (output / "summary.json").write_text(result, encoding="utf-8")
    return result


if __name__ == "__main__":
    main()
