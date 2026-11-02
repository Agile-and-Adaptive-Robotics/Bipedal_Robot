"""Static anatomical torque-span diagnostic; no training or time integration.

This tests the frozen NMF activation basis, not the complete nonlinear SNS.
Additional reflex/posture inputs may supply authority outside that basis.
"""
import os
for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_key] = "1"
import argparse
import hashlib
import json
import sys
from pathlib import Path

import mujoco
import numpy as np

HERE = Path(__file__).resolve().parent
SPINAL = HERE.parent
REPO = SPINAL.parents[2]
sys.path.insert(0, str(SPINAL / "functional_gait2392"))
from mechanics import CoupledArms, JOINTS

MODEL = REPO / "Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml"
BASIS = SPINAL / "synergy_basis.npz"
REFERENCE = SPINAL / "fsa_results/gait_joint_angles_phase.npz"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def spectrum(matrix):
    values = np.linalg.svd(matrix, compute_uv=False)
    tolerance = values[0] * max(matrix.shape) * np.finfo(float).eps if values.size else 0.
    return {"shape": list(matrix.shape), "rank": int(np.count_nonzero(values > tolerance)),
            "singular_values_Nm_per_unit_input": values.tolist(), "numerical_tolerance": float(tolerance)}


def audit():
    before = {str(p): sha(p) for p in (MODEL, BASIS, REFERENCE)}
    m = mujoco.MjModel.from_xml_path(str(MODEL))
    d = mujoco.MjData(m)
    arms = CoupledArms(m)
    with np.load(BASIS, allow_pickle=False) as basis:
        W = np.zeros((m.nu, 12))
        for side, offset in (("r", 0), ("l", 6)):
            for name, row in zip(basis["muscle_names_" + side], basis["W_" + side]):
                aid = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, str(name))
                if aid < 0:
                    raise ValueError("Basis actuator missing: " + str(name))
                W[aid, offset:offset + 6] = row
    with np.load(REFERENCE, allow_pickle=False) as ref:
        names = [str(x) for x in ref["coordinate_names"]]
        coordinates = ref["r_mean"].copy()
        phase = ref["phase"].copy()
    # Source .mot rotations in this cached extraction are degrees; translations
    # are metres. Internal rotations alone are used, and pelvis remains qpos0.
    columns = [names.index(name) for name in JOINTS]
    samples = []
    for percent in (0., 25., 50., 75.):
        index = int(np.argmin(np.abs(phase - percent)))
        d.qpos[:] = m.qpos0
        d.qvel[:] = 0.
        q = np.radians(coordinates[index, columns])
        for name, angle in zip(JOINTS, q):
            if name.startswith("knee_angle") and angle > 0.:
                raise ValueError("Reference posture hyperextends knee; do not silently clip")
        d.qpos[arms.addresses] = q
        arms.project(d)
        mujoco.mj_forward(m, d)
        R = arms.update(d).copy()
        finite_difference_error = None
        if not samples:
            finite_difference_error = float(np.max(np.abs(R - arms.finite_difference(d))))
        # dF/da for MuJoCo's native muscle at this pose and zero velocity.
        # Activation dynamics and passive forces are deliberately outside this
        # instantaneous active-force derivative, as are BPA pressure maps.
        gain = np.zeros(m.nu)
        for aid in range(m.nu):
            if m.actuator_gaintype[aid] != mujoco.mjtGain.mjGAIN_MUSCLE:
                raise ValueError("Diagnostic requires the native muscle force law")
            gain[aid] = -mujoco.mju_muscleGain(
                d.actuator_length[aid], 0., m.actuator_lengthrange[aid],
                m.actuator_acc0[aid], m.actuator_gainprm[aid, :9])
        full = R * gain[None, :]
        projected = full @ W
        represented = np.flatnonzero(np.any(W != 0., axis=1))
        samples.append({"phase_percent": float(phase[index]),
            "internal_angles_deg": dict(zip(JOINTS, np.degrees(q).tolist())),
            "full_native_muscle_activation_map": spectrum(full),
            "represented_muscles_activation_map": spectrum(full[:, represented]),
            "frozen_bilateral_K6_basis_activation_map": spectrum(projected),
            "active_muscle_count_at_pose": int(np.count_nonzero(gain)),
            "equality_aware_moment_arm_max_finite_difference_error_m": finite_difference_error,
            "row_norm_full_Nm": dict(zip(JOINTS, np.linalg.norm(full, axis=1).tolist())),
            "row_norm_K6_Nm": dict(zip(JOINTS, np.linalg.norm(projected, axis=1).tolist()))})
    after = {str(p): sha(p) for p in (MODEL, BASIS, REFERENCE)}
    if before != after:
        raise RuntimeError("A frozen input changed during the static diagnostic")
    return {"machine": "easteregg2", "operation": "static_design_diagnostic",
        "no_time_integration_or_training": True, "source_sha256": before,
        "anatomical_internal_coordinates": list(JOINTS),
        "floating_base_coordinates": ["pelvis_" + s for s in ("tx", "ty", "tz", "tilt", "list", "rotation")],
        "software_nq": m.nq, "software_nv": m.nv, "native_actuators": m.nu,
        "represented_basis_muscles": int(np.count_nonzero(np.any(W != 0., axis=1))),
        "basis_channels": 12, "reference_scope": "One old walking cycle, four static postures; not running/generalization validation",
        "map": "B=R(q)*diag(native_dF_da(q,velocity=0)); B_basis=B*W",
        "limitations": [
            "W is the fixed NMF activation reconstruction, not the derivative of the compiled SNS conductance circuit.",
            "Native force curves and asymmetrical source fit parameters are retained, not refitted or certified here.",
            "No contact equilibrium, harness removal, pressure limits, muscle activation dynamics, or closed-loop stability is tested.",
            "Linear span/rank does not enforce tension-only cones, activation bounds, or sufficient torque margins.",
            "Reflex and posture pathways can add independent control beyond the fixed primitive basis.",
            "A periodic gait need not request arbitrary independent torques at every internal coordinate."],
        "samples": samples}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({"report": str(args.output), "ranks": [
        [row["phase_percent"], row["full_native_muscle_activation_map"]["rank"],
         row["represented_muscles_activation_map"]["rank"],
         row["frozen_bilateral_K6_basis_activation_map"]["rank"]] for row in result["samples"]]}))
