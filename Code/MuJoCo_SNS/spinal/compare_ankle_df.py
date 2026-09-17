"""Compare OpenSim and converted-MuJoCo right ankle dorsiflexor capacity.

The supplied OpenSim tables sweep ankle_angle_r from -90 to +90 with the
four dorsiflexors activated. MuJoCo is evaluated at the same coordinate
values, zero velocity, and activation state 1. Muscle tension is reported as
-actuator_force. Torque uses a central-difference tendon-length derivative
with respect to the OpenSim ankle coordinate, so its sign is directly
comparable and does not depend on MuJoCo's raw transmission Jacobian.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mujoco
import numpy as np

from bsolve_ik import apply_eq_followers, build_air_model
from gait_phase import IK_MOT, read_storage


HERE = Path(__file__).parent
ROOT = HERE.parent
OUT = HERE / "ankle_df_results"
FORCE_FILE = ROOT / "Opensim_R_ankle_DF_Force"
TORQUE_FILE = ROOT / "Opensim_R_ankle_DF_Torque"
MUSCLES = ("ext_dig_r", "ext_hal_r", "per_tert_r", "tib_ant_r")
COLORS = ("#0072B2", "#009E73", "#E69F00", "#D55E00")


def activate(model, data, actuator_ids):
    data.act[:] = 0.0
    for actuator_id in actuator_ids:
        address = int(model.actuator_actadr[actuator_id])
        if address < 0:
            raise RuntimeError(
                f"actuator {actuator_id} has no activation state")
        data.act[address] = 1.0


def metrics(reference, prediction):
    reference = np.asarray(reference)
    prediction = np.asarray(prediction)
    scale = max(float(np.max(np.abs(reference))), 1e-12)
    valid = np.abs(reference) > 0.05 * scale
    return {
        "peak_abs_opensim": float(np.max(np.abs(reference))),
        "peak_abs_mujoco": float(np.max(np.abs(prediction))),
        "peak_ratio_mujoco_over_opensim":
            float(np.max(np.abs(prediction)) / scale),
        "median_abs_ratio_on_nonzero_reference":
            float(np.median(np.abs(prediction[valid] / reference[valid]))),
        "normalized_rmse":
            float(np.sqrt(np.mean((prediction - reference) ** 2)) / scale),
        "correlation": float(np.corrcoef(reference, prediction)[0, 1]),
        "sign_agreement_on_nonzero_reference":
            float(np.mean(np.sign(reference[valid]) ==
                          np.sign(prediction[valid]))),
    }


def main():
    OUT.mkdir(exist_ok=True)
    _, force_names, force_values = read_storage(FORCE_FILE)
    _, torque_names, torque_values = read_storage(TORQUE_FILE)
    angle_name = "/jointset/ankle_r/ankle_angle_r/value"
    assert force_names[0] == angle_name and torque_names[0] == angle_name
    angles = force_values[:, 0]
    assert np.allclose(angles, torque_values[:, 0])
    osim_force = np.column_stack([
        force_values[:, force_names.index(name)] for name in MUSCLES])
    osim_torque = np.column_stack([
        torque_values[:, torque_names.index(name)] for name in MUSCLES])

    model, data = build_air_model(id_mode=False)
    actuator_ids = np.asarray([
        mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_ACTUATOR, name)
        for name in MUSCLES], dtype=int)
    if np.any(actuator_ids < 0):
        raise RuntimeError("one or more dorsiflexor actuators are absent")
    joint_id = mujoco.mj_name2id(
        model, mujoco.mjtObj.mjOBJ_JOINT, "ankle_angle_r")
    dof = int(model.jnt_dofadr[joint_id])
    qpos = int(model.jnt_qposadr[joint_id])
    saved = np.load(HERE / "bsolve_out.npz", allow_pickle=True)
    signs = json.loads(str(saved["signs"].item()))
    coordinate_sign = float(signs.get("ankle_angle_r", 1.0))

    mujoco.mj_resetDataKeyframe(model, data, 0)
    base_qpos = data.qpos.copy()
    mj_force = np.zeros_like(osim_force)
    mj_torque = np.zeros_like(osim_torque)
    mj_torque_raw = np.zeros_like(osim_torque)
    delta = 1.0e-5
    for row, angle in enumerate(angles):
        data.qpos[:] = base_qpos
        data.qpos[qpos] = coordinate_sign * np.deg2rad(angle)
        apply_eq_followers(model, data)
        center = data.qpos.copy()

        lengths = []
        for direction in (1.0, -1.0):
            data.qpos[:] = center
            data.qpos[qpos] += direction * coordinate_sign * delta
            apply_eq_followers(model, data)
            data.qvel[:] = 0.0
            data.act[:] = 0.0
            mujoco.mj_forward(model, data)
            lengths.append(data.actuator_length[actuator_ids].copy())
        arm_os = (lengths[0] - lengths[1]) / (2.0 * delta)

        data.qpos[:] = center
        data.qvel[:] = 0.0
        activate(model, data, actuator_ids)
        mujoco.mj_forward(model, data)
        raw_force = data.actuator_force[actuator_ids].copy()
        mj_force[row] = -raw_force
        mj_torque[row] = arm_os * raw_force
        # Diagnostic only. In this converted/patched model MuJoCo's raw
        # actuator_moment is zero for these paths even though tendon length
        # changes with ankle angle. The FD derivative above is therefore the
        # authoritative OpenSim-compatible moment arm, as in bsolve_ik.py.
        mj_torque_raw[row] = (
            coordinate_sign * data.actuator_moment[dof, actuator_ids]
            * raw_force)

    ik_time, ik_names, ik_values = read_storage(IK_MOT)
    del ik_time
    gait_ankle = ik_values[:, ik_names.index("ankle_angle_r")]
    gait_range = [float(gait_ankle.min()), float(gait_ankle.max())]
    joint_range = np.sort(
        coordinate_sign * np.rad2deg(model.jnt_range[joint_id])).tolist()
    in_gait = (angles >= gait_range[0]) & (angles <= gait_range[1])

    result = {
        "method": {
            "activation": 1.0,
            "velocity": 0.0,
            "force": "MuJoCo tension = -data.actuator_force",
            "torque":
                "central-difference d(tendon length)/d(OpenSim ankle angle) "
                "times MuJoCo actuator_force",
            "raw_transmission_crosscheck":
                "diagnostic only; actuator_moment is zero for these "
                "converted paths, so the FD tendon derivative is required",
            "coordinate_sign_mujoco_per_opensim": coordinate_sign,
            "input_header_note":
                "OpenSim tables say inDegrees=no, but the -90..+90 values "
                "are interpreted as degrees by magnitude and sweep intent",
        },
        "angle_deg": [float(angles[0]), float(angles[-1])],
        "opensim_ik_gait_range_deg": gait_range,
        "mujoco_joint_range_deg": joint_range,
        "raw_vs_fd_torque_max_abs_difference_nm":
            float(np.max(np.abs(mj_torque_raw - mj_torque))),
        "raw_actuator_torque_max_abs_nm":
            float(np.max(np.abs(mj_torque_raw))),
        "muscles": {},
    }
    for column, muscle in enumerate(MUSCLES):
        result["muscles"][muscle] = {
            "mujoco_Fmax_N":
                float(model.actuator_gainprm[actuator_ids[column], 2]),
            "force_full_sweep": metrics(
                osim_force[:, column], mj_force[:, column]),
            "torque_full_sweep": metrics(
                osim_torque[:, column], mj_torque[:, column]),
            "force_gait_range": metrics(
                osim_force[in_gait, column], mj_force[in_gait, column]),
            "torque_gait_range": metrics(
                osim_torque[in_gait, column], mj_torque[in_gait, column]),
        }
    result["combined_torque"] = {
        "full_sweep": metrics(osim_torque.sum(axis=1),
                              mj_torque.sum(axis=1)),
        "gait_range": metrics(osim_torque[in_gait].sum(axis=1),
                              mj_torque[in_gait].sum(axis=1)),
    }

    fig, axes = plt.subplots(2, len(MUSCLES), figsize=(15, 6.8),
                             sharex=True)
    for column, (muscle, color) in enumerate(zip(MUSCLES, COLORS)):
        force_ax, torque_ax = axes[:, column]
        for ax in (force_ax, torque_ax):
            ax.axvspan(*gait_range, color="0.90", zorder=-10,
                       label="IK gait range" if column == 0 else None)
            ax.axvline(joint_range[0], color="0.55", ls=":", lw=0.8)
            ax.axvline(joint_range[1], color="0.55", ls=":", lw=0.8)
            ax.grid(alpha=0.20)
        force_ax.plot(angles, osim_force[:, column], color="black", lw=1.5,
                      label="OpenSim")
        force_ax.plot(angles, mj_force[:, column], color=color, lw=1.5,
                      label="MuJoCo")
        torque_ax.plot(angles, osim_torque[:, column], color="black", lw=1.5)
        torque_ax.plot(angles, mj_torque[:, column], color=color, lw=1.5)
        torque_ax.axhline(0, color="0.6", lw=0.7)
        force_ax.set_title(muscle)
        torque_ax.set_xlabel("OpenSim ankle_angle_r (deg)")
    axes[0, 0].set_ylabel("muscle tension at a=1 (N)")
    axes[1, 0].set_ylabel("ankle generalized torque (N m)")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle(
        "Right ankle dorsiflexor capacity: OpenSim reference vs converted "
        "MuJoCo\nshading = measured IK gait range; dotted = MuJoCo joint limits")
    fig.tight_layout()
    fig.savefig(OUT / "ankle_df_opensim_vs_mujoco.png", dpi=220)
    fig.savefig(OUT / "ankle_df_opensim_vs_mujoco.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.axvspan(*gait_range, color="0.90", zorder=-10, label="IK gait range")
    ax.plot(angles, osim_torque.sum(axis=1), color="black", lw=1.8,
            label="OpenSim sum")
    ax.plot(angles, mj_torque.sum(axis=1), color="#0072B2", lw=1.8,
            label="MuJoCo sum")
    ax.axhline(0, color="0.6", lw=0.7)
    ax.set(xlabel="OpenSim ankle_angle_r (deg)",
           ylabel="summed dorsiflexor torque (N m)",
           title="Combined right dorsiflexor torque capacity at activation 1")
    ax.grid(alpha=0.20)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "ankle_df_total_torque.png", dpi=220)
    fig.savefig(OUT / "ankle_df_total_torque.pdf")
    plt.close(fig)

    np.savez_compressed(
        OUT / "ankle_df_comparison.npz", angle_deg=angles,
        muscle_names=np.asarray(MUSCLES), opensim_force=osim_force,
        mujoco_force=mj_force, opensim_torque=osim_torque,
        mujoco_torque=mj_torque, mujoco_torque_raw=mj_torque_raw)
    (OUT / "ankle_df_comparison.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")
    lines = [
        "# OpenSim versus MuJoCo right ankle dorsiflexor capacity",
        "",
        f"IK gait ankle range: {gait_range[0]:.2f} to "
        f"{gait_range[1]:.2f} deg.",
        "",
        "| muscle | force ratio MJ/OS (gait) | torque ratio MJ/OS "
        "(gait) | torque sign agreement |",
        "|---|---:|---:|---:|",
    ]
    for muscle in MUSCLES:
        block = result["muscles"][muscle]
        lines.append(
            f"| {muscle} | "
            f"{block['force_gait_range']['median_abs_ratio_on_nonzero_reference']:.3f} | "
            f"{block['torque_gait_range']['median_abs_ratio_on_nonzero_reference']:.3f} | "
            f"{block['torque_gait_range']['sign_agreement_on_nonzero_reference']:.3f} |")
    total = result["combined_torque"]["gait_range"]
    lines.extend([
        "",
        f"Combined gait-range torque magnitude ratio: "
        f"{total['median_abs_ratio_on_nonzero_reference']:.3f}.",
        "",
        "Ratios are capacity comparisons at activation state 1 and zero "
        "velocity, not dynamic gait-force predictions.",
    ])
    (OUT / "ankle_df_comparison_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"raw-vs-FD torque max difference: "
          f"{result['raw_vs_fd_torque_max_abs_difference_nm']:.3g} N m")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
