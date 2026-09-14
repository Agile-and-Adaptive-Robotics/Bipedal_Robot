"""Torque budget: muscle torque SUPPLY vs DEMAND (Ben 2026-09-13: "do
the math to make sure we have sufficient force to produce sufficient
torque to stabilize it").

At the standing keyframe, with the TRUNK FMAX FIX applied:
  1. model mass + weight + total force capacity,
  2. gravity torque demand per relevant dof (qfrc_bias at rest; NOTE:
     for a stance leg the big demands in WALKING come from the GRF
     moment, not gravity of the distal segment - the GRF section below
     covers the ankle),
  3. SUPPLY per functional group = sum(Fmax * |fd moment arm|) at the
     keyframe. fd = central-difference dl/dtheta with the equality
     followers re-projected (= OpenSim's arm; actuator_moment is WRONG
     at the knee through the couplers).
  4. measured realized force of the fixed trunk actuators at act=1
     (their lengthrange is converter-garbage, so the fraction of Fmax
     realizable at this pose is empirical),
  5. walking peak ankle external moment from subject01_walk1 GRF
     (vy x CoP-to-ankle lever; ankle center taken from the standing
     pose, walk pose differs by ~cm).

Usage: python _torque_budget.py
"""
import io
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import mujoco

import bsolve_ik as B
import runner as R
from muscle_map import classify

G = 9.81
SIDE_JOINTS = {"hip": "hip_flexion", "knee": "knee_angle",
               "ankle": "ankle_angle"}
GROUP_JOINT = {"hip_ext": "hip", "hip_flex": "hip", "knee_ext": "knee",
               "knee_flex": "knee", "ankle_pf": "ankle",
               "ankle_df": "ankle", "trunk_ext": "lumbar",
               "trunk_flex": "lumbar"}
GROUPS = ("hip_ext", "hip_flex", "knee_ext", "knee_flex", "ankle_pf",
          "ankle_df", "trunk_ext", "trunk_flex")


def main():
    model = mujoco.MjModel.from_xml_path(str(R.MODEL))
    data = mujoco.MjData(model)
    model = R.apply_harness(model, data)          # patches incl. Fmax fix
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    mujoco.mj_forward(model, data)
    acts = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(model.nu)]
    fmax = np.maximum(model.actuator_gainprm[:, 2], 0.0)
    mass = float(model.body_mass.sum())

    print("== torque budget (trunk Fmax fix APPLIED) ==")
    print(f"mass {mass:.1f} kg -> weight {mass * G:.0f} N | "
          f"sum Fmax {fmax.sum():,.0f} N ({fmax.sum() / mass:.1f} N/kg)")

    dofs = {f"hip_flexion_{s}": model.joint(f"hip_flexion_{s}").dofadr[0]
            for s in ("r", "l")}
    dofs.update({f"knee_angle_{s}": model.joint(f"knee_angle_{s}").dofadr[0]
                 for s in ("r", "l")})
    dofs.update({f"ankle_angle_{s}": model.joint(f"ankle_angle_{s}").dofadr[0]
                 for s in ("r", "l")})
    dofs["lumbar_extension"] = model.joint("lumbar_extension").dofadr[0]
    print("\n-- gravity torque demand, qfrc_bias at the keyframe (N*m) --")
    for jn, adr in dofs.items():
        print(f"  {jn:20s} {data.qfrc_bias[adr]:+8.2f}")

    print("\n-- measured actuator_force at act=1 (fixed trunk set) --")
    for a in ("ercspn_r", "ercspn_l", "intobl_r", "extobl_r", "ext_hal_r"):
        i = acts.index(a)
        data.act[:] = 0.0
        data.act[i] = 1.0
        mujoco.mj_forward(model, data)
        F1 = float(data.actuator_force[i])
        lr = model.actuator_lengthrange[i]
        print(f"  {a:10s} Fmax {fmax[i]:6.0f} N -> act=1 force {F1:+8.1f} N "
              f"({100 * abs(F1) / fmax[i]:5.1f}% of Fmax; L "
              f"{data.actuator_length[i]:.3f} m, range "
              f"[{lr[0]:.2f},{lr[1]:.2f}])")
    data.act[:] = 0.0
    mujoco.mj_forward(model, data)

    # fd moment arms at the keyframe
    rows = list(dofs.values())
    A = B.fd_moments(model, data, data.qpos[None, :].copy(), rows)[0]
    row_of = {jn: i for i, jn in enumerate(dofs)}
    trunk_row = row_of["lumbar_extension"]

    print("\n-- torque SUPPLY: sum(Fmax * |fd arm|) at keyframe (N*m) --")
    print(f"  {'group':10s} {'right':>8s} {'left':>8s} {'both':>8s}")
    supply = {}
    for g in GROUPS:
        vals = []
        for side in ("r", "l"):
            s = 0.0
            for i, a in enumerate(acts):
                if not a.endswith("_" + side):
                    continue
                mi = classify(a)
                if mi is None or g not in mi.groups:
                    continue
                if g.startswith("trunk"):
                    arm = abs(float(A[trunk_row, i]))
                else:
                    jn = SIDE_JOINTS[GROUP_JOINT[g]] + "_" + side
                    arm = abs(float(A[row_of[jn], i]))
                s += fmax[i] * arm
            vals.append(s)
        supply[g] = vals
        print(f"  {g:10s} {vals[0]:8.1f} {vals[1]:8.1f} "
              f"{vals[0] + vals[1]:8.1f}")

    # walking peak joint demand from the bsolve ID residuals (the RIGHT
    # demand numbers: inverse dynamics with the measured GRF along
    # subject01_walk1; a naive vy x CoP lever here mixes walk-frame CoP
    # with the model origin - 8 N*m/kg nonsense, caught immediately)
    print("\n-- walking peak joint demand (bsolve ID residuals, "
          "subject01+GRF) vs supply --")
    z = np.load(R.HERE / "bsolve_out.npz", allow_pickle=True)
    fit = [n for n in z["fit_joints"]]
    tau = z["tau_fit"]
    pair = {"hip_flexion_r": "hip_flex", "knee_angle_r": "knee_ext",
            "ankle_angle_r": "ankle_pf", "lumbar_extension": "trunk_ext"}
    for jn, g in pair.items():
        i = fit.index(jn)
        col = tau[:, i]
        dem = float(np.max(np.abs(col)))
        sup = supply[g][0] + (supply[g][1] if g.startswith("trunk") else 0)
        print(f"  {jn:20s} demand {dem:6.1f} N*m | supply({g}) {sup:6.1f} "
              f"N*m | margin x{sup / max(dem, 1e-9):.2f}")
    print("\n(caveats: arms at the standing pose only - the ankle PF arm "
          "shrinks in late stance plantarflexion, so the 1.4x ankle margin "
          "is the optimistic end; ankle demand 2.7 N*m/kg runs high vs "
          "literature ~1.5, likely the known subtalar/CoP residuals)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
