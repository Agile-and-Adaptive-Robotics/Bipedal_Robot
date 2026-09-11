"""A/B stability diagnostic for the gait2392 MJCF + pelvis rig (NaN blocker).

Pure MuJoCo (no SNS network): isolates whether the PLANT diverges before the
closed loop is even considered. One suspect at a time:

  --pp-contact      keep pathpoint sphere collisions (old behavior; default
                    in patch_xml is now to disable them)
  --leg-damping X   leg hinge damping (converted default 0.05)
  --boundmass X     pathpoint body boundmass (converted default 0.01)
  --no-ground       suspended in air (ground plane collision off)
  --rigid           rigid y/z rig defaults instead of the compliant x tether
  --act X           constant activation on all muscles (default 0 = passive)
  --seconds T       sim duration (default 3 s)

Reports: keyframe contact census (which geoms touch, forces), per-1s
max |qacc| / max |qvel|, first non-finite time + offending joints, verdict.
"""
from __future__ import annotations

import sys

import numpy as np
import mujoco

import runner as R

BODY = mujoco.mjtObj.mjOBJ_BODY


def geom_name(model, gid):
    return mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_GEOM, gid)


def dof_joint_names(model, dofs):
    # build dof -> joint map once (hinges/slides = 1 dof each in this model)
    out = []
    starts = [(int(model.jnt_dofadr[j]), j) for j in range(model.njnt)]
    starts.sort()
    for b in dofs:
        name = None
        for da, j in starts:
            if da <= b:
                name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, j)
            else:
                break
        out.append(name)
    return out


def main(argv):
    opts = dict(pp_contact=False, leg_damping=None, boundmass=0.1,
                no_ground=False, rigid=False, act=0.0, seconds=3.0)
    args = list(argv)
    while args:
        a = args.pop(0)
        if a == "--pp-contact":
            opts["pp_contact"] = True
        elif a == "--leg-damping":
            opts["leg_damping"] = float(args.pop(0))
        elif a == "--boundmass":
            opts["boundmass"] = float(args.pop(0))
        elif a == "--no-ground":
            opts["no_ground"] = True
        elif a == "--rigid":
            opts["rigid"] = True
        elif a == "--act":
            opts["act"] = float(args.pop(0))
        elif a == "--seconds":
            opts["seconds"] = float(args.pop(0))

    tag = " ".join(argv) if argv else "(defaults: pp-contact off)"
    model = mujoco.MjModel.from_xml_path(str(R.MODEL))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    key_pose = R.capture_pose(model)
    model = R.apply_harness(
        model, data,
        kxy=2.0e5 if opts["rigid"] else 1500.0,
        no_ground=opts["no_ground"],
        leg_damping=opts["leg_damping"],
        boundmass=opts["boundmass"],
        pp_contact=opts["pp_contact"])
    data = mujoco.MjData(model)
    R.seed_pose(model, data, key_pose)
    mujoco.mj_forward(model, data)

    print(f"\n=== diag_stab: {tag}")
    print(f"nq={model.nq} nv={model.nv} nu={model.nu} dt={model.opt.timestep} "
          f"integrator={model.opt.integrator}")

    # ---- contact census at the keyframe ----
    print(f"contacts at keyframe: {data.ncon}")
    for c in range(min(data.ncon, 12)):
        con = data.contact[c]
        f6 = np.zeros(6)
        mujoco.mj_contactForce(model, data, c, f6)
        print(f"  {geom_name(model, con.geom1):24s} <-> "
              f"{geom_name(model, con.geom2):24s} dist={con.dist:+.5f} "
              f"fn={abs(f6[0]):8.2f} N")

    # ---- simulate ----
    nsteps = int(opts["seconds"] / model.opt.timestep)
    t_nan = None
    worst_acc, worst_vel = 0.0, 0.0
    acc_hist = []
    for k in range(nsteps):
        if opts["act"] > 0:
            data.ctrl[:] = opts["act"]
        try:
            mujoco.mj_step(model, data)
        except Exception as e:  # mujoco raises on some instabilities
            print(f"mj_step exception at t={k * model.opt.timestep:.3f}: {e}")
            break
        t = k * model.opt.timestep
        finite = np.all(np.isfinite(data.qacc)) and np.all(np.isfinite(data.qvel))
        if not finite:
            t_nan = t
            break
        ma = float(np.max(np.abs(data.qacc)))
        mv = float(np.max(np.abs(data.qvel)))
        worst_acc = max(worst_acc, ma)
        worst_vel = max(worst_vel, mv)
        acc_hist.append(ma)
        if (k + 1) % int(0.5 / model.opt.timestep) == 0:
            qkey = {jn: np.degrees(
                data.qpos[model.joint(jn).qposadr[0]])
                for jn in ("hip_flexion_r", "knee_angle_r", "ankle_angle_r")}
            print(f"t={t:5.2f}s max|qacc|={ma:10.1f} max|qvel|={mv:7.2f} "
                  f"ncon={data.ncon:3d} "
                  f"knee_r={qkey['knee_angle_r']:7.1f}deg "
                  f"ankle_r={qkey['ankle_angle_r']:6.1f}deg")

    if t_nan is not None:
        bad = np.flatnonzero(~np.isfinite(data.qacc) | ~np.isfinite(data.qvel))
        print(f"NaN/Inf at t={t_nan:.3f}s in dofs: "
              f"{dof_joint_names(model, bad[:10])}")
        print("VERDICT: UNSTABLE")
    else:
        print(f"stable {opts['seconds']}s | worst max|qacc|={worst_acc:.1f} "
              f"worst max|qvel|={worst_vel:.2f} "
              f"final knee_r={np.degrees(data.qpos[model.joint('knee_angle_r').qposadr[0]]):.1f}deg")
        print("VERDICT: OK")


if __name__ == "__main__":
    main(sys.argv[1:])
