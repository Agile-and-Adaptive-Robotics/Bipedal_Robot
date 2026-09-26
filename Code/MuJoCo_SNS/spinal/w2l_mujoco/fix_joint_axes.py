"""M3 body fix: correct the 8 hinge axes in the M1 MJCF (v2, MuJoCo-native).

DEFECT (found 2026-09-25, milestone 3): make_w2l_mjcf.py emit_body composed
the joint frame with the body's PARENT-RELATIVE rotation (b.R_al) and wrote
the resulting WORLD-frame axis into MuJoCo's body-LOCAL `axis` attribute:

    jR = mat_mul(b.R_al, j["R_local_al"])     # wrong base: needs R_world
    axis_mj = al2mj_vec(mat_vec(jR, [1,0,0])) # world vector -> body-local slot

For femur/toe the error is small; for tibia/foot (~90 deg rotations) the knee
and ankle shipped as VERTICAL-axis (yaw) joints with ~zero sagittal muscle
moment arm (kinematic probe: knee tendon arms 0.00-0.36 mm per +10 deg).
Consequences: passive "crumple" is actually a leg twirl; knee/ankle muscles
produce no sagittal torque; M2's closed loop could never step.

FIX (this file; aproj and M1 artifacts untouched): the desired WORLD axis per
joint comes from the aproj ground truth computed with the generator's own
(correct) dump formula  world_AL = R_world(child) @ R_joint(local) @ e_x
(the same numbers as M1 report section 4b), mapped AL->MJ. The body-LOCAL
axis MuJoCo needs is then read off MuJoCo's own body rotation:
    axis_local = R_worldbody_MJ^T @ world_MJ.
Write w2l_mjcf_fixed.xml; verify every hinge is lateral in world.

Run:  C:\\Users\\Ben Bolen\\.conda\\envs\\myo\\python.exe fix_joint_axes.py
"""
from __future__ import annotations

import io
import os
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import make_w2l_mjcf as G  # noqa: E402  (parses the aproj read-only; main-guarded)

SRC = os.path.join(HERE, "w2l_mjcf.xml")
DST = os.path.join(HERE, "w2l_mjcf_fixed.xml")


def main() -> int:
    m = G.parse_aproj()
    G.world_transforms(m["root"])
    # name -> aproj body WITH its joint (segment names are shadowed by
    # LinearHillMuscle overlay bodies of the same name; keep the jointed one)
    jbody = {}
    def rec(b):
        if b.joint is not None:
            jbody[joint_key(b)] = b
        for c in b.children:
            rec(c)
    def joint_key(b):
        return b.joint["name"]
    rec(m["root"])
    assert len(jbody) == 8, f"expected 8 jointed bodies, got {sorted(jbody)}"

    txt = open(SRC, encoding="utf-8").read()

    # ---- knee RANGE sign reconciliation (M3, documented) ----
    # Kinematic evidence (probe5b on the axis-fixed body): the knee FLXOR
    # muscles (femur-back -> tibia-back) lengthen per +1 deg of qpos, i.e.
    # knee FLEXION = -qpos, while the shipped limit is [0, +60] (positive
    # only) -> flexion blocked at 0. The aproj LowerLimit=0 is "straight"
    # and UpperLimit=60 the flexion extreme in APROJ coordinates, so the
    # range must be negated in MuJoCo coordinates: [-60, 0]. Hip/ankle ranges
    # are left as transported (hip straddles rest 0; ankle keeps its DF arc).
    for jn in ("knee_L", "knee_R"):
        pat = re.compile(r'(<joint name="%s"[^>]*?range=")([^"]+)(")' % jn)
        hits = pat.findall(txt)
        assert len(hits) == 1, f"{jn}: expected 1 range, got {len(hits)}"
        lo, hi = (float(x) for x in hits[0][1].split())
        assert lo == 0.0 and hi > 0, f"{jn}: unexpected range {lo},{hi}"
        txt = pat.sub(lambda mo, lo=lo, hi=hi:
                      mo.group(1) + "%.10g %.10g" % (-hi, -lo) + mo.group(3),
                      txt, count=1)
        print(f"  knee range {jn}: [{lo:g},{hi:g}] -> [{-hi:g},{-lo:g}] "
              "(flexion = -qpos; see comment above)")

    os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
    import numpy as np
    import mujoco

    # body world rotations from the CURRENT xml (axes do not affect frames)
    m0 = mujoco.MjModel.from_xml_path(SRC)
    d0 = mujoco.MjData(m0)
    mujoco.mj_forward(m0, d0)

    print("joint     world_AL (aproj)        world_MJ desired      body-local axis written")
    new_axes = {}
    for jname in sorted(jbody):
        b = jbody[jname]
        jR_world = G.mat_mul(b.R_world, b.joint["R_local_al"])
        ax_al = G.mat_vec(jR_world, [1.0, 0.0, 0.0])       # world AL axis (M1 table 4b)
        ax_mj = G.al2mj_vec(ax_al)                          # desired world MJ axis
        jid = mujoco.mj_name2id(m0, mujoco.mjtObj.mjOBJ_JOINT, jname)
        bid = m0.jnt_bodyid[jid]
        R = d0.xmat[bid].reshape(3, 3)
        loc = R.T @ np.array(ax_mj)
        new_axes[jname] = "%.9g %.9g %.9g" % tuple(loc)
        print(f"{jname:<9} ({ax_al[0]:+.4f},{ax_al[1]:+.4f},{ax_al[2]:+.4f})"
              f"  ({ax_mj[0]:+.4f},{ax_mj[1]:+.4f},{ax_mj[2]:+.4f})"
              f"  ->  ({loc[0]:+.6f}, {loc[1]:+.6f}, {loc[2]:+.6f})")
        pat = re.compile(r'(<joint name="%s"[^>]*?axis=")([^"]+)(")' % re.escape(jname))
        hits = pat.findall(txt)
        assert len(hits) == 1, f"{jname}: expected 1 joint tag, got {len(hits)}"
        txt = pat.sub(lambda mo: mo.group(1) + new_axes[jname] + mo.group(3), txt, count=1)

    hdr = ("<!-- M3 AXIS FIX (fix_joint_axes.py): identical to w2l_mjcf.xml except\n"
           "     the 8 hinge axes are corrected to the aproj ground truth\n"
           "     (world AL axis = R_world @ R_joint @ e_x per make_w2l_mjcf.py's own\n"
           "     dump formula; emitted in MuJoCo body-local coords). The generator's\n"
           "     emit_body composed the joint frame with b.R_al (parent-relative)\n"
           "     and wrote a world-frame vector into the body-local axis slot, so\n"
           "     knee/ankle shipped as vertical-axis yaw joints with ~zero sagittal\n"
           "     muscle moment arm. M1 report section 4b world-axis table is the\n"
           "     reference this patch reproduces. -->\n")
    open(DST, "w", encoding="utf-8").write(hdr + txt)
    print("wrote", DST)

    # ---------------- verification pass ----------------
    mm = mujoco.MjModel.from_xml_path(DST)
    dd = mujoco.MjData(mm)
    mujoco.mj_forward(mm, dd)
    print("\nverified world axes (MJ): all 8 must be ~lateral (|y| > 0.99)")
    ok = True
    for i in range(mm.njnt):
        n = mujoco.mj_id2name(mm, mujoco.mjtObj.mjOBJ_JOINT, i)
        ax = dd.xmat[mm.jnt_bodyid[i]].reshape(3, 3) @ mm.jnt_axis[i]
        lat = abs(ax[1]) > 0.99
        ok &= lat
        print(f"  {n:<8} ({ax[0]:+.4f}, {ax[1]:+.4f}, {ax[2]:+.4f})"
              f"  {'lateral OK' if lat else 'NOT LATERAL -- FAIL'}")
    print("AXIS FIX:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
