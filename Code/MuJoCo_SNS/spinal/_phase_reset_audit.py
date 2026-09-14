"""Dynamic sign audit for the v5 sensory phase-reset pathways.

Two independent checks (run from Code/MuJoCo_SNS/spinal, myo env):

A. SNS pathway signs (dynamic, real compiled dynamics): drive the half-
   centers at DRIVE=4 and feed a tonic 1 nA signal into each new port;
   compare mean RG-E/RG-F potentials against the gain-0 baseline.
   Expect: HIP_EXT_SIG excites RG-E / inhibits RG-F; HIP_FLEX_SIG excites
   RG-F / inhibits RG-E. (No static actuator_moment involved - that trap
   does not apply here, this is the compiled network itself.)

B. Muscle-length sign check (static MuJoCo kinematics - lengths are pure
   transmission kinematics, valid static; the 'static actuator_moment
   lies' trap concerns moments, not lengths): sweep hip_flexion_r and
   verify the runner's ext_sig rises with hip EXTENSION and flex_sig rises
   with hip FLEXION velocity (finite difference).
"""
from __future__ import annotations

import numpy as np

import build_network as bn
import mujoco
from muscle_map import classify
from params import DT, E_HI, G

import runner as R

HERE = __import__("pathlib").Path(__file__).parent

AUDIT_ACTS = [
    # right side: representative pools incl. both signal groups
    "glut_max1_r", "semimem_r", "psoas_r", "iliacus_r", "sar_r", "tfl_r",
    "rect_fem_r", "vas_lat_r", "soleus_r", "tib_ant_r",
    # left side minimal
    "glut_max1_l", "psoas_l", "vas_lat_l", "soleus_l",
]


def run_case(g_e: float, g_f: float, ext_port: float = 0.0,
             flex_port: float = 0.0, dur: float = 12.0):
    G["phase_reset_e"] = g_e
    G["phase_reset_f"] = g_f
    net = bn.build(AUDIT_ACTS, dt=DT, interleg=True)
    u = net.make_inputs()
    u[net.input_index("DRIVE")] = 4.0
    u[net.input_index("POSTURE")] = 1.0
    if ext_port:
        u[net.input_index("HIP_EXT_SIG_r")] = ext_port
    if flex_port:
        u[net.input_index("HIP_FLEX_SIG_r")] = flex_port
    nsteps = int(dur / DT)
    rec = np.zeros((nsteps, 4))   # RG_E_r, RG_F_r, RG_E_l, RG_F_l
    for k in range(nsteps):
        v = net.step(u)
        rec[k] = (v[net.idx["RG_E_r"]], v[net.idx["RG_F_r"]],
                  v[net.idx["RG_E_l"]], v[net.idx["RG_F_l"]])
    tail = rec[nsteps // 2:]
    # period from RG_E_r threshold crossings
    on = tail[:, 0] > 0.5 * E_HI
    rises = np.flatnonzero(np.diff(on.astype(int)) == 1)
    period = float(np.mean(np.diff(rises))) * DT if len(rises) >= 2 else float("nan")
    return dict(e=float(tail[:, 0].mean()), f=float(tail[:, 1].mean()),
                e_amp=float(np.ptp(tail[:, 0])), period=period)


def part_a():
    print("== A. SNS pathway-sign audit (dynamic) ==")
    base = run_case(0.0, 0.0)
    print(f"  baseline       : RG-E {base['e']:.3f} mV  RG-F {base['f']:.3f} mV"
          f"  E-amp {base['e_amp']:.2f}  period {base['period']:.3f} s")
    ok = True
    for g_e, g_f, ext, flex, tag in (
            (1.0, 0.0, 1.0, 0.0, "EXT_SIG=1 (g_e=1)"),
            (0.0, 1.0, 0.0, 1.0, "FLEX_SIG=1 (g_f=1)"),
            (2.0, 2.0, 1.0, 0.5, "both g=2, EXT=1 FLEX=0.5")):
        r = run_case(g_e, g_f, ext, flex)
        de, df = r["e"] - base["e"], r["f"] - base["f"]
        which = "EXT" if ext else "FLEX"
        print(f"  {tag:24s}: RG-E {r['e']:.3f} ({de:+.3f})  RG-F {r['f']:.3f} "
              f"({df:+.3f})  period {r['period']:.3f} s")
        if which == "EXT":
            good = de > 0 and df < 0
        else:
            good = df > 0 and de < 0
        ok &= good
        print(f"    -> {'OK' if good else 'SIGN WRONG'} "
              f"(expect {which}: E up / other down)")
    print("== A result:", "PASS" if ok else "FAIL", "==")
    return ok


def part_b():
    print("== B. muscle-length sign audit (static kinematics) ==")
    model = mujoco.MjModel.from_xml_path(str(R.MODEL))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    Lrange = model.actuator_lengthrange
    Lmid = Lrange.mean(axis=1)
    Lhalf = np.maximum((Lrange[:, 1] - Lrange[:, 0]) / 2, 1e-3)
    acts = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(model.nu)]
    aid = {n: i for i, n in enumerate(acts)}
    ie = np.array([aid[a] for a in acts if a.endswith("_r")
                   and "hip_ext" in classify(a).groups])
    ifl = np.array([aid[a] for a in acts if a.endswith("_r")
                    and "hip_flex" in classify(a).groups])
    print(f"  group sizes r-side: hip_ext {ie.size}, hip_flex {ifl.size}")
    jadr = model.joint("hip_flexion_r").qposadr[0]
    rows = []
    for deg in (-30, -15, 0, 15, 30, 45, 60):
        mujoco.mj_resetDataKeyframe(model, data, 0)
        data.qpos[jadr] = np.radians(deg)
        mujoco.mj_forward(model, data)
        L = data.actuator_length
        ln = (L - Lmid) / Lhalf
        ext_raw = -float(np.mean(ln[ie]))
        # finite-difference flexion velocity (+deg/s): d(flex_sig)/dt
        data2_q = data.qpos[jadr] + np.radians(1.0)
        mujoco.mj_resetDataKeyframe(model, data, 0)
        data.qpos[jadr] = data2_q
        mujoco.mj_forward(model, data)
        ln2 = (data.actuator_length - Lmid) / Lhalf
        vnorm = (ln2 - ln) / DT   # +1 deg/s of hip flexion
        flex_raw = -float(np.mean(vnorm[ifl]))
        rows.append((deg, ext_raw, flex_raw))
        print(f"  hip={deg:+3d} deg: ext_sig {ext_raw:+.3f}  "
              f"flex_vel_sig@+1deg/s {flex_raw:+.3f}")
    ok = (rows[0][1] > rows[-1][1]        # extension increases ext_sig
          and all(r[2] > 0 for r in rows))  # flexion velocity -> positive
    # signal at every angle (the magnitude may vary with length - that is
    # length-dependent spindle sensitivity, not a sign error)
    print("== B result:", "PASS" if ok else "FAIL", "==")
    return ok


if __name__ == "__main__":
    a = part_a()
    b = part_b()
    print("OVERALL:", "PASS" if (a and b) else "FAIL")
