"""M3 closed-loop air smoke: single-LH-RG original CPG -> M1-fixed body.
Writes w2l_air.xml from the FIXED xml (root lifted +0.30 m), runs the loop,
prints joint ranges / ncon / simple alternation. Knobs via args te tf tau
ctrl cap damp dur.
"""
import io, os, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import numpy as np
import mujoco

W2L = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
sys.path.insert(0, W2L)

kv = {}
for arg in sys.argv[1:]:
    if arg.startswith("--") and "=" in arg:
        k, v = arg[2:].split("=" , 1)
        kv[k] = float(v)
TE = kv.get("te", 3.0)
TF = kv.get("tf", 4.0)
TAU = kv.get("tau", 0.25)
CAP = kv.get("cap", 1.0)
DAMP = kv.get("damp", 1.5)
DUR = kv.get("dur", 6.0)

FIXED = os.path.join(W2L, "w2l_mjcf_fixed.xml")
AIR = os.path.join(W2L, "w2l_air.xml")
txt = open(FIXED, encoding="utf-8").read()
old = 'pos="-3.454 0 0.99298"'
assert old in txt
hdr = ("<!-- AIR VARIANT for milestone 3 (test_w2l_air.py): identical to\n"
       "     w2l_mjcf_fixed.xml (axis-fixed body, see fix_joint_axes.py) except\n"
       "     the welded Root is lifted +0.30 m so NO leg pose can reach the\n"
       "     ground plane (leg vertical reach ~0.95 m). Makes 'no ground\n"
       "     contact' provable (ncon == 0 all run). -->\n")
open(AIR, "w", encoding="utf-8").write(hdr + txt.replace(old, 'pos="-3.454 0 1.29298"'))

import build_w2l_orig_net as B
B.NAP["tau_max_h"] = TAU
net = B.build()
m = mujoco.MjModel.from_xml_path(AIR)
d = mujoco.MjData(m)
m.dof_damping[:] = DAMP
if kv.get("stiff", 1.0) > 0.0:
    # stiffer joint-limit constraints (runtime-only; defaults are too soft for
    # 1000-1500 N muscles on 0.4-3.5 kg segments - measured 20 deg violations)
    m.jnt_solimp[:] = np.array([0.9, 0.99, 0.001, 0.5, 2.0])
    m.jnt_solref[:] = np.array([0.006, 1.0])
acts = {a: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, a)
        for a in net.muscle_outputs}
JNTS = ("hip_L", "knee_L", "ankle_L", "hip_R", "knee_R", "ankle_R")
qadr = {j: m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)] for j in JNTS}

u = net.make_inputs()
iE = net.input_index("TONIC L RG ext")
iF = net.input_index("TONIC L RG flx")
iS = net.input_index("Stimulus_1")

NPHY = int(DUR / 0.001)
NSUB = 2                      # physics steps per 2 ms net step
qlog = np.zeros((NPHY, 6))
ncon_max = 0
for i in range(NPHY):
    t = i * 0.001
    if i % NSUB == 0:
        u[iE] = TE; u[iF] = TF
        u[iS] = 10.0 if t < 0.01 else 0.0
        V = net.step(u)
        for a, c in net.muscle_ctrl(V).items():
            d.ctrl[acts[a]] = min(c, CAP)
    mujoco.mj_step(m, d)
    ncon_max = max(ncon_max, d.ncon)
    qlog[i] = [np.degrees(d.qpos[qadr[j]]) for j in JNTS]
    if not np.isfinite(d.qpos).all():
        print(f"[FAIL] non-finite qpos at t={t:.3f}")
        sys.exit(1)

an = qlog[int(1.0 * 1000):]     # drop the first second (kickoff transient)
flex = -an                       # flexion-positive convention
print(f"te={TE} tf={TF} tau={TAU} cap={CAP} damp={DAMP} dur={DUR}"
      f" | ncon_max={ncon_max} (0 = airborne)")
names = ("hipL", "kneeL", "anklL", "hipR", "kneeR", "anklR")
for k, n in enumerate(names):
    q = flex[:, k]
    print(f"  {n}: {q.min():+7.1f}..{q.max():+7.1f} deg (range {q.max()-q.min():5.1f})"
          f"  mean {q.mean():+6.1f}")
hl = flex[:, 0] - flex[:, 0].mean()
hr = flex[:, 3] - flex[:, 3].mean()
r = float(np.corrcoef(hl, hr)[0, 1])
print(f"  hip L/R flexion Pearson r = {r:+.3f} (antiphase want < -0.5)")
# period from L hip flexion peaks
from numpy import diff, flatnonzero
sig = flex[:, 0]
on = sig > sig.mean() + 0.3 * sig.std()
pk = flatnonzero(on[1:] & ~on[:-1]) + 1
if len(pk) >= 3:
    print(f"  L hip flexion excursions: {len(pk)}, mean interval "
          f"{diff(pk).mean():.0f} ms")
else:
    print(f"  L hip flexion excursions: {len(pk)} (too few)")
