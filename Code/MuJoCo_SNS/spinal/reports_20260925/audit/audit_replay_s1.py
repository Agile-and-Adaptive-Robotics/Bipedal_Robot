"""AUDIT (goal4) - spot-reproduce the STAGE-1 winner of one variant.

Replicates the exact objective code path of _curriculum_{variant}.py
stage 1 (env pins, BASE_MUL merge from best_walk_params_v10.json,
set_stage, runner args, npz re-score), but writes the npz to an
AUDIT-PRIVATE scratch name so no campaign artifact is touched.

Usage: python audit_replay_s1.py w2lvar|syn6
"""
import importlib
import json
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SPINAL = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, SPINAL)
os.chdir(SPINAL)

variant = sys.argv[1]
assert variant in ("w2lvar", "syn6")
scratch_npz = f"audit_replay_{variant}_s1.npz"
if os.path.exists(scratch_npz):
    os.remove(scratch_npz)

# ---- replicate main()'s setup (env pins + BASE_MUL) --------------------
os.environ["AARL_NET"] = variant
os.environ["AARL_NPZ"] = scratch_npz
prev = json.loads(open("best_walk_params_v10.json",
                       encoding="utf-8").read())
BASE_MUL = dict(prev["multipliers"])
BASE_MUL["renshaw"] = 0.5 if variant == "w2lvar" else 0.0
if variant == "syn6":
    BASE_MUL["syn6"] = 1.0
    BASE_MUL["syn6_brainstem"] = 0.0

cur = importlib.import_module(f"_curriculum_{variant}")
J = json.loads(open(f"curriculum_{variant}_stage1.json",
                    encoding="utf-8").read())
p = {**BASE_MUL, **J["params"]}
print(f"[audit:{variant}] json score={J['score']!r} trial={J['trial']}")
print(f"[audit:{variant}] merged p ({len(p)} keys): {p}")

cur.set_stage(1, p)

# ---- stage-1 runner invocation (exact stage-1 args) --------------------
args = ["--no-ground", "--no-afferents", "--no-interleg",
        "--time", "14", "--drive", repr(p["drive"])]
print(f"[audit:{variant}] runner args: {args}", flush=True)
m = cur.R.main(args)
print(f"[audit:{variant}] runner args: {args}", flush=True)
if m is None:
    print(f"[audit:{variant}] runner main() returned None (air route "
          f"returns no metrics dict) - scoring from the npz as the "
          f"objective does")
else:
    print(f"[audit:{variant}] runner summary: {m}")

# ---- exact stage-1 objective recompute (from the AUDIT npz) ------------
import numpy as np
z = np.load(scratch_npz, allow_pickle=True)
t, q, neuro = z["t"], z["q"], z["neuro"]
names = [str(x) for x in z["neuro_names"]]
joints = [str(x) for x in z["key_joints"]]
i_rge, i_knee = names.index("RG_E_r"), joints.index("knee_angle_r")
print(f"[audit:{variant}] cfg={z['cfg']} npz columns: RG_E_r idx "
      f"{i_rge}; knee_angle_r idx {i_knee}")
mask = (t >= 5.0) & (t <= 17.0)
finite = bool(np.all(np.isfinite(q[mask])) and
              np.all(np.isfinite(neuro[mask])))
knee = q[mask, i_knee]
rom_ok = (-360.0 < float(knee.min()) < 360.0) and \
         (-360.0 < float(knee.max()) < 360.0)
rge = neuro[mask, i_rge]
on = rge > 0.5 * max(rge.max(), 1e-9)
rises = int(np.sum(np.diff(on.astype(int)) == 1))
span = float(rge.max()) - float(rge.min())
if not finite:
    score = -200.0
elif not rom_ok:
    score = -200.0
elif rises > 30:
    score = -200.0
elif rises < 3 or span < 1.0:
    score = -10.0 + 0.05 * (-float(knee.min()))
else:
    score = 3.0 * rises + 0.5 * (-float(knee.min()))
jval = float(J["score"])
delta = score - jval
pct = abs(delta) / max(abs(jval), 1e-9) * 100.0
print(f"[audit:{variant}] finite={finite} rom_ok={rom_ok} rises={rises} "
      f"span={span:.4f} knee_min={float(knee.min()):.2f} "
      f"knee_max={float(knee.max()):.2f}")
print(f"[audit:{variant}] RECOMPUTED {score!r} vs json {jval!r} "
      f"delta {delta:+.6f} ({pct:.4f}%)")
print(f"[audit:{variant}] REPLAY:",
      "PASS" if pct <= 1.0 else "FAIL")
# cleanup: keep the evidence npz in the audit folder
import shutil
shutil.move(scratch_npz,
            os.path.join("reports_20260925", "audit", scratch_npz))
sys.exit(0 if pct <= 1.0 else 1)
