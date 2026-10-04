"""GOAL-1 baseline eval driver (spiking campaign 2026-10-02, easteregg2).

Runs the current production set through the documented winner-loading
flows and prints every headline metric:

  s3k    : AARL_NET unset (default 410-n build), curriculum_stage3.json
           (curr_s3k_nocross t34) via _curriculum.set_stage(4, ..., full_rules=0)
           == the documented rescore flow (reports_20260925/tmp/
           rescore_s3k_fixed_ref.py). Expected anchor -161.56754173676563.
  w2lvar : AARL_NET=w2lvar, curriculum_w2lvar_stage5.json (t68 winner)
           via _curriculum_w2lvar.set_stage(5, ...) == the supgate replay flow.
  syn6   : AARL_NET=syn6, curriculum_syn6_stage5.json (t0 seed, study argmax)
           via _curriculum_syn6.set_stage(5, ...).

NOTE (stdout wrap trap, learned from _supgate_replay_w2lvar_s5_t4.py):
no sys.stdout wrapper here - the curriculum modules wrap stdout at
import; an extra wrapper orphans the previous one.

Usage: D:\\Anaconda\\envs\\myo\\python.exe baseline_eval.py s3k|w2lvar|syn6
"""
import importlib
import json
import os
import sys
from pathlib import Path

SPINAL = Path(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
os.chdir(SPINAL)
sys.path.insert(0, str(SPINAL))

which = sys.argv[1] if len(sys.argv) > 1 else "s3k"

if which == "spiking":
    # stretch (goal-1 D): ONE ground eval attempt of the spiking mirror
    # at the s3k winner params (the default-architecture mirror; the
    # spiking net is NOT tuned - mapped parameters only)
    os.environ["AARL_NET"] = "spiking"
    os.environ["AARL_NPZ"] = "spinal_run_spkbase_spiking.npz"
    import _curriculum as C
    import params as P
    import runner as R

    mul = json.loads((SPINAL / "best_walk_params_v10.json")
                     .read_text(encoding="utf-8"))["multipliers"]
    st3 = json.loads((SPINAL / "curriculum_stage3.json")
                     .read_text(encoding="utf-8"))["params"]
    C.BASE_MUL = dict(mul)
    C.BASE_MUL["renshaw"] = 0.5
    C.set_stage(4, {**mul, **st3, "full_rules": 0.0})
    P.G["full_rules"] = 0.0
    drive = st3["drive"]
elif which == "s3k":
    assert "AARL_NET" not in os.environ, "s3k baseline: AARL_NET must be unset"
    os.environ["AARL_NPZ"] = "spinal_run_spkbase_s3k.npz"
    import _curriculum as C
    import runner as R

    mul = json.loads((SPINAL / "best_walk_params_v10.json")
                     .read_text(encoding="utf-8"))["multipliers"]
    st3 = json.loads((SPINAL / "curriculum_stage3.json")
                     .read_text(encoding="utf-8"))["params"]
    C.BASE_MUL = dict(mul)
    C.BASE_MUL["renshaw"] = 0.5
    C.set_stage(4, {**mul, **st3, "full_rules": 0.0})
    import params as P
    P.G["full_rules"] = 0.0
    drive = st3["drive"]
elif which in ("w2lvar", "syn6"):
    os.environ["AARL_NET"] = which
    os.environ["AARL_NPZ"] = f"spinal_run_spkbase_{which}.npz"
    mod = importlib.import_module(
        "_curriculum_w2lvar" if which == "w2lvar" else "_curriculum_syn6")
    import runner as R

    win = json.loads((SPINAL / f"curriculum_{which}_stage5.json")
                     .read_text(encoding="utf-8"))
    prev = json.loads((SPINAL / "best_walk_params_v10.json")
                      .read_text(encoding="utf-8"))
    BASE_MUL = dict(prev["multipliers"])
    BASE_MUL["renshaw"] = 0.5
    p = {**BASE_MUL, **win["params"]}
    mod.set_stage(5, p)
    drive = p["drive"]
    print(f"{which} stage-5 winner: study={win.get('study')} "
          f"trial={win.get('trial')} recorded_score={win.get('score')!r}",
          flush=True)
else:
    raise SystemExit(f"unknown variant {which!r}")

print(f"BASELINE EVAL {which}: AARL_NPZ={os.environ['AARL_NPZ']} "
      f"drive={drive!r}", flush=True)
m = R.main(["--eval", "--drive", repr(drive)])
k = m.get("kine") or {}
print(f"== baseline {which} ==", flush=True)
print(f"kine_score = {m['kine_score']!r}", flush=True)
print(f"nan = {m['nan']}  t_end = {m['t_end']:.2f}  dx = {m['dx']:.3f} m",
      flush=True)
print(f"kz = {m['kz']:.4f}  tilt_max = {m['tilt_max']:.2f} deg  "
      f"knee_min = {m['knee_min']:.2f} deg  hip_amp = {m['hip_amp']:.2f} deg",
      flush=True)
print(f"duty(RG_E_r) = {m['duty']:.3f}  burst_r = {m['burst_r']}", flush=True)
print(f"bal_fell = {m['bal_fell']}  bal_tilt_max = {m['bal_tilt_max']:.2f}",
      flush=True)
sel = {}
for kk in ("contact_frac_r", "contact_frac_l", "n_cycles_r", "n_cycles_l",
           "bilateral", "knee_min", "hip_amp", "duty_r", "duty_l",
           "cadence", "kine_score"):
    if kk in k:
        sel[kk] = k[kk]
print(f"kine fields: {json.dumps(sel)}", flush=True)

# per-joint ranges from the npz the eval wrote
import numpy as np
z = np.load(SPINAL / os.environ["AARL_NPZ"], allow_pickle=True)
q, names = z["q"], list(z["key_joints"])
print("joint ranges over the full run [deg]:", flush=True)
for jn in ("hip_flexion_r", "knee_angle_r", "ankle_angle_r",
           "hip_flexion_l", "knee_angle_l", "ankle_angle_l",
           "pelvis_tilt"):
    if jn in names:
        col = q[:, names.index(jn)]
        print(f"  {jn:16s} {col.min():8.2f} .. {col.max():8.2f}", flush=True)
print("BASELINE_DONE", flush=True)
