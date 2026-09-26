"""Harvest check for w2lvar stage 4: trial distribution + winner repro
through the exact stage-4 objective path (--no-ground --eval)."""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
os.chdir(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, os.getcwd())

import optuna

st = optuna.load_study(study_name="curr_w2lvar_s4_walk_nocontact",
                       storage="sqlite:///optuna_w2lvar.db")
vals = [t.value for t in st.trials]
print(f"n_trials={len(st.trials)} "
      f"COMPLETE={[t.state.name for t in st.trials].count('COMPLETE')}")
print("values:", [None if v is None else round(v, 3) for v in vals])
nan400 = sum(1 for v in vals if v is not None and v <= -399.5)
frozen = sum(1 for v in vals if v is not None and -320.5 <= v <= -319.5)
clip = sum(1 for v in vals if v is not None and -315.5 <= v <= -314.5)
above = sum(1 for v in vals if v is not None and v > -315.0)
print(f"nan(-400)={nan400}  frozen(-320)={frozen}  clip(-315)={clip}  "
      f"real-kine(>-315)={above}")

win = json.loads(open("curriculum_w2lvar_stage4.json",
                      encoding="utf-8").read())
p = win["params"]
print("winner params:", p)

# EXACT stage-4 objective reproduction (code path from
# _curriculum_w2lvar.py set_stage(4, p) + objective() stages-4/5 branch).
import optuna_walk_v10 as OW
import params as P
import runner as R

os.environ["AARL_NET"] = "w2lvar"
os.environ["AARL_NPZ"] = "spinal_run_w2lvar.npz"
os.environ.pop("AARL_KY", None)
os.environ.pop("AARL_PELVIS_TY", None)
BASE = dict(json.loads(open("best_walk_params_v10.json",
                            encoding="utf-8").read())["multipliers"])
BASE["renshaw"] = 0.5
full = {**BASE, **p}
OW.load_fitted_baseline()
OW.set_params(full)
P.TAU["rg_nap_h"] = float(p["rg_nap_h"])
P.G["ia_to_mn"] = float(p["ia_to_mn"])
P.G["ia_to_antagonist"] = float(p["ia_to_antagonist"])
P.G["ii_to_mn"] = float(p["ii_to_mn"])
P.G["ib_to_mn_inh"] = float(p["ib_to_mn_inh"])
P.G["rg_mutual_inh"] = float(p["rg_mutual_inh"])
P.G["rg_weak_exc"] = float(p["rg_weak_exc"])
P.G["pf_to_mn"] = float(p["pf_to_mn"])

m = R.main(["--no-ground", "--eval", "--drive", repr(p["drive"])])
if m["nan"]:
    score, extra = -400.0, "NAN"
elif m.get("kine") is None:
    score, extra = -320.0, "NO CYCLES (frozen sentinel)"
else:
    score = max(float(m["kine_score"]), -315.0)
    if float(m["tilt_max"]) > 40.0:
        score -= 10.0
    extra = (f"kine_score={float(m['kine_score']):.3f} "
             f"tilt_max={float(m['tilt_max']):.2f} "
             f"duty={m.get('duty')} cadence={m.get('cadence')} "
             f"cycles={m.get('cycles')}")
print(f"REPRO: {extra} recomputed_score={score:.6f} "
      f"json_score={win['score']:.6f} "
      f"delta={abs(score - win['score']):.6f}")
print("exploit gate:",
      "PASS (real cycle pattern, not on a sentinel floor)"
      if score > -320.0 and abs(score - win["score"]) < 0.5 else "FAIL")
