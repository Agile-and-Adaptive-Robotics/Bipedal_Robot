"""Harvest check for w2lvar stage 3: trial distribution + winner repro
through the exact stage-3 objective path (--stand-eval 8 --rig-scale S)."""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
os.chdir(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, os.getcwd())

import optuna

st = optuna.load_study(study_name="curr_w2lvar_s3_balance",
                       storage="sqlite:///optuna_w2lvar.db")
vals = [t.value for t in st.trials]
print(f"n_trials={len(st.trials)} "
      f"COMPLETE={[t.state.name for t in st.trials].count('COMPLETE')}")
print("values:", [None if v is None else round(v, 3) for v in vals])
nan = sum(1 for v in vals if v is not None and v <= -199.5)
fall = sum(1 for v in vals if v is not None and -199.5 < v <= -149.5)
stander = sum(1 for v in vals if v is not None and v > -149.5)
print(f"nan(-200)={nan}  fall(-150)={fall}  standers(>-149.5)={stander}")

win = json.loads(open("curriculum_w2lvar_stage3.json",
                      encoding="utf-8").read())
p = win["params"]
print("winner params:", p)

# EXACT stage-3 objective reproduction (code path from
# _curriculum_w2lvar.py set_stage(3, p) + objective() stage-3 branch).
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
P.G["vest_prop"] = float(p["vest_prop"])

m = R.main(["--stand-eval", "8", "--rig-scale",
            repr(float(p["rig_scale"]))])
if m["nan"]:
    score = -200.0
    extra = "NAN"
elif m["bal_fell"]:
    score = -150.0
    extra = "FELL"
else:
    score = (100.0
             - 400.0 * float(m["bal_sway"])
             - 1.0 * float(m["bal_tilt_max"])
             - 40.0 * abs(0.5 - float(m["bal_contact_sym"])))
    extra = (f"sway={float(m['bal_sway']):.4f} "
             f"tilt_max={float(m['bal_tilt_max']):.3f} "
             f"contact_sym={float(m['bal_contact_sym']):.4f} "
             f"com_floor_min={m.get('kz')}")
print(f"REPRO: {extra} recomputed_score={score:.6f} "
      f"json_score={win['score']:.6f} "
      f"delta={abs(score - win['score']):.6f}")
print("exploit gate:", "PASS (genuine stander, not on a sentinel floor)"
      if score > -149.5 and abs(score - win["score"]) < 0.5
      else "FAIL")
