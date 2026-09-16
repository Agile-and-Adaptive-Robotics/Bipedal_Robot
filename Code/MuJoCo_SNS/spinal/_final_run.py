"""Final v11 run: apply curriculum stage-3 winner inline (no curriculum
import — avoids the closed-stdout issue), full 22 s."""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import params as P
import runner as R

st3 = json.loads(open("curriculum_stage3.json", encoding="utf-8").read())
prev = json.loads(open("best_walk_params_v10.json",
                       encoding="utf-8").read())
fit = json.loads(open("fitted_walk_params.json", encoding="utf-8").read())
BASE = dict(
    e2_pf=P.W_PF_MN["E2"]["ankle_pf"], f1_df=P.W_PF_MN["F1"]["ankle_df"],
    f1_kf=P.W_PF_MN["F1"]["knee_flex"],
    post_kneext=P.W_POSTURE["knee_ext"], post_hipext=P.W_POSTURE["hip_ext"],
    desc_f=P.G["descend_to_rg_f"], e2_adapt=P.PF_SHAPE["E2"][1])

for ph, tbl in fit["W_PF_MN"].items():
    for g, w in tbl.items():
        P.W_PF_MN[ph][g] = float(w)
for g, w in fit["W_POSTURE"].items():
    P.W_POSTURE[g] = float(w)

mul = prev["multipliers"]
gain = float(prev["pf_gain"])
for ph in P.W_PF_MN:
    for g in P.W_PF_MN[ph]:
        P.W_PF_MN[ph][g] *= gain
for g in P.W_POSTURE:
    P.W_POSTURE[g] *= gain
P.TAU["rg_adapt"] = st3["params"]["rg_adapt"]
P.G["descend_to_rg_e"] = st3["params"]["desc_e"]
P.G["descend_to_rg_f"] = st3["params"]["desc_f"]
P.G["rg_to_pf"] = st3["params"]["rg_to_pf"]
P.PF_SHAPE["E2"] = (P.PF_SHAPE["E2"][0], 1.7)  # fitted default
P.G["phase_reset_e"] = float(st3["params"]["phase_reset_e"])
P.G["phase_reset_f"] = float(st3["params"]["phase_reset_f"])
P.G["heel_rge"] = float(st3["params"]["heel_rge"])
P.G["toe_rge"] = float(st3["params"]["toe_rge"])
P.G["ib_rge"] = float(st3["params"]["ib_rge"])
P.G["ia_in"] = float(st3["params"]["ia_in"])
P.G["renshaw"] = 0.5
P.G["ankle_post_walk_trim"] = float(st3["params"].get(
    "ankle_post_walk_trim", 0.63))
P.BAL["kx"] = 150.0
drive = float(st3["params"]["drive"])

print(f"v11 curriculum final: drive {drive:.3f}, heel {P.G['heel_rge']:.2f}, "
      f"toe {P.G['toe_rge']:.2f}, ib {P.G['ib_rge']:.2f}, "
      f"ia_in {P.G['ia_in']:.2f}, renshaw 0.5, trim "
      f"{P.G['ankle_post_walk_trim']:.2f}", flush=True)
R.main(["--drive", repr(drive), "--time", "22"])
import shutil
shutil.copy("spinal_run.npz", "curriculum_final.npz")
print("saved curriculum_final.npz")
