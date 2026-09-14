"""Ankle dorsiflexion experiment (air-stepping): strengthen swing DF
drive and trim push-off PF, report ankle range. Usage: _ankle_test.py.
"""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

# load fitted + best6 exactly like the flags do, so mutations start from
# the tuned config
import runner as R
import params as P

with open(R.HERE / "fitted_walk_params.json", encoding="utf-8") as f:
    fit = json.load(f)
for ph, tbl in fit["W_PF_MN"].items():
    for g, w in tbl.items():
        P.W_PF_MN[ph][g] = float(w)
for g, w in fit["W_POSTURE"].items():
    P.W_POSTURE[g] = float(w)
best = json.loads((R.HERE / "best_walk_params_v6.json").read_text("utf-8"))
bp = best["params"]
gain = float(best["pf_gain"])
for ph in P.W_PF_MN:
    for g in P.W_PF_MN[ph]:
        P.W_PF_MN[ph][g] *= gain
for g in ("e2_pf", "f1_df", "f1_kf"):
    key = {"e2_pf": ("E2", "ankle_pf"), "f1_df": ("F1", "ankle_df"),
           "f1_kf": ("F1", "knee_flex")}[g]
    P.W_PF_MN[key[0]][key[1]] = bp[g]
P.W_POSTURE["knee_ext"] = bp["post_kneext"]
P.W_POSTURE["hip_ext"] = bp["post_hipext"]
P.TAU["rg_adapt"] = bp["rg_adapt"]
P.G["descend_to_rg_e"] = bp["desc_e"]
P.G["descend_to_rg_f"] = bp["desc_f"]
P.G["rg_to_pf"] = bp["rg_to_pf"]
P.PF_SHAPE["E2"] = (P.PF_SHAPE["E2"][0], bp["e2_adapt"])
P.G["phase_reset_e"] = bp.get("phase_reset_e", 0.0)
P.G["phase_reset_f"] = bp.get("phase_reset_f", 0.0)
P.G["f1_kneext_inh"] = bp.get("f1_kneext_inh", 0.0)
P.G["renshaw"] = 0.5
P.BAL["kx"] = bp["kx"]
DRIVE = bp["drive"]

CASES = [
    ("baseline", {}),
    ("df_boost", {"F1_ankle_df": 3.0, "F2_ankle_df": 0.30,
                  "E2_ankle_pf": 0.7}),
    ("df_boost2", {"F1_ankle_df": 5.0, "F2_ankle_df": 0.45,
                   "E2_ankle_pf": 0.5, "tib_ant_ii": 1.0}),
]

from muscle_map import classify  # noqa: E402


def apply(muts):
    # re-pin the three tuned knobs, then apply mutations
    P.W_PF_MN["E2"]["ankle_pf"] = bp["e2_pf"]
    P.W_PF_MN["F1"]["ankle_df"] = bp["f1_df"]
    P.W_PF_MN["F1"]["knee_flex"] = bp["f1_kf"]
    for k, v in muts.items():
        if k == "F1_ankle_df":
            P.W_PF_MN["F1"]["ankle_df"] = bp["f1_df"] * v
        elif k == "F2_ankle_df":
            P.W_PF_MN["F2"]["ankle_df"] = v
        elif k == "E2_ankle_pf":
            P.W_PF_MN["E2"]["ankle_pf"] = bp["e2_pf"] * v
        elif k == "tib_ant_ii":
            P.AFF["ii_gain"] = P.AFF["ii_gain"]  # unused marker


for name, muts in CASES:
    apply(muts)
    R.main(["--no-ground", "--no-afferents", "--no-interleg",
            "--time", "14", "--drive", repr(DRIVE)])
    z = np.load(R.HERE / "spinal_run.npz", allow_pickle=True)
    t, q = z["t"], z["q"]
    ji = [str(a) for a in z["key_joints"]].index("ankle_angle_r")
    jk = [str(a) for a in z["key_joints"]].index("knee_angle_r")
    m = (t >= 5.0) & (t <= 17.0)
    ar = float(np.max(q[m, ji]) - np.min(q[m, ji]))
    amx = float(np.max(q[m, ji]))
    kr = float(np.max(q[m, jk]) - np.min(q[m, jk]))
    print(f"RESULT {name}: ankle range {ar:.1f} deg (max {amx:+.1f}, "
          f"ref range 23.1, dorsi = positive) | knee range {kr:.1f} deg",
          flush=True)
