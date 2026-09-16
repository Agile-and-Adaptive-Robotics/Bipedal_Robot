"""Test swing-phase PF suppression (KINH->ankle_pf, f1_anklepf_inh) for
air-stepping dorsiflexion. Reuses the exact v6+limits+Renshaw config.
"""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

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

CASES = [("anh0", 0.0), ("anh05", 0.5), ("anh10", 1.0)]
for name, anh in CASES:
    P.G["f1_anklepf_inh"] = anh
    R.main(["--no-ground", "--no-afferents", "--no-interleg",
            "--time", "14", "--drive", repr(DRIVE)])
    z = np.load(R.HERE / "spinal_run.npz", allow_pickle=True)
    t, q = z["t"], z["q"]
    jj = [str(a) for a in z["key_joints"]]
    ji, jk = jj.index("ankle_angle_r"), jj.index("knee_angle_r")
    m = (t >= 5.0) & (t <= 17.0)
    print(f"RESULT {name}: ankle {q[m, ji].min():+.1f}..{q[m, ji].max():+.1f} "
          f"deg (range {q[m, ji].max() - q[m, ji].min():.1f}, ref 23.1; "
          f"dorsi = positive) | knee range "
          f"{q[m, jk].max() - q[m, jk].min():.1f} deg", flush=True)
