"""Rig-share analysis of the stage-3 balance winners (2026-10-21).

Answers Ben's 'stands or leans on the rig' question with numbers: from
each saved standing-eval npz, the post-warmup mean TOTAL foot normal
load (log_contact/2 - heel and toe each accumulate every contact of
their foot, so the per-foot column is 2x the true normal force), the
rig-carried share = weight - foot load, pelvis height, and knee/tilt
excursions. Nominal model weight 72 kg x 9.81 = 706 N (the converted
Gait2392-robotbody mass).
"""
import io
import json
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

W_NOM = 72.0 * 9.81
CASES = [
    ("stock quiet", "s3verify_stock_20261020.npz", 4.0),
    ("stock Ib-abl", "s3ibablate_stock_20261020.npz", 4.0),
    ("w2lvar quiet", "s3verify_w2l_bal_20261020.npz", 4.0),
    ("w2lvar Ib-abl", "s3ibablate_w2l_bal_20261020.npz", 4.0),
    ("syn6 quiet", "s3verify_syn6_bal_20261020.npz", 4.0),
    ("syn6 Ib-abl", "s3ibablate_syn6_bal_20261020.npz", 4.0),
    ("stock s3p last-trial quiet (npz=final eval, t19)",
     "spinal_run_curr_s3p_20261020.npz", 4.0),
]

out = {}
for name, fnz, warm in CASES:
    try:
        z = np.load(fnz, allow_pickle=True)
    except FileNotFoundError:
        print(f"{name}: MISSING {fnz}")
        continue
    t = z["t"]
    n = int(np.sum(t > 0)) + 1
    n = min(n, len(t))
    tt = t[:n]
    i0 = int(np.searchsorted(tt, warm))
    ct = z["contact"][i0:n] if "contact" in z.files else None
    q = z["q"][i0:n]          # already degrees: [tilt, list, rot, hipf_r,
                              #  knee_r, ankle_r, ..., knee_l(9), ...]
    com = z["com"][i0:n]
    rec = {}
    if ct is not None and ct.shape[0]:
        foot = 0.5 * (ct[:, 0] + ct[:, 1])   # true total normal force (N)
        rec["foot_load_mean_N"] = float(np.mean(foot))
        rec["foot_load_pct_weight"] = float(100.0 * np.mean(foot) / W_NOM)
        rec["rig_share_pct"] = float(
            100.0 * (W_NOM - np.mean(foot)) / W_NOM)
        rec["contact_duty_R_L"] = (float((ct[:, 0] > 20).mean()),
                                   float((ct[:, 1] > 20).mean()))
    rec["pelvis_tilt_max_deg"] = float(np.max(np.abs(q[:, 0])))
    rec["knee_flex_max_deg_R_L"] = (float(-np.min(q[:, 4])),
                                    float(-np.min(q[:, 9])))
    rec["com_z_mean"] = float(np.mean(com[:, 2]))
    rec["com_z_min"] = float(np.min(com[:, 2]))
    out[name] = rec
    print(f"\n== {name} ({fnz})")
    for k, v in rec.items():
        if isinstance(v, tuple):
            print(f"   {k:22s} = {v[0]:.3f} / {v[1]:.3f}")
        else:
            print(f"   {k:22s} = {v:.2f}")

with open("_s3_rigshare_20261021.json", "w", encoding="utf-8") as f:
    json.dump(dict(weight_nominal_N=W_NOM, cases=out), f, indent=1)
print("\nwrote _s3_rigshare_20261021.json")
