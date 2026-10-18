"""Write the 2026-10-16 drive-sweep winner param json (JSON-RULE
conformant: the full param set the run used, so a future loader branch
`P.G[k] = float(p.get(k, 0.0))` consumes it directly).

Winner: g1.5_a1_d0 (gamma_l=gamma_r=1.5, alpha_l=alpha_r=1.0) on the
anchor-A base (s3k trial-34 winner + full_rules=0).
Also writes the runner-up g2.5_a1_d0 for reference.
"""
import json

v10 = json.loads(open("best_walk_params_v10.json", encoding="utf-8").read())
base = dict(v10["multipliers"])
base["renshaw"] = 0.5
s3k = json.loads(open("reports_20260923/s3k_trial34_full_params.json",
                      encoding="utf-8").read())["params"]
base.update(s3k)
base["full_rules"] = 0.0            # anchor A recipe (sweep base)

for ident, gl, gr, al, ar in (
        ("winner_g1.5_a1_d0", 1.5, 1.5, 1.0, 1.0),
        ("runnerup_g2.5_a1_d0", 2.5, 2.5, 1.0, 1.0)):
    p = dict(base)
    p["gamma_l"], p["gamma_r"] = gl, gr
    p["alpha_l"], p["alpha_r"] = al, ar
    p["alpha_to_mn"] = 0.5          # fixed ALPHA->MN conductance (params.G)
    doc = dict(
        id=ident,
        note=("2026-10-16 human L/R alpha+gamma drive sweep winner; "
              "base = anchor A (s3k trial-34 + full_rules=0); loader: "
              "set_stage(4, params) then P.G[k]=float(p[k]) for "
              "gamma_l/gamma_r/alpha_l/alpha_r before build; "
              "reproduce: drive_sweep_full_20261016.py"),
        params=p,
        knobs=dict(gamma_l=gl, gamma_r=gr, alpha_l=al, alpha_r=ar),
    )
    out = f"reports_20261016/drive_sweep_human/{ident}.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=1)
    print("wrote", out)
