"""Assemble the per-model standing-balance verdict table (2026-10-21,
resume session) from the on-disk artifacts. Numbers are read from the
JSONs written by the studies/verifiers - no hand-typed values except
file labels. Emits verdict_table_20261021.json + .md in the balance
report folder."""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
B = os.path.join("reports_20261020", "balance")
S = "."


def jload(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)


out = {}

# --- quiet stage-3 winners: DB pull + fresh verify + Ib ablation ------
quiet = jload(os.path.join(S, "_s3_pull_20261020.json"))
stats = jload(os.path.join(S, "_s3_stats_20261021.json"))
rig = jload(os.path.join(S, "_s3_rigshare_20261021.json"))
ver = {c: jload(os.path.join(S, f"_s3_verify_{c}_20261020.json"))
       for c in ("stock", "syn6_bal", "w2l_bal", "w2l_fix")}
abl = {c: jload(os.path.join(S, f"_s3_ibablate_{c}_20261020.json"))
       for c in ("stock", "syn6_bal", "w2l_bal", "w2l_fix")}

# --- s3p perturbation winners: saved winner jsons + fresh verify -----
s3p = {}
for fam in ("stock", "syn6", "w2lvar"):
    s3p[fam] = jload(os.path.join(S, f"curriculum_s3p_{fam}.json"))
    try:
        s3p[fam]["fresh"] = jload(
            os.path.join(S, f"_s3p_verify_{fam}_20261021.json"))
    except FileNotFoundError:
        s3p[fam]["fresh"] = None

out["quiet_winners"] = {
    "stock": dict(db=quiet["sqlite:///optuna_walk.db::curr_s3_balance"],
                  fresh_score=ver["stock"]["recomputed_score"],
                  db_value=ver["stock"]["db_value"]),
    "w2lvar": dict(db=quiet["sqlite:///optuna_w2lvar.db::curr_w2lvar_s3_balance"],
                   fresh_score=ver["w2l_bal"]["recomputed_score"],
                   db_value=ver["w2l_bal"]["db_value"]),
    "syn6": dict(db=quiet["sqlite:///optuna_syn6.db::curr_syn6_s3_balance"],
                 fresh_score=ver["syn6_bal"]["recomputed_score"],
                 db_value=ver["syn6_bal"]["db_value"]),
}
out["s3p_winners"] = s3p
out["rig_share"] = rig
out["ib_ablation"] = {c: dict(fell=a["metrics"]["bal_fell"],
                              tilt=a["metrics"]["bal_tilt_max"],
                              sway=a["metrics"]["bal_sway"])
                      for c, a in abl.items()}
out["s3p_study_counts"] = {
    k: dict(n_total=stats[k]["n_total"], n_complete=stats[k]["n_complete"],
            hist=stats[k]["hist"])
    for k in stats if "s3p" in k or k.endswith("s3_balance")}

with open(os.path.join(B, "verdict_table_20261021.json"), "w",
          encoding="utf-8") as f:
    json.dump(out, f, indent=1, default=str)
print("wrote verdict_table_20261021.json")
