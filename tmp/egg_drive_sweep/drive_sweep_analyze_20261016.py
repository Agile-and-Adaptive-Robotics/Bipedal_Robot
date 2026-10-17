"""Analyze the 2026-10-16 human L/R alpha+gamma drive sweep.

Reads reports_20261016/drive_sweep_human/results.jsonl, prints a table,
writes summary CSV next to it, and applies the campaign walking bar:
  WALKS iff not nan AND (cycles_r + cycles_l) >= 6 AND contact duty
  (mean of per-leg duties where defined) in [0.4, 0.75] AND not bal_fell.
Per-leg detail is reported alongside (one-leg-dominant gaits are
explicitly flagged, not hidden in the total).
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "reports_20261016", "drive_sweep_human",
                   "results.jsonl")
rows = []
for ln in open(SRC, encoding="utf-8"):
    try:
        rows.append(json.loads(ln))
    except Exception:
        pass
# keep the LAST record per id (re-runs append)
last = {}
for r in rows:
    last[r["spec"]["id"]] = r
rows = [last[k] for k in sorted(last)]

out_rows = []
for r in rows:
    s = r["spec"]
    duties = [d for d in (r.get("duty_r"), r.get("duty_l"))
              if d is not None]
    duty = sum(duties) / len(duties) if duties else None
    cyc = (r.get("cycles_r", 0) or 0) + (r.get("cycles_l", 0) or 0)
    walks = (not r.get("nan", True)
             and cyc >= 6
             and duty is not None and 0.4 <= duty <= 0.75
             and not r.get("bal_fell", True))
    out_rows.append(dict(
        id=s["id"], gamma_l=s["gamma_l"], gamma_r=s["gamma_r"],
        alpha_l=s["alpha_l"], alpha_r=s["alpha_r"],
        kine_score=r.get("kine_score"), nan=r.get("nan"),
        cycles_r=r.get("cycles_r"), cycles_l=r.get("cycles_l"),
        duty_r=r.get("duty_r"), duty_l=r.get("duty_l"),
        duty_mean=duty, ds=r.get("ds"),
        T_r=r.get("T_r"), T_l=r.get("T_l"),
        contact_frac_r=r.get("contact_frac_r"),
        contact_frac_l=r.get("contact_frac_l"),
        tilt_max=r.get("tilt_max"), knee_min=r.get("knee_min"),
        dx=r.get("dx"), bal_com_z_min=r.get("bal_com_z_min"),
        bal_fell=r.get("bal_fell"), WALKS=walks))

hdr = ("id", "gamma_l", "gamma_r", "alpha_l", "alpha_r", "kine_score",
       "cyc_r", "cyc_l", "duty_r", "duty_l", "duty_mean", "nan", "fell",
       "tilt", "knee_min", "dx", "WALKS")
print(("%-16s %5s %5s %5s %5s %12s %5s %5s %6s %6s %8s %5s %5s "
       "%6s %8s %7s %5s") % hdr)
for o in out_rows:
    print(("%-16s %5.2f %5.2f %5.2f %5.2f %12.4f %5d %5d %6s %6s %8s "
           "%5s %5s %6.1f %8.1f %7s %5s") % (
        o["id"], o["gamma_l"], o["gamma_r"], o["alpha_l"], o["alpha_r"],
        o["kine_score"] if o["kine_score"] is not None else float("nan"),
        o["cycles_r"] or 0, o["cycles_l"] or 0,
        "%.2f" % o["duty_r"] if o["duty_r"] is not None else "-",
        "%.2f" % o["duty_l"] if o["duty_l"] is not None else "-",
        "%.2f" % o["duty_mean"] if o["duty_mean"] is not None else "-",
        o["nan"], o["bal_fell"], o["tilt_max"], o["knee_min"],
        "%.2f" % o["dx"] if o["dx"] is not None else "-",
        o["WALKS"]))

csvp = SRC.replace("results.jsonl", "summary.csv")
with open(csvp, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    w.writerows(out_rows)
print(f"\nwrote {csvp}  ({len(out_rows)} rows)")

walkers = [o for o in out_rows if o["WALKS"]]
print(f"\nWALKERS (>=6 cycles, duty 0.4-0.75, no NaN, no fall): "
      f"{len(walkers)}")
for o in sorted(walkers, key=lambda x: -(x["kine_score"] or -1e9)):
    print(f"  {o['id']}: kine {o['kine_score']:.2f}, cycles "
          f"{o['cycles_r']}r+{o['cycles_l']}l, duty {o['duty_mean']:.2f}")
