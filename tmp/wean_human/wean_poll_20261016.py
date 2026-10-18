"""Poll the wean campaign rounds: rows, bilateral count, top rows."""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reports_20261016", "wean_human")

for rd in [int(a) for a in (sys.argv[1:] or ["2", "4"])]:
    rows = []
    for p in glob.glob(os.path.join(OUT, f"round{rd}_w*.jsonl")):
        rows.extend(json.loads(l) for l in open(p, encoding="utf-8"))
    bil = [r for r in rows
           if (r["metrics"].get("k_n_cycles_r") or 0) >= 3
           and (r["metrics"].get("k_n_cycles_l") or 0) >= 3]
    print(f"r{rd} n {len(rows)} bil {len(bil)} best_bil "
          f"{round(max((r['value'] for r in bil), default=-999), 2)}")
    top = sorted(rows, key=lambda x: -x["value"])[:4]
    for r in top:
        m = r["metrics"]
        print(f"  {r['value']:.2f} kine {round(m.get('kine_score') or -999, 2)}"
              f" cyc {m.get('k_n_cycles_r')}/{m.get('k_n_cycles_l')}"
              f" duty {round(m['k_duty_r'], 2) if m.get('k_duty_r') else '-'}"
              f"/{round(m['k_duty_l'], 2) if m.get('k_duty_l') else '-'}"
              f" kz {round(m.get('kz') or 0, 2)}"
              f" scale {round(r['params'].get('drive_scale', 1.0), 2)}"
              f" t{r['number']}")
