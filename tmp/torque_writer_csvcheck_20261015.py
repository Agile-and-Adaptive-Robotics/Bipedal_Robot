import csv, os
from collections import OrderedDict

root = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization"

for tag, p in [("v1", os.path.join(root, "Routing27", "Results", "routing27_summary_full.csv")),
               ("v2", os.path.join(root, "Routing27_v2", "Results", "routing27_summary_full.csv"))]:
    with open(p, newline='') as f:
        rows = list(csv.DictReader(f))
    ids = OrderedDict()
    for r in rows:
        ids.setdefault(r["actuatorId"], []).append(r)
    chosen = [(r["actuatorId"], r["masterName"], r["config"], r["feasible"], r["marginWorstRel"])
              for r in rows if r["chosen"] == "1"]
    errored = sorted({k for k, v in ids.items()
                      if all(x["exitflagS"] == "error" for x in v)}, key=int)
    chosen_feasible = [c for c in chosen if c[3] == "1"]
    chosen_posmargin = [c for c in chosen if c[3] == "1" and float(c[4]) > 0]
    print(f"--- {tag} ({os.path.basename(os.path.dirname(os.path.dirname(p)))}) ---")
    print(f"total rows: {len(rows)}; distinct actuators present: {len(ids)}; ids: {sorted(ids, key=int)}")
    print(f"actuators with a chosen=1 row: {len(chosen)}")
    print(f"error actuators (all rows error): {len(errored)} -> ids {errored}")
    print(f"chosen rows with feasible=1: {len(chosen_feasible)}")
    for c in chosen_feasible:
        print("   feasible-chosen:", c)
    print(f"chosen+feasible+positive margin: {len(chosen_posmargin)}")
    for c in chosen_posmargin:
        print("   posmargin-chosen:", c)
    print()
