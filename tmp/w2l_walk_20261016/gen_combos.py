import json
import sys

stage = sys.argv[1]
combos = []

if stage == "a":      # symmetric strong-drive grid (protocol defaults)
    i = 0
    for te in (2.0, 3.0, 4.0, 5.0, 6.0):
        for tf in (3.0, 4.0, 5.0, 6.0, 7.0):
            i += 1
            combos.append({"id": f"a{i:02d}", "args": {"te": te, "tf": tf}})
elif stage == "b":    # asymmetric GAMMA around the symmetric best
    # filled in after reading stage A results
    raise SystemExit("edit me after stage A")
else:
    raise SystemExit(f"unknown stage {stage}")

json.dump(combos, open(f"combos_{stage}.json", "w"), indent=1)
print(len(combos), "combos ->", f"combos_{stage}.json")
