import csv

path = "tmp/paper_bridge/capability_table_20261022/per_muscle_capability_20261023.csv"
rows = list(csv.DictReader(open(path, newline="", encoding="utf-8")))
better = []
worse = []
equal = []
for r in rows:
    hand = float(r["hand_meet_pct"])
    opt = float(r["opt_meet_pct"])
    name = r["muscle"]
    if opt > hand:
        better.append((name, hand, opt))
    elif opt < hand:
        worse.append((name, hand, opt))
    else:
        equal.append(name)
print(f"optimized strictly better (meet%): {len(better)} of 27")
for n, h, o in sorted(better, key=lambda t: t[2] - t[1], reverse=True):
    print(f"  {n}: hand {h:.0f} -> opt {o:.0f}")
print(f"optimized strictly worse: {len(worse)} of 27")
for n, h, o in sorted(worse, key=lambda t: t[1] - t[2], reverse=True):
    print(f"  {n}: hand {h:.0f} -> opt {o:.0f}")
print(f"ties: {len(equal)} {equal}")
