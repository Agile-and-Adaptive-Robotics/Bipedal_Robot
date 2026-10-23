import csv

base = r"D:\GitHub\Bipedal_Robot\tmp\paper_bridge\capability_table_20261022"
rows = list(csv.DictReader(open(base + r"\per_muscle_capability_20261023.csv", encoding="utf-8")))
inc = [r for r in rows if r["excluded"].strip().lower() == "no"]
print("total rows:", len(rows), " included:", len(inc))


def mean(col):
    return sum(float(r[col]) for r in inc) / len(inc)


for c in ["hand_meet_pct", "hand_active_pct", "opt_meet_pct", "opt_active_pct"]:
    print(c, round(mean(c), 2))

meet_better = sum(1 for r in inc if float(r["opt_meet_pct"]) > float(r["hand_meet_pct"]))
meet_tie = sum(1 for r in inc if float(r["opt_meet_pct"]) == float(r["hand_meet_pct"]))
act_better = sum(1 for r in inc if float(r["opt_active_pct"]) > float(r["hand_active_pct"]))
both = sum(
    1
    for r in inc
    if float(r["opt_meet_pct"]) > float(r["hand_meet_pct"])
    and float(r["opt_active_pct"]) > float(r["hand_active_pct"])
)
print("meet strictly larger:", meet_better, " ties:", meet_tie, " lower:", len(inc) - meet_better - meet_tie)
print("active strictly larger:", act_better)
print("improved on BOTH measures:", both)
