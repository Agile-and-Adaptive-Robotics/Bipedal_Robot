import csv

rows = list(csv.DictReader(open("tmp/paper_bridge/capability_table_20261022/per_muscle_capability_20261023.csv", newline="", encoding="utf-8")))
kept = [r for r in rows if r.get("excluded") != "yes"]
print(f"kept {len(kept)}")
meet_better = [r for r in kept if float(r["opt_meet_pct"]) > float(r["hand_meet_pct"])]
act_better = [r for r in kept if float(r["opt_active_pct"]) > float(r["hand_active_pct"])]
both = [r for r in kept if float(r["opt_meet_pct"]) > float(r["hand_meet_pct"]) and float(r["opt_active_pct"]) > float(r["hand_active_pct"])]
either = [r for r in kept if float(r["opt_meet_pct"]) > float(r["hand_meet_pct"]) or float(r["opt_active_pct"]) > float(r["hand_active_pct"])]
meet_tie = [r for r in kept if float(r["opt_meet_pct"]) == float(r["hand_meet_pct"])]
print(f"torque measure (meet share): better {len(meet_better)}, tie {len(meet_tie)}")
print(f"RoM measure (active share): better {len(act_better)}")
print(f"better on BOTH: {len(both)} -> {[r['muscle'] for r in both]}")
print(f"better on EITHER: {len(either)}")
print(f"active-detail rows (opt_act vs hand_act, kept):")
for r in kept:
    print(f"  {r['muscle']}: meet {float(r['hand_meet_pct']):.0f}->{float(r['opt_meet_pct']):.0f}  active {float(r['hand_active_pct']):.0f}->{float(r['opt_active_pct']):.0f}")
