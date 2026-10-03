import csv
from collections import defaultdict

p = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\gait_validation_20260930.csv"
rows = list(csv.DictReader(open(p)))
grp = defaultdict(list)
for r in rows:
    grp[r["variant"]].append((r["reference"], float(r["kine_score"])))
for v, rs in grp.items():
    self_score = [s for n, s in rs if n == "kine_score"][0]
    refs = [s for n, s in rs if n != "kine_score"]
    print("%s: self=%.2f  n_refs=%d  mean=%.1f  min=%.1f (%s)  max=%.1f (%s)" % (
        v, self_score, len(refs), sum(refs) / len(refs),
        min(refs), min(rs, key=lambda x: x[1])[0],
        max(refs), max(rs, key=lambda x: x[1])[0]))
    ong = [s for n, s in rs if n.startswith("ong_speed")]
    run = [s for n, s in rs if "_Run_" in n]
    print("   ong_speed family mean=%.1f (n=%d)  arnold Run family mean=%.1f (n=%d)" % (
        sum(ong) / len(ong), len(ong), sum(run) / len(run), len(run)))
