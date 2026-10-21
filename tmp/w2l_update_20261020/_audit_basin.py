import csv
from collections import Counter

p = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261020\w2l_contact2\verify_trials.csv"
rows = list(csv.DictReader(open(p, encoding="utf-8")))
print("total rows:", len(rows))
blocks = [
    ("early(04:0x-04:40)", 0, 9),
    ("r030_basin(te5.27)", 9, 21),
    ("r100_basin(te4.35)", 21, 33),
    ("ben_basin", 33, 45),
    ("robust_basin", 45, 57),
    ("comm2_basin", 58 - 1, 70 - 1),
    ("comm4_basin", 70 - 1, 82 - 1),
    ("contra05_basin", 82 - 1, 94 - 1),
]
for name, a, b in blocks:
    blk = rows[a:b]
    c = Counter(r["survive"] for r in blk)
    print(f"{name:22s} n={len(blk):2d} survive_True={c.get('True',0)} survive_False={c.get('False',0)}")
