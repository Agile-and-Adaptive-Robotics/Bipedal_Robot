import csv
import json
import os

ROOT = r"D:\Github\Bipedal_Robot\SADb_audit"

# 1) chapter arithmetic: post-prune 957 + removed 436
post = {"papers": 561, "pathways": 39, "models": 244, "reviews": 109, "recorders": 4}
removed = {"papers": 382, "orphans": 54}
post_total = sum(post.values())
print("post-prune base:", post_total, "(chapter: 957)")
print("pre-prune base = 957 + 436 =", post_total + removed["papers"] + removed["orphans"])
print("pre-prune papers = 561 + 382 =", 561 + 382, "(audited corpus 943)")

# 2) deletion manifest
mpath = os.path.join(ROOT, "deletion_manifest.json")
if os.path.exists(mpath):
    d = json.load(open(mpath, encoding="utf-8"))
    if isinstance(d, list):
        print("deletion_manifest entries:", len(d))
    else:
        print("deletion_manifest keys:", list(d.keys())[:10])
        for k in ("total", "count", "records"):
            if k in d:
                print(" ", k, "=", d[k])

# 3) deletion log rows
lpath = os.path.join(ROOT, "deletion_log.csv")
if os.path.exists(lpath):
    rows = list(csv.DictReader(open(lpath, encoding="utf-8-sig")))
    print("deletion_log.csv rows:", len(rows))
    if rows:
        print("  columns:", list(rows[0].keys())[:8])

# 4) task1 purge log
tpath = os.path.join(ROOT, "task1_purge_log.csv")
if os.path.exists(tpath):
    rows = list(csv.DictReader(open(tpath, encoding="utf-8-sig")))
    print("task1_purge_log.csv rows:", len(rows))

# 5) prune review docs mention of totals
for f in ("prune_review_20261001.md", "prune_review_round2.md"):
    p = os.path.join(ROOT, f)
    if os.path.exists(p):
        txt = open(p, encoding="utf-8", errors="replace").read()
        for kw in ("1393", "1,393", "1357", "1,357", "1300", "1,300", "pre-prune", "before the prune", "957", "1,000"):
            if kw in txt:
                i = txt.find(kw)
                print(f, "::", kw, "::", txt[max(0, i - 80) : i + 80].replace("\n", " "))
