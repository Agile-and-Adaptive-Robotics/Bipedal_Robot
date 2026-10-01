"""Dump per-Louvain-cluster composition so clusters can be named (2026-09-28).

Writes export/cluster_profiles.json: per cluster — size, year span, top
primary authors, top title keywords (stopworded), top animals, and 12 sample
titles. A human (or the supervisor session) reads this and writes
export/cluster_labels.json {clusterIndex: "Label"}, consumed by build_app.py.
"""
import json, os, re
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
EXPORT = os.path.join(HERE, "export")
records = json.load(open(os.path.join(EXPORT, "sadb_export.json"), encoding="utf-8"))
layout = json.load(open(os.path.join(EXPORT, "sadb_layout.json"), encoding="utf-8"))

STOP = set("""a an and the of in on to for with without by from during via their its his her
our between within into over under both either neither or nor is are was were be been
using used use how what which who whom whose that this these those it its as at than
then so such can could may might will would shall should must not no nor only also more
most less least very much many few new novel role roles effect effects study studies
review approaches approach mechanism mechanisms control controls system systems
""".split())

byid = {r["id"]: r for r in records}
clusters = {}
for rid, L in layout.items():
    cl = L.get("cl", -1)
    if cl >= 0:
        clusters.setdefault(cl, []).append(rid)

prof = {}
for cl, ids in sorted(clusters.items()):
    rs = [byid[i] for i in ids if i in byid]
    yrs = [r["year"] for r in rs if r["year"]]
    authors = Counter(r["primary"] for r in rs if r["primary"])
    animals = Counter(a for r in rs for a in r["animals"])
    kw = Counter()
    for r in rs:
        for w in re.findall(r"[a-z]{4,}", r["title"].lower()):
            if w not in STOP:
                kw[w] += 1
    prof[cl] = {
        "size": len(rs),
        "years": [min(yrs), max(yrs)] if yrs else [],
        "top_authors": authors.most_common(8),
        "top_keywords": kw.most_common(15),
        "top_animals": animals.most_common(6),
        "sample_titles": [r["title"][:90] for r in rs[:12]],
    }

json.dump(prof, open(os.path.join(EXPORT, "cluster_profiles.json"), "w"),
          ensure_ascii=False, indent=1)
print("clusters:", len(prof), "| sizes:", {k: v["size"] for k, v in prof.items()})
