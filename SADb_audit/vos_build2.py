"""Rebuild VOSviewer map + network for the FULL Airtable corpus (2026-09-18).

Replaces the 2026-09-14 build (which covered 552 papers from local CSVs):
now reads export/sadb_export.json (all 943 records, one Airtable fetch of truth)
and re-resolves every DOI against OpenAlex. Also writes export/sadb_layout.json
(x, y, cluster, citations per record id) consumed by app/build_app.py.

VOSviewer outputs keep the established formats:
  vosviewer/sadb_map.txt       (id, label, doi, year, first_author, weight, source)
  vosviewer/sadb_network.txt   (source, target, weight)
  vosviewer/openalex_enrichment.csv / openalex_misses.txt / sadb_preview.png

Run with the myo env python (needs networkx + matplotlib).
"""
import csv, json, os, time, urllib.parse, urllib.request
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
VOS = os.path.join(HERE, "vosviewer")
EXPORT = os.path.join(HERE, "export")
MAILTO = "benjamin.bolen@pdx.edu"

records = json.load(open(os.path.join(EXPORT, "sadb_export.json"), encoding="utf-8"))
print("records:", len(records))


def norm_doi(d):
    d = (d or "").strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d


# --- source classification ---
def rec_ids_from(csv_path):
    ids = set()
    with open(csv_path, encoding="utf-8-sig") as f:
        for row in csv.reader(f):
            for cell in row:
                if cell.startswith("rec") and len(cell) == 17:
                    ids.add(cell)
    return ids


orig_ids = rec_ids_from(os.path.join(HERE, "airtable_papers_slim.csv"))
demo_ids = rec_ids_from(os.path.join(HERE, "airtable_created50_ids.csv"))
digest_ids = {"rechmhCcKA0ILQWKG", "recw5L1NKiDuwcdQY", "recp5dgFg2CMMDQAz",
              "rec8MVmUOx1oDk02j", "recB7n3thv2EEgswz", "rec1Q4YaKntrUybgY",
              "recCjTyA6Mz9tHw0G", "recZJU3Tv1NeOUDiG", "recEBU9QYruCYKihU",
              "recIjP2Mi5teZyaov", "recMZoiH13ytYevS7", "recp7yrg2a8E7DgyR",
              "recRmBChe7Lolwn0F", "recHoAYp5AUdgYEZA", "reck2RL4OlhuakV9q",
              "recafnPacssmJFhBt", "recenxdRrnk3ydqyg", "recoNlbHxaPTmt2Rs",
              "recpHeyTiV5C1veLV", "recl7iChHRGAfNn1c"}
rest_dois = set()
with open(os.path.join(HERE, "airtable_rest_import_clean.csv"), encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        d = norm_doi(row.get("DOI") or "")
        if d:
            rest_dois.add(d)


def classify(rec):
    rid = rec["id"]
    if rid in orig_ids:
        return "originals99"
    if rid in demo_ids:
        return "demo50"
    if rid in digest_ids:
        return "digest20"
    if norm_doi(rec["doi"]) in rest_dois:
        return "rest383"
    return "task5import"


for r in records:
    r["src"] = classify(r)
print("sources:", Counter(r["src"] for r in records))

# --- OpenAlex enrichment (cached in export/openalex_raw.json — only missing DOIs are fetched) ---
CACHE = os.path.join(EXPORT, "openalex_raw.json")
cache = {}
if os.path.exists(CACHE):
    cache = json.load(open(CACHE, encoding="utf-8"))
    print("cache loaded:", len(cache), "DOIs")
dois = sorted({norm_doi(r["doi"]) for r in records if norm_doi(r["doi"])})
todo = [d for d in dois if d not in cache]
print("to fetch:", len(todo))
B = 50
for i in range(0, len(todo), B):
    chunk = todo[i:i + B]
    filt = "doi:" + "|".join(chunk)
    url = ("https://api.openalex.org/works?filter=" + urllib.parse.quote(filt, safe="|")
           + "&select=id,doi,cited_by_count,referenced_works&per-page=50&mailto=" + MAILTO)
    got = 0
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(
                    url, headers={"User-Agent": "sadb-vosviewer/0.2 (mailto:%s)" % MAILTO}),
                    timeout=60) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            break
        except Exception as e:
            print("  retry", i // B, e)
            time.sleep(3 + 3 * attempt)
    else:
        print("  FAILED batch", i // B)
        continue
    for w in data.get("results", []):
        d = norm_doi(w.get("doi") or "")
        if d:
            cache[d] = {"oa_id": w["id"], "cited_by_count": w.get("cited_by_count", 0),
                        "refs": w.get("referenced_works", [])}
            got += 1
    print("openalex batch %d: %d/%d" % (i // B, got, len(chunk)))
    time.sleep(0.6)
json.dump(cache, open(CACHE, "w"))
oa = cache

resolved = {r["id"]: oa[norm_doi(r["doi"])] for r in records if norm_doi(r["doi"]) in oa}
misses = [norm_doi(r["doi"]) for r in records if norm_doi(r["doi"]) and norm_doi(r["doi"]) not in oa]
print("openalex resolved:", len(resolved), "misses:", len(misses))
oa2doi = {v["oa_id"]: k for k, v in oa.items()}

doi2rec = {}
for r in records:
    d = norm_doi(r["doi"])
    if d:
        doi2rec.setdefault(d, r["id"])


def doi2rec_get(d):
    return doi2rec.get(d)


edge_count = {}
for d, v in oa.items():
    for ref in v["refs"]:
        if ref in oa2doi and oa2doi[ref] != d:
            a, b = sorted((d, oa2doi[ref]))
            edge_count[(a, b)] = edge_count.get((a, b), 0) + 1
print("internal citation edges:", len(edge_count))

# directed citation adjacency by Airtable record id: {citerId: [citedIds]} (app "References"/"Cited by" views)
oaid2rec = {oa[d]["oa_id"]: rid for d, rid in doi2rec.items() if d in oa}
cites = {}
for d, v in oa.items():
    ra = doi2rec_get(d)
    if not ra:
        continue
    for ref in v["refs"]:
        rb = oaid2rec.get(ref)
        if rb and rb != ra:
            cites.setdefault(ra, [])
            if rb not in cites[ra]:
                cites[ra].append(rb)
json.dump(cites, open(os.path.join(EXPORT, "sadb_cites.json"), "w"))
n_edges = sum(len(v) for v in cites.values())
print("directed corpus citation pairs:", n_edges)

# --- VOSviewer map + network (formats unchanged) ---
os.makedirs(VOS, exist_ok=True)


def first_surname(s):
    a = " ".join((s or "").split())
    return a.split(" et al")[0].split(" and ")[0].split(",")[0].split(";")[0].strip()


with open(os.path.join(VOS, "sadb_map.txt"), "w", encoding="utf-8", newline="\n") as f:
    f.write("\t".join(["id", "label", "doi", "year", "first_author", "weight", "source"]) + "\n")
    for r in records:
        d = norm_doi(r["doi"])
        if d and d in oa:
            nid = oa[d]["oa_id"].rsplit("/", 1)[-1]
        elif d:
            nid = "doi:" + d
        else:
            nid = "nodoi:" + r["id"]
        f.write("\t".join([nid, (r["title"] or "")[:120], d, str(r["year"] or ""),
                           first_surname(r["primary"]),
                           str(oa[d]["cited_by_count"]) if d and d in oa else "0",
                           r["src"]]) + "\n")

with open(os.path.join(VOS, "sadb_network.txt"), "w", encoding="utf-8", newline="\n") as f:
    f.write("source\ttarget\tweight\n")
    for (a, b), c in sorted(edge_count.items(), key=lambda kv: -kv[1]):
        f.write("%s\t%s\t%d\n" % (oa[a]["oa_id"].rsplit("/", 1)[-1],
                                  oa[b]["oa_id"].rsplit("/", 1)[-1], c))

with open(os.path.join(VOS, "openalex_enrichment.csv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["doi", "oa_id", "cited_by_count", "n_refs", "n_in_corpus_refs"])
    for d, v in sorted(oa.items(), key=lambda kv: -kv[1]["cited_by_count"]):
        w.writerow([d, v["oa_id"], v["cited_by_count"], len(v["refs"]),
                    sum(1 for ref in v["refs"] if ref in oa2doi)])
with open(os.path.join(VOS, "openalex_misses.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(misses))
print("wrote vosviewer files")

# --- layout for the HTML app ---
import networkx as nx

G = nx.Graph()
doi2rec = {}
for r in records:
    d = norm_doi(r["doi"])
    if d:
        doi2rec.setdefault(d, r["id"])
for (a, b), c in edge_count.items():
    ra, rb = doi2rec.get(a), doi2rec.get(b)
    if ra and rb:
        G.add_edge(ra, rb, weight=c)
iso = [r["id"] for r in records if r["id"] not in G]
print("graph:", G.number_of_nodes(), "nodes,", G.number_of_edges(), "edges,", len(iso), "isolated")

try:
    from networkx.algorithms.community import louvain_communities as community_fn
    comms = community_fn(G, weight="weight", seed=42)
    algo = "louvain"
except Exception:
    try:
        from networkx.algorithms.community import asyn_lpa_communities as community_fn
        comms = list(community_fn(G, weight="weight", seed=42))
        algo = "asyn_lpa"
    except Exception:
        from networkx.algorithms.community import greedy_modularity_communities as community_fn
        comms = list(community_fn(G, weight="weight"))
        algo = "greedy_modularity"
print("clusters:", len(comms), "via", algo)

cl_of = {}
for ci, com in enumerate(sorted(comms, key=len, reverse=True)):
    for n in com:
        cl_of[n] = ci

pos = nx.spring_layout(G, weight="weight", seed=42, iterations=120)
if iso:
    import math
    xs = [p[0] for p in pos.values()] or [0]
    ys = [p[1] for p in pos.values()] or [0]
    cx, cy = sum(xs) / len(xs), sum(ys) / len(ys)
    rad = max(max(abs(x - cx) for x in xs), max(abs(y - cy) for y in ys)) * 1.35 + 1.0
    for k, rid in enumerate(iso):
        ang = 2 * math.pi * k / len(iso)
        pos[rid] = (cx + rad * math.cos(ang), cy + rad * math.sin(ang))

layout = {}
for r in records:
    rid = r["id"]
    d = norm_doi(r["doi"])
    x, y = pos.get(rid, (0.0, 0.0))
    layout[rid] = {
        "x": round(float(x), 5), "y": round(float(y), 5),
        "cl": cl_of.get(rid, -1),
        "c": oa[d]["cited_by_count"] if d and d in oa else 0,
    }
json.dump(layout, open(os.path.join(EXPORT, "sadb_layout.json"), "w"))
print("wrote export/sadb_layout.json")

# --- static preview PNG ---
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(16, 12))
    cmap = plt.get_cmap("tab20")
    for rid, p in pos.items():
        ci = cl_of.get(rid, -1)
        ax.scatter(p[0], p[1], s=4 + 0.6 * layout[rid]["c"] ** 0.5,
                   c=[cmap(ci % 20) if ci >= 0 else "#bbbbbb"])
    ax.set_title("SADb corpus — %d papers, %d citation edges (%s clusters)" %
                 (len(records), G.number_of_edges(), algo))
    ax.set_xticks([]), ax.set_yticks([])
    fig.savefig(os.path.join(VOS, "sadb_preview.png"), dpi=110, bbox_inches="tight")
    print("wrote sadb_preview.png")
except Exception as e:
    print("preview skipped:", e)
