"""Build the citation graph data for the SADb HTML app (2026-09-22).

(VOSviewer export retired 2026-09-22 — generic viewer, no drill-through; the
vosviewer/ folder and its map/network files were deleted. The app is the viz.)

Reads export/sadb_export.json, resolves every DOI against OpenAlex (cached in
export/openalex_raw.json — only missing DOIs are fetched), and writes:
  export/sadb_cites.json   DIRECTED citation adjacency {citerId: [citedIds]}
                           (drives the app's References / Cited-by views and
                           the excitatory/inhibitory synapse links)
  export/sadb_layout.json  {recordId: {x, y, cl, c}} — spring layout of the
                           in-corpus citation network + Louvain cluster +
                           citation count (drives the Topic-landscape map mode)

Run:  myo python build_citation_graph.py      (after export_corpus.py)
"""
import csv, json, math, os, time, urllib.parse, urllib.request
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
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


# --- OpenAlex enrichment (cached) ---
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
                    url, headers={"User-Agent": "sadb-cites/1.0 (mailto:%s)" % MAILTO}),
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

resolved = sum(1 for r in records if norm_doi(r["doi"]) in oa)
misses = [norm_doi(r["doi"]) for r in records if norm_doi(r["doi"]) and norm_doi(r["doi"]) not in oa]
print("openalex resolved:", resolved, "misses:", len(misses))
oa2doi = {v["oa_id"]: k for k, v in oa.items()}

doi2rec = {}
for r in records:
    d = norm_doi(r["doi"])
    if d:
        doi2rec.setdefault(d, r["id"])

edge_count = {}
for d, v in oa.items():
    for ref in v["refs"]:
        if ref in oa2doi and oa2doi[ref] != d:
            a, b = sorted((d, oa2doi[ref]))
            edge_count[(a, b)] = edge_count.get((a, b), 0) + 1
print("internal citation edges:", len(edge_count))

# directed citation adjacency by Airtable record id (referenced_works are full
# OpenAlex URLs, so the lookup goes oa_id -> record id)
oaid2rec = {oa[d]["oa_id"]: rid for d, rid in doi2rec.items() if d in oa}
cites = {}
for d, v in oa.items():
    ra = doi2rec.get(d)
    if not ra:
        continue
    for ref in v["refs"]:
        rb = oaid2rec.get(ref)
        if rb and rb != ra:
            lst = cites.setdefault(ra, [])
            if rb not in lst:
                lst.append(rb)
json.dump(cites, open(os.path.join(EXPORT, "sadb_cites.json"), "w"))
print("directed corpus citation pairs:", sum(len(v) for v in cites.values()))

# --- layout + clusters for the app's Topic-landscape mode ---
import networkx as nx

G = nx.Graph()
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
