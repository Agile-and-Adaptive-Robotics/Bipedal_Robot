"""Build VOSviewer map + network files for the SADb corpus (2026-09-14, EB475WS4).

Sources (552 papers total):
  - SADb_audit/airtable_papers_slim.csv        (99 originals: id/title/author/year/animals, DOI)
  - SADb_audit/airtable_created50_ids.csv      (50 demo-50: zotero_key -> airtable_id)
  - SADb_audit/personal_SAD_items_slim.csv     (zotero metadata for the 50; keyed by key)
  - SADb_audit/airtable_rest_import_clean.csv  (383 rest-import: zotero_key/title/author/year/DOI)
  - digest20 inline below                      (created in Airtable 2026-09-14)

Outputs in SADb_audit/vosviewer/:
  sadb_map.txt               VOSviewer map (id, label, doi, year, first_author, weight, source)
  sadb_network.txt           citation edges within the corpus (source, target, weight)
  openalex_enrichment.csv    per-DOI OpenAlex metadata (oa_id, cited_by_count, n_refs, in_corpus_refs)
  openalex_misses.txt        DOIs OpenAlex could not resolve
  sadb_preview.png           quick networkx spring preview (not a VOSviewer product)

OpenAlex: batched doi| filters (50/call), polite pool.
"""
import csv, json, os, sys, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)
OUT = HERE
MAILTO = "benjamin.bolen@pdx.edu"

DIGEST20 = [
    ("10.7554/elife.98841", "Operation regimes of spinal circuits controlling locomotion and the role of supraspinal drives and sensory feedback", "Rybak et al.", 2024),
    ("10.7554/elife.103504", "Operation of spinal sensorimotor circuits controlling phase durations during tied-belt and split-belt locomotion after a lateral thoracic hemisection", "Rybak et al.", 2025),
    ("10.1152/jn.00104.2024", "Forelimb movements contribute to hindlimb cutaneous reflexes during locomotion in cats", "Harnie et al.", 2024),
    ("10.1113/jp286151", "Changes in intra- and interlimb reflexes from hindlimb cutaneous afferents after staggered thoracic lateral hemisections during locomotion in cats", "Mari et al.", 2024),
    ("10.1113/jp286808", "Changes in intra- and interlimb reflexes from forelimb cutaneous afferents after staggered thoracic lateral hemisections during locomotion in cats", "Mari et al.", 2024),
    ("10.1016/j.neunet.2024.106422", "A spinal circuit model with asymmetric cervical-lumbar layout controls backward locomotion and scratching in quadrupeds", "Zhu et al.", 2024),
    ("10.1098/rsos.240207", "Sensory feedback and central neuronal interactions in mouse locomotion", "Molkov et al.", 2024),
    ("10.1152/jn.00248.2023", "Dynamic spinal reflex adaptation during locomotor adaptation", "Refy et al.", 2023),
    ("10.3389/fncir.2023.1235181", "Lumbar V3 interneurons provide direct excitatory synaptic input onto thoracic sympathetic preganglionic neurons, linking locomotor, and autonomic spinal systems", "Chacon et al.", 2023),
    ("10.1016/j.cub.2023.07.014", "Distinct roles of spinal commissural interneurons in transmission of contralateral sensory information", "Laflamme et al.", 2023),
    ("10.1371/journal.pcbi.1012101", "Balancing central control and sensory feedback produces adaptable and robust locomotor patterns in a spiking, neuromechanical model of the salamander spinal cord", "Pazzaglia et al.", 2025),
    ("10.1371/journal.pcbi.1013494", "A physiologically inspired hybrid CPG/Reflex controller for cycling simulations that generalizes to walking", "Severini et al.", 2025),
    ("10.1016/j.cub.2025.09.030", "A spinal circuit for skilled locomotion", "Toscano et al.", 2025),
    ("10.1016/j.expneurol.2023.114496", "Spinal control of locomotion before and after spinal cord injury", "Danner et al.", 2023),
    ("10.1523/jneurosci.2015-22.2023", "Excitatory and Inhibitory Descending Commissural Interneurons Differentially Integrate Supraspinal and Segmental Sensory Signals", "Giorgi et al.", 2023),
    ("10.1007/s00422-023-00970-z", "The Bcm rule allows a spinal cord model to learn rhythmic movements", "Kohler et al.", 2023),
    ("10.1038/s42003-024-06843-w", "Interlimb coordination is not strictly controlled during walking", "Arai et al.", 2024),
    ("10.1152/jn.00331.2025", "Speed-dependent locomotor adjustments following staggered thoracic lateral hemisections in adult cats", "Yassine et al.", 2025),
    ("10.1101/2025.11.11.687930", "Mechanisms of adaptive interlimb coordination to sudden ground loss: a neuromusculoskeletal modeling study", "Shinohara et al.", 2025),
    ("10.7554/elife.107480", "Linking spinal circuit reorganization to recovery after thoracic spinal cord injury", "Shevtsova et al.", 2025),
]


def norm_doi(d):
    if not d:
        return ""
    d = d.strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d


def clean_ws(s):
    return " ".join((s or "").split())


def first_surname(author_field):
    a = clean_ws(author_field or "")
    if not a:
        return ""
    return a.split(" et al")[0].split(" and ")[0].split(",")[0].split(";")[0].strip()


nodes = []  # dicts: doi, label, author, year, source

# 1) originals 99
with open(os.path.join(SAD, "airtable_papers_slim.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        nodes.append({
            "doi": norm_doi(r.get("DOI") or r.get("doi") or ""),
            "label": clean_ws(r.get("title") or r.get("Title") or ""),
            "author": clean_ws(r.get("author") or r.get("Author") or ""),
            "year": (r.get("year") or r.get("Year") or "").strip(),
            "source": "originals99",
        })

# 2) demo-50 via zotero metadata
keys50 = {}
with open(os.path.join(SAD, "airtable_created50_ids.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        k = (r.get("zotero_key") or "").strip()
        if k:
            keys50[k] = True
if keys50:
    with open(os.path.join(SAD, "personal_SAD_items_slim.csv"), encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            k = (r.get("key") or r.get("zotero_key") or "").strip()
            if k in keys50:
                nodes.append({
                    "doi": norm_doi(r.get("DOI") or r.get("doi") or ""),
                    "label": clean_ws(r.get("title") or ""),
                    "author": clean_ws(r.get("author") or r.get("creators") or ""),
                    "year": (r.get("year") or "").strip(),
                    "source": "demo50",
                })

# 3) rest-import 383
with open(os.path.join(SAD, "airtable_rest_import_clean.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        nodes.append({
            "doi": norm_doi(r.get("DOI") or ""),
            "label": clean_ws(r.get("title") or ""),
            "author": clean_ws(r.get("author") or ""),
            "year": (r.get("year") or "").strip(),
            "source": "rest383",
        })

# 4) digest 20
for doi, title, author, year in DIGEST20:
    nodes.append({"doi": norm_doi(doi), "label": clean_ws(title), "author": author,
                  "year": str(year), "source": "digest20"})

# dedupe by doi (keep first), keep no-DOI rows keyed by label
seen = set()
uniq = []
for n in nodes:
    key = n["doi"] if n["doi"] else "nodoi::" + n["label"].lower()
    if key in seen:
        continue
    seen.add(key)
    uniq.append(n)
nodes = uniq
print("nodes:", len(nodes), {s: sum(1 for n in nodes if n["source"] == s) for s in set(n["source"] for n in nodes)})

# --- OpenAlex enrichment (batched) ---
dois = [n["doi"] for n in nodes if n["doi"]]
oa = {}          # doi -> dict(oa_id, cited_by_count, refs)
misses = []
B = 50
for i in range(0, len(dois), B):
    chunk = dois[i:i + B]
    filt = "doi:" + "|".join(chunk)
    url = ("https://api.openalex.org/works?filter=" + urllib.parse.quote(filt, safe="|")
           + "&select=id,doi,cited_by_count,referenced_works&per-page=50&mailto=" + MAILTO)
    got = 0
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "sadb-vosviewer/0.1 (mailto:%s)" % MAILTO}), timeout=60) as resp:
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
            oa[d] = {"oa_id": w["id"], "cited_by_count": w.get("cited_by_count", 0),
                     "refs": w.get("referenced_works", [])}
            got += 1
    print("openalex batch %d: %d/%d resolved" % (i // B, got, len(chunk)))
    time.sleep(0.6)

for n in nodes:
    if n["doi"] and n["doi"] not in oa:
        misses.append(n["doi"])
print("openalex resolved:", len(oa), "misses:", len(misses))

# oa_id -> doi for corpus-edge mapping
oa2doi = {v["oa_id"]: k for k, v in oa.items()}

# edges: citation pairs within corpus
edge_count = {}
for d, v in oa.items():
    for r in v["refs"]:
        if r in oa2doi and oa2doi[r] != d:
            a, b = sorted((d, oa2doi[r]))
            edge_count[(a, b)] = edge_count.get((a, b), 0) + 1
print("internal citation edges:", len(edge_count))

# --- VOSviewer map file ---
os.makedirs(OUT, exist_ok=True)
map_path = os.path.join(OUT, "sadb_map.txt")
with open(map_path, "w", encoding="utf-8", newline="\n") as f:
    f.write("\t".join(["id", "label", "doi", "year", "first_author", "weight", "source"]) + "\n")
    for n in nodes:
        wid = ""
        if n["doi"] and n["doi"] in oa:
            wid = oa[n["doi"]]["oa_id"].rsplit("/", 1)[-1]
        nid = wid if wid else ("doi:" + n["doi"] if n["doi"] else "nodoi:" + n["label"][:60])
        f.write("\t".join([
            nid, n["label"][:120], n["doi"], n["year"], first_surname(n["author"]),
            str(oa[n["doi"]]["cited_by_count"]) if n["doi"] and n["doi"] in oa else "0",
            n["source"],
        ]) + "\n")
print("wrote", map_path)

net_path = os.path.join(OUT, "sadb_network.txt")
with open(net_path, "w", encoding="utf-8", newline="\n") as f:
    f.write("source\ttarget\tweight\n")
    for (a, b), c in sorted(edge_count.items(), key=lambda kv: -kv[1]):
        wa = oa[a]["oa_id"].rsplit("/", 1)[-1]
        wb = oa[b]["oa_id"].rsplit("/", 1)[-1]
        f.write("%s\t%s\t%d\n" % (wa, wb, c))
print("wrote", net_path)

enr_path = os.path.join(OUT, "openalex_enrichment.csv")
with open(enr_path, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["doi", "oa_id", "cited_by_count", "n_refs", "n_in_corpus_refs"])
    for d, v in sorted(oa.items(), key=lambda kv: -kv[1]["cited_by_count"]):
        w.writerow([d, v["oa_id"], v["cited_by_count"], len(v["refs"]),
                    sum(1 for r in v["refs"] if r in oa2doi)])
print("wrote", enr_path)

with open(os.path.join(OUT, "openalex_misses.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(misses))
print("wrote misses")
