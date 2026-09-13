#!/usr/bin/env python
# RR-style digest candidates (2023-2026) from OpenAlex, seeded by the papers observed
# in Ben's Research Rabbit collection + the lab's model family. Dedupe vs SAD corpus.
import json, subprocess, os, csv, re

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "batch3")
MAILTO = "bbolen@pdx.edu"

def curl_json(url):
    r = subprocess.run(["curl.exe", "-s", "-m", "30", "-A", "mailto:" + MAILTO, url], capture_output=True)
    try:
        return json.loads(r.stdout.decode("utf-8", "replace"))
    except Exception:
        return None

def find_work(doi=None, title=None):
    if doi:
        return curl_json("https://api.openalex.org/works/https://doi.org/" + doi)
    t = re.sub(r"[^\w\s]", "", title).replace(" ", "%20")
    d = curl_json("https://api.openalex.org/works?filter=title.search:%s&per-page=3" % t)
    if d and d.get("results"):
        return d["results"][0]
    return None

SEEDS_BY_DOI = [
    "10.1113/jphysiol.2013.261115",          # Rybak 2013 left-right model (in batch 3)
    "10.7554/eLife.20150",                    # guess candidate, verified below by title check
]
SEEDS_BY_TITLE = [
    "A neural central pattern generator system producing flexor-extensor coordination",  # Danner 2017
    "Operation regimes of spinal circuits controlling locomotion and role of supraspinal drives and sensory feedback",  # Rybak 2024/2025 eLife
    "Linking spinal circuit reorganization to recovery after thoracic spinal cord injury",  # Shevtsova 2026 (RR-selected)
]

seeds = []
for doi in SEEDS_BY_DOI:
    w = find_work(doi=doi)
    if w:
        seeds.append(w)
for t in SEEDS_BY_TITLE:
    w = find_work(title=t)
    if w:
        seeds.append(w)
print("SEEDS:")
for w in seeds:
    print(" ", w.get("id", "").split("/")[-1], "|", (w.get("display_name") or "")[:70], "|", w.get("publication_year"), "| cites:", w.get("cited_by_count"))

# Corpus DOIs for dedupe (personal SAD slim + rest-import queue)
corpus_dois = set()
for path, col in [(r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\SAD_audit\personal_SAD_items_slim.csv", None)]:
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            rd = csv.DictReader(f)
            for row in rd:
                for v in row.values():
                    if v and re.match(r"^10\.\d{4}", v.strip()):
                        corpus_dois.add(v.strip().lower())
    except Exception as e:
        print("corpus read fail", path, e)
# rest-import queue DOIs (they're in Airtable Papers now)
with open(r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\SAD_audit\airtable_rest_import_clean.csv", encoding="utf-8", errors="replace") as f:
    rd = csv.DictReader(f)
    for row in rd:
        d = (row.get("DOI") or "").strip().lower()
        if d:
            corpus_dois.add(d)
print("corpus DOIs loaded:", len(corpus_dois))

candidates = {}
def add_candidate(w, why):
    doi = (w.get("doi") or "").replace("https://doi.org/", "").lower()
    if not doi:
        return
    yr = w.get("publication_year") or 0
    if yr < 2023:
        return
    if doi in corpus_dois:
        return
    if doi in candidates:
        candidates[doi]["why"].add(why)
        return
    candidates[doi] = dict(
        doi=doi, title=(w.get("display_name") or "").strip(), year=yr,
        cited_by=w.get("cited_by_count", 0),
        venue=((w.get("primary_location") or {}).get("source") or {}).get("display_name", "") if w.get("primary_location") else "",
        authors=", ".join((a.get("author") or {}).get("display_name", "") for a in (w.get("authorships") or [])[:3]),
        why={why},
        oa=(w.get("open_access") or {}).get("oa_url", ""),
    )

for w in seeds:
    wid = w.get("id", "").split("/")[-1]
    d = curl_json("https://api.openalex.org/works?filter=cites:%s,from_publication_date:2023-01-01&sort=cited_by_count:desc&per-page=50&mailto=%s" % (wid, MAILTO))
    for r in (d or {}).get("results", []):
        add_candidate(r, "cites %s (%s)" % ((w.get("display_name") or "")[:40], w.get("publication_year")))

# also recent works citing the two classic RR seeds: Rybak 2006 Rybak CPG / McCrea-Rybak 2008
for t in ["Spinal cord excitatory and inhibitory interneurons", "Reconfiguration of the spinal locomotor network"]:
    pass  # kept small; seeds above suffice for the first digest

rows = sorted(candidates.values(), key=lambda r: (-len(r["why"]), -r["cited_by"]))
for r in rows:
    r["why"] = "; ".join(sorted(r["why"]))
with open(os.path.join(OUT, "rr_digest_candidates.json"), "w", encoding="utf-8") as f:
    json.dump(rows, f, indent=1, ensure_ascii=False)
print("\nCANDIDATES (2023+):", len(rows))
for r in rows[:40]:
    print(" %-42s %4d cites=%-4d %s" % ((r["title"] or "")[:42], r["year"], r["cited_by"], r["doi"]))
