#!/usr/bin/env python
# RR digest round 2: better seeds = the actual RR recommendations + Danner 2017 (found properly).
import json, subprocess, os, csv, re

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "batch3")
MAILTO = "bbolen@pdx.edu"

def curl_json(url):
    r = subprocess.run(["curl.exe", "-s", "-m", "30", "-A", "mailto:" + MAILTO, url], capture_output=True)
    try:
        return json.loads(r.stdout.decode("utf-8", "replace"))
    except Exception:
        return None

TITLE_QUERIES = [
    "Simultaneous control of forward and backward locomotion by spinal sensorimotor circuits",   # Audet 2023
    "A sensory signal related to left-right symmetry modulates intra- and interlimb cutaneous reflexes during locomotion in intact cats",  # Mari 2023
    "Locomotor speed control circuits in the caudal brainstem",                                   # Nature tab (seed family)
    "The role of V3 neurons in speed-dependent interlimb coordination during locomotion in mice", # Zhang (published)
    "A neural central pattern generator producing flexor-extensor coordination",                  # Danner 2017 variants
    "Neuromusculoskeletal model for the investigation of adaptive locomotion" ,                   # Danner 2017 alternate title guess
]
seeds = []
for t in TITLE_QUERIES:
    q = re.sub(r"[^\w\s]", " ", t).replace(" ", "%20")
    d = curl_json("https://api.openalex.org/works?filter=title.search:%s&per-page=3&mailto=%s" % (q, MAILTO))
    for r in (d or {}).get("results", [])[:1]:
        seeds.append(r)

seeds_by_doi = [
    "10.1113/jphysiol.2013.261115",
    "10.7554/elife.107480",
]
for doi in seeds_by_doi:
    w = curl_json("https://api.openalex.org/works/https://doi.org/" + doi + "?mailto=" + MAILTO)
    if w and w.get("id"):
        seeds.append(w)

print("SEEDS:")
seen_ids = set()
for w in seeds:
    wid = w.get("id", "").split("/")[-1]
    if wid in seen_ids:
        continue
    seen_ids.add(wid)
    print(" ", wid, "|", (w.get("display_name") or "")[:70], "|", w.get("publication_year"), "| cites:", w.get("cited_by_count"))

corpus_dois = set()
with open(r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\SAD_audit\personal_SAD_items_slim.csv", encoding="utf-8", errors="replace") as f:
    for row in csv.DictReader(f):
        for v in row.values():
            if v and re.match(r"^10\.\d{4}", v.strip()):
                corpus_dois.add(v.strip().lower())
with open(r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\SAD_audit\airtable_rest_import_clean.csv", encoding="utf-8", errors="replace") as f:
    for row in csv.DictReader(f):
        d = (row.get("DOI") or "").strip().lower()
        if d:
            corpus_dois.add(d)

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
    src = ((w.get("primary_location") or {}).get("source") or {})
    candidates[doi] = dict(
        doi=doi, title=(w.get("display_name") or "").strip(), year=yr,
        cited_by=w.get("cited_by_count", 0), venue=src.get("display_name", ""),
        authors=", ".join((a.get("author") or {}).get("display_name", "") for a in (w.get("authorships") or [])[:3]),
        why={why}, oa=(w.get("open_access") or {}).get("oa_url", ""),
        oa_status=(w.get("open_access") or {}).get("oa_status", ""),
    )

for w in seeds:
    wid = w.get("id", "").split("/")[-1]
    label = (w.get("display_name") or "")[:38]
    d = curl_json("https://api.openalex.org/works?filter=cites:%s,from_publication_date:2023-01-01&sort=cited_by_count:desc&per-page=50&mailto=%s" % (wid, MAILTO))
    for r in (d or {}).get("results", []):
        add_candidate(r, "cites [%s %s]" % (label, w.get("publication_year")))
    # shared-reference route: referenced_works of 2023+ papers are heavy; instead pull
    # this seed's references that are themselves 2023+ (recent work the seed builds on)
    for refid in (w.get("referenced_works") or [])[:200]:
        pass  # per-ref fetch too chatty; skip round 2

rows = sorted(candidates.values(), key=lambda r: (-len(r["why"]), -r["cited_by"]))
for r in rows:
    r["why"] = "; ".join(sorted(r["why"]))
with open(os.path.join(OUT, "rr_digest_candidates.json"), "w", encoding="utf-8") as f:
    json.dump(rows, f, indent=1, ensure_ascii=False)
print("\nCANDIDATES (2023+):", len(rows))
for r in rows:
    print(" %-58s %4d cites=%-4d %s" % ((r["title"] or "")[:58], r["year"], r["cited_by"], r["doi"]))
