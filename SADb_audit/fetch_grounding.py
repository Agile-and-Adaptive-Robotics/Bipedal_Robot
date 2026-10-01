"""Fetch grounding text (abstracts) for every un-curated Papers record (2026-09-28).

Sources, in order: OpenAlex (abstract_inverted_index, batched by DOI), then
Europe PMC (abstractText, single-DOI queries, cached). Writes:
  export/grounding.json      {recordId: {doi,title,year,abstract,source,oa_cited,topic}}
  curation_queue/queue.json  ordered list of curatable papers (id,title,primary,
                             year,doi,abstract,source) — the subagent work queue
  curation_queue/notext.json bare papers with no findable text (log as no-text)

Stdlib only. Run: myo python fetch_grounding.py
"""
import json, os, re, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
EXPORT = os.path.join(HERE, "export")
QUEUE = os.path.join(HERE, "curation_queue")
os.makedirs(QUEUE, exist_ok=True)
MAILTO = "benjamin.bolen@pdx.edu"


def norm_doi(d):
    d = (d or "").strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d


def get_json(url, headers=None, timeout=60):
    req = urllib.request.Request(url, headers=headers or {"User-Agent": "sadb-grounding/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


records = json.load(open(os.path.join(EXPORT, "sadb_export.json"), encoding="utf-8"))
bare = [r for r in records if not r["has_notes"]]
print(f"bare records: {len(bare)}; with DOI: {sum(1 for r in bare if norm_doi(r['doi']))}")

ground_path = os.path.join(EXPORT, "grounding.json")
ground = {}
if os.path.exists(ground_path):
    ground = json.load(open(ground_path, encoding="utf-8"))
epmc_cache_path = os.path.join(EXPORT, "epmc_cache.json")
epmc_cache = {}
if os.path.exists(epmc_cache_path):
    epmc_cache = json.load(open(epmc_cache_path, encoding="utf-8"))

doi2rec = {}
for r in bare:
    d = norm_doi(r["doi"])
    if d:
        doi2rec.setdefault(d, r["id"])

todo = [d for d in sorted(doi2rec) if not (ground.get(doi2rec[d], {}).get("abstract") or "").strip()]
print(f"DOIs needing abstract fetch: {len(todo)}")

# --- pass 1: OpenAlex abstracts, 50 DOIs per request ---
B = 50
for i in range(0, len(todo), B):
    chunk = todo[i:i + B]
    filt = "doi:" + "|".join(chunk)
    url = ("https://api.openalex.org/works?filter=" + urllib.parse.quote(filt, safe="|")
           + "&select=id,doi,title,publication_year,abstract_inverted_index,cited_by_count,primary_topic"
           + "&per-page=50&mailto=" + MAILTO)
    data = None
    for attempt in range(3):
        try:
            data = get_json(url)
            break
        except Exception as e:
            print("  oa retry", i // B, e)
            time.sleep(3 + 3 * attempt)
    if data is None:
        print("  OA FAILED batch", i // B)
        continue
    for w in data.get("results", []):
        d = norm_doi(w.get("doi") or "")
        rid = doi2rec.get(d)
        if not rid:
            continue
        inv = w.get("abstract_inverted_index")
        abstract = ""
        if inv:
            pos = {}
            for word, idxs in inv.items():
                for ix in idxs:
                    pos[ix] = word
            abstract = " ".join(pos[k] for k in sorted(pos))
        rec = {"doi": d, "title": w.get("title") or "", "year": w.get("publication_year"),
               "abstract": abstract.strip(), "source": "openalex",
               "oa_cited": w.get("cited_by_count", 0),
               "topic": ((w.get("primary_topic") or {}).get("display_name") or "")}
        old = ground.get(rid, {})
        if not abstract.strip():
            rec["abstract"] = old.get("abstract", "")
            rec["source"] = old.get("source", "openalex")
        else:
            rec["source"] = "openalex"
        ground[rid] = rec
    if i % (B * 4) == 0:
        print(f"  openalex batch {i//B}/{(len(todo)+B-1)//B}")
    time.sleep(0.55)

# --- pass 2: Europe PMC fallback for still-empty abstracts ---
need_epmc = [d for d in todo
             if not (ground.get(doi2rec.get(d, ""), {}).get("abstract") or "").strip()
             and d not in epmc_cache]
print(f"Europe PMC fallback: {len(need_epmc)} DOIs")
for k, d in enumerate(need_epmc):
    rid = doi2rec[d]
    try:
        q = urllib.parse.quote(f'DOI:"{d}"')
        j = get_json(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}"
                     f"&resultType=core&format=json&pageSize=1")
        res = (j.get("resultList") or {}).get("result") or []
        ab = (res[0].get("abstractText") or "") if res else ""
        epmc_cache[d] = ab
        if ab.strip():
            rec = ground.setdefault(rid, {"doi": d})
            rec.update({"abstract": re.sub(r"<[^>]+>", " ", ab).strip(), "source": "epmc"})
    except Exception as e:
        epmc_cache[d] = ""
        if (k % 25) == 0:
            print("  epmc", k, e)
    time.sleep(0.25)
    if (k % 40) == 39:
        json.dump(epmc_cache, open(epmc_cache_path, "w"))
json.dump(epmc_cache, open(epmc_cache_path, "w"))
for d, ab in epmc_cache.items():
    rid = doi2rec.get(d)
    if rid and ab.strip() and not (ground.get(rid, {}).get("abstract") or "").strip():
        rec = ground.setdefault(rid, {"doi": d})
        rec.update({"abstract": re.sub(r"<[^>]+>", " ", ab).strip(), "source": "epmc"})

# also keep no-DOI bare papers listed (never groundable via these routes)
json.dump(ground, open(ground_path, "w"), ensure_ascii=False)

by_rec = {r["id"]: r for r in records}
queue, notext = [], []
for r in bare:
    g = ground.get(r["id"], {})
    entry = {"id": r["id"], "title": r["title"], "primary": r["primary"],
             "year": r["year"], "doi": norm_doi(r["doi"]),
             "abstract": g.get("abstract", ""), "source": g.get("source", ""),
             "topic": g.get("topic", ""), "oa_cited": g.get("oa_cited", 0)}
    if entry["abstract"].strip():
        queue.append(entry)
    else:
        entry["reason"] = "no-doi" if not entry["doi"] else "no-abstract-found"
        notext.append(entry)

json.dump(queue, open(os.path.join(QUEUE, "queue.json"), "w"), ensure_ascii=False, indent=1)
json.dump(notext, open(os.path.join(QUEUE, "notext.json"), "w"), ensure_ascii=False, indent=1)
print(f"queue: {len(queue)} curatable | notext: {len(notext)} "
      f"({sum(1 for x in notext if x['reason']=='no-doi')} no-DOI)")
