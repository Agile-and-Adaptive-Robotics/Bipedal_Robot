"""PDF hunt pass 2 (2026-09-20): find OA PDF URLs for remaining_no_pdf.csv records.

Verify-before-attach (the Sept-11/15 lessons): every URL is range-GET'd and must
start with the %PDF magic BEFORE it is attached. Attaches via URL only (the PAT
cannot upload binaries). Results -> author_fix/pdf_hunt_pass2_20260920.csv.
Run: myo python pdf_hunt2.py
"""
import csv, json, os, re, time, urllib.parse, urllib.request

SAD = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(SAD, "author_fix", "pdf_hunt_pass2_20260920.csv")
MAILTO = "benjamin.bolen@pdx.edu"
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 SADb/1.0 (mailto:%s)" % MAILTO}

txt = open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()
PAT = re.search(r"PAT\s*:\s*(pat[^\s]+)", txt).group(1)
BASE = "appMQTnobUNRytIp7"


def air(method, path, body=None):
    req = urllib.request.Request("https://api.airtable.com/v0/" + path,
                                 data=json.dumps(body).encode() if body is not None else None,
                                 headers={"Authorization": "Bearer " + PAT, "Content-Type": "application/json"},
                                 method=method)
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def norm(d):
    d = (d or "").strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d


# hunt list: ALL DOI-bearing records without a PDF (live export), minus already-hunted
verdicted = set()
try:
    with open(OUT, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            verdicted.add(norm(row.get("doi")))
except FileNotFoundError:
    pass
export = json.load(open(os.path.join(SAD, "export", "sadb_export.json"), encoding="utf-8"))
doi2rec = {}
todo = []
for r in export:
    d = norm(r["doi"])
    if d:
        doi2rec[d] = r["id"]
        if not r["has_pdf"] and d not in verdicted:
            todo.append(d)
print(f"corpus-no-pdf-with-doi to-hunt now: {len(todo)} (already verdicted: {len(verdicted)})")

# OpenAlex candidates (batched)
cand = {}
B = 50
for i in range(0, len(todo), B):
    chunk = todo[i:i + B]
    filt = "doi:" + "|".join(chunk)
    url = ("https://api.openalex.org/works?filter=" + urllib.parse.quote(filt, safe="|")
           + "&select=doi,best_oa_location,locations&per-page=50&mailto=" + MAILTO)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "sadb-pdfhunt/1.0 (mailto:%s)" % MAILTO})
        data = json.loads(urllib.request.urlopen(req, timeout=60).read().decode())
    except Exception as e:
        print("OA batch fail", i // B, e)
        time.sleep(3)
        continue
    for w in data.get("results", []):
        d = norm(w.get("doi") or "")
        urls = []
        locs = [w.get("best_oa_location")] + (w.get("locations") or [])
        for loc in locs:
            u = (loc or {}).get("pdf_url")
            if u and u not in urls:
                urls.append(u)
        if urls:
            cand[d] = urls
    print(f"OA batch {i//B}: {len(cand)} with candidates so far")
    time.sleep(0.5)
print("records with OA candidates:", len(cand))


def verify(url):
    """True if a range GET returns bytes starting with %PDF."""
    try:
        req = urllib.request.Request(url, headers={**UA, "Range": "bytes=0-63"})
        with urllib.request.urlopen(req, timeout=30) as r:
            head = r.read(64)
        return head.startswith(b"%PDF")
    except Exception:
        return False


results = []
n_att = 0
for i, d in enumerate(todo):
    rid = doi2rec[d]
    urls = cand.get(d, [])
    status, good = "no-oa-candidate", ""
    if urls:
        status, good = "no-verified-pdf", ""
        for u in urls:
            if verify(u):
                good = u
                status = "verified"
                break
        if good:
            try:
                rec = air("GET", f"{BASE}/Papers/{rid}?fields%5B%5D=Attachments")
                existing = [{"id": a["id"]} for a in rec.get("fields", {}).get("Attachments", [])]
                air("PATCH", f"{BASE}/Papers/{rid}", {"fields": {"Attachments": existing + [{"url": good}]}})
                status = "ATTACHED"
                n_att += 1
            except Exception as e:
                status = f"attach-failed: {str(e)[:80]}"
                good = ""
    results.append([d, rid, len(urls), status, good])
    if i % 25 == 0:
        print(f"{i}/{len(todo)} checked, {n_att} attached so far")
    time.sleep(0.4)

with open(OUT, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["doi", "record_id", "n_candidates", "status", "verified_url"])
    w.writerows(results)
print(f"DONE: {len(results)} records checked, {n_att} PDFs attached -> {OUT}")
from collections import Counter
print(Counter(r[3].split(':')[0] for r in results))
