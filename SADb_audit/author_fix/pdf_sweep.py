"""PDF sweep: for all 552 corpus DOIs, find anonymously-fetchable OA PDF URLs
(OpenAlex best_oa_location), verify with a range GET (%PDF magic), and emit:
  - pdfs_to_attach.json  chunked Airtable update payloads (upsert on DOI,
    Attachments = [{url}]) for verified-fetchable PDFs
  - pdfs_paywalled.csv   DOIs with no fetchable PDF (need PSU-browser upload)
Run with any python (stdlib only).
"""
import csv, json, os, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)
MAILTO = "benjamin.bolen@pdx.edu"

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

# corpus DOIs from all sources
dois = set()
with open(os.path.join(SAD, "airtable_papers_slim.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r.get("DOI"):
            dois.add((r.get("DOI") or "").strip())
with open(os.path.join(SAD, "airtable_rest_import_clean.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r.get("DOI"):
            dois.add(r["DOI"].strip())
with open(os.path.join(SAD, "batch3", "rr_digest_shortlist.json"), encoding="utf-8") as f:
    for s in json.load(f):
        dois.add(s["doi"])
dois.add("10.7554/elife.107480")
dois = sorted({norm_doi(d) for d in dois if d})
print("corpus DOIs:", len(dois))

# OpenAlex: best_oa_location.pdf_url
oa = {}
B = 50
for i in range(0, len(dois), B):
    chunk = dois[i:i + B]
    filt = "doi:" + "|".join(chunk)
    url = ("https://api.openalex.org/works?filter=" + urllib.parse.quote(filt, safe="|")
           + "&select=doi,best_oa_location&per-page=50&mailto=" + MAILTO)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "sadb-pdf/0.1 (mailto:%s)" % MAILTO}), timeout=60) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            break
        except Exception as e:
            print("  retry", i // B, e)
            time.sleep(3)
    else:
        continue
    for w in data.get("results", []):
        d = norm_doi(w.get("doi") or "")
        loc = (w.get("best_oa_location") or {})
        pdf = loc.get("pdf_url") or ""
        if d and pdf:
            oa[d] = pdf
    print("openalex batch %d done (cumulative pdf urls: %d)" % (i // B, len(oa)))
    time.sleep(0.5)

print("OA pdf urls found:", len(oa))

# verify: anonymous GET, first bytes = %PDF
def fetchable(pdf_url):
    try:
        req = urllib.request.Request(pdf_url, headers={
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)",
            "Range": "bytes=0-1023"})
        with urllib.request.urlopen(req, timeout=25) as r:
            head = r.read(1024)
            return head.startswith(b"%PDF")
    except Exception:
        return False

ok = {}
for d, u in sorted(oa.items()):
    if fetchable(u):
        ok[d] = u
    time.sleep(0.2)
print("verified fetchable:", len(ok), "| blocked/paywalled:", len(oa) - len(ok))

os.makedirs(HERE, exist_ok=True)
with open(os.path.join(HERE, "pdfs_to_attach.json"), "w", encoding="utf-8") as f:
    recs = [{"fields": {"DOI": d, "Attachments": [{"url": u}]}} for d, u in ok.items()]
    chunks = [recs[i:i + 50] for i in range(0, len(recs), 50)]
    json.dump({"chunks": chunks}, f, ensure_ascii=False)
    print("attach chunks:", len(chunks))

# paywalled / blocked: DOIs with an OA url that failed verification, plus DOIs with no OA url at all
with open(os.path.join(HERE, "pdfs_paywalled.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["doi", "reason"])
    for d in dois:
        if d in ok:
            continue
        w.writerow([d, "oa-url-blocked" if d in oa else "no-oa-url"])
print("wrote pdfs_paywalled.csv")
