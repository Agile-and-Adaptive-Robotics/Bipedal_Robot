"""Retry the 48 attach-failed records from pass 2, capture Airtable's error body,
and stage verified PDFs to D:\\sadb_pdf_staging for drag-drop."""
import csv, json, os, re, time, urllib.parse, urllib.request

SAD = os.path.dirname(os.path.abspath(__file__))
STAGE = r"D:\sadb_pdf_staging"
MAILTO = "benjamin.bolen@pdx.edu"
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 SADb/1.0 (mailto:%s)" % MAILTO}
os.makedirs(STAGE, exist_ok=True)

txt = open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()
PAT = re.search(r"PAT\s*:\s*(pat[^\s]+)", txt).group(1)
BASE = "appMQTnobUNRytIp7"


def air(method, path, body=None):
    req = urllib.request.Request("https://api.airtable.com/v0/" + path,
                                 data=json.dumps(body).encode() if body is not None else None,
                                 headers={"Authorization": "Bearer " + PAT, "Content-Type": "application/json"},
                                 method=method)
    try:
        with urllib.request.urlopen(req, timeout=90) as r:
            return json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        return {"_err": e.code, "_body": e.read().decode("utf-8", "replace")[:300]}


def norm(d):
    d = (d or "").strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d


# the 48 attach-failed + export metadata
failed = []
with open(os.path.join(SAD, "author_fix", "pdf_hunt_pass2_20260920.csv"), encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        if row["status"].startswith("attach-failed"):
            failed.append(norm(row["doi"]))
exp = json.load(open(os.path.join(SAD, "export", "sadb_export.json"), encoding="utf-8"))
meta = {}
for r in exp:
    d = norm(r["doi"])
    if d:
        meta[d] = r
print("retrying:", len(failed))

# candidates in one OpenAlex call (48 < 50)
filt = "doi:" + "|".join(failed)
url = ("https://api.openalex.org/works?filter=" + urllib.parse.quote(filt, safe="|")
       + "&select=doi,best_oa_location,locations&per-page=50&mailto=" + MAILTO)
req = urllib.request.Request(url, headers={"User-Agent": "sadb-pdfhunt/1.0 (mailto:%s)" % MAILTO})
data = json.loads(urllib.request.urlopen(req, timeout=60).read().decode())
cand = {}
for w in data.get("results", []):
    d = norm(w.get("doi") or "")
    urls = []
    for loc in [w.get("best_oa_location")] + (w.get("locations") or []):
        u = (loc or {}).get("pdf_url")
        if u and u not in urls:
            urls.append(u)
    if urls:
        cand[d] = urls


def verify(url):
    try:
        req = urllib.request.Request(url, headers={**UA, "Range": "bytes=0-63"})
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.read(64).startswith(b"%PDF")
    except Exception:
        return False


results = []
n_att = n_staged = 0
for d in failed:
    rid = meta.get(d, {}).get("id", "")
    surname = re.sub(r"[^\w]", "", (meta.get(d, {}).get("primary") or "Unknown").split()[-1])
    year = meta.get(d, {}).get("year") or "nd"
    stage_path = os.path.join(STAGE, f"{surname}_{year}__{d.replace('/', '_')}.pdf")
    status, err, staged = "", "", False
    good = ""
    for u in cand.get(d, []):
        if not verify(u):
            continue
        good = u
        # all hunt-list records have ZERO attachments (has_pdf=False), so no read-before-append
        res = air("PATCH", f"{BASE}/Papers/{rid}", {"fields": {"Attachments": [{"url": u}]}})
        if "_err" in res:
            status = f"attach-422: {res['_body'][:200]}"
        else:
            status = "ATTACHED"
            n_att += 1
        break
    if not good:
        status = status or "no-verified-pdf"
    # stage a local copy regardless (drag-drop path) if one isn't there
    if good and not os.path.exists(stage_path):
        try:
            req = urllib.request.Request(good, headers=UA)
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
            if body.startswith(b"%PDF"):
                open(stage_path, "wb").write(body)
                staged = True
                n_staged += 1
        except Exception as e:
            status += f" | stage-failed: {str(e)[:60]}"
    elif os.path.exists(stage_path):
        staged = True
    results.append([d, rid, status, good, stage_path if staged else ""])
    time.sleep(0.4)

out = os.path.join(SAD, "author_fix", "pdf_hunt_retry_20260920.csv")
with open(out, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["doi", "record_id", "status", "verified_url", "staged_path"])
    w.writerows(results)
print(f"DONE: {n_att} attached, {n_staged} staged locally -> {out}")
from collections import Counter
print(Counter(r[2].split(':')[0].split(' ')[0] for r in results))
