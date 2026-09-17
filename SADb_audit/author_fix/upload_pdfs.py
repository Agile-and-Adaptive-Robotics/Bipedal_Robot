"""Upload staged Zotero PDFs to Airtable attachments — full auto.

Flow per record (Airtable Attachment Upload API):
  POST /v0/{base}/attachments            -> {id, uploadUrl}
  PUT  uploadUrl (raw bytes)             -> 200
  PATCH /v0/{base}/Papers/{recordId}     -> Attachments = [{"id": id}]

Safety:
  - only corpus DOIs from pdfs_from_zotero.csv (sensory afferent / motor control scope)
  - only records whose Attachments field is currently EMPTY (no overwrites)
  - case-insensitive DOI match against a fresh full-table pull
Token comes from env AT_PAT; never written to disk.
"""
import csv, json, os, sys, time, urllib.request

BASE = "appMQTnobUNRytIp7"
TABLE = "Papers"
PAT = os.environ.get("AT_PAT", "")
HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)
STAGING = os.path.join(SAD, "pdf_staging")

def api(method, path, body=None, raw=None, raw_ctype=None):
    url = path if path.startswith("http") else "https://api.airtable.com/v0/" + path
    data = None
    headers = {"Authorization": "Bearer " + PAT}
    if body is not None:
        data = json.dumps(body).encode()
        headers["Content-Type"] = "application/json"
    if raw is not None:
        data = raw
        headers["Content-Type"] = raw_ctype or "application/pdf"
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    with urllib.request.urlopen(req, timeout=60) as r:
        payload = r.read()
    return json.loads(payload) if (payload and not raw) else (payload or b"")

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

# 1) validate token
me = api("GET", "meta/whoami")
print("token ok, user:", me.get("id"))

# 2) corpus list (279)
targets = []
with open(os.path.join(HERE, "pdfs_from_zotero.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        targets.append(norm_doi(r["doi"]))
print("corpus targets with zotero pdf:", len(targets))

# 3) full table pull: id, DOI, has-attachments?
recs = {}
offset = ""
while True:
    q = "?pageSize=100&fields%5B%5D=" + "DOI" + "&fields%5B%5D=" + "Attachments"
    if offset:
        q += "&offset=" + offset
    page = api("GET", BASE + "/Papers" + q)
    for r in page.get("records", []):
        recs[r["id"]] = {"doi": norm_doi((r.get("fields") or {}).get("DOI")),
                         "atts": (r.get("fields") or {}).get("Attachments") or []}
    offset = page.get("offset")
    if not offset:
        break
print("table records pulled:", len(recs))

# 4) resolve: corpus target + empty attachments
todo = []
skipped_hasfile = 0
doi_to_ids = {}
for rid, v in recs.items():
    if v["doi"]:
        doi_to_ids.setdefault(v["doi"], []).append(rid)
for d in targets:
    ids = doi_to_ids.get(d, [])
    if not ids:
        print("  WARN no record for", d)
        continue
    rid = ids[0]
    if recs[rid]["atts"]:
        skipped_hasfile += 1
        continue
    todo.append((rid, d))
print("to upload:", len(todo), "| skipped (already has files):", skipped_hasfile)

# 5) staged-file map: doi -> local path (prefer personal library copies first)
inv = {}
with open(os.path.join(SAD, "pdf_inventory.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r["file_exists"].lower() == "true" and r["doi"]:
            inv.setdefault(r["doi"], []).append((r["library"], r["local_path"], r["surname"], r["year"]))

def staged_file(d):
    cands = inv.get(d) or []
    if not cands:
        return None
    pref = [c for c in cands if c[0] == "users/0"] or cands
    return pref[0][1]

# 6) upload loop
status = []
done = 0
t0 = time.time()
for rid, d in todo:
    path = staged_file(d)
    if not path or not os.path.exists(path):
        status.append((d, rid, "no-local-file", ""))
        continue
    surname_year = os.path.splitext(os.path.basename(path))[0]
    surname_year = surname_year.split("__")[0] + ".pdf"
    try:
        up = api("POST", BASE + "/attachments",
                 body={"contentType": "application/pdf", "filename": surname_year})
        api("PUT", up["uploadUrl"], raw=open(path, "rb").read())
        api("PATCH", BASE + "/Papers/" + rid,
            body={"fields": {"Attachments": [{"id": up["id"]}]}})
        status.append((d, rid, "ok", surname_year))
        done += 1
        if done % 20 == 0:
            print("uploaded %d/%d (%.0fs)" % (done, len(todo), time.time() - t0))
    except Exception as e:
        status.append((d, rid, "FAIL", str(e)[:200]))
    time.sleep(0.4)

with open(os.path.join(HERE, "upload_status.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["doi", "record_id", "status", "filename_or_error"])
    w.writerows(status)
print("DONE: uploaded %d/%d in %.0fs -> upload_status.csv" % (done, len(todo), time.time() - t0))
