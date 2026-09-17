"""TASK 3 v2: attach staged PDFs to Airtable via trycloudflare tunnel URL.
Airtable fetches from the public tunnel -> stores its own copy.
Deletes the local file after a verified attach. Env: AT_PAT. Env LIMIT: max records.
"""
import csv, json, os, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)
STAGING = r"D:\sadb_pdf_staging"
BASE = "appMQTnobUNRytIp7"
AT_PAT = os.environ["AT_PAT"]
TUNNEL = os.environ.get("TUNNEL", "").rstrip("/")
LIMIT = int(os.environ.get("LIMIT", "999999"))

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

def air(method, path, body=None):
    req = urllib.request.Request(
        ("https://api.airtable.com/v0/" + path) if not path.startswith("http") else path,
        data=json.dumps(body).encode() if body is not None else None,
        headers={"Authorization": "Bearer " + AT_PAT, "Content-Type": "application/json"},
        method=method)
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())

targets = []
with open(os.path.join(HERE, "pdfs_from_zotero.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        targets.append(norm_doi(r["doi"]))

recs = {}
offset = ""
while True:
    q = f"{BASE}/Papers?pageSize=100&fields%5B%5D=DOI&fields%5B%5D=Attachments"
    if offset:
        q += "&offset=" + offset
    page = air("GET", q)
    for r in page.get("records", []):
        recs[norm_doi((r.get("fields") or {}).get("DOI"))] = (r["id"], (r.get("fields") or {}).get("Attachments") or [])
    offset = page.get("offset")
    if not offset:
        break

inv = {}
with open(os.path.join(SAD, "pdf_inventory.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r["file_exists"].lower() == "true" and r["doi"]:
            inv.setdefault(r["doi"], []).append((r["library"], r["local_path"], r["surname"], r["year"]))

todo = []
for d in targets:
    if d not in recs or d not in inv:
        continue
    rid, atts = recs[d]
    if any((a.get("size") or 0) > 1000 for a in atts):
        continue  # already has a real file
    pref = [c for c in inv[d] if c[0] == "users/0"] or inv[d]
    local = pref[0][1]
    if not os.path.exists(local):
        continue
    todo.append((rid, d, local))
todo = todo[:LIMIT]
print("todo:", len(todo))

# staged copies are named Surname_Year__ATTKEY.pdf — map att_key -> staging path
staged_by_key = {}
for f in os.listdir(STAGING):
    if "__" in f and f.endswith(".pdf"):
        staged_by_key[f.rsplit("__", 1)[1][:-4]] = os.path.join(STAGING, f)
# doi -> att_key from the inventory
doi_to_key = {}
with open(os.path.join(SAD, "pdf_inventory.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r["file_exists"].lower() == "true" and r["doi"]:
            doi_to_key[r["doi"]] = r["att_key"]

log = []
ok = 0
t0 = time.time()
for rid, d, local in todo:
    attkey = doi_to_key.get(d, "")
    staged = staged_by_key.get(attkey)
    if not staged or not os.path.exists(staged):
        log.append((d, rid, "NO-STAGING-COPY", attkey))
        continue
    fname = os.path.basename(staged)
    url = TUNNEL + "/" + urllib.parse.quote(fname)
    try:
        resp = air("PATCH", f"{BASE}/Papers/{rid}",
                   body={"fields": {"Attachments": [{"url": url}]}})
        # size/filename populate asynchronously after Airtable's fetch — poll
        size = 0
        for _ in range(6):
            time.sleep(2)
            chk = air("GET", f"{BASE}/Papers/{rid}")
            atts = (chk.get("fields") or {}).get("Attachments") or []
            size = atts[0].get("size", 0) if atts else 0
            if size and size > 1000:
                break
        if size and size > 1000:
            ok += 1
            try:
                os.remove(local)
            except OSError:
                pass
            log.append((d, rid, "ok", f"size={size}"))
        else:
            log.append((d, rid, "NO-FETCH", f"size={size}"))
    except urllib.error.HTTPError as e:
        log.append((d, rid, "HTTP-FAIL", str(e.code) + " " + e.read().decode()[:150]))
    except Exception as e:
        log.append((d, rid, "FAIL", str(e)[:150]))
    time.sleep(0.3)

with open(os.path.join(HERE, "tunnel_status.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["doi", "record_id", "status", "detail"])
    w.writerows(log)
print("DONE: attached %d/%d in %.0fs -> tunnel_status.csv" % (ok, len(todo), time.time() - t0))
