"""TASK 3: file.io bridge — upload staged PDFs for unattached corpus records,
attach the expiring link to the Airtable record (Airtable fetches its own copy),
then delete the staged file from the repository.
Env: AT_PAT. Derivation mirrors upload_pdfs.py (same skip-empty/skip-hasfile rules).
"""
import csv, json, os, time, urllib.request, urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)
STAGING = os.path.join(SAD, "pdf_staging")
BASE = "appMQTnobUNRytIp7"
AT_PAT = os.environ["AT_PAT"]

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
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())

# corpus targets with staged files
targets = []
with open(os.path.join(HERE, "pdfs_from_zotero.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        targets.append(norm_doi(r["doi"]))
inv = {}
with open(os.path.join(SAD, "pdf_inventory.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r["file_exists"].lower() == "true" and r["doi"]:
            inv.setdefault(r["doi"], []).append((r["library"], r["local_path"], r["surname"], r["year"]))

# table state
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

todo = []
reasons = {}
for d in targets:
    if d not in recs:
        reasons["no-record"] = reasons.get("no-record", 0) + 1
        continue
    if recs[d][1]:
        reasons["has-attachments"] = reasons.get("has-attachments", 0) + 1
        continue
    if d not in inv:
        reasons["no-staged-file"] = reasons.get("no-staged-file", 0) + 1
        continue
    reasons["todo"] = reasons.get("todo", 0) + 1
    pref = [c for c in inv[d] if c[0] == "users/0"] or inv[d]
    todo.append((recs[d][0], d, pref[0][1]))
print("DEBUG targets:", len(targets), "recs:", len(recs), "inv:", len(inv))
print("reasons:", reasons)

def file_io_upload(path):
    # multipart POST via curl (python multipart is messy on 3.14 stdlib)
    import subprocess
    p = subprocess.run(["curl.exe", "-s", "-m", "120", "-F", f"file=@{path}", "https://file.io"],
                       capture_output=True)
    j = json.loads(p.stdout.decode())
    return j.get("link"), j.get("status"), j.get("message", "")

log = []
ok = 0
t0 = time.time()
for i, (rid, d, path) in enumerate(todo):
    name = os.path.basename(path).split("__")[0] + ".pdf"
    try:
        link, st, msg = file_io_upload(path)
        if not link:
            log.append((d, rid, "FILEIO-FAIL", (msg or st or "")[:120]))
            time.sleep(2)
            continue
        air("PATCH", f"{BASE}/Papers/{rid}", body={"fields": {"Attachments": [{"url": link}]}})
        log.append((d, rid, "ok", name))
        ok += 1
        os.remove(path)  # delete from repository once Airtable accepted it
        if ok % 15 == 0:
            print("attached %d/%d (%.0fs)" % (ok, len(todo), time.time() - t0))
    except urllib.error.HTTPError as e:
        log.append((d, rid, "HTTP-FAIL", str(e.code) + " " + e.read().decode()[:150]))
    except Exception as e:
        log.append((d, rid, "FAIL", str(e)[:150]))
    time.sleep(0.6)

with open(os.path.join(HERE, "fileio_status.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["doi", "record_id", "status", "detail"])
    w.writerows(log)
print("DONE: attached %d/%d in %.0fs -> fileio_status.csv" % (ok, len(todo), time.time() - t0))
