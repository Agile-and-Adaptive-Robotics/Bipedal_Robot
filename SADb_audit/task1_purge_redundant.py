# TASK 1: Delete PDF attachments in personal SAD collections (RBG8ZVYP + NZBY3DE7
# + F6RKEZNY + 2HJUIDGU) that are redundant because a copy exists in pdf_staging,
# the AARL group library, or Airtable. Uses Zotero web API (permanent delete).
# Keys via env ZOT_KEY (user 631450) + AT_PAT. Every deletion logged.
import csv, json, os, time, urllib.request, urllib.error

SAD = r"D:\Github\Bipedal_Robot\SADb_audit"
ZKEY = os.environ["ZOT_KEY"]
UID = "631450"
AT_PAT = os.environ.get("AT_PAT", "")
BASE = "appMQTnobUNRytIp7"
COLLS = ["RBG8ZVYP", "NZBY3DE7", "F6RKEZNY", "2HJUIDGU"]
STAGING = os.path.join(SAD, "pdf_staging")

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

def zreq(method, path, headers=None):
    url = f"https://api.zotero.org/users/{UID}{path}"
    h = {"Authorization": f"Bearer {ZKEY}"}
    h.update(headers or {})
    req = urllib.request.Request(url, headers=h, method=method)
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.status, r.read()

# 1) staged att keys (redundancy anchor 1: staging copy exists)
staged_keys = set()
for f in os.listdir(STAGING):
    if f.endswith(".pdf") and "__" in f:
        staged_keys.add(f.rsplit("__", 1)[1][:-4])

# 2) AARL group DOI->has-pdf map from local inventory (redundancy anchor 2)
aarl_doi = set()
with open(os.path.join(SAD, "pdf_inventory.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r["library"] == "groups/735051" and r["file_exists"].lower() == "true" and r["doi"]:
            aarl_doi.add(r["doi"])

# 3) Airtable DOI->has-file map (redundancy anchor 3)
at_has = set()
if AT_PAT:
    offset = ""
    hdr = {"Authorization": "Bearer " + AT_PAT}
    while True:
        q = BASEQ = f"https://api.airtable.com/v0/{BASE}/Papers?pageSize=100&fields%5B%5D=DOI&fields%5B%5D=Attachments"
        if offset:
            q += "&offset=" + offset
        req = urllib.request.Request(q, headers=hdr)
        with urllib.request.urlopen(req, timeout=40) as r:
            page = json.loads(r.read().decode())
        for rec in page.get("records", []):
            d = norm_doi((rec.get("fields") or {}).get("DOI"))
            if d and (rec.get("fields") or {}).get("Attachments"):
                at_has.add(d)
        offset = page.get("offset")
        if not offset:
            break
print("anchors: staged=%d aarl_dois=%d airtable_dois=%d" % (len(staged_keys), len(aarl_doi), len(at_has)))

# 4) enumerate personal SAD attachments (local API has full data incl. version)
cands = []
for coll in COLLS:
    start = 0
    while True:
        import subprocess
        p = subprocess.run(["curl.exe", "-s", "-m", "60",
            f"http://localhost:23119/api/users/0/collections/{coll}/items?format=json&limit=100&start={start}&itemType=attachment"],
            capture_output=True)
        items = json.loads(p.stdout.decode("utf-8", "replace")) if p.stdout else None
        if not items or not isinstance(items, list) or len(items) == 0:
            break
        for it in items:
            if it["data"].get("contentType") != "application/pdf":
                continue
            cands.append({"key": it["key"], "version": it["version"],
                          "parent": it["data"].get("parentItem", ""),
                          "filename": it["data"].get("filename", ""),
                          "coll": coll})
        start += 100
        if start > 3000:
            break
    time.sleep(0.2)
print("personal SAD pdf attachments:", len(cands))

# parent DOIs: paginate ALL personal top-level items once (reliable)
pmap = {}
start = 0
while True:
    p = subprocess.run(["curl.exe", "-s", "-m", "60",
        f"http://localhost:23119/api/users/0/items?format=json&limit=100&start={start}&itemType=-attachment"],
        capture_output=True)
    items = json.loads(p.stdout.decode("utf-8", "replace")) if p.stdout else None
    if not items or not isinstance(items, list) or len(items) == 0:
        break
    for it in items:
        if it["data"].get("parentItem"):
            continue
        pmap[it["key"]] = norm_doi(it["data"].get("DOI"))
    start += 100
    if start > 5000:
        break
print("personal top-level items mapped:", len(pmap), "| parents resolved:", sum(1 for c in cands if c["parent"] in pmap))

# 5) decide + delete
log = []
deleted = 0
kept = 0
for c in cands:
    d = pmap.get(c["parent"], "")
    in_staging = c["key"] in staged_keys
    in_aarl = d in aarl_doi and bool(d)
    in_air = d in at_has and bool(d)
    if not (in_staging or in_aarl or in_air):
        kept += 1
        log.append((c["key"], c["coll"], d, "KEPT (no redundant copy)"))
        continue
    try:
        status, _ = zreq("DELETE", f"/items/{c['key']}",
                         headers={"If-Unmodified-Since-Version": str(c["version"])})
        deleted += 1
        log.append((c["key"], c["coll"], d, f"DELETED via {in_staging and 'staging' or ''}{in_aarl and '+AARL' or ''}{in_air and '+Airtable' or ''}"))
    except urllib.error.HTTPError as e:
        log.append((c["key"], c["coll"], d, f"FAIL HTTP {e.code} {e.read().decode()[:120]}"))
    except Exception as e:
        log.append((c["key"], c["coll"], d, f"FAIL {str(e)[:120]}"))
    time.sleep(0.15)

with open(os.path.join(SAD, "task1_purge_log.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["att_key", "collection", "doi", "action"])
    w.writerows(log)
print(f"DONE: deleted {deleted}, kept {kept} -> task1_purge_log.csv")
