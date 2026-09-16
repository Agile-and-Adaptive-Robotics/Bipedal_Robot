"""Attach Ben's 17 found PDFs to their Airtable records via the tunnel bridge,
move consumed PDFs to pdfs_consumed, and mark paperhunt_status.csv.
"""
import csv, json, os, shutil, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
PDFS = r"D:\sadb_pdf_staging\pdfs"
CONSUMED = r"D:\sadb_pdf_staging\pdfs_consumed"
BASE = "appMQTnobUNRytIp7"
AT_PAT = os.environ["AT_PAT"]
TUNNEL = os.environ["TUNNEL"].rstrip("/")

def air(method, path, body=None):
    req = urllib.request.Request(
        ("https://api.airtable.com/v0/" + path) if not path.startswith("http") else path,
        data=json.dumps(body).encode() if body is not None else None,
        headers={"Authorization": "Bearer " + AT_PAT, "Content-Type": "application/json"},
        method=method)
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())

os.makedirs(CONSUMED, exist_ok=True)
pairs = []
with open(os.path.join(HERE, "pdf_match_pairs.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        pairs.append(r)

log = []
ok = 0
for p in pairs:
    fname = p["pdf_file"]
    src = os.path.join(PDFS, fname)
    url = TUNNEL + "/pdfs/" + urllib.parse.quote(fname)
    try:
        air("PATCH", f"{BASE}/Papers/{p['doi']}", body=None) if False else None
        # need record id: paperhunt rows carry airtable_id; pairs CSV doesn't. Look up by DOI.
        q = f"{BASE}/Papers?pageSize=5&filterByFormula=" + urllib.parse.quote("{" + "DOI" + "}='" + p["doi"] + "'")
        recs = air("GET", q).get("records", [])
        if not recs:
            log.append((p["doi"], "NO-RECORD", fname))
            continue
        rid = recs[0]["id"]
        air("PATCH", f"{BASE}/Papers/{rid}", body={"fields": {"Attachments": [{"url": url}]}})
        size = 0
        for _ in range(8):
            time.sleep(2)
            chk = air("GET", f"{BASE}/Papers/{rid}")
            atts = (chk.get("fields") or {}).get("Attachments") or []
            size = max((a.get("size") or 0) for a in atts) if atts else 0
            if size > 1000:
                break
        if size > 1000:
            shutil.move(src, os.path.join(CONSUMED, fname))
            ok += 1
            log.append((p["doi"], rid, "ok size=%d" % size, fname))
        else:
            log.append((p["doi"], rid, "NO-FETCH size=%d" % size, fname))
    except Exception as e:
        log.append((p["doi"], "?", "FAIL", str(e)[:150]))
    time.sleep(0.4)

print("attached %d/%d" % (ok, len(pairs)))
for l in log:
    print("  ", l)

# mark paperhunt_status.csv
path = os.path.join(HERE, "paperhunt_status.csv")
rows = list(csv.DictReader(open(path, encoding="utf-8-sig")))
okdois = {p["doi"] for p in pairs}
for r in rows:
    if r["doi"] in okdois:
        r["Downloaded"] = "x -> attached"
fieldnames = list(rows[0].keys())
if "Attached" not in fieldnames:
    fieldnames.append("Attached")
for r in rows:
    if r["doi"] in okdois:
        r["Attached"] = "y 2026-09-16"
with open(path, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)
print("paperhunt_status.csv updated")
