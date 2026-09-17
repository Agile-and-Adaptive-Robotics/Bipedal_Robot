"""Attach local Zotero PDFs of the 19 new AARL items to their Airtable records.
Copy child PDF from Zotero storage into D:\sadb_pdf_staging\tmpattach (the served
root), attach via the live tunnel, verify, then delete the temp copy.
Env: AT_PAT, TUNNEL.
"""
import csv, json, os, shutil, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
AT_PAT = os.environ["AT_PAT"]
TUNNEL = os.environ["TUNNEL"].rstrip("/")
STAGING = r"D:\sadb_pdf_staging"
TMP = os.path.join(STAGING, "tmpattach")
STORAGE = r"C:\Users\Ben Bolen\Zotero\storage"
BASE = "appMQTnobUNRytIp7"

fresh = json.load(open(os.path.join(HERE, "..", "new_refs_fresh.json"), encoding="utf-8"))
RECORD_IDS = {
    "JMJIQ9ZZ": "rec7FBx5woaWolX0i", "6CEVM6YE": "rec00Fih3nO6FmfNV",
    "6NIJ7VIE": "rec3fLgtC69pKIljl", "8YCE28D9": "recDoc59f9eBJy5Af",
    "4HP4LJDQ": "rec54bBKJSjsVkg60", "6Z7RBVRC": "recVOknfb17E9AxjJ",
    "PXUSZ7NR": "recRtx7hS4xPDLv4P", "ZYA3GGJF": "recJhNrjI9iPIJB1E",
    "SUEWLVPW": "rec4QLIO1aXeQFced", "E3U29CXK": "recwy6OuthAArEsQ9",
    "R29RM7JI": "recx0HOyhL3E7SkL2", "WGZSJTQH": "rec00oIejxpXTzNt9",
    "TQYBQQUG": "recKcZ0SBPTxJu485", "TBHKUUAK": "recRqtgGMEHCYFJ4t",
    "FNQCAAU7": "recPnKgpWB1oLiyWJ", "4S96Y3FG": "recUDAo88KmPkHNdT",
    "BLL4S52K": "recFsVAudLXpREicY", "IZ2IA2ML": "recWgf8LP7dtq1UPZ",
    "7ZJV4GTK": "rec92HUnN53mhFrCc",
}
for o in fresh:
    o["record_id"] = RECORD_IDS.get(o["key"], "")

def air(method, path, body=None):
    req = urllib.request.Request(
        ("https://api.airtable.com/v0/" + path) if not path.startswith("http") else path,
        data=json.dumps(body).encode() if body is not None else None,
        headers={"Authorization": "Bearer " + AT_PAT, "Content-Type": "application/json"},
        method=method)
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())

os.makedirs(TMP, exist_ok=True)
log = []
ok = 0
for o in fresh:
    # find record id by DOI (or title match fallback handled manually later)
    rid = o.get("record_id", "")
    if not rid:
        log.append((o["key"], "NO-RECORD-ID", ""))
        continue
    # find local pdf child
    try:
        with urllib.request.urlopen(f"http://localhost:23119/api/groups/735051/items/{o['key']}/children?format=json", timeout=30) as r:
            kids = json.loads(r.read().decode())
    except Exception as e:
        log.append((o["key"], "CHILD-FAIL", str(e)[:100]))
        continue
    local = None
    for k in kids:
        if k["data"].get("contentType") == "application/pdf" and k["data"].get("linkMode") in ("imported_file", "imported_url"):
            cand = os.path.join(STORAGE, k["key"], k["data"].get("filename") or "")
            if os.path.exists(cand):
                local = cand
                break
    if not local:
        log.append((o["key"], rid, "NO-LOCAL-PDF", ""))
        continue
    tmpname = o["key"] + ".pdf"
    tmp = os.path.join(TMP, tmpname)
    shutil.copy(local, tmp)
    url = TUNNEL + "/tmpattach/" + tmpname
    try:
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
            ok += 1
            os.remove(tmp)
            log.append((o["key"], rid, "ok size=%d" % size, fname_tmp := tmpname))
        else:
            log.append((o["key"], rid, "NO-FETCH size=%d" % size, tmpname))
    except Exception as e:
        log.append((o["key"], rid, "FAIL", str(e)[:150]))
    time.sleep(0.4)

with open(os.path.join(HERE, "newref_attach_status.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["key", "record_id", "status", "detail"])
    w.writerows(log)
print("attached %d/%d" % (ok, len(fresh)))
for l in log:
    print("  ", l)
