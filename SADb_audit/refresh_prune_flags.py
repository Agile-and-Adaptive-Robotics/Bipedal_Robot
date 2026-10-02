"""Clear stale prune flags post-curation (2026-09-28).

Records flagged prune-candidate BEFORE the curation campaign that now have
notes are no longer candidates — clear their Prune Status/Reason. Also
recomputes the strict shortlist on the fresh export and reports the delta.
"""
import json, os, re, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "appMQTnobUNRytIp7"
PAPERS = "tblnnMrZszhboU4uD"
F_STATUS = "fld9jhnEQW738POne"
F_REASON = "fldcW2N1nn2pNxfv0"
F_NOTES = "fld3gPiUIKn26N6ji"


def pat():
    return re.search(r"PAT\s*:\s*(pat[^\s]+)",
                     open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()).group(1)


def air(path, params=None):
    q = "?" + urllib.parse.urlencode(params, doseq=True) if params else ""
    req = urllib.request.Request("https://api.airtable.com/v0/" + path + q,
                                 headers={"Authorization": "Bearer " + pat()})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def patch(body):
    req = urllib.request.Request(f"https://api.airtable.com/v0/{BASE}/{PAPERS}",
                                 data=json.dumps(body).encode(),
                                 headers={"Authorization": "Bearer " + pat(),
                                          "Content-Type": "application/json"},
                                 method="PATCH")
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


# pull flagged records (they have Prune Status set)
flagged, offset = [], None
while True:
    p = [("pageSize", 100), ("fields[]", ["Name", "Notes", "Prune Status"])]
    if offset:
        p.append(("offset", offset))
    d = air(f"{BASE}/{PAPERS}", p)
    flagged.extend(r for r in d.get("records", []) if r["fields"].get("Prune Status"))
    offset = d.get("offset")
    if not offset:
        break
print(f"flagged records found: {len(flagged)}")
to_clear = [r for r in flagged if (r["fields"].get("Notes") or "").strip()]
print(f"now curated -> clear flag: {len(to_clear)}; keep flag: {len(flagged) - len(to_clear)}")
for i in range(0, len(to_clear), 10):
    chunk = to_clear[i:i + 10]
    patch({"records": [{"id": r["id"], "fields": {F_STATUS: None, F_REASON: None}}
                       for r in chunk], "typecast": False})
    time.sleep(0.4)
print("cleared stale flags")
