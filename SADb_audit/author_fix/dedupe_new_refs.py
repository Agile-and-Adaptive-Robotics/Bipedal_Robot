"""Dedupe the 33 new refs against Airtable DOIs; keep genuinely-new non-webpage ones."""
import json, os, urllib.request

PAT = os.environ["AT_PAT"]
new = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\new_refs_full.json", encoding="utf-8"))

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

recs = {}
offset = ""
hdr = {"Authorization": "Bearer " + PAT}
while True:
    q = "appMQTnobUNRytIp7/Papers?pageSize=100&fields%5B%5D=DOI"
    if offset:
        q += "&offset=" + offset
    req = urllib.request.Request("https://api.airtable.com/v0/" + q, headers=hdr)
    with urllib.request.urlopen(req, timeout=40) as r:
        page = json.loads(r.read().decode())
    for r in page.get("records", []):
        recs[norm_doi((r.get("fields") or {}).get("DOI"))] = r["id"]
    offset = page.get("offset")
    if not offset:
        break

fresh, dupe, webpage = [], [], []
for o in new:
    if o["itemType"] == "webpage":
        webpage.append(o)
        continue
    if o["doi"] and norm_doi(o["doi"]) in recs:
        dupe.append(o)
        continue
    fresh.append(o)

print("webpage (skip):", len(webpage), [w["key"] for w in webpage])
print("already in Airtable:", len(dupe), [d["key"] for d in dupe])
print("NEW to import:", len(fresh), [f["key"] for f in fresh])
json.dump(fresh, open(r"D:\Github\Bipedal_Robot\SADb_audit\new_refs_fresh.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
