import csv, json, os, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)
AT_PAT = os.environ.get("AT_PAT", "")

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

targets = []
with open(os.path.join(HERE, "pdfs_from_zotero.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        targets.append(norm_doi(r["doi"]))
print("targets:", len(targets), "sample:", targets[:3])

inv = {}
with open(os.path.join(SAD, "pdf_inventory.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r["file_exists"].lower() == "true" and r["doi"]:
            inv.setdefault(r["doi"], []).append(r["local_path"])
print("inventory dois:", len(inv), "targets with staged file:", sum(1 for d in targets if d in inv))

hdr = {"Authorization": "Bearer " + AT_PAT}
recs = {}
offset = ""
while True:
    q = f"appMQTnobUNRytIp7/Papers?pageSize=100&fields%5B%5D=DOI&fields%5B%5D=Attachments"
    if offset:
        q += "&offset=" + offset
    req = urllib.request.Request("https://api.airtable.com/v0/" + q, headers=hdr)
    with urllib.request.urlopen(req, timeout=40) as r:
        page = json.loads(r.read().decode())
    for r in page.get("records", []):
        recs[norm_doi((r.get("fields") or {}).get("DOI"))] = (r["id"], (r.get("fields") or {}).get("Attachments") or [])
    offset = page.get("offset")
    if not offset:
        break
print("table dois:", len(recs))

n_ids = sum(1 for d in targets if d in recs)
n_att = sum(1 for d in targets if d in recs and recs[d][1])
n_staged = sum(1 for d in targets if d in inv)
print("targets with record:", n_ids, "| with existing atts:", n_att, "| with staged:", n_staged)
todo = [d for d in targets if d in recs and not recs[d][1] and d in inv]
print("actual todo:", len(todo))
