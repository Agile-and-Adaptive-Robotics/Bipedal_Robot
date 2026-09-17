import csv, json, os, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)
PAT = os.environ["AT_PAT"]

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

corpus = {}
for f in (os.path.join(SAD, "airtable_papers_slim.csv"),
          os.path.join(SAD, "airtable_rest_import_clean.csv")):
    with open(f, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            d = norm_doi(r.get("DOI"))
            if d:
                corpus[d] = True
with open(os.path.join(SAD, "batch3", "rr_digest_shortlist.json"), encoding="utf-8") as f:
    for s in json.load(f):
        corpus[norm_doi(s["doi"])] = True
corpus[norm_doi("10.7554/elife.107480")] = True

recs = {}
offset = ""
hdr = {"Authorization": "Bearer " + PAT}
while True:
    q = "appMQTnobUNRytIp7/Papers?pageSize=100&fields%5B%5D=DOI&fields%5B%5D=Attachments"
    if offset:
        q += "&offset=" + offset
    req = urllib.request.Request("https://api.airtable.com/v0/" + q, headers=hdr)
    with urllib.request.urlopen(req, timeout=40) as r:
        page = json.loads(r.read().decode())
    for r in page.get("records", []):
        recs[norm_doi((r.get("fields") or {}).get("DOI"))] = (r.get("fields") or {}).get("Attachments") or []
    offset = page.get("offset")
    if not offset:
        break

with_att = sum(1 for d in corpus if d in recs and any((a.get("size") or 0) > 1000 for a in recs[d]))
stubs = [d for d in corpus if d in recs and recs[d] and not any((a.get("size") or 0) > 1000 for a in recs[d])]
none_at_all = [d for d in corpus if d not in recs or not recs[d]]
print("corpus records:", len(corpus))
print("with real attachment (>1KB):", with_att)
print("with stub/empty-size attachment:", len(stubs))
print("with no attachment field at all:", len(none_at_all))
with open(os.path.join(HERE, "remaining_no_pdf.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["doi", "state"])
    for d in stubs:
        w.writerow([d, "stub"])
    for d in none_at_all:
        w.writerow([d, "none"])
print("wrote remaining_no_pdf.csv")
