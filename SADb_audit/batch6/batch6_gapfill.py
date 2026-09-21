"""Batch 6 gap-fill + vocabulary dump.
(a) Crossref abstracts for the two Elsevier chapters; Europe PMC full abstract for Hatz 2012.
(b) Models / Review Papers name lists (dedupe check) + Feedback id->name table.
"""
import json, os, re, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))

def get_json(url, headers=None):
    req = urllib.request.Request(url, headers=headers or {"User-Agent": "SADb-curation/1.0 (mailto:research@example.org)"})
    with urllib.request.urlopen(req, timeout=45) as r:
        return json.loads(r.read().decode("utf-8", "replace"))

for label, doi in [("Thiry 2020", "10.1016/b978-0-12-816477-8.00011-9"),
                   ("Prochazka 2017", "10.1016/b978-0-12-803766-9.00008-7")]:
    try:
        d = get_json("https://api.crossref.org/works/" + doi)["message"]
        ab = d.get("abstract", "(no abstract in Crossref)")
        ab = re.sub(r"<[^>]+>", " ", ab)
        print(f"=== CROSSREF {label} | container={d.get('container-title','')} ===\n{ab}\n")
    except Exception as e:
        print(f"=== CROSSREF {label} FAILED: {e} ===")

try:
    q = urllib.parse.quote('DOI:"10.1152/jn.00944.2011"')
    d = get_json(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json&resultType=core")
    hit = d["resultList"]["result"][0]
    print("=== EPMC Hatz 2012 full abstract ===")
    print(hit.get("abstractText", "(none)"))
except Exception as e:
    print(f"EPMC Hatz failed: {e}")

pat = re.search(r"PAT\s*:\s*(pat[^\s]+)", open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()).group(1)

def air_list(table):
    out, off = {}, None
    while True:
        url = f"https://api.airtable.com/v0/appMQTnobUNRytIp7/{table}?pageSize=100"
        if off:
            url += "&offset=" + off
        req = urllib.request.Request(url, headers={"Authorization": "Bearer " + pat})
        d = json.loads(urllib.request.urlopen(req, timeout=60).read().decode())
        for rec in d["records"]:
            out[rec["id"]] = rec["fields"].get("Name", "")
        off = d.get("offset")
        if not off:
            return out

fb = air_list("Feedback")
print("\n=== FEEDBACK id -> name ===")
for i, n in fb.items():
    print(f"{i}  {n}")

for table, label in [("Models", "MODELS"), ("Review%20Papers", "REVIEW PAPERS")]:
    names = air_list(table)
    hits = [n for n in names.values()
            if re.search(r"ekeberg|fujiki|c[oô]t[eé]|prochazka|proch|thiry|nichols|hatz|cazalets|boulland", n, re.I)]
    print(f"\n=== {label}: {len(names)} records; relevant existing: {sorted(hits)}")
