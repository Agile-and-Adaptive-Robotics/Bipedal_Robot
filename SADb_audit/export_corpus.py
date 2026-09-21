"""Export the full Airtable Papers corpus to export/sadb_export.json + .csv.

Reads the PAT from D:\Github\api_credentials_local.txt (never stores it).
Stdlib only. Run: python export_corpus.py
Regenerate after any curation batch — the HTML app + pivot consume these files.

NOTE: Airtable's REST API keys response `fields` by FIELD NAME (even when you
request by field id). "Models copy" has two same-named fields; fldh983rtt2YtMZQX
is empty everywhere (2026-09-11 audit) so any value seen belongs to the live
Review Papers link fldBLowhKJcjuFMSK.
"""
import csv, json, os, re, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "export")
os.makedirs(OUT, exist_ok=True)
BASE = "appMQTnobUNRytIp7"
PAPERS = "tblnnMrZszhboU4uD"

FETCH = ["Name", "Notes", "Attachments", "Feedback", "Primary Author", "Animal",
         "Year", "Models", "Models 2", "Models copy", "DOI", "Secondary Authors"]


def pat():
    txt = open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()
    m = re.search(r"PAT\s*:\s*(pat[^\s]+)", txt)
    if not m:
        raise SystemExit("PAT not found in credentials file")
    return m.group(1)


_PAT = pat()


def air(path, params=None):
    q = ""
    if params:
        q = "?" + urllib.parse.urlencode(params, doseq=True, quote_via=urllib.parse.quote)
    req = urllib.request.Request("https://api.airtable.com/v0/" + path + q,
                                 headers={"Authorization": "Bearer " + _PAT})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def names_for(table_id, label):
    out, offset = {}, None
    while True:
        p = [("pageSize", 100)]
        if offset:
            p.append(("offset", offset))
        d = air(f"{BASE}/{table_id}", p)
        for rec in d.get("records", []):
            out[rec["id"]] = rec["fields"].get("Name", "")
        offset = d.get("offset")
        if not offset:
            break
    print(f"{label}: {len(out)} records")
    return out


def air_post(path, body):
    req = urllib.request.Request("https://api.airtable.com/v0/" + path,
                                 data=json.dumps(body).encode(),
                                 headers={"Authorization": "Bearer " + _PAT,
                                          "Content-Type": "application/json"},
                                 method="POST")
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


papers_raw, offset = [], None
while True:
    p = [("pageSize", 100)]
    if offset:
        p.append(("offset", offset))
    d = air(f"{BASE}/{PAPERS}", p)  # full fetch: a fields[] filter 422s here
    papers_raw.extend(d.get("records", []))  # because two fields share the
    offset = d.get("offset")                 # name "Models copy"
    if not offset:
        break
print(f"Papers: {len(papers_raw)} records")

fb_names = names_for("tblot5mo4s5KgN5le", "Feedback")
md_names = names_for("tblsBq9IEv7dZe6fn", "Models")
rv_names = names_for("tblSEubKcRId4wYMK", "Review Papers")


def sel(v):
    if isinstance(v, list):
        return [x.get("name", "") if isinstance(x, dict) else str(x) for x in v]
    if isinstance(v, dict):
        return v.get("name", "")
    return v


out = []
for rec in papers_raw:
    f = rec.get("fields", {})
    notes = f.get("Notes", "") or ""
    out.append({
        "id": rec["id"],
        "title": f.get("Name", ""),
        "primary": f.get("Primary Author", ""),
        "secondary": sel(f.get("Secondary Authors", [])) if isinstance(f.get("Secondary Authors"), list) else [],
        "year": f.get("Year", ""),
        "doi": (f.get("DOI", "") or "").strip(),
        "animals": sel(f.get("Animal", [])) if isinstance(f.get("Animal"), list) else [],
        "feedback": [fb_names.get(r, r) for r in f.get("Feedback", [])],
        "models_ref": [md_names.get(r, r) for r in f.get("Models", [])],
        "models2": [md_names.get(r, r) for r in f.get("Models 2", [])],
        "reviews": [rv_names.get(r, r) for r in f.get("Models copy", [])],
        "has_pdf": bool(f.get("Attachments")),
        "has_notes": bool(notes.strip()),
        "notes": notes.strip()[:1500],
    })

with open(os.path.join(OUT, "sadb_export.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False)

cols = ["id", "title", "primary", "secondary", "year", "doi", "animals", "feedback",
        "models_ref", "models2", "reviews", "has_pdf", "has_notes", "notes"]
with open(os.path.join(OUT, "sadb_export.csv"), "w", encoding="utf-8-sig", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    for r in out:
        w.writerow({c: (("; ".join(r[c]) if isinstance(r[c], list) else r[c])) for c in cols})

n_notes = sum(1 for r in out if r["has_notes"])
n_pdf = sum(1 for r in out if r["has_pdf"])
n_rev = sum(1 for r in out if r["reviews"])
n_an = sum(1 for r in out if r["animals"])
n_fb = sum(1 for r in out if r["feedback"])
print(f"wrote {len(out)} records | notes {n_notes} | pdf {n_pdf} | animals {n_an} | feedback {n_fb} | review-links {n_rev}")
