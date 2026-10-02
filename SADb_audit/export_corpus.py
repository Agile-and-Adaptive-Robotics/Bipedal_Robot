"""Export the full Airtable Papers corpus to export/sadb_export.json + .csv.

Reads the PAT from D:\Github\api_credentials_local.txt (never stores it).
Stdlib only. Run: python export_corpus.py
Regenerate after any curation batch — the HTML app + pivot consume these files.

2026-09-28 schema v2: the duplicate "Models copy" fields were renamed
(fldBLowhKJcjuFMSK is now "Review Papers"; the empty twin is "(unused) …"),
which UNIQUE-IFIES all field names — a fields[] filter no longer 422s, so this
fetch is now projection-filtered. New curation-layer fields exported:
Afferent Types, Animal Study Potential, Robot/Sim Translation, Prune Status.
"""
import csv, json, os, re, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "export")
os.makedirs(OUT, exist_ok=True)
BASE = "appMQTnobUNRytIp7"
PAPERS = "tblnnMrZszhboU4uD"

FETCH = ["Name", "Notes", "Attachments", "Feedback", "Primary Author", "Animal",
         "Year", "Models referenced", "Is the model paper", "Review Papers", "DOI",
         "Secondary Authors", "Afferent Types", "Animal Study Potential",
         "Robot/Sim Translation", "Prune Status"]


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
        q = "?" + urllib.parse.urlencode(params, doseq=True)
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


papers_raw, offset = [], None
while True:
    p = [("pageSize", 100), ("fields[]", FETCH)]
    if offset:
        p.append(("offset", offset))
    d = air(f"{BASE}/{PAPERS}", p)
    papers_raw.extend(d.get("records", []))
    offset = d.get("offset")
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
        "afferents": sel(f.get("Afferent Types", [])) if isinstance(f.get("Afferent Types"), list) else [],
        "feedback": [fb_names.get(r, r) for r in f.get("Feedback", [])],
        "models_ref": [md_names.get(r, r) for r in f.get("Models referenced", [])],
        "models2": [md_names.get(r, r) for r in f.get("Is the model paper", [])],
        "reviews": [rv_names.get(r, r) for r in f.get("Review Papers", [])],
        "has_pdf": bool(f.get("Attachments")),
        "has_notes": bool(notes.strip()),
        "notes": notes.strip()[:1500],
        "animal_study": (f.get("Animal Study Potential", "") or "").strip()[:1000],
        "robot_sim": (f.get("Robot/Sim Translation", "") or "").strip()[:1000],
        "prune": sel(f.get("Prune Status", "")) or "",
    })

with open(os.path.join(OUT, "sadb_export.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False)

cols = ["id", "title", "primary", "secondary", "year", "doi", "animals", "afferents",
        "feedback", "models_ref", "models2", "reviews", "has_pdf", "has_notes", "notes",
        "animal_study", "robot_sim", "prune"]
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
n_af = sum(1 for r in out if r["afferents"])
print(f"wrote {len(out)} records | notes {n_notes} | pdf {n_pdf} | animals {n_an} | "
      f"feedback {n_fb} | afferents {n_af} | review-links {n_rev} | bare {len(out)-n_notes}")

# --- archive merge (Ben, 2026-10-02: "we'll keep them in the app") ---
# Records deleted from Airtable (over the free-plan record cap) are re-added
# from their full-row snapshots in archive/deleted_records_*.json, flagged
# archived=True, so the HTML app keeps their rules/insight forever.
import glob as _glob
live_ids = {r["id"] for r in out}
for path in sorted(_glob.glob(os.path.join(OUT, "..", "archive", "deleted_records_*.json"))):
    for snap in json.load(open(path, encoding="utf-8")):
        if snap["id"] not in live_ids:
            snap["archived"] = True
            out.append(snap)
n_arch = sum(1 for r in out if r.get("archived"))
with open(os.path.join(OUT, "sadb_export.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False)
with open(os.path.join(OUT, "sadb_export.csv"), "w", encoding="utf-8-sig", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols + ["archived"])
    w.writeheader()
    for r in out:
        w.writerow({c: (("; ".join(r[c]) if isinstance(r[c], list) else r[c])) for c in cols}
                   | {"archived": bool(r.get("archived"))})
print(f"after archive merge: {len(out)} records ({n_arch} app-only / removed from Airtable)")
