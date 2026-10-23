"""Export the full Airtable Papers corpus to export/sadb_export.json + .csv.

Reads the PAT from D:\Github\api_credentials_local.txt (never stores it).
Stdlib only. Run: python export_corpus.py
Regenerate after any curation batch — the HTML app + pivot consume these files.

2026-09-28 schema v2: the duplicate "Models copy" fields were renamed
(fldBLowhKJcjuFMSK is now "Review Papers"; the empty twin is "(unused) …"),
which UNIQUE-IFIES all field names — a fields[] filter no longer 422s, so this
fetch is now projection-filtered. New curation-layer fields exported:
Afferent Types, Animal Study Potential, Robot/Sim Translation, Prune Status.

2026-10-23 FULL PRE-PRUNE BASE: the export also carries the non-paper records
— Models, Review Papers, Feedback pathways, Recorders — tagged "recordType".
Retained records come from the Airtable working base via the SAME access path
and credential source as the papers. The 54 deleted orphan non-paper records
(Models 3, Review Papers 47+1, Feedback 3, Recorders 0) are reconstructed from
the local deletion snapshots (deletion_manifest.json satellite ids, the
prune_shortlist.csv names, review_orphans_20261002.json) and flagged
archived=True with the same convention as the archived papers. Everything is
assembled in memory and written ONCE at the end: a failed Airtable pull never
clobbers the previous export.
"""
import csv, json, os, re, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "export")
os.makedirs(OUT, exist_ok=True)
BASE = "appMQTnobUNRytIp7"
PAPERS = "tblnnMrZszhboU4uD"
MODELS = "tblsBq9IEv7dZe6fn"
REVIEWS = "tblSEubKcRId4wYMK"
FEEDBACK = "tblot5mo4s5KgN5le"
RECORDERS = "tbloZUIALmRULPBFY"

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


def table_records(table_id, label):
    out, offset = [], None
    while True:
        p = [("pageSize", 100)]
        if offset:
            p.append(("offset", offset))
        d = air(f"{BASE}/{table_id}", p)
        out.extend(d.get("records", []))
        offset = d.get("offset")
        if not offset:
            break
    print(f"{label}: {len(out)} records")
    return out


def names_for(records):
    return {rec["id"]: rec["fields"].get("Name", "") for rec in records}


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

fb_recs = table_records(FEEDBACK, "Feedback")
md_recs = table_records(MODELS, "Models")
rv_recs = table_records(REVIEWS, "Review Papers")
rc_recs = table_records(RECORDERS, "Recorders")
fb_names = names_for(fb_recs)
md_names = names_for(md_recs)
rv_names = names_for(rv_recs)


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
        "recordType": "paper",
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

n_notes = sum(1 for r in out if r["has_notes"])
n_pdf = sum(1 for r in out if r["has_pdf"])
n_rev = sum(1 for r in out if r["reviews"])
n_an = sum(1 for r in out if r["animals"])
n_fb = sum(1 for r in out if r["feedback"])
n_af = sum(1 for r in out if r["afferents"])
print(f"live papers {len(out)} | notes {n_notes} | pdf {n_pdf} | animals {n_an} | "
      f"feedback {n_fb} | afferents {n_af} | review-links {n_rev} | bare {len(out)-n_notes}")

# --- retained non-paper records (Review Papers / Models / Feedback / Recorders) ---
def nonpaper(rec, rt):
    f = rec.get("fields", {})
    notes = f.get("Notes", "") or ""
    return {
        "id": rec["id"],
        "recordType": rt,
        "title": f.get("Name", "") or "",
        "primary": "", "secondary": [], "year": "", "doi": "",
        "animals": [], "afferents": [], "feedback": [], "models_ref": [],
        "models2": [], "reviews": [],
        "has_pdf": bool(f.get("Attachments")),
        "has_notes": bool(notes.strip()),
        "notes": notes.strip()[:1500],
        "animal_study": "", "robot_sim": "", "prune": "",
    }


non_paper = ([nonpaper(r, "model") for r in md_recs]
             + [nonpaper(r, "review") for r in rv_recs]
             + [nonpaper(r, "feedback") for r in fb_recs]
             + [nonpaper(r, "recorder") for r in rc_recs])
print(f"retained non-paper records: {len(non_paper)}")

# --- archive merge (Ben, 2026-10-02: "we'll keep them in the app") ---
# Records deleted from Airtable (over the free-plan record cap) are re-added
# from their full-row snapshots in archive/deleted_records_*.json, flagged
# archived=True, so the HTML app keeps their rules/insight forever.
import glob as _glob
live_ids = {r["id"] for r in out}
for path in sorted(_glob.glob(os.path.join(OUT, "..", "archive", "deleted_records_*.json"))):
    for snap in json.load(open(path, encoding="utf-8")):
        if snap["id"] not in live_ids:
            snap["recordType"] = "paper"
            snap["archived"] = True
            out.append(snap)
n_arch = sum(1 for r in out if r.get("archived"))
print(f"papers after archive merge: {len(out)} records ({n_arch} app-only / removed from Airtable)")

# --- archived non-paper orphans, reconstructed from the deletion snapshots ---
# Sources of record: deletion_manifest.json satellite ids (what was deleted
# 2026-10-01), prune_shortlist.csv (their Names), review_orphans_20261002.json
# (the 47 Review Papers records orphaned+deleted 2026-10-02, id->Name).
RT_BY_TABLE = {"Models": "model", "Review Papers": "review", "Feedback": "feedback"}
orphans = {}
with open(os.path.join(HERE, "prune_shortlist.csv"), encoding="utf-8-sig") as fh:
    for row in csv.DictReader(fh):
        if row["table"] in RT_BY_TABLE:
            orphans[row["record_id"]] = {"id": row["record_id"],
                                         "recordType": RT_BY_TABLE[row["table"]],
                                         "title": (row["name/title"] or "").strip()}
rev_orphans = json.load(open(os.path.join(HERE, "review_orphans_20261002.json"),
                             encoding="utf-8"))
for rid, name in rev_orphans.items():
    assert rid not in orphans, f"review orphan {rid} already snapshotted"
    orphans[rid] = {"id": rid, "recordType": "review", "title": name}
manifest = json.load(open(os.path.join(HERE, "deletion_manifest.json"), encoding="utf-8"))
sat_ids = {}
for tbl, rt in RT_BY_TABLE.items():
    for rid in manifest["satellite"][tbl]:
        sat_ids[rid] = rt
missing = {rid for rid, rt in sat_ids.items()
           if rid not in orphans or orphans[rid]["recordType"] != rt}
extra = set(orphans) - set(sat_ids) - set(rev_orphans)
assert not missing, f"orphan snapshots missing manifest satellite ids: {sorted(missing)}"
assert not extra, f"orphan ids beyond manifest satellites + review orphans: {sorted(extra)}"

live_all = live_ids | {r["id"] for r in non_paper}
archived_non = []
for o in orphans.values():
    if o["id"] in live_all:
        continue  # a live Airtable copy would win; never duplicate
    o.update({"primary": "", "secondary": [], "year": "", "doi": "", "animals": [],
              "afferents": [], "feedback": [], "models_ref": [], "models2": [],
              "reviews": [], "has_pdf": False,
              "has_notes": bool(o["title"]), "notes": "",
              "animal_study": "", "robot_sim": "", "prune": "",
              "archived": True})
    archived_non.append(o)
orph_by_type = {}
for o in archived_non:
    orph_by_type[o["recordType"]] = orph_by_type.get(o["recordType"], 0) + 1
print(f"archived non-paper orphans: {len(archived_non)} {orph_by_type}")

out += non_paper + archived_non

# --- single write: assembled fully in memory first, so a failed Airtable pull
# above never clobbers the previous export ---
with open(os.path.join(OUT, "sadb_export.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False)

cols = ["id", "title", "primary", "secondary", "year", "doi", "animals", "afferents",
        "feedback", "models_ref", "models2", "reviews", "has_pdf", "has_notes", "notes",
        "animal_study", "robot_sim", "prune"]
with open(os.path.join(OUT, "sadb_export.csv"), "w", encoding="utf-8-sig", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols + ["record_type", "archived"])
    w.writeheader()
    for r in out:
        w.writerow({c: (("; ".join(r[c]) if isinstance(r[c], list) else r[c])) for c in cols}
                   | {"record_type": r.get("recordType", "paper"),
                      "archived": bool(r.get("archived"))})

by_type, by_type_arch = {}, {}
for r in out:
    by_type[r.get("recordType", "paper")] = by_type.get(r.get("recordType", "paper"), 0) + 1
    if r.get("archived"):
        by_type_arch[r.get("recordType", "paper")] = by_type_arch.get(r.get("recordType", "paper"), 0) + 1
n_arch = sum(1 for r in out if r.get("archived"))
print(f"wrote {len(out)} records by type: {by_type}")
print(f"archived {n_arch} by type: {by_type_arch} | live {len(out) - n_arch}")
