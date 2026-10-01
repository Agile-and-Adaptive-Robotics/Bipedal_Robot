"""Goal 1: compute the prune shortlist for the over-cap Airtable base (2026-09-28).

Base census: Papers 943 + Models 172 + Review Papers 75 + Feedback 38 +
Recorders 4 = 1232 records (free-tier cap = 1000; need ~232 gone to fit).

Scoring a Papers record as prune-candidate requires ALL of:
  - no Notes, no Feedback link, no Models 2 link, no Review Papers link
  - isolated in the corpus citation network (no cluster, 0 in/degree)
  - OpenAlex cited_by_count <= CITE_MAX (default 2)
Candidates get Prune Status = "prune-candidate" + a Prune Reason written to
Airtable (flags only — NEVER deletes; Ben approves every deletion).

Also audits the satellite tables for orphan records (Models with no Paper
link, Review Papers with no Paper link, Feedback with no Papers link) into
prune_shortlist.csv — those tables are the cheapest place to recover records.

Stdlib only. Run: myo python prune_shortlist.py
"""
import csv, json, os, re, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "appMQTnobUNRytIp7"
PAPERS = "tblnnMrZszhboU4uD"
FLD_STATUS = "fld9jhnEQW738POne"    # Prune Status (singleSelect)
FLD_REASON = "fldcW2N1nn2pNxfv0"    # Prune Reason (singleLineText)
CITE_MAX = 2


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


def air_patch(body):
    req = urllib.request.Request(
        f"https://api.airtable.com/v0/{BASE}/{PAPERS}",
        data=json.dumps(body).encode(),
        headers={"Authorization": "Bearer " + _PAT, "Content-Type": "application/json"},
        method="PATCH")
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


records = json.load(open(os.path.join(HERE, "export", "sadb_export.json"), encoding="utf-8"))
layout = json.load(open(os.path.join(HERE, "export", "sadb_layout.json"), encoding="utf-8"))
cites = json.load(open(os.path.join(HERE, "export", "sadb_cites.json"), encoding="utf-8"))
cited_by = {}
for ra, rbs in cites.items():
    for rb in rbs:
        cited_by.setdefault(rb, []).append(ra)

cands, rows = [], []
for r in records:
    rid = r["id"]
    L = layout.get(rid, {})
    indeg, outdeg = len(cited_by.get(rid, [])), len(cites.get(rid, []))
    oa_c = L.get("c", 0)
    rich = (r["has_notes"] or r["feedback"] or r["models2"] or r["reviews"])
    connected = (L.get("cl", -1) >= 0 or indeg > 0 or outdeg > 0)
    if (not rich) and (not connected) and oa_c <= CITE_MAX:
        reason = (f"no notes/pathway/model/review links; isolated in corpus network; "
                  f"OpenAlex cites {oa_c}; {'PDF' if r['has_pdf'] else 'no PDF'}")
        cands.append((rid, reason))
    rows.append({
        "record_id": rid, "title": r["title"][:100], "year": r["year"],
        "doi": r["doi"], "has_notes": int(r["has_notes"]),
        "feedback_links": len(r["feedback"]), "oa_cited": oa_c,
        "cluster": L.get("cl", -1), "in_deg": indeg, "out_deg": outdeg,
        "has_pdf": int(r["has_pdf"]), "verdict": "",
    })

# --- satellite-table orphan audit (cheap record recovery) ---
def fetch_table(tid, label, fields):
    out, offset = [], None
    while True:
        p = [("pageSize", 100), ("fields[]", fields)]
        if offset:
            p.append(("offset", offset))
        d = air(f"{BASE}/{tid}", p)
        out.extend(d.get("records", []))
        offset = d.get("offset")
        if not offset:
            break
    print(f"{label}: {len(out)} records")
    return out


orphans = []
models = fetch_table("tblsBq9IEv7dZe6fn", "Models", ["Name", "Paper"])
for m in models:
    if not m["fields"].get("Paper"):
        orphans.append(("Models", m["id"], m["fields"].get("Name", ""), "no Paper link"))
reviews = fetch_table("tblSEubKcRId4wYMK", "Review Papers", ["Name", "Paper"])
for m in reviews:
    if not m["fields"].get("Paper"):
        orphans.append(("Review Papers", m["id"], m["fields"].get("Name", ""), "no Paper link"))
fbs = fetch_table("tblot5mo4s5KgN5le", "Feedback", ["Name", "Papers"])
for m in fbs:
    if not m["fields"].get("Papers"):
        orphans.append(("Feedback", m["id"], m["fields"].get("Name", ""), "no Papers link"))

out_csv = os.path.join(HERE, "prune_shortlist.csv")
with open(out_csv, "w", encoding="utf-8-sig", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["table", "record_id", "name/title", "year", "reason/verdict"])
    for rid, reason in cands:
        r = next(x for x in records if x["id"] == rid)
        w.writerow(["Papers", rid, r["title"][:100], r["year"], reason])
    for t, rid, name, why in orphans:
        w.writerow([t, rid, name, "", why])
n_paper_c, n_orph = len(cands), len(orphans)
print(f"\nPapers prune-candidates: {n_paper_c}")
print(f"satellite orphans: {n_orph} "
      f"(Models {sum(1 for o in orphans if o[0]=='Models')}, "
      f"Review {sum(1 for o in orphans if o[0]=='Review Papers')}, "
      f"Feedback {sum(1 for o in orphans if o[0]=='Feedback')})")
print(f"removing ALL candidates+orphans would free {n_paper_c + n_orph} of the "
      f"232 needed to reach the 1000-record cap -> {1232 - n_paper_c - n_orph} remaining")

# --- write flags to Airtable (Papers only; satellites just listed in CSV) ---
if cands:
    items = list(cands)
    done = 0
    for i in range(0, len(items), 10):
        chunk = items[i:i + 10]
        body = {"records": [{"id": rid,
                             "fields": {FLD_STATUS: "prune-candidate",
                                        FLD_REASON: reason[:200]}}
                            for rid, reason in chunk], "typecast": False}
        for attempt in range(3):
            try:
                air_patch(body)
                done += len(chunk)
                break
            except Exception as e:
                print(f"  batch {i//10} retry {attempt}: {e}")
                time.sleep(3 + 3 * attempt)
        else:
            print(f"  FAILED batch at {i}")
        time.sleep(0.35)
    print(f"flagged {done} Papers records as prune-candidate in Airtable")
print("wrote", out_csv)
