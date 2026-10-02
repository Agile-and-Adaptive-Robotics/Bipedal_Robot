"""Build the audited deletion manifest for "the 90" (Ben's go, 2026-10-01).

Groups (all from local files — zero API):
  A  general NN/ML theory (clusters 9,10,12)      B  aerodynamics (cluster 11)
  C1 flow/flight sensing, no neuroscience signal  D  bio-inspired mechanisms (14,15)
  E  fish swimming biomechanics (8, borderline — included per "the 90")
  P  prior strict prune-candidate flags (export prune field)
  +  satellite orphans (Models 3 / Review 1 / Feedback 3, from prune_shortlist.csv)

Safety checks before writing the manifest:
  - dedupe across groups
  - list any delete-set paper carrying Models/Review links (links will empty)
  - snapshot every deleted paper's full export row for the audit trail
Output: deletion_manifest.json + deleted_records_20261001.json (full snapshots).
"""
import csv, json, os, re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
records = json.load(open(os.path.join(HERE, "export", "sadb_export.json"), encoding="utf-8"))
layout = json.load(open(os.path.join(HERE, "export", "sadb_layout.json"), encoding="utf-8"))

NEURO = re.compile(r"neur|afferent|reflex|sensorimotor|gangli|cpg|central pattern|"
                   r"motor control|electrophysiolog|spike|synap|muscle spindle|"
                   r"propriorecept|propriocept|sensill|eme", re.I)

byid = {r["id"]: r for r in records}
groups = defaultdict(list)

for r in records:
    cl = layout.get(r["id"], {}).get("cl", -1)
    if cl in (9, 10, 12):
        groups["A_nn_ml_theory"].append(r["id"])
    elif cl == 11:
        groups["B_aerodynamics"].append(r["id"])
    elif cl == 7 and not NEURO.search((r["title"] or "") + " " + (r.get("notes") or "")):
        groups["C1_flow_sensing_no_neuro"].append(r["id"])
    elif cl in (14, 15):
        groups["D_bioinspired_mechanisms"].append(r["id"])
    elif cl == 8:
        groups["E_fish_biomechanics"].append(r["id"])
    if r.get("prune") == "prune-candidate":
        groups["P_prior_flags"].append(r["id"])

sat = defaultdict(list)
with open(os.path.join(HERE, "prune_shortlist.csv"), encoding="utf-8-sig") as fh:
    for row in csv.DictReader(fh):
        if row["table"] == "Papers":
            continue  # papers already covered above via P group? No — CSV rows were
            # the ORIGINAL strict flags; the export prune field is the live truth.
        sat[row["table"]].append(row["record_id"])

papers = sorted({rid for g in groups.values() for rid in g})
overlap_note = []
for rid in papers:
    r = byid.get(rid, {})
    if r.get("models2") or r.get("reviews"):
        overlap_note.append((rid, r.get("models2"), r.get("reviews")))

manifest = {"papers": papers,
            "paper_group_counts": {k: len(v) for k, v in sorted(groups.items())},
            "papers_total": len(papers),
            "satellite": {"Models": sat.get("Models", []),
                          "Review Papers": sat.get("Review Papers", []),
                          "Feedback": sat.get("Feedback", [])},
            "links_that_will_empty": overlap_note}

json.dump(manifest, open(os.path.join(HERE, "deletion_manifest.json"), "w"), indent=1)
snapshots = [byid.get(rid) for rid in papers if rid in byid]
json.dump(snapshots, open(os.path.join(HERE, "deleted_records_20261001.json"), "w"),
          ensure_ascii=False, indent=1)

print("groups:", manifest["paper_group_counts"])
print("papers to delete:", len(papers))
print("satellite orphans:", {k: len(v) for k, v in manifest["satellite"].items()})
print("papers with model/review links that will empty:", len(overlap_note))
for rid, m2, rv in overlap_note[:10]:
    print("   ", rid, byid[rid]["title"][:60], "| models2:", m2, "| reviews:", rv)
