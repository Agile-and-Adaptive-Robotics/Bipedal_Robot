"""Triage staged PDFs by folder: which are NOT in the Airtable corpus yet,
grouped by Zotero folder, split into in-scope (sensory afferent / motor control)
vs out-of-scope. Outputs author_fix/aarl_gap_by_folder.csv (in-scope, review list).
"""
import csv, os, json

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

# corpus DOI set (all Airtable papers)
corpus = set()
for f, col in ((os.path.join(SAD, "airtable_papers_slim.csv"), "DOI"),
               (os.path.join(SAD, "airtable_rest_import_clean.csv"), "DOI")):
    with open(f, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            d = norm_doi(r.get(col))
            if d:
                corpus.add(d)
with open(os.path.join(SAD, "batch3", "rr_digest_shortlist.json"), encoding="utf-8") as f:
    for s in json.load(f):
        corpus.add(norm_doi(s["doi"]))
corpus.add(norm_doi("10.7554/elife.107480"))
print("corpus DOIs:", len(corpus))

# in-scope folders (Ben's criterion: sensory afferent feedback + motor control)
IN_SCOPE = {
    "Sensory Feedback", "Proprioceptive Feedback", "Afferent", "Ia Afferent Cite",
    "Local Feedback", "CPGs", "CPG PDC", "CPG", "Dual-Layered CPGs",
    "Neural Oscillators", "Neuron Modeling", "Models", "Sensors",
    "Biped", "Quadruped", "Biomechanics", "Muscles", "Synergies",
    "Animal Data", "Swimming", "Balance Control", "Balance",
    "Bipedal balance control", "Stability in Legged Locomotion",
    "Damping in Locomotion", "Spikes to Muscle Activation",
    "Spiking/NonSpiking Networks", "Functional Subnetwork Approach",
    "Synchronization", "Mammals", "Insects", "Fish", "Crustacean",
    "Mollusk", "Amphibians", "Respiration", "Behavioral Examples",
    "Plasticity", "Foundational",
}

# staged + foldered rows
rows = []
with open(os.path.join(SAD, "staged_with_folders.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        rows.append(r)

seen_doi = set()
gap = {}      # (folder) -> list of unique papers
in_air = 0
for r in rows:
    d = norm_doi(r["doi"])
    if not d or d in corpus:
        in_air += 1
        continue
    if d in seen_doi:
        continue
    seen_doi.add(d)
    folders = [x for x in (r["folders"] or "").split(";") if x]
    scope = sorted({f for f in folders if f in IN_SCOPE})
    others = sorted({f for f in folders if f not in IN_SCOPE})
    bucket = ";".join(scope) if scope else "OUT-OF-SCOPE (" + ";".join(others[:3]) + ")"
    gap.setdefault(bucket, []).append({
        "doi": d, "surname": r["surname"], "year": r["year"],
        "title": r["title"], "folders": r["folders"],
        "att_key": r["att_key"], "library": r["library"],
    })

print("staged rows:", len(rows), "| already in Airtable (rows):", in_air)
print("unique papers NOT in Airtable:", len(seen_doi))
print()
print("=== in-scope folders (gap counts, unique papers) ===")
tot = 0
for b in sorted(gap):
    if b.startswith("OUT-OF-SCOPE"):
        continue
    print("%6d  %s" % (len(gap[b]), b))
    tot += len(gap[b])
print("in-scope gap total:", tot)
out_tot = sum(len(v) for k, v in gap.items() if k.startswith("OUT-OF-SCOPE"))
print("out-of-scope total:", out_tot)

with open(os.path.join(HERE, "aarl_gap_by_folder.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["folder_bucket", "doi", "surname", "year", "title", "folders", "att_key", "library"])
    w.writeheader()
    for b, items in sorted(gap.items()):
        for it in items:
            it2 = dict(it)
            it2["folder_bucket"] = b
            w.writerow(it2)
print("wrote aarl_gap_by_folder.csv")
