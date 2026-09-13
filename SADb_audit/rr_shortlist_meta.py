#!/usr/bin/env python
# Fetch Crossref metadata for the RR-digest shortlist -> rr_digest_shortlist.json
import json, subprocess, os, sys

MAILTO = "bbolen@pdx.edu"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "batch3")

SHORTLIST = [
    ("10.7554/elife.98841", "Rybak lab: operation regimes of spinal circuits (eLife final version) - RR recommended"),
    ("10.7554/elife.103504", "Operation of spinal sensorimotor circuits controlling phase transitions (eLife) - RR family"),
    ("10.1152/jn.00104.2024", "Forelimb movements contribute to hindlimb cutaneous reflexes - cutaneous locomotor reflexes"),
    ("10.1113/jp286151", "Changes in intra/interlimb reflexes from hindlimb cutaneous afferents (J Physiol) - load/cutaneous gating"),
    ("10.1113/jp286808", "Sister paper: forelimb cutaneous reflex changes (J Physiol)"),
    ("10.1016/j.neunet.2024.106422", "Spinal circuit model, asymmetric cervical-lumbar layout - model paper"),
    ("10.1098/rsos.240207", "Sensory feedback and central neuronal interactions - sensory+CPG integration"),
    ("10.1152/jn.00248.2023", "Dynamic spinal reflex adaptation during locomotor adaptation - reflex plasticity"),
    ("10.3389/fncir.2023.1235181", "Lumbar V3 interneurons direct excitatory input - extends Zhang V3 line"),
    ("10.1016/j.cub.2023.07.014", "Distinct roles of spinal commissural interneurons - left-right circuitry (Rybak 2013 experimental kin)"),
    ("10.1371/journal.pcbi.1012101", "Balancing central control and sensory feedback - CPG+feedback model"),
    ("10.1371/journal.pcbi.1013494", "Physiologically inspired hybrid CPG/reflex controller - robot control, directly lab-relevant"),
    ("10.1016/j.cub.2025.09.030", "A spinal circuit for skilled locomotion"),
    ("10.1016/j.expneurol.2023.114496", "Spinal control of locomotion before and after SCI"),
    ("10.1523/jneurosci.2015-22.2023", "Excitatory and inhibitory descending commissural interneurons"),
    ("10.1007/s00422-023-00970-z", "BCM rule lets spinal cord model learn rhythmic - learning CPG model"),
    ("10.1038/s42003-024-06843-w", "Interlimb coordination is not strictly controlled during walking"),
    ("10.1152/jn.00331.2025", "Speed-dependent locomotor adjustments following staggered hemisection (Frigon family)"),
    ("10.1101/2025.11.11.687930", "Adaptive interlimb coordination to sudden ground loss, neuromusculoskeletal cat CPG model - ALREADY IN BEN'S ZOTERO (C5LYIBNK), needs Airtable only"),
]

rows = []
for doi, why in SHORTLIST:
    r = subprocess.run(["curl.exe", "-s", "-m", "30", "-A", "mailto:" + MAILTO,
                        "https://api.crossref.org/works/" + doi], capture_output=True)
    try:
        m = json.loads(r.stdout.decode("utf-8", "replace"))["message"]
    except Exception:
        rows.append(dict(doi=doi, why=why, fetch="FAILED"))
        continue
    auth = m.get("author", [])
    year = None
    for k in ("published-print", "published-online", "issued", "created"):
        if m.get(k, {}).get("date-parts"):
            year = m[k]["date-parts"][0][0]
            if k != "created":
                break
    rows.append(dict(
        doi=doi, why=why,
        title=(m.get("title") or [""])[0],
        container=(m.get("container-title") or [""])[0],
        year=year,
        authors="; ".join("%s %s" % (a.get("family", ""), a.get("given", ""))[:40] for a in auth[:4]) if False else "; ".join(((a.get("family") or "") + " " + (a.get("given") or "")).strip() for a in auth[:4]),
        first_author=(auth[0].get("family") if auth else ""),
    ))
    print("%-28s %4s %-52s %s" % (doi, year, ((m.get("title") or [""])[0])[:52], rows[-1]["first_author"]))

with open(os.path.join(OUT, "rr_digest_shortlist.json"), "w", encoding="utf-8") as f:
    json.dump(rows, f, indent=1, ensure_ascii=False)
print("\nwrote rr_digest_shortlist.json with", len(rows), "rows")
