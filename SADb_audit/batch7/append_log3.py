"""Append the VOSviewer-retirement WORKFLOW row to curation_log.csv."""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(os.path.dirname(HERE), "curation_log.csv")

rows = [
    ["2026-09-22", "VIZ", "WORKFLOW", "ben",
     "VOSviewer RETIRED as a dead end (Ben ruling): generic viewer, no drill-through, "
     "mega-cited bubbles dominate. vosviewer/ folder (map/network/README) DELETED; "
     "vos_build2.py replaced by build_citation_graph.py (app-only outputs: "
     "export/sadb_cites.json directed adjacency + export/sadb_layout.json "
     "layout/clusters; OpenAlex cache unchanged). App upgraded: bubble map gained a "
     "Layout selector (Year x citations | Topic landscape = the citation-network "
     "layout) and a Neuron style - papers drawn as neurons (soma/dendrites/axon, "
     "deterministic per id); focusing a paper draws its synapses per Ben's convention: "
     "excitatory open triangles FROM citing papers, inhibitory filled circles ONTO "
     "cited papers; spotlight now refits to the neighborhood; verified in-browser. "
     "Roadmap for online mode + Web of Science/Google Scholar enrichment + a "
     "Zotero-like PDF-capture extension + multi-lab deployment (Quinn / Chiel / "
     "Bueschges / Szczecinski / Webster-Wood) written in app/README.md. Airtable "
     "remains the primary tool",
     "", "", "", ""],
]

with open(LOG, "a", encoding="utf-8", newline="") as fh:
    csv.writer(fh).writerows(rows)
print("logged")
