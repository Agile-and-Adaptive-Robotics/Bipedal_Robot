"""Append batch-7 rows + WORKFLOW milestone to curation_log.csv."""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(os.path.dirname(HERE), "curation_log.csv")
DATE = "2026-09-20"

rows = [
    [DATE, "7", "L8WLNNF3", "recDNNmx4fL6zAhv1", "Yakovenko and Prochazka 2001 chapter", "no-text", "", "",
     "NO TEXT FOUND (Zotero local/web, Europe PMC, Crossref, OpenAlex) — chapter in 'Biomechanics and Neural Control of Posture and Movement'. Left bare. BEN: library pass"],
    [DATE, "7", "6V4RWQNW", "recwIiaV3T7VT6UtQ", "Nourse 2019 Drosophila inter-leg CPG", "yes", "Insects; Arthropods", "",
     "Grounded: Zotero local abstract. FLAG: abstract says 'intra-leg load feedback' without naming campaniform sensilla — no Feedback link until full text confirms (Ben's call)"],
    [DATE, "7", "I7WF27CS", "reczvpAn0Yl9BlIQZ", "Ziskind-Conhaim and Hochman 2017 review", "yes", "Mice; Mammals", "",
     "Grounded: Zotero local abstract. Review Papers twin 'Ziskind-Conhaim and Hochman 2017' = reckp7LCl56L0mh0a"],
    [DATE, "7", "KMATJMI4", "recZNGVV5CmEzj8R4", "Andersson and Grillner 1983 hip entrainment", "yes", "Cat; Mammals", "",
     "Grounded: Zotero local abstract. FLAG: no Feedback link — the paper demonstrates hip afferent entrainment but does NOT identify the afferent class, so 'Ia or II stance to swing' would overclaim. Ben may want a dedicated record 'hip afferent entrainment of CPG'"],
    [DATE, "7", "GGHBEBW4", "recHzZPAbILVjbMvT", "Quinlan and Kiehn 2007 commissural INs", "yes", "Mice; Mammals", "",
     "Grounded: Zotero local abstract"],
    [DATE, "7", "WJ4ZTW5P", "recbEFSrkhmkQiXhA", "Lewinger 2006 SCASM robots", "no-text", "", "",
     "NO TEXT FOUND — no DOI; OpenAlex has the work (cited_by 3) but no abstract. Robot paper. BEN: library pass if wanted"],
    [DATE, "7", "XIJH8N3A", "reci1mqLtEvfwGu9Z", "Bouyer and Rossignol 2003 II spinal cats", "yes", "Cat; Mammals", "Cutaneous stance modification",
     "Grounded: Zotero local abstract. YEAR FILLED 2003 (was blank; per DOI/journal). Note: part II of a companion pair"],
    [DATE, "7", "3URZQV5H", "recXVuB69yaIn6437", "Duysens and Forner-Cordero 2018 review", "yes", "Human", "",
     "Grounded: Zotero local abstract. Review Papers twin 'Duysens and Forner-Cordero 2018' = recvisSCDWFyYXAfJ. FLAG: near-twin name vs existing 'Duysens and Forner-Cordero 2019' (recbuZPxtSFoexUn3's review, a different paper) — Ben may want to unify or differentiate"],
    [DATE, "7", "8TLXCQR4", "recH9rSuUtGykbuWE", "Bondy 2016 multifunctional CPG", "yes", "Cat; Mammals", "",
     "Grounded: Zotero local abstract. Models twin 'Bondy 2016' = recVUX3okwyChn7ew. FLAG: paw-skin afferent regime switching not linked to a Feedback record (cutaneous flexor excitation would overclaim)"],
    [DATE, "7", "BU4APFEI", "recSEg1aVNs16JtjP", "Kjaerulff and Kiehn 1996 lesion study", "yes", "Rat; Mammals", "",
     "Grounded: Zotero local abstract"],
    [DATE, "7", "WORKFLOW", "verify", "Batch 7 COMPLETE + verified: queue rows 52-61 (8 curated, 2 no-text: Yakovenko 2001 chapter, Lewinger 2006). 3 twin records (Reviews: Ziskind-Conhaim and Hochman 2017 reckp7LCl56L0mh0a, Duysens and Forner-Cordero 2018 recvisSCDWFyYXAfJ; Models: Bondy 2016 recVUX3okwyChn7ew), paper-side auto-links verified live. Year fix: Bouyer/Rossignol Part II -> 2003 (was blank). Re-pull of all 10 ids: every field echo matched. RUNNING TOTAL: 66/383 queue papers curated (batches 1-7); 4 no-text total awaiting Ben library pass. NEXT = batch 8 = rows 62-70 (grounding already in batch6/ground_*.txt); audit subagent due after batch 10", "", "", "", ""],
]

with open(LOG, "a", encoding="utf-8", newline="") as fh:
    csv.writer(fh).writerows(rows)
print(f"appended {len(rows)} rows to curation_log.csv")
