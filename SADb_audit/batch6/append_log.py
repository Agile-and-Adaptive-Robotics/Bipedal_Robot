"""Append batch-6 rows + WORKFLOW milestone to curation_log.csv (UTF-8, quoted)."""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(os.path.dirname(HERE), "curation_log.csv")
DATE = "2026-09-18"

rows = [
    [DATE, "6", "8QVW6ZRY", "recBMnoy8KuuLLSdG", "Fujiki 2018 split-belt phase reset", "yes", "Rat; Mammals", "",
     "Grounded: Zotero local abstract. Models twin 'Fujiki 2018' = recxgsj4gWFoEA3Jz (auto-linked Models 2). No Feedback link: abstract pins no specific pathway (flag: full text would detail the hip-flexor afferent class)"],
    [DATE, "6", "D8PKQH5X", "recLkJ2b8tJBQhhQO", "Thiry 2020 mouse genetics review chapter", "no-text", "", "",
     "NO TEXT FOUND anywhere (Zotero local/web abstract, Europe PMC, Crossref, OpenAlex) — Elsevier chapter in 'The Neural Control of Movement'. Left bare per no-grounding-no-note rule. BEN: library pass needed"],
    [DATE, "6", "JV5ZJCH9", "recMIXEbZvmgbZrNx", "Cote 2018 spinal interneurons review", "yes", "Cat; Mice; Rat; Human; Mammals", "",
     "Grounded: Zotero local abstract. Review Papers twin 'Cote 2018' = rectKTUSAaUJIw5tG (auto-linked via Models copy field)"],
    [DATE, "6", "L2RUU8KR", "recYLE3O2RRPv3yTp", "Boulland 2012 MPIO cell tracking", "yes", "", "",
     "Grounded: Zotero local abstract (truncated at in-vivo results — note covers in-vitro findings only). Animals blank: cell-line methods paper, no whole-animal preparation in abstract. Peripheral to locomotor core (kept: Perreault co-author)"],
    [DATE, "6", "XT86CDBR", "recvuLXfmP2ns7vAO", "Nichols and Ross 2009 lambda-model force feedback", "yes", "Cat; Mammals", "Ib excitatory; Ib inhibition",
     "Grounded: Zotero local abstract"],
    [DATE, "6", "X24SE498", "rec2vnLEhwtLZDpnJ", "Ekeberg 2004 stick insect walking sim", "yes", "Stick Insect; Insects; Arthropods", "Chordotonal organ multi-synaptic excitation; Chordotonal organ multi-synaptic inhibition",
     "Grounded: Europe PMC abstract. Models twin 'Ekeberg 2004' = receLzSMQA9RYo4t3"],
    [DATE, "6", "BKERX7AM", "reccbV3L4NG7Kq8ZQ", "Prochazka 2017 neuromuscular models chapter", "no-text", "", "",
     "NO TEXT FOUND anywhere — Elsevier chapter in 'Bioinspired Legged Locomotion' (co-authors Gosgnach, Capaday, Geyer). Left bare. BEN: library pass needed"],
    [DATE, "6", "U4HBKUSE", "recXdsqbLWdkE9OL0", "Cazalets 1992 5-HT/EAA CPG activation", "yes", "Rat; Mammals", "Fictive locomotion without sensory feedback",
     "Grounded: Zotero local abstract. Fictive-without-feedback link: isolated-spinal pharmacological preparation, same usage as Grillner & Zangger 1975"],
    [DATE, "6", "VMKUQ2IQ", "reci22Nke6ZTYEXCJ", "Hatz 2012 ankle extensor feedback gains", "yes", "Cat; Mammals", "Ib excitatory; type II excitatory",
     "Grounded: Zotero local abstract + Europe PMC FULL abstract (Zotero's was truncated): Ib force feedback primary modulator, group II small tonic, Ia none"],
    [DATE, "6", "LQPY5BP2", "recJ3ZtAYEvlGiThj", "Prochazka and Ellaway 2012 sensory review", "yes", "Mammals; Human", "",
     "Grounded: Zotero local abstract. Review Papers twin 'Prochazka and Ellaway 2012' = reckH8W6X4AxkLFeH"],
    [DATE, "6", "WORKFLOW", "verify", "Batch 6 COMPLETE + verified: 10 queue rows covered (41-45,47-51; row 46 McCrea 1980 already pilot-curated). 8 curated, 2 no-text (Thiry 2020, Prochazka 2017 — need Ben library pass). 4 twin records created (2 Models, 2 Review Papers), all paper-side auto-links verified live. Re-pull of all 10 ids: every field echo matched. task5 auto-curation did NOT cover the rest-import queue (only the pilot overlap) — reconcile_queue.py is the reusable reconciliation tool. Corpus export built: SADb_audit/export_corpus.py -> export/sadb_export.{json,csv} (943 records). NEXT = batch 7 = rows 52-60 (grounding for 52-60 already fetched in batch6/ground_*.txt). Papers table = 943 records, 490 with empty Notes", "", "", "", ""],
]

with open(LOG, "a", encoding="utf-8", newline="") as fh:
    csv.writer(fh).writerows(rows)
print(f"appended {len(rows)} rows to curation_log.csv")
