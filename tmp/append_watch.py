add = """

## Deliverable 2: the monthly skimmer (the dissertation's own promise)

Chapter 7 (70-overall.tex) already promises, verbatim: "the monthly
skimmer--curator workflow sweeps the recent literature for candidate records
and routes every one of them to a human ruling", and names it one of the
corpus's two growth engines. This task makes that sentence true in software.

### Spec

- New script `SADb_audit/watch/monthly_watch.py` (generalizes the existing
  `rr_digest_*` pattern — read `rr_digest_2023plus.md` and
  `rr_digest_round2.py` first; the mechanism is proven: OpenAlex candidate
  pull, dedupe against the corpus, shortlist, Crossref verification,
  import-ready payloads).
- **Trigger A — prominent authors**: config file `watch/watch_config.json`
  with an `authors` list. Seed it two ways and let Ben edit afterwards:
  (1) derive the top ~20 primary authors by retained-paper count from
  `export/sadb_export.json`; (2) a hand-pinned core (Rybak, Shevtsova,
  McCrea, Deng, Grillner, Pearson, Prochazka, Jankowska at minimum).
- **Trigger B — keyword collection**: `keywords` list in the same config,
  seeded from the corpus vocabulary (Ia / Ib / group II afferent, Renshaw,
  spinal locomotion / CPG, proprioceptive feedback, and the like). Ben owns
  this list.
- **Sweep window**: publications from the first day of the previous month to
  the run date (one-month lookback plus a small overlap; the script is
  idempotent — re-running a month never duplicates a report row).
- **Sources**: OpenAlex and Crossref ONLY (both open, no key, already the
  sanctioned route in the app's Search tab and the rr digest). No scraping of
  Scholar/WoS — that constraint is standing (app README).
- **Curation pipeline, exactly the established one**: candidates deduped
  against the corpus DOI set (943) and against previous shortlists
  (`watch/state.json`), collapsed preprint/version twins, then Crossref
  verification of every shortlist entry, then a one-line "why relevant"
  (which trigger fired: author / keyword hit).
- **Delivery**: a dated report `watch/reports/YYYY-MM.md` (table: DOI, year,
  first author, short title, trigger, why relevant) plus
  `watch/reports/YYYY-MM_shortlist.json` in the batch3 pattern (Name /
  Author / Year / DOI) so the existing import paths (`zotero_add_digest.py`,
  Airtable import files) work unchanged. Every entry is a CANDIDATE routed to
  a human ruling — nothing is written to Airtable or Zotero automatically
  ("the machine drafts, and the human rules").
- **Scheduling**: provide (do not auto-install) a Windows Task Scheduler
  command in the report, e.g. `schtasks /Create /SC MONTHLY /D 1 /ST 09:00
  /TN SADbMonthlyWatch /TR "<myo python path> D:\\Github\\Bipedal_Robot\\SADb_audit\\watch\\monthly_watch.py"`;
  Ben decides whether to activate it. A ZCode-side cron is an acceptable
  alternative but the script must work standalone.
- Constraints that stay true: text-only artifacts, no credentials in the repo
  (any key the import step needs is passed at run time, never stored), growth
  policy respected (substantive records or clearly flagged stub twins, never
  bare rows), and the run must work from the `myo` env
  (`C:\\Users\\Ben Bolen\\.conda\\envs\\myo\\python.exe`).

### Definition of done (deliverable 2)

- One real run completes and writes a dated report (an empty shortlist is a
  valid result — say so in the report); re-running the same month changes
  nothing; a config edit (add one author, one keyword) changes the next run.
- The two dissertation sentences above could be pointed at this script
  verbatim. Note (do not perform) any prose tweak that would then be apt in
  `1393_REPORT.md`.
"""
with open(
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\APP_1393_TASK.md",
    "a",
    encoding="utf-8",
) as f:
    f.write(add)
print("appended monthly-skimmer spec")
