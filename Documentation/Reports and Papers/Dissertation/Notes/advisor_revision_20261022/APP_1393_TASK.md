# Side task (2 deliverables): (1) full 1,393-record pre-prune base in the packaged app; (2) the monthly skimmer

Written 2026-10-22 for a fresh session. Owner: Ben Bolen. Do not commit/push;
Ben commits. The repo is PUBLIC — never print or commit Airtable tokens or any
credentials.

## Verified facts (do not re-litigate; re-verify cheaply if needed)

- The pre-prune Airtable base held **1,393 records** = **943 papers** (the full
  audited corpus) + **450 linked non-paper records** (Review Papers, Models,
  Feedback pathways, Recorders).
- Post-prune base: 957 records (561 papers + 39 pathway + 244 model + 109 review
  + 4 recorder). Prune removed 436 records (382 papers + 54 orphaned non-paper
  records).
- The packaged app (`SADb_audit/app/sadb_app.html`) already embeds **all 943
  papers**, including the 382 pruned ones under an `archived` flag
  (`SADb_audit/export/sadb_export.json`: 943 records, `archived` True on 382;
  `SADb_audit/app/check_app.py` prints "papers embedded: 943").
- The app does **NOT** embed the non-paper records: the export payload has no
  Feedback / Models / Review Papers / Recorders record sets (verified 2026-10-22;
  the words appear only in descriptive text; `knowledge_base/` has derived
  pathway/afferent/animal/cluster pages, not the records).
- So "the app has the 1,393" means: add the 450 non-paper records, 396 retained
  (live in Airtable) + 54 deleted-orphaned (recoverable from the deletion
  snapshots), each under the same archived convention as the papers.

## Goal

Extend the export + app build chain so the single-file app carries the **full
1,393-record pre-prune base**, with papers and non-paper records distinguishable,
pruned/archived records flagged (papers: existing `archived` flag; non-paper:
same convention), and the Table view able to filter them. Update `check_app.py`
to assert the new counts (943 papers of which 382 archived; 450 non-paper of
which 54 archived; total 1,393).

## Where things live

- Rebuild chain (from `SADb_audit/app/README.md`):
  1. `myo python SADb_audit/export_corpus.py`   # Airtable -> export/sadb_export.{json,csv}
  2. `myo python SADb_audit/build_citation_graph.py`
  3. `myo python SADb_audit/app/build_app.py`   # -> app/sadb_app.html
  - Python: `C:\Users\Ben Bolen\.conda\envs\myo\python.exe`
- How the 382 archived papers already make it into the export: read
  `export_corpus.py` first — it merges the Airtable working base (561 papers)
  with the deletion snapshot to produce 943. Extend the SAME mechanism for the
  non-paper tables.
- Deleted non-paper record IDs: `SADb_audit/deletion_manifest.json`
  (`satellite`: Models 3, Review Papers 1, Feedback 3) and
  `SADb_audit/deletion_log.csv` (`record_ids` column; the 47 review records of
  the final sweep). Also see `restore2_items.json`, `deletion_manifest`-adjacent
  snapshots, and `corrected_purge_audit.csv` / `purge_audit_corrected.csv`.
- Retained non-paper records (396): in the Airtable working base; whatever
  access path `export_corpus.py` already uses, reuse it (do not add new
  credentials anywhere).
- App internals: `SADb_audit/app/_app_main.js` (logic, extracted for lint),
  `build_app.py` (embedding), `check_app.py` (payload sanity). The payload rides
  in a `sadb-data` script tag plus `CITES`.

## Requirements

1. Do not regress the current guarantees: offline single file, no network except
   the existing Search (online) tab, small text-only artifacts (no PDFs/binaries),
   existing filters/views keep working, citation graph unchanged (papers only).
2. Non-paper records are a NEW dimension (record type), not new tables in the
   UI unless minimal: at minimum they must be present in the payload, countable,
   and filterable by record type and archived status. A simple approach: extend
   the paper-record schema with a `recordType` field on new entries and tag
   papers `recordType: "paper"`, or add a second payload array — pick whichever
   is smaller in `build_app.py`/`_app_main.js`.
3. After building, run `check_app.py`; it must print the new counts and OK.
4. Verify totals: 1,393 records (943 + 450); 436 archived (382 papers + 54
   non-paper); 957 not archived.
5. Constraint from the 2026-10-22 session: Airtable is the editable working
   base; the app is the primary deliverable. Keep the rebuild chain one-command
   per step.

## Dissertation cross-reference (for Ben, not for this session to edit)

Chapter 5 currently says the explorer "embeds the corpus whole, including the
papers the prune archived from the working base" (true today) and states the
1,393-record pre-prune size. If this task lands, the prose can be upgraded from
"full audited corpus" to "full pre-prune base" in a later pass; the .tex files
under `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/proposals/`
are a separate pending patch set — do not edit them here.

## Definition of done

- `app/sadb_app.html` rebuilt; `check_app.py` asserts and prints 1,393 total
  (943 papers / 450 non-paper; 436 archived); Table view filters by record type
  and archived flag; everything still works offline.
- A short report file `SADb_audit/app/1393_REPORT.md` with what changed, counts,
  and the check output. Leave the working tree clean for Ben to commit.


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
  /TN SADbMonthlyWatch /TR "<myo python path> D:\Github\Bipedal_Robot\SADb_audit\watch\monthly_watch.py"`;
  Ben decides whether to activate it. A ZCode-side cron is an acceptable
  alternative but the script must work standalone.
- Constraints that stay true: text-only artifacts, no credentials in the repo
  (any key the import step needs is passed at run time, never stored), growth
  policy respected (substantive records or clearly flagged stub twins, never
  bare rows), and the run must work from the `myo` env
  (`C:\Users\Ben Bolen\.conda\envs\myo\python.exe`).

### Definition of done (deliverable 2)

- One real run completes and writes a dated report (an empty shortlist is a
  valid result — say so in the report); re-running the same month changes
  nothing; a config edit (add one author, one keyword) changes the next run.
- The two dissertation sentences above could be pointed at this script
  verbatim. Note (do not perform) any prose tweak that would then be apt in
  `1393_REPORT.md`.
