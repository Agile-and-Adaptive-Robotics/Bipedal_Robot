# 1393_REPORT — SADb app now carries the full pre-prune base (2026-10-23)

## What changed

The single-file app ([sadb_app.html](sadb_app.html)) previously embedded only the
943 papers. It now embeds the **full pre-prune base: 1,393 records** — the 943
papers plus the 450 non-paper records (Review Papers, Models, Feedback pathways,
Recorders), every record tagged with a record type and carrying the archived
flag under the same convention as the papers.

Changed files:

- `SADb_audit\export_corpus.py` — exports the retained non-paper records from the
  Airtable working base using the SAME access path and credential source as the
  papers (no new credential anywhere); reconstructs the 54 deleted orphan
  non-paper records from the local deletion snapshots and flags them
  `archived=True`; every record gets `recordType`; the export is assembled in
  memory and written ONCE at the end, so a failed Airtable pull can never clobber
  the previous export. The Airtable base was only read, never modified.
- `SADb_audit\build_citation_graph.py` — one added line: the graph/layout input
  is filtered to `recordType == "paper"`, so the citation graph stays papers-only.
- `SADb_audit\app\build_app.py` — payload now carries all records (compact key
  `rt` = record type; papers keep every existing field and flag). The Table view
  gains two filters — **Record type** (paper / review / model / feedback pathway /
  recorder) and **Archive** (live / archived app-only). Pivot, Bubble map, and the
  citation graph stay papers-only; drill-through now matches exactly the drilled
  paper cell; the detail drawer shows the record type and hides paper-only
  citation actions on non-paper records.
- `SADb_audit\app\check_app.py` — asserts the full pre-prune counts (below) and
  still extracts `_app_main.js` for `node --check`.
- `SADb_audit\app\README.md` — counts and filter documentation updated.

Papers were not touched: same 943 records with the same 382 archived flags as
before this change.

## Orphan reconstruction (54 archived non-paper records)

Retained non-paper records (396) come live from Airtable. The deleted orphans are
rebuilt from the deletion snapshots of record:

- **Models 3** — `deletion_manifest.json` satellite ids (Shevtsova 2016, Geyer and
  Herr 2010, Rybak 2006), names from `prune_shortlist.csv`.
- **Review Papers 48** — 1 manifest satellite (Prochazka 2007, name from
  `prune_shortlist.csv`) + the 47 orphaned Review Papers deleted 2026-10-02,
  id→name from `review_orphans_20261002.json`.
- **Feedback pathways 3** — manifest satellite ids, names from
  `prune_shortlist.csv` (one has an empty Name in Airtable; carried as empty).
- **Recorders 0** — all 4 Recorders are live; nothing to reconstruct.

`export_corpus.py` cross-checks the reconstructed id set against the manifest
satellite ids (asserts nothing missing and nothing beyond manifest + review
orphans) and never duplicates a record that is live in Airtable.

## Counts (read back from the final rebuilt `export\sadb_export.json`)

| | total | archived | live |
|---|---|---|---|
| Papers | 943 | 382 | 561 |
| Models | 247 | 3 | 244 |
| Review Papers | 157 | 48 | 109 |
| Feedback pathways | 42 | 3 | 39 |
| Recorders | 4 | 0 | 4 |
| **All records** | **1393** | **436** | **957** |

## check_app.py output (run 2026-10-23)

```
records embedded: 1393
papers: 943 (382 archived)
non-paper: 450 (54 archived) by type: {'model': 247, 'review': 157, 'feedback': 42, 'recorder': 4} | archived by type: {'model': 3, 'review': 48, 'feedback': 3}
archived total: 436 | live: 957
with notes: 897
with afferents: 95
clusters: 18
citing papers: 730 | pairs: 11130
cluster labels: 18
extracted 36 KB of JS; placeholders left: 0
OK
```

Rebuild chain (one command per step, myo python):

```
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" SADb_audit\export_corpus.py
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" SADb_audit\build_citation_graph.py
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" SADb_audit\app\build_app.py
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" SADb_audit\app\check_app.py
```

`node --check app\_app_main.js` on the extracted JS: clean. The app remains one
offline HTML file; text-only artifacts only; the Search (online) tab is still the
only network surface.

## NOTE (not performed here)

The dissertation's Chapter 5 can later be upgraded from "full audited corpus" to
"full pre-prune base" now that the app genuinely carries all 1,393 pre-prune
records. That wording change belongs to the separate .tex patch set and was
deliberately NOT made in this change (nothing under `Documentation\` was
touched).
