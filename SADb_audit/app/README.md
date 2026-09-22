# SADb Explorer — `sadb_app.html`

Single-file, fully OFFLINE explorer for the Sensory Afferent Database corpus
(943 papers). Double-click it — no server, no network, no install. Data is
embedded at build time from the Airtable export; rebuild after any curation
batch.

## Rebuild chain

```
myo python SADb_audit/export_corpus.py          # Airtable -> export/sadb_export.{json,csv}
myo python SADb_audit/build_citation_graph.py   # OpenAlex (cached) -> export/sadb_cites.json + sadb_layout.json
myo python SADb_audit/app/build_app.py          # -> app/sadb_app.html
```

`export\sadb_export.csv` is also the Excel-ready flat corpus (the pivot-table
source independent of Airtable).

## What's inside

- **Table** — search (title/author/year/DOI), column sort, filters (animal,
  feedback pathway, import source, year range, has-notes, has-PDF); click a row
  for the detail drawer (curation note, DOI link, tags).
- **Pivot** — count papers by any two dimensions (cluster × decade default;
  author, animal, pathway, source, notes, PDF...). Click any cell to drill
  through to the filtered table.
- **Bubble map** — two Layouts:
  - *Year × citations* (publication year vs log citations),
  - *Topic landscape* (spring layout of the in-corpus citation network, 18
    Louvain clusters — proximity ≈ citation similarity).
  - **Style: Neurons** — each paper is a stylized neuron (soma sized by
    citations, dendrites, axon stub). Focusing/spotlighting a paper draws its
    synapses: **excitatory (open triangles) from the papers that cite it**,
    **inhibitory (filled circles) onto the papers it cites** (Ben's semantic
    assignment, 2026-09-22; shapes follow his SNS diagram convention).
  - *View*: All papers / References / Cited by / Research review / Animal
    studies / Models. Colors = Ben's accessible 7-color set (`Code\Matlab\Colors.m`).
  - Interactions: left-click = focus (swappable with details via "Left click
    action"), right-click / Ctrl+click (Mac) / Shift+Enter = details, Tab +
    arrows + Enter keyboard path, Esc = back, wheel = zoom, drag = pan.

## Roadmap (Ben, 2026-09-22)

Airtable remains the primary data tool; this app is the corpus browser.

1. **v1 — OFFLINE (done)**: this file. Snapshot-based; rebuild after batches.
2. **v2 — ONLINE mode**: the same UI fetching live data — Airtable REST (or a
   small hosted JSON) for records + OpenAlex for fresh citation counts, with
   Web of Science (needs institutional subscription; lab has PSU access) and
   Google Scholar (no official API — only via sanctioned tools, ToS-sensitive)
   as enrichment sources. OpenAlex remains the default citation source.
3. **v3 — Browser extension** (Zotero-connector-like, Manifest V3): on a
   publisher page, capture DOI/metadata, resolve the PDF through the user's
   university-library access, and push to Zotero + Airtable. This turns every
   lab member with library access into a PDF-hunting node (kills the
   remaining_no_pdf backlog organically).
4. **v4 — Multi-lab deployment**: user base = AARL plus the collaborating labs
   on the grant (Quinn, Chiel, Büschges, Szczecinski, Webster-Wood). Hosted
   read-only snapshots (GitHub Pages works for the offline file) + shared
   curation backends (Airtable base shared with editors per lab); curation
   conventions per the SADb spec so cross-lab notes stay compatible.

Constraints that stay true at every stage: small text-only artifacts in the
repo (no PDFs/binaries — repo-size campaign), no proxied URLs, Ben approves
any Airtable schema change.
