# SADb Explorer — `sadb_app.html`

Single-file explorer for the Sensory Afferent Database corpus. The core
(Table / Pivot / Map) is fully OFFLINE — double-click it, no server, no
network, no install. Data is embedded at build time from the Airtable
export; rebuild after any curation batch. The **Search (online)** tab is the
only part that touches the network.

## Rebuild chain

```
myo python SADb_audit/export_corpus.py          # Airtable -> export/sadb_export.{json,csv}
myo python SADb_audit/build_citation_graph.py   # OpenAlex (cached) -> sadb_cites.json + sadb_layout.json
myo python SADb_audit/app/build_app.py          # -> app/sadb_app.html
```

`export\sadb_export.csv` is also the Excel-ready flat corpus (the pivot-table
source independent of Airtable). Checks: `app\check_app.py` (payload sanity +
extracts the JS for `node --check app\_app_main.js`).

## What's inside

- **Table** — search (title/author/year/DOI/notes), column sort, filters
  (animal, feedback pathway, afferent type, import source, year range,
  has-notes, has-PDF); click a row for the detail drawer (curation note,
  afferents, animal-study potential, robot/sim translation, DOI link, tags).
- **Pivot** — count papers by any two dimensions. Topic clusters are NAMED
  (18 labels from `export\cluster_labels.json`, e.g. "Afferent control of
  locomotion (classics)"). Rows and columns each sort by count, A→Z, or Z→A
  (Primary Author alphabetical finally possible — Rybak sits as secondary
  author on many papers, count-order hides that). "Import source (legacy)" is
  demoted to the end of the dimension list. Click any cell to drill through;
  the Table view then shows a **"◀ Back to pivot"** button so you can always
  dive back out.
- **Bubble map** — two layouts: *Year × citations* and *Topic landscape*
  (spring layout of the in-corpus citation network; proximity ≈ citation
  similarity). Style: Neurons (soma sized by citations, dendrites, axon stub).
  - **Research-Rabbit-style FOCUS MODE (2026-09-28)**: click a paper → the
    map becomes its citation neighborhood out to 1 / 2 / 3 degrees
    (selector). **Direct connections are colored; 2nd-degree nodes are
    grayed**; everything else drops out. Green open-triangle synapses =
    excitatory input from papers that CITE the focus; red filled-circle
    synapses = the papers the focus CITES (Ben's semantic assignment,
    2026-09-22). **Left panel** = the focus paper (details, actions);
    **right panel** = the connection list (cited-by / references / 2nd
    degree), click any row to hop the focus there. Esc = back to all papers.
  - *View*: All papers / References / Cited by / Research review / Animal
    studies / Models. Colors = Ben's accessible 7-color set
    (`Code\Matlab\Colors.m`).
  - Interactions: wheel = zoom, drag = pan, right-click/Ctrl+click/Shift+Enter
    = details (swappable via "Left click action"), Tab + arrows + Enter
    keyboard path, Esc = back.
- **Search (online)** — live OpenAlex query (CORS-open, no key, needs
  internet) with abstracts, citation counts, and a badge when a DOI is
  already in the corpus. Each result carries launch buttons that open the
  query in YOUR browser session: Google Scholar, PubMed, Web of Science
  (opens the WoS basic-search page with the query copied to the clipboard —
  PSU login applies there), plus doi.org and OpenAlex links and a
  "⤓ suggest for corpus" download that drops a small JSON the curator can
  import later. Scraping Scholar/WoS from inside the app is deliberately NOT
  done (their terms of service; institutional auth lives in the browser).

## Roadmap (Ben, 2026-09-22)

Airtable remains the primary data tool; this app is the corpus browser.

1. **v1 — OFFLINE (done)**: snapshot browser.
2. **v2 — ONLINE mode (done 2026-09-28 in-app)**: live OpenAlex search tab;
   WoS + Scholar via sanctioned launch buttons.
3. **v3 — Browser extension** (Zotero-connector-like, Manifest V3): on a
   publisher page, capture DOI/metadata, resolve the PDF through the user's
   university-library access, and push to Zotero + Airtable.
4. **v4 — Multi-lab deployment**: hosted read-only snapshots (GitHub Pages
   works for the offline file) + shared curation backends.

Constraints that stay true at every stage: small text-only artifacts in the
repo (no PDFs/binaries — repo-size campaign), no proxied URLs, Ben approves
any Airtable schema change.
