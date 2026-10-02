# Sensory Afferent Database (SADb) — Explorer Package

A self-contained snapshot of the AARL literature corpus on sensory feedback
for locomotion, plus the tools to browse it. Nothing to install.

## Start here

1. **Double-click `sadb_app.html`** (any modern browser — Chrome/Edge/Firefox).
   The Table / Pivot / Map tabs work fully offline.
2. Start on the **Table** tab: search, filter by animal, feedback pathway, or
   afferent type; click a row for full notes.
3. The **Pivot** tab counts papers by any two dimensions (e.g. topic cluster ×
   decade, or feedback pathway × feedback pathway). Rows/columns sort by count
   or alphabetically. Click any cell to drill through; the "◀ Back to pivot"
   button brings you back.
4. The **Bubble map** shows every paper as a neuron (soma size = citations).
   Click any paper to enter the **focus view**: its citation network out to
   1–3 degrees of separation. Direct connections are colored — green
   open-triangle synapses mark papers citing the focus (excitatory input),
   red filled circles mark papers the focus cites — and second-degree papers
   are grayed. The left panel holds the focus paper; the right panel lists
   its connections (click to hop). Layouts: year × citations, or the topic
   landscape where proximity ≈ citation similarity.
5. The **Search (online)** tab queries OpenAlex live (needs internet) and
   badged results that are already in the corpus. Each result has buttons to
   open the same query in Google Scholar, PubMed, or Web of Science — the
   institutional-subscription sources open in your own browser session.

## Also in this package

- `knowledge_base\` — the same corpus as browsable markdown: an index plus
  one page per topic cluster, feedback pathway, afferent type, and animal,
  each listing the papers with their distilled notes.
- `sadb_export.csv` — the flat corpus (open in Excel; pivot-table source).
- `CURATION_NOTES.md` — how the corpus is curated, and what is still open.

## Data provenance (honest accounting)

- Corpus: 943 papers curated from the Airtable "Sensory Feedback" base
  (Zotero personal + AARL group libraries reconciled against it).
- Citation counts and the citation network come from OpenAlex; the in-corpus
  graph is restricted to corpus DOIs.
- Curation notes are written from each paper's abstract (OpenAlex / Europe
  PMC), not from memory; papers whose text could not be found remain
  uncured rather than guessed.
- Generated 2026-09-30 by Ben Bolen's research tooling (PSU, Agile and Adaptive
  Robotics Lab).
