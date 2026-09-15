# VOSviewer integration — SADb corpus (built 2026-09-14, EB475WS4)

Research-Rabbit-style bubble map of the whole Papers corpus, built from
OpenAlex citation data. Files here (regenerable, all small):

| file | what it is |
|---|---|
| `sadb_map.txt` | VOSviewer MAP: 552 papers — id (OpenAlex W-id), label, doi, year, first_author, weight = cited_by_count, source (originals99/demo50/rest383/digest20) |
| `sadb_network.txt` | VOSviewer NETWORK: 9,400 citation edges between corpus members (OpenAlex referenced_works, deduped, weight = reciprocal-citation count) |
| `openalex_enrichment.csv` | per-DOI: oa_id, cited_by_count, n_refs, in-corpus refs (sorted by citations) |
| `openalex_misses.txt` | DOIs OpenAlex could not resolve (0 on first build — all 502 DOIs resolved) |
| `sadb_preview.png` | static matplotlib preview (NOT a VOSviewer product) |
| `vos_build.py` | rebuilds all of the above from the local CSVs + OpenAlex |
| `vos_preview.py` | rebuilds the PNG (networkx spring layout) |

## How to open in VOSviewer

Desktop (vosviewer.com, free, needs Java) **or** the web app `https://app.vosviewer.com`:

1. **File → Open → VOSviewer map file…** → pick `sadb_map.txt`.
   If it asks, choose to create the map based on a network.
2. **File → Open → VOSviewer network file…** → pick `sadb_network.txt`.
3. In the right panel: **Size by = weight (citations)**; **Color by = cluster**
   (VOSviewer computes its own topic clusters from the citation network) or
   by weight for a citation-density heat view. Zoom/pan freely; search box
   matches labels/authors.

The web app also accepts both files dragged in together.

## Reading it

- The big central component (495/552 papers) is the citation-connected literature;
  57 singletons are mostly 2023-2026 papers too new to have accumulated in-corpus
  citing links (orange digest20 dots sit near the core they were mined from).
- Bubble size = OpenAlex cited_by_count (all-citations, not in-corpus).

## Regenerating after future curation batches

1. `vos_build.py` reads four local CSVs (originals99, demo-50 via personal_SAD_items_slim,
   rest-import 383) plus **the inline `DIGEST20` list** — extend that list (or replace it
   with an Airtable export) as new papers are created. Then rerun:
   `"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" vos_build.py`
2. `python vos_preview.py` for the PNG.

Known cosmetic issue: a few `originals99` labels carry cp437 mojibake from the legacy
slim CSV (e.g. "BA¼schges"); Airtable holds the correct characters — re-export labels
from Airtable if it bothers you.
