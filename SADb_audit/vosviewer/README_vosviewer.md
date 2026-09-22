# VOSviewer integration — SADb corpus (full-corpus rebuild 2026-09-18; self-contained map 2026-09-21)

Research-Rabbit-style bubble map of the whole Papers corpus, built from
OpenAlex citation data. Rebuilt by `SADb_audit\vos_build2.py` (the older
`vos_build.py` covered only 552 papers and is superseded). Files here:

| file | what it is |
|---|---|
| `sadb_map.txt` | VOSviewer MAP: 943 papers — id (OpenAlex W-id), label, doi, year, first_author, weight = cited_by_count, **x, y, cluster** (built-in layout + 18 Louvain clusters, same as the HTML app), source. SELF-CONTAINED: loads with no prompts. |
| `sadb_network.txt` | VOSviewer NETWORK: 11,108 citation edges (optional — only draws links on top of the map layout) |
| `openalex_enrichment.csv` | per-DOI: oa_id, cited_by_count, n_refs, in-corpus refs |
| `openalex_misses.txt` | DOIs OpenAlex could not resolve (6 of 892) |
| `sadb_preview.png` | static matplotlib preview (NOT a VOSviewer product) |

## How to open in VOSviewer

Desktop (vosviewer.com, free, needs Java) **or** the web app `https://app.vosviewer.com`:

**Easiest: drag `sadb_map.txt` onto the page by itself.** It carries its own
positions and clusters, so it renders immediately — no "create map based on a
network" prompt, no x/y warning. Then set **Size by = weight** (citations) and
**Color by = cluster** (our 18 topic clusters, matching the HTML app).

**Optional: add `sadb_network.txt`** to draw the citation links on top. Drop it
onto the page second, or in the Open dialog put the MAP in the map-file slot
and the NETWORK in the network-file slot — swapping them throws
"line 1 missing or incorrect" (the network header is not a map header).

## Reading it

- 816/943 papers are in the citation-connected component; 127 singletons are
  mostly 2023-2026 papers too new to have in-corpus citing links (they sit on
  the outer ring).
- Bubble size = OpenAlex cited_by_count (all-citations, not in-corpus).

## Regenerating

Run the standard chain (reads a fresh `export/sadb_export.json`, OpenAlex
results are cached in `export/openalex_raw.json`):

```
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" SADb_audit\vos_build2.py
```
