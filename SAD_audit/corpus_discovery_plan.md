# Corpus discovery tool plan — a private "Research Rabbit" for the SAD corpus

Written 2026-09-12 (laptop session). Goal: RR-style discovery (similar work + newest
articles + bubble graph) over OUR corpus (Airtable Papers 532 / personal Zotero SAD 538),
tied to the SAD curation tags instead of locked in someone else's web app.

## Why not just keep using Research Rabbit
RR is good for interactive triage, but: its collections are separate from SAD (no
Animals/Feedback tags), no programmatic access, and recommendations can't be batch-audited
or preserved. The plan below reproduces the two features Ben actually uses — "similar
work" graphing and "newest articles" feeds — on top of the corpus we already curate.

## Data spine (no new infrastructure)
- **OpenAlex** (free, no key; polite pool with `?mailto=`): DOI → work record gives
  `referenced_works` (its reference list), `cited_by_count`, `cited_by_api_url`,
  `topics`/`concepts`, year, venue. Citation edges between corpus members and
  recent outside papers are all derivable from this one endpoint.
- **Semantic Scholar Graph API** as fallback/supplement (abstracts + TLDRs; rate-limited,
  needs free key for real volume).
- **Corpus inventory**: Airtable Papers (has DOI field now, 93% coverage) +
  `personal_SAD_items_slim.csv` for Zotero-side metadata.
- Runs on any machine's existing python env (laptop `myoconv`, EB475WS4 `myo`,
  easteregg2 `myo`) — only needs `requests`; pyvis/networkx optional for phase 3.

## Phases

### Phase 1 — Corpus enrichment (build once, refresh on demand) — ~half a day
Script `oa_enrich.py`: pull every Papers DOI through OpenAlex; cache
`work_id, referenced_works, cited_by_count, topics, year` to one JSON/CSV in
`SAD_audit\`. Output: coverage report (papers OpenAlex can't find; missing DOIs to backfill).

### Phase 2 — Bubble graph, zero-code route (evaluate first — Ben's VOSviewer call)
Export the enriched corpus as a VOSviewer "network/map" CSV (label, year, citations,
+ optional `references` for bibliographic coupling). VOSviewer (free desktop) renders
the citation-network bubble map in one click: bubble size = cited_by_count, cluster =
co-citation community, labels on zoom. If that answers the "how do we see our corpus"
question, stop here for graphing.

### Phase 3 — Interactive HTML map tied to SAD tags (scriptable route)
`corpus_map.py` → pyvis/networkx HTML in `SAD_audit\maps\`:
- nodes = corpus papers (+1-hop cited works, greyed)
- edges = citation links (OpenAlex `referenced_works` ∩ corpus, both directions)
- node size = cited_by_count; color = Airtable **Animals** (or Feedback pathway)
- hover = title/year/note snippet; filters: year slider, tag picker.
This is the shareable lab version of RR's graph, colored by OUR vocabulary.

### Phase 4 — Recommendations engine (the real RR replacement)
`recommend.py` (weekly, cron-able):
1. Seeds = corpus papers (optionally: a pinned subset, e.g. CPG/afferent cores).
2. Candidate pool = OpenAlex `cited_by_api_url` of seeds + topic filters, restricted to
   `from_publication_date` = last 2–3 years; score = Σ(cites our seeds) + topic overlap
   + recency + venue; drop anything already in Papers (DOI match).
3. Output a **Markdown digest** (`SAD_audit\digests\YYYY-WW.md`): ranked table
   (title/authors/DOI/why-recommended/abstract link). Human triage: mark
   import/skip. `--import` creates Airtable Papers rows flagged for the normal
   curation pipeline (notes grounded before any write — same rules as batches).
This is exactly the RR "Similar Work"/"Newest Discoveries" loop, but the candidates,
scores, and imports stay ours.

## This session's head start (2026-09-12)
Ben's live RR collection was inspected and its recent recommendations enumerated
(the 2023–2026 items) as the first digest cohort — see `curation_log.csv` and the
session summary. The RR collection itself can be exported (all citations → BibTeX)
and used as RR's opinion of "similar", to benchmark phase 4's scoring.

## Decisions for Ben
1. VOSviewer (GUI, instant) vs pyvis HTML (scriptable, tag-colored) — or both.
2. Digest cadence (weekly?) and whether imports auto-create Airtable rows or wait
   for his triage click.
3. Scope of "corpus" for graphing: all 532 Papers, or the curated/Feedback-linked
   subset only.
