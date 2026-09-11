# SAD Audit — Sensory Afferent Database three-way reconciliation

Purpose: compare the Sensory Afferent Database across (1) Ben's personal Zotero
library, (2) the AARL Zotero group library (group id **735051**), and
(3) the Airtable base that is supposed to host it (https://airtable.com).

State as of 2026-09-10 (EB475WS4 session):

## Done — Zotero legs captured (files in this folder)

Fetched from Zotero 9 **local API on port 23119** (`http://localhost:23119/api/...`,
read-only, no key needed; NOT 23127 — that port doesn't serve the API here).
Fetch script: `fetch_zotero.ps1` (paginate limit=100; PowerShell's
Invoke-RestMethod FAILS against this local server — the script shells out to
curl.exe instead; reuse that pattern).

- `personal_all_collections.json` — all personal collections
- `personal_SAD_items_full.json` / `_slim.csv` — **538** unique top-level items in the
  personal "Sensory Afferent Database" collection (RBG8ZVYP, 3 subcollections; 483 with DOI)
- `group_collections.json` — **256** AARL group collections
- `group_items_full.json` (21 MB) / `group_items_slim.csv` — **3739** unique top-level
  group items (2269 with DOI)

## Preview finding (already computed)

**Personal SAD ⊆ AARL group, essentially exactly:**
- 483/483 personal DOIs are present in the group (0 missing)
- 531/538 personal titles match a group item by normalized title
- Only 1 personal item matches neither: `DLHNVKD7` (bookSection, EMPTY title/DOI, 2015)
  — a placeholder entry, safe to ignore or delete
- `preview_personal_not_in_group.csv` — the DOI-level "personal not in group" list (empty)

## The AARL "SAD" corpus = the Biology folder subtree (Ben, 2026-09-10)

SAD material is NOT one collection in AARL — it lives under the **Biology**
collection (root `Z4NMMVVJ`) and its subfolders. Analysis in
`analyze_biology.ps1`, outputs `group_Biology_tree.csv`,
`group_Biology_items_slim.csv`, `preview_personal_not_in_Biology.csv`:

- **Biology subtree = 27 collections, 911 unique top-level items**
- Core: **Sensory Feedback** (`N9R9QB36`) with **532 items** ≈ the personal SAD set
- Rest spread over CPGs 66, Flight 53, Proprioceptive Feedback 48, Fly Eye 36,
  Neuron Modeling 36, Sensors 30, Animal Data 27, Biomechanics 25, Swimming 24,
  Neural Oscillators 23, Biped 22, Models 20, Muscles 14, and smaller ones
  (full tree with counts in `group_Biology_tree.csv`)
- Personal SAD vs Biology subtree: **537/538 found** (miss = the empty
  placeholder `DLHNVKD7`) → personal SAD ⊆ AARL Biology, confirmed

## Done — Airtable leg (2026-09-11, EB475WS4 session)

MCP tools connected (user-scope `mcp.airtable.com/mcp`). Exactly one base visible:
**"Sensory Feedback" (`appMQTnobUNRytIp7`)**, permissionLevel create. Tables:
Feedback, Papers, Models, Recorders, Review Papers. **Papers is the corpus table,
99 records total.** No DOI field anywhere ("Field 10"/"Field 11" are empty leftover
columns); titles are line-wrapped pastes with occasional typos. Review Papers is a
4-row tagging table (Prochazka 2007, Buschges 2011, Kiehn 2016, Büschges 2008), not
an extra corpus.

- Capture: `airtable_papers_slim.csv` (id/title/author/year/animals, whitespace collapsed)
- Diff: `diff_airtable.ps1` (strict norm = trim+lower+collapse-ws, loose = +strip
  non-alphanumerics — same convention as analyze_biology.ps1, plus whitespace handling
  for the line-wrapped Airtable titles)
- Adjudication + corrected gaps: `make_corrected_gaps.ps1`

### Results (three-way reconciliation)

- **Airtable (99) → AARL Biology (911):** 88 auto-matched (79 strict + 9 loose) +
  7 manual adjudications = **95/99 in Biology**. One more (Geyer 2003 "Positive force
  feedback in bouncing gaits?", 10.1098/rspb.2003.2454, INYJVQYV) is in the AARL
  group but OUTSIDE the Biology subtree. **3 Airtable papers are nowhere in Zotero**:
  Hultborn 1971 (recurrent inhibition / Ia pathway), Zill 2015 (force feedback
  insect legs), Dietz 2003 (spinal cord pattern generator). Details:
  `airtable_not_in_Biology_ADJUDICATED.csv`.
- **Biology (911) → Airtable:** 96 auto-matched rows + 10 adjudicated rows =
  **106/911 present → 805 missing** (`Biology_not_in_airtable_corrected.csv`; raw
  auto-match file `Biology_not_in_airtable.csv` beside it). 106 rows ↔ 95 Airtable
  records ⇒ the Biology subtree carries **~11 internal near-duplicate rows**
  (McVea ×3 — one with wrong year 2015, Akay ×2, Human-CPG ×3, etc.).
- **Personal SAD (538) → Airtable:** **95/538 present → 443 missing**
  (`personalSAD_not_in_airtable_corrected.csv`).
- Per-collection coverage: Sensory Feedback 72/532 on Airtable; Models 14/20;
  Animal Data 11/27; Proprioceptive Feedback 2/48; CPGs, Flight, Fly Eye, Swimming,
  Biped, Muscles, Synergies etc. 0 — Airtable is a small curated subset, not a
  mirror.

### Decisions parked for Ben (nothing written to Airtable)

1. Bulk-import gap into Papers? Options: all 805 Biology items, the 443
   personal-core items, or the 460 Sensory-Feedback-only items. MCP create works in
   50-record batches; would want a DOI field added first (base has none).
2. Add the 3 Zotero-missing papers to Zotero by DOI (keeps Zotero the master).
3. Fix 4 Airtable title typos (stick inspect→insect, walkde→walker, truncated
   Hiebert title, Akay dropped words) — one-line edits, need his OK.
4. Zotero-side cleanup of the ~11 Biology internal duplicates — his call.

## Second correction round (2026-09-11 evening, same session)

**Zotero purge criterion — Ben's rule was FILE-level, I applied ITEM-level.**
Corrected audit (`corrected_purge_audit.ps1` → `purge_audit_corrected.csv`):
Airtable has files on 98/149 records (essentially all originals; the 50 new
records are bare). Under the strict rule (file in BOTH AARL and Airtable):
**78/80 purges were correct; only 2 over-purged** — AnimatLab (Cofer 2010,
recf6gvkwsPmIHY3N) and Büschges 1995 pilocarpine (rec0k5AJjSxNrEexP) — and
BOTH still have their PDFs in Airtable (recoverable by download; AARL lacks
their files). Web-API delete remains permanent — restore path = download from
Airtable attachment → drag onto the Zotero item (GUI) or API file upload.

**Airtable curation pass on the 50 (the "missing notes/attachments" complaint):**
- **14 Review Papers records created + linked** to their Papers entries
  (Prochazka 1999, Duysens 2004, Kiehn 2006, Brownstone and Bui 2010, Hultborn
  1998, Nishikawa 2007, Arber 2012, Grillner and Kozlov 2021, Perreault and
  Giorgi 2019, Li 2023, Sengupta and Bagnall 2023, Chiel and Beer 1997,
  Grillner 1995, Mori 1987).
- **9 Models records created + linked** ("this paper is this model" via
  Models.Paper): Deng 2022, Szczecinski 2015, Yakovenko 2018, Boxerbaum 2011,
  Deng 2018, Lyttle 2017, Pickard 2020, Bui and Brownstone 2015, Danner 2016.
- **Notes written on all 50 papers** (draft quality, Alex-template style —
  distilled findings, NOT abstracts; flagged for Ben's review pass), **Animal
  tags filled on 44**, **Feedback links added on 9 clear cases** (cutaneous
  stance modification, Ia reciprocal inhibition, Ia or II stance-to-swing,
  Ia monosynaptic excitation, fictive-without-feedback).
- Still open on the 50: PDF attachments (files live in Zotero/AARL, not
  uploaded to Airtable — bulk upload is a Ben decision), deeper Feedback-type
  elaboration per Ben's diagram-style naming (needs his insight pass).

## Ben's rulings + follow-up session (2026-09-11, same EB475WS4 session)

### Rulings
- **AARL Zotero duplicates: DO NOT clean up.** They may feed separate .bib files.
- Approved: add 3 missing papers to Zotero, fix Airtable typos + years, DOI field
  (done, see below), and the personal-SAD attachment purge with two-part criterion:
  remove attachment ONLY if (not dissertation-cited) AND (already in both AARL
  Zotero and Airtable).
- Ben wants batch writes verified ("show me a 50-record batch done right") and
  GLM-5.3 review — this session runs GLM-5.3; an independent audit subagent also
  re-verified all writes. **Audit verdict: OVERALL PASS** (all 4 checks: 90 DOI
  backfills exact, 21-record correction batch exact, 6 renames exact incl. α/ü,
  6 expected blanks + total DOI count 93). Post-audit fill: Poppele & Bosco Year
  → 2003 (DOI 10.1016/S0166-2236(03)00073-0; audit flagged it as the one missing
  year).

### Base schema semantics (decoded from link data — answers Ben's questions)
- **Papers."Models"** = models the paper references/was influenced by (symmetric
  twin of Models."Papers Cited"). The "was influenced by models" field Ben wanted
  ALREADY EXISTS — it's this column.
- **Papers."Models 2"** = "this paper IS this model" (twin of Models."Paper", the
  model's source paper). Ben's reading confirmed.
- **Papers."Models copy" (fldBLowhKJcjuFMSK, last one)** = link to **Review Papers**
  entries ("the review that covers this paper") — NOT models. Only 3 papers linked.
  The other "Models copy" (fldh983rtt2YtMZQX) is empty everywhere — dead field.
- **Feedback.Papers** (fldcSKk0TWYimoWy) is empty everywhere; the live link is
  Papers."Feedback" from the paper side.
- **Models tab** = conceptual models of sensory afferents, named "Author Year"
  (19 records); "Rules Used" links the feedback pathways the model implements.
- **Review Papers tab** = literature reviews; "Paper" links the review's own Papers
  entry; "Papers Cited" (currently empty) is where in-library cited papers should go.
- **Recorders**: 4 people (Alexander Hunt, Sabian Limones, Ben Bolen, Kaiyu Deng).
  **All Feedback entries without a recorder are Ben's** (his rule). Alex's entries
  are the oldest (2020-12) and are the curation template.
- **Curation template (from Alex's exemplars)**: paper first in Papers; Notes =
  distilled insight sentences usable in our writing (NOT abstract pastes — the
  Geyer&Herr 2010 note is the counter-example); Feedback names should be
  pathway+function+phase in the style of Ben's reflex diagrams ("Ib stance to
  swing", "type II contralateral stance-to-swing"), not bare "type Ia monosynaptic"
  (acceptable only for pathway-existence experiments); Animals = all relevant
  animals.

### Airtable writes DONE this session (all verified + audited)
1. **DOI field created** at end of Papers (`flddp7XZI4RklGWzq`), then backfilled
   from Zotero matches: **93/99 papers have DOIs** (2 batches of 50+40, then 3 more).
   6 blanks: Prilutsky book, Deng 2019, Herr 2002, Brown 1914, Ijspeert 2007,
   Perreault 1995 — Zotero rows for those lack DOIs too.
   One conflict resolved by hand: Latash 2020 filed twice in Zotero (bioRxiv +
   published) → wrote the published DOI 10.3389/fnins.2020.598888.
   Mapping file: `airtable_doi_backfill.csv`; generator `gen_doi_backfill.ps1`.
2. **Correction batch (21 Papers records)**: 7 title fixes (McVea word order,
   stick inspect→insect, α-motoneuron, Hiebert truncation completed, walkde→walker,
   Akay dropped "Proprioceptive and", Dietz singular→"Generators"), 7 year fixes
   (Ivashko→2003, Shevtsova V1→2022, Nichols→2018, Danner eLife→2017, Shevtsova
   two-joint→2016, Frigon split-belt→2016, Perreault→1995), 6 author typo fixes
   (Ryback→Rybak, Shevstova→Shevtsova ×3, trailing comma/period/space ×3),
   3 new DOIs (Hultborn/Zill/Dietz). Szczecinski stays 2014 (volume year; Zotero
   holds online-first 2013). Year checker: `check_years.ps1`.
3. **Renames**: Feedback "la or ll"→"Ia or II swing to stance", "Ipsilatoral"→
   "Ipsilateral", "Chordontonal Organ"→"Chordotonal organ"; Models "Shevstova
   2015/2016"→"Shevtsova", "Durr 2004"→"Dürr 2004".
   `airtable_papers_slim.csv` re-patched to post-fix state (now has a DOI column).
- Flagged for Ben, NOT changed: Models-name years "Van Der Noot 2017" (paper is
  2018) and "Schilling and Cruse 2015" (linked paper is the 2020 one) — these are
  lab identifiers; rename only on his word. Select option "Chordontonal organ"
  (Sensory System) also misspelled — option rename not exposed via current tools.

### Zotero 9.0.6 local API is READ-ONLY → SOLVED via Ben's web-API key (2026-09-11)
POST /api/users/0/items → "Endpoint does not support method" (local). Ben supplied
a web-API key (userID 631450, username bbolen, user-library write + files — the key
is passed as a script param and NOT stored in any file). Via api.zotero.org:

- **3 papers ADDED to personal SAD (RBG8ZVYP)**, 200/3-successful/0-failed:
  Hultborn 1971 = `6MNCI8I6`, Zill 2015 = `8H7VACP4`, Dietz 2003 = `928TAKM7`.
  Script: `zotero_web_add3.ps1` (payload `add3_items.json`).
- **Attachment purge EXECUTED: 107 attachments deleted from 80 personal-SAD items**
  (all matched Ben's criterion: not dissertation-cited, on both AARL + Airtable;
  re-verified current parent before deleting — 107/107, 0 drift). Script:
  `zotero_purge.ps1`. **IMPORTANT: Zotero web-API DELETE is PERMANENT — it bypasses
  the client trash** (confirmed: deleted keys now 404 "Item not found"). Files stop
  counting against storage immediately.
- PS 5.1 gotchas hit (do not relearn): `@($json | ConvertFrom-Json)` collapses a
  JSON array to ONE object — assign first (`$items = $json | ConvertFrom-Json`)
  then `@($items)`; inline `&`-containing URLs get eaten by cmd quoting — write
  .ps1 files; curl.exe output must be read via `-o file` + `Get-Content -Raw`.

### First 50-record CREATE batch (2026-09-11, after Ben's "I don't see any additions")
Ben's read was right: everything above was edits, no new rows. Executed the create
demo: **50 new Papers records** from `personalSAD_not_in_airtable_corrected.csv`
(his personal-SAD gap list; first 50 clean items). Fields: Name/Primary Author/
Year/DOI. **Table now 149 records.** Files: `airtable_create50.csv` (payload),
`airtable_created50_ids.csv` (zotero_key → airtable_id mapping), builders
`build_create50.ps1` + `fix_create50.ps1`. Gotchas hit: personal full JSON has
DOUBLE creators (single-"name" variant + firstName/lastName variant — dedupe by
surname); two mangled Zotero items EXCLUDED from import for manual curation:
`C6252AEY` (title is a PDF filename "Horcholle-Bossavit - Encoding of muscle
contractile tension...") and `QDN8VDSK` ("catsof ankle extensor muscle activity
in walking" — fused/truncated title); `MCU88JT6` Zotero author mangled ("Ml") —
imported as "Shik et al." (Shik/Severin/Orlovsky 1966). Two Zotero-side record
fixes suggested for Ben eventually. Remaining gap after this batch: 393 personal
items (443 - 50) not yet on Airtable. Independent audit subagent re-verified
(separate from the earlier edit-batch audit). **Create-batch audit verdict:**
CHECK 2 (count = 149) PASS, CHECK 3 (no duplicate DOIs/titles table-wide) PASS,
CHECK 1 (50 records vs CSV) 47/50 exact — the 3 "failures" were mojibake in the
LOCAL CSV artifact (cp437 double-encoding of Procházka / Lafrenière-Roula /
en-dash; Airtable held the correct characters, proven by codepoint round-trip),
plus one deliberate deviation: 3SB4QSIF DOI written uppercase JNEUROSCI to match
the base's house convention (DOIs are case-insensitive). CSV repaired in place
(2 passes). **PS 5.1 trap discovered: .ps1 files saved UTF-8 WITHOUT BOM get
their non-ASCII literals (Γ etc.) misread as ANSI — use codepoint literals
([char]0x0393) in scripts, never paste mojibake-lookup regexes as literals.**
The original CSV mojibake source: console round-trips of curl/PS output —
always author Airtable payloads from -Raw file reads, not console echoes.

### Attachment purge eligibility (computed, THEN EXECUTED — see Zotero section)
`attachment_eligibility.ps1` scanned the dissertation (126 .tex under
Dissertation\, 132 cited keys, 99 cited DOIs + 129 titles), fetched 607 personal
attachments via local API (567 file attachments), and applied Ben's criterion:
**107 attachments on 80 items were eligible** (`personal_SAD_attachment_purge_eligible.csv`);
17 dissertation-cited items with attachments were protected. Executed permanently
via web API (see above). Personal SAD ∩ Airtable = 96/538 (post-fix).

## Notes

- Zotero storage cap: metadata sync is free; only attachment files count. This audit
  uses metadata only — no storage impact. Attachment reconciliation (if wanted later)
  needs the web API with an API key.
- `group_items_full.json` is 21 MB — keep this folder out of git, or gitignore the
  big JSONs (Ben commits via GitHub Desktop; folder is currently untracked).
