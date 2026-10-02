# Column renames — applied vs. proposed (2026-09-28)

Ben asked for "the columns renamed like we talked about." That conversation is
not in this repo, so this file records exactly what was renamed (all
field-ID-stable, scripts keep working) and what still needs his word.

## Applied (2026-09-28, field IDs unchanged)

| Field ID | Old name | New name | Why |
|---|---|---|---|
| fldBLowhKJcjuFMSK | Models copy | **Review Papers** | it IS the live reverse link to the Review Papers table (this also kills the duplicate-name GET-422 bug) |
| fldh983rtt2YtMZQX | Models copy | Models copy (unused, empty) | empty everywhere; renamed to uniquify — delete on Ben's word |
| fldnoM8tlsUo7SPJA | Field 10 | Field 10 (unused) | junk name; contains one stray `pat-write-probe` artifact |
| fldXcGzCtnfISOcEn | Field 11 | Field 11 (unused) | junk name, empty everywhere |
| fldeiA2v7KuPeR0b6 | (auto) | Cited by (in corpus) | inverse of the new Cites self-link |

## New fields added (2026-09-28)

| Field | Type | Purpose |
|---|---|---|
| Afferent Types (fldcWSx0wyFJHjYZc) | multi-select | Type 1 (legacy), Ia, Ib, II, III/IV, Mechanosensory, Cutaneous, Heat, Nociceptive, Flexor reflex afferents |
| Animal Study Potential (fldYtDbTshGMbRK6w) | long text | how the paper's info could seed / be tested in an animal study |
| Robot/Sim Translation (fld5PGyWPjLYONLTD) | long text | which pathway / synapse / loss-of-function fact enables sim or robotic testing |
| Cites (in corpus) (fldqQUpWF6Lbi6hdp) | self-link | in-corpus citation edges (730 papers, 11,130 directed pairs) |
| Prune Status (fld9jhnEQW738POne) | single-select | prune-candidate / keep / pruned — flags only, deletions are Ben's |
| Prune Reason (fldcW2N1nn2pNxfv0) | text | one-line justification |

## Proposed but NOT applied (need Ben's OK)

1. **"Models 2" → "Is the model paper"** — that field (reverse of the Models
   table's Paper link) is non-empty exactly when the paper IS a model; the
   current name reads like a duplicate of "Models".
2. **"Models" → "Models referenced"** — the forward link, models the paper
   uses (its description already says this).
3. Delete the three "(unused)" fields outright.
4. Feedback-table "Models copy" (fldjcTBtIC1xHzl7j) — same duplicate-name
   family, unused; rename or delete.

Renaming any of these is one MCP call plus a one-line FETCH update in
`export_corpus.py`; say the word.

---
Ben's ruling (file comment, 2026-10-01): "Everything here looks good. Go for it."

## Status after execution (2026-10-01)

- ✅ "Models 2" → **Is the model paper** (fldbIkufcnBjpnUln) — applied;
  `export_corpus.py` + app labels updated (rebuild pending API budget).
- ✅ "Models" → **Models referenced** (fldR641pV7rYou7jA) — applied; same.
- ✅ Feedback "Models copy" → **Models copy (unused, empty)** (fldjcTBtIC1xHzl7j)
  — verified used on 0 of 38 records, renamed.
- ❌ The three Papers junk-field DELETIONS (Field 10/11, empty Models copy):
  the REST + metadata field-DELETE endpoints both 404 with this PAT (no
  schema scope) — **Ben: 2 clicks each in the Airtable UI** (field menu →
  Delete field). They're all empty, nothing can be lost.

## Twin-record audit (same session; twins_profile.json)

Models 247 = 230 stubs + 14 rich (Notes/Rules Used/Papers Cited) + 3 orphans.
Review Papers 157 = 156 stubs + 0 rich + 1 orphan. The 386 stubs carry NOTHING
but a name + link to their paper — they are the "repeats" to delete for the
record cap, NOT the papers (see prune_review_20261001.md / Ben's 2026-10-01
question).
Everything here looks good. Go for it.