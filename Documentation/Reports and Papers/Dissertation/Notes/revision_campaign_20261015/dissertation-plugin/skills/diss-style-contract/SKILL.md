---
name: diss-style-contract
description: Enforces Ben's dissertation writing style for any new or revised prose in Documentation/Reports and Papers/Dissertation. Use before writing or editing dissertation chapter text — voice, punctuation, spelling, sourcing, and banned AI-slop patterns.
---

# Dissertation style contract (enforcement skill)

The canonical spec lives in the repo and is the authority — read it before writing:

- `Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/style-contract.md`
  (voice, person/tense, punctuation, spelling, numbers, citations, slop bans, rhythm,
  reviewer-comment etiquette, change discipline)

Related specs in the same folder:

- `supervisor-rubric.md` — the four gates; Gate 2 holds the slop regex table.
- `mr-nomenclature.md` — required for any math edit.

## Hard rules (never violate; details in the contract)

1. Professional ME-committee academic prose. Calibrate against the exemplars cited in the
   contract (`10-introduction.tex:6-24`, `20-methods.tex:161-163`, `40-discussion.tex:9-17`).
2. Oxford comma, always. American spelling.
3. Every number has a named source (measurement file, computed result, or cited reference
   at the cited page). No source, no number. Never silently harmonize disagreeing numbers.
4. Never paste reviewer/committee text into the dissertation; no meta-commentary about
   revisions ("as requested…").
5. Banned outright: delve; "it is important to note"; "plays a crucial role"; metaphorical
   "landscape"; moreover/furthermore/additionally chains; em-dash padding (≤1 per ~500
   words, never clause-splicing); uniform paragraph rhythm. Full tiered list with regexes:
   rubric Gate 2.
6. Units through siunitx (`\qty`, `\qtylist`, `\SIrange`). No bare "10mm".
7. An edit answers its comment and stops — no drive-by rewording of untouched sentences.

## Procedure

1. Read the style contract (and `mr-nomenclature.md` if math is touched).
2. Write/revise the passage.
3. Self-check the changed paragraphs against the rules above and the Gate 2 regex table.
4. Hand the batch to the diss-supervisor agent for the four-gate review.
