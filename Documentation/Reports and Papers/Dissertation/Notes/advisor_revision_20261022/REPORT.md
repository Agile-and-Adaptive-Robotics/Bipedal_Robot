# Advisor-review revision session — 2026-10-22 (EB475WS4), round 2

Round 1 delivered the restructure; the principal review ruled it "would not apply
as-is": verified content was treated as must-stay, contributions were repeated
instead of selected, SADb was product-first in heading only, two contradictions
survived, and the delivery was too broad. Round 2 (this document's state) applies
selection, compression, and explanation — not another expansion. Everything is
LOCAL: nothing uploaded to Overleaf, no commits or pushes, no comments resolved,
no tracked edits accepted or rejected.

## Final measurements

| Chapter | Baseline narrative | Final narrative | Δ narrative | Captions | Tables |
|---|---:|---:|---:|---:|---:|
| 10-introduction | 2,093 | 2,227 | +134 | 84→84 | 0 |
| 15-background | 5,217 | 5,387 | +170 | 529→500 | 46→46 |
| 30-torque | 8,993 | 9,154 | +161 | 1,517→1,517 | 459→459 |
| 50-database | 2,588 | 2,085 | **−19.4%** | 338→338 | 519→519 |
| 60-sns | 10,851 | 6,794 | **−37.4%** | 1,069→1,096 | 316→463 |
| 70/80/94 | — | person sweep only | — | — | — |

SNS total reading matter (narrative + captions + table bodies):
12,236 → 8,353 (**−31.7%**) — reported this way because the round-1 headline
narrative number hid the caption/table burden the review flagged.
SADb final = 2,085 narrative words as measured against the final files
(the earlier 2,063 predated the added 1,393-record sentence).
PDF: 294 → **280 pages**. Verification build of record: `tmp/build_live8.log`
— 0 errors, 0 undefined references/citations, 0 missing figures; overfull-hbox
multiset vs. baseline: 0 added, 2 removed. Capability-table page fits the text
block exactly (ink-bounds check, `renders/v4_p164.png`).

## What round 2 changed

1. **SNS: the capability table now REPLACES the tool prose.** The four-family
   paragraph is a one-sentence pointer; per-tool subsections carry only what the
   table cannot (constructions and their relevant limits, stated precisely).
   Conversion narration, secondary demonstration numbers, unit-level spiking
   residuals, per-perturbation scores, and duplicate status accounts were cut
   (itemized in `CUT_MAP.md`). `tab:spike_beer` deleted (numbers preserved in one
   sentence). The closing section went from six subsections to three; the
   tools-summary subsection that repeated the table is gone; the walking status is
   stated once. Kept per the review: the Deng–Li opening, the completed/proposed
   distinction, the Hunt/Young lineage paragraph.
2. **SNS contradiction fixed (the flagged scoring issue).** The curriculum
   paragraph claimed the −25 fallback "sits below every genuine gait"; it sits
   numerically ABOVE every genuine gait score (−25 > −144…−267), which is exactly
   why it is a sentinel excluded from score comparisons. The text now says so and
   matches the table captions.
3. **SADb round 2.** Database-first opening (SADb is sentence 1). Import-batch,
   recovery, and per-stage prune counts cut; decision authority stated once in
   polished prose ("I held every decision that required judgment or carried
   risk") instead of "approved by me" stage annotations. DOI contradiction fixed:
   "unique identifier wherever a record carries one" (550 of 561 carry DOIs).
4. **App-first ruling applied.** The packaged explorer is the database's primary
   form; Airtable is the editable working base. Verified before writing: the app
   embeds all 943 audited papers, including the 382 pruned from Airtable, each
   under an archived flag (`export/sadb_export.json`, `app/check_app.py`: "papers
   embedded: 943"). `SADb_audit/app/README.md` roadmap line updated to match.
5. **First person everywhere.** "the author" swept to I/my/me in chapters 1, 2,
   5, 6 (proposals) and in new minimal patches for 70-overall, 80-conclusion,
   94-AppendixC (chapters 3–4 contain none). `style-contract.md` §2 superseded
   with the 2026-10-22 rulings (first person; self-plagiarism gating removed;
   own papers cited wherever one exists).
6. **bolen_design_2020 citation reverted** from 60-sns per Ben: it is a poster
   schematic, not a paper, and its connectome is not the current one.
7. **Self-plagiarism gating removed everywhere** (Ben's ruling). The Gate-1
   rewording from round 1 stays as prose (it reads well); the gating itself no
   longer constrains edits, recorded in the style contract.
8. **Torque patch kept as "candidate pending provenance checks"** per the review.

## Delivery format (round 2)

- **Discrete, independently applicable sub-patches** for 10, 15, 30 under
  `patches/<chapter>/` — each verified to apply cleanly to the pristine baseline
  (`git apply --unidiff-zero`), grouped by fix family; applying all of a chapter's
  sub-patches reproduces the proposal exactly. See `patches/README.md`.
- **Whole-chapter patches** only where the advisor-mandated restructure cannot be
  split: 50-database, 60-sns (plus the three minimal person-sweep patches for
  70/80/94). Full-chapter context patches for every chapter remain alongside.
- `proposals/` still carries the drop-in files for reference.

## Open decisions for Ben

1. **BG-24** — the advisor's "Delete?" targets this sentence in the background
   (master's-study paragraph): "PAM routings were then placed by hand in OpenSim
   beside the Gait2392 paths they imitate, moved to respect Festo BPA contraction
   limits and to keep actuators from interfering with one another, and the result
   was saved as the gait2392_robotbody model that the later chapters of this
   dissertation continue to refine." The proposal keeps it, reframed. Delete or
   keep is yours. (BG numbers are comment-thread IDs in CHATGPT_REPORT.md, in
   review order — not line numbers.)
2. **SNS-14** (advisor: "that is all you are saying about the other tools?") vs.
   the reduce directive: current answer is the Hunt/Young lineage paragraph plus
   precisely-stated limits. Expand or answer in the reply thread.
3. **04-remove-graphical-abstract.patch** (BG-04) still needs your explicit
   sign-off (it deletes a figure).
4. **Figure-level comments** outside .tex scope: BG-03 (motor drawing position),
   SNS-11 (overlapping figure text) — graphics pass needed.

## Provenance and constraints

- One recorded correction (round 1, unchanged): the first three torque points were
  on the earlier direct route, not the rebuilt route. Manuscript staleness flags
  (Main.tex lines 556/560; line 377 `\newton\cdot\meter` inconsistency in live
  text) still stand. No numbers were changed silently; every round-2 removal is
  itemized in `CUT_MAP.md`. Bibliography bytes untouched.
- Throwaway build trees (`ProofFinal/livebuild_20261022/`, `livebase_20261022/`,
  `D:\temp\diss_build_{live,base}_20261022/`) are safe to delete. `tmp/` scratch
  is disposable.

## Files to commit

The whole folder `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/` (proposals, patches incl.
sub-patch folders and README, baseline anchor, CUT_MAP, SNS_MOVE_MAP,
COMPLIANCE_TABLE, renders, scripts, this report), plus the two edited ruling
documents outside it: `Notes/revision_campaign_20261015/style-contract.md` (§2
supersession) and `SADb_audit/app/README.md` (app-primary roadmap line).
