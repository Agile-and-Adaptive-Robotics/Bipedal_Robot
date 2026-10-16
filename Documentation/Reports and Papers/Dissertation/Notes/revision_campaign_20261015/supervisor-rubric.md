# Supervisor rubric — four gates for the dissertation revision campaign

Campaign: post-defense revision, 2026-10-15. Paths relative to the repo root
(`D:\Github\Bipedal_Robot`). Run all four gates over every revision batch **before** the
batch is declared done. Each gate returns PASS, FLAG (with dispositions required), or FAIL
(with file:line). Overall verdict: GO only if no gate failed and every flag is dispositioned
in the report; otherwise BLOCK with the specific failures.

Scans run over **changed paragraphs only** (diff hunks against the pre-batch state), not
whole chapters — existing grandfathered text is out of scope unless the batch touches it.

---

## Gate 1 — Plagiarism and self-plagiarism

**1a. External sources.** No verbatim run of 7+ words from any source Ben does not own,
except: quoted with quotation marks and a citation, standard technical phrases (equipment
names, model names), or equation terms. Equations reproduced from a paper carry a citation;
re-set them in dissertation notation, do not screenshot.

**1b. Ben's own two papers — adaptation vs. duplication.** The dissertation's publication
note (`ProofFinal/chapters/40-discussion.tex:5`) is the authority on overlap:

- `bolen_2026` — Bolen, Elzein, Pang, and Hunt, *Actuators*, 2026 ("Isometric Force
  Characterization of Braided Pneumatic Actuators", bib entry at
  `ProofFinal/thesis.bib:1086`). **Published.** Dissertation sections
  `sec:disc_effects`–`sec:disc_bio`, with their figures and tables, "appeared in" it.
  Rules for prose taken from it: adapted, not duplicated — restructured for the
  dissertation's frame, dissertation-specific connective text, no verbatim runs of 10+
  words beyond technical necessity. Figures/tables may be reused (own work) with the
  publication note intact and captions re-checked against the dissertation's context.
- `bolen_isometric_2025` — Bolen, Burgess, Morrow, Elzein, Pang, Scharzenberger, and Hunt,
  "Isometric Torque Values About Robotic Knee Using Braided Pneumatic Actuators",
  unpublished manuscript in preparation (bib entry at `ProofFinal/thesis.bib:1113`);
  `sec:disc_error` is its subject. **Unpublished and unsubmitted.** More latitude applies,
  but the prose is still reworded so the dissertation is not a straight dump of the
  manuscript, and the note's statement ("has not been submitted") must stay true.

**Procedure.** For each changed hunk that overlaps either paper's subject matter: locate the
paper's manuscript text (Ben maintains these; do not guess a path), extract word n-grams
(7-word window for external, 10-word for Ben's own), and diff. Report every hit with
file:line on both sides plus a disposition (technical phrase / reworded / quoted).
The publication note itself is checked for continued accuracy.

**Fail conditions.** Unquoted external verbatim ≥7 words; published-paper (`bolen_2026`)
verbatim ≥10 words without disposition; publication note rendered false by the edit.

---

## Gate 2 — AI-slop patterns

Run case-insensitive regex over the changed text (LaTeX source, comments stripped).

**Tier A — hard fail, no legitimate use in Ben's voice:**

| Pattern (PCRE, /i) | Notes |
|---|---|
| `\bdelv(e|es|ed|ing)\b` | delve |
| `\bit is important to note\b` | |
| `\bit is worth noting that\b` | |
| `\bit should be noted that\b` | |
| `\bplays? (a|an) (crucial|pivotal|vital|key) role\b` | in this table the inner pipes are shown escaped; the working regex uses plain `|` |
| `\blandscape\b` | metaphorical use; literal terrain sense judged |
| `\btapestry\b` | |
| `\b(a|the) testament to\b` | |
| `\bcutting-edge\b` | |
| `\bseamless(ly)?\b` | |
| `\bmyriad\b` \| `\bplethora\b` | |

**Tier B — flag, judged in context:**

`moreover`/`furthermore`/`additionally` as sentence openers; `holistic`; `leverage` as a
verb; `underscore(s/d)`; `showcase(s/d)`; `paradigm`; `comprehensive`; `multifaceted`;
`in the realm of`; `navigate the complexities`; `not only … but also` chains.

**Structural checks (not regex):**

- **Connective chains:** ≥2 sentence-opening `Moreover,/Furthermore,/Additionally,` in one
  paragraph, or ≥3 in one revised passage → flag.
- **Em-dash padding:** count `---` in changed prose. More than one per ~500 words, or any
  habitual use splicing independent clauses → flag. En-dash `--` ranges are exempt and
  correct (`20-methods.tex:80`).
- **Uniform paragraph rhythm:** in a revised section, ≥70% of paragraphs within ±1 sentence
  of the same length AND ≥60% opening with the same pattern (e.g., "The + noun + verb") →
  flag; a human confirms it reads mechanical.
- **Reviewer-text leak:** any occurrence of committee/comment language in the dissertation
  ("as requested", "the committee suggested", quoted comment fragments) → hard fail; also
  check the diff for meta-apologies ("This section was previously unclear").

---

## Gate 3 — Factual drift

Every number in a changed paragraph is spot-checked against the cited source file, at the
scale of the whole batch.

**Procedure.** Extract every numeric token from the hunks (including siunitx arguments).
For each, record: value + unit, kind, and the exact source it was checked against:

- Measured values → `Testing_Data/` files (e.g., `Results_table_10mm_pinned.xlsx` tabs,
  `Results_table_20mm.xlsx` tabs) at the specific row/tab.
- Computed/model values → the generating script and its saved results
  (`Code/Matlab/Mesh_Optimization/Results/`, `2022_02_Festo/Dig_out/` mats, spinal reports).
- Literature values → the cited document itself; for datasheets, the cited page
  (the existing text already does page-level citation, `40-discussion.tex:13`).
- Dissertation-internal claims → the target label actually says it (open the referenced
  equation/table/section and compare).

**Fail conditions.** Any number with no named source; any number that changed without the
source changing; any silent harmonization of disagreeing numbers (style contract §5);
a value whose unit class changed in the edit.

---

## Gate 4 — Comment-intent

Did the edit actually answer the comment? Judged **against the comment text itself**, not
against the editor's summary of what was done.

**Procedure.** For each comment in the batch, produce a table row:

1. **Comment (verbatim).** Quoted exactly, including any reviewer emphasis.
2. **Claimed edits.** File:line ranges, before → after (hunk summaries).
3. **Judgment.** Read the comment as written and the post-edit text as written. Does the
   change do the specific thing asked — not a nearby simplification, not a partial that
   reads as full? A comment asking for analysis gets analysis, not restatement; one asking
   for deletion is answered by deletion, not softening; one asking for a citation gets a
   real citation, not a vaguer claim.
4. **Verdict:** ADDRESSED / PARTIALLY ADDRESSED (state what remains) / NOT ADDRESSED /
   DECLINED (only with a stated reason that Ben can overrule).

**Fail conditions.** Any comment marked ADDRESSED that is not; any comment silently
dropped from the batch report; a DECLINED without a reason.

---

## Output contract

The supervisor returns, per batch: per-gate verdicts with file:line findings, the Gate 4
table in full, and an overall **GO** or **BLOCK**. Findings quote the offending text.
The supervisor edits nothing itself; discrepant judgments escalate to Ben rather than being
labeled away.
