---
name: diss-revise
description: End-to-end workflow for revising Bolen dissertation chapters against committee comments — from comment triage through a four-gate supervisor review. Use when asked to revise, respond to comments on, or edit Dissertation chapters.
---

# Dissertation revision workflow

Canonical specs (same folder, repo-relative to the repo root):

- `Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/style-contract.md`
- `Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/supervisor-rubric.md`
- `Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/mr-nomenclature.md`

Working set: `Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/*.tex`
(build: `ProofFinal/main.tex`; the `Notes/*.organized.tex` files are organization copies,
not the build set).

## Steps

1. **Triage.** Take one comment at a time. Quote it verbatim into the batch report first —
   Gate 4 judges the edit against the comment text, so the exact wording must be on record.
2. **Locate.** Find every passage the comment implicates (search all chapters; a comment
   about one section often implicates its echoes elsewhere — e.g., the same claim restated
   in the introduction, discussion, and conclusion).
3. **Plan the edit** against the style contract and, for math, the nomenclature mapping
   table. Decide: rewrite / add analysis / delete / add citation / decline (with reason).
4. **Edit.** Minimum churn that answers the comment. Touching a paragraph subjects the
   whole paragraph to the style contract.
5. **Numbers.** Any number in changed text gets a source line in the report
   (file/tab/row/page) — Gate 3 will re-check a sample of every one.
6. **Self-check** Gate 2 regexes over the diff hunks.
7. **Supervise.** Hand the batch to the diss-supervisor agent for the four-gate review
   (plagiarism/self-plagiarism incl. bolen_2026 and bolen_isometric_2025 adaptation rules,
   slop patterns, factual drift, comment-intent). BLOCK verdicts go back to step 3;
   scope questions and DECLINED calls go to Ben, not to reinterpretation.
8. **Record.** Batch report: comment verbatim, edits with file:line before → after,
   number sources, gate verdicts, dispositions.

## Standing rules

- Never paste comment/reviewer text into the dissertation itself.
- `\label{}`s and equation numbers are preserved by math edits (nomenclature spec §4).
- The publication note at `40-discussion.tex:5` must remain true after every batch.
- Recompile through the existing build path after math-touching batches and check the log.
