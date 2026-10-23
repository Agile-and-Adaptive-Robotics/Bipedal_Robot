# Round-4 briefing (for the next session of this revision line)

Status: round 2 is landed (REPORT.md in this folder). Round 3 is being produced
by ChatGPT on Ben's instruction. Round 4 is ZCode's pass over the round-3
output. This file carries the checklist and the verified facts so round 4 does
not need this session's context.

## Protocol for round 4

1. Load the `dissertation` skill FIRST (AGENTS.md rule for dissertation
   sessions), then read this file, REPORT.md, CUT_MAP.md, and the round-3
   output wherever Ben puts it.
2. Ben requires the **plus/minus scoring system**. It lives in
   `Documentation/Reports and Papers/Dissertation/summary.md` — ground rule 3
   and the "Machine-teaching ledger" section. The rules: every entry written or
   revised carries a score marker (**add score** / **minus score**) plus, where
   possible, the specific words that earned it — the labeled example is the
   part that teaches; numeric totals are bookkeeping. Praise earns add score;
   tracked deletions and critical comments earn minus score; unchanged regions
   earn a slight add score, discounted for read depth. Read the ledger and the
   9-item learning log in that file before scoring round 4. Ben's 2026-10-22
   clarification: NO add score for regions without comments (overall section,
   conclusion, appendices) — absence of comments is not evidence of review.
3. Run the diss-supervisor gates on the round-3 diff (they caught real
   failures in round 1). Self-plagiarism gating is REMOVED per Ben's
   2026-10-22 ruling (style-contract.md §2 and §6 carry it).
4. The pass is: relocate supporting results, delete unnecessary statements,
   narrow unsupported claims. NOT a general rephrasing sweep.
5. Appendix moves must preserve provenance: every relocated number keeps its
   table/figure/label, cross-references from the main chapter point to the
   appendix, and the measurement report gains a plus/minus accounting of what
   moved where (main-chapter words vs. appendix words reported separately —
   moving is not cutting, and Ben's instruction says to report it as such).
6. End with a clean verification build (build_one.bat in this folder; MiKTeX
   full sequence, fresh build dir) and re-measure with measure_words.py.

## The five round-3 problems (verify ChatGPT's round 3 actually did them)

1. **Appendix relocation (the standing instruction round 2 missed).** The SNS
   main chapter should keep a few representative demonstrations that establish
   the tools' capabilities, and MOVE the wider comparisons and failure
   catalogues to an appendix (revised 93-AppendixB proposal expected; round 2
   produced none). Candidates named by Ben: the gait-reference sweep
   (sec:gait_library_results + tab:gait_sweep), the pathway-ablation table
   (tab:walker_ablations + its two paragraphs), the multiple spiking
   comparisons (Simulink twins detail, AnimatLab conversion catalogue +
   tab:spike_animatlab), the conversion-failure material, and the detailed
   tuning discussion. References to appendix material do NOT count as moving.
   Preserve labels that are externally referenced (check with check_labels.py
   across ALL chapters; sec:gait_library_results and sec:walker_variant_results
   were externally referenced as of round 2 — either keep pointers in the main
   text or update the referrers).
2. **Overclaim fix (60-sns, sec:disc_spiking first paragraph, round-2 line ~370).**
   "What spiking does not change is the behavior of the controlled model" is
   too strong — the next paragraph documents timing changes and different gait
   behavior. Defensible claim: keeping the motor-neuron layer graded preserves
   the continuous activation interface; it does not guarantee unchanged
   movement. Also fix the arithmetic phrasing in the same section: 0.399 s to
   0.621 s is a 55.6% LONGER period (equivalently ~35.6% lower frequency);
   "runs 55.6% slower" reads as a 55.6% frequency reduction and is wrong.
3. **Huh?-candidate narrowing/deletion (do not always rewrite — some should go).**
   - Gap 3 in sec:sns_gaps: "no graphical circuit construction on maintained
     simulation tools" is broader than the evidence; the defensible claim is
     that the program's existing workflow lacked the needed integration.
   - 50-database: "The database answers queries that no bibliography can" —
     drop the universal ("answers queries a bibliography cannot" or just show
     the queries).
   - The discarded-bibliometric-viewer history in the SADb opening (the
     parenthetical about the off-the-shelf viewer retired within weeks) —
     Ben says failed iterations like this can simply go unsaid. DELETE, don't
     rewrite.
4. **App vs. database, stated plainly (50-database).** The five linked tables
   are the Airtable base. The app currently distributes: all 943 paper records
   (561 retained + 382 archived), the citation graph, topic clusters, and
   derived knowledge-base pages. It does NOT currently carry the non-paper
   record sets (Feedback/Models/Review Papers/Recorders) — that is the open
   APP_1393_TASK.md side task. Prose must not describe planned app capability
   as completed. Numbers: pre-prune base 1,393 records = 943 papers + 450
   non-paper; PDFs are attachments (318 flagged corpus-wide), never in the app.
5. **Patch README correctness — FIXED 2026-10-22** (application from
   ProofFinal/, baseline pre-flight check, whole-chapter patches as review
   diffs). Verify round 3 didn't reintroduce stale instructions.

## Fixed on 2026-10-22 after the round-2 review (already done; don't redo)

- REPORT.md SADb number corrected to 2,085 narrative words (the 2,063
  predated the 1,393-record sentence; SNS confirmed 6,794).
- patches/README.md rewritten: sub-patches apply from ProofFinal/ (headers
  target chapters/...), baseline pre-flight check added, whole-chapter patches
  documented as review diffs (copy from proposals/ instead of git apply).
- style-contract.md §6 reconciled with the 2026-10-22 reuse ruling (former
  Gate 1 adaptation rules superseded; cite own papers wherever one exists).

## Standing constraints (unchanged)

No commits/pushes (Ben commits), nothing uploaded to Overleaf, no comment
resolutions or tracked-edit accept/rejects, bibliography bytes untouched,
no silent number changes, protected Chapter 3 figures untouched, repo is
public (no credentials anywhere), "model" not "plant", Oxford commas.

Also standing from summary.md ground rules and learning log: no em-dashes in
new prose (split into two sentences); a section's last sentence bridges to the
next and states the dissertation's stake; completed-work framing only (never
"planned/future work" phrasing); quantitative over qualitative; problem and
solution in the same entry, never making Ben scroll to match them.
