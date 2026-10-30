# Assistant review rubric and failure costs

30 October 2026. This rubric evaluates the assistant's editing/review failures. It does not grade Ben, assign a PhD score, score his AI-tool comparison, or treat a reviewer's tone as evidence of scientific misconduct. It is a persistent error-analysis/checklist system; no model weights were retrained.

## Weighted penalties

The objective is zero avoidable misses. Report the number of confirmed unique failures and their weighted cost, not a spurious percentage probability of acceptance. Costs can exceed 100 across a campaign. Count one underlying miss once at its highest applicable cost; record associated comments as evidence rather than adding duplicate penalties. Do not count a comment as unresolved merely because its review thread is still open.

| Cost | Assistant failure | Trigger and required response |
|---:|---|---|
| 40 | Wording/evidence failure that needlessly exposes the author to an integrity challenge | Verify the actual trigger. Remove unnecessary self-defense, ungrounded claims, or contradictory selection rules; preserve accurate methods and the author's reply. A hostile comment alone proves no scientific fault. |
| 30 | Fabricated source, result, completion claim, or unauthorized scientific change | Stop clearance; recover original evidence and disclose the exact error. No silent harmonization. |
| 25 | Misreading or overriding Ben's directive; presenting an assistant inference as an advisor requirement | Record instruction provenance verbatim. Withdraw the invented task; do not replace it with another unsolicited task. |
| 20 | Missed substantive tracked edit, resurrected deletion, or failure to read an author's plain-text changes | Capture author, changed text, anchor, and date; reconcile against the effective source. Metadata absence is not evidence of absent edits. |
| 12 | Unanswered 'huh?', 'which?', 'what do you mean?' | Locate the actual anchor and repair the concept/referent. If the anchor is uncertain, say so instead of declaring success on a nearby paragraph. |
| 10 | Figure/table/claim mismatch, illegible evidence, wrong caption, or broken cross-reference | Inspect the final compiled size and content; check units, letters, denominators, and consistency with the claim. |
| 8 | A substantive advisor tracked replacement that an adequate editing pass should have caught | Treat as a labeled example of a recurring reader problem; learn the principle, not merely the replacement string. Mere use of track changes is not an automatic defect. |
| 6 | Repeated result inventories, generic caveats, operational walkthroughs, or long unfocused paragraphs | Prefer deletion, consolidation, or an existing display. Do not add a decorative figure to satisfy a quota. |
| 2 | Copyediting-only tracked correction | Record the recurring grammar/notation rule and apply it consistently without broad rewriting. |

## Historical calibration and regression rules

| Evidence | Assistant failure / cost | Prevention rule |
|---|---|---|
| Chapter 5 comment5.13, 29 Oct09:40; Ben reply 14:57 | Integrity-inviting curation narration, category40. The original wording invited a gratuitous inference about the author rather than describing selection. | State the true active curation rules; preserve Ben's reply; check cross-section denominators/status definitions. No accusation is inferred from tone. |
| Ben's saved Word anchors, 29 Oct16:14 | Treating absent Word comments/tracking as absent new text, category20. | Read the complete current document text and compare content, including headings and Disposition fields. |
| Advisor tracked deletions 08:47/08:49, Chapter 4 | Missing effective deletions in a source capture, category20. | Read tracked edits separately from comments; check that deleted passages stay deleted in effective output. |
| Chapter 5 comment5.23 screenshot | Calling the earlier polish complete while the criticized screenshot remained, category10. | Verify the included asset itself at final size, not just revised prose. |
| This audit's initial 'documented assessment method' recommendation, corrected by Ben 30 Oct | Misattributing an inferred task to the advisor, category25. | Exclude the AI-comparison issue entirely. It was Ben's judgment of accuracy and speed relative to his own work. No assessment-method, attribution-check, benchmark, or wording task remains. |

These five documented calibration events have weighted cost 115, counted once each. This is a historical sample used to set priorities, not a score for this manuscript or an exhaustive failure count. Current live tracked-change totals were not obtained and are not invented.

## Acceptance checks before saying an edit is done

1. Identify authority: newest supplied ZIP, advisor thread/anchor, Ben reply, Ben directive, or new recommendation. Keep these distinct.
2. Read the effective source and actual included figure. A saved comment or old report is not current-state proof.
3. Trace each changed number, equation, model identity, and claim to the existing record. Resolve contradictions before stylistic rewrites.
4. State the engineering result before its display when appropriate; define the target and the conditions under which it holds.
5. Assign every paragraph one purpose. Remove duplicate recap and operational detail before adding graphics.
6. Check final-size axes, units, legends, letters, captions, table denominators, and first callouts. Equations belong in text unless naturally part of a schematic.
7. Compare article/dissertation/supplement carryovers explicitly. Never assume a change to one changed the others.
8. Build the reviewed snapshot and inspect affected pages. A clean build is necessary, but does not clear scientific or layout findings.
9. Report addressed / partial / pending / deferred-by-Ben separately, with evidence. Do not accept tracked changes, resolve threads, or upload during a findings-only audit.

## Boundaries retained

FVU is settled by Ben's reply: no refit, metric replacement, or optimization rerun. The database classification exercise remains deferred. Published Chapter 20 is protected. The waived individual vastus/reverse-pulley worked example is not revived. True assistance conditions and negative results remain visible; administrative debugging history need not remain in the scientific narrative.

Use `SOURCE_MANIFEST.json`, the three subagent reports, and `APPENDICES_A_B_C_AUDIT.md` as the labeled examples for the next pass. Add a regression rule only when an observed failure justifies it; do not manufacture new requirements from every suggestion.
