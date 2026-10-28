## Measured test dimensions and target comparison (2026-10-28)

Latest measured-dimensions/target clarification: MEASURED_DIMENSIONS_AND_TARGET_20261028.md. Table C7 now records individual resting/tendon lengths: flexors 561/28 and 560/30 mm; extensors 828/92 and 826/85 mm. Shared prediction inputs remain separate. Both manuscripts explicitly state that the -86.5-degree measurement exceeds the human target; the 2.06 N m shortfall is against the pressure-specific prediction only. All nine extensor measurements retain their meet/exceed conclusion. No data, coefficients, fit values, or curves changed. Local dissertation 321 pages/live 322; paper 33 + supplement 10. No commits/pushes.

## Appendix C corrections (2026-10-28)

Latest Appendix C pass: APPENDIX_C_CORRECTIONS_20261028.md. Registered bone drawings, G travel/force distinction, biomimetic evaluator tables, optimizer/physical bounds, Training/Validation wording, and shared reverse-pulley schematic integrated. Final-test table explicitly lists model inputs; actual per-actuator measurements and the noted model/prose metadata differences remain unconfirmed. Local dissertation 321 pages/live 322 (0 errors/0 warnings); paper 33 + supplement 10, no errors or undefined references, existing paper font warnings. Results and coefficients unchanged. No commits/pushes.

## SNS table rewrite (2026-10-28)

Table 6.1 now uses the earlier seven-item contribution list in plain Model or tool / What it is used for columns. Former curriculum Figure 6.4 is Table 6.2. Knee/arm figures renumber automatically to 6.4/6.5. Live 308 pages / local 307, zero errors or undefined references. See committee revision SNS_TABLES_20261028.md.

## Individual pulley example retired by Ben (2026-10-28)

Paper and dissertation transition directly to the Bolen (2019) whole-leg comparison table. Individual vastus sizing/example removed from active prose and archived; no longer pending for submission. Table/results unchanged. Both live builds have zero errors or undefined references. See committee revision WHOLELEG_TRANSITION_20261028.md.

## Knee comparison interpretation corrected (2026-10-28)

Figure 6.5(A) is a startup response under constant 0.5 N m load and sustained drive, not a transient reflex/recovery test. Title, caption, and interpretation corrected; traces and values unchanged. Details in committee revision KNEE_STARTUP_FRAMING_20261028.md.

## Pipeline update (2026-10-28)

Full pipeline moved to Figure 7.4; new Figure 6.2 covers implemented construction and neural-controller tools. Uploaded and verified: 309 live / 308 local pages, zero errors or undefined references. Origin-only vastus instruction corrected. See committee revision PIPELINE_REVISION_20261028.md.

## Reflex figure accuracy update (2026-10-28)

Figures 6.5/6.6 now show all six implemented connections; Figure 6.6 panel letters and caption aligned. Biceps stretch excites biceps and inhibits triceps; Ib remains autogenic inhibition. Results unchanged. Uploaded to Overleaf: 308 pages, zero errors; local 307. See REFLEX_FIGURE_ACCURACY_20261028.md in committee revision notes.

## Latest Overleaf publication update (2026-10-28)

Ben authorized live upload. Current revised dissertation: 307 pages locally, 309 on Overleaf; torque paper: 32 pages plus 10-page supplement. Live builds have zero errors and no undefined references. Recovered W2L feedback comparison and perturbation figures now replace the earlier main comparison; Di Russo rejected adaptation reported in Appendix D. Repeated support caveats removed. Protected muscle chapter and biomimetic results frozen. Details: Documentation/Reports and Papers/Dissertation/Notes/committee_revision_20261027/OVERLEAF_UPLOAD_20261028.md.

## SADb/SNS narrative revision — October 28, 2026

Current local dissertation delivery: 304 pages, zero build errors or unresolved references. SADb leads with the 1,393-record literature collection and identifies the Papers category separately for scoped app statistics. SNS now connects the AnimatLab/MATLAB attempts, knowledge gap, conversion repairs, block/interface tools, alternative models, curriculum, and representative functional demonstrations. Overall discussion and conclusion transitions were revised. No torque/biomimetic results, torque graphics, or published Chapter 20 changes in this pass. Six chapter/appendix sources remain in ProofFinal; the unchanged curriculum figure now belongs in Figures/60-sns. See Notes/committee_revision_20261027/SADB_SNS_NARRATIVE_20261028.md and Archive/narrative_tools_20261028 for line evidence, source-count provenance, hashes, and visual checks. Root PDF/ZIP updated; no upload, commit, or push.

## Local publication pass — October 28, 2026

Current local delivery: dissertation 302 pages; torque paper 32 pages plus 10-page supplement. Builds have no errors or unresolved references. Published Chapter 20 and its supplied graphics remain byte-identical to the authoritative download. All 28 new paper comments and 25 tracked-change entries were reviewed; 27 comments are addressed locally, while the verified vastus-medialis attachment-placement example and reverse-pulley schematic await ZCode. Live comments remain open and tracked changes remain undispositioned; no upload occurred. Numerical additions, the three corrected dissertation prose GoF values, changed Figure 10 curves, and gold flexor pressure-specific predictions are explicitly disclosed in Notes/committee_revision_20261027/PUBLICATION_VERIFICATION_20261028.md and Archive/publication_pass_20261028. Root delivery PDFs/ZIPs are updated; originals are preserved. No commits/pushes.

# Reconciliation checkpoint — October 28, 2026

Current dissertation: 299 pages, zero live Overleaf errors or warnings. Root PDF/ZIP and ProofFinal sources are updated; build files are in Archive/builds/reconciliation_20261028. The delivery PDF is a local compilation of the synchronized source. Published Chapter 20 is preserved. The revision-folder and historical CHATGPT_REPORT reconciliation is complete for the recovered directives; C31/C32 remain unidentifiable question-mark comments without anchors. Exact historical saved-run evaluator provenance remains unestablished. See committee_revision_20261027/REVIEW_REPORT.md and NARRATIVE_REVISION_20261028.md for the completed work, evidence, and limits. Earlier page counts and unfinished-handoff descriptions below are historical.

Torque paper: optimized pinned-flexor Figure 8 restored; full extensor panels in Figure 9; physical validation combined in Figure 12; biomimetic results in main Table 2; former Table 1 moved to Supplementary Table S1. Abstract bracket projection restored in both documents. No numerical results changed in this narrative pass. Current count is 12 figures plus four tables; no further result removal is authorized solely to force the 15-item limit.

# DISS_STATUS.md — living dissertation status table

Maintained by /diss and /digest runs. Rows updated in place; never rewrite Ben's
prose. Last full sweep: 2026-10-27 (digest_2026-10-27.md).

## Current committee copy (supersedes historical status below)

Deadline: **October 28, 2026, 9 a.m. Pacific**, advisor/committee review.

- Working source: ProofFinal root support files and chapters/*.tex only; figures in sibling Figures. Root Bolen_Dissertation.pdf/.zip are the revised delivery, based on Ben's supplied current ZIP.
- Build: 306 pages, five successful passes, no unresolved references/citations or errors; one 0.669 pt overfull line. Build logs preserved in Archive/ProofFinal/build_sadb_final.
- Requested appendix/figure revisions installed. Thelen equations 1–7, current SNS evidence/connectomes, full five-flexor/nine-extensor before/after tables, shorter captions, current knee illustration, and bone-backed routes included. Table 4.3 unchanged.
- History audit: supplied ZIP matches Git commit 1cb0304b; 19 current live source files matched the previous delivery. Removed rejected graphical abstract, combined old Figures 4.10/4.11, exported current twelve-block SNS library, and shortened 26 additional captions. See the review report for the limits of historical verification.
- Chapter 50 now uses the cropped sadb_bubble.pdf with six callouts and expanded app functionality, verified against the supplied local HTML. Pivot scope corrected to all 943 snapshot papers. The two BPA schematic PDFs are directly in Figures/70-overall and remain Figure 7.2.
- Published Chapter 20 restored byte-for-byte to the supplied ZIP; original captions and prose retained. Its five figures use eight active graphic files, all now directly in Figures/20-muscles and byte-identical to the journal final-proof assets.
- Figures now use each chapter's own folder, including Appendices D/E, with Components and Deprecated. The old inactive KneeReflex_circuit.pdf is now in Deprecated; no file-lock cleanup remains.
- Older source trees and build/upload staging moved to one Archive. Notes contains review records and evidence, not alternate working source copies.
- Remaining: exact saved-run mechanics provenance; the proposed measurement-uncertainty task was rejected and removed. Bounded force-vector audit is complete and has sub-0.001 N m torque effects for the tested pinned configurations. No broad refit required by that check alone.
- Read committee_revision_20261027/REVIEW_REPORT.md and EASTEREGG2_HANDOFF.md. Earlier status tables and handoffs below are historical, not the current submission checklist.

## 2026-10-27 addendum (submission-framed; supersedes the 10-15 rows' currency)

**Ben's correction (2026-10-27): the dissertation is NOT submitted. Polish and
submit by 2026-10-28, 10:00 am Pacific (13:00 Eastern).** Defense stays 11-06.

- Build of record: `ProofFinal\build_r12_torque_review_20261027` (320 pp, 0
  errors, 0 undefined refs, 30 overfull). `Bolen_Dissertation.pdf` (12.4 MB) +
  `.zip` packaged 10-27 14:40 from the local tree.
- **TOP RISK — tree divergence:** the 10-26 A00–A23 advisor-facing rewrites live
  ONLY on the Overleaf tree (268 pp); the torque-physical-audit patches P01/P02/
  P11/P12 live ONLY in local ProofFinal. No single PDF carries both. Reconcile
  before submitting (`LIVE_ASTRA_20261026\after\*_LIVE.tex` are the exact
  post-edit sources).
- Unapplied diss-side: P08 (Ext20mm physical fig + wording), P09 (FlxGrp 2×2),
  E6 (finish "Then we don't need to mention…"), Appendix C.8 stale 40.2 N·m @
  804 N clause vs 27.23 N·m record of record.
- 0 placeholders / FLAG-BEN / TODO markers in ProofFinal. Flx20mm pickup
  satisfied (30-torque.tex:577).
- Post-addendum pull (10-27 evening, `76e70659` + merge): the AH001
  robot-balance wave from easteregg2 (55 files, `Code\MuJoCo_SNS\robot_torso\`)
  — 20 s supported crouched stance VERIFIED at 94.68% rig support, bit-identical
  backends, locomotion still unverified — plus the robot anatomical-frame ruling
  (first bullet in AGENTS.md Standing conventions). **No dissertation-source
  files in the pull**; the CHATGPT_REPORT.md delta is AH001 status text, not
  advisor-thread changes, so the ledger counts above stand. The AGENTS.md pull
  conflict with this session's deadline correction is resolved — both edits
  live. Defense-relevant only; does not affect the submission punch list.
- Full plan + disposition table: digest_2026-10-27.md (f).

| Chapter / section | Status | Blocking item | Last touched |
|---|---|---|---|
| main.tex (+graphicspath) | Stable; `\graphicspath{{../Figures/}}` single entry | Canonical build of record FAILED (main.pdf write-lock); clean rebuild outstanding — CAMPAIGN_REPORT §6 | 2026-10-16 |
| 10-introduction (Ch.1) | Revision campaign ACCEPTED (round 0; C01, C02); Organization rewritten to new structure | Ben: reword ~27–40-word bolen_2026 reuse or extend publication note; afferent-generalization wording (CAMPAIGN_REPORT §8.2) | 2026-10-15 |
| 15-background (Ch.2) | Revision campaign ACCEPTED (round 1; 15 ids); all changed-text numbers verified | 5 bib keys pending steward merge (needed_bib_keys.md); C10 figure picks; "faster than real time" claim contradicts 60-sns.tex:71 (§8.2) | 2026-10-15 |
| 20-muscles (Ch.3; was 20-methods, split) | Revision campaign ACCEPTED (round 0; C28, T10); mm-unit caption fix verified | Ben ruling: adaptation pass for sec:disc_* verbatim runs; text-vs-table RMSE preserved = matches published source (§8.2) | 2026-10-15 |
| 30-torque (Ch.4; was 30-results) | Revision campaign ACCEPTED (round 1; 17 ids); round-1 cam FAIL resolved against design-of-record mat | 2 minor precision flags: cam surrogate budget (default 1,000 vs winning 2,000), seeding clause flexor-only (§8.2) | 2026-10-15 |
| 40-ankle (Ch.5; new) | Revision campaign ACCEPTED (round 1; C12, C13, C14, C44, T06); no-claimed-values posture restored | Ben's standing confirmation of checkpoint-10 testbed claims (sole unattested source, ideally video locations); peterka_2002 bib key (§8) | 2026-10-15 |
| 50-database (Ch.6; new) | Revision campaign ACCEPTED (round 1; no comment ids claimed) | Deliberate disposition: ch.6 vs ch.8 coverage-stats coexistence + triple-"95" (§8.3) | 2026-10-15 |
| 60-sns (Ch.7) | Revision campaign ACCEPTED (round 1; 16 ids); 410-neuron, block-census 12/15, 2 ms-step fixes verified | Prior-walker placement (checkpoint 3a) confirmation + ledger annotation; minor wording notes (§8.2/§8.3) | 2026-10-15 |
| 70-overall (Ch.8; new whole-work chapter) | Revision campaign ACCEPTED (round 0; T11); retitle "The Work as a Whole and Future Work" applied | C12/R20 claims unresolvable vs ledger — ratify retitle or take offered revert; sec:methods_overview removal ratification (§8.1) | 2026-10-15 |
| 80-conclusion + 92/93/94 appendices (back matter) | Revision campaign ACCEPTED (round 0; C46, T12–T15) | T14 disposition (target sentence still present in both copies, lines 205/223) | 2026-10-15 |
| Figures tree | 52/52 targets exist + cam-winner/route-check/GIF candidates staged; QC verdicts written (1 PASS, 6 FAIL) | Ben's fix decisions per figure_qc.md — no figure was edited/regenerated (CAMPAIGN_REPORT §5) | 2026-10-16 |
| Overleaf zip/pdf | Synced 2026-10-15 (zip-is-canonical ports landed in ProofFinal); upload bundle READY | Upload ONLY on Ben's go (latex-overleaf skill); T01–T15 pending in review UI — re-run sync after any reject (§10) | 2026-10-16 |
| Defense deck | pptx untouched since 9-24; slideshow GIFs/PNGs landing | Wire 15+2 new assets into deck; visual QC | 2026-10-15 |

Advisor comment ledger: `CHATGPT_REPORT.md` (2026-10-15): 46 open threads /
50 messages; 15 pending tracked changes T01–T15; 32 resolved R01–R32;
commenters ajh26 ×63, mcneal ×11, bbolen83 ×4. Restructure directive = C12.
Reward signal of record: "things are looking pretty good" (ajh26, 2026-10-15).

Revision campaign 2026-10-15 (closed 2026-10-16): all 9 chapter slices accepted
(46 C-ids + 12 T-ids claimed addressed; C12/R20 in the whole-work slice unresolvable vs
ledger); reward 108 (+1 baseline = the ajh26 signal above; weights per the ledger).
Full record — verdicts, checks run/not run, figure QC, Ben queue, Overleaf gate:
[CAMPAIGN_REPORT.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/CAMPAIGN_REPORT.md>).
