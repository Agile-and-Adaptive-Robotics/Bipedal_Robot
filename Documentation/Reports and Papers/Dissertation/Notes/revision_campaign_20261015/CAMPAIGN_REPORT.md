# CAMPAIGN REPORT — Dissertation Revision Campaign (revision_campaign_20261015)

Written 2026-10-16 by the campaign reporter. Sources = this campaign's own inputs: the
nine per-chapter supervisor outcomes (quoted throughout as "supervisor notes"), the
figure-QC verdicts ([figure_qc.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/figure_qc.md>)), the final build result, the reward score, and the flagged
items from the sync/restructure/audit phases. **Nothing was re-measured for this report**;
every number below is carried from those inputs, and every "verified" claim was verified
by the chapter supervisors in their sessions, not by the reporter.

Campaign artifacts referenced (all verified present by directory listing this session):
[comment_ledger.json](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/comment_ledger.json>),
[needed_bib_keys.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/needed_bib_keys.md>),
[source_sync.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/source_sync.md>),
[migration_map.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/migration_map.md>),
[restructure_log.txt](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/restructure_log.txt>),
[style-contract.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/style-contract.md>),
[supervisor-rubric.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/supervisor-rubric.md>).

---

## 1. Verdict

**All nine chapter slices were ACCEPTED (Gate 4).** Four on the first pass (rounds 0):
whole-work discussion, introduction, muscles, back matter. Five after one revision round:
ankle, database, background, SNS, torque. No slice sits in FAIL or pending state.

The campaign is **not** fully closed. Outstanding: the canonical compile of record FAILED
(§6); 24 comment anchors are inferred rather than ledger-anchored and await Ben (§7); a
queue of deliberate dispositions needs Ben's rulings (§8); the Overleaf upload is staged
but gated on Ben's go (§10).

## 2. Per-chapter comment status

| File (chapter) | Slice | Status | Rounds | Comments addressed |
|---|---|---|---|---|
| 70-overall.tex (Ch.8) | Whole-work discussion | accepted | 0 | 3 claimed: C12\*, R20\*, T11 |
| 10-introduction.tex (Ch.1) | Introduction | accepted | 0 | 2: C01, C02 |
| 20-muscles.tex (Ch.3) | Muscles chapter | accepted | 0 | 2: C28, T10 |
| 40-ankle.tex (Ch.5) | Ankle chapter | accepted | 1 | 5: C12, C13, C14, C44, T06 |
| 50-database.tex (Ch.6) | Database chapter | accepted | 1 | none claimed |
| 15-background.tex (Ch.2) | Background | accepted | 1 | 15: C03–C11, C15, C35, C36, T02, T05, T07 |
| 60-sns.tex (Ch.7) | SNS chapter | accepted | 1 | 16: C05–C08, C15–C23, C37, T03, T04 |
| 80-conclusion + 92/93/94 appendices | Back matter | accepted | 0 | 5: C46, T12, T13, T14, T15 |
| 30-torque.tex (Ch.4) | Torque chapter | accepted | 1 | 17: C24–C27, C29–C34, C38–C43, C45 |

\* **C12 and R20 as claimed by the whole-work writer could not be reconciled with the
ledger** (§8.1): the ledger's actual C12 (comment_ledger.json:113-121) is ajh26's
balance-testbeds MOVE comment homed to "ankle" — claimed and addressed there. "R20"
exists nowhere in the ledger or campaign notes. The whole-work chapter's retitle
("The Work as a Whole and Future Work") therefore rests on an unverified comment mandate
and needs Ben's ratification (the writer offered revert).

**Coverage arithmetic** (from the table above): 65 id claims across slices resolve to all
**46 C-series ids (C01–C46)** plus **12 distinct T-series ids** (T02–T07, T10–T15) and R20
(unresolvable). C05–C08 and C15 are claimed by BOTH background and SNS (multi-home
comments); C12 is claimed by both overall and ankle (only the ankle home matches the
ledger). **T01, T08, and T09 are claimed by no writer batch** — the sync notes state
T01–T15 content is already applied in the working tree (2026-10-15 zip → ProofFinal), so
those three carried via the sync, not via a writer.

## 3. Math / numeric change list — pointer

The **writers' `mathChanges` lists live in the run results of the writer batches** of this
campaign (each writer returned its math-change list alongside its outcome); they are not
reproduced in this folder and not copied into this report. What this report carries is the
**supervisors' verification** of the numbers those edits touched, condensed per chapter in
§4; the full verification detail is in the supervisor notes of each chapter's outcome
quoted in the campaign's run record.

## 4. What was verified, and what was not run (per chapter)

Condensed from the supervisor notes. "Verified" = the supervisor executed the check in
their session, as described; "NOT RUN" = stated as such by the supervisor.

**70-overall (rounds 0).** Numbers verified: 943/833/95/363/18/11,130 against
SADb_audit/README.md:587-588 + knowledge_base/INDEX.md:3 (values unchanged from
pre-campaign; in-sentence scope matches source); 11,130/18 also at README.md:569,546;
0.9999 and exact-period agreement at 60-sns.tex:142; 60% pressure ceiling at
60-sns.tex:157; 0.1 ms and 0.9 s at 60-sns.tex:73; 4 Hz / ~1 Hz against the
bolen_design_2020 course-report PDF; "three diameter-dependent coefficients" at
20-muscles.tex:318,50. Ref integrity by label grep (inbound refs resolve in 15-background,
30-torque, 40-ankle, 60-sns, 80-conclusion, 93-AppendixB, 10-introduction:30,41,47).
NOT RUN: compile (concurrent writer batches in the tree). Hazards carried to §8: the
ch.6↔ch.8 coverage-stats coexistence and triple-"95"; C12/R20 unresolvable; the writer's
`sec:methods_overview` removal (no ledger mandate; verified zero dangling references —
findstr over chapters/*.tex + main.tex returned no `sec:methods_overview` references);
stale writer flags about figure_qc.md (it exists) and the intro (now describes the NEW
structure).

**10-introduction (rounds 0).** N-gram check vs the published bolen_2026 manuscript run
(10-gram/7-gram sweep): 36/43 hits, all inside sentences identical on both sides of the
git diff; the batch's new sentences have zero overlap. All referenced labels exist; zero
dangling `sec:problem` refs document-wide; all six citation keys present in thesis.bib.
NOT RUN: n-gram check vs bolen_isometric_2025 (manuscript not in the repo); full compile
(rubric names none; verification was grep-based). Carried to §8: the pre-existing ~27–40-
word reuse from bolen_2026 (grandfathered, disposition = Ben's); the
"afferent feedback did not generalize to weight-bearing motion" wording, which asserts a
failure mode 60-sns.tex:18,20 does not state (chapter records suspended-only
demonstration, no ground-contact signals, support required).

**20-muscles (rounds 0).** No new plagiarism/self-plagiarism: no ≥7-word external runs
and no ≥10-word runs from bolen_2026 in the new paragraphs (longest overlap = 6-word
technical phrase). Unit fix verified: fig:ModelCompare caption 0.257/0.233/0.300/0.450 mm
→ 257/233/300/450 mm against tab:model_compare rows, body text, and article line 507 (a
unit-magnitude error). Cross-chapter claims verified (30-torque.tex:75, 210, 259-267;
60-sns.tex:10). The writer's claim that figure_qc.md is absent is FALSE (process error,
flagged; outcome unaffected — QC covers only torque figures). Preserved correctly (per
style-contract §5, no silent harmonization): the text-vs-table RMSE discrepancy at
20-muscles.tex:291 vs tab:Fstar — identical in the published source of record, so the
dissertation matches its source verbatim. Carried to §8: pre-existing ≥10-word verbatim
runs from bolen_2026 in mechanically-touched paragraphs (disclosed by the publication
note; "adapted, not duplicated" standard unmet for sec:disc_* prose — Ben's ruling).

**40-ankle (rounds 1).** Round-2 mitigation verified in text: the contested framing is
deleted; the section now states the recordings "exist as the author's records of his own
experiments" and that timings/set angles "have not been extracted, and no such values are
presented" (40-ankle.tex:29) — no claimed values. Independent scratch rebuild of the FULL
document by the supervisor: pdflatex ×3 + bibtex, exit 0, 250-page PDF (9,712,497 bytes),
zero errors, zero undefined references; undefined citations = exactly the 6 pending
steward-merge keys, of which peterka_2002 (page 111) is the only one from this chapter.
**UNVERIFIED-ATTESTED (carryover):** "both self-balanced and walked" (IP 2.0;
pre-campaign 30-results.tex:346 claimed balance only), "within about 5°" (no repo
artifact exists — searched), "Most testbed configurations held balance indefinitely", and
the peterka_2002 protocol-match claim — sole named source is the writer-reported Ben
checkpoint (10) of 2026-10-15, recorded nowhere in the repo. Ben's standing confirmation
(ideally with video file locations) is the outstanding disposition. thesis.bib still lacks
peterka_2002 (verified, findstr 0 hits); verified entry queued in needed_bib_keys.md.

**50-database (rounds 1).** Accepted with no comment ids claimed and no supervisor notes
carried in the outcome. Note the phase-A/phase-B history in §9 (skeleton only during
restructure; the database writer filled it in the same working tree — see the §8.3
stats-coexistence item).

**15-background (rounds 1).** Every number in the changed text checks against its source
(round-1 verifications stand: Aschenbeck greyhound details and Steele 0.079/151° from
thesis.bib abstracts; Liang 2020 Table 1 TFA 22 MPa and its ">100 times" sentence;
Kengoro 1670 mm/56.5 kg/174 DoF; 50-min film = 50:13 per transcript note); the round-2
edits REMOVED the one unverified number (the "within roughly ten percent" tolerance,
deleted) and ADDED transection-level verification (needed_bib_keys.md:35-39 quotes the
film transcript). Dependency open: 5 cited keys absent from ProofFinal/thesis.bib
(utherst_decerebrate_2025, stuart_reflections_2007, severin_1967, geijtenbeek_scone_2019,
shinohara_mechanisms_2025 — regex check of all 106 cited keys); full entries queued in
needed_bib_keys.md; severin_1967 carries the verify-against-Stuart-2007 caveat. C10
(borrowed background figures) still open in the document: no figures added; deferral
states a candidate list with licensing, which the rubric permits. T03/T04 wording
conflict carried to §8.2.

**60-sns (rounds 1).** Round-2 fixes verified FIXED: (1) "406-neuron" → "a 410-neuron
two-level spinal network" (60-sns.tex:71), matching DESIGN.md:1928/1932; (2) library
block counts rewritten and verified by block census of the build scripts (sns_build_library.m
creates 8 = 7 core + SynSum; +4 actuator blocks = 12; +3 spiking blocks = 15; chapter
arithmetic 7+1+3+1=12 and 12+3=15 internally consistent, matches
goal4_dissertation_spiking.md:118; the ">6 synapses" claim verified at
sns_build_library.m:20 and :207); (3) "0.1 ms neural step" → "a 2 ms neural step, matching
the physics-engine timestep", verified against runner.py (one net.step :1535 + one
mj_step :1595 per loop iteration :1243; nsteps :1074; MJCF timestep 0.002 implicitfast
:250). Round-1 Gate-1 duplication RESOLVED: intro compressed to pointer paragraphs; n-gram
overlap vs 15-background re-run = 0 seven-gram / 0 ten-gram (was 60/40). Gate-2 slop scan
re-run over ~4,499 added words: Tier A NONE, Tier B NONE. Compile check: pdflatex ×2 from
ProofFinal → output written, 0 errors, 0 undefined references, 0 multiply-defined labels
(check artifacts deleted). Still open (carried to §8): background:183 "faster than real
time" contradiction with 60-sns.tex:71 ("Full-body closed-loop co-simulation is not yet
running"); prior-walker placement rests on the writer-attested "Ben checkpoint 3a
override" (ledger C15 still says newHome=background); minor notes (40 mm element's Fmax
is a fixed class constant, so "implementing the fitted force surfaces" is slightly
generous; stale fprintf block-counts in the build scripts; "air-stepping gate" wording;
three grandfathered "% ALT TEXT:" comments).

**Back matter (rounds 0).** Accepted; C46, T12–T15 claimed; no supervisor notes carried
in the outcome. T14's live-file mismatch is in §9 (sync flags).

**30-torque (rounds 1).** Round-1 FAIL RESOLVED: 30-torque.tex:433 now quotes the
design-of-record mat's values (+30.83% vs required, +39.77% vs plain, 269–507 kPa, 96.9 mm
bend radius, moment arm 59.4–102.8 mm mean 72.4, clearance 43.1 mm) — each matching the
supervisor's MATLAB recompute of Vas_Pam_20mm_Result_cam_20261014_0416.mat — and
attributes +31.30%/+42.34% explicitly to the pre-adjustment run-8 winner
(tmp\optcam_run8_20261013.log). "7,000 evaluations" replaced by per-configuration budgets
(extensor 3,000 = Opt_run_Ext.m:29 + mat funccount; flexor 1,000 = Opt_run.m:79 + mat
funccount); stage settings split verified (flexor 150/600/4e-3 at
minimizeFlxPin10mm_2brk.m:39-40,102-105; extensor 50/150/1e-3 at minimizeExt10mmX3.m:111-116).
Residual minor flags (dispositioned, non-failing; §8): the cam-variant parenthetical says
1,000 surrogate evaluations but the winning run used 2,000 via env override (log:
"exceeded the function evaluation limit"); the "both stages ... seeded" clause is
flexor-only (extensor InitialPopulationRange commented out). NOT RUN: bolen_isometric_2025
n-gram diff (Gate 1b) — the unpublished manuscript is not in the repo; the fix-delta
itself contains no external verbatim and no new citations.

## 5. Figure QC — verdict table (fixes are Ben's; no figure was edited, regenerated, or restyled)

| Figure | Verdict | Reason (one line) |
|---|---|---|
| Vas_Pam_20mm_cam_winner_1.png | **FAIL** | orange vs yellow "Required (+5%)" curves near-overlap in hue and position; dash pattern is the only separator; yellow-on-white reads weakly |
| Vas_Pam_20mm_cam_winner_2.png | **FAIL** | gray dashed `zero` vs gray dotted `required +5%` differ only by dash (indistinguishable in grayscale); dotted line merges with the x-axis baseline |
| Vas_Pam_20mm_cam_winner_3.png | **FAIL** | `0 (no stretch)` and `KMAX` use identical gray dotted styling; the "0" line barely visible |
| Vas_Pam_20mm_cam_winner_4.png | **PASS** | single curve, unit-bearing axes, ~8 pt effective ticks; nothing to regenerate |
| route_check_run8_20261014.png | **FAIL** | 8-panel 1817×900 px native scales to ~0.38× at a ~6.2 in print width → ~2.5–4 pt tick text, sub-legible short of a full landscape page |
| route_check_run8b_20261014.png | **FAIL** | same ~0.37× print collapse; p3–p8 route cluster an unreadable blob; p4/p6 orange-vs-yellow-only separation; θk only in tiny titles |
| Knee_Extensor_20mm.gif (first frame) | **FAIL** | frame 0 is 155×807 px, aspect destroyed, no legible text/labels/legend; unusable as a print still |

Cross-figure note: all plotted knee-angle ranges (−113..+10°) are consistent with the
dissertation convention (negative = flexion), but no figure states it — a caption-level
fix, not a regeneration defect; it drove no verdict. Method (per figure_qc.md):
dimensions/DPI + GIF frame-0 via PIL; content by AI vision with tick-glyph pixel heights
cross-checked; print scale grounded in ProofFinal's actual includegraphics widths.
**Not checked:** compiled-PDF rendering (no \includegraphics exists for these candidates
yet), grayscale conversion by a specific print provider, GIF frames beyond the first.
Full detail: [figure_qc.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/figure_qc.md>).

## 6. Build of record — FAILED

The final build agent returned `ok:false` with exactly:

```
! I can't write on file `main.pdf'.
! Emergency stop.
!  ==> Fatal error occurred, no output PDF file produced!
```

`needsRerun:false` as delivered by the build agent; recorded here as-is. This is
consistent with the ankle supervisor's finding that the canonical build environment is
blocked by a viewer lock / invocation issue: build_wf/main.log (07:32) is a fatal abort at
02-abstract, and Dissertation/build_wf_anklecheck/texput.log (07:46) an emergency-stop
abort with no output — while 02-abstract.tex on disk is clean ASCII (both aborts look
transient/environmental). The one successful full compile this campaign is the ankle
supervisor's independent scratch rebuild (§4, 40-ankle: 250 pages, exit 0, only the 6
pending steward-merge bib keys undefined). **A clean canonical build remains
outstanding** and should be re-run once the main.pdf lock is cleared. Also note the
restructure-phase build gotcha: bibtex must run from build_wf with ABSOLUTE
BIBINPUTS/BSTINPUTS pointing at ProofFinal — a relative path silently opened the wrong
thesis.bib and produced 130 missing-entry warnings.

## 7. Reward ledger

**Score: 108.** The **+1 baseline** component records the advisor's 2026-10-15 assessment
— "things are looking pretty good" (ajh26; already the reward signal of record in
[DISS_STATUS.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/DISS_STATUS.md>)).
The remaining weighting is **per the campaign ledger**
([comment_ledger.json](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/comment_ledger.json>)).
The decomposition beyond the +1 baseline was not part of this report's inputs and is not
inferred here. Context for the score: the ledger opened with 46 C-threads, 15 tracked
changes (T01–T15), and 32 resolved R-items; §2 accounts where the 46 C and 12 T ids
landed.

## 8. Needs Ben — consolidated queue

### 8.1 Comment-id integrity (highest priority — affects recorded mandates)
1. **C12-as-quoted and R20 (whole-work chapter)**: the writer's claimed mandates are
   UNRESOLVABLE — quoted text appears nowhere in comment_ledger.json, migration_map.md,
   source_sync.md, summary files, or CHATGPT_HANDOFF.md (all searched by the supervisor).
   The ledger's real C12 is the ankle balance-testbeds comment (addressed in 40-ankle).
   Action: ratify the retitle "The Work as a Whole and Future Work" or take the writer's
   offered revert (edit is safe either way — labels ch:overall/ch:futurework preserved,
   inbound refs verified); and reconcile the writer's comment-id source against the
   ledger.
2. **24 inferred anchors awaiting Ben**: C10, C11, C24–C34, C35, C36, C37, C38–C43, C45,
   C46 (full list as delivered: C10, C11, C24, C25, C26, C27, C28, C29, C30, C31, C32,
   C33, C34, C35, C36, C37, C38, C39, C40, C41, C42, C43, C45, C46). These were addressed
   under anchors the writers inferred, not anchors the ledger pins — Ben's confirmation
   retires the inference.
3. **Prior-walker placement in 60-sns** rests on the writer-attested "Ben checkpoint 3a
   override" (ledger C15 says newHome=background; no campaign file records the override).
   Document state is coherent; needs confirmation + a ledger annotation.

### 8.2 Text rulings
4. **T03/T04 accepted-TC wording**: 15-background.tex:183 now reads "these tools were
   used to convert..." (passive) while the advisor-accepted TC state was "I used these
   tools to convert" — the accepted wording exists nowhere in the document (grep = 0
   hits). Restore the TC form (here or in 60-sns) or ratify the passive and reconcile the
   Overleaf review state.
5. **Background:183 vs 60-sns:71 contradiction**: "faster than real time" vs "Full-body
   closed-loop co-simulation is not yet running". The SNS statement is the honest one;
   fix on the background side.
6. **Intro wording**: (a) pre-existing ~27–40-word reuse from published bolen_2026
   (grandfathered; not covered by the 20-muscles publication note) — reword or
   consciously extend the note's umbrella; (b) "afferent feedback did not generalize to
   weight-bearing motion" asserts a failure mode ch:sns does not record — recommended
   rewording supplied by the supervisor; (c) optional pathway-notes wording alignment
   with 50-database.tex:18.
7. **Muscles**: ruling on a dedicated adaptation pass for sec:disc_* prose (pre-existing
   ≥10-word verbatim runs from bolen_2026, disclosed by the publication note); the
   text-vs-table RMSE discrepancy is preserved deliberately (matches the published source
   of record).
8. **Torque residuals (minor)**: cam-variant budget precision ("1,000 by default; the
   winning run used 2,000"); scope the "seeded" clause to the flexor stage.
9. **SNS minors**: 40 mm Fmax wording; "air-stepping gate" → gait/configuration; three
   grandfathered "% ALT TEXT:" comments.

### 8.3 Structural/content dispositions
10. **Ch.6 ↔ Ch.8 database stats**: 50-database.tex now carries retained-corpus stats
    (561/523/86/258 at 50-database.tex:18) while 70-overall.tex:6 carries full-audited-
    corpus stats (833/95/363 of 943). Both sourced, both explicitly scoped — they do not
    contradict, but there is no cross-reference, and "95" now names three distinct
    quantities across chapters (95 afferent-classified; 95-of-99 Airtable resolutions; 95
    pruned records). Decide deliberately whether ch.8 keeps the full-corpus stats or
    defers to ch.6; do not silently harmonize (style-contract §5).
11. **sec:methods_overview removal** (writer-initiated, no ledger mandate; verified safe
    — zero dangling refs, content duplicated by chapter openings). Ratify or overrule.
12. **C10 background figures**: candidate list with licensing is in the chapter's
    deferral (Deng 2019 Figs 1/2/4 MDPI OA; Rybak/McCrea schematic; Shevtsova 2016; Di
    Russo 2023 Fig 5/2 JNE OA; Shinohara 2025; Stuart 2007 Fig 7) — Ben picks, then
    insert with captions/license lines.
13. **Ankle attested claims**: Ben's standing confirmation of the checkpoint-10 testbed
    claims (balance+walk, ~5°, indefinite hold, peterka protocol match), ideally with
    video/record locations. Until then the section honestly claims no values.
14. **Grandfathered**: 30-torque.tex:472 FLAG-BEN (measured-torque paragraph pending test
    results); the synthesis "meet or exceed" claim is supported at the redesign/model
    level only (30-torque.tex:407).

### 8.4 Dependencies
15. **Bib merge (owner: steward)**: 6 keys queued in
    [needed_bib_keys.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/needed_bib_keys.md>)
    (5 from 15-background + peterka_2002 from 40-ankle; severin_1967 carries its
    bibliographic caveat). Until merged, the document builds with exactly those 6
    undefined citations.
16. **Canonical build**: clear the main.pdf viewer lock, then re-run the clean build (§6).
17. **Figure decisions**: accept/regenerate per the §5 verdicts — all fixes are Ben's;
    nothing has been changed.

## 9. Flagged items from the sync / restructure phases (carried verbatim in substance)

**Sync (2026-10-15 zip → ProofFinal; zip-is-canonical):**
- T01–T15 are still PENDING in Overleaf's review UI, but the zip (and now ProofFinal)
  already carry their content APPLIED. **If Ben rejects any tracked change on Overleaf,
  the Overleaf file reverts and re-diverges from ProofFinal — re-run the sync after any
  reject.**
- T14 (delete "Regenerating the results files updates this table." from 94-AppendixC)
  does NOT match the live file: the sentence is still present in BOTH copies (lines 205
  and 223, grep-verified). Disposition = Ben's ruling; nothing changed either way.
- The "28-line introduction" premise did not hold: ProofFinal's intro measures 50 lines /
  14,529 B and was already byte-equal to the Overleaf intro modulo the figs-path line and
  EOF newline. If a 28-line working copy exists elsewhere, Ben should point at it.
- No upload/ mirror exists anywhere under Documentation/.../Dissertation/ (searched) —
  the do-not-touch instruction had no target; the next Overleaf push recreates the
  mirror.
- Three 15-background figure resizes and 92-AppendixA's [h]→[htbp] were untracked direct
  Overleaf edits (no T-id, author not in the ledger); ported as-is under the
  zip-is-canonical rule.

**Restructure (Phase A):**
- 50-database.tex received a skeleton only — no verbatim section moved in (Phase B
  placeholders); the substantive SADb description stays in 15-background sec:bg_sns for
  now. (The database writer subsequently filled the chapter — see §8.10.)
- sec:foundation (old 30-results lines 7-9) was placed whole into 70-overall because it
  spans four topics; Phase B may redistribute it.
- Label trade-off: ch:methods → 20-muscles (Ch.3), ch:results → 30-torque (Ch.4) keeps
  "Chapters 3-4" ranges monotone, but two refs render one chapter off semantically and
  need Phase B re-pointing: 15-background.tex:86 and 15-background.tex:134.
- 10-introduction.tex:36-50 (Organization) hard-coded the OLD structure at restructure
  time; the intro writer has since rewritten it to the new structure (verified by the
  overall/intro supervisors), but Phase-B consistency checks should re-confirm.
- T02–T07 were already applied by the sync before the restructure; each was re-verified
  by exact-string search in its new home (T06 now in 40-ankle; T07 in 15-background's
  prior-walker section).
- bibtex absolute-path gotcha — see §6.
- Figures/60-sns/CPG_airstepping_figs/ (49 files, unreferenced by ProofFinal) moved
  wholesale from 30-results; external drafts referencing the old path (e.g.
  Dissertation/CPG_spinal_section_draft.tex) need updating.

## 10. Next step — Overleaf upload (gated on Ben)

The upload bundle is **ready**. **Upload only on Ben's go**, via the `latex-overleaf`
skill. At upload time:
1. Push the ProofFinal state to Overleaf (this also recreates the absent upload/ mirror
   per §9).
2. After upload, Ben should disposition **T01–T15 in the Overleaf review UI** — any
   reject reverts the Overleaf side and requires re-running the sync (§9, first bullet).
   T14's target sentence is still live in both copies (§9).
3. Re-run a clean canonical build after the main.pdf lock is cleared (§6), then the
   Overleaf PDF refresh.

---

*Reporter's note: this report touched exactly two files — this CAMPAIGN_REPORT.md (new)
and the chapter-file rows of
[DISS_STATUS.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/DISS_STATUS.md>)
(updated in place). No dissertation source, figure, or campaign artifact was modified.*
