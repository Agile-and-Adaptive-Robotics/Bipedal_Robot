# Rulings and session record — 2026-10-17/18 (writing-quality campaign)

Ben's directives from his 2026-10-17 handoff message (machining-day session), as
executed by the writing-quality campaign (workflow dwfrun-89431ed1). Submission
2026-10-21; defense 2026-11-06.

## Rulings captured (near-verbatim where quoted)

1. **Journal standard for chapters 4, 6, and 7** — "the whole torque chapter should
   be written for a journal. So should the sns chapter and the final discussion
   chapter before conclusion, as best as possible. At least make sure it is
   dissertation worthy."
2. **BPA-actuated robot belongs in the SNS chapter's DISCUSSION** — his explicit
   placement (already reflected in SNS_DISCUSSION_REDESIGN_20261017.md §B).
3. **Wean-campaign fold-in: YES** — Ben's message asks for "the scoring results of
   the robot physical models using different neural controllers" in the SNS
   discussion and states he is "not clear if any have successfully walked without a
   harness." Resolution after checking the committed easteregg2 reports
   (reports_20261016, now local and tracked): harness-free walking remains
   UNDEMONSTRATED; folded in exactly that way — human model: feedback-sustained
   bilateral stepping at zero hand-drive (−145.07 at drive 0.381, ladder to scale 0),
   still in place under the support rig; W2L robot: antiphase progression at 60%
   weaned drive under its harness. The 60-sns "robot stands without walking"
   negative was re-scoped to the untuned-default-gains case it described.
4. **χ4 named** — the Discussion's candidate-fourth-term sentence now names it
   $\chi_{4}$ (Ben speaks of "design factors Xi0-Xi4").
5. **"Pool" word** — Appendix C's table marker and the "test/validation pool"
   phrasing replaced with "fit set" / "validation set" wording (both Appendix C and
   Chapter 4, consistently).
6. **SADb tool-name alignment** — full name + (SADb) at first use per chapter,
   "the SADb"/"the database" afterward, everywhere.
7. **Stale commented numbers** — the stale pre-optimization GoF comment block in
   30-torque deleted; stray authoring comments cleaned; flagged in the writer
   report.
8. **OpenSim isometric render** — interpreted as: a print-crisp isometric render of
   the OpenSim-derived robot model (the converted Gait2392 robot-body model,
   rendered through its MuJoCo conversion) staged in Figures/30-torque/ as a
   candidate replacement/companion for the intrinsically low-res master's-era leg
   figure (BEN_TODO item 8). NOT wired into the .tex — Ben dispositions placement.
9. **Voice ("I")** — Ben predicts the advisor will want more first-person "I" to
   humanize the document. NOT mass-applied (style contract §2 says never "I";
   4 days pre-submission a full-voice change is high-risk). The advisor's own
   ACCEPTED tracked change T03/T04 ("I used these tools to convert...") was restored
   in 15-background as the single sanctioned "I". A light first-person pass in the
   personal chapters remains available on Ben's nod after the advisor sees it.
10. **AGENTS.md cleaned 2026-10-18** (this session, Ben's mandate, no deferral):
    stale results/mat-file narratives, HX711 detail, and the neuromechanics campaign
    blow-by-blow removed; contracts, paths, rulings, and current results-of-record
    preserved; pre-clean file recoverable from git history.
11. **r6-vs-r5** ("in some ways r6 is worse than r5"): the regression candidates
    this campaign repaired are the r6-era gaps the FLOW review catalogued — no
    torque chapter opener, the heavy SNS preamble, the pulley-results duplication,
    the stale ch6↔ch7 state mismatch, and the campaign-language residue. r6's page
    growth (264→282) traced to that residue, not to new content.
12. **Torque publication note re-scoped (Ben, 2026-10-18):** "the whole chapter 4 is
    an unpublished manuscript in preparation, because it is" — the note now says
    the whole chapter is the manuscript (section list removed, orphaned "also"
    gone) and the redundant "This manuscript has not been submitted for
    publication" sentence deleted. CRediT contributions kept. Applied directly in
    30-torque.tex while the workflow ran (the torque editor's item-2 "drop also"
    may therefore report a benign conflict).
13. **Tables 4.2/4.3 training-vs-validation labels WRONG (Ben, 2026-10-18):** "No T3
    means the 47 cm BPA was not evaluated. Filtered results columns 2 and 3 mean
    that tests 2 and 4 were the validation, and therefore tests 1 and 5 were the
    training." Verified against the mats (extensor `holdCols=[2 3 4]`, three held
    per fold over the 5-test pool {40, 40-tendon, 43, 43-tendon, 47}; the 46 cm
    test is not in the flxr77 pool). The tables' CELL VALUES reproduce the mats
    exactly (re-eval 2026-10-18: flexor row 77, extensor flxr77 row 32, rigid
    escape paths for the zero-correction columns, `baselineScores` for the
    extensor) — the errors are the labels and the means: flexor "47.9 =
    validation" (it was not evaluated in the noT3 identification; validation =
    45.7 + 40.6-tendon, training = 48.5 + 41.0) and extensor "45.7 = validation /
    mean(fitting) over six" (row 32's fold: fitting = 43.6 + 43.2-tendon,
    validation = 40.0 + 40.0-tendon + 46.5; 45.7 outside the pool). Corrections
    applied to tab:FlxPinGof, tab:ExtPinGof, the §sec:ongoing summary sentences,
    and Appendix C's kf(3) note once the in-flight torque editor settled.
    Provenance: values from `Testing_Data\2022_02_Festo\minimizeFlxPin10_results.mat`
    (front row 77) and `minimizeExt10mmX3_results.mat` (front row 32, pair locked to
    row 77 via LOCKSRC) — the CLEAN-NAMED records of record per Ben (2026-10-18;
    verified to carry the same fronts, picks, LOCKSRC, and baselineScores as the
    dated duplicates `..._20260908_2brkt_2trans_noT3.mat` / `..._20260920_flxr77.mat`,
    which Ben will DELETE along with the other pinned-knee result mats).
    Re-evaluated with the evaluators of record (`reeval_gof_20260926.m`, itself
    repointed to the clean names, and the rigid escape paths).
    **Ben's follow-up rulings (same day):** (a) flexor test 3 (47 cm) is OUT of
    the main body — appendix-only, noted as "improved but not used for training
    or evaluation"; (b) the extensor identification is 3 training / 3 validation
    / 3 not-used over the nine tests: training {43, 43-tendon, 46 cm},
    validation {40, 40-tendon, 47 cm}, not used {42, 42-tendon, 48 cm}. Applied:
    Table 4.2 reduced to the four identification tests (means over the correct
    pairs: training 8.656→1.624, validation 7.588→1.490 RMSE), Table 4.3
    relabeled (means: fitting 2.498→0.997, validation 2.554→1.098 RMSE), the
    §sec:identification and §sec:ongoing sentences corrected, Appendix C Use
    columns + prose + tab:app_gof tags fixed.

## Where Ben's other open items stand after this session

- IK-references question ("where do you have those saved?"):
  `Code\MuJoCo_SNS\spinal\gait_refs\*.npz` — 43 cycles (Arnold running 8×4, Ong
  walking 0.5–2.0 m/s + self-selected, Falisse predicted). Staging and intake
  history: reports_20260923/goal5_thumbdrive_campaign.md + AGENTS.md.
- The prior session's six handoff-named rulings (χ4 sentence, SADb alignment,
  pool word, stale commented numbers, wean fold-in, isometric render) were NOT
  previously recorded in BEN_TODO — they are itemized above (items 4–8) because
  the prior session's final report was lost with its context.
- Execution status of every edit: see the workflow's revision report artifact and
  the per-writer EditResults (workflow dwfrun-89431ed1); verification build of
  record this pass: ProofFinal build_r7 (fresh directory).
