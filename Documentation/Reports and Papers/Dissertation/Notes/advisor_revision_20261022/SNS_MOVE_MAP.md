# SNS chapter (60-sns.tex) move map — baseline → proposal, 2026-10-22

Baseline = fresh live Overleaf source (sha256-prefix 9d1c581a, 444 lines, 10,851
narrative words). Proposal = `proposals/60-sns.tex` (sha256-prefix 0f638ef9, 420
lines, 7,868 narrative words; net cut 2,983 words, 27.5%). All 42 baseline labels
kept (2 added: `sec:sns_gaps`, `tab:sns_capability`); all 19 figure files kept;
all 6 tables kept verbatim; every recorded number preserved. Paragraph
conservation accounting, baseline line numbers:

| Baseline ¶ (line) | Disposition in proposal |
|---|---|
| Preamble ¶1 (4) mechanics/validation + SNS definition | KEPT, compressed ~15%, now opener ¶1 |
| Preamble ¶2 (6) "Two definitions" vocabulary | REMOVED — Background Sections 2.5–2.8 already define it; one pointer sentence retained in opener ¶1 ("defined in Sections 2.5–2.8 and used here without redefinition") |
| Preamble ¶3 (8) lineage/SNS-line/ancestry | COMPRESSED to one clause in opener ¶1 (citations to rybak/mccrea/deng kept there); nishizaki/fukuoka sentence CUT (peripheral; claim lives in Background) |
| Preamble ¶4 (10) "This chapter concerns tooling" laboratory ¶ | REPLACED — its new-combination list became §6.3 New Tools ¶1; its engine-reuse rationale became the final ¶ of §6.2 (the gaps section); its "re-running adds nothing" argument CUT (one trace survives in the Hunt/Fletcher ¶) |
| Preamble ¶5 (12) roadmap | REPLACED by 3-sentence roadmap ¶ in the new chapter order (opener ¶2) |
| §prior_walker ¶1 (16) | KEPT with the combination named explicitly (Li biomechanics + Deng controller, first person "Morrow and I") |
| §prior_walker ¶2 (18) biomechanics | KEPT ~95% (trimmed "standing in for" clause) |
| §prior_walker ¶3 (20) controller + hand tuning + NIH | KEPT ~90%; ending sharpened: air-stepping stated plainly as "never produced walking … never exercised under load" (answers SNS-08/"Preliminary Results:" removal by the advisor's applied edit) |
| NEW §6.2 What the AnimatLab Attempt Exposed (4 gap ¶s + engines ¶) | NEW (~470 words) — replaces the removed preamble; the four gaps are Ben's comment-reply list (SNS-01): Gait2392+CPG and robotbody+CPG couplings (gap 1), Li walker + Deng controller with contact feedback locked in AnimatLab (gap 2), graphical block programming on maintained tools incl. the Hunt/Young predecessor libraries (gap 3), CAD + BPA + controller in one loop (gap 4); the several-engines rationale is Ben's SNS-06 reply |
| §existing_tools ¶1 (35) | KEPT, compressed ~35% into one ¶ (SCONE detail kept: headless practice baseline + balance tutorial anchor) |
| §existing_tools ¶2 (37) SADb segments | REMOVED per advisor SNS-13 ("This whole paragraph seems unnecessary") — the SADb linkage survives once at the L94 wiring point and once in the new ¶1 |
| §new_tools ¶1 (41) four families | COMPRESSED + Table tab:sns_capability added (capability/evidence/status: complete, partial import, provisional, proposed — Ben's status classes) |
| §new_tools ¶2 (43) re-implementation rationale | Cut here; rationale lives in §crossengine_methods; last sentence (standard of proof) kept as closing ¶ of §new_tools |
| NEW Hunt/Fletcher ¶ | NEW (~180 words) answering SNS-15 with actual sources: hunt_development_2017 (Hunt's hindlimb-walking controller in the lineage's block-based Simulink libraries) and young_analyzing_2019 (Fletcher Young's full-muscle rat hindlimb model); states what the derivative adds (Ia/Ib encoders, diameter-specific BPA blocks, CAD coupling, spiking additions, regeneration) |
| §AnimatLab org (47, 49) | KEPT, two ¶s merged and compressed ~25%; figures (model architecture, motifs, GUI) unchanged |
| §MuJoCo (90, 92, 94, 96) | Four ¶s → two, compressed ~20%; every number kept (92 actuators, dynamic sign check, 2 ms, 410-neuron census, 0.9 s antiphase, SADb wiring, harness statement) |
| §Simulink library (114, 116) | Compressed ~15%; block inventory, reflex demos, CAD defect (welded hinge + re-mate fix), provisional-parameter statement kept |
| §ports (136, 138) | Two ¶s → one, compressed ~20%; Deng files, burst-termination pointer, export script, two import routes, bridge kept |
| §crossengine_methods (142) | Compressed ~15%; harness definition, W2L core, protocol, headless reruns, scope limit kept |
| §spiking_methods (146, 148, 150, 152, 154, 156) | Six ¶s → one long ¶ + gate list, compressed ~35%; design rule, selector + bit-exact proof, adapting-if rationale, calibration summary (3e-4 mV), 3 new blocks, AnimatLab conversion, 5 gates all kept; unit-test numbers left to Appendix B where already tabulated |
| §tuning (160, 162, 164) | Compressed ~15%; TPE, no-hand-tuning rule, gates, bifurcation rationale, 43-ref library, structure-as-data kept |
| §Results intro (168) | Replaced by 2-sentence verification-framing intro (expected/test/observed/limit) |
| §phase-1 (172) | Reframed in expected/observed/limit form, compressed ~20% |
| §simulink_demos (183, 185, 187, 189) | Compressed ~15%; all numbers kept (15.79°, 43.5°, 2.00/7.85°, 9.7–48.0°, 1.938/1.94 s, r=−0.86, 4.2e-06 mV, >40 kN vs 469 N, 96%, 1,446 blocks/149 joints) |
| §crossengine_results (193, 195) | Reframed expected/observed; latch finding kept; compressed ~10% |
| §animatlab_validation (206, 208) | Reframed; all metrics kept (0.450/0.467 s, 0.529 cycle, r=0.976, 2.25 Hz, −56.4/−55 mV, r=−0.996, 0.637 m/s, 1.299±0.024/1.333±0.080 s, duty 0.443/0.555); un-built-replications honesty kept |
| §gait_library (212, 230) | Compressed ~25%; table kept; −8816/−8761 means kept; negative reading kept |
| §walker_variants (234, 236, 265, 267) | Compressed ~20%; full ablation table kept; all scores kept; reading ¶ merged with drive-perturbation results |
| §spiking_results 271–390 | Compressed ~25%; every gate result, table, and figure kept; the figure ref `fig:spike_animatlab` (dangling in the baseline — no such figure block exists) was REMOVED from the text (bug fix, logged) |
| §discussion opener (394) | Replaced by 1-sentence opener |
| §disc_spiking (398–408) | Four ¶s → two, compressed ~35%; all four findings kept (interface unchanged, timing +robustness, naive-copy failure, bilateral-gait lead) |
| §sns_curriculum (412, 421) | Two ¶s → one, compressed ~25%; five stages, standing-before-walking, gated scoring kept; figure kept |
| §sns_bpa_integration (425) | Compressed ~20%; one-source statement, 60% ceiling, quasi-static limit kept |
| §bpa_walker (429, 431) | Compressed ~15%; air-stepping numbers (70% of 411 N, 1.542 s, 63.5°, 880/1372 N, 372 kPa, 7+5 contacts) and all negatives kept verbatim in substance |
| §sns_controller_comparison (435, 437) | Compressed ~15%; family-level reading, wean ladder (−161.57→−145.07, duty 0.71/0.65, 8/8 basin, drive to 0.0, dies below 85%), and the exact "Unsupported, harness-free walking has not been demonstrated" statement KEPT |
| §sns_tools_summary (441, 443) | Compressed ~15%; handoff to Chapter 7 kept |

## SADb chapter (50-database.tex) accounting, baseline → proposal

Baseline 2,588 → proposal 2,201 narrative words (cut 387, 15%). Labels: all 10
kept + `sec:sadb_methods` added. Figures/tables: all kept verbatim.

| Baseline ¶ (line) | Disposition |
|---|---|
| Resource ¶1 (6) | Reworked: plain "requires numbers" framing (DB-01/02), ipsi/contra/both (DB-03), softened "established experimentally" with the unknowns (DB-04), readers lead (DB-07) |
| Reviews ¶ (8) | MOVED into Curation Methods ¶1 (DB-08: "goes in your methods section"), compressed ~30% |
| Product ¶ (10) | REPLACED by two finished-product ¶s: five tables + relationships + notes/pathway-records defined (DB-10/11/12/21/22); "first-class" gone (DB-11); "canonical" → unique (DB-14); Airtable working home + explorer packaging stated together (DB-13); both presentations verified from bib and stated (DB-09); "I built it" credit (DB-25) |
| One-program ¶ (12) | MOVED to chapter end (DB-16) |
| Bare "Methods" line (14) | Replaced by real \section{Curation Methods}\label{sec:sadb_methods} (DB-17/18) |
| Three-way reconciliation ¶ (16) | Compressed ~40%: sources + verified import + growth + audited 943 kept; batch chronology numbers (483 DOIs, 95/99, 443, 532) dropped (DB-20) |
| Note rules ¶ (18) | Merged into methods ¶3 (grounding rule, vocabularies, flagged rulings) |
| Final-wave ¶ (20) | Compressed ~30%; AI gates + author-as-decision-authority (DB-23/25/26) |
| Prune ¶ (22) | Compressed ~35%; criteria first (inclusion scope, rules, rubric, protected classes — DB-24), counts once, archived flag |
| Explorer ¶ (24) | Merged into product ¶3, compressed ~50% |
| Corpus ¶1 (28) | KEPT minus the operating-cap clause |
| Corpus ¶2 (30) | DEDUPLICATED (the "382 records … archived flag" sentence appeared twice in one ¶ in the baseline); counts unchanged |
| Triage ¶s (48, 50) | Compressed to one ¶; both presentations stated once at chapter top and once here |
| Use ¶1 (54) | KEPT, ~10% trim |
| Extension ¶s (99, 101) | Merged and trimmed ~30% |
| Honest-accounting ¶ (103) | KEPT minus one clause; one-program mapping appended (from line 12) |


## Round-2 (2026-10-22) paragraph accounting — SNS

Deleted outright (content absorbed by tab:sns_capability or the merged discussion):
- "The standard of proof for each tool..." (1 sentence).
- tab:spike_beer (table; numbers preserved in one prose sentence).
- sec:sns_tools_summary subsection (2 paragraphs) — its unique surviving content is
  the closing handoff sentence pair ("What the tools are for is the program they
  define next..."; "Every tool is reusable...").
- The controller-comparison opening paragraph and BPA-integration subsection were
  MOVED (not deleted) into the merged "BPA-Actuated Robot and the Walking Status"
  subsection (sec:sns_controller_comparison), each compressed.
Labels retired: sec:sns_bpa_integration, sec:bpa_walker, sec:sns_tools_summary
(no external references; verified by grep across all chapters). All externally
referenced labels preserved: sec:prior_walker, sec:preliminary_sim_methods,
sec:spiking_methods, sec:sns_tuning_methods, sec:preliminary_sim_results,
sec:gait_library_results, sec:walker_variant_results, sec:disc_spiking,
sec:sns_controller_comparison, sec:sns_gaps, sec:sns_curriculum.
Every verified number removed from prose in round 2 is listed in CUT_MAP round-2
section; none were changed, and each remains traceable to its source table, figure,
or Appendix tabulation.
