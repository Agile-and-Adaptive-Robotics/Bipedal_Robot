# Cut map — advisor-driven reduction pass, 2026-10-22

Baseline: fresh live Overleaf source (captured 2026-10-22), copied to
`baseline/chapters/` in this folder. Baseline narrative/caption word counts from
`measure_words.py` (narrative = prose outside captions and table bodies):

| Chapter | Narrative | Captions | Status |
|---|---:|---:|---|
| 10-introduction | 2,093 | 84 | targeted repairs only |
| 15-background | 5,217 | 529 | targeted repairs + comp-neuro de-duplication |
| 30-torque | 8,993 | 1,517 | reconcile with Oct-22 test records (bounded) |
| 50-database (SADb) | 2,588 | 338 | restructure, target ~25-35% narrative cut |
| 60-sns | 10,851 | 1,069 | restructure per Ben's 5-part order, target ~40-50% cut |
| 70-overall | 3,383 | 240 | mirror-language touch-ups only |

## 60-sns (SNS) — restructure to Ben's five-part order

New section order (all existing \label commands retained; every externally
referenced label preserved: sec:prior_walker, sec:preliminary_sim_methods,
sec:spiking_methods, sec:sns_tuning_methods, sec:preliminary_sim_results,
sec:gait_library_results, sec:walker_variant_results, sec:disc_spiking,
sec:sns_controller_comparison):

| Baseline block (live lines) | Action | Purpose | Comments addressed | Est. words saved |
|---|---|---|---|---|
| Preamble ¶1-2 (L4-6) | replace | fold into short opener + AnimatLab section | SNS-01, SNS-07 | 480 |
| Preamble ¶3 (L8) | move+compress | one lineage sentence into opener; rest cut (duplicates Background) | SNS-01 | 180 |
| Preamble ¶4 (L10) | replace | "laboratory" ¶ cut; its new-combination list moves to New Tools; engine-reuse rationale moves to gaps section | SNS-03, SNS-06, SNS-07 | 340 |
| Preamble ¶5 roadmap (L12) | replace | 3-sentence roadmap in the new order | SNS-01 | 150 |
| §prior_walker (L14-20) | keep+reframe | PART 1: the Deng+Li combination, air stepping, hand-tuning failure; first-person; name the combination | SNS-08, Ben structure | -60 (net add) |
| NEW §What the attempt exposed | new | PART 2: the four gaps from Ben's replies, each naming the tool limitation | SNS-01, SNS-06, SNS-14, Ben structure | -480 (net add, replaces preamble) |
| §existing_tools (L33-37) | compress | one ¶; SADb segment-list ¶ cut (kept at L94 wiring point) | SNS-13, SNS-14 | 300 |
| §new_tools (L39-43) | replace | inherited vs built; capability/evidence table; status classes (complete/partial/provisional/proposed); Hunt + Fletcher Young comparison from actual sources (hunt_development_2017, young_analyzing_2019) | SNS-15, Ben structure | 220 |
| §AnimatLab org (L45-51) | compress | keep lineage-archive role + reorganized-model status; keep 3 figures | SNS-11 (figure already fixed) | 260 |
| §MuJoCo conversion (L88-96) | compress | keep census, repair-scripts, BPA interface, harness honesty | — | 420 |
| §Simulink library (L112-116) | compress | keep block inventory, reflex demos, CAD defect | SNS-15 | 400 |
| §ports (L134-138) | compress | keep | — | 120 |
| §crossengine_methods (L140-142) | compress | protocol detail to one ¶ | — | 170 |
| §spiking_methods (L144-156) | compress | calibration/gates detail already in App. B; keep design rule + gate list | — | 560 |
| §tuning (L158-164) | compress | keep optimizer, gates, reference-library rationale | — | 260 |
| §Results intro (L166-168) | replace | verification framing: expected/test/observed/limit | Ben structure | 90 |
| §phase-1 (L170-172) | compress | keep honesty limits | — | 130 |
| §simulink_demos (L181-189) | compress | keep all verified numbers, tighter prose | — | 380 |
| §crossengine_results (L191-195) | compress | keep exact-match numbers + latch finding | — | 130 |
| §animatlab_validation (L204-208) | compress | keep metrics; keep un-built-replications honesty | — | 160 |
| §gait_library (L210-230) | compress | keep table + negative reading | — | 200 |
| §walker_variants (L232-267) | compress | keep sweep table + reading + drive tests | — | 420 |
| §spiking_results (L269-390) | compress | keep headline outcomes + tables/figures | — | 640 |
| §Discussion (L392-443) | compress | PART 5: short handoff; spiking discussion de-duplicated against results; curriculum, BPA integration, BPA walker, walking status kept precise but tighter | Ben structure | 1,100 |
| **Total** | | | | **~6,700 cut, ~540 added ≈ net -5,600 to -6,100 narrative words (52-56%)** |

Working target 40-50% (net ~4,300-5,400); if the executed pass lands above 50%
with all verified numbers, figures, tables, and labels preserved, the excess is
reported as justified (the preamble and discussion carried the duplication).
All 10 figures, all 5 tables, all recorded numbers retained.

## 50-database (SADb) — finished-product-first restructure

| Baseline block | Action | Purpose | Comments addressed | Est. words |
|---|---|---|---|---|
| §1 ¶1 (L6) | rework | readers-and-questions lead; plain wording | DB-01, DB-02, DB-07 | -80 |
| §1 ¶2 reviews (L8) | move | to Curation Methods as starting point | DB-08 | 0 |
| §1 ¶3 (L10) | replace | the finished product: five tables, relationships, notes/pathway records defined | DB-10, DB-11, DB-12, DB-13, DB-14, DB-21, DB-22 | -60 |
| §1 ¶4 one-program (L12) | move | Ch6 mapping to chapter end | DB-15, DB-16 | 0 |
| Bare "Methods" line (L14) | replace | real \section{Curation Methods} | DB-17, DB-18 | 0 |
| Three-way reconciliation ¶ (L16) | compress | sources + verified import in 4 sentences; drop batch chronology | DB-20 | -160 |
| Note rules ¶ (L18) | compress | keep grounding rule, vocabularies, grammar | DB-21 | -60 |
| Final wave ¶ (L20) | compress | AI-assist + gates + human rulings; author as decision authority | DB-23, DB-24, DB-25, DB-26 | -90 |
| Prune ¶ (L22) | compress | criteria-first (rules/rubric/protected), counts once, archived flag | DB-24 | -120 |
| Explorer ¶ (L24) | compress | merge into product section (why own explorer) | DB-13 | -110 |
| §2 corpus counts (L26-30) | compress | counts once each (retained vs audited scoped) | DB-20 | -90 |
| §3 triage (L46-50) | compress | two presentations stated; 2024 study | DB-09 | -60 |
| §4 use (L52-103) | keep+trim | tables stay; end with Ch6 mapping | DB-16 | -80 |
| **Total** | | | | **~-910 narrative (35%)** |

Verified this pass: two presentations (Neuroscience 2023 poster bolen_sensory_2023
= the database; Neuroscience 2024 poster bolen_using_2024 = the AI comparative
study — advisor's DB-09 memory is correct); decision authority = the author
(SADb_audit README: pilot batches reviewed by Ben before continuing; prune
review docs prune_review_20261001.md / round2; every deletion approved by the
author). No invented decisions.

## 10-introduction — targeted repairs (no length change target)

| Block | Action | Comments |
|---|---|---|
| Gap ¶s (L8-16) | reorder afferent gap before pipeline gap; expand "assembly" sentence; fix "their"; reword "fourth gap" sentence; exact count replaces "hundreds if not thousands"; four-gaps sentence in chapter order | INTRO-07, -08, -10, -11, -13, -15, -16, -17 |
| L22 Steele/Muscle Mutt ¶ | connect to the Biomimetic Humanoid Project explicitly | INTRO-14 |
| Roadmap ¶s (L38-47) | drop "Gait2392" (too specific); define spiking-neuron mirrors at first use; reword afferent-feedback clause; fix "Its also scopes" typo | INTRO-21, -23, -26, -30, -32/-33 |

## 15-background — targeted repairs + comp-neuro de-duplication

| Block | Action | Comments |
|---|---|---|
| Gait2392 name (L111) | state 23 DoF / 92 muscles | BG-05 |
| Moment-arm ¶ (L113) | reword "first-class"; scope "used throughout" | BG-06, -07, -09 |
| Knee definitions ¶ (L115) | replace "stressed"; split run-on; reword "written against"; state Gait2392's knee method | BG-10, -11, -13, -14 |
| Four-bar ¶ (L119) | split run-on; fix non sequitur | BG-15, -16 |
| Master's-model ¶s (L131-152) | re-sequence (what compared stated; no saved-file narration; contraction-limit clause explained; "followed"→calculated) | BG-19, -24, -25, -26, -27, -28, -29, -30 |
| fig:bpa_graphical_abstract (L100) | PROPOSE removal (flagged for Ben; advisor says confusing/not needed) | BG-04 |
| sec:bg_compneuro (L210-214) | move CPG-specific model ¶s into the CPG section; keep section scope general | BG-43, -44, -45 |
| Motor-center figure (L36) | figure item, not text — reported to Ben (figures already re-uploaded this pass) | BG-03 |

## 30-torque — bounded reconciliation (no length target)

Inspect Knee_Torque_Test/ProofFinal/Main.tex at fd1c13f4/7862e4e2 + test records;
resolve the six-point-sweep vs "full sweep pending" contradiction; separate raw
measurements / provisional reductions / finalized comparisons / optimizer
predictions; use current reaction-point constants and LC-angle convention
(Ben 2026-10-04: torque = F·sin(LC+0.83°)·d). Adapt only dissertation-relevant
changes with exact provenance.

## Conservation rules for all batches

Every removed or moved paragraph, number, figure block, and label is logged in
the patch files. No number changes without provenance. Labels preserved.
Bibliography untouched (comparison uses existing keys only).


## Round-2 cuts (2026-10-22, after principal review: "selection, compression, and explanation — not another expansion")

60-sns (10,851 -> 6,794 narrative; captions 1,069 -> 1,096; table bodies 316 -> 463;
combined reading matter 12,236 -> 8,353, -31.7%):
- Four-tool-families paragraph reduced to a one-sentence pointer; tab:sns_capability
  is now the primary tool description (round 1 had table + full prose restatement).
- "The standard of proof..." one-liner deleted (duplicated the verification framing).
- AnimatLab-organization paragraph: instrumentation/bilateral-rebuild detail folded to
  three sentences; figure set unchanged.
- MuJoCo conversion: repair-script narration compressed; "Python reproduces MATLAB"
  parity sentence, 0.9 s rhythm-generator detail, and the Selionov tonic-drive aside
  CUT (logged; the parity fact remains implied by the audited-pipeline claim).
- Simulink library: block inventory to one sentence; CAD hinge-mate defect kept
  verbatim-precise (relevant limitation), URDF-trial detail cut.
- Ports, cross-engine methods, spiking methods: each ~40-55% shorter; protocol
  details retained in results or Appendix; sub-stepping/residual numbers moved out
  (Appendix tabulations remain).
- Results: Simulink demos + ports compressed ~45% (10 nA/20 ms pulse, 96%-recruited,
  block-count breakdown cut); AnimatLab validation per-model SDs/duty/threshold
  detail cut; gait-sweep sub-artifact sentence (-8816/-8761) cut; deletion-stimulus
  per-perturbation scores and pelvis-tilt numbers cut (table + prose summary kept).
- Spiking results: unit-level residuals moved to "tabulated residuals"; tab:spike_beer
  DELETED (six rows folded into one prose sentence: sag 2.02 vs 2.00, off 7.85,
  modulation within 0.02 deg); hybrid-subnetwork paragraph halved.
- Discussion: SIX subsections -> THREE. sns_bpa_integration + bpa_walker +
  sns_controller_comparison merged into "The BPA-Actuated Robot and the Walking
  Status" (label sec:sns_controller_comparison preserved; sec:bpa_walker and
  sec:sns_bpa_integration labels retired, no external references — verified);
  sns_tools_summary subsection DELETED (repeated the capability table; replaced by a
  two-sentence closing handoff). Walking status now stated once.
- CONTRADICTION FIXED: curriculum paragraph said the -25 fallback "sits below every
  genuine gait" — wrong direction; it sits numerically ABOVE every genuine gait
  score (genuine gaits -144 to -267; -25 > -144), which is exactly why it is a
  sentinel excluded from score comparisons, and the text now says so. Matches
  tab:walker_ablations/tab:gait_sweep captions.

50-database (2,588 -> 2,063 narrative):
- Opening restructured DATABASE-FIRST: the SADb is the subject of sentence 1; the
  controller-requirements motivation compressed from ~200 to ~130 words.
- Methods: import-batch/recovery counts cut (374-in-fifteen-batches, 47-of-55 PDF
  recovery, 24-of-70 abstract recovery -> one clause); prune three-stage enumeration
  merged into one sentence (95/290/4+47 stage counts cut, totals 436/382/54 and the
  351-protected clause kept); "each approved by me before execution" dropped
  (decision authority already stated once, polished, in the division-of-labor
  paragraph).
- CONTRADICTION FIXED: "DOI as the unique identifier of every record" vs 550-of-561
  with DOIs -> "unique identifier wherever a record carries one".

App-first ruling applied (Ben, 2026-10-22): the packaged explorer is the database's
primary form; Airtable is the editable working base; every explorer build embeds the
full 943-record audited corpus including the 382 papers pruned from the working base
(verified against export/sadb_export.json archived flags and app check_app.py:
"papers embedded: 943"). SADb_audit/app/README.md roadmap line updated to match.

bolen_design_2020 citation REVERTED from 60-sns (Ben: it is a poster schematic, not
a paper, and its connectome is not the current one). No new citation added there.
