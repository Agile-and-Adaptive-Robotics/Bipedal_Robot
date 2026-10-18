# Bridge Plan — Frontiers manuscript "Design Factors for Isometric Robotic Knee Torque Using Braided Pneumatic Actuators"

Written by the bridge planner from: `tmp/paper_bridge/BRIEF.md` (contract), the four mining digests
(`mining/chapter_digest.md`, `mining/results_ledger.md`, `mining/paper_state.md`,
`mining/actuators_style.md`), direct reads of `Dissertation/ProofFinal/chapters/30-torque.tex`
(lines 332–547 this session) and `tmp/paper_bridge/PlotPrompt_utf8.txt` (Ben's plotting
instructions, read in full), plus direct disk checks of `ProofFinal/Figures/` (14 files) and
`Knee_Torque_Test/Tables/` (holds only Ben's `CoefficientTable.xlsx`, 2026-05-27).

**How to use this plan.** Section writers write fragments per §2; the tables builder builds §5;
the figures regenerator executes §6; the assembler merges per §2 and obeys §8 (numbers of record)
and §7 (citations). Nothing outside `ProofFinal/` + `Tables/` (xlsx only) + `tmp/paper_bridge/`
is ever written. Ben-only decisions live in §9 (OPEN QUESTIONS) — never guess them.

**Display-item budget (Frontiers cap 15):** 9 figures + 5 tables = 14. One spare slot.

---

## 1. Target skeleton and word budgets

Calibrated against the Actuators style card §10 (body ~5,500–7,000 words, Discussion-dominant):

| # | Paper section | Words (target) | Display items |
|---|---|---|---|
| — | Title / Abstract (~150 w, 5 sentences) / Keywords (8) | — | — |
| 1 | Introduction | ~800 | none |
| 2 | Methods and Materials (2.1–2.9) | ~1,300 | testJigs, KneeICR, xiFrameTransform, bktFrame |
| 3 | Results (3.1–3.4) | ~1,100 | FlxPin_group, ExtPin_group, Ext10mm_52cm, FlxGrp |
| 4 | Discussion (4.1–4.5, LONG) | ~2,900 (≥ 2× Results) | Route figure + tables T1–T5 |
| 5 | Conclusion | ~280 | none |

Deviation from the style card's "Results ≥ 1.5× Methods" ratio is deliberate: this paper's Methods
carries the model + identification + redesign procedures (the chapter's own arrangement), so the
preserved invariant is Discussion-dominance (Discussion ≥ 2× Results), which holds.

**Architecture decision (BRIEF-mandated, mirrors the chapter):** Results reports the *initial
single-bracket identification* story; the LONG Discussion carries the *two-stage re-identification
(adopted χ0–χ3, two-bracket)* and the *route redesign* as ADDITIONAL RESULTS, with the
fictive-strain (ε620) framing used wherever correction terms are interpreted.

---

## 2. Section-by-section mapping (chapter → paper)

Action codes: **BRIDGE** = adapt dissertation text (Ben's own unpublished material; proper
adaptation, not verbatim dump), **REWRITE** = rework existing Main.tex prose to journal register,
**NEW** = written fresh for the paper, **CUT** = not carried.

| Paper element | Source | Action | Content and instructions |
|---|---|---|---|
| Title block (Main.tex:64) | BRIEF working title | REWRITE | Long: "Design Factors for Isometric Robotic Knee Torque Using Braided Pneumatic Actuators"; short: "Design Factors for Isometric Robotic Knee Torque". |
| Abstract (Main.tex:72–82) | Actuators style §3 template + this plan §4 | REWRITE | 5 sentences, ~150 words, no citations, "(1)…(2)…(3)…" contribution sentence. Draft skeleton in §4 below. Kill the enumeration-style "1)…2)" and future-tense promise of the current draft. |
| Keywords (Main.tex:80) | — | REWRITE | Replace all 8 (incl. the odd "Function Fit"); list in §4. |
| **1 Introduction** | Main.tex:84–95 (611 w) + chapter context | REWRITE + EXPAND | Actuators funnel, ~800 w, citation-dense (~20 calls): (P1) biomimetic/musculoskeletal robots need muscle-like actuators — deploy uncited bib keys (kurumaya_musculoskeletal_2016, hitzmann_anthropomorphic_2018, asano_musculoskeletal_2019, liu_robotic_2018, ijspeert_amphibious_2020); (P2) one-clause tour of alternatives → BPAs (aschenbeck_design_2006, veale_towards_2016, madden_artificial_2004, suzumori_trends_2018); (P3) force is characterized (bolen_2026, bolen_determination_2019) but joint torque depends on moment arm, routing, and fixture compliance (menegaldo_moment_2004, young_analyzing_2019, suderman_moving_2012, brand_model_1982); (P4) model-by-model narrative \citet review of BPA/torque modeling (hunt_modeling_2017, hunt_development_2017, morrow_optimization_2020, martens_modeling_2017, tondu_static_2022, sarosi_new_2012, chou_measurement_1996) landing on the gap: rigid-body torque models neglect fixture compliance, length offsets, and wrapping losses, and none is validated across joint configurations with one consistent parameter set; (P5) "The present work…" — contributions: (1) hybrid torque method isolating error sources + correction terms incl. fictive-strain wrapping loss expressed through ε620; (2) two-stage identification with a cross-configuration-consistent stiffness pair, validated on held-out tests and the biomimetic knee; (3) attachment-route redesign meeting human isometric torque targets, with first hardware points. Mention the biomimetic knee with steele_biomimetic_2018 early (P1 or P5). |
| **2 Methods and Materials** | | RESTRUCTURE | The model equations currently living in Main's Results (eqs 7–25) migrate here (the chapter's own arrangement: Methods 3.1.4). |
| 2.1 Overview | Main.tex:100 | REWRITE | Add the design-factor framing: torque = f(force, moment arm, effective length); each identified factor is a correctable design quantity. |
| 2.2 Actuator force model | Main.tex:104–158 + Appendix A eqs | COMPRESS + BRIDGE | ~150–200 w: recap F*(ε*,P*) closed form (eq:app_fstar pattern), Fmax, composition; **define ε, ε620, ε\* here** (eq:app_strain / 20-muscles.tex definitions; ε620 = (l0−l620)/l0, ε\* = ε/ε620) — this is the fictive-strain basis used everywhere later. All characterization values cited to bolen_2026 (do not re-tabulate). |
| 2.3 Robot leg model and knee test stands | Main.tex:162–195 + chapter sec:robot_model (:9–50) | BRIDGE | Leg geometry from bone scans of a 1.83 m person; pinned revolute stand and biomimetic knee (steele_biomimetic_2018, steele_development_2017); ICR migration with KneeICR figure; testJigs figure. |
| 2.4 Equipment and data collection | Main.tex:197–206 + chapter :52–61 | BRIDGE | Equipment list verbatim-in-substance (80/20 1515, crane scale, CALT DYLY-103 + HX711, MPX5700, Parker+Husky regulators, VTUG-10 manifold, BlackBoard C, 3/16 in Kevlar rope) — all through siunitx; drop the Arduino/Windows note (stays commented in Main). |
| 2.5 Hybrid torque calculation | Main.tex:207–221 (thin) + chapter sec:hybrid_method (:63–79) | EXPAND | eq:momentarm_def, eq:project_ma (Young projected moment arm, cite young_analyzing_2019), eq:M_vec, eq:M_scalar; state the hybrid's role: measured lm and rk isolate kinematic error from force error. |
| 2.6 Improved torque model: length offset and bracket compliance | Main Results eqs 7–22 + chapter sec:offset_stiffness (:86–274) | BRIDGE (compress) | Keep ~9 equations: eq:Lm (χ0), eq:bracket_frame_transform (two-rotation frame; NEW fig xiFrameTransform), eq:bktstiffness (y-symmetric diag(χ1,χ2,χ1)), the projected compliance eq:bktprojectedcompliance, eq:keff, eq:equivalentstiffness (bracket+tendon series), eq:activelengthcorrected, eq:complianceforcebalance (fzero balance), eq:bracket_series_projection (two-bracket generalization). Compress the intermediate component forms into prose. One sentence on frame-convention invariance under y-symmetric K. bktFrame figure (caption panel-D fix, see §6). |
| 2.7 Fictive strain and wrapping loss | chapter sec:fictive_contraction (:276–297) | BRIDGE | Define fictive strain (excess contraction the force model must see beyond the actuator's physical contraction, in ε620-normalized units); ζ = 1 − ε\*; eq:deltaL (Δl = χ3·R·\|θ_wrap−θ_k\|·ζ²), eq:x3 (corrected length). State the physical picture: the wrapped arc stops contributing usable contraction as the BPA bends over the joint. |
| 2.8 Identification procedure | chapter sec:identification (:299–311) + 94-AppendixC tab:app_opt_bounds + ledger §7 | BRIDGE | Multiobjective genetic algorithm (NSGA-II-type, gamultiobj) settings: flexor population 150, ≤600 generations, spread tolerance 4e-3; extensor refit population 50, ≤150 generations, tolerance 1e-3; normalized search space; bounds (flexor χ0 ∈ [0,20] mm, χ1 ∈ [3e4,1e6] N/m, χ2 ∈ [5e3,2e4] N/m; extensor χ0 ∈ [−20,+5] mm, χ3 ∈ [0,1]); two-stage design (flexor two-bracket front first; extensor refit with stiffness pair locked to the flexor front); train/validation splits (flexor: 4 training tests + 47.9 cm held out; extensor pool {40.0, 40.0-tendon, 43.6, 43.2-tendon, 46.5} fitting + 45.7 cm held out; three tests excluded from the pool); GoF metrics RMSE/FVU/Max defined. |
| 2.9 Route-redesign procedure | chapter sec:optimization (:313–330) + sec:route_redesign_results setup (:506–508) + Appendix A sec:app_opensim_targets | NEW (method; results stay in Discussion) | Gait2392-derived robot-body model; human isometric torque targets = OpenSim tendon-excursion moment arms × max isometric force (vastus medialis extensor, biceps femoris short head flexor) — revive the commented human-benchmark citations at Main.tex:225 (hoy_musculoskeletal_1990, yamaguchi_planar_1989, delp_interactive_1990, seth_opensim_2011/2018, sherman_what_2013, millard_flexing_2013, thelen_adjustment_2003 — all keys exist). eq:opt_cost (J = w_w·J_worst + w_t·J_torque + …, w_w = 1e5, w_t = 1e3); surrogate + pattern-search pipeline (budgets: extensor 3,000 surrogate evaluations, flexor 1,000, pattern-search refinement 15,000; do NOT state a cam budget — FLAG-5); hard feasibility via prohibitive penalties; required 5% torque margin above the human envelope; actuator plan: two ⌀20 mm BPAs per route, 620 kPa ceiling, maximum usable contraction 0.255; route-to-route clearance requirement 25 mm; cam variant and mono-/bi-pulley reeving variants defined in one short paragraph each. |
| **3 Results** | chapter 3.2 (:332–476) + Main.tex:228–455 | REWRITE (register; equations moved to Methods) | Report, don't interpret. All numbers at the INITIAL identification (values of record per §8). |
| 3.0 lead | chapter :334 | BRIDGE | Three predictions compared per configuration: original rigid-body, improved BPA characterization alone, corrected model; hybrid used to isolate error sources. |
| 3.1 Flexor on pin joint | chapter :336–362 + Main.tex:231–378 | REWRITE | Overprediction of original + char-only; hybrid matches → kinematics mispredict lm/rk. Initial coefficients in running text (Actuators Pattern A): χ0 = 11.9 mm, χ1 = 102 kN/m, χ2 = 10.6 kN/m. Initial GoF: 48.5 cm fit RMSE 1.616 / FVU 0.0958 / Max 3.2527; 45.7 cm validation 2.1214 / 0.0878 / 5.1690 (chapter :348; ledger §2, DR). Figure: regenerated FlxPin_group (2×2, spec in §6). |
| 3.2 Extensor on pin joint | chapter :385–398 + Main.tex:380–418 | REWRITE | Zero-point shift → χ0 opposite sign; load-dependent divergence at higher torque → χ3. Initial extensor terms in text: χ0 = −10.9 mm, χ3 = 0.123 (stiffness pair shared). **Figure-led: no numeric initial-GoF exists in the record for these tests — do not invent any; the adopted-terms numbers arrive in Discussion 4.2 (T2).** Figure: regenerated ExtPin_group. |
| 3.3 Extensor on biomimetic knee | chapter :423–444 + Main.tex:420–442 | REWRITE | l0 = 51.8 cm everywhere (fix Main's 52.0); measured below prediction; initial-terms GoF RMSE 0.8920 / FVU 0.3396 / Max 1.4188 (chapter :430; ledger §2). Figure: regenerated Ext10mm_52cm. Pass-through difficulty references the SUPPLEMENTARY figure (see §6, passthru decision). |
| 3.4 Flexor on biomimetic knee | chapter :448–476 + Main.tex:443–455 | REWRITE | 10 mm, l0 41.5 cm + 21 mm tendon @ 604 kPa: RMSE 0.3007 / FVU 0.0102 / Max 0.4647. 20 mm @ 620 kPa (label per §8): RMSE 1.8024 / FVU 0.0499 / Max 3.5678. Human comparison: single-BPA torque does not exceed human magnitude ca. 60–100° flexion; outside ±20% of the human mean 78–103°; femur ultimate strength caps pressure (sub-620 points). Figure: regenerated FlxGrp. Close by setting up the redesign motivation (one sentence, no interpretation). |
| **4 Discussion (LONG)** | chapter 3.3 + appendices | see §3 below | ~2,900 w, subsections 4.1–4.5; carries ADDITIONAL RESULTS. |
| **5 Conclusion** | — | NEW | One paragraph ~280 w, no citations: recap hybrid method → correction terms (initial and adopted, one consistent pair) → fictive-strain expression → redesign meeting human targets → first hardware points → the design-workflow payoff (first-iteration torque design for BPA-driven joints). |
| Back matter (RII, CRediT, funding, acknowledgments, data availability, supplement) | Main.tex:476–507 | KEEP | Do not touch the author block. CRediT defect ("Investigation, B.B., C.M., and C.M.") is flagged (§9), not fixed. |
| Supplementary | frontiers_SupplementaryMaterial.tex | UNTOUCHED | Workflow rule. Main text stops referencing supp tab:TorqueFit (its numbers disagree with the record — §8); passthru photo stays there and is cited from Discussion 4.1. |

Cleanup that rides along (assembler; from paper_state §7/§10): add `\usepackage{mathabx}` (flagged
missing) and convert all 8 `$\phi$` diameter uses to `$\diameter$` (supp untouched); all bare units
→ siunitx; the typo ledger (L88 "robot robot", L373 "the the", L384 "offet", L415 "gamultobj",
L445 "acheivable" ×2, garbled L462 splice, panel-D caption "(A)"→"(D)", missing "Fig." before
\ref at L445) — all fixed or deleted by rewriting the containing text. Remove the dangling scratch
fragment at Main.tex:469–474 (its content is properly written in Discussion 4.2f).

---

## 3. Long-Discussion outline (ADDITIONAL RESULTS, fictive-strain framing throughout)

Discussion word budget ~2,900; ~26 citation calls; interpretation-level math allowed (Actuators
precedent). Every number below is DR or DR-verified in `mining/results_ledger.md` (§ refs given).

**4.1 Sources of error in isometric torque generation** (~500 w; BRIDGE chapter sec:disc_error
:480–498 — overlap section, adapt Ben's own text)
- Compliance in winch cable, artificial tendons, and brackets → the χ0/χ1/χ2 family (chapter :484).
- Wrapped-path modeling error; short moment arms amplify path error (young_analyzing_2019);
  the deliberate no-patella simplification (chapter :486).
- BPA stretch/compression at high flexion and ±20 mm Z-slip off the femur bolt (chapter :486;
  testJigs C).
- The ExtPin load-dependent divergence → loss of functional lm → why χ3 was created; the
  fictive-strain reading introduced here in one sentence (the wrapped arc stops contributing
  usable contraction; formalized in 4.2). Reference the supplementary pass-through figure.
- Close on the simple→complex progression (chapter :491).

**4.2 Follow-up identification: two-stage re-identification and the adopted correction terms**
(~700 w; NEW prose from chapter sec:ongoing :500–504 + Appendix C; ADDITIONAL RESULTS)
- (a) The coefficients reported with the results figures are the initial single-bracket
  identification; re-identifying over the full test sets with the two-stage procedure (Methods
  2.8) yields the adopted values.
- (b) **Table T3 (adopted terms + bounds):** flexor χ0 = +3.94 mm, χ1 = 3.998e4 N/m,
  χ2 = 1.473e4 N/m; extensor χ0 = −6.38 mm, χ3 = 0.158 on the same locked pair (ledger §1,
  DR-verified). Fictive-strain reading of each term: χ0 a constant effective length offset;
  χ1/χ2 the lumped effective stiffness of bracket, fixtures, and winch cable (chapter :484
  supports the lumping); χ3 the scale of the wrapping-length loss Δl = χ3·R·\|θ_wrap−θ_k\|·ζ² with
  ζ = 1 − ε\*: the loss is expressed in ε620-normalized (fictive) contraction units, vanishes as
  the actuator approaches full contraction (ζ → 0 at ε\* → 1), and grows quadratically toward
  resting length.
- (c) Non-identifiability: the pinned-flexor data admit a family of near-equivalent stiffness
  pairs; the adopted pair was selected by a requirement external to the fit — both redesign
  configurations must meet the required 5% torque-overshoot target simultaneously (flexor
  +5.00%, extensor +5.043% at −20.2°; ledger §1).
- (d) What the corrections buy — **Tables T1 (pinned flexor) and T2 (pinned extensor)**: training
  mean RMSE 8.12 → 1.56 N·m (4-test mean, state the set); held-out 47.9 cm 2.29 N·m; extensor
  fitting mean 2.60 → 1.16 N·m; held-out 45.7 cm 0.49 N·m; FVU < 1 everywhere except the 43.2 cm
  tendon test; biomimetic extensor RMSE 3.03 → 1.17 N·m, FVU 4.49 → 0.67 (ledger §3).
- (e) Solver assignment rationale: genetic multiobjective retained where the Pareto front is
  needed for holdout validation; surrogate + pattern search for the redesign (chapter :502).
- (f) Two-to-three-sentence interpretation (replaces the dangling Main.tex:469–474 scratch): the
  corrected force balance reads as an effective spring–damper in ε\* with k(ε\*,P\*) =
  ∂F/∂ε\*\|_{P\*}; setting mass and external force to zero gives the effective damping
  c = k·ε\*/ε̇\* (present tense, no working-math residue; chapter :493–498).

**4.3 Route redesign meeting human isometric torque targets** (~800 w; NEW from chapter
sec:route_redesign_results :506–538 + sec:disc_pulley :540–544; ADDITIONAL RESULTS)
- (a) Motivation (from 3.4's shortfall) + shared actuator plan (two ⌀20 mm BPAs, 620 kPa, max
  usable contraction 0.255); targets = OpenSim isometric profiles (vastus medialis, biceps femoris
  short head) + 5% overshoot envelope.
- (b) Direct extensor route: predicted two-BPA torque 35.45–109.4 N·m vs requirement
  20.25–50.32 N·m; zero shortfall angles; thinnest margin +5.02% at −20.2° (ledger §6.1).
- (c) Cam redesign (current extensor design of record): min margin +30.83% vs required envelope,
  +39.77% vs plain human target (winner pre-adjustment +31.30%/+42.34%); envelope pressure span
  269–507 kPa vs the 150 kPa regulated supply floor; wrap bend radius 96.9 mm; condyle clearance
  43.1 mm at 10° extension; moment arm 59.4–102.8 mm, mean 72.4 mm (ledger §6.2). NEW figure 9
  (renamed copy of the dissertation route-geometry figure — §6).
- (d) Direct flexor: −81.23 to −3.864 N·m vs requirement −23.16 to −2.461 N·m (peak 23.16 N·m);
  min margin fraction 0.0501 (+5.01% at −107°); stated assumption: both actuators share one
  attachment route; mirrored pair = next refinement (ledger §6.3).
- (e) Mono-pulley variant: peak 85.74 N·m (85.56 insertion-side) vs 81.23 direct → +4.51 N·m
  (5.6%); off-axis ≤ 0.803 N·m (<1% of peak); 1:1 redirect not a multiplier; objectives not
  commensurable across variants (0.08465 vs 74.47) — keep the honest caveat (ledger §6.4).
- (f) Bi-pulley batch — **Table T5**: all five muscle sites nonpositive margin; soleus 0.9968 /
  tibialis anterior 0.9288 (converged, flag 1), rectus femoris 3.423 (worst constraint 0.2423),
  medial gastrocnemius + biceps femoris short head at the infeasibility sentinel (flag −2);
  finite clearances 4.26–19.14 mm vs 25 mm required; BFSH target 40.2 N·m at 804 N; screening
  conclusion: reeved routes do not fit these corridors; redirect effort to clearance-relaxed
  sites (ledger §6.5; chapter tab:bipulley_batch).
- (g) First physical measurements of the redesigned ⌀20 mm extensor route near full extension:
  1.163 / 7.849 / 23.437 N·m at ~10–11°, hung loads 3.7/25.2/76 N, 388.5/388/592 kPa; the first
  and third points span a twentyfold torque range at essentially one angle; full sweep and flexor
  measurements pending (chapter :544; ledger §5 DR-verified — subject to OPEN QUESTION 3).
- (h) Honest-status paragraph, whole-body extension (chapter :538; ledger §6.6, chapter-only):
  first pass 19/27 with all-but-one routes infeasible (BFLH +0.03%); restructured pass 12/27
  terminated on a route-builder error, 15 finalized, 8 constraint-feasible, 6 meeting torque at
  every angle (worst-angle margins +0.25% to +250.7%). Status only, no table, no artifacts quoted
  beyond the chapter. **CUT: adductor-magnus restoration** (robot-body model bookkeeping,
  off-scope; §9).
- (i) Scope caveat: redesigns are optimizer-stage results, not yet propagated into the simulation
  model trees (chapter :536, :542).

**4.4 Design factors for isometric torque** (~450 w; NEW synthesis; no new numbers)
- The design-factor triangle: moment arm r, normalized contraction ε\* (usable force), and path +
  fixture compliance (the χ terms). Practical instructions in the Actuators direct-address
  register: budget stiffness when moment arms are short; route to keep ε\* high through the
  functional range because the wrapping loss scales with (1 − ε\*)²; moment-arm vs
  contraction-budget trade-off (the cam buys 59.4–102.8 mm of arm with margin and clearance);
  pressure headroom as a design resource (269–507 kPa vs the 620 kPa ceiling).
- Comparison-with-literature paragraph deploying uncited bib keys: musculoskeletal legs and
  tendon-routing work (kurumaya_musculoskeletal_2016, hitzmann_anthropomorphic_2018,
  asano_musculoskeletal_2019, liu_robotic_2018, he_mechanism_2020, romero_four-bar_2020,
  chen_realizing_2019, leibach_development_2020, schilling_adaptive_2022); credit-and-compare
  tone (§7 of the style card).

**4.5 Limitations and future work** (~350 w)
- Shared-route assumption (mirrored pair pending); optimizer-stage routes not yet in the
  simulation trees; physical validation begun but incomplete; the stiffness pair selected by an
  external design requirement rather than uniquely identified — an independent stiffness
  experiment would close the non-identifiability; the fictive-strain term identified on one
  diameter family — transfer to other diameters untested; per-actuator ε620 recording practice
  (continuity with the force paper's conclusion); route-builder concatenation fix and
  larger-actuator/relaxed-corridor options for the short whole-body routes.

---

## 4. Title, abstract, and keywords strategy

**Title (long):** Design Factors for Isometric Robotic Knee Torque Using Braided Pneumatic
Actuators. **Short:** Design Factors for Isometric Robotic Knee Torque.

**Abstract** (~150 words, 5 sentences, zero citations — Actuators template sentence by sentence):
1. *Motivation + asset:* BPAs reproduce muscle-like force and contraction behavior, making them
   attractive for musculoskeletal robot legs.
2. *Gap:* predicting the isometric torque a BPA delivers about a joint remains poor because
   attachment compliance, length offsets, and routing losses over the joint are absent from
   rigid-body torque models.
3. *What we did:* we measured isometric knee torque on pinned and biomimetic knee test stands
   across actuator resting lengths, isolated error sources with a hybrid torque calculation, and
   identified correction terms by multiobjective optimization.
4. *Enumerated contributions:* "(1) a corrected torque model whose length offset, effective
   stiffness, and fictive-strain wrapping loss — expressed through the unloaded 620 kPa
   contraction — halve prediction error on held-out tests across configurations; (2) a
   cross-configuration-consistent stiffness pair selected by a simultaneous design requirement;
   and (3) redesigned attachment routes that meet or exceed human isometric knee-torque targets
   with quantified margins, confirmed by first hardware measurements."
5. *Implication:* the corrected model and redesign pipeline support first-iteration torque design
   of BPA-actuated robot joints.

**Keywords (8):** braided pneumatic actuator; artificial muscle; isometric torque; moment arm;
system identification; actuator routing; biomimetic robotics; knee joint.

---

## 5. Tables list (LaTeX float in Main + .xlsx export into `Knee_Torque_Test/Tables/`)

All numbers from `mining/results_ledger.md` (confidence DR/DR-verified) or the chapter as noted.
Booktabs, caption ABOVE the table, full-sentence captions defining abbreviations and units
("RMSE and maximum absolute residual (Max.) are in N·m"; "FVU denotes the fraction of variance
unexplained"). Coefficients carry `\qty`. The existing `tab:OptCoeff` table is RETIRED — its
five initial coefficients move to Results running text (Actuators Pattern A); the adopted-values
table T3 supersedes it. (Ben may override — §9.)

| # | Label | Caption / name | Columns | Rows | Ledger source | .xlsx filename |
|---|---|---|---|---|---|---|
| T1 | `tab:FlxPinGof` | Pinned-knee flexor torque-prediction goodness of fit per test | Test (l0 + role); RMSE/FVU/Max under "Improved BPA characterization alone"; RMSE/FVU/Max under "Corrected model" | 5 tests + mean (training) + mean (validation) | §3 flexor table (= chapter tab:FlxPinGof :373–380) | `PinnedFlexor_GoF.xlsx` |
| T2 | `tab:ExtPinGof` | Pinned-knee extensor torque-prediction goodness of fit per test (stiffness pair locked to flexor) | same six-metric layout | 6 pool tests + mean (fitting) + mean (validation) | §3 extensor table (= chapter tab:ExtPinGof :409–417) | `PinnedExtensor_GoF.xlsx` |
| T3 | `tab:AdoptedXi` | Adopted correction terms of the two-stage identification | Parameter; value (`\qty`); search bounds; configuration/role | χ0/χ1/χ2 (flexor), χ0/χ3 (extensor) — 5 rows; caption carries the locked-pair + external-selection statements | §1 (adopted + bounds rows) | `AdoptedCorrectionTerms.xlsx` |
| T4 | `tab:Redesign` | Route-redesign outcomes against human isometric torque targets | Route variant; predicted torque span; human-scaled requirement; minimum margin (%); worst angle; distinguishing note | direct extensor; cam extensor; direct flexor; mono-pulley flexor | §6.1–6.4 | `RouteRedesign_Summary.xlsx` |
| T5 | `tab:BiPulley` | Bi-pulley batch on five muscle sites | Muscle site; objective; outcome (flag); minimum route clearance (mm) | 5 sites; caption states the 25 mm clearance requirement and the 101 sentinel | §6.5 + chapter tab:bipulley_batch :520–534 | `BiPulley_Batch.xlsx` |

Placement: T1–T3 in Discussion 4.2; T4–T5 in Discussion 4.3. Export script lives in
`tmp/paper_bridge/tabscripts/` (openpyxl via `D:/Anaconda/python.exe`); ONLY the .xlsx files are
created in `Tables/` (Ben's `CoefficientTable.xlsx` untouched).

---

## 6. Figures list (keep / regenerate / new)

Master conventions (BRIEF + `Plot_KneeFlx_all.m`): Arial, FontSize 12 bold axes, LineWidth 2,
minor ticks both axes, TickLength [0.025 0.05], legend 8 pt, GridLineStyle none, panels (A)–(D)
12 pt bold Arial EdgeColor none, shared axes labels, EPS into `ProofFinal/Figures/` under the
SAME filename when replacing. Palette (`PlotPrompt_utf8.txt` hex list): measured = indigo
`#0000FF` filled scatter (MarkerFaceAlpha 0.75, sz 60); hybrid = gold `#FFD700`; original
uncorrected = gray `[0.4 0.4 0.4]` dashed LineWidth 2; corrected/optimized prediction =
magenta `#CD34B5` LineWidth 2.5; improved-BPA-characterization-only = light orange `#FA8775`
(chapter FlxGrp convention). QC: every result-figure caption states the sign convention
("negative values denote knee flexion" on flexor figures; the positive-extension analogue on
extensor figures); no hue-only legend separation; legible at ~6.2 in width.

**KEEP (3):**
| File | Where | Notes |
|---|---|---|
| `Figures/testJigs.pdf` | Methods 2.3 | unchanged; caption gains nothing (photo panel labels already there) |
| `Figures/KneeICR.eps` | Methods 2.3 | unchanged |
| `Figures/bktFrame.pdf` | Methods 2.6 | unchanged file; caption panel-D mislabel "(A)"→"(D)" fixed in .tex |

**COPY, NEW TO PAPER (2):** (owner: figures agent; fallback assembler)
| File (new) | Source (read-only copy) | Where | Notes |
|---|---|---|---|
| `Figures/xiFrameTransform.pdf` | `Dissertation/Figures/30-torque/xiFrameTransform.pdf` | Methods 2.6 | two-rotation bracket-frame construction; keep its alt text spirit in the caption |
| `Figures/ExtensorRoute_redesigned.pdf` | `Dissertation/Figures/30-torque/route_geometry_ext20mm_20260925.pdf` | Discussion 4.3 | RENAMED copy (no dates in paper artifacts); eight-panel femur-frame view of the optimized extensor route; text must not reference the source filename |

**REGENERATE (4)** (stale result figures; same filenames; data from `Testing_Data/2022_02_Festo`
with `Code/Matlab/Functions` + `Code/Matlab/Robot_Data` on path; scripts in
`tmp/paper_bridge/figscripts/`; `matlab -batch`):
| File | Spec | Data / model source | Risk |
|---|---|---|---|
| `FlxPin_group.eps` | 2×2 tiles; columns titled "Optimization" (48.5 cm) / "Validation" (45.7 cm); row labels "Pre-optimization" / "Optimized"; legends in tiles 2 and 4; row 1 = Hybrid (gold), Measured (indigo), Original (gray dashed), Improved BPA model (light orange); row 2 = Measured, Improved BPA model, Optimized prediction (magenta 2.5); XLim [-120 20] (PlotPrompt final paragraph overrides the [-120 15] master default here), YLim [-35 0], XTick [-120 -90 -60 -30 0] | run `minimizeFlxPin(0.0119, 1.02e5, 1.06e4)` for the corrected series + zero-correction + hybrid series per the PlotPrompt's bpa-structure usage; raw tests via `FlxPinBPASet.mat`/raw tabs | moderate: initial-era bpa structs easiest obtained by calling the evaluator at the Table-1 coefficient values (PlotPrompt pattern) |
| `ExtPin_group.eps` | 3 panels (A)(B)(C) per PlotPrompt ext layout (tileIdxs [1 5 10], tileOrder [4 3 1]); same palette; YLim [0 15]; panel titles use the TRUE resting lengths of the plotted tests (ledger §4 values; not the legacy rounded strings) | run `minimizeExtX3(-0.0109, 1.02e5, 1.06e4, 0.123, idx)` + zero-correction + hybrid | moderate |
| `Ext10mm_52cm.eps` | single panel, same four series (measured indigo, hybrid gold, original gray dashed, optimized magenta); positive torque; caption label l0 = 51.8 cm | bio-extensor prediction path in `Testing_Data/2022_02_Festo` (locate the evaluator section that produced the chapter's fig:Ext10mm_52cm); raw tab `ExtTest10mm_*` era data | if the producing code cannot be located unambiguously, SKIP per the mission (old figure stays; assembler informed) |
| `FlxGrp.eps` | panel A: ⌀10 mm, l0 41.5 cm + 21 mm tendon @ 604 kPa — Measured (indigo), corrected prediction (magenta 2.5), characterization-only (light orange dashed); panel B: ⌀20 mm @ 620 kPa — Measured (indigo), corrected (magenta), human maximum-isometric line (black dotted), sub-620 measured points (pink dots, 560/421/325/281 kPa); sign-convention sentence in caption; length label per §8 policy (default 42.0 cm + 15 mm) | `Results_table_20mm.xlsx` tab `FlxTest20mm_42cm (3)` (RML 415 mm — FLAG-3) + `Knee_Flexor_Data_20mm.m`-era prediction machinery + OpenSim human-target export | highest: mixed provenance (tab vs chapter labels) — if ambiguous, SKIP and flag (Ben decides label first) |

**DROP from main text (1):** `KneeExt_10mm_passthru.png` — remains in the supplement (same image,
supp fig `fig:passThru`); Discussion 4.1 cites it as the supplementary figure. Frees one display
slot for `xiFrameTransform`. (If Ben prefers it in the main text, the budget is exactly 15 — §9.)

**OUT OF SCOPE:** the PlotPrompt's muscle-length / moment-arm / relative-strain working figures
(Ben's analysis plots, not paper display items); orphaned files (`F-P_contraction.eps`,
`FlexText4_pinned_FishScale_leak.png`, `FlxPin_angleMisalign.png`, `Flx_brktBend.png`,
`Gait2392_rendered.PNG`) stay unreferenced; `Figures/Ext10mm_52cm.eps` KEEPS its filename
(52 cm legacy name) even though text says 51.8 cm — filename stability over label accuracy.

---

## 7. Citation strategy and gaps

**Existing keys to deploy** (all 62 keys verified present; 35 currently uncited — headroom):
Intro tour + musculoskeletal legs (kurumaya_musculoskeletal_2016, hitzmann_anthropomorphic_2018,
asano_musculoskeletal_2019, liu_robotic_2018, he_mechanism_2020, chen_realizing_2019,
leibach_development_2020, mahapatra_multi-legged_2020, romero_four-bar_2020, veale_towards_2016,
madden_artificial_2004, suzumori_trends_2018); BPA modeling comparisons in 4.4
(martens_modeling_2017, tondu_static_2022, sarosi_new_2012, sarosi_comparative_2017,
sarosi_investigation_2013, chou_measurement_1996, aksoz_design_2019, al-rawi_joint_1985,
andrews_general_1985, scharzenberger_design_2019, li_influence_2020, liang_comparative_2020);
moment-arm methodology (menegaldo_moment_2004, brand_model_1982, suderman_moving_2012); human
targets (revive the commented Main.tex:225 block: hoy_musculoskeletal_1990, yamaguchi_planar_1989,
sherman_what_2013, millard_flexing_2013, thelen_adjustment_2003, seth_opensim_2011);
muscle-model context where needed (arnold_model_2010, zajac_determining_1989); manufacturer data
with page pointers where pressures/limits are invoked (festo_2026 — Actuators house precedent).

**Citation GAPS (report via missingCitations; bibliographer verifies before adding; never
invent; fallback = generic algorithm naming + the existing MATLAB RRIDs):**
1. **NSGA-II** — claim: the identification uses an elitist nondominated-sorting genetic algorithm
   (gamultiobj). Suggested: Deb, Pratap, Agarwal, Meyarivan (2002), IEEE Trans. Evol. Comput.
   6(2):182–197, doi:10.1109/4235.996017.
2. **Surrogate-based optimization** — claim: the redesign uses radial-basis-function surrogate
   optimization. Suggested: Jones, Schonlau, Welch (1998), J. Global Optim. 13:455–492,
   doi:10.1023/A:1008306431147.
3. **Pattern search / mesh-adaptive direct search** — claim: pattern-search refinement of the
   surrogate winner. Suggested: Audet, Dennis (2006), SIAM J. Optim. 17(1):188–217,
   doi:10.1137/040610330.
4. (OPTIONAL) **Human isometric torque literature context** for §3.4/4.3 — the OpenSim-derived
   targets already stand on delp_interactive_1990 + seth_opensim_2018; add only if Ben wants an
   experimental-torque cross-reference.

**Bib hygiene:** `bolen_2026` is an `@unpublished` entry (bib L1085) — it must be updated to the
published Actuators 15(5):230 form before submission (flagged, §9). Both bib copies must stay in
sync (currently identical modulo line endings).

---

## 8. Number-of-record policy (binding for all writers; never harmonize silently)

1. **Results quotes the INITIAL identification** (chapter :348, :396, :430, :451 = ledger §2);
   **Discussion quotes the ADOPTED two-stage set** (ledger §1/§3). The biomimetic-extensor test
   legitimately carries BOTH GoF pairs (initial 0.8920/0.3396/1.4188 in Results; adopted
   1.17/0.67/2.21 in Discussion 4.2d) — label each by its parameter set.
2. Superseded exploration sets (pk107lock, row-107 alternates — ledger FLAG-6) are never quoted.
3. Training-mean flexor RMSE is the 4-test mean 8.122 → 1.557 (FLAG-4); state the set averaged.
4. Extensor resting length: 51.8 cm everywhere in text/captions (Main's 52.0 fixed; the FILE
   `Ext10mm_52cm.eps` keeps its name).
5. 20 mm biomimetic flexor label: default = chapter text of record 42.0 cm + 15 mm tendon
   (open question 4; the plotted tab records RML 415 mm/14 mm — unresolved).
6. Cam-run surrogate budget is NOT stated (FLAG-5); budgets only for extensor 3,000 / flexor
   1,000 / pattern search 15,000.
7. First hardware points (1.163/7.849/23.437 N·m) quote the chapter values; the lab tab's row-16
   direct-formula variant (~12% lower) is NOT quoted (FLAG-1; open question 3).
8. Whole-body 27-actuator numbers are chapter-record only (ledger §6.6) — quoted as honest
   status with no artifact-level precision beyond the chapter.
9. Main.tex stops referencing the supplementary `tab:TorqueFit` (its numbers disagree with the
   record — paper_state §3.7); the new main-text tables T1/T2 + Results text carry the GoF.
10. Every number in every fragment needs its named source (ledger row or chapter line) in the
    writer's returned numberSources — Gate 3 samples ten of them.

---

## 9. OPEN QUESTIONS for Ben (decisions only Ben can make; defaults stated)

1. **Author list:** Main.tex lists 4 authors; the dissertation publication note
   (30-torque.tex:5) names 7 contributors for bolen_isometric_2025 (Elzein, Pang, Scharzenberger
   missing). Do not change until Ben rules. (BRIEF flag; standing.)
2. **Author Contributions defect** (Main.tex:490): "Investigation, B.B., C.M., and C.M." — C.M.
   duplicated, L.B. absent. Needs Ben's intent for the CRediT split; left untouched.
3. **First hardware torque points (FLAG-1):** the chapter's 1.163/7.849/23.437 N·m come from the
   Adjoint-transform section whose rig constants (d = 292.9 mm, angle −90.83°) are inherited from
   the 40 mm test with an in-code "UPDATE after rig measurement" note; the lab tab's row-16
   variant reads ~12% lower and row 17 is empty. Confirm the 20 mm rig constants (or have them
   measured) before submission. Default: quote the chapter values.
4. **20 mm biomimetic flexor label (FLAG-3):** chapter text/figure say l0 = 42.0 cm + 15 mm
   tendon; `tab:TorqueFit` row says 41.5 cm; the plotted lab tab records RML 415 mm + 14 mm.
   Which label is the value of record? Default: 42.0 cm + 15 mm (chapter text). The FlxGrp
   regeneration needs this before it runs.
5. **Supplementary GoF table refresh:** the untouched supp `tab:TorqueFit` disagrees with the
   values of record. Options: refresh it before submission, or retire it (Main no longer cites
   it). Default: leave untouched now; Ben decides pre-submission.
6. **New algorithm citations** (Deb 2002 / Jones 1998 / Audet & Dennis 2006 — §7): add via the
   bibliographer, or keep generic algorithm naming + RRIDs? Default: request all three.
7. **Retire `tab:OptCoeff`?** Plan folds the initial coefficients into Results running text
   (Actuators Pattern A) and lets T3 carry parameter values. Keep the table instead if Ben
   prefers (budget goes to exactly 15/15).
8. **Pass-through figure location:** plan moves `KneeExt_10mm_passthru.png` out of the main text
   (it remains as the supplementary figure) to make room for `xiFrameTransform` within the
   15-item cap. Confirm.
9. **CUT content:** adductor-magnus P2 restoration (robot-body bookkeeping) is cut from the
   paper; the whole-body 27-actuator extension stays as one honest-status paragraph. Confirm.
10. **`bolen_2026` bib entry** is `@unpublished`; update to the published Actuators 15(5):230
    form before submission (key stays).

---

## 10. Risks and fallbacks (for downstream agents)

- **TeX missing** (`D:\MiKTeX\miktex\bin\x64\pdflatex.exe` not present at mining time): every
  compile is PENDING per the BRIEF; never treat as manuscript failure. Compile the supplementary
  first if xr complains.
- **siunitx v3 vs the `\DeclareUSUnit/\US` v2 shims** (Main.tex:21–24): verify on first compile;
  the class already loads times/caption/etc. — add packages with NO option clashes; never add
  geometry; mathabx is the one sanctioned addition.
- **Figure regeneration skips:** if a figure's data source is ambiguous (esp. FlxGrp, FLAG-3;
  Ext10mm_52cm producer code), SKIP per mission and keep the old file — the assembler keeps the
  existing `\includegraphics` and says so.
- **Evaluator scripts may write files despite no save calls** (ledger session-integrity note:
  `Knee_Extensor_20mm.m` regenerated CSV/GIF via writetable/imwrite) — the figures agent greps
  for writers in any script it runs and restores tracked byproducts (`git checkout --`) if hit.
- **Label collisions:** fragments bring chapter-style labels; the assembler must dedupe against
  surviving Main labels (keep one label per object; check refs resolve).
- **Line numbers** (`\linenumbers`) stay on; Frontiers trailing Figures/Tables stub sections stay;
  tables are real LaTeX, never images.
