# Dissertation evidence inventory from the 2026-10-24/25 desktop sessions (ZCode, EB475WS4)

## Rulings recorded 2026-10-26 (Ben)

- **Terminology**: integrator/solver steps are **time steps** (never
  "steps" alone); gait cycles use **stride steps** or **foot steps**.
  Applies to all prose, figures, decks.
- **Figure rejections + v3 spec** (all three demos): knee reflex = needs
  kinematic-response overlay + comparison (ON vs OFF), muscles routed
  clear of the thigh; knee CPG = kinematics + neural architecture +
  four-corner half-center diagram; beer = humerus vertical, forearm
  horizontal, natural cup-hold posture, reflex-pathway panel. Implemented
  in `demos\sns_demo_figures_v3.m` → staged in
  `Figures\60-sns\components\fig_{kneereflex,kneecpg,beer}_v3.png`;
  animations regenerated to the same anatomy.
- **Advisor name**: decks carried "M. Anthony Hunt" — correct is
  **Alexander J Hunt** (fixed where not file-locked).
- **Deng 2023 mined** (kaiyu2023_fulltext.txt; findings
  `dev\deng2023_mine_findings.txt`): 8 muscles/limb with 3 biarticulars
  (BFP/RF/GA focus), AnimatLab 2, conductance-based non-spiking neurons
  with persistent Na, piecewise-linear synapses (per-synapse g, E_lo,
  E_hi), weak cross-excitatory RG connectivity, two-layer synapse
  strength trade-off, 2 nA perturbation protocol — the paper whose
  perturbation panels easteregg2 reproduced in `sec:w2l_perturbation`.

STATUS: evidence + figures STAGED; ZERO prose edits to Revision_D or any
.tex. This file is the compiled handoff — each item carries its verified
numbers, its evidence path, and a PROPOSED dissertation placement with a
draft two-sentence result. Folding any of these into Revision_D goes
through the campaign's supervisor gates (3-gate contract, Ben 2026-10-24)
and the no-silent-number-changes rule. Chapter numbering per the current
structure (60-sns etc.).

## A. Platform-implementation evidence (SNS ch. 60)

**A1. Deng 2019 CPG rebuilt + verified in base-Simulink (R2025a).**
Numbers: 10 bursts from one 10 nA/1 ms pulse; period 1.938 s (numpy
reference 1.94 s); ext-vs-flx correlation −0.857 (documented −0.86).
Evidence: demos\sns_run_deng_demo.m run log 2026-10-24; rebuilt trio
.slsx + rewritten sns_build_deng_cpg.m.
Draft result: "The Deng two-layer CPG was rebuilt in Simulink with the
persistent-sodium half-centers (fixed τh = 350 ms) reproducing the
self-sustained 1.94 s rhythm from a single current pulse, matching the
numpy reference period to within 0.01 s and the antagonist correlation
to 0.003."

**A2. Li two-layer contact-driven CPG built from the AnimatLab mine.**
All 58 nodes / 98 connections, every synapse carrying its own per-synapse
conductance (1–8 µS). Smoke: both stance CPGs oscillate ~2.0 s period,
inter-leg correlation −0.271 (antiphase tendency) under scripted antiphase
heel contacts. LABELED adaptation: graded stance cells + slow adaptation
(tau = documented tau_Ca 7.5 s), not spiking I&F.
Evidence: demos\SNS_Li_CPG.slx, sns_build_li.m, dev\run_li_smoke_20261024.m.
Draft result: needs the adaptation labeled, so: "A Simulink transcription
of the Li contact-driven CPG preserved the complete connectome and
per-synapse conductances; with graded stance cells substituting the
original spiking dynamics (documented substitution), the network produces
rhythmic per-leg stepping with contralateral antiphase under alternating
contact input."

**A3. Toolbox tutorial twins (implementation-equivalence evidence).**
Tutorial 2 (5-neuron exc/inh/modulatory): all five finals match the
sns-toolbox numpy backend to ≤1.1e-05 mV. Tutorial 8 (persistent-Na pair):
winner/loser steady states exact (20.000 / −0.010 mV), loser pointwise
0.595 mV, winner excursion within 7% and 5 ms.
Evidence: dev\sns_tutorial2_twin_20261025.m, dev\sns_tutorial8_twin_20261024.m,
results\tutorial2_ref.mat, tutorial8_ref.mat (+probes).
Draft result: "Simulink twins of the toolbox tutorial circuits agree with
the reference numpy backend to microvolt precision (five-neuron network)
and reproduce the persistent-sodium pair's winner-take-all steady states
exactly."

**A4. Full-network verification (existing, reproduced on desktop).**
410-neuron production export: max deviation 4.201e-06 mV vs numpy, all
neurons at t = 0.3 s; winner slice all 4 checks. Evidence: capability
report + logs (2026-10-24).

## B. Per-target conductance evidence (SNS ch. 60 + methods)

**B1. Deng 2019's own AnimatLab project (ModelDB 264844) mined.**
90 non-spiking neurons, 196 synapses, 58 named types with distinct
per-path conductances: PF→MN per joint = hip flx 3.6325 / hip ext 2.6563
/ knee flx 1.5165 / knee ext 4.9296 / ankle flx 4.5222 / ankle ext
4.0538 µS. 114/196 edges carry the untouched 0.5 µS per-instance GUI
default → the effective strength in AnimatLab projects is the TYPE
conductance; mines now record both.
Evidence: dev\deng2019_mine_20261025.json, mine + check scripts.
Draft result (methods-flavored): "The original Deng project assigns each
pathway a dedicated synapse type with its own conductance; per-instance
values retain the uniform editor default, establishing type-level
conductances as the published model's effective strengths."

**B2. Twin updated to the paper's conductances and re-verified** (A1
numbers unchanged; hip PF→MN now 3.6325/2.6563 µS from the paper, was
Nourse 3.632/2.565).

## C. Curriculum evidence (SNS ch. 60, the five-part curriculum)

**C1. Air stepping (stages 1–2), first Simulink instance: PASS.**
Two BPA-actuated planar legs (8× BPA_20mm, Festo surfaces), Deng RG + PF
layers, pelvis fixed, no contact: hips ±47°, knees ±46°, L-vs-R hip
correlation −0.411 (antiphase). Labels: single-shared-CPG v1;
decoupled-joint mechanics; swing frequency currently pendulum-dominated
(~0.31 s) vs RG 1.938 s (tuning open).
Evidence: demos\SNS_AirStep.slx, sns_build_air_step.m,
dev\run_air_step_20261025.m.
Draft result: "Under the staged validation curriculum, the first two
stages (air stepping) were demonstrated in Simulink: a synthetic
rhythm-generator network drove both BPA-actuated legs in sustained
alternating swing with contralateral antiphase (hip correlation −0.41)
while unloaded."

**C2. Literature anchors staged (NOT citation-ready).** Reflex stiffness
share 30–60% (Nichols & Houk 1976; Sinkjaer 1988; Mirbagheri 2000) vs
our beer demo's ~74%; Thelen 1984/1986 air-stepping; Lokomat/BWSTT
30–50% initial support. All web-sourced, full texts UNVERIFIED.
Evidence: reports_20261024\literature_values_reflex_curriculum.md.

## D. Figures (staged, NOT placed into the dissertation tree)

**D1. v2 leg + cup animations/stills** — anatomically shaped limbs,
braided McKibben BPA bellies, level beer surface in a tilted pint glass.
Files: results\animations\v2\*.gif, results\pictures\v2\*.png (4 stills).
PENDING: Ben's visual approval, then copy into Figures\60-sns\components\
+ captions, under the style contract (Tol palette, Arial, 10 pt floor,
alt text).

**D2. Deng demo figure** regenerated (sns_deng_demo.png, 1.938 s run).

## E. Infrastructure results (SNS ch. 60 cross-platform subsection)

**E1. MuJoCo↔Simulink bridge built on EB475WS4** (R2025a + MinGW64, 4/4
mex targets; bidirectional S-function). Bit-exact tests a/b/c NOT yet
re-run on this machine — required before any claim.
**E2. Timing:** pure-Python MuJoCo 2.3.7 steps a 28-DoF model at
~107,000 steps/s (3,400 steps/s at 211 DoF); the Simulink bridge's 2 s
example run costs 46.2 s wall (compile-dominated). Bridge = integration
vehicle, not throughput.
Evidence: bridge build log + timing scripts in mujoco_bridge\matlab\.

## Not done / blockers

- **Table 6.1 RULING (Ben, 2026-10-25): the table "sucks and is
  confusing" — the section should be restructured around KNOWLEDGE GAPS
  and what this work did to FILL them, not around the five-task matrix.**
  easteregg2's campaign owns the restructure (its
  `B2_table61_explainer.md` draft); this inventory's items A–E slot into
  that gaps-and-fills frame as the fill evidence.
- NO prose in any chapter; NO compile (no source changed).
- Table 6.1 disambiguation still awaits Ben (his "five-part curriculum"
  definition is recorded; the table-vs-curriculum mapping question is
  open).
- Deng 2023 paper (Ben's links) not yet downloaded — needs a
  logged-in browser; then mine like B1.
- Bridge bit-exact tests; AnimatLab headless run of the paper aproj;
- easteregg2 campaign (arbitration, BW, BC stages) running separately —
  its own report carries its dissertation rows.
