# Spinal network for Gait2392 in MuJoCo (SNS-Toolbox) — design & status

Built 2026-09-09. Files in `Code\MuJoCo_SNS\spinal\`.

## 2026-09-17 VERIFICATION (EB475WS4: model, figures, and draft audited)

The documented local interpreter works:
`C:\Users\Ben Bolen\.conda\envs\myo\python.exe` (Python 3.10.21,
NumPy 1.22.4, SciPy 1.9.3, MuJoCo 2.3.7, SNS-Toolbox 1.5.2).
`_fix_check.py`, `_panels_check.py`, `audit_signs.py`, and the new
`_live_plant_audit.py` pass. `_smoke_nap.py` ran the full patched
MuJoCo--SNS loop at 2 ms and passed: finite dynamics, five RG bursts per
leg in the 5--15 s walk window, 2.11 s period, E-duty 0.67, right knee
-106.3..+11.6 deg, right hip -20.1..+50.9 deg. `spinal_run.npz/.png`
now contain that verified air smoke run.

Figure audit found and fixed four concrete problems:

1. `runner.py` plotted radian joint data under a degree label and used
   pre-NaP hard-coded neural column positions (including PF columns for
   the left RG summary). Plots and summaries now resolve channels by
   `NEURO_NAMES`, plot degrees, and restrict rhythm metrics to the true
   5--15 s walk window.
2. `_render_panels.py` still drew the retired tonic DRIVE->PF edge and
   called a stage-3 `score=-100` sentinel a tuned winner. The edge is
   removed; failed stages are ignored and representative gains are
   labelled as such.
3. `draw_circuit.py` checked only compiled-but-undrawn edge groups, so
   extra fictional arrows could pass. The contract is now two-way. It
   immediately caught the V3 target drawn/labeled as contra RG-E; the
   compiled circuit targets contra InE. The figure and conductance label
   are corrected, and the missing/extra edge-group contract passes.
4. `_figure_hindlimb_style.py` now resolves RG/PF columns from the NPZ's
   `neuro_names` instead of positional constants.

The live Fmax repair is confirmed: bilateral ercspn/intobl/extobl/ext_hal
are exactly 2500/900/900/162 N after `apply_harness`. The broader draft
claim was too strong: among 86 non-pruned active actuators, 78 are within
15% of stock OpenSim and eight exceed 15% (worst gluteal paths about
31--34%). `CPG_spinal_section_draft.tex` now reports that result, the
verified air metrics (four complete intervals/five onsets, 2.108 s,
cycle duty 0.691), and the honest curriculum state: stage 1 valid;
stages 2--3 remain `-100` sentinel diagnostics, not winners.

Regenerated and synchronized: `circuit_dengstyle.{pdf,svg,png}`,
`circuit_literature.{pdf,svg,png}`, `sns_diagram_panels.{pdf,png}`, and
`hindlimb_style_nap_air.{pdf,png}` in the
spinal `figures/` and dissertation `CPG_airstepping_figs/` folders.

Visual follow-up against the local historical AnimatLab figures and the
Shinohara-2025, Shevtsova-2026, Rybak-2024/2025, Deng-2019, Di-Russo-2023,
Jankowska, Rahmati, and Klishko PDFs: the detailed circuit is retained as the
zoomable compiled-edge audit figure; the toolbox-native figure is now a
page-readable A/B (RG/PF) over C (representative knee motor/reflex circuit)
composition with duplicate rendering ports removed. Air traces show only the
settled 6--15 s commanded-walk interval. All current figures now follow the
recent literature convention extensor=blue, flexor=vermillion/orange. This is
presentation-only. Final source/dissertation copies are hash-identical.
After Ben's annotated visual review, the preliminary dissertation draft now
uses `circuit_literature.pdf` first: a bilateral RG--PF--MN--muscle hierarchy
with centered c1/V3 relays, plus one complete antagonist/reflex knee motif.
The toolbox panels remain a software-native artifact. `circuit_dengstyle.pdf`
remains on a dedicated page as the compiled-edge audit.

The annotated review exposed genuine audit-figure defects: orphaned
hard-coded MN-to-activation starts, no activation-map-to-flexor-muscle edge,
and no flexor muscle-to-Ia/II/Ib encoder paths. Heel/toe-to-InE and
PF-F1-to-IaIN also used unconnected raw start coordinates. These are fixed.
The dense audit now shows complete extensor and flexor interfaces, omits
partial hip glyphs, and separately asserts the plant-interface paths because
they do not exist in SNS-Toolbox's `net.connections`. The new reader figure
asserts every displayed solid neural edge against a freshly compiled network.
Its final visual pass also restored RC->IaIN recurrent disinhibition, RG-E
stance gating of IB-EXC, and aggregate sensory projections to both RG and PF.

### 2026-09-17 directional commissural + activation/FSA correction

Ben caught one remaining literature-figure abstraction error. The compiled
network already had four distinct directional commissural cells, but the
reader figure collapsed them into a shared C1 and a shared V3. The corrected
paths are:

- RG-E_l -> V3_l-to-r -> InE_r and RG-E_r -> V3_r-to-l -> InE_l
  (both synapses excitatory);
- RG-F_l -> C1_l-to-r -| RG-F_r and RG-F_r -> C1_r-to-l -| RG-F_l.

`draw_literature_circuit.py` now draws the four cells separately and
`_fix_check.py` asserts every compiled source/sign/destination. The
representative diagram build sets `rg_weak_exc=0`; direct ipsilateral
RG-E<->RG-F excitation is absent from the audited graph.

Activation provenance was also corrected. The surviving
`ResultsBSolve/zz_bsolve_StaticOptimization_activation.sto` does not match
`bsolve_out.npz['acts']` (RMSE 0.3924, correlation -0.0119 on common
frames). Therefore the circuit backsolve uses the unambiguous converted-
MuJoCo ridge/NNLS `acts` array. Earlier four/five-component NMF values from
the STO do not validate the converted-model PF count.

`fsa_backsolve.py` implements the nonspiking leaky-integrator Function
Subnetwork Approach target mapping (`V_MN=5 mV*activation`) and exports the
rank scan, PF time courses/muscle weights, MN reconstructions, and inverse PF
current requirements to `fsa_results/`. Held-out interleaved-frame VAF at
four PFs is only 0.852 right/0.823 left. Six is the smallest shared count
above 0.90 (0.924 right/0.912 left); in-sample centered VAF is 0.926/0.913.
Ten-seed robustness is tight in VAF, but one right sixth component is
spatially unstable and the record contains only about 1.6 cycles. Treat six
as the engineering target for this dataset, not a universal biological
primitive count, until longer/multiple trials confirm it.

The dynamic excitatory PF->MN fit reaches centered VAF 0.863/0.835. About
16--17% of MN samples require net negative current, so higher-fidelity closure
needs explicit phase-specific inhibition (or an equivalent negative-current
pathway), not only excitatory PF conductances plus membrane leak.

The OpenSim-4.6 `normal.mot` experiment is quarantined from this inference.
It required an analysis-only model derivative disabling failed-equilibrium
`lat_gas_r` plus six pelvis residuals; the vertical residual averaged
739.4 N. It is a reproducible kinematics-only, residual-supported solve in
`ResultsNormalSO/`, not a ground-contact walking activation target.

### 2026-09-17 phase/PF semantics + ankle-capacity audit

Do not conflate the existing E1/E2/F1/F2 early/late stance/swing
pattern-formation channels with the six NMF synergies. Figures call the
latter S1--S6 only. A synergy is an observed temporal coefficient plus muscle
weight vector; it is not a PF neuron or layer. A joint/functional PF layer
may contain extensor and flexor half-centers whose outputs express two
synergy-like patterns, while one synergy may combine multiple PF layers.
MN pools may receive convergent drive from several layers; this is expected
to matter most for biarticular hip+knee and knee+ankle muscles. The current
FSA fit already permits all six unconstrained sources to project to every MN,
but it does not identify an anatomical PF-layer decomposition.

Their stance-rescaled peaks are approximately S1=96%, S2=40%, S3=10%,
S4=82%, while S5/S6 are not bilaterally phase-stable. Next interpretation
work should test S3<->S6 weight reassignment with a PF-informed constrained
factorization, not hand edits, and overlay synergy/PF contributions on both
the early-stance and swing-phase knee-flexion episodes (the classic
double-knee pattern). Preserve monoarticular/biarticular identity so the same
MN pool can receive appropriate hip+knee or knee+ankle PF drive.

`gait_phase.py` uses measured GRF heel strike/toe-off and maps stance to
0--50%, swing to 50--100%. The current 2-s activation record contains only
one complete stride per side (measured duty 0.622/0.615), so these are
cycle-normalized traces, not a multi-cycle statistical average.
`plot_gait_joint_angles.py` adds the requested 5x2 OpenSim-coordinate
kinematic figure (sagittal column; YZ-plane column with intentionally blank
knee/MTP panels).

`compare_ankle_df.py` compares Ben's OpenSim right-dorsiflexor force/torque
sweeps with the patched converted MuJoCo model at activation 1 and zero
velocity. In the IK gait range (-8.84..+16.02 deg), MJ/OS force ratios are
1.007 ext_dig, 0.966 ext_hal, 0.945 per_tert, 1.033 tib_ant; torque ratios
are 1.021, 0.974, 0.973, and 1.041. Summed torque is 1.023x OpenSim.
Therefore static force/moment-arm capacity is not the dorsiflexion deficit;
inspect recruitment timing/magnitude, antagonist activity, and dynamic
force-velocity effects. As at the knee, use equality-aware FD tendon moments:
raw MuJoCo `actuator_moment` is zero for these converted ankle paths.

## 2026-09-16 NIGHT (CURRICULUM STATE + INTERLEG LATCH DIAGNOSIS —
## READ FIRST if resuming; handed to ChatGPT for off-peak work)

**Stage 1 DONE: best 137.254 @ trial 78** (winner in
curriculum_stage1.json; 100 trials; ADAP-era all-time best was 57.4).
Winner is a FAST air rhythm (29 cycles / 0.37 s / duty 0.19 at
rg_nap_h=0.17, drive 3.15) — the air objective rewards rise count;
stage 2/3 must slow it toward the 1.2 s human cycle (rg_nap_h is
searched; stage 3's kine objective enforces cadence).

**Stages 2-3 RAN AND FLATLINED: best -100 @ trial 0 in BOTH (every
trial -100 = no-countable-cycles sentinel).** Two root causes, both
diagnosed empirically (isolation matrix below):

1. (FIXED) heel/toe->PF_E edges were wired in _build_rg BEFORE the PF
   cells exist -> hard crash on the first stage-2 trial with
   ib_e_central > 0. Moved into _build_pf. _fix_check.py now asserts
   every new central pathway at the compiled-net level (PASSES).
2. (KNOBS ADDED, TUNING PENDING) **interleg latch**: with interleg ON,
   the network bilateral-E-latches on ground (E-duty 0.90-1.00, knees
   pinned -4..+11, RG_F 0.11 s chatter, tilt 31) — IDENTICAL with
   v3_gain 0 and 0.5, and identical with afferents ON/OFF, so the
   driver is the c1 cross-side F-F inhibition at full rg_mutual_inh
   strength against NaP plateau neurons. Isolation matrix
   (_diag_gait.py, stage-1 winner config):
     air, no-interleg, deaff   : 29 cycles, 0.37 s, duty 0.19  (alive)
     air, no-interleg, afferents: 29 cycles, 0.37 s             (alive)
     ground, NO-interleg        : 29 cycles, 0.37 s, knee -85   (alive)
     ground, interleg (c1 1.0)  : LATCHED (v3 0 and 0.5 identical)
   FIXES SHIPPED: G["c1_gain"] (default 1.0) and G["v3_gain"]
   (default 0.0 = pathway absent) multiply rg_mutual_inh; V3 now
   targets the CONTRALATERAL InE (per Ben's Shinohara reading:
   "v3 -> excite -> contralateral RG_E IN" — V3->contra-RG_E direct
   was my error and is also a bilateral E-E positive loop); both are
   searched in stages 2-3 (c1 [0.1,1.5], v3 [0,0.5]). ALSO: the
   stage-2 objective was mislabeled from the start — it ran the GROUND
   eval; now stage 2 is truly air-afferented (--no-ground, interleg
   ON, air rises+knee objective) and only stage 3 runs the ground
   kine eval.

**RESUME RECIPE (next session, off-peak):**
1. Optional 30-s characterization: air + interleg c1=1.0 (the cancelled
   test) — tells you whether stage-2's seed region latches in air too.
2. Purge the garbage studies: curr_s2_air_aff, curr_s3_ground
   (all -100; 0 informative trials). KEEP curr_s1_air_deaff (complete).
3. Relaunch: `resume_curriculum_lam.bat` (stage 1 no-ops through with
   its saved winner, stages 2-3 fresh at 100 trials each).
4. On completion: `post_curriculum_deliverables.bat` (stage-3 winner
   22 s run + GIF + hindlimb + overlay + winner-gain figures).
5. Fill the `\fillme{}` slots in
   `Documentation\...\Dissertation\CPG_spinal_section_draft.tex`
   (copy-paste map in CPG_DISSERTATION_UPDATE_NOTES.md; the stage-1
   air numbers in PART B are already final).
6. Second commit (results), then push both commits (1fb515f = the
   architecture commit, local-only as of tonight).

## 2026-09-16 EVENING (NaP ARCHITECTURE OVERHAUL — LIVES NOW; answers
## the toolbox-correction section directly below)

The RG is **real persistent-Na conditional bursters in production**
(this section resolves the caveat in the "TOOLBOX CORRECTION" section
below — the real class WAS tested, in the real half-center):

- **Measurements** (`_nap_test.py`, `_nap_fixed_test.py`): the stock
  voltage-dependent tau_h(V) oscillates WITHOUT ADAP but ~10x too fast
  (period 0.05-0.12 s in every swept region); **FIXED tau_h** gives
  period 1.17 s at tau 0.35 s (~3.2x tau scaling; 0.70 -> 2.2-2.4 s).
- **Implementation**: `SNS_NumpyFixedTau(SNS_Numpy)` in build_network.py
  — verbatim forward() copy with ONE changed line (tau_b = tau_max_b);
  `compile()` swaps the class in. TRAP: sns_toolbox 1.5.2 puts
  everything in `SNS_Numpy.forward` (no `__forward_pass__`); a subclass
  override under the mangled name silently never runs (bit me once).
  Re-diff the forward body on any toolbox upgrade.
- **Params**: `params.NAP` (rescaled to our 0..5 mV range: e_ion 8 =
  plateau ceiling, e_m 2, s_h -2 / e_h 3.5) + `TAU["rg_nap_h"]` (0.35 s)
  = THE period knob, searched by the curriculum.

**Architecture changes shipped the same day (build_network.py):**
ADAP retired; PRESET/PREA pathway REMOVED (Ben circled it; replaced by
per-muscle afferent->central wiring); DRIVE->PF tonic edge removed
(audit #23); G_W weak mutual excitation RG-E<->RG-F (Deng 2022,
conditional); commissurals per Shinohara 2025 (RG_F->c1 IN->INHIBIT
contra RG_F; RG_E->V3 IN->EXCITE contra RG_E — sign flipped — plus
V3->contra IBEXC extensor MNs = audit P2a, G v3_to_ibexc); per-muscle
afferent->CENTRAL same-group excitation (Rybak 2025 SF-E1/SF-E2 +
Ben's spec): extensor Ib+II -> PF_E/RG_E/InE (ib_e_central,
ii_e_central), flexor Ia+II -> PF_F/RG_F/InF (ia_f_central,
ii_f_central), flexor Ia also INHIBITS contra RG_F (ia_f_contra_f),
heel/toe ride the extensor pathway. Runner: HIP_*_SIG ports removed
(NEURO_NAMES lost 4 channels — mind old npz neuro column indices).

**Untuned results** (drive 2.5, tau_h 0.35, `nap_air_walk.npz` /
air_nap.gif): 5 cycles, 2.11 s period (0.47 Hz), **E-duty 0.86-0.90** —
duty never exceeded 0.34 across all ADAP-era tuning. NaP plateau at
zero drive + POSTURE bias = tonic extensor stand state (Ben OK'd).
**Curriculum relaunched fresh on this architecture**: stage 1 search =
drive / rg_nap_h / desc_e / desc_f / rg_to_pf; 21 trials in, best
109.8 (ADAP-era stage-1 all-time best: 57.4).

**Figures regenerated + copied to Dissertation\CPG_airstepping_figs:**
circuit_dengstyle.{pdf,svg,png} (NaP RG, G_W, c1/V3 excitatory + P2a
arm, afferent->central arrows, HIP COLUMN added, no ADAP/PRESET/PREA/
DRIVE->PF; edge-contract assert passes), sns_diagram_panels.png
(spinal_layers.py rebuilt: NaP RG + InE/InF + G_W only),
hindlimb_style_nap_air.png, air_nap.gif. `_panels_check.py` = the
NaP-wiring assert suite. NOTE `_fix_check.py` predates the ADAP/PRESET
removal — refresh before rerunning it.

## 2026-09-16 (TOOLBOX CORRECTION, Ben caught it): persistent-Na class EXISTS

**`sns_toolbox` 1.5.2 ships `NonSpikingNeuronWithPersistentSodiumChannel`**
(Tutorial 8; constructor: membrane_capacitance, membrane_conductance,
g_ion, e_ion, k_m/slope_m/e_m, k_h/slope_h/e_h, tau_max_h, name, color —
Tutorial 8 was executed in `sns_tutorials` 2026-09-14). Earlier claims in
this file's history that "the toolbox can't express persistent-Na / NaP
dynamics, so the ADAP loop is the only burst-termination substitute" were
WRONG. Consequences:

- The RG half-centers CAN be rebuilt as literal Deng-2022-style HC neurons
  with intrinsic NaP burst termination (m/h gates, tau_max_h) instead of
  the ADAP-loop workaround — the cross-platform-consistency route, since
  the Simulink `SNS_Deng_Library.slx` (fixed tau_h 350 ms) and AnimatLab
  LinearHill (tau_h.max fixed) ports already work that way.
- CAVEAT before trusting the class: the quenching result (toolbox tau_h(V)
  collapses at depolarized V and kills the Deng oscillator;
  `deng_cpg_ode.py`) was measured on the HAND-CODED formula, not the
  class. Test the real class in the actual Deng circuit first — if its
  tau_h(V) collapses the same way, check whether fixed-tau_h semantics are
  expressible via tau_max_h before assuming the port transfers.
- When porting values: SNS_Library vs sns_toolbox use OPPOSITE synapse
  saturation conventions (ThrPre/Elo) — standing trap, still applies.
- Also logged in `Code\Matlab\SNS_Simscape\README_SNS_Simscape.md`
  (Deng section correction note) and AGENTS.md (spinal insights (5) +
  AnimatLab latch note).

## 2026-09-16 (CONNECTION FIXES — Ben: "the architecture is not correct,
## modify it" — READ THIS FIRST)

Two wiring bugs confirmed by the compiled-network audit (found while
verifying the diagrams against build_network.py) and **FIXED in
build_network.py**:

1. **HEEL pathway was wired TWICE** — `heel_in -> RG-E` (exc) and
   `heel_in -> RG-F` (inh) were each added twice, verbatim (duplicate
   lines 270-275), doubling the effective heel contact conductance.
   Now wired ONCE per target (`_fix_check.py` asserts count == 1).
2. **IaIN had NO Ia afferent input** — the docstring claimed
   "Ia -> IaIN -> antagonist MN" but the built IaIN received only the
   PF_F1 phase gate (0.5 exc) and RC->IaIN disinhibition; the antagonist
   pathway was swing-gated but not stretch-driven. FIX: Ia -> IaIN exc
   at G["ia_to_mn"] (same conductance as the homonymous Ia->MN arc; no
   new tunable knob). All 4 representative pools verified.

**Contamination consequences (important for re-picking):**
- Stage 1 (deafferented air) is UNAFFECTED: both fixed blocks live
  behind `stance_fb` / `ia_in` conditionals that are False at stage-1
  gains — the compiled stage-1 network is byte-identical pre/post fix
  (`_fix_check.py` section 3). The in-flight 100-trial laminated stage-1
  study stayed valid and was NOT re-run.
- Stages 2-3 imported the FIXED code (heel now 1x, IaIN afferent live),
  so their winners are the first tuned on correct wiring.
- ALL pre-2026-09-16 stage-2/3-equivalent tunings (v7-v11 studies,
  deleted 2026-09-16) carried the heel-2x bug + afferentless IaIN.

**Diagrams rebuilt to match the fixed architecture** (Ben's staleness
complaint confirmed: the Dissertation copies were from 09-12/09-13/09-15
and showed DIRECT RG<->RG + PF<->PF inhibition with RC/IaIN as
"NOT IMPLEMENTED" ghosts — the laminated inverse of the real circuit):
- `spinal_layers.py` (panels source) + `draw_circuit.py fig_deng` now
  draw the true laminated architecture with REAL InE/InF, PF_IN_E/F,
  mutual Renshaw, IaIN (afferent + gate + disinhibition), heel/toe/LBIN,
  PREA onset loops; fig_deng's edge contract (every compiled group must
  be drawn) PASSES; gain labels load live from the winner jsons.
  Network-size metadata corrected 410 -> 418 (laminated INs).
- OLD FIGURES/GIFS DELETED (Ben directive): all pre-fix run outputs from
  spinal root (7 gifs, fig1-6, spinal_circuit, check_rhythm,
  walk-phases pngs) and ALL stale SNS figures in
  Dissertation\CPG_airstepping_figs (deng/vclasses/weights/panels/
  spinal/hindlimb/overlay/gif sets) — only architecture-independent
  muscle_force_compare.png kept. New post-curriculum set regenerates
  into the same filenames (circuit_dengstyle, sns_diagram_panels,
  ground_v11_curriculum.gif, hindlimb_style_v11_curriculum,
  opensim_overlay_gait_cycles). NOT auto-regenerated (need fig_full
  laminated surgery if wanted back): circuit_full_vclasses.*,
  circuit_weights.*, sns_diagram_spinal.png, hindlimb_style_{air,ground,
  ground_pose}.
- Laminated curriculum re-run (Ben 2026-09-16): contaminated studies
  deleted (curr_s1/s2/s3 mixed + v7/v8/v9/v10; v1-v6 kept as pure
  pre-lamination history; per-trial histories dumped to
  `curriculum_prelam_history_20260916.json`), 100 trials/stage via
  `_curriculum.py` chain (`run_curriculum_laminated.bat`, log
  `curriculum_lam_20260916.log`). Stage-1 best reproduced the 25-trial
  laminated run exactly (57.436 @ trial 22, seeded determinism) and
  extended it. Final deliverables chain: `post_curriculum_deliverables.bat`
  (_final_run_lam.py applies the stage-3 winner via _curriculum.set_stage
  — bit-faithful, unlike old _final_run.py which skipped v10 multiplier
  overrides).

## 2026-09-16 later (Ben: "there are incorrect connections at the RG and
## PF layers" — PFA removed, commissurals laminated, ADAP proven required)

**RG/PF layer corrections (build_network.py):**
- **PFA self-adaptation loops REMOVED** from the PF layer (Ben flagged
  09-15 "don't exist in Shevtsova"; kept then; flagged again 09-16 —
  now gone). PF cells keep their PF_SHAPE tau_m; the tau_a half of
  PF_SHAPE is INERT (kept in params for the record). Consequence: E2's
  burst is no longer hard-truncated (its tail decays with the 72 ms
  membrane tau instead) — stage-1 winner config re-scored 57.4 -> 43.9
  (rhythm intact, 2 rises, knee -75.8); full re-tune required.
- **Cross-side commissurals LAMINATED**: the direct RG-F_r->RG-F_l and
  RG-E_r->RG-E_l inhibitions are now routed through per-side commissural
  INs CIN-F_side / CIN-E_side (tau = TAU.rg; same effective gains
  rg_mutual_inh and 0.5x; Shevtsova V0 analog — no direct synapse
  between pattern-generating pools across the midline either).
  Network counts: 203/side + 8 shared = 414 neurons (was 418).
- **ADAP loops KEPT — proven load-bearing** (`_repro_s1.py`): with
  rg_adapt_inh=0 the rhythm dies completely (0 bursts, objective
  57.4 -> -4.2, joints freeze tonic). The plain NonSpikingNeuron
  half-center has no intrinsic burst termination; ADAP is the
  toolbox-native substitute for Deng's persistent-Na h-gate (which
  sns_toolbox 1.5.2 neurons cannot express; consistent with
  _tau_h_check.py: toolbox tau_h(V) quenches Deng oscillators). If Ben
  wants literal Deng neurons, that is the Simulink SNS_Deng_Library
  route — his call.
- End-to-end proof recorded en route: the heel/IaIN fixes left stage-1
  dynamics BIT-EXACT (trial-22 re-evaluated 57.436 on the fixed code,
  pre-PFA-removal).

**Per-muscle reflex audit vs Ben's spec** ("each flexor should have Ia
and II feedback; each extensor type II and mechanosensory feedback") —
ALREADY SATISFIED, no change needed: every muscle has its own MN + Ia +
II + Ib neurons (x92) in build_network; the runner's presynaptic gates
differentiate: flexors get FULL Ia gain (swing stretch-velocity reflex)
+ 0.5x II; extensors get stance-BOOSTED II (0.3+0.7*stance) and their
autogenic Ib is stance-SUPPRESSED (1-0.7*stance) so extensor force
routes through the IBEXC stance load-sharing group + LBIN instead
(mechanosensory), plus heel/toe contact mechanosensors at the RG.

**Synergy -> PF-layer study** (`_synergy_pf_study.py`, NMF on the
subject01 SO backsolve): explained var 0.72@3 / 0.89@4 / 0.92@5 /
0.95@6 — no sharp elbow, BIC keeps improving to 8. Caveat: the SO file
covers only ~1.6 gait cycles (126 frames, 2 s), so phase-folds are
coarse and the leading variance direction is LEFT-vs-RIGHT split (one
stance synergy per leg) — classic with short records. The windows at
4-5 synergies: push-off (plantarflexor+hip-ext burst, wrapping the
cycle boundary), contralateral stance, hip-abductor tone spanning
stance (10-55%), late-stance/knee-ext prep (60-85%), early swing
(5-15%). ANSWER: the data supports 4-5 phase windows per side — the
CURRENT 4-cell PF (E1/E2/F1/F2) has the RIGHT COUNT as a phase
decomposition, but the data reshuffles composition: hip ABDUCTORS
belong to EARLY/mid-stance (weight acceptance) while the fitted table
puts them mostly in F2 (0.301!) — that column is contradicted by the
data and is a retune target. A full synergy-based PF (task-defined
cells crossing joints) remains the bigger redesign option.

**Biarticular hip+knee muscles — which PF layer?** Currently BOTH:
W_PF_MN weight = primary group FULL + secondary group 0.5x, per phase.
From the fitted table: rect_fem (knee_ext+hip_flex) = 0.5x of each;
hamstrings bifemlh/semiten/semimem (knee_flex+hip_ext) E1 = 0.129 +
0.5x0.083 = 0.170; gastrocs (ankle_pf+knee_flex) = 0.068+0.5x0.129.
Physiologically defensible per the synergy data (hamstrings appear in
both the push-off AND swing synergies) — keep, it is a tunable table.

**Curriculum on the corrected architecture:** curr_s1_air_deaff
PURGED (its 86 trials were tuned with PFA present — unconditional
topology, so PFA removal contaminates ALL stages, unlike the heel/IaIN
fixes). Full re-run relaunched 2026-09-16 evening:
stage1->stage2->stage3 x 100 trials (`resume_curriculum_lam.bat`, log
`curriculum_lam_20260916.log`). Diagrams (panels + deng) regenerated to
match; stale old figures/gifs already deleted (earlier 2026-09-16
section).


## MORNING REPORT — overnight 2026-09-12/13 (v5 phase-reset + Shevtsova
## figure + v6 quad suppression) — READ THIS FIRST

**TL;DR: the swing-knee quad suppression (lever #2) is the lever that
works: best kine_score improved −61.18 → −59.43 (eval) and −62.7 →
−60.6 (full 22 s), all runs stayed up. The tonic phase-reset (lever #1
as formulated) is effectively inert — the optimizer avoids it. Also:
the regression gate caught a REAL pre-existing bug — `--best` never
applied pf_gain (silent since v3).**

1. **BEST CONFIG = the v6 winner** (study ground_walk_v6_kneext, trial
   106): `python runner.py --fitted --best6` (full 22 s) or add
   `--eval` (12 s, scores −59.426). Reproduction verified:
   −59.4xx bit-consistent. v5 winner (no suppression): `--fitted
   --best5`. Old v4b: `--fitted --best`.
   Key values: f1_kneext_inh 0.590, phase_reset_e 0.388,
   phase_reset_f 0.035, drive 2.85, pf_gain 2.13, rg_adapt 1.10,
   desc_e 1.75, e2_adapt 1.90, kx 121.
2. **Numbers (full 22 s, winner)**: kine −60.591, duty 0.473 (ref 0.61),
   knee_min +10.1 deg (ref −69.7 — still extension-dominant mean
   cycle), cadence 1.4 Hz (ref 0.81), tilt 35.3 deg, stayed up.
   Eval-window: kine −59.426, n_cycles 8, rmse_knee ~19.5.
3. **WHAT MOVED THE NEEDLE**: v6's F1→KINH→knee_ext-MN inhibition
   (phase-gated by F1 = swing-only). TPE engaged it immediately (top
   trials f1_kneext_inh 0.4-0.6), buying knee-shape RMSE (20.9 → ~19.5)
   and a large full-run jump. v5's HIP_EXT_SIG/HIP_FLEX_SIG →
   PRESET → RG tonic gated inputs: sign-audited PASS but the optimizer
   drives those gains toward 0-0.25 — a ≤1 nA gated tonic input cannot
   move an RG running on ~4 nA of DRIVE.
4. **WHAT FAILED / SURPRISES**:
   - **`--best` pf_gain bug (PRE-EXISTING, FIXED)**: `if "pf_gain" in
     best["params"]` was always False (pf_gain lives at the JSON top
     level) → the gain-scaling branch NEVER ran → every `--fitted
     --best` "reproduction" since v3 silently evaluated a pf_gain=1.0
     config. The recorded v4b/v3 full-run numbers (E-duty 0.27, knee
     −14..+26, −63.3) describe gain-1.0 configs, not the true winners.
     Fixed (scale only fitted-file entries; pf_gain from doc root);
     gate now BIT-EXACT: `--fitted --best --eval --drive 2.5868` =
     −61.1741850610298 on the v5 code at gains 0.
   - phase_reset_e ≈ 0.6-0.75 + weak drive LATCHES RG-E on (no rhythm;
     the stance-gated ext-signal is positive feedback) — stability
     constraint for any future sensory→RG pathway.
   - drive sensitivity is small (4-dp rounding ≈ 0.005 kine) — v5/v6
     evaluate with repr(drive) so their jsons reproduce exactly.
5. **DELIVERABLES**: `draw_circuit.py --vclasses` →
   `circuit_full_vclasses.{pdf,svg,png}` in Dissertation\
   CPG_airstepping_figs (Shevtsova-2026/Rybak-2015 V-class strip + tags;
   class-by-class mapping table in the section below); circuit_*
   figures REGENERATED (their pf_gain composite was fixed too — weights
   figure numbers changed slightly); v5/v6 sweeps + trials CSVs
   (v5_sweep.csv, v5_results.csv, v6_results.csv); per-global-best full
   runs v5_best_trial*.npz (9) and v6_best_trial*.npz (10); audits
   _phase_reset_audit.py PASS, _kinh_audit.py PASS; RUNNER_DUMP_STATE
   env hook in runner.py for config diffing.
6. **NEXT SESSION (priority order)**:
   a. The RG duty/cadence structure is THE remaining bottleneck (duty
      0.47 vs 0.61, cadence 1.4 vs 0.81 — suppression sped the cycle
      up). Try: transient phase-ONSET reset (flexion-velocity PULSE
      into PRESET_F, not tonic), stance-biased PF windows, or
      weaning-relevant pelvis-balance work per the v4 list.
   b. knee_min still +10: consider a stronger KINH range (searched 0-2,
      winner 0.59 — not boundary-pinned; maybe search 0-4 with E-duty
      shaped reward) or ankle-PF override retune (v4 item 3).
   c. plot_run.py on v6_best_trial106.npz + refresh ground figs/gifs in
      the Dissertation folder (NOT done tonight).
   d. ~~The trunk Fmax fix (8 one-newton muscles) still needs Ben's
      go~~ — DONE same day, see the postscript below.
   Nothing committed (per standing rule); sqlite has both studies for
   resume (`python optuna_walk_v5.py 50` continues v5).

## 2026-09-13 (day): trunk Fmax fix (Ben's go received), torque budget,
## v6 render + dissertation draft

- **Fmax audit** (`_fmax_audit.py`): every actuator's gainprm[2] vs
  stock `gait2392_simbody.osim` max_isometric_force. The 84 leg/active
  muscles all match within the converter's RoM normalization
  (≤15%, both directions — Ben: "use OpenSim muscle max force over
  human RoM as the benchmark", which is what the converter does).
  **8 trunk actuators shipped at 1 N active** vs stock ercspn 2500,
  intobl/extobl 900, ext_hal 162 (the passive slot already carried the
  stock value; lengthrange "0.01 1" is also converter garbage —
  flagged, not fixed). `_fmax_audit.py` re-runs the check any time.
- **Fix** (patch_xml repair 2c, build-asserted to exactly 8 tags):
  gainprm[2] := stock values. **Post-fix the v6 winner IMPROVES**:
  eval kine −59.241 (was −59.426), tilt_max 27.4 (was 33.3), kz 0.909,
  dx halved — the BAL_TRK trunk controller finally has real muscles
  (measured ercspn act=1 force 1637 N = 65% of Fmax at the standing
  pose despite the bad lengthrange). Standing solve now engages 16
  muscles >0.1 (was 8). No retune required; old `--best`/`--best5`
  numbers were pre-fix and remain valid AS RECORDED.
- **Torque budget** (`_torque_budget.py` + `_torque_demand_walk.py`):
  mass 78.5 kg (770 N); fd moment arms at the standing pose × stock
  Fmax: hip_ext 251 N·m/side, knee_ext 277, ankle_pf 296, trunk_ext
  213. Demand from the bsolve ID residuals (subject01 + measured GRF):
  hip 74, knee 98, ankle 212, lumbar 100 N·m. MARGINS: hip 3.4x,
  knee 2.8x, trunk 2.1x, ankle 1.4x (tightest; optimistic end — the
  PF arm shrinks in late-stance dorsiflexion; the 2.7 N·m/kg ankle
  demand also runs hot vs ~1.5 literature, likely the known subtalar/
  CoP residuals). Pre-fix trunk supply was 0.09 N·m — the rig was
  doing 100% of the trunk work.
- **v6 render**: ground_walk_v6_trial106.gif + fig1..5 PNGs rendered
  from v6_best_trial106.npz (spinal_run.npz now holds that run);
  copies in Dissertation\CPG_airstepping_figs.
- **Dissertation draft**:
  `Documentation\Reports and Papers\Dissertation\CPG_spinal_section_draft.tex`
  — architecture, V-class correspondence (table), v5/v6 refinement,
  plant audit + torque budget, honest limitations (rigid tendon,
  no V3). Ready to adapt; figures referenced from
  CPG_airstepping_figs/.

## 2026-09-13 (Ben's critique): Deng-style circuit figure, drawn FROM
## the compiled network (circuit_dengstyle)

Ben rejected the hand-laid-out circuit_full ("pile of garbage") and set
the standard: the Deng/Nourse Fig-6A schematic — layered RG/PF/motor
circuit with EVERY IN and connection explicit, red conversion-map
circles at the muscle interface, Renshaw + IaIN + IbIN visible.

**The reference, exactly** (Nourse et al. 2023 Fig 6A + Table A6 — the
image Ben sent; paper PDF pulled from MDPI, full synapse table
extracted; also _nourse2023.txt kept here):

RG layer: 2 half-center neurons (voltage-gated fast-Na, endogenous
bursting; Table A5) that **mutually inhibit via TWO dedicated
non-spiking INs** — HC_ext→(exc 2.749 µS) Ext-IN→(inh 2.749) HC_flx and
mirror. NO direct HC↔HC synapses and NO mutually excitatory RG synapses
anywhere in this circuit (if you remember "mutually excitatory", that
is not this figure — the mutual coupling is entirely inhibitory, and
speed is set by its level). PF layer: one HC pair PER JOINT (hip; knee+
ankle share one), "constructed in a similar manner" (same IN-laminated
mutual inhibition); RG HC→ same-half PF HC, weak exc 0.1 µS. Motor
circuit per joint: PF→MN exc (1.5–4.9 µS per muscle); Ia afferent
(formatted from muscle TENSION) → IaIN (PF-gated 0.5 exc — phase
control of reciprocal inhibition; IaIN→antagonist MN 2.0 inh;
IaIN↔IaIN 0.5 inh); MN→RC 0.5 exc, RC→MN 0.5 inh, RC↔RC 0.5 inh,
RC→IaIN 0.5 inh (disinhibition); Ib afferent→IbIN→MN 0.59 µS
**EXCITATORY** (positive force feedback), PF→Ib 2.0 µS near-rest shunt
(graded presynaptic control). Interface: act = 1/(1+e^{s(x0−V)})+y0,
s=0.153, x0=−70 mV, MN Vrest −100 mV (Fig 6B); feedback = muscle
tension formatted as Ia and Ib input currents. No II afferents, no
interleg, no balance cells in Deng.

**Ours, exactly** (_net_edges.py dumps every connection from the
COMPILED network object — 118 edges in the 4-muscle representative
build; the full net scales to 92 muscles):

- RG: RG-E↔RG-F mutual inhibition DIRECT g=4.0 (Deng's IN lamination
  LUMPED — drawn as ghost Ext-IN/Flx-IN in the figure); ADAP-E/F
  adaptation INs (RG→ADAP exc 2.5, ADAP→RG inh 2.5) — burst
  termination setting cycle period, the non-spiking substitute for
  Deng's persistent-Na; DRIVE→E 1.7/F 1.4, POSTURE→E 0.8; commissural
  RG-F_r↔RG-F_l inh 4.0, RG-E_r↔RG-E_l inh 2.0 (Deng has none);
  v5 PRESET_E/F INs (hip ext/flex signals; E↑F↓ / F↑E↓).
- PF: 4 phase cells E1/E2/F1/F2 + PFA self-adaptation INs (Deng: HC
  pairs per joint); RG→PF exc 2.4; DRIVE→PF 0.05; PF↔PF reciprocal
  inh 4.0 (4 pairs); v6 KINH (F1→KINH 1.5, KINH→knee-ext MNs 0.59,
  swing-gated).
- Motor circuit per muscle (×92): Ia→MN homonymous EXC 0.6 (Deng has
  NO monosynaptic Ia→MN — ours is closer to biology); Ia→antagonist MN
  INH 0.4 DIRECT (Deng: via IaIN with PF phase-gating + RC
  disinhibition — ghost in the figure); II→MN exc 0.4 (not in Deng);
  Ib→MN autogenic INH 0.35 (Deng: IbIN→MN EXCITATORY — opposite sign;
  our positive-force pathway is the stance-gated IB-EXC reversal:
  Ib→IBEXC 0.5, RG-E gate 1.0, IBEXC→extensor MNs 0.6);
  PF→MN × W_PF_MN; POSTURE/POST_i; BAL family (full net only).
- **Renshaw: NOT IMPLEMENTED** (ghost in the figure, Deng's exact
  wiring quoted). One-line add in _wire_muscle if wanted
  (MN→RC→same-pool MN, conditional gain like KINH).
- Interface maps (runner.py, now drawn as red circles): a =
  clip(V/5 mV, 0, 1) linear-saturating (Deng: sigmoid); feedback
  len/vel/force-normalized encoders (L−Lmid)/Lhalf, L̇/0.6 m/s, F/Fmax
  with presynaptic phase/speed gating (Deng: tension-based Ia/Ib only).

**Figure**: `draw_circuit.py --which deng` → circuit_dengstyle.{pdf,
svg,png} (figures/ + Dissertation\CPG_airstepping_figs). STRUCTURE-
DRIVEN: it builds the representative network, reads every connection
from the compiled SNS object, and ASSERTS each (src,dst,sign) group is
drawn — the figure cannot drift from the code again. Deng's Table A6
is quoted in the figure footer for side-by-side comparison. The old
circuit_full/_vclasses stay (both-sides + commissural view + V-class
strip) but the dengstyle figure is the one that shows the full
connection grammar.

## 2026-09-14: RoM limits + Renshaw (Ben's go), toolbox-native diagram,
## tutorials + skill

- **JOINT RoM LIMITS ON (patch_xml repair 2f2, unconditional)**: all 14
  driver hinges shipped limited="false" (ranges inert). Flipped to
  limited="true" with the CONVERTED/stock ranges (knee [-120,+10] deg
  is the true stock RoM; hips +-120; ankle/subtalar/mtp +-90 - stock
  gait2392 itself uses wide ranges, verified by
  _parse_osim_ranges.py). Followers stay unlimited (repair 2f).
  **Result: knee hyperextension GONE** (ground walk knee now
  -0.3..+11.4 deg vs -14..+26 before; the +10 cap visibly enforced).
- **RENSHAW CELLS IMPLEMENTED** (Ben's go): per-pool RC with MN->RC exc
  1.0, RC->homonymous MN inh G["renshaw"], RC<->RC inh (each pair
  once); conditional topology, default 0, `--renshaw X` / json key
  "renshaw". Deng A6: 0.5/0.5/0.5 - we ran 0.5.
- **Ground walk (limits + Renshaw 0.5, `runner --fitted --best6
  --renshaw 0.5`)**: stayed up, COM min 0.90 (best yet), tilt 30.5,
  hip -0.9..+18.3, knee -0.3..+11.4, ankle -40..0, 15 cyc at 1.32 Hz,
  E-duty 0.15. HONEST: gait re-timed around the limits (the old
  hyperextended knee was load-bearing) - kinematics retune (v7) needed
  for duty/cadence; the ARTIFACT is fixed. Run saved v6_rom_ground.npz.
- **AIR-STEPPING TRANSFORMED** (`--no-ground --no-afferents
  --no-interleg`, v6_rom_air.npz): knee -97..+10 deg TRUE deep swing
  flexion (cap respected), hip -21..+50, 0.41 Hz (5 cycles in window)
  - right at the Ivanenko air-stepping regime (~0.3 Hz preferred; we
  were 0.58 Hz with a pinned knee). Pelvis tilt +-2.5 deg in air.
  THIS is the dissertation air-stepping result.
- Figures: `hindlimb_style_{ground,air}.png/.pdf` (_figure_hindlimb_
  style.py; SNS-Toolbox figure_hindlimb layout: RG / PF / MN rows over
  knee + hip angle rows) + ground_rom.gif / air_rom.gif (all in
  Dissertation\CPG_airstepping_figs).
- **Toolbox-native diagram** (`_render_diagram.py` + `spinal_layers.py`):
  the circuit rebuilt as Tutorial-4-style Network SUBCLASSES
  (RhythmGeneratorNetwork / PatternFormationNetwork / MotorColumnNetwork,
  conductances live from params) composed with add_network and rendered
  by the OFFICIAL sns_toolbox.renderer -> sns_diagram_spinal.png
  (Dissertation folder too). This is the sns_diagram.png style: colored
  layer groups, invhouse inputs, house outputs. 1.5.2 has NO RG/PF
  presets (only the 5 arithmetic ones) - our layer classes are the
  reusable pattern. NOTE: add_network flattens; per-layer colors carry
  through.
- **test_max_realistic_net_size.py** (Ben Q): upstream benchmark,
  master API, default Torch - does NOT run unmodified on 1.5.2
  (backends take the compiled params dict via net.compile there).
  Bounded numpy version `_test_max_net_numpy.py`: 8000-neuron spiking
  net compiles + steps in 4.4 s (~N^2; memory wall ~25-30k on 32 GB).
  Our 410-neuron net is nowhere near limits.
- **TUTORIALS DONE + SKILL WRITTEN**: all 9 official tutorial notebooks
  fetched to spinal\sns_tutorials\; 1/2/4/8 executed clean after
  installing the graphviz BINARY (conda install -n myo graphviz; pip
  graphviz is bindings-only; dot.exe lives in <env>\Library\bin, not on
  PATH unless activated). Skill at
  C:\Users\Ben Bolen\.zcode\skills\sns-toolbox\SKILL.md (junction into
  ZCode_Skills repo - Ben commits): env, API traps (shape arg, compiled
  params, two class trees, connections-dict schema, flattening
  add_network), Tutorial-4 pattern, rendering, net sizing, our MuJoCo
  interface conventions.
- Trunk lean (Ben: "fix the femurs, trunk supporting itself"): ground
  tilt 30.5 deg post-limits (was 33-35). The lean is hip-extensor
  saturation pitching the pelvis backward (v4 diagnosis); with real
  trunk muscles (Fmax fix) BAL_TRK now has authority but its gains are
  pre-fix. Next: retune kp_trk/kd_trk + hip-ext posture bias with
  limits+Renshaw in place, THEN wean the pelvis-pitch rig spring.

## 2026-09-14 (Simulink session cross-report): Deng CPG verified, tau_h
## collapse CONFIRMED, constant-DRIVE validated

The MuJoCo-Simulink bridge session (see README_SNS_Simscape.md + the
AGENTS.md bridge section for their full write-up) built Deng's RG and
PF layers as Simulink .slx (persistent-Na HCs, IN-laminated mutual
inhibition, params from Nourse Tables A4-A7) and verified them against
an independent numpy ODE reference (`spinal\deng_cpg_ode.py`): ONE
10 nA / 20 ms pulse -> self-sustained 1.938 s alternation (numpy:
1.94 s), RG ext/flx correlation -0.86. It runs indefinitely on constant
DRIVE.

- **tau_h(V) collapse CONFIRMED numerically** (_tau_h_check.py, this
  session): with Deng's h-gate (S=-0.6, E=-60, tau_max=350 ms),
  sns_toolbox's tau formula tau(V) = tau_max * z_inf(V) *
  sqrt(K*exp(S*(E-V))) peaks at only ~140 ms near -50 mV (K=0.01) and
  collapses to <0.5 ms at V >= -30 mV, -> 0 at 0 mV, for every K. The
  slow burst-terminating process is destroyed exactly where the
  neuron spends its burst -> the oscillator quenches. Fixed
  tau_h = 350 ms restores clean bursting (their verified fix).
  IMPLICATION FOR OUR NETWORK: this is WHY the ADAP adaptation loop is
  mandatory in our non-spiking RG - the toolbox has no working slow
  intrinsic process in either its non-spiking neurons or its
  persistent-Na implementation. Any future toolbox port of
  Deng/persistent-Na cells needs the fixed-tau patch (or a custom
  neuron). Also from that session: SNS_Library vs toolbox use OPPOSITE
  synapse-saturation conventions (ThrPre/Elo) - keep straight when
  porting values.
- **Constant DRIVE validated**: the tuned 410-neuron net self-sustains
  at constant DRIVE = 2.5 nA for 20 s (1.21 s alternation, constant
  amplitude). The Simulink copy's tonic settling was an integrator
  summation-order artifact (same chaos family as our fp-drift findings)
  - the network sits near a bifurcation. Constant descending drive is
  the correct operating mode; and for the future basal-ganglia layer,
  descending control should be modeled as a bias RELEASED
  (disinhibition), not commanded - Ben + session agreed.
- **Animatlab latch cross-pointer**: the ~0.1 ms tau_h collapse is the
  prime suspect for the Animatlab RG latch (Animatlab treats
  tau_h.max as a fixed constant, which is why the Simulink port
  works). First thing to check in the laptop .aproj: how the
  LinearHill Na h-gate time constant is implemented.

## 2026-09-14: v7 retune under corrected physics (RoM + Renshaw + ANH)

- **Ankle dorsiflexion (Ben's complaint) — solved mechanism**: boosting
  swing DF drive x3-x5 does NOTHING (max stays -14 deg: 9 PF muscles
  ~10 kN vs 3 DF ~1.6 kN + ligament spring). What works is SWING-PHASE
  PF SUPPRESSION: `params.G["f1_anklepf_inh"]` routes the F1-driven
  KINH IN onto the ankle_pf MNs (same pattern as the knee). At 1.0 the
  air-stepping ankle sweeps **-55..+7.1 deg — true dorsiflexion**
  (range 62 vs ref 23; overshoot = tunable). Tests:
  _ankle_test.py / _ankle_anh_test.py.
- **v7 study** (`ground_walk_v7_rom`, optuna_walk_v7.py, 140 trials,
  seeded v6 winner, RENSHAW 0.5 FIXED, f1_anklepf_inh searched): best
  **-62.523 (trial 100)**, plateaued from ~trial 100 (100-137 flat).
  Reproduces BIT-EXACT via `runner --fitted --best7` after a json fix
  (below). Winner full-22s: stayed up, tilt 28.6 (best yet), duty
  0.18-0.20, knee pinned AT the +10 cap (mean-cycle knee range 0.3
  deg!), cadence ~0.4 Hz. HONEST: under corrected physics the walk is
  a slow stiff-kneed shuffle — the unlimited-model crutch (-59.4) cost
  ~3 points and 140 trials recovered half. Scalar tuning has CONVERGED
  again; duty/cadence/knee-flexion are architecture-level (transient
  phase reset, FSA-analytic seeding of PF/RG timing).
- **Interesting**: under the capped knee, phase_reset_f IS engaged
  (0.5-0.8 in the top trials; it was avoided in v5) — swing-trigger
  feedback matters once flexion is physically cheap.
- **NEW JSON LESSON (same family as pf_gain)**: the generated v7 set
  renshaw=0.5 in set_params but never wrote it into the json params ->
  `--best7` silently ran Renshaw-less (kine off by 0.14, caught by the
  reproduction check). Fixed: script writes renshaw into eff dict +
  json patched. RULE going forward: EVERY new params.G knob a study
  uses must appear in the saved json AND have an `if key in best`
  branch in the --best loader — the repro check is what proves it.

## 2026-09-14 (Ben's start-pose directive): THE POSE WAS THE UNLOCK

Ben: start the model leaning FORWARD ~5 deg, knees 5-15 deg flexed,
ankles dorsiflexed, legs in different stance phases (L hip flexed, R
extended) - the old keyframe started dead-straight, trunk already back.

- `apply_start_pose()` in runner.py (default on; `--straight-start`
  reverts): pelvis_tilt -5 (sign auto-checked against the torso
  up-vector: reports "torso lean -5.0 deg"), knees -10, ankles +5,
  hip_flexion_l +25 / hip_flexion_r -15, pelvis_ty 0.90, followers
  re-projected (local copy of bsolve's apply_eq_followers; runner
  imports bsolve circularly). Applied BEFORE capture_pose (now reads
  data.qpos, not model.key_qpos), the standing solve, and the rig (so
  springref holds the new pose). Standing solve now engages 21 muscles.
- **v8 study** (`ground_walk_v8_pose`, seeded v7, 60 trials): best
  **-47.889 (trial 55)** - a 14.6-point jump over the v7 plateau
  (-62.5) and far past every earlier result. The tuned winner leans on
  the new machinery: f1_kneext_inh 1.00, phase_reset_f 1.72, drive 3.0.
  First full-22s capture at trial 24 already -49.6.
- Winner full-22s: stayed up, knee_min -92.9 deg (deep flexion IN RoM),
  tilt 16.9 deg (best trunk yet), duty 0.24, cadence 0.3 Hz, full-run
  kine -62.5. CAVEAT: eval (12 s window) vs full (22 s incl. ramp-down)
  diverge - the slow 0.3 Hz cadence means few cycles late; extending
  the study (+80 trials, running) to shape cadence/duty.
- Panel diagram: `_render_panels.py` renders each spinal_layers
  subnetwork separately via the official renderer and composes ->
  sns_diagram_panels.png (RG / PF / motor panels, tuned gains labeled).
  THE tidy layered look Ben asked for, still 100% code-driven.

## 2026-09-14 (Ben): start pose replaced by the model's own normal.mot

Ben supplied the gait2392 Coordinates-panel values for the "normal"
pose (normal.mot; degrees) and asked why we don't use THAT as the start
pose. Correct - it replaces my hand-invented numbers. Now in
START_POSE_DEG (runner.py): pelvis_tilt -1.87, list -0.47, rotation
+2.43, hips R +24.6/-1.35 add/rot +1.09, knee_r -3.94, ankle_r -1.7,
hips L -16.6/+2.68/+1.28, knee_l -8.2, ankle_l +9.8 (DORSIFLEXED),
lumbar +1.87/+0.47/-2.43, pelvis_ty 0.96. Sign mapping is direct: our
repair-1 hip axis flips exist precisely so OpenSim +flexion/+adduction
== MuJoCo +qpos. **AUTO-SIGN-FLIP REMOVED**: my IMU-lean heuristic
flipped -1.87 -> +1.87 (measured IMU lean at the OpenSim values is
-0.0 deg = upright, which is what "normal" means) - the pose is now
applied AS GIVEN, IMU lean only reported. 20 coordinates applied.
- v8 (hand-pose) continuation STOPPED mid-run; fresh study
  **ground_walk_v8b_normal** (optuna_walk_v8b.py, seeded v7 winner)
  running on the canonical pose; writes best_walk_params_v8b.json
  (`--best8b`), v8b_results.csv, v8b_best_trial*.npz.
- The v8 hand-pose breakthrough stands as the finding (pose quality is
  a first-order lever: -62.5 plateau -> -47.9 in 60 trials); v8b gives
  the canonical-pose numbers. Hand pose vs normal.mot are qualitatively
  close (forward lean, knees slightly flexed, split hips) so the v8
  tuned direction carries as the v8b seed.
- **v8b first launch COLLAPSED onto the frozen plateau (-65.000 best,
  every trial) — SECOND sentinel-calibration lesson** (the first was
  v4's -25): the normal.mot pose is nearly upright, its genuine walkers
  score **-76..-86** on the 12-s eval, so the -65 no-rhythm sentinel
  and -80 NaN gate sat ABOVE real walking -> frozen models outscored
  walkers and TPE had nothing to learn. Fixes: (1) eval schedule
  lengthened (1 s stand / 1 s ramp / 11 s walk / wind-up = 16 s total)
  so 0.3-0.9 Hz gaits yield >= 3 countable cycles; (2) gates lowered to
  -100 (frozen) and -110 + t_end (NaN); (3) seed = the v8 hand-pose
  winner (a demonstrated stepper). Fresh study running
  (ground_walk_v8b_normal; the poisoned one deleted from the db).
  RULE: whenever the plant, pose, or objective changes, re-check that
  the no-rhythm/NaN sentinels sit BELOW the worst genuine walker.

## 2026-09-14 (Ben): TIPTOE ROOT CAUSE + amplitude-honest objective

Ben (looking at opensim_overlay_gait_cycles.png): the sim curves are
flat lines (hip ~25, knee ~+10 pinned, ankle ~-33) vs OpenSim's full
sinusoids — "you're fitting a flat line to a cosine curve and calling
it better", and the model walks on TIPTOES; check the talus/calcaneus/
MTP angles in OpenSim.

- **TIPTOE ROOT CAUSE (measured, _foot_flat_check.py)**: the FOOT
  angles I took from normal.mot were correct (ankle -1.7, subtalar 0,
  mtp 0 = flat), but I also took pelvis_ty = 0.96 from that table.
  With MESH-VERTEX foot bottoms measured in OUR converted model, the
  soles sit 6.3 cm (right) / 3.0 cm (left) ABOVE ground at ty 0.96 —
  the rig held the pelvis there, the feet dangled into passive PF
  (ankle -33 constant = the tiptoe). OpenSim's 0.96 does not match our
  converted foot geometry. Right sole touches at ty 0.903, left (back
  leg) at ~0.925 -> START_PELVIS_HEIGHT = 0.92 (contact compromise the
  knees + soft contact absorb). Foot angles STAY from normal.mot.
- **OBJECTIVE FIX (kine_ref.py)**: amplitude terms strengthened — range
  errors now cost 0.5 (hip) / 0.3 (knee, NEW - knee range previously
  only entered via knee_min so a knee pinned at the RoM cap scored the
  same as real flexion) / 0.5 (ankle), up from 0.10/none/0.10. A
  flat-line gait now totals ~-121 while real steppers land
  -45..-75. The old v8b winner re-scored -74.4 under the honest
  objective+height (and its knee range went 0.3 -> 41 deg at the
  flat-foot height - the tiptoe support was loading the walk).
- Both changes invalidate prior tuning -> fresh study
  **ground_walk_v9_flat** (optuna_walk_v9.py, seeded v8b winner,
  `--best9`, v9_results.csv, 60 trials running). Sentinel gates
  (-100/-110) re-checked: still below the -95-band of genuine walkers.
- Air-stepping "stationary bicycle" (Ben): hip/knee are already
  roughly there (hindlimb_style_air: hip -21..+50, knee -97..+10
  quasi-sinusoidal); the ankle PF hang was the passive-bias issue -
  f1_anklepf_inh now gives dorsiflexion, and amplitude tuning should
  finish the bicycle look.
- **v9 RESULT (60 trials)**: best **-65.402 (trial 36)**, reproduced
  bit-exact (`runner --fitted --best9`). AMPLITUDE SOLVED — the
  winner's mean cycle: range_hip 46.4 (ref 43.3), range_knee 66.9
  (ref 70.5), range_ankle 24.1 (ref 23.1), knee_min -78.1 (ref
  -69.7), tilt 15.7. What REMAINS is phase/timing, not amplitude:
  rmse_hip 30.6 with full excursion (hip phase roughly inverted vs
  the RG-anchored cycle), stance starts crouched (knee -67 at
  contact), ankle oscillates around a -45 deg PF OFFSET (set-point,
  from POSTURE_OVERRIDE soleus/tib_post + balance PD), duty 0.18,
  cadence 0.36 Hz. Overlay regenerated from v9_best_trial36.npz
  (Dissertation folder). NEXT: the remaining gaps are exactly the
  transient phase-reset / FSA-analytic-seeding items + an ankle
  set-point trim; extension +60 trials running (optuna_v9b.log).
- **v9 extension (+60, 120 total)**: best **-62.243 (trial 119)**,
  still climbing but decelerating (~3 pts / 60 trials). Winner full-22s
  (v9_best_trial119.npz, also copied to spinal_run.npz + Dissertation
  folder): stayed up, tilt 16.6, knee -81.4 (RoM-respecting), duty
  0.167, cadence 0.3 Hz. 4-cycle overlay regenerated: knee/hip
  waveforms real, ankle rides a -50 deg PF OFFSET with small
  oscillation. The study is resumable (`python optuna_walk_v9.py 60`);
  stopped extending here - the remaining deficits (duty, ankle
  set-point, hip phase vs RG anchor) are the architecture items, not
## 2026-09-15 (Ben's architecture critique): LAMINATED IN-mediated mutual
## inhibition — the REAL Shevtsova/Deng architecture

Ben's critique (with Shevtsova Fig 2 and his own Animatlab diagram as
evidence): our RG had DIRECT RG-E↔RG-F inhibition (no INs), and our PF
had DIRECT PF↔PF reciprocal inhibition with PFA SELF-adaptation loops.
Shevtsova/Deng use IN-LAMINATED architecture: RG-E excites InE, InE
inhibits RG-F (never direct); PF-E excites PF_IN_E, PF_IN_E inhibits
PF-F (never direct). PFA self-loops don't exist in Shevtsova.

**FIXED in build_network.py:**
- _build_rg: added InE_<side> and InF_<side> INs (tau = TAU["rg"]).
  RG-E → InE exc (g=rg_mutual_inh) → RG-F inh (g=rg_mutual_inh);
  RG-F → InF exc → RG-E inh. Direct RG↔RG edges REMOVED.
- _build_pf: added PF_IN_E_<side> and PF_IN_F_<side> INs (tau = TAU["pf"]).
  E1+E2 → PF_IN_E exc → F1+F2 inh; F1+F2 → PF_IN_F exc → E1+E2 inh.
  Direct PF↔PF inhibition edges REMOVED.
- PFA kept (documented: our window-shaping addition, not in Shevtsova).
- ADAP kept (burst termination, separate from mutual-inhibition routing).
- Mutual Renshaw fix applied (both directions between distinct RCs).
- IaIN population added (Ia → IaIN → antagonist MN, PF_F1 phase-gated,
  RC→IaIN inh) — replaces direct Ia→antagonist when G["ia_in"] > 0.
- AFF_E/AFF_F semi-closed sensory loop relays added (Shevtsova principle:
  active phase's afferents excite that phase's PF and RG through relay
  INs — three-layer loop: muscle → afferent → PF → RG → MN → muscle).
  Gains: aff_e_rg, aff_f_rg, aff_e_pf, aff_f_pf (all default 0).

**Architecture is now CORRECT per Shevtsova/Deng.** All prior tuning
invalidated (different wiring = different dynamics). Curriculum stages
1-3 must be re-run. Stage 1 laminated: best 57.4 (trial 22) — lower
score expected because the laminated path has two synapses per
inhibition vs one direct, halving effective inhibition at same g; the
optimizer compensates with higher rg_adapt and desc_f.

**NMF synergy analysis** (_synergy_nmf.py on backsolved activations):
3 synergies = 82% var (left stance, right stance, transition — classic
Ivanenko). 5 synergies = 92% (separates hip abductors, knee extensors,
plantarflexors). KEY FINDING: synergies CROSS joint boundaries — no
clean "hip PF" vs "knee PF" decomposition. The natural grouping is
task-based (stance push-off, swing initiation), not joint-based. The
current 4-cell phase-based PF is a simplification; a synergy-based PF
would be data-driven but requires full redesign.

## 2026-09-15 (Ben's curriculum directive): staged tuning + mechanosensory
## stance feedback + IaIN — the honest state

Ben's question ("do you start fully deafferented in air, tune, then do
supported walking and tuning cycles while reintroducing connections?")
exposed that we didn't. His directive became the new methodology, and
the following were implemented and tuned in a staged curriculum:

### What was implemented (all off-by-default conditional topology)

1. **Heel/toe contact mechanosensors** (audit P1a): per-foot normal
   contact forces from MuJoCo contacts on calcn/toes geoms, normalized
   to BW fractions, fed as HEEL_c / TOE_c / LOAD_c input ports.
   HEEL_IN → RG-E exc + RG-F inh (S2W trigger, Conway/Hultborn 1987);
   TOE_IN → RG-E exc (late-stance prolongation). Gains: heel_rge,
   toe_rge (both default 0).
2. **LBIN stance-Ib group IN** (audit #15/P1a): per side, receives the
   stance-group IBEXC outputs (RG-E-gated Ib) and contact load → RG-E
   excitation (G["ib_rge"], default 0). Dominguez 2020 full text places
   these INs INSIDE the rhythm-generating layer.
3. **IaIN population** (audit #7/P1b): when G["ia_in"] > 0, Ia → IaIN →
   antagonist MN replaces the direct edge; PF_F1 → IaIN phase gate
   (Deng A6 0.5); RC → IaIN inh (Hultborn recurrent disinhibition)
   when Renshaw is on.
4. **Mutual Renshaw fix** (audit #12): RC↔RC between distinct pools,
   both directions (was one-directional). No self-synapse.

### The curriculum (staged tuning, Ben's methodology)

- **Stage 1 (air, deafferented)**: rhythm core only (drive, rg_adapt,
  desc_e, desc_f, rg_to_pf). Objective: 3×cycle-rises + swing knee
  flexion depth. 25 trials → best 104.6.
- **Stage 2 (air, afferented)**: + phase_reset_e/f, heel_rge, toe_rge.
  Seeded stage 1. 25 trials → best −68.1 kine. Heel 0.63, toe 0.38
  engaged.
- **Stage 3 (ground, full P1a/P1b)**: + ib_rge, ia_in,
  ankle_post_walk_trim. Seeded stage 2. 30 trials → best −72.5 kine
  (trial 16). Heel 0.91, toe 0.63, ib 0.68, ia_in 0.62, trim 0.63 —
  ALL new pathways tuned to nonzero.
- Final v11 run (full 22 s): stayed up, COM min 0.87, tilt −2..+18.6,
  hip −23.6..+58.3, knee −78.8..+12.7 (INSIDE RoM, deep swing
  flexion), ankle −89.9..−0.8, no hyperextension. All artifacts
  regenerated: ground_v11_curriculum.gif,
  hindlimb_style_v11_curriculum.png, opensim_overlay_gait_cycles.png
  (14 cycles), curriculum_final.npz.

### HONEST state vs OpenSim

Amplitudes are now in the right ballpark (hip 58 vs 43, knee 79 vs 70,
ankle range within bounds). The three named gaps are:
1. **Duty** ~0.14–0.2 vs 0.61 — the RG E-burst is still too short;
   P1a's ib_rge tuned to 0.68 but the prolongation mechanism needs
   more stance-Ib gain or a different integrator.
2. **Cadence** 0.46 Hz vs 0.81 — the curriculum model walks slower
   than OpenSim; rg_adapt tuned to 1.27 (was 0.88 pre-curriculum).
3. **Ankle** −90..−0.8 vs −10..+14 — still PF-dominant; the ankle trim
   (0.63) helped but tib_ant swing drive still loses to 9 PF muscles.
These are the architecture items, not scalar-tuning items. The next
levers are: (a) FSA-analytic seeding of PF/RG timing from the
backsolve, (b) transient onset-triggered reset INSIDE the RG layer
(now correctly placed per Rybak 2015 + Dominguez 2020), (c) MyoSim-
style compliant tendons.

## 2026-09-15 (LAMINATED architecture + curriculum results — HONEST)

The laminated IN-mediated architecture (InE/InF in RG, PF_IN_E/F in PF,
matching Shevtsova/Deng) was implemented and tuned through the
curriculum. RESULTS ARE MIXED:

- **Stage 1 (air, deafferented, 25 laminated trials):** best 57.4.
  Lower than pre-lamination 104.6 — expected: two synapses per
  inhibition path halve effective inhibition at the same g. The
  circuit oscillates but needs re-tuning for the new architecture.
- **Stage 2 (air, afferented, 25+25=50 mixed trials):** best −68.1
  (trial 16, pre-lamination). Laminated trials (25–49) did not beat
  it. The study mixes architectures — TPE learning is confounded.
- **Stage 3 (ground, 30+30=60 mixed trials):** best −72.53 (trial 16,
  pre-lamination). Same issue. The laminated trials scored worse,
  expected from the two-synapse inhibition path.

**HONEST CONCLUSION:** the laminated architecture is CORRECT per
Shevtsova/Deng and is the right long-term architecture. But the 30
laminated trials were insufficient for convergence on the new
architecture — the study is confounded by mixing architectures in one
TPE history. TO PROPERLY TUNE: delete the contaminated studies, re-run
the full curriculum (stage 1→2→3) with laminated-only trials and
ENOUGH budget (100+ trials per stage, since the two-synapse path
changes the effective gain scale). The alternative — reverting to
direct synapses — would be architecturally wrong per Shevtsova 2026
and would defeat the purpose of the exercise.

**NMF synergy analysis** (_synergy_nmf.py on backsolved activations
from the SO file): 3 synergies = 82% var (left stance, right stance,
transition — classic Ivanenko). 5 synergies = 92% (separates hip
abductors, knee extensors, plantarflexors). KEY FINDING: synergies
CROSS joint boundaries — no clean "hip PF" vs "knee PF" decomposition.
The natural grouping is task-based (stance push-off, swing
initiation), not joint-based. The current 4-cell phase-based PF is a
simplification; a synergy-based PF would be data-driven but requires
full redesign.

### Env incident (repaired)
conda graphviz install clobbered myo\python.exe — repaired via
--force-reinstall python=3.10.21; all pip pins survived; dot.exe OK.
WARNING: verify python.exe after any conda transaction in this env.

## 2026-09-14 (Ben Q): MuJoCo vs OpenSim MUSCLE FORCE comparison
## (answers "what spring" + "does force/torque match, R^2 + phase")

- `_muscle_force_compare.py` (+ `muscle_force_compare.csv`,
  `figures\muscle_force_compare.png`, dissertation copies +
  `muscle_model_appendix_draft.tex` in the Dissertation folder):
  MuJoCo (rigid-tendon actuator) vs OpenSim Thelen2003 (compliant
  tendon) forces along the subject01 walk, driven with IDENTICAL
  kinematics + IDENTICAL OpenSim SO activations (first attempt compared
  act=1 vs SO-driven — apples/oranges, fixed).
- RESULT: matched-activation R^2 **0.87-0.98 on 10/12 muscles**
  (tib_ant 0.95, per_brev 0.98, vas_lat 0.94, psoas 0.93, glut_max2
  0.93, tfl 0.93, med/lat_gas 0.89/0.88, soleus 0.74, semimem 0.57;
  rect_fem -1.15 and bifemsh -0.26 the outliers — patella reroute +
  scale offset), mean-force ratios within +/-22%, and **NO phase
  shift** (best lag 0-2 frames = 0-34 ms) -> Ben's tendon-slack
  hypothesis is NOT supported as a timing issue: the rigid-vs-compliant
  signature is SHAPE not phase (OpenSim compliant tendon carries soleus
  force through mid-stance while rigid MuJoCo follows activation dips
  — the missing series-elastic smoothing).
- The "spring" in my earlier note = patch_xml ligament surrogates
  (ankle/subtalar/mtp joint stiffness 10 N-m/rad added because the
  conversion lost OpenSim coordinate stops) — a JOINT spring, not a
  tendon element; it does not give the ankle a muscle-force path.
- MJCF muscle class (for the appendix): dyntype=muscle,
  dynprm 0.01/0.04 s (activation rise/fall), force = Fmax*[a*FL(Ln)*FV
  - FP(Ln)], Ln normalized to the actuator lengthrange, gainprm[2] =
  Fmax (RoM-peak-normalized to OpenSim, <=15% dev).

## 2026-09-14 (late): npz deletion (other session), RC wiring fixes,
## v10 overlay verdict on hip phase + ankle offset

- **npz deletion (other chat, NOT ours)**: 83 npz (986 MB) removed from
  the working tree incl. all *_best_trial*.npz captures; figures/gifs/
  jsons intact; blobs remain in git history (commit a8746f9). NOTHING
  of ours is blocked: every capture is REGENERABLE from its winner
  json via runner --fitted --bestN (fresh v10 winner run regenerated;
  spinal_run.npz now = v10 trial 56 run). Cleanup of OUR dead files:
  404-stub tutorial notebooks (bad cmd download), one-off probe/generator
  scripts, state dumps — done; utilities/logs kept.
- **RENSHAW WIRING FIXED (Ben caught both)**: (a) the DIAGRAM network
  (spinal_layers.py) had RC->RC SELF-loops — an autapse is not standard
  connectivity; (b) the real build wired each RC pair ONE-DIRECTIONALLY
  (order-condition in the pair loop) — lopsided. Both now MUTUAL between
  different pools (Hultborn's cat data; Deng A6 "RC->RC" is between the
  ext/flx RCs), no self-synapse. Panels diagram regenerated. NOTE:
  v10's numbers were tuned with the lopsided wiring — the mutual fix
  lands in the next retune (v11).
- **v10 overlay verdict** (regenerated from the fresh winner run):
  KNEE now close (swing flexion -75 at ~62 pct vs OpenSim -70 at 60
  pct; our stance starts -12 crouched and lacks the stance-extension
  wave). HIP still phase-inverted (ours flexes to +40 mid-"cycle"
  while OpenSim extends) — with E-duty at 0.08-0.2, the RG-E-anchored
  cycle is dominated by F-phase motion, so the inversion is largely a
  DUTY/ANCHOR symptom: until stance is >40 pct of the cycle, the
  anchored mean cannot look like OpenSim's. ANKLE: the trim tuned to
  its floor (0.097) yet the mean ankle still rides -52 deg — the
  offset is NOT primarily the standing POST tone; remaining carriers =
  passive ankle spring equilibrium + E2 PF drive + PF-side reflex
  load + contact. Next diagnostic: a stance-phase DF hold (angle-PD
  or DF-biased reflex gains) rather than more tone trimming.

## 2026-09-14 (evening): ankle set-point trim + TRANSIENT phase reset
## + ENV INCIDENT (repaired)

- **ENV INCIDENT (self-inflicted, repaired)**: the `conda install -n
  myo graphviz` transaction CLOBBERED `<env>\python.exe` (only
  pythonw.exe + python310.dll survived; conda-meta still listed
  python-3.10.21). Repaired with `conda install -n myo --force-reinstall
  -y python=3.10.21`; numpy 1.22.4 / scipy 1.9.3 / pip packages all
  intact, dot.exe survived. LESSON for the skill: on this conda stack,
  install graphviz with care and verify python.exe after ANY conda
  transaction in this env.
- **Ankle set-point trim**: `params.G["ankle_post_walk_trim"]`
  (default 1.0 = v9-identical) scales the POST bias of ankle_pf-group
  muscles toward 0 as drive rises — the soleus/tib_post standing tone
  held the -45 deg PF offset through gait; real soleus tonic EMG drops
  with locomotor drive. Search range 0.05-1.0.
- **TRANSIENT phase reset (the lever, implemented)**: PRESET_E/F now
  carry a FAST self-adaptation loop (PRESET -> PREA(tau 0.08) ->
  PRESET, gain PHASE_RESET["adapt_g"]=1.5), making them high-pass
  ONSET detectors: a sustained hip signal emits a brief pulse at its
  onset instead of a tonic bias (tonic <= 1 nA was proven inert vs
  ~4 nA DRIVE). Rectifying synapses pass only the onset pulse. Active
  whenever phase_reset gains > 0 (v9 winner reload now uses it by
  design). Audit (_phase_reset_audit.py re-run): signs PASS, mean
  effect +0.34 -> +0.10 mV as expected for high-pass.
- **v10 study** (`ground_walk_v10_transient`, optuna_walk_v10.py,
  seeded v9 winner, +ankle_post_walk_trim searched, `--best10`,
  v10_results.csv, 60 trials RUNNING).
- **v10 RESULT (60 trials)**: best **-65.375 (trial 56)**, reproduced
  BIT-EXACT (`runner --fitted --best10`) after TWO loader catches:
  (1) the `--best10` flag branch was LOST in my successive same-anchor
  edits (the chain silently ended at --best9 — a variation on the
  parallel-edit hazard; the state-dump diff + the missing
  "loaded best_walk_params_v10.json" print caught it), and (2) the
  ankle_post_walk_trim loader branch was missing (the JSON RULE
  again). Winner full-22s: knee -75, tilt 15.7, duty 0.19, stayed up;
  eval amplitudes hold (hip 45.4/43.3, knee 62.0/70.5, ankle
  25.4/23.1, knee_min -77.1/-69.7). The transient reset + trim
  re-tuned to parity with v9's tonic plateau; the trim winner value
  0.097 confirms the PF standing-tone diagnosis (near-zero standing
  soleus tone wanted during gait). Remaining: duty/cadence/hip-phase
  (architecture). Study resumable
  (`python optuna_walk_v10.py 60`).

## Goal

Two-level spinal cord network (McCrea–Rybak RG + PF) with proprioceptors for
each of Gait2392's 92 muscles, driving the MyoConverter MJCF of
`gait2392_simbody` in MuJoCo: standing → walking → standing, speed via
descending drive AND reflex-gain modulation, with PF→MN weights back-solvable
from OpenSim IK/SO activation patterns.

## Architecture (per side; 410 neurons total = 201/side + 8 shared, 376 inputs)

| layer | cells | notes |
|---|---|---|
| descending | DRIVE, POSTURE, BAL_PF, BAL_DF, BAL_TRK_EXT, BAL_TRK_FLX, BAL_LAT_R, BAL_LAT_L | MLR/postural surrogates, external inputs (the 4 later balance cells are why this is 410, not the stale 406) |
| RG | RG-E, RG-F + ADAP-E, ADAP-F | half-center: mutual inhib. + slow (0.9 s) self-adaptation; frequency rises with DRIVE; stance-biased drive split (E stronger) sets duty |
| PF | PF-E1, PF-E2, PF-F1, PF-F2 + PFA-* | forced mode (no own rhythm); (tau, adapt) shape multipliers stagger their windows; reciprocal inhib. on conflicting pairs |
| MN+afferents | per muscle: MN, Ia, II, Ib | Ia ~ tendon velocity (pure-signal: no resting tone — baseline drove constant reciprocal inhibition and crushed flexors); II ~ length (small baseline); Ib ~ force |
| Ib load sharing | IBEXC per stance group | Ib afferents → group IN, gated by RG-E → excites that group's MNs (stance reflex reversal / load sharing) |
| balance | BAL_PF/BAL_DF | pelvis-COM PD → ankle PF/DF + hip flex/ext MNs (ankle+hip strategy), fades with DRIVE |

Standing uses a **solved posture pattern**: static optimization at the
keyframe (torque rows = gravity − contact − passive, plus co-contraction
preload; activations regularized toward a physiological prior) injected as
per-MN bias currents (`POST_<muscle>`). This is the same machinery as the
planned IK/SO back-solve, applied to standing.

Ia/II/Ib gains are speed-modulated presynaptically (`params.MOD`, scaled by
DRIVE) and phase-gated from RG potentials — the "reflex modulation = speed
control" knob.

## Literature grounding (verified; mostly in Ben's Zotero)

- **RG+PF two-level**: Rybak et al. 2006 J Physiol (PF = multiple
  populations exciting synergist MN pools, forced by RG, low persistent Na);
  McCrea & Rybak 2008; Rybak et al. 2015 eNeuro review. The 2006 model
  explicitly anticipates multiple phase-shifted PF units for bifunctional
  muscles → Ben's "similar but phase-shifted PF" observation.
- **Speed**: increased drive shortens stance, swing ~constant (Danner 2017,
  Rybak 2024 eLife RP98841 — both RG-only, no PF). **Reflex gain modulation
  as speed control**: Bunz, Ijspeert, Schmitt 2026 Sci Rep
  10.1038/s41598-026-48509-z (offline+online, walk↔run); Ross/Rybak-lab mouse
  spindle-gain work (2018, in Zotero).
- **Precedent for exactly this stack**: SNS-Toolbox paper (Nourse et al.
  2023, Biomimetic 8:247) MuJoCo case study = rat hindlimbs, two-layer RG+PF
  CPG (Deng-style), PF shapes joint commands, MNs adjusted by Ia/Ib, speed
  set by RG inhibition level. Our network is the human-92-muscle analogue.
- **Synergies**: Ivanenko et al. 2004 (5 modules, speed-invariant);
  Kibushi et al. 2018 (speed carried by timing/amplitude); Avaltroni 2024
  spinal maps. → back-solve design in `fit_synapses.py`.
- **Closest competitor**: Di Russo, Ijspeert, Bouri 2023 J Neural Eng
  (CPG+reflex human neuromechanical model) — read before claiming novelty.
- Also in Zotero Sensory Afferent Database: Markin et al. 2016
  (large-scale PF→MN neuromechanical cat model), load-receptor review 2000,
  stick-insect multisensory control, air-stepping by tonic drive (2009),
  human balance NN papers (2017/2000).

## Verified status (this session)

- `check_rhythm.py`: half-centers oscillate; DRIVE=4 nA → period 0.844 s,
  E duty 0.29-0.33, left–right antiphase corr −0.74; PF windows alternate
  (E cluster ~0.08, F cluster ~0.72 of the cycle).
- Network builds/compiles/steps; runner closes the loop on the full
  92-muscle model with logging/summary/NaN forensics.

## BREAKTHROUGH 2026-09-10 (machine EB475WS4 session): the NaN blocker is DEAD

Root causes, in the order they were found (all fixed in `runner.patch_xml`
+ `apply_harness`; every fix is asserted at build time):

1. **boundmass 0.01 is too small at human muscle forces.** With 92 muscles
   pulling on the equality-coupled pathpoint DoFs the mass matrix goes
   singular ("Inertia matrix too close to singular" → NaN) at as little as
   0.3 whole-body co-activation. → `boundmass=0.1` (survives 1.0
   co-contraction for 6 s in diag_stab.py).
2. **The equality-coupled "conditional pathpoint" slide joints carried
   range limits sampled only around the keyframe** (MyoConverter artifact).
   The equality drives them exactly along their polycoef curves, which
   leave those boxes within 10-30 deg of knee flexion (14 of 36 coupled
   dofs outside their range at knee=10 deg, 31 at 30 deg) → range-vs-
   equality fight = huge constraint forces + jammed joints. → strip the
   ranges, `limited="false"` (NOTE: default class sets limited="true"
   explicitly, so dropping `range` alone errors out at compile) and add
   `armature="0.5"` to the followers (0.1 was not enough: rect_fem P3_y
   went singular at t=8.9 s with the CPG co-contraction).
3. **knee_angle driver range: NO FLIP — Ben was right (2026-09-10 PM).**
   The converter PRESERVED OpenSim's flexion-negative convention: real-
   actuator tests (unjammed followers, `_muscle_direction_test.py`) show
   gravity buckling the standing knee NEGATIVE (flexion), semimem/
   bifemsh/med_gas driving NEGATIVE (to -77 deg), vas_lat/rect_fem
   driving POSITIVE (extension). The converted range [-2.094, 0.1745]
   (= -120 deg flexion / +10 deg extension) is CORRECT as shipped.
   Moreover knee_angle ships `limited="false"` — the range is INERT, so
   an earlier "range flip" (2g, applied and then REVERTED the same
   evening) changed nothing dynamically. CAUTION: a kinematics sweep that
   rotates knee_angle WITHOUT letting the coupled translation dofs follow
   (the tibia's rolling center) mirrors the apparent foot path and
   "proves" the wrong convention — don't repeat that. Related: the
   audit_signs moment-injection method reports moment-arm signs that are
   SIGN-INVERTED vs real activation (docstring caveat added); the hip
   "axis flips" from the previous session were physical no-ops (axis
   negation only relabels qpos signs).
4. **collision="predefined" + 19 explicit ground <pair> lines**: predefined
   pairs IGNORE contype/conaffinity, so muting the ground geom does
   nothing (the first "air" runs dragged on an invisible floor). → patch
   surgery on the <contact> section: remove it entirely for suspended
   tests; keep only the 6 foot pairs on the ground (OpenSim gait2392 has
   no self-collision either).
5. **The pelvis rig springs were landing** (they always had been in the
   original code — a mid-session refactor briefly broke them, caught by
   _rig_check). Rig now also writes its damping (was accepted but never
   written before: undamped rig spring) and the 3 pelvis ROTATION hinges
   get springs too (150/25 ground, 400/40 air): with zero yaw damping the
   pelvis freely spun about vertical at 40-68 rad/s on the ground
   (foot-friction spin-up from left-right standing-solve asymmetries).
6. **Standing solve fixes**: always solve on a GROUND-ON model copy (in
   air the solve degenerates: nothing carries gravity → all activations
   <=0.2, wrong muscles); fit only DoF rows with nonzero muscle moment
   (pelvis translation rows carry the full body weight but no muscle can
   act on them → the ridge smeared everything to mush); stronger
   co-contraction preloads (hip -40, knee -60, ankle +40 N·m), lam 0.25.

Results after the fix stack (myo env on EB475WS4):
- **Deafferented + suspended in air: 22 s, zero warnings, 11 alternating
  RG_F bursts, COM stable on the rig.** knee range still only ~-1..1 deg
  in air: quads (9.2 kN total Fmax) vs hamstrings (5.5 kN) co-contract —
  at PF-cell current levels the net knee torque never crosses zero. The
  tonic DRIVE->PF term was halved (drive_to_pf 0.4→0.2) and PF reciprocal
  inhibition raised (2.0→3.0) to narrow the windows; W_PF_MN retuned
  (F1 knee_flex 1.2, F2 knee_ext 0, E1 knee_ext 0.25/ankle_df 0.10).
- **Ground walk with afferents (default runner config): 22 s, stayed up,
  11 bursts, contacts 3-11.** After gating the II-afferent BASELINE with
  the stance gate (runner.py; the ungated i0_ii=1.0 nA was a global ~0.2
  co-contraction floor on every extensor MN - found via diag_phase.py),
  joint motion appears: knee range 0..+17 deg, ankle -87..0 deg. NOTE ON
  SIGNS: with the true flexion-NEGATIVE knee convention, that +17 deg is
  EXTENSION/hyperextension-side motion, not flexion - the swing knee
  still does not flex in the walk; the ankle sweep is dominated by
  plantarflexion (toe-pointing). Same observation, correct labels.
- NEXT (in order): (1) swing-knee FLEXION (negative excursion) - the
  flexors have full range now; the co-contraction standoff is the
  remaining foe (phase table: quads 0.25 vs hams 0.27 in swing); consider
  stance-gating the quad drive harder + the PFA clamp 1.5->2.5;
  (2) ankle plantarflexion overshoot (-61 deg after the rig fix, was
  -87 = toe-pointing; raise the gated II floor / tib_ant swing drive or
  lower E2 ankle_pf);
  (3) COM bounce/drift - BAL gains (kx=150 saturates max_current=6 at
  com_x -0.3) and stance duty 0.29 vs human 0.6; (4) wean rig springs
  (tx tether 1500 N/m, pelvis rotations 150 N·m/rad, NEW 2026-09-10 PM
  after Ben watched the viewer: lumbar_extension/bending/rotation
  150 N·m/rad + hip_rotation 50 N·m/rad - the torso was flipping upside
  down over the free lumbar hinge and the legs free-spinning on their
  long axes; bounce min-COM improved 0.62->0.71) toward true balance;
  (5) then vestibular/ocular (Ben) and cerebellum/BG layers. Tools:
  diag_stab.py (plant A/B), diag_phase.py (phase-aligned activations),
  _muscle_direction_test.py (per-muscle real-activation direction test).
  runner.py --view now plays in REAL TIME with an initial side-view
  camera (orbit: left-drag, zoom: scroll, pan: right-drag, double-click
  a body to track it).

## 2026-09-10/11 night session (Ben staged tuning plan)
1. Literature (Ivanenko 2002 JN air-stepping): preferred cadence at 100 percent BWS = 36+-8 steps/min ~ 0.3 Hz cycles - air stepping is ~3x SLOWER than walking, not 5x faster. Extensor duty 53+-3 percent, joints near-sinusoidal.
2. Deafferented air, no interleg (--no-afferents --no-ground --no-interleg): WORKING. 21 s clean, 0.58 Hz, knee -75..+15 deg true swing flexion, hip -10..+60. Keys: decouple amplitude from frequency (low DRIVE shrank all amplitudes to sub-mV - doubled rg_to_pf/pf_to_mn/descend conductances, slow ADAP tau 1.9 s); subtalar/mtp ligament surrogate springs 10 N-m/rad (OpenSim stops lost in conversion - foot flopped to 128 deg); follower armature 1.0. Remaining: ankle PF-heavy, E-duty 0.16 vs 0.53.
3. Afferented ground + interleg (--drive 2.5): 21 s, stayed up, 6 bursts, COM y +-2-3 cm, adduction +-8 deg (BAL_LAT abductors engaging under load). PELVIS LIMBO: pelvis_tilt still 36 deg with 400 N-m/rad rig + IMU trunk controller - the IMU levels the TORSO (lumbar counter-tilts) but pelvis pitch is driven by SATURATED hip extensors (glut_max 0.94, semimem 1.0) pitching the planted-leg pelvis backward. Trunk muscles cannot fight hip torques - needs the IK back-solve.
4. New pieces: interleg toggle (--no-interleg); BAL_LAT_R/L stance-gated abductor strategy (mujoco -y = opensim +z); BAL_TRK_EXT/FLX IMU vestibular surrogate (torso up-vector PD -> ercspn/obliques). W_PF_MN: E1 knee_ext 0.10, E2 ankle_pf 0.35, F1 knee_flex 1.80, F2 ankle_df 0.45.
5. NEXT - IK/SO back-solve (Ben plan): OpenSim IK walking kinematics -> map to MuJoCo joints -> per-timestep NNLS back-solve along trajectory -> fit W_PF_MN, extract human hip/pelvis torque balance (pelvis_tilt +-5 deg), EMG ordering. Then ground duty 0.6, ankle balance, wean rig springs.

## 2026-09-11 optimization (Ben asked for fast + efficient tuning on ground)
SNS-Toolbox ships no optimizer; literature standard for expensive gait sims = Bayesian optimization (Calandra 2014 biped BO, Antonova 2016, Ryu/Geyer 2021 CPG optimality). Built optuna_walk.py: Optuna TPE, sqlite persistence (optuna_walk.db, resumable), runner --eval mode (12 s schedule, quiet, metrics dict: nan/dx/kz/tilt_max/knee_min/hip_amp/burst_r).
v1 (60 trials, 8 params incl. walk_drive/rg_adapt/desc_e/rg_to_pf/e2_pf/f1_df/f1_kf/post_kneext/kx): best 1.071 - STABLE upright shuffle: 26 s up at tether 1500/300/80 N/m (identical), 0.82 Hz, COM y +-5 mm, adduction +-9 deg, but knees 0..+13 ext-side (stiff-legged). Lesson: stability penalties dominate motion rewards - shuffle wins.
v2 (40 trials, added post_hipext param + motion-demanding objective: knee reward saturating at 35 deg): the v1 winner ENQUEUED as seed stayed best (1.091) - deep-flexion attempts all scored worse on ground. Honest result: at current fidelity, ground stability costs swing flexion; the missing piece is balanced activations (IK back-solve), not more parameter search. v2 objective shaping is ready for a rerun after that.
runner --best loads best_walk_params.json (drive included). runner npz now saves qfull (full qpos) + cfg air/ground; render_frames/render_video are npz-driven (render exactly the recorded run, config auto-detected, mj_fwdPosition so muscles follow bones). Ground media: ground_walk.gif, ground_walk_phases.png, ground_fig*.png in Dissertation folder. SELF-SUPPORTED: not yet - vertical rig still carries weight; BAL cannot hold sagittal/lateral alone; next lever = IK back-solve then wean rig.

## Notes for next session (Ben, 2026-09-11 evening)
1. fig4_gait_cycles.png OpenSim overlay: Ben suspects DOUBLE degree-conversion (our load_benchmark does np.degrees() on the .mot values - verify the freshly generated subject01_walk1_ik.mot units first: hip_flexion_r ~ +-0.3 means radians (keep conversion), ~ +-20 means already degrees (drop np.degrees)).
2. Rename the sim legend entry in fig4 to SNS sim mean (distinct from OpenSim IK).
3. ground_fig4_gait_cycles.png (dissertation folder) should ALSO carry the OpenSim overlay (plot_run was run before the benchmark file existed for that variant).
4. OpenSim 4.3 CLI works locally (Scale+IK in 90 s, artifacts committed in Gait2392_Robotbody); stay on 4.3 unless Python-API scripting is needed (then 4.6 into a py3.11 env, side-by-side, nothing on PATH).

RESOLVED same evening: (1) subject01_walk1_ik.mot says inDegrees=yes - values
are degrees, the np.degrees() in load_benchmark was the double conversion
(plot_run.py now honors the header flag); (2) legend renamed "SNS sim mean";
(3) rerun plot_run.py on a ground run and refresh the ground_fig* copies in
Dissertation\CPG_airstepping_figs (NOT yet done - do it on the post-v3 run).

## 2026-09-11 night session (EB475WS4): the IK/NNLS back-solve chain

Ben's staged chain EXECUTED: MuJoCo-vs-OpenSim muscle validation ->
per-timestep NNLS activation back-solve along the IK trajectory ->
W_PF_MN/W_POSTURE refit (limbo fix) -> v3 optimizer rerun with the
balanced pattern -> rig weaning.

New tools (spinal/): `bsolve_ik.py` (validation + back-solve, writes
bsolve_out.npz/.png/report), `fit_pf.py` (refit -> fitted_walk_params.json),
`wean_rig.py` (--rig-scale ladder), `optuna_walk.py` v3 (study
ground_walk_v3), probes `diag_force/diag_frontal/diag_knee/diag_trans.py`
(keep - they are the regression tests for the traps below). runner.py
gained --fitted (load fitted_walk_params.json), --rig-scale S (all RIG
stiffness*S, damping*sqrt(S); ligament surrogates untouched), --best
pf_gain support, and an eval `duty` metric.

VALIDATION (subject01_walk1_ik.mot: 121 frames 0.5-2.5 s, cycle 1.23 s,
duty 0.70; OpenSim reference via opensim-cmd 4.3 AnalyzeTool on
subject01_simbody.osim - NOTE 4.x has NO standalone StaticOptimization
tool, SO runs as an ANALYSIS inside Analyze, and the lengths file is
`*_MuscleAnalysis_Length.sto`):
- Muscle lengths: median r 0.89, all 78 muscles r>=0.81, 45/78 r>0.9.
  Large RMSE on thigh biarticulars (~5.8 cm at r=0.96) is the
  subject01-vs-generic-gait2392 SCALE offset, not shape error - r is the
  metric. Moment arms r 0.74-0.95 with 100% sign agreement after folding
  in MuJoCo's transmission sign (qfrc = -F dl/dtheta). Conversion sound.
- NNLS activations vs OpenSim SO: glut max/med 0.7-0.93, tib_ant 0.70,
  lat_gas 0.73; peronei/tib_post negative (subtalar geometry + scale);
  median r 0.22, mean torque residual 0.33. Group profiles:
  bsolve_groups.png.

THREE SILENT TRAPS (each verified by a probe script; do not relearn):
1. PELVIS SLIDES LOAD IDENTITY. The converter preserved OpenSim's
   coordinate VALUES (keyframe qpos[pelvis_ty]=0.95 <-> pelvis world z
   0.95; qpos[pelvis_tz]=+0.1 moves the pelvis to world y=-0.1). The z-up
   remap lives in the BODY FRAMES/AXES. A y/z swap in the loading code
   puts the pelvis 2 cm above ground and 1 m lateral (diag_trans.py).
   GRF VECTORS still remap os(x,y,z) -> mj(x,-z,y).
2. EQUALITY COUPLERS POISON INVERSE DYNAMICS. At deep knee flexion the
   polyfit pathpoint followers (fitted near the straight keyframe)
   generate ~890 N*m of spurious constraint wrench at the knee row
   (diag_knee.py). The bsolve ID model zeroes follower armature AND
   disables equalities (mjDSBL_EQUALITY): clean tree ID along the
   measured trajectory (follower dofs are excluded rows anyway; their
   boundmass inertia is virtual). FORWARD SIM STILL NEEDS the
   equalities - ID only.
3. ACTUATOR FORCE CHANNEL: the actuators are dyntype=muscle -
   data.act IS the activation and drives actuator_force under
   mj_forward; data.ctrl is the excitation target and is INERT under
   mj_forward (only feeds act dynamics in mj_step). ctrl=1 force 0,
   act=1 force -2655 N on soleus_r (diag_force.py).

Back-solve formulation: per IK frame, tau = mj_inverse(q, v, a) with
measured GRF subtracted via mj_jac at the CoP on calcn_r/l (forces
remapped, see trap 1); solve lsq_linear with the SO-style ridge
(lam = 0.05 * median column norm) and bounds [0,1]. Plain NNLS is
FORBIDDEN here - it dumps a~1000 into muscles whose FL~0 (near-zero
columns soak residual). 6 Hz zero-phase filtering of the kinematics
matches SO's lowpass_cutoff_frequency_for_coordinates=6. Frontal
signs (adduction/subtalar/list/rotation) are UNSEEABLE by muscle-length
matching - bsolve_ik sweeps them against ID+GRF consistency (frontal
residual 627 -> 30 N*m; found hip_adduction_l -1, subtalar_l -1,
hip_rotation_r -1, hip_rotation_l +1; signs stored in bsolve_out.npz).

REFIT (fit_pf.py): per functional group, NNLS of the back-solved
activation profile over the four recorded PF-cell phase windows (from
spinal_run.npz; only right-side PF cells are logged - windows are
identical per side). THE LIMBO FIX, at the source: hip_ext stance drive
stack 1.17 -> 0.27 vs human peak 0.34 (the old W_POSTURE hip_ext 0.22
alone was ~2/3 of the human PEAK, and E1+E2 stacked on top -> glut_max
0.94/semimem 1.0 saturation -> pelvis limbo). Output:
fitted_walk_params.json (full W_PF_MN/W_POSTURE + pf_gain=1.0).

V3 OPTIMIZER (optuna_walk.py, study ground_walk_v3): new pf_gain
dimension (log 0.5-8) - the back-solved weights carry honest human
amplitudes ~10x SMALLER than the hand-tuned table the network's gain
structure was tuned against (old F1 knee_flex 1.80 vs fitted 0.148); at
gain 1 the sim barely moves (hip_amp 0.85 deg). Plus the 5 phase
multipliers [0.5,1.8], drive/rg/kx ranges, and a duty-0.60 reward.
CAUTION: do NOT pass --fitted inside the optuna eval call - the reload
wipes the trial's mutations (verified: identical metrics with and
without gain). The final config reproduces with `runner --fitted --best`
(--best applies pf_gain to the whole table; the 5 knob values in the
json are already effective and are written AFTER the gain to avoid
double-apply).

WEANING (wean_rig.py): runs the ladder rig-scale 1.0 -> 0.6 -> 0.4 ->
0.25 -> 0.15 with --fitted --best; a stage passes when the full schedule
completes, COM height > 0.62, pelvis tilt < 35 deg; stops at first fail,
keeps wean_stage*_S*.npz per stage.

## 2026-09-11 night RESULTS (v3 winner + weaning verdict)

- v3 study ground_walk_v3, 40 trials: best score 1.521 (trial 37) vs
  v2's 1.091. Winner: pf_gain 0.52 (the optimizer went DOWN from the
  seed - human-shaped patterns need LESS brute-force gain, the opposite
  of the hand-tuned table's direction), drive 1.49, rg_adapt 0.83,
  desc_e 1.33, rg_to_pf 1.98, kx 220. Saved in best_walk_params.json
  (effective values + pf_gain + multipliers; reproduce with
  `runner --fitted --best`).
- v3 winner 12 s eval: no NaN, kz 0.877, tilt 27.1 deg, dx 0.151 m,
  9 bursts - but knee_min -0.6/hip_amp 0.9 deg on the SHORT window.
- FULL 22 s ground run at rig-scale 1.0 (wean_stage0_S1.npz): PASS -
  knee -22.4 deg REAL swing flexion, hip amp 37.7 deg, tilt 29 deg,
  kz 0.87, dx 0.15 m. The v2-era "knees 0..+13 ext-side stiff shuffle"
  is GONE at full rig.
- Weaning ladder: S=0.8 completes 22 s and keeps stepping (knee -22.7,
  hip 37.4) but tilt 36.7 deg = progressive lean; S=0.6 tilt 46.3 deg.
  VERDICT: support boundary S~0.8-1.0. The pelvis rotation assist is
  what the weaned rig misses - the BAL_TRK IMU levels the TORSO, the
  (now human-scaled) hip extensors hold the legs, but nothing yet holds
  the PELVIS pitch in the frontal-sagittal sense. Next lever per Ben's
  plan: ground duty 0.6 + ankle balance, then the pelvis-balance piece,
  THEN wean below 0.8.
- Figures: plot_run.py rerun on wean_stage0_S1.npz (spinal_run.npz
  currently holds that run) with the inDegrees fix - fig4 now carries
  the OpenSim IK overlay with the "SNS sim mean" legend; ground_fig1/3/4
  + ground_walk.gif refreshed in Dissertation\CPG_airstepping_figs.

## 2026-09-11 late night: publication circuit figures (draw_circuit.py rewrite)

draw_circuit.py is now a figure SUITE (Ben asked for paper figures of the
circuit): `python draw_circuit.py [--which core|full|weights|all]
[--source params|fitted|best] [--fmt pdf,svg,png]` -> figures/.
  - circuit_core: one side, DRIVE/POSTURE/POST_i -> RG+ADAP -> PF(+PFA)
    -> MN ellipses; reflex block; each MN pool carries ONE label (its
    dominant W_PF_MN entry); edges drawn for w >= 0.03; full table =
    the weights figure.
  - circuit_full: both sides (left ghost) + interleg commissurals as
    arcs over the top + BAL family; rig marked external.
  - circuit_weights: W_PF_MN x groups heatmap + W_POSTURE column.
CRITICAL PROPERTY: numbers are never hardcoded - --source composites
params.py -> fitted_walk_params.json -> (best) pf_gain x knobs exactly
like `runner --fitted --best` (verified by spot-check vs the jsons), and
each figure prints a gray source note. RERUN THIS AFTER ANY RETUNE
(e.g. the v4 trunk-Fmax fix) so figures never drift from the sim.
Gotcha recorded the hard way: matplotlib arc3 with a +x chord bulges
DOWN for positive rad - over-the-top arcs need negative rad. Old
spinal_circuit.png is superseded (its weights were the stale pre-refit
hand-tuned table). Copies of all 9 files (pdf/svg/png x 3) live in
Dissertation\CPG_airstepping_figs\ as circuit_*.

## 2026-09-11 night, addendum: Ben's knee/patella + muscle-count notes

- MOMENT-ARM DEFINITION (Ben): OpenSim's arm = dl/dtheta with the whole
  geometry following. MuJoCo's data.actuator_moment is the RAW Jacobian
  and does NOT propagate through the eq couplers - at the knee it misses
  the 36 vastii pathpoint followers, at hip_flexion 5 more. bsolve_ik.py
  now computes moment arms by CENTRAL DIFFERENCE dl/dtheta with the
  followers re-projected at every perturbed pose (fd_moments), which is
  the matching quantity AND the right B matrix for the back-solve.
  Keep this in mind anywhere knee moments matter.
- VASTII -> TIBIA VERIFIED (Ben: "the opensim model has no patella"):
  all 36 knee eq couplers drive the vas_med/vas_int/vas_lat pathpoint
  bodies (P3-P5) on both sides - the vastii actuate the knee/tibia
  through them (rect_fem rides vas_med's P4 since the patella reroute).
  Quantitatively (fd_moments vs OpenSim MuscleAnalysis knee arms):
  vas_med/lat/int mean arm -0.047 m vs OpenSim -0.045 m, same sign, but
  the small AC fluctuations anti-correlate (r -0.67..-0.81) -> the knee
  coupler polys deviate from OpenSim's true pathpoint trajectories away
  from the keyframe (polyfit artifact, same family as the range-limit
  and ID-wrench issues). Same likely explains the residual peronei/
  tib_post SO mismatches (subtalar couplers) and knee median arm r 0.31.
- "ALL 78 MUSCLES" RESOLVED: the comparison table = 92 actuators MINUS
  14 with Fmax (gainprm[2]) <= 5 N. Only 6 are the intentional prunes
  (quad_fem/gem/peri r/l). THE OTHER 8 SHIP AT Fmax = 1 N FROM THE
  CONVERTER: ercspn r/l, intobl r/l, extobl r/l, ext_hal r/l. That
  means the IMU trunk controller (BAL_TRK_FLX/EXT -> ercspn/obliques)
  has been driving muscles with 1 newton of capacity - the torso is
  held by the rig springs alone. diag_fmax.py prints the audit.
  RECOMMENDED FIX (needs Ben's go - it invalidates the current v3
  tuning): set the 8 Fmax values from stock gait2392_thelen2003 in
  patch_xml (like the prune patch), re-run the standing solve, refit,
  v4. Until then, trunk-control claims about ercspn/obliques are
  placeholders.

## 2026-09-12: v4 kinematics-matching campaign (Ben: "fine-tune until the
kinematics are similar to OpenSim")

New tool: kine_ref.py - reference cycle from subject01_walk1_ik.mot
phased by measured GRF onsets (right cycle 0.63-1.86 s = 1.23 s /
0.81 Hz, stance duty 0.61, knee_min -69.7 deg, hip range 43.3, ankle
range 23.1; careful: the mot parser index excludes the time column -
bit twice in one day). runner --eval now reports metrics["kine"] +
kine_score: cycle-normalized hip/knee/ankle shape RMSE (mean RG-E cycle
vs reference, offset removed) + peak-knee + range + duty errors, single
number, 0 = perfect match.

Campaign: optuna_walk v4 (study ground_walk_v4_kine) had a BROKEN
landscape - the -25 no-rhythm sentinel beat every genuine walker
(-40..-64), TPE collapsed onto that plateau in 18 trials; killed. v4b
(study ground_walk_v4b_kine) rescaled: no-rhythm = -65 (what a frozen
model truly scores), NaN -80+t_end, falls -10, tilt>40 -5. 60 trials ->
best kine_score -61.17 (trial 49; barely above baseline -63.6). Winner
params (saved in best_walk_params.json, study v4b): drive 2.59,
rg_adapt 0.88, desc_e 1.52, desc_f 1.08, rg_to_pf 1.96, pf_gain 2.18,
e2_pf 0.49, f1_df 1.78, f1_kf 1.79, e2_adapt 1.38, post_kneext 0.61,
post_hipext 0.79, kx 106. Full 22 s run: stayed up, 13 cycles at
1.18 Hz, hip -9..+25 deg, knee -14..+26, ankle -43..0, E-duty 0.27.
runner --best now also applies desc_f and e2_adapt when present.

HONEST DIAGNOSIS: scalar knobs on the current architecture CANNOT reach
OpenSim kinematics - the optimizer converged to ~-61 twice. The gap is
architectural, in priority order:
1. E-DUTY 0.27 vs 0.61 - the half-center + adaptation tops out ~0.3.
   Levers: hip-extension/loaded sensory phase-reset into the RG
   (Rybak-style Ia/II hip signals switching E->F), stance-biased PF
   windows, Ib load-sharing gain as stance prolonger.
2. SWING KNEE still extension-dominant (+26/-14 vs -70/+1) - the quad
   co-contraction standoff again; F1 knee_flex gains saturate. Needs
   phase-specific quad SUPPRESSION (F1 -> reciprocal inhibition onto
   knee_ext MNs), not more flexor drive.
3. ANKLE PF-dominant (-43 deg vs -9..+16) - POSTURE_OVERRIDE keeps
   soleus/tib_post tonically high (0.55/0.35) and the stand_frac floor
   carries it into gait; tib_ant swing drive (f1_df 1.78) loses.
4. CADENCE 1.18 vs 0.81 Hz - rg_adapt already at the range edge (0.88).
5. Cycle-to-cycle phase jitter: mean-cycle ranges collapse (hip 2.6 deg
   of a 34 deg raw range) - the rhythm is not stereotyped; sensory
   phase resetting (1) is also the fix here.
6. The 8 one-Newton trunk muscles (night addendum above) make BAL_TRK
   cosmetic - the Fmax fix is on the model list (needs Ben's go).
BEST CURRENT CONFIG = v4b winner (reproduce: runner --fitted --best);
the v3 winner remains better on the STABILITY-shaped objective (its
study is untouched in the db). kine_ref.py is the acceptance metric for
any future kinematics-matching campaign.

## 2026-09-12/13 night: Shevtsova-comparable circuit figure (--vclasses)

Source papers (Ben's Zotero, read from the local full-text caches):
**Shevtsova, Lockhart, Rybak, Magnuson, Danner, Smith, Poirazi 2026,
"Linking spinal circuit reorganization to recovery after thoracic spinal
cord injury", eLife 14:RP107480** (item TABLDVVC; the model = Danner-2017
four-RG architecture, Fig. 2 schematic + Table 1 weights), and **Rybak,
Shevtsova, Kiehn 2015 eNeuro review** (item BMB659JG; the V0/V1/V2a/V2b/
V3/dI6 class definitions). NOTE: the 2026 paper itself draws only
InF/InE (RG mutual-inhibition populations IniF/IniE), V2a, V0V, V0D,
V3-E/V3-F, Ini, InE1 + the LPN interlimb set — V1/V2b/dI6 come from the
2015 review framework; the 2026 connectivity is: RG-F→InF(0.4)→RG-E(−1),
RG-E→InE(0.4)→RG-F(−0.1), RG-F→V0D(0.7)⇢c-RG-F(−0.07), RG-F→V2a(1)→
V0V(1)⇢c-Ini(0.6)→i-RG-F(−0.04/−0.07), RG-F→V3-F(0.4)→c-RG-F(+0.03),
RG-E→V3-E(0.35)→c-RG-E(+0.02), V3-E→InE1(1)→c-RG-E(−0.045).

`draw_circuit.py --vclasses` (new flag) renders **circuit_full** with the
V-class annotation strip (4 columns: V0D(V0c)/dI6, V1/V2b, V2a, V3) plus
in-figure tags at the annotated elements. Glyphs unchanged (Ben's
conventions: open circle, open triangle = exc, filled dot = inh, ellipse
= muscle, Okabe-Ito tints). Output: figures/circuit_full.{pdf,svg,png},
copied to Dissertation\CPG_airstepping_figs as
**circuit_full_vclasses.\*** (the un-annotated circuit_full.\* stays).
Rerun with `--source best` after any retune, like the rest of the suite.

CLASS-BY-CLASS MAPPING (our element <-> V-class <-> Shevtsova 2026
element) — feeds the dissertation discussion:

| our circuit (gait2392 spinal SNS) | V-class (Rybak 2015) | Shevtsova 2026 element (weight) | comment |
|---|---|---|---|
| cross-side RG inhibition: RG_F↔RG_F (g 4.0), RG_E↔RG_E (g 2.0), all INHIBITORY | **V0D** (V0c in Shevtsova-2015 naming; dI6 is the 3rd inhibitory CIN class) | RG-F→V0D (0.7) ⇢ contra RG-F (−0.07) | alternation; ours is a direct lumped synapse, theirs via a CIN population |
| RG-E↔RG-F half-center mutual inhibition (g 4.0, direct) | **V1/V2b** (ipsilateral inhibitory; also Ia-IN/Renshaw family) | RG-F→InF (0.4)→RG-E (−1); RG-E→InE (0.4)→RG-F (−0.1) | theirs is routed through explicit IniF/IniE populations, ours is a direct lumped synapse |
| PF reciprocal inhibition E2↔F1, F2↔E1, E1↔F1, E2↔F2 (g 4.0) | **V1/V2b** | (same InF/InE family, pattern layer) | keeps PF windows mutually exclusive |
| Ia reciprocal inhibition onto antagonist MN pools (g 0.4) | **V1** (classic Ia inhibitory interneuron) | (not drawn in 2026 fig) | textbook reciprocal pathway |
| RG→PF excitation (g rg_to_pf) and PF→MN excitation (g pf_to_mn × W_PF_MN) | **V2a** (ipsilateral excitatory, Chx10) | RG-F→V2a (1)→V0V (1); V2a rhythmically recruits with frequency | in the full 2015 models V2a also excites MN pools — our PF→MN fan is the analogous excitatory relay |
| **NO analog — every cross-side connection we have is inhibitory** | **V3** (commissural excitatory, synchrony) | RG-F→V3-F (0.4)→c-RG-F (+0.03); RG-E→V3-E (0.35)→c-RG-E (+0.02); V3-E→InE1 (1)→c-RG-E (−0.045) | **we cannot express left–right SYNCHRONY gaits**; adding an excitatory commissural class is the architectural answer if ever needed |
| ADAP-E/F burst-termination loop | (none — intrinsic INaP adaptation in Shevtsova) | slow inactivation of the centers | ours is an explicit slow inhibitory interneuron |
| v5 PRESET_E/F phase-reset interneurons (hip ext/flex signals → RG) | (none of the named classes; group I/II position/velocity pathways) | afferent phase reset (Rybak 2006b; hip-reset experiments) | new in v5, gains default 0 |

Also noted for the dissertation text: in the 2026 model the FLEXOR center
is the driven/bursting one and frequency rises by SHORTENING the
extensor phase — the opposite drive allocation of our RG (our DRIVE
biases E, and our problem is a too-SHORT E phase, duty 0.27 vs 0.61).

## 2026-09-12/13 overnight: v5 sensory phase-reset (+ the --best pf_gain
## bug it caught) and phase-3 swing-knee quad suppression

Implements the v4 diagnosis's lever #1 (hip afferent phase-reset into the
RG: duty 0.27 vs 0.61, cadence 1.18 vs 0.81 Hz, cycle jitter) and, if
time allowed, lever #2 (swing-knee quad suppression).

### Implementation (all off-by-default; v4b behavior preserved exactly)

- `params.G` gained `phase_reset_e`, `phase_reset_f` (default 0.0) +
  `PHASE_RESET = dict(inh=1.0, stance_gate=(0.3,0.7))`; `params.TAU`
  gained `preset=0.04`. Phase 3 added `G["f1_kneext_inh"]` (default 0.0).
- `build_network.py`: per side, input ports HIP_EXT_SIG/HIP_FLEX_SIG drive
  PRESET_E/PRESET_F interneurons; PRESET_E excites RG-E and inhibits
  RG-F (prolongs stance), PRESET_F excites RG-F and inhibits RG-E
  (triggers swing). Phase 3: PF_F1 -> KINH inhibitory IN -> every primary
  knee_ext MN of the side (phase-gated by F1 itself = swing-only).
  **CONDITIONAL TOPOLOGY**: the new neurons/synapses are built ONLY when
  the corresponding gain > 0. At gains 0 the compiled network is
  byte-identical to v4b. REASON (learned the hard way, reggate_v5_0/2):
  zero-conductance synapses still change the dense-matrix BLAS summation
  order, and the chaotic 12 s sim amplifies the 1e-16 drift into visible
  kine_score shifts (~0.004-0.6 depending on what else changed).
- `runner.py`: computes both signals from EXISTING afferent data —
  `ext_sig = clip(-mean(len_norm[hip_ext group r/l]), 0, ...) * (0.3 +
  0.7*stance)` (hip-extensor group mean normalized length INVERTED: it
  rises as the hip extends; stance-gated II-style) and `flex_sig =
  clip(-mean(vel_norm[hip_flex group]), 0, ...)` (hip-flexor group mean
  shortening velocity — positive in flexion). Group membership is
  primary OR secondary (hamstrings' hip_ext arm, rect_fem's hip_flex arm
  count). Feeds the ports (gain * signal), logs RAW signals to
  log_neuro columns 14-17 (HIP_{EXT,FLEX}_SIG_{r,l} — appended, old
  column indices unchanged for kine_ref/plot_run). CLI: `--phase-reset
  E F`, `--kneext-inh X`, `--best5`/`--best6` (load best_walk_params_
  v5/v6.json through the same loader as --best).
- v5 signal magnitudes in gait: ext_sig ~0.1-0.4 nA, flex_sig ~0.3-1 nA
  at swing onset; RG-E/F receive ~1.5-4 nA from DRIVE — a weak
  perturbation at gain <= 1 (matters for reading the sweep below).

### THE REGRESSION-GATE SAGA — a real pre-existing bug found (fix first,
### as Ben's orders said)

Gate requirement: `runner --fitted --best --eval` must reproduce v4b
trial 49 (−61.1741850610298). It initially returned −61.82.

1. First suspect (fp drift from added topology) — fixed by conditional
   topology, but the gap remained.
2. State-dump hook added (RUNNER_DUMP_STATE=<file> env var dumps every
   mutable param right after arg parsing) and the study path
   (optuna_walk set_params) vs the --best path diffed numerically.
3. ROOT CAUSE: `--best` tested `if "pf_gain" in best` where
   `best = json["params"]` — but pf_gain lives at the JSON TOP LEVEL.
   The gain-scaling branch NEVER RAN: **every `--fitted --best`
   "reproduction" since the v3 pf_gain feature silently evaluated a
   pf_gain=1.0 config.** The recorded v4b full-22 s numbers (12 bursts,
   E-duty 0.27, knee −14..+26, kine −63.3) and the v3 weaning runs
   describe gain-1.0 configs, NOT the true study winners. Historical v3
   "verified: identical metrics with and without gain" was wrong.
4. Fix: read `doc.get("pf_gain")` from the document root; scale ONLY the
   fitted-file entries (the fitted JSON lacks the 4 trunk keys — E1/E2
   trunk_ext, F2 trunk_flex, W_POSTURE trunk_ext — scaling params
   defaults for those was also wrong, though inert: the 8 trunk muscles
   ship at 1 N and never move the sim). Same composite fix applied in
   draw_circuit.effective_tables (weights figure numbers change slightly
   — regenerate+recopy done 2026-09-13 night).
5. GATE NOW PASSES BIT-EXACT: `runner --fitted --best --eval --drive
   2.5868` (the study eval's 4-dp-rounded drive) = −61.1741850610298 on
   the v5 code at gains 0. With the json's full-precision drive the same
   config gives −61.181 (Δ0.007 = the :.4f rounding in optuna_walk.py's
   eval call — NOT chaos; drive sensitivity is ~0.005 per 3e-5 nA).
   v5/v6 evaluate trials with repr(drive) so their jsons reproduce
   bit-exactly via --best5/--best6 with no --drive override.
6. DETERMINISM verified: identical reruns are bit-identical
   (−61.817683226763805 twice, det_a/det_b).

### Sign audits (dynamic, no static actuator_moment anywhere)

- `_phase_reset_audit.py` = PASS. Part A: tonic 1 nA into each new port
  on the compiled network — EXT_SIG: RG-E +0.34 mV / RG-F −0.36 mV;
  FLEX_SIG: RG-F +0.21 / RG-E −0.10; signs correct at gains 1 and 2.
  Part B (static MuJoCo kinematics — lengths are pure transmission
  kinematics, the static-moment trap does not apply): ext_sig peaks at
  hip extension (+0.302 at −30 deg, monotonic to −0.492 at +60 deg
  flexion); flex-velocity signal positive at EVERY angle (4-8 nA per
  1 deg/s of flexion — length-dependent sensitivity, not a sign error;
  an earlier FAIL here was a bug in the audit's own condition).
- `_kinh_audit.py` = PASS: with f1_kneext_inh 1.5, KINH_r peaks 3.3 mV
  during PF-F1 bursts and pulls MN_vas_lat_r from 1.02 mV (gain 0) to
  −1.30 mV during swing (drop 0.99 -> 2.59). Real suppression, bounded.

### Hand sweep (v5_sweep.csv, 16 rows, 4x4 grid at the frozen v4b point)

At gains (0, 0.25, 0.5, 1.0)² with ALL other knobs held at the v4b
winner: kine −61.18 -> −61.62 (no improvement), duty 0.467 -> 0.38-0.46
(flatt-to-down, AWAY from 0.61), knee_min pinned ~+11 deg, cadence
1.067 unchanged, all stayed up. HONEST READING: a <= 1 nA gated input
cannot move a RG running on ~4 nA of saturated DRIVE — the sweep is the
LEAST favorable testbed (all other knobs frozen); the joint search is
the real test. Brief's rule stands: if the study plateaus at the 2.0
gain boundary, widen ONLY the phase-reset ranges and continue.

### Study ground_walk_v5_phase (optuna_walk_v5.py; sqlite, resumable)

Search space = v4b's 13 dims + phase_reset_e/f in [0, 2]; seed = v4b
winner (trial 0, gains 0, scored −61.181 eval = matches the sweep and
the fixed --best path bit-for-bit); per-trial rows -> v5_results.csv
(trial, gains, kine_score, duty, knee_min, rmse_hip/knee/ankle, cadence,
tilt_max, stayed_up); every NEW GLOBAL BEST triggers the full 22 s run ->
v5_best_trial<N>.npz + a run=full22 CSV row. walk_drive is runner-local
— full22 passes repr(drive) explicitly (bug caught before it bit).
Startup crashes fixed: study.best_value raises on queued-only studies
(guard: completed trials only) and the enum is optuna.trial.TrialState.

RESULTS (150 trials, COMPLETE): best kine_score **−60.987 (trial 137**,
phase_reset_e 0.222, phase_reset_f 0.219, drive 2.52, rg_adapt 1.28,
desc_e 1.64, desc_f 0.86, pf_gain 2.40, e2_adapt 1.61, post_hipext
0.99). Saved best_walk_params_v5.json; `runner --fitted --best5 --eval`
reproduces −60.98655543758134 (verified). Gain trajectory across the
top trials: TPE first latched ~0 gains, ended accepting ~0.22 — the
+0.19 over the v4b seed came mostly from drive/rg_adapt/e2_adapt
fine-tuning, NOT from the phase-reset. 9 new-global-best full-22 s
captures (v5_best_trial{0,20,42,59,63,105,124,131,137}.npz), all stayed
up; full-22s kine −62.6..−63.1, duty 0.44-0.48, knee_min ~+10 (still
extension-dominant), cadence 1.1-1.2 Hz — the duty/cadence/knee PLATEAU
DID NOT MOVE. HONEST VERDICT: tonic stance-gated position + velocity
inputs into the RG (this formulation) is effectively INERT at gains the
optimizer will accept — lever #1 needs a different formulation
(transient phase-ONSET reset — a flexion-velocity PULSE triggering the
F transition, not tonic excitation — or the v4 list's alternative:
stance-biased PF windows / Ib load-sharing as stance prolonger). The
no-rhythm trials clustered at phase_reset_e ~0.6-0.75 with weak random
drive: the stance-gated ext-signal is positive feedback that can LATCH
RG-E on (worth remembering as a stability constraint on any future
sensory-RG pathway). Per Ben's widen rule: gains never pinned at the
2.0 boundary — they were AVOIDED, so widening upward was not justified
and the ranges were left as briefed.

### Phase 3 (lever #2): ground_walk_v6_kneext (optuna_walk_v6.py)

v5 space + f1_kneext_inh in [0, 2]; KINH topology conditional (above);
seeded with the v5 winner (falls back to v4b defaults if the v5 json is
absent); writes best_walk_params_v6.json (runner --best6), CSV
v6_results.csv, full22 captures v6_best_trial<N>.npz. Launch ONLY after
v5 finishes (seed needs the v5 json; and both studies write
spinal_run.npz — never run concurrently).

RESULTS (110 trials, COMPLETE): best kine_score **−59.426 (trial 106**,
f1_kneext_inh 0.590, phase_reset_e 0.388, phase_reset_f 0.035, drive
2.85, rg_adapt 1.10, desc_e 1.75, desc_f 0.83, pf_gain 2.13, e2_adapt
1.90). Saved best_walk_params_v6.json (`runner --fitted --best6`
reproduces; verified). 10 full-22 s captures
(v6_best_trial{0,15,16,17,21,54,86,101,102,106}.npz), all stayed up;
winner full22 kine −60.591, duty 0.473, knee_min +10.1, cadence 1.4 Hz,
tilt 35.3. THE LEVER THAT WORKED: unlike the phase-reset gains (v5),
TPE ENGAGED the suppression dimension immediately and monotonically —
top trials all carry f1_kneext_inh 0.4-0.6, and the improvement
(−60.99 -> −59.43 eval) is carried by the knee-shape RMSE term
(20.9 -> ~19.5) with the full-22 s score jumping to −60.4..−60.8 (v5's
best full run was −62.6..−63.1). Swing-knee quad suppression is the
real architectural lever so far; knee_min is still ~+10 deg
(extension-dominant mean cycle) and cadence rose to ~1.4 Hz (ref 0.81)
— the NEXT bottleneck is now clearly the RG duty/cadence structure, not
the knee musculature.
