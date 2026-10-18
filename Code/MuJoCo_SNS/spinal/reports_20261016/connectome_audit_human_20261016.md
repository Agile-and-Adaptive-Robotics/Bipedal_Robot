# Connectome audit — HUMAN-body walker nets (campaign 2026-10-16, stage 2+3)

Auditor stage: "connectome auditor (human)". Everything below was RUN on
easteregg2 (`D:\Anaconda\envs\myo\python.exe`, cwd
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal`) unless marked LOCAL
(laptop myoconv env, files hash-identical to egg — SHA1-verified for all
13 audited sources before any edit).

## What ran (exact commands)

1. Edge inventory at the s3k winner config + default build + Ben's three
   exports:
   `D:\Anaconda\envs\myo\python.exe _audit_human_20261016.py`
   → `reports_20261016/audit_human_edges_20261016.txt` (1512 lines; full
   dump of every compiled connection, both builds, grouped + rule-filtered).
   The winner config is applied exactly as production does it
   (`_curriculum.set_stage(4, v10 multipliers + renshaw 0.5 +
   reports_20260923/s3k_trial34_full_params.json == curriculum_stage3.json
   params`, verified equal), so the audited net is the one the runner builds.
   `no_cross=1.0` in the winner therefore APPLIES ITS OVERRIDES:
   c1_gain→0.1, v3_gain→0.0, ia_f_contra_f→0, contra_swing→0, contra_kinh→0,
   and `set_stage` forces `renshaw=0.5`.
2. Duplicate-edge probe: `_dup_probe_20261016.py` (output in this report).
3. Anchor evals pre/post change:
   `_s3k_anchor_eval_20261016.py A|B` → logs
   `anchor_pre_change.log`, `anchor_post_change_A.log`,
   `anchor_post_change_B.log`.
4. Defaults + variant gates:
   `reports_20260925/audit/audit_defaults_gate.py` (the repo's canonical
   gate) + default-G AARL_NET=w2lvar/syn6 builds, before and after the
   drive change.
5. L/R drive topology probe at the winner config:
   `_lrdrive_probe_20261016.py` → `lrdrive_probe_20261016.txt`.

## Build counts (measured, not read)

| build | counts (n/i/syn) | note |
|---|---|---|
| DEFAULT (module-default G) | (410, 376, 1186) | canonical gate PASS (pre and post drive change) |
| s3k WINNER, full_rules=1 (production topology) | (878, 382, 8396) | incl. 92 RC cells (renshaw 0.5 forced) + 2 KINH cells |
| s3k winner + 4 drive knobs ON | (882, 386, 8492) | +4n/+4p/+96s, drive wiring only |
| w2lvar default-G | (888, 382, 7866) | == pre-change == 2026-10-03 rulings record |
| syn6 default-G | (794, 382, 3033) | == pre-change == rulings record |

## Rule verdicts (s3k production config = curr_s3k_nocross t34, full_rules=1)

### (a) RG+PF feedback PER PF LAYER — **NOT MET in production**
- The winner runs `joint_pf=0` → the PHASE-cell PF: 12 PF neurons = 8
  phase cells (PF_E1/E2/F1/F2 × r,l) + 2×(PF_IN_E, PF_IN_F). One
  ipsilateral RG pair drives its 4 phase cells at rg_to_pf=2.6005
  (8 RG→PF edges total; dump §(a)).
- Ben's drawing (ben_rules_20260924.json): HC-RG-E-ipsilateral →
  HC-PF-E_hip / HC-PF-E_Knee / HC-PF-E_Ankle each at **0.1** with a
  SEPARATE IN-PF lamination per micro-layer (2.749). I.e. three
  joint-specific PF micro-layer pairs per side — none exist in the
  production build. Joint differentiation in production comes only from
  per-muscle W_PF_MN weights on the shared phase cells (340 PF→MN edges
  broadcast to all joints) — exactly the "per joint broadcast" the rule
  rejects.
- The compliant topology EXISTS but is not what production runs:
  `G["joint_pf"]` (6 joint HCs, default 0, off in the winner), w2lvar
  (HIP PF + merged KNEE+ANK PF pair per side — a documented lumping of
  two micro-layers), syn6 (6 synergy PF layers, Eq-18 weights).

### (b) heel = stance reset AT the PF layer; toe → TOEDF only — **NOT MET in production**
- Production heel edges (dump §(b)): HEEL_s→InE_s exc 0.3753
  (heel_rge) + HEEL_s→InF_s **inh** 0.3753 + HEEL_s→PF_E1/E2_s exc
  0.2435 (the ib_e_central shared pathway, because foot mechanosensors
  ride the extensor Ib route). NO heel→PF_IN edge (heel_pf_layer=0 in
  the winner; the key exists, default 0).
- Production toe edges: TOE_s→InE_s exc 0.1457 + 0.2435 (double
  synapse, see finding N4) + TOE_s→PF_E1/E2_s exc 0.2435. **TOEDF does
  not exist** (toe_df_inh=0; 0 edges into or out of TOEDF) — toe is NOT
  dorsiflexion-inhibition-only; it is an RG-lamination extensor
  exciter, the opposite specialization.
- Ben's drawing: heel → ipsi InE 0.5 + ipsi PF E-lamination INs
  (hip+Knee+Ankle, 0.5 each) + CONTRA InF 0.5; toe → toe IN 1.0 →
  IN-PF_dorsiflexion_inhibit 5.0 → HC-PF-Dorsiflexion (inh 2.749), and
  toe touches NOTHING else.
- Compliant implementations exist off-by-default (heel_pf_layer,
  toe_df_inh, heel_in_f_exc, ia_pf_f, ii_pf_f — all 0 in the winner);
  w2lvar implements the drawing reading (ipsi InE + contra InF + both
  PF E-INs 0.5; TOEDF 5.0) with the documented KNEE+ANK merge side
  effect. The s3k stack cannot express the contra-InF heel edge without
  a new key (2026-10-03 ruling F5).

### (c) proprioception AUTOGENIC per muscle — **MET (reflex arcs); group terms are the drawn motifs**
- Direct afferent→MN edges at winner: autogenic 92 (Ia→own MN 0.6),
  non-autogenic **0**. The default build still carries 436 direct
  Ia→antagonist-MN inhibitory edges (ia_in=0 path); production
  (ia_in=0.825) routes reciprocal inhibition through per-muscle IaIN
  (Ia→IaIN 0.6; IaIN→antagonist MNs 0.4; IaIN↔IaIN 0.5, 436 directed
  edges post-2026-10-03-fix).
- II: II→IIX(1.0)→own MN 0.4; II→IIIN(1.0)→antagonist MNs 0.4 — the
  drawn ii_exc/ii_inh motifs (gains differ from drawing: F6 ruling).
- Ib: Ib→IBIN(1.0)→own MN 0.35 inh (autogenic); Ib→IBEXC(0.5, RG_E-gated
  at 1.0)→group extensor MNs 0.6 — the drawn ib_rev stance reversal.
- The central projections (ib_e/ia_f/ii_f/ii_e_central) send per-muscle
  afferents onto the side's SHARED central cells (PF_E1/E2, RG_E/F,
  InE/InF) — same shape as Ben's drawing (SN-Ia/SN-II → shared
  HC-RG-F/HC-PF-F nodes), not per-joint.

### (d) Ib load handling "PF-E 0.5, RG-E gate 0.1" — **NOT LITERAL in production**
- Ben's drawing ib_load motif: Ib grp → RG-E 0.5, Ib grp → IN-InE 0.5,
  Ib grp → HC-PF-E 0.5; and the per-layer RG-E→HC-PF-E drive gate 0.1.
- Production instead: per-muscle Ib→{PF_E1,E2, RG_E, InE} all at
  ib_e_central=**0.2435** (winner-tuned, not 0.5); RG→PF drive at
  rg_to_pf=**2.6005** (not 0.1); the stance reversal is a separate
  IBEXC IN (Ib→IBEXC 0.5; RG_E→IBEXC gate at **1.0**, not 0.1;
  IBEXC→MN 0.6); plus LBIN→RG_E at ib_rge=0.6741 (Dominguez pathway,
  not in the drawing).
- The literal dress exists in the variants: w2lvar family-sums the
  drawing's 0.5 (0.5/27 per-muscle onto the four E targets) and syn6
  uses 0.5 verbatim (both per WIRING_RULINGS_20261003 + its probes).
- Reading: rule (d)'s "unchanged" held for the 09-24 ruling round; the
  s3k production stack keeps its own tuned conductances and adds the
  IBEXC/LBIN machinery. Flagged as a doc-vs-code gap, not a bug.

### (e) afferents wired DIRECTLY into RG half-centers — **PRESENT (4 knob families, 184 edges), but the "original W2L has zero" premise is CORRECTED**
- Direct afferent→RG edges in the production build (dump §(e)): Ib→RG-E
  ×54 @0.2435 (ib_e_central); Ia→RG-F ×38 @0.3986 (ia_f_central);
  II→RG-F ×38 @0.3839 (ii_f_central); II→RG-E ×54 @0.093
  (ii_e_central). Contralateral afferent→RG is ABSENT in production
  (ia_f_contra_f zeroed by no_cross). Heel/toe reach the RG only
  through the InE/InF lamination (no direct heel/toe→RG edges at
  full_rules=1 — the laminated ruling IS honored); LBIN→RG-E (0.6741)
  is disynaptic-plus from Ib via IBEXC.
- DOCTRINE STATUS (important): the premise "in the original W2L ZERO
  afferents reach RG or PF half-centers" was an audit error CORRECTED
  2026-10-03 — WIRING_RULINGS_20261003.md ledger row 1 records Ben's
  screenshot verification that afferents DO reach the RG/PF
  half-centers in the original Deng model, and Ben's own master drawing
  carries SN-Ia/SN-II → HC-RG-F / HC-PF-F / IN-InF at 0.5 (edges
  re-dumped this round). So the four central knobs are known,
  Ben-sanctioned deviations from the OLD reading only. They remain
  flagged here because the ask names them; Shevtsova's abstraction
  (ben_shevtsova export) indeed has NO afferents at all.
- Cross-check per export (re-dumped this run, on egg):
  - ben_rules_20260924.json: SN-Ia_2843/SN-II_2844 → {HC-RG-F_2811,
    IN-InF_2810, HC-PF-F_2817} ×0.5 (direct afferent→HC edges present);
    heel→{IN-InE_2806, IN-PF_2810/2812/2814, contra IN-InF_2813} ×0.5;
    toe→IN-PF_dorsiflexion_inhibit 5.0; Ib grp→{RG-E, IN-InE, HC-PF-E}
    ×0.5; RG-E→HC-PF-E 0.1 per micro-layer.
  - ben_shevtsova_20260924.json: zero afferents; "brainstem gamma"→
    RG-E_l/r 0.5, "brainstem alpha"→RG-F_l/r 0.5 + alpha inhibits
    V0V/V0D 0.5 (the alpha/gamma naming precedent for the drive split).
  - ben_shinohara_20260924.json: flexor ii/ia →{RG-F, IN-F, PF-F} 0.53
    (eq10); extensor Ib →{RG-E, IN-E, PF-E} 0.16 (eq11); supraspinal→
    RG-E 0.15 / RG-F 0.02.

## New findings this round (not in WIRING_RULINGS_20261003.md)

- **N1 (HIGH, verified by probe): stacked duplicate IBEXC→LBIN
  synapses.** `build_network.py` `_wire_muscle` adds
  IBEXC_<group>_<side>→LBIN_<side> (0.5) once per MUSCLE in the group
  with no dedup guard (the 2026-10-03 audit's lazy-creation bug family,
  opposite direction). Probe output: 54 synapses on 8 distinct pairs —
  ankle_pf ×9, hip_ext ×7, hip_abd ×6, knee_ext ×4 per side — so the
  effective LBIN drive conductance is group-size-weighted (4.5 / 3.5 /
  3.0 / 2.0 vs the intended flat 0.5; trunk_ext ×1 = 0.5). The stance
  load signal into RG-E is therefore unintentionally scaled by extensor
  group size. NOT fixed here (it would change the tuned winner's
  dynamics); fix = add the edge once per group in the `gname not in
  idx` creation branch. Same code path serves w2lvar (unconditional
  LBIN edge there is single per group per its own builder — not
  re-verified this round; flagged unconfirmed).
- **N2 (MEDIUM, verified): the 2026-10-03 rulings doc's s3k winner
  topology count (786, 382, 3846) under-counts the PRODUCTION path.**
  `set_stage` (used by the curriculum objective AND _replay_best)
  forces `P.G["renshaw"] = 0.5` unconditionally and the v10 base
  multipliers bring `f1_kneext_inh=1.584` — so the real production
  build is (878, 382, 8396): +92 RC cells, +4416 RC-family synapses
  (incl. 4140 RC↔RC), +2 KINH cells with KINH→26 MNs at 1.584/0.929.
  The rulings' "Renshaw cells are absent from the s3k production
  winner" is wrong for every production entry point.
- **N3 (MEDIUM, verified by eval): the classic s3k t34 config NO LONGER
  WALKS at full_rules=1 under current physics.** Flavor-B anchor
  (production topology, winner gains/drive) scores **kine -25.0**
  (pre-change log == post-change log, bit-exact), vs its 2026-09 era
  db value -159.45 and vs the mutfix retune -192.21. The -161.56754173676563
  number that reproduces bit-exact is the **full_rules=0** flavor (the
  AGENTS 2026-09-26 gate recipe), NOT the production full_rules=1
  build — the "classic s3k anchors bit-exact at -161.567" shorthand
  conflates the two configs. Phase-4 bootstrap should tune from the
  mutfix winner, not classic t34.
- **N4 (LOW, by construction): parallel same-pair synapses from
  stacking knob pathways.** HEEL_s→InE_s exists twice (heel_rge 0.3753
  + ib_e_central 0.2435) and TOE_s→InE_s twice (toe_rge 0.1457 +
  0.2435) — conductances add; each is a separately-ruled pathway, so
  this is documented behavior, not a bug, but effective heel drive
  into InE is the sum.

## L/R alpha/gamma hand drive — IMPLEMENTED, default-off bit-identical

Four default-0 knobs per model (exactly the campaign spec): `gamma_l`,
`gamma_r`, `alpha_l`, `alpha_r` (+ conductance key `alpha_to_mn`,
default 0.5). Per side s, knob > 0 builds ONLY that side's cells:
- GAMMA_s cell + GAMMA_c_s port → RG_E_s (descend_to_rg_e) + RG_F_s
  (descend_to_rg_f) — rhythm-level; runner feeds
  `u[GAMMA_c_s] = drive * G["gamma_s"]`, so side s's RG sees
  drive·(1+gamma_s) through the tuned conductances.
- ALPHA_s cell + ALPHA_c_s port → EVERY ipsilateral MN pool at
  alpha_to_mn — tonic MN-pool drive; runner feeds
  `u[ALPHA_c_s] = drive * G["alpha_s"]`.
Implemented in build_network.py (parent; s3k/human default plant),
build_network_w2lvar.py, build_network_syn6.py (both inherit the
ALPHA→MN fan-out from the parent `_wire_balance`), params.py (knobs),
runner.py (the two-line port feed; port-name-guarded so absent ports
leave `u` untouched). Files: +133 insertions total, no existing line
changed semantically (git diff --stat above).

Gates (all RUN on egg post-change):
- defaults gate (410, 376, 1186) PASS;
- w2lvar (888, 382, 7866) / syn6 (794, 382, 3033) unchanged;
- anchor A (full_rules=0) BIT-EXACT -161.56754173676563;
- anchor B (full_rules=1) bit-exact -25.0 (its own pre-change value);
- winner+knobs topology probe (882, 386, 8492) with exactly the drawn
  wiring (lrdrive_probe_20261016.txt).

JSON RULE recipe for phase 4: add the four keys to the study's KEYS +
a `P.G[k] = float(p.get(k, 0.0))` loader branch (set_stage-style); the
runner needs nothing further (it reads G at feed time). Wean ladder =
scale gamma_l/r down toward 0 in the winner json.

## Not run / not verified this round
- No full 22 s walk with knobs ON (behavioral effect of the drive is
  phase 4's bootstrap sweep, out of this stage's scope; only build
  topology + bit-identity were gated here).
- w2lvar's own LBIN stacking was not re-probed (N1's w2lvar note is
  unconfirmed).
- build_network_spiking.py (BPA plants) not audited/modified — outside
  this ask's human-body scope.
- Ben's exports were parsed on the laptop AND re-dumped by the egg
  audit driver (identical bytes, SHA1-verified); no JSON was edited.
