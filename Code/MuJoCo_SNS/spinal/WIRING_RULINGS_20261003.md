# Wiring rulings: s3k, syn6, w2lvar (2026-10-03)

Audit of the three walker networks against the rules we actually trust.
Per model: current wiring, rules checked against, conflicts found, choice
made and why, remaining uncertainty. Two provable wiring bugs were fixed
(gated on the counts gate, which passes before and after); everything
else is a documented ruling, not a code change. No gain value was
changed anywhere. No ZCODE-fenced block or adopted Section 4.4 value was
touched.

## Sources audited

- Ben's master rules drawing: [Circuit_rules_CONNECTOME_md__connectome.json](D:/Github/Bipedal_Robot/Neuromechanical_Models/Mujoco_SNS_models/Circuit_rules_CONNECTOME_md__connectome.json), byte-identical to the tracked copy [ben_rules_20260924.json](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/ben_rules_20260924.json) (verified with `fc /b`: IDENTICAL).
- Ben's Shevtsova and Shinohara exports: [ben_shevtsova_20260924.json](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/ben_shevtsova_20260924.json) and [ben_shinohara_20260924.json](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/ben_shinohara_20260924.json).
- Replication drafts in [replication](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/replication): `rybak_rules.json`, `shevtsova_rules.json`, `shinohara_rules.json` (each synapse carries a citation and an exists/new status).
- The Deng wiring audit: [deng_feedback_wiring_audit.md](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/reports_20260924/deng_feedback_wiring_audit.md).
- SADb knowledge base: [knowledge_base](D:/Github/Bipedal_Robot/SADb_audit/knowledge_base), pages read in full this round: `pathways/ia-monosynaptic-excitation.md`, `pathways/ia-reciprocal-inhibition.md`, `pathways/ib-excitatory.md`, `pathways/ib-inhibition.md`, `pathways/ib-stance-to-swing.md`, `pathways/type-ii-excitatory.md`, `pathways/ii-inhibitory.md`, `pathways/cutaneous-stance-modification.md`.
- Code: [build_network.py](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/build_network.py), [build_network_syn6.py](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/build_network_syn6.py), [build_network_w2lvar.py](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/build_network_w2lvar.py), [params.py](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/params.py), [synergy_model.py](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/synergy_model.py), [fsa_backsolve.py](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_backsolve.py), [runner.py](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/runner.py), plus the s3k winner record `reports_20260923/s3k_trial34_full_params.json`.

Method: I built each network (default gains, the s3k trial-34 winner gains, and the variant selectors) and dumped every compiled connection with sign and conductance, then compared edge-by-edge against the rule files. Probe scripts and their outputs live in [wiring_audit_20260903](D:/Github/Bipedal_Robot/tmp/wiring_audit_20261003).

## Gates (exact commands and outcomes)

The ask names `tmp/digest_20261002/spinal_counts_gate.py` with the myo
python. Run exactly as shipped, before any change:

```
C:\Users\Ben Bolen\.conda\envs\myo\python.exe tmp\digest_20261002\spinal_counts_gate.py
TypeError: build() missing 1 required positional argument: 'model_actuators'
```

The shipped script cannot pass: it calls `bn.build()` with no arguments
and unpacks `net, meta = ...`, while `build_network.build(model_actuators, dt, interleg)`
returns a single `SpinalNetwork`. I therefore ran the repo's own
equivalent defaults gate (same canonical counts, real API) before making
changes:

```
C:\Users\Ben Bolen\.conda\envs\myo\python.exe reports_20260925\audit\audit_defaults_gate.py
counts neurons/inputs/synapses = (410, 376, 1186)
GATE: PASS
```

I then repaired the digest gate script to use the real API (same check,
same canonical counts; the repair is documented in its header), and
applied the wiring fixes. After the fixes, both gates pass:

```
C:\Users\Ben Bolen\.conda\envs\myo\python.exe tmp\digest_20261002\spinal_counts_gate.py
counts: (410, 376, 1186)
COUNTS GATE PASS
C:\Users\Ben Bolen\.conda\envs\myo\python.exe reports_20260925\audit\audit_defaults_gate.py
counts neurons/inputs/synapses = (410, 376, 1186)
GATE: PASS
```

Note on the standing rule "deterministic gates are run by the script":
the ask explicitly instructs running this named check with my own tools
before and after any fix, so I ran it; the repair also means the
workflow's own later run of that file will work instead of erroring.

## 0. Changes applied this round

1. **s3k/w2lvar mutual-inhibition completion (the one provable wiring
   bug).** Under `full_rules > 0`, the IaIN-to-IaIN and IBIN-to-IBIN
   mutual inhibition specified by Deng Table A6 and Rybak 2006 Table 2
   was one-directional: the INs are created lazily inside the per-muscle
   wiring loop, so the earlier-wired muscle's antagonist loop ran before
   the later muscle's IN existed, and the `in self.idx` guards silently
   skipped the edge. Measured pre-fix (probe `probe_mutual_dir.py`):
   218 of 436 directed edges per family in the stock full_rules build,
   every present edge pointing later-to-earlier in wiring order; w2lvar
   identical (probe `probe_w2l_mutual.py`); syn6 unaffected (436/436, it
   pre-creates all four motif INs). The in-code comment "both directions
   emerge from the two muscles' wiring loops" is true for RC-to-RC (RCs
   are created in the neuron pass) and was false for IaIN/IBIN. Fix: a
   completion pass (`_complete_in_mutual`, build_network.py:928) that
   adds only the missing directed edges at the same 0.5 conductance. No
   existing edge, neuron, or index changes; the default build (both
   gains 0) never enters the pass, so the 410/376/1186 gate is
   structurally unaffected, and it passes. Post-fix measurements:
   436/436 in both families in both builders; s3k winner-topology build
   went from (786, 382, 3410) to (786, 382, 3846); w2lvar from
   (888, 382, 7430) to (888, 382, 7866); syn6 unchanged (794, 382,
   3033). Consequence, stated loudly per the no-silent-change rule: all
   recorded `full_rules > 0` eval scores (the s3k family and the w2lvar
   curriculum) are PRE-FIX values produced with half the mutual
   inhibition; a winner replay today will score differently. No gain
   number was altered.
2. **Digest gate repair** described above.

---

## 1. s3k (the stock [build_network.py](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/build_network.py) build)

### Current wiring (verified by building, not just reading)

- Default build (all conditional gains 0): 410 neurons, 376 inputs, 1186
  synapses, counts gate PASS.
- Production configuration is the stage-3/4 winner
  (`s3k_trial34_full_params.json`): full_rules 1.0, ia_in 0.825,
  heel_rge 0.375, toe_rge 0.146, ib_rge 0.674, ib_e_central 0.244,
  ia_f_central 0.399, ii_f_central 0.384, ii_e_central 0.093, c1_gain
  0.498, v3_gain 0.095, v3_to_ibexc 0.878, ia_f_contra_f 0.797,
  f1_anklepf_inh 0.929, contra_kinh 1.425, contact_onset 0.814. Built at
  those gains the network is (786, 382, 3846) after the fix; every edge
  below was dumped and checked (probe `probe_s3k_winner.py`).
- RG: persistent-Na half-centers with FIXED tau_h
  (`SNS_NumpyFixedTau`, params rg_nap_h 0.35), laminated mutual
  inhibition through InE/InF, weak mutual excitation conditional.
- PF: four phase cells PF_E1/E2/F1/F2 with PF_IN_E/F lamination
  (s3k production runs the phase-cell PF; joint_pf is 0 in the winner,
  matching the Deng audit note).
- Ia: monosynaptic Ia to MN (0.6), Ia to IaIN (0.6) to antagonist MN
  (0.4), PF_F1 phase gate onto IaIN (0.5), RC-to-IaIN disinhibition when
  Renshaw is on.
- II (full_rules): II to IIX (1.0) to MN (0.4); II to IIIN (1.0) to
  antagonist MNs (0.4); direct II-to-MN edge correctly absent.
- Ib (full_rules): Ib to IBIN (1.0) to MN (0.35, inhibitory);
  stance-gated reversal Ib to IBEXC (0.5, gated by RG_E at 1.0) to MN
  (0.6); IBEXC to LBIN (0.5), LBIN to RG_E (0.674).
- Heel/toe (full_rules branch, build_network.py:458): heel to ipsi InE
  EXC and heel to ipsi InF INH; toe to InE EXC; no direct heel/toe edges
  onto the half-centers.
- Central afferent knobs: extensor Ib to {PF_E1, PF_E2, RG_E, InE};
  flexor Ia and II to {PF_F1, PF_F2, RG_F, InF}; flexor Ia additionally
  inhibits the contralateral RG_F (ia_f_contra_f).
- Commissurals: RG_F to CIN_F (inhibits contra RG_F); RG_E to CIN_E
  (excites contra InE, plus contra IBEXCs via v3_to_ibexc).

### Rules checked against, and verdicts

- **F1 (pass). Ia motif.** Monosynaptic homonymous excitation plus
  reciprocal inhibition through an IaIN population with phase gating
  matches the SADb corpus (Talpalar 2011, Feldman and Orlovsky 1975,
  Capaday 1986 for phase gating, Pratt and Jordan 1987 for
  RC-to-IaIN), Deng A6, and Rybak 2006 Table 2 (PF-to-IaIN 0.4 exists;
  IaIN mutual now bidirectional after the fix).
- **F2 (pass). Ib stance reversal.** Autogenic inhibition at rest plus
  stance-gated group excitation is exactly the reversal doctrine of the
  SADb Ib pages (Pearson and Collins 1993; Procházka 1997; Grey 2007;
  Dietz and Duysens 2000; Donelan and Pearson 2004 loop gains), and
  LBIN-to-RG_E realizes the rhythm-layer access of Dominguez 2020,
  Conway 1987, and Gossard 1994.
- **F3 (fixed). IaIN/IBIN mutual inhibition was one-directional.** See
  change 1 above.
- **F4 (conflict, documented choice: keep). Direct afferent-to-HC
  central knobs.** The Deng wiring audit flags ia_f_central,
  ii_f_central, ii_e_central, ib_e_central as the clearest deviation
  from the original W2L (which had zero afferents onto half-centers).
  But Ben's own master-rules drawing contains SN-Ia and SN-II to
  HC-RG-F EXC 0.5 and Ib grp to HC-RG-E EXC 0.5 as direct edges, and
  Shinohara 2025 eq 10 and 11 target j in {RG, IN, PF} directly, while
  Shevtsova has no afferents at all. Choice: Ben's drawing and Shinohara
  (the two sources that speak to afferent targets at this level of
  detail) support the direct edges; the Deng-original doctrine is the
  minority position here and is already recorded as such in the Deng
  audit. Keep, gains stay as tuned.
- **F5 (conflict, documented choice: keep). Heel routing side and
  sign.** Ben's drawing: heel to IPSI InE EXC 0.5, heel to IPSI PF
  E-lamination INs EXC 0.5, and heel to CONTRA InF EXC 0.5 (edge heel to
  IN-InF_2813, a node in the drawing's contralateral RG block). The
  stock full_rules branch instead wires heel to ipsi InE EXC plus heel
  to ipsi InF INH (both push the ipsilateral side toward extension, a
  coherent mechanism, but neither edge matches the drawing's contra-InF
  excitation), and the per-PF-layer keys (heel_pf_layer, heel_in_f_exc)
  target ipsilateral cells only, so the drawing's contra edge is not
  expressible in the stock stack today. Choice: keep s3k as tuned (it
  predates the drawing and its numbers are records); the drawing reading
  is implemented and verified in w2lvar; if Ben wants the contra edge in
  the stock stack it needs a new gain key, which I did not add silently.
- **F6 (deviation, documented choice: keep). II relay gains.** II to
  IIIN uses 1.0 (the generic relay constant) where the master rules draw
  0.5, and IIIN-to-antagonist-MN reuses G["ia_to_antagonist"] (0.4)
  where the drawing has a dedicated 0.5 (syn6 implements both verbatim).
  This is knob reuse, not a transcription error; changing it would shift
  tuned behavior for no rule gain. Document and keep.
- **F7 (pass). DRIVE-to-PF removed.** Conflicts with Rybak 2006 (MLR to
  PF 0.5) but matches Deng and Shevtsova; the removal is recorded in
  `rybak_rules.json` itself as the resolved status. Keep.
- **F8 (pass). V3 target.** V3 excites the contralateral InE rather than
  the contralateral RG_E, matching Ben's rules drawing and Shinohara's
  crossed V3-to-IN-E row, and avoiding the bilateral E-to-E latch the
  code documents. v3_to_ibexc adds the crossed-extensor reinforcement
  (Rybak 2025 SF-E2 analog).

### Remaining uncertainty (s3k)

- Shevtsova's asymmetric lamination weights (IniF to RG_E −1 vs IniE to
  RG_F −0.1) are collapsed to one symmetric rg_mutual_inh; all paper
  weights are dimensionless and scaled to the stack's 0-to-5 mV regime,
  so only relative structure, not absolute conductance, is comparable.
- Renshaw cells are absent from the s3k production winner (no renshaw
  key in the trial-34 params, so G remains 0) although Deng A6 carries
  them; the s3k stack ran that way and I did not change it.

---

## 2. syn6 ([build_network_syn6.py](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/build_network_syn6.py))

### Synergy implementation verdict: conceptually right, verified end-to-end

Checked against [synergy_basis.npz](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/synergy_basis.npz),
[synergy_model.py](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/synergy_model.py), and the FSA back-solve
([fsa_backsolve.py](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_backsolve.py) plus
`fsa_results/fsa_backsolve.npz`). All numbers below are from probes run
this round, not from documentation.

- **Weight construction.** W_r and W_l are exactly [43 muscles x 6];
  `W == gain.T` verbatim from the fsa npz `{side}_gain` for both sides
  (np.array_equal), names identical to the fsa names, every name
  classifiable, and the runner's six pruned actuators absent from the
  basis. This is the layout `synergy_model.stage_basis` writes, nothing
  re-factorized.
- **Rank.** rank(W) = 6 on both sides (full rank). The fsa report
  confirms 6 as the smallest shared rank with centered VAF at or above
  0.90 (right 5, left 6), and the stored-H replay reproduces the
  recorded basis quality: R2 0.9261 (r) and 0.9126 (l).
- **Eq-18 mapping.** All 311 PF_S-to-MN edges equal
  `fsa_backsolve.analytical_conductance(W[act, k])` exactly (atol 1e-12,
  probe `probe_syn6_build.py`); the builder imports the audited function
  rather than reimplementing it, so the math cannot drift. Hand check on
  five sample weights matches (for example k = 1.008493 gives
  g = 1.704957 by hand and by function). Validity bound: max W is 1.167
  (r) and 1.235 (l), below dE/R = 8/5 = 1.6, so zero entries are
  Eq-18-invalid, matching the module docstring claim.
- **Per-synergy MN mapping.** Edges exist only where W > 0; pruned and
  out-of-basis muscles receive no PF drive while keeping their MN and
  afferent arc, mirroring the runner's force-zeroing of the same list
  (runner.PRUNE_MUSCLES equals the builder's set, verified by import).
  This is exactly the NOTE point-1 design in synergy_model.py: each basis
  column becomes one PF-population-to-MN-pool projection.
- **Family and phase assignment.** The fsa npz carries
  `{side}_pf_phase_mean` and `{side}_phase_grid` (the builder's
  commented warning about the grid key is correct as written; no silent
  fallback fired). Symmetrized stance fractions give families F, E, E,
  F, F, E; the S5 bilateral instability claimed in the docstring is real
  (stance_frac S5 = 0.581 r vs 0.239 l) and the symmetrization is the
  documented resolution of synergy_model NOTE open item (b). Only S5
  falls inside the mixed band, and it receives the split RG drive
  0.984155 from RG_E and 1.415845 from RG_F, which is exactly
  rg_to_pf(2.4) times (0.410065, 0.589935), the unrounded symmetrized
  stance fraction; the committed channels are single-family at 2.4 and
  excluded-channel lamination behaves as documented. The dorsiflexion
  channel is S5 on both sides (largest ankle_df W mass among F-family
  channels), so TOEDF inhibits PF_S5.
- **Open items (a) and (c) of the synergy_model NOTE** are resolved by
  design in this builder: channels map 1:1 to PF populations (a), and
  phase-windowing is realized through RG-family membership plus measured
  peak-phase tau staggering rather than replayed H (c).

### Rules dress (verified by probe)

Heel and toe dress, commissurals, and the Shinohara autogenic motif all
match their stated sources at the stated gains: heel to InE 0.5, heel to
PF_IN_E 0.5, toe to TOEDF 5.0, TOEDF to PF_S5 at 2.749; V0V to contra
InE 0.6, V0D to contra RG_F 0.07, V3E to contra RG_E 0.02 and contra InE
1.0 with drive edges 1.0/0.7/0.35/1.0; per-muscle Ia to MN 2.0, Ia to
IaIN 1.0, IaIN to antagonist MN 0.5, II to IIX 1.0 to MN 0.5, II to
IIIN 0.5 to antagonist MN 0.5, Ib to IBIN 1.0 to MN 0.5; mutual
IaIN-to-IaIN and IBIN-to-IBIN at 0.5 with 436 of 436 directed edges
(both directions; this builder pre-creates its INs, so the s3k/w2lvar
bug never existed here). IBEXC wiring uses the rule gains 0.5/0.5 for
all five extensor stance groups.

### Conflicts found

- **S1 (conflict, documented choice: do not change). Heel-to-InF side.**
  syn6 wires heel to the IPSILATERAL InF (build_network_syn6.py:487).
  Ben's master-rules drawing has exactly one heel-to-InF-family edge,
  heel to IN-InF_2813, and that node belongs to the drawing's
  contralateral RG block (with its own lamination), so the drawing reads
  heel to CONTRA InF. w2lvar implements the contra reading
  (build_network_w2lvar.py:265); the AGENTS.md prose summary
  ("heel_IN to InE/InF + the three PF-layer INs") drops the side, which
  is where the ipsi reading came from. Mechanically the two ipsilateral
  edges partially cancel (InE-exc inhibits RG_F while InF-exc inhibits
  RG_E), whereas the drawing's ipsi-InE-exc plus contra-InF-exc do not
  cancel. Ruling: the drawing is the precise artifact and wins; the
  syn6 heel_to_inf edge should be contralateral. Not changed: syn6's
  dress is unconditional (no default-0 knob), the curr_syn6 studies were
  tuned on the current wiring, and re-routing would invalidate that
  tuning without a rerun. For the next syn6 retune, move heel_to_inf to
  the contralateral InF, or make both sides gain-gated.
- **S2 (deviation, documented choice: keep). IBEXC has no stance gate.**
  The stock and w2lvar IBEXC is gated by RG_E; syn6's is not. That
  matches Ben's drawing literally (the ib_rev motif is purely
  afferent-to-IN-to-MN) but drops the SADb's phase-dependent reversal
  gating (the Procházka 1997 switch from negative in posture to positive
  in locomotion); in syn6 the reversal is structural (extensor groups
  only). The module does not list this among deviations D1 to D3, so it
  is recorded here.
- **S3 (documented behavior, not a bug). LBIN load dress is inert in the
  recorded winner.** The load dress (Ib grp to RG-E, IN-InE, HC-PF-E,
  all 0.5) is gated on G["ib_rge"], and the stage-5 winner has
  ib_rge = 0.0, so those master-rules edges carried no signal in the
  recorded runs. The code documents why (unconditional force-proportional
  tonic Ib E-latches the RG in air, gate c2 measurement). Verified: at
  ib_rge = 0.5 the ten LBIN out-edges and 54 Ib-to-LBIN edges appear.

### Remaining uncertainty (syn6)

- **Analytic Eq-18 versus dynamic fit.** syn6 uses the analytic
  conductances; the fsa back-solve itself quantifies the
  overlapping-input correction: analytic-only centered VAF is 0.837 (r)
  and 0.795 (l) versus 0.863 and 0.835 for the dynamically fitted
  conductances (g_fit mean roughly half of g_analytic). The builder's D1
  rationale (the PF source itself reaches E_HI under S3K drive, fsa
  analytic peak 1.27) is the documented justification; the 3-to-4 VAF
  point gap is the honest cost. A g_fit-based dress is the available
  upgrade if Ben wants it.
- The tau staggering interpolates the stock PF_SHAPE constants over
  within-family peak-phase rank, so per-side assignments differ slightly
  (for example S6 gets 0.6 on the right and S3 gets 0.6 on the left);
  this follows the measured peak phases and is intended, not a bug.

---

## 3. w2lvar ([build_network_w2lvar.py](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/build_network_w2lvar.py))

### Current wiring (verified by probe `probe_w2lvar_build.py`)

Build at the VARIANT_G overlay: (888, 382, 7866) after the fix. Two PF
pairs per side (HIP and the merged KNEE+ANKLE cell) with per-pair IN
lamination at the drawing's 2.749 in both directions; RG drive at 2.4
(calibrated, D-noted). PF-to-MN weights are 2.0 times the fitted
[joint_pf_weights.json](D:/Github/Bipedal_Robot/Code/MuJoCo_SNS/spinal/joint_pf_weights.json)
rows summed over each cell's source rows (verified numerically for
vas_lat, soleus, glut_max1, iliacus, semimem, med_gas, tib_ant; for
example med_gas gets KNEE-E + ANK-E = 1.9997, semimem gets HIP-E 0.965
plus the biarticular KNEE-F cross 1.461), with the anatomical 1.0
fallback and trunk deliberately absent (documented: the W2L walker has
no trunk). Heel dress matches the drawing: ipsi InE EXC 0.5, CONTRA InF
EXC 0.5, and both PF E-lamination INs EXC 0.5. Toe dress: TOEDF EXC 5.0
then inhibition of the merged KNEE-F cell at 2.749 (the documented
side effect that knee-flex drive dips with dorsiflexion inhibition).
Shevtsova commissurals with the INI relay retained (V2a 4.0, V0V 4.0,
V0V-to-INI 2.4, INI-to-contra-RG_F 2.0, V0D 2.8/1.2, V3E 1.4/0.04/1.0,
DRIVE inhibits V0V and V0D at 2.0). Autogenic central projections
divide the drawing's aggregate node gains by the per-side family count
(Ib 0.5/27 = 0.018519 onto the four E targets; flexor Ia and II
0.53/19 = 0.027895 onto the four F targets), so the family sum
reproduces the one-node conductance; measured and verified. Renshaw
always built (MN-to-RC 1.0, RC-to-MN 0.5, RC-to-IaIN 0.5, 4140 RC
mutual edges = full same-side broadcast). IBEXC keeps the parent's
RG_E gate. KINH absent at defaults. The PF_ANK names are idx aliases
onto the merged KNEE cells for the runner's watch columns: zero edges
source from them, true neuron count 888 confirmed.

### Conflicts found

- **W1 (fixed). Mutual inhibition one-directional**, identical bug and
  fix as s3k (change 1 above); 436/436 after.
- **W2 (documented behavior). Zero-gain LBIN-to-RG_E edge always
  built** (VARIANT_G sets ib_rge 0.0, and the variant wires the edge
  unconditionally). It is inert but present in the counts, and the
  in-code comment says why it is kept (runner contract). The variant is
  not under the default bit-exactness contract, so this is accepted.
- **W3 (doc drift, no change). heel_rge does double duty**: the runner
  scales the HEEL_c port current by it (runner.py:1407) while the
  builder also uses it as the conductance of the ipsi-InE and contra-InF
  edges, so tuning it changes both the sensor gain and two synapses at
  once (effective squared scaling). The VARIANT_G comment ("port-current
  scale; edges are wired variant-side") reads as if the edges used
  variant-local gains, which they do not. Recorded here rather than
  edited, since a comment rewrite in a tuned variant file buys nothing.
- **W4 (documented deviation, keep). Commissural gains are calibrated,
  not literal.** The ×4.0 scale is the w2l_cpg calibration, and two
  crossed-F legs were raised for the measured antiphase lock (INI-to-RG_F
  0.075 to 2.0; V0D-to-RG_F 0.07 to 1.2); all four Shevtsova cell
  classes are kept, and the raises are annotated in GAINS.
- **W5 (open polarity note).** Exciting the contra InF suppresses the
  partner's RG_F (keeps the partner extensor). At ipsilateral heel
  strike the partner is at terminal stance, so the drawing's edge
  delays partner swing onset; the code comment calls it crossed stance
  support. The wiring follows the drawing either way; only the
  commentary's functional story is uncertain.

### Remaining uncertainty (w2lvar)

- The merged KNEE+ANKLE cell is a documented lumping of two drawing
  micro-layers (one E-IN instead of separate knee and ankle INs, and
  TOEDF hitting knee-flex drive with dorsiflexion drive); Ben's
  per-micro-layer lamination is only fully realized if the pairs are
  split again.
- Trunk muscles carry POSTURE and BAL drive but no PF drive, so trunk
  control in this variant is purely the balance layer.

---

## 4. Source-conflict ledger (choices in one place)

| # | Conflict | Sources | Choice |
|---|----------|---------|--------|
| 1 | Afferents onto RG/PF half-centers | Deng original: none; Ben master rules + Shinohara eq10/11: direct; Shevtsova: no afferents | Ben's drawing + Shinohara (kept, F4) |
| 2 | Heel-to-InF side and sign | Drawing: ipsi InE exc + contra InF exc; stock full_rules: ipsi InE exc + ipsi InF inh; AGENTS prose: unspecified side | Drawing wins; stock kept as tuned, contra edge realized in w2lvar, syn6 flagged S1 |
| 3 | II relay gains | Drawing: II-to-IIIN 0.5, IIIN-to-ant-MN 0.5; stock: 1.0 and reuse of ia_to_antagonist 0.4 | Keep stock calibration, document (F6); syn6 uses verbatim |
| 4 | DRIVE-to-PF | Rybak: MLR-to-PF 0.5; Deng/Shevtsova: RG only | Already removed 2026-09-16; keep (F7) |
| 5 | IBEXC stance gate | Drawing: no gate; stock/w2lvar: RG_E-gated; SADb: phase-dependent reversal | syn6 literal (S2); stock/w2lvar gated; both defensible, both documented |
| 6 | V3-E crossed gain asymmetry | Ben's file: 0.02 one direction, 0.5 the other; task text: small exc | 0.02 both directions (both variants, D2-documented) |
| 7 | Ini relay folding | Shevtsova: separate Ini (V0V-to-Ini-to-RG_F) | syn6 folds into InE (D2); w2lvar keeps a real INI cell |

## 5. What was not run or not verified

- No MuJoCo runner evaluation was run (no ground/air walk replay), so
  the behavioral effect of the mutual-inhibition fix on the recorded
  s3k and w2lvar winner scores is quantified only at the topology level
  (+436 directed edges each); replaying the winners is Ben's call.
- The SADb corpus pages beyond the nine listed above were not read this
  round; cluster pages and animal pages were not used.
- `build_network_spiking.py` (the 2026-10-02 goal-1 spiking mirror) is
  outside this ask's three models and was not audited.
- Ben's Shinohara drawing's duplicated-gain rows (Ia vs II subpopulation
  duplicates) were taken as the drawing states them; the variant's
  choice to reuse the stack's ia_to_mn and ii_to_mn calibrations instead
  is documented in the w2lvar docstring and was not re-litigated.

Probes: [probe_syn_basis.py](D:/Github/Bipedal_Robot/tmp/wiring_audit_20261003/probe_syn_basis.py), [probe_syn6_build.py](D:/Github/Bipedal_Robot/tmp/wiring_audit_20261003/probe_syn6_build.py), [probe_census.py](D:/Github/Bipedal_Robot/tmp/wiring_audit_20261003/probe_census.py), [probe_w2lvar_build.py](D:/Github/Bipedal_Robot/tmp/wiring_audit_20261003/probe_w2lvar_build.py), [probe_s3k_winner.py](D:/Github/Bipedal_Robot/tmp/wiring_audit_20261003/probe_s3k_winner.py), [probe_mutual_dir.py](D:/Github/Bipedal_Robot/tmp/wiring_audit_20261003/probe_mutual_dir.py), [probe_w2l_mutual.py](D:/Github/Bipedal_Robot/tmp/wiring_audit_20261003/probe_w2l_mutual.py). Run with
`C:\Users\Ben Bolen\.conda\envs\myo\python.exe` from
`D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal`.
