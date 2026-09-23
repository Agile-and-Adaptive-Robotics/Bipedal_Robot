# Shinohara replication draft — FULL EXTRACTION 2026-09-23 (ZCode)

Source (RESOLVED — the old skeleton said "locate exact paper"): Shinohara,
Ambe, Kim, Mano, Tata Ramalingasetty, Lockhart, Markin, Ausborn, Rybak,
Danner & Aoi, "Mechanisms of adaptive interlimb coordination to sudden
ground loss: a neuromusculoskeletal modeling study", bioRxiv
2025.11.11.687930 (posted November 12, 2025, CC-BY). Full text held
locally in TWO copies with identical length (91,491 chars):
`spinal\shinohara_2025_biophiv687930_fulltext.txt` and
`spinal\lit_shinohara2025.txt`. Everything below was read from that text
THIS session. A cat hindlimb neuromusculoskeletal model — same
Rybak/Danner CPG lineage as ours, PLUS per-muscle afferent feedback INTO
the CPG centers. This is the paper that directly specs two of our named
gaps (sensory phase reset into the RG; crossed support for the partner
leg).

FIRST PASS by ZCode — Ben edits before build. Machine-readable edge list:
`shinohara_rules.json`.

## 1. Neuron populations (Fig 1 + model description)

Per side (left/right CPG), activity-based (non-spiking):
- **RG-F, RG-E** — rhythm generator half-centers; both have persistent Na;
  when UNCOPLED, RG-E sits TONIC under supraspinal drive and the rhythm
  is defined by RG-F ("Rhythmic oscillation of the RG is defined by the
  RG-F, which provides rhythmic inhibition of the RG-E through the IN-F").
- **IN-F, IN-E** — inhibitory interneurons of the RG mutual inhibition.
- **CINs: C1, V3** — commissural. **C1**: ipsilateral RG-F → inhibits
  CONTRAlateral RG-F. **V3**: ipsilateral RG-E → (crossed) excites
  contralateral IN-E → inhibits contralateral RG-F (verbatim: "the
  population of V3 neurons, which mediate the inhibition of the
  contralateral RG-F by the ipsilateral RG-E through the contralateral
  IN-E").
- **PF-F, PF-E** — pattern formation centers (NaP-bearing but
  sub-rhythmogenic; ḡNaP = 0.5 nS vs RG 4.5 nS — Table A.2).
- **MN-m ×7**: IP (hip flexor), GM (hip ext), VL (knee ext), TA (ankle
  flexor), SO (ankle ext), BF (hip ext/knee flex biart), GA (knee flex/
  ankle ext biart). Flexors {IP, TA, BF} ← PF-F; extensors {GM, VL, SO,
  GA} ← PF-E.
- Musculoskeletal side: 2D cat hindlimbs on a treadmill (forelimbs fixed
  to a stationary platform), 7 Hill-type muscles/leg, constant moment
  arms, RK4 @ 0.04 ms; treadmill + ground contact via viscoelastic
  elements; muscle activation through a low-pass filter.

## 2. Full connection table (Table A.1 + Table B.3 + eqs 10-11)

Central weights are dimensionless α/β multiplying gSynE/gSynI = 10 nS —
NOT µS. Citation per edge = Table A.1/B.3 unless noted.

| # | source | target | sign | gain | note |
|---|---|---|---|---|---|
| 1 | supraspinal drive (d=1.0, γ) | RG-F | exc | γ_F = 0.02 | flexor-defined rhythm |
| 2 | supraspinal drive | RG-E | exc | γ_E = 0.15 | E center tonic when uncoupled |
| 3 | RG-F | IN-F | exc | 0.4 | |
| 4 | RG-F | C1 | exc | 1 | |
| 5 | RG-F | PF-F | exc | 0.7 | |
| 6 | RG-E | IN-E | exc | 0.4 | |
| 7 | RG-E | V3 | exc | 0.5 | |
| 8 | RG-E | PF-E | exc | 0.7 | **EXTRACTION AMBIGUITY**: the plain-text table renders this row "- - - 0.4 - 0.5 0.7 -" which column-aligns the 0.7 under PF-F, but the text ("PF-F and PF-E ... receive the excitatory input from the RG neurons") and the E/F symmetry require RG-E→PF-E. Default reading = RG-E→PF-E 0.7; VERIFY against the PDF table before wiring |
| 9 | V3 (contralateral) | IN-E | exc | 0.15 | the crossed extensor-support edge |
| 10 | IN-F | RG-E | inh | 0.7 | |
| 11 | IN-F | PF-E | inh | 2.1 | strong PF gating |
| 12 | IN-E | RG-F | inh | 0.6 | |
| 13 | IN-E | PF-F | inh | 0.3 | |
| 14 | C1 (contralateral) | RG-F | inh | 1 | crossed flexor inhibition |
| 15 | PF-F | MN-IP / MN-TA / MN-BF | exc | 0.09 / 0.10 / 0.08 | Table B.3, CMA-ES optimized |
| 16 | PF-E | MN-GM / MN-VL / MN-SO / MN-GA | exc | 0.07 / 0.12 / 0.08 / 0.12 | Table B.3 |

**Afferent feedback (eqs 9-11) — ALL of it enters the EXCITATORY synaptic
current (s_j sits inside the I_SynE bracket of eq 7):**

| # | source signal | target(s) | gain | encoding |
|---|---|---|---|---|
| 17 | flexor muscle velocity (IP/TA/BF) | RG-F, IN-F, PF-F (all three) | k^v_F = 0.0007 | (v/lmax)^0.6, lengthening-only (positive v only) |
| 18 | flexor muscle length (IP/TA/BF) | RG-F, IN-F, PF-F (all three) | k^l_F = 0.53 | (l − 0.9·lmax)/lmax, only beyond 90% of max length |
| 19 | flexor muscle velocity+length | OWN MN (autogenic) | ×τv=1.0 / ×τl=1.1 on the same gains | eq 10 |
| 20 | extensor muscle force (GM/VL/SO/GA) | RG-E, IN-E, PF-E (all three) | k^f_E = 0.16 | F/Fmax, positive-only |
| 21 | extensor muscle force | OWN MN (autogenic) | ×τf=1.1 → 0.176 | eq 11: POSITIVE force feedback, UNGATED |

No IaIN, no Renshaw, no Ib autogenic inhibition, no cutaneous sensors —
their motor level is deliberately minimal; the afferent action is entirely
RG/PF/MN excitation. (Their discussion cites cat Ia/II flexor stretch and
Ib extensor unload physiology as the justification.)

## 3. Dynamics notes (eqs 3-8 + Table A.2, exact)

- C = 20 pF; ENa = 55 mV; ESynE = −10 mV; gSynE = gSynI = 10 nS;
  Vth = −50 mV; Vmax = 0 mV; drive d = 1.0.
- Table A.2: gLeak = 4.5 (RG) / 2.8 (IN, CIN) / 1.6 (PF, MN) nS;
  ḡNaP = 4.5 (RG) / 0.5 (PF) / 0.3 (MN) nS; ESynI = −75 (IN, CIN) /
  −70 (PF, MN) mV; ELeak = −62.5 / −60 / −64 / −64 mV.
- NaP kinetics: m∞ = 1/(1+exp(−(V+40)/6)) instantaneous; h∞ =
  1/(1+exp((V+45)/4)); **τh(V) = 320 + 320/cosh((V+35)/15) ms** → range
  320–640 ms. Same trap lesson as always: this bell's FLOOR (320 ms) is
  what makes NaP bursting survive — the sns_toolbox schedule that
  collapses toward ~0.1 ms quenches it. Our fixed τh = 350 ms sits just
  above their floor; keep fixed (`params.py:34`).
- PF and MN carry NaP but are NOT intrinsically rhythmogenic (low ḡNaP);
  they oscillate only when driven by the RG — a different PF philosophy
  from our phase-window cells.

## 4. The ground-loss result we want to replicate (their Figs 5/7/8 + Discussion)

1. Foot enters a hole → the leg hyperextends → flexors stretch beyond
   0.9·lmax → length feedback (k^l_F = 0.53) fires into RG-F/IN-F/PF-F →
   flexor burst starts EARLY → foot is pulled out ("increasing the length
   feedback from the flexor muscles caused the flexor activity to begin
   just after the foot entered the hole", Discussion).
2. Contralateral leg: extensor force feedback (k^f_E = 0.16) into RG-E/
   IN-E/PF-E maintains extensor activity (weight support), reinforced by
   the crossed V3 pathway from the ipsilateral (holed) RG-E → contra
   IN-E → suppression of contra RG-F ("along with the reciprocal
   contribution from the ipsilateral flexor center (Fig. 8) maintained
   the contralateral extensor activity").
3. WITHOUT these afferent loops the model fails to continue walking after
   ground loss (their Appendix C, Movies 6/7). On flat ground all three
   feedbacks matter little (their Appendix C Fig C.9) — the loops are
   specifically PERTURBATION circuitry.
4. Nullcline analysis (Fig 7, Appendix D): the adaptation is an afferent-
   driven change in NETWORK DYNAMICS (nullcline geometry), not just added
   input amplitude — their stated mechanism-level claim.
5. Parameter determination: 13 parameters by CMA-ES (7 PF→MN weights +
   6 afferent gains/multipliers) on randomized ±2 cm belt levels; cost
   favors long walking with low muscle activity + falling penalty (cost
   formula's exact terms not extracted — not verified beyond this).

## 5. Mapping onto our SNS classes

| Shinohara element | our element | status |
|---|---|---|
| RG-F/RG-E NaP, flexor-defined rhythm | RG_E/RG_F toolbox NaP fixed τh | EXISTS (drive allocation differs: ours biases E) |
| IN-F/IN-E lamination | RG-F→InF→RG-E etc. (g 4.0) | EXISTS |
| C1 (RG-F→c-RG-F inh, 1) | RG-F→CIN_F (4.0)→contra RG-F (4.0) | EXISTS |
| V3 (RG-E→V3 0.5→c-IN-E 0.15→c-RG-F inh 0.6) | RG-E→CIN_E (1.2)→contra InE (1.2)→contra RG-F (4.0) | EXISTS-STRUCTURAL (our crossed V3 lands on contra InE — same disynaptic shape; gains differ; conditional v3_gain, default 0) |
| PF-F/PF-E per side | PF_E1/E2/F1/F2 phase windows (+ joint_pf variant) | EXISTS-DEVIANT (4 windows vs 2 centers; ours lacks the IN→PF cross inhibition at 2.1/0.3 strength) |
| PF→MN weights (Table B.3) | W_PF_MN fitted table + pf_gain | EXISTS |
| flexor v/l → RG-F, IN-F, PF-F | AFF_F→RG-F (0.4) + AFF_F→PF-F (0.3) semi-closed loops, v11b | EXISTS-CLASS, default 0 (`build_network.py:490-491,545-553`) — but ours fans from ONE lumped runner-side signal, not per-muscle, no 0.9·lmax threshold, no v^0.6 compression, and does NOT hit the IN |
| extensor force → RG-E, IN-E, PF-E | Ib-aff→RG-E (0.5), →InE (0.5), →PF (0.5) via ib_e_central + LBIN→RG-E (0.68) | EXISTS-CLASS (conditional; encoding differs: our per-muscle Ib encoders vs their summed k^f_E) |
| extensor force → own MN (positive, ungated) | Ib→IBEXC (0.5)→MN (0.6), RG-E-gated | EXISTS-DEVIANT (ours is stance-gated + via IBEXC; theirs ungated direct) |
| (none) | IaIN, RC, KINH, HEEL/TOE, II afferents | OURS BEYOND the paper |

## 6. EXISTS vs NEW (cross-checked against the compiled net THIS session)

Ground truth: `_net_edges.py` default (82 edges, exit 0) + all-conditionals
probe (150 edges, exit 0). Commands in rybak_draft.md §6.

**NEW (not in our net at any gain):**
1. **Per-muscle afferent fan-in to the centers**: their k^l_F/k^v_F/k^f_E
   sum SEVEN per-muscle signals into RG/IN/PF of the matching phase. Our
   AFF_E/AFF_F are single lumped input ports driven by runner-side
   formulas. The per-muscle encoder fan-in (we already have 92 per-muscle
   Ia/II/Ib encoders!) is the missing wiring: flexor-group encoders →
   RG-F/IN-F/PF-F, extensor-group force encoders → RG-E/IN-E/PF-E.
   NOT IN PAPER - inference (beyond the paper): our version could reuse
   the existing per-muscle encoders rather than new receptors.
2. **Flexor afferent → IN-F** (their eq 10 includes j ∈ {IN-F}): our AFF_F
   does not target InF. NEW edge.
3. **Length threshold 0.9·lmax + v^0.6 compression + positive-only
   rectification** as the afferent ENCODING — our encoders use (L−Lmid)/
   Lhalf linear, L̇/0.6, F/Fmax. An encoding change, not an edge.
4. **Ungated autogenic extensor force→MN excitation** — ours exists only
   stance-gated via IBEXC.

**EXISTS already** (conditional, some default 0): AFF_E/AFF_F→RG-E/RG-F
and →PF; Ib→RG-E/InE/PF (ib_e_central); LBIN→RG-E; C1/CIN_F; CIN_E→contra
InE (their V3 chain); PF→MN tables.

## 7. Which Shinohara circuits address our known gaps

- **Frozen left leg / crossed support**: their answer is TWO-SIDED —
  (a) the affected leg's own swing is triggered by its flexor LENGTH
  feedback into RG-F (edge #18), and (b) the partner leg's stance is
  protected by crossed V3→contra-IN-E→contra-RG-F suppression plus the
  partner's own extensor force feedback (edges #9, #20). We have (b)'s
  crossing (CIN_E, default 0) and the lumped (a) port; we lack per-muscle
  (a). The paper does NOT contain a crossed flexor-EXCITATORY edge — that
  edge comes from Shevtsova 2026 (V3-F→c-RG-F 0.03) / Rybak 2015 Fig 7A
  (CINe-F). Corrected against the earlier skeleton's assumption.
- **Sensory phase reset into the RG**: this paper is the most concrete
  spec — flexor length/velocity → RG-F excitation terminates stance
  early (exactly our duty/frozen-leg failure mode); extensor force → RG-E
  prolongs stance. Placement at the RG matches the Rybak 2015 reset-level
  rule (audit §6). Our HEEL/TOE/LBIN cover the extensor side; the
  flexor-side per-muscle fan-in is the NEW piece.
- **IaIN↔IaIN / Ib-IN mutual inhibition**: NOT in this paper (no IaIN,
  no Ib INs). Both edge classes ALREADY exist in our compiled net under
  the `full_rules` conditional (IaIN↔IaIN g 0.5, build_network.py:795-800
  Deng A6 rule-1 addendum; IBIN↔IBIN g 0.5 in the same 2026-09-21 block —
  both verified in the corrected 170-edge probe dump). Literature sources
  remain Rybak 2006a Table 2 (−0.1) and Di Russo rule 3. Nothing to build;
  provenance merge only (see rybak_draft.md §6/§8, supervisor block
  2026-09-23).

## 8. Replication plan (Ben edits before build)

1. Unit mapping: their α/β dimensionless × gSyn 10 nS vs our µS regime —
   scale by behavior (same as the Deng/Nourse port; SNS_Units doc).
2. Biped adaptation: their cat 7-muscle/leg list → our 92-muscle gait2392
   groups (flexor group = hip/knee/ankle flexor pools; extensor force
   group = EXTENSOR_STANCE_GROUPS in muscle_map.py — the mapping choice
   is Ben's).
3. Afferent ports: extend the v11b AFF_E/AFF_F machinery to per-muscle
   fan-in with their encoding (0.9·lmax threshold, v^0.6) as conditional
   code, default 0, byte-identity preserved.
4. MSK side: our converted gait2392 + hole/ground-loss scenario needs a
   runner-side "hole" environment (treadmill with a dropping belt) — new
   scenario code, not in the runner today.
5. SUCCESS CRITERION: their Fig 5/7 behavior — with afferents ON, walking
   continues after sudden ground loss under one foot (early flexor burst
   on the affected side, maintained extensor activity on the support
   side); with the loops OFF, gait collapses (their Movies 6/7).

## Verified / not verified

- Verified (read from the bioRxiv full text this session): population
  inventory + C1/V3 verbatim description; Table A.1/A.2/B.3 values;
  eqs 3-11 including the afferent target sets and encodings; τh(V)
  formula; the Discussion/Fig 5/7/8 mechanism statements; Appendix C/D
  existence claims; CMA-ES 13-parameter setup.
- Not verified: the exact cost-function terms (formula not fully
  extracted); the RG-E→PF-E column alignment (extraction ambiguity,
  flagged in §2 edge #8); figure-number references (Fig 5/7/8 cited from
  the Discussion text, not opened); no code was built or run for this
  draft.
