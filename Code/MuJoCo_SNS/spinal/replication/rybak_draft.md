# Rybak replication draft — FIRST PASS (ZCode 2026-09-23)

Papers covered (the "Rybak two-level RG/PF rhythm model" family):

1. **Rybak 2006a** — Rybak, Shevtsova, Lafreniere-Roula & McCrea, "Modelling
   spinal circuitry involved in locomotor pattern generation: insights from
   deletions during fictive locomotion", J Physiol 577(Pt 2):617–639,
   DOI 10.1113/jphysiol.2006.118703, PMC1890439. Full text downloaded THIS
   session to `spinal\lit_rybak2006_fulltext.txt` (open access).
2. **Rybak 2006b** — Rybak, Stecina, Shevtsova & McCrea, "...insights from
   the effects of afferent stimulation", J Physiol 577(Pt 2):641–658,
   PMID 17008375. **NOT open access** (EuropePMC `isOpenAccess: N`, no PMCID)
   — no full text obtained; afferent-stimulation specifics below come from
   2006a's framework + the 2015 review + Rybak 2024 eLife instead. Stated
   plainly: this paper was NOT read. NOTE 2026-09-23 (supervisor block):
   an earlier fetch attempt with a guessed PMCID (PMC1890436) returned an
   UNRELATED article (Peng et al., HIF-1α/carotid body, J Physiol
   577:705–716); the stray file `lit_rybak2006b_fulltext.txt` it produced
   was DELETED — do not trust any copy you may have seen; no Rybak-2006b
   text exists locally.
3. **Rybak 2015** — Rybak, Dougherty & Shevtsova, "Organization of the
   Mammalian Locomotor CPG: Review ...", eNeuro 2(5):e0069-15.2015. Local
   full text `spinal\_rybak2015_plain.txt`.
4. **Rybak 2024 eLife** (local cache named lit_rybak2025.txt; the file
   header reads "Rybak et al. eLife 2024;13:RP103504") — "Operation of
   spinal sensorimotor circuits controlling phase durations during tied-belt
   and split-belt locomotion after a lateral thoracic hemisection". Local
   full text `spinal\lit_rybak2025.txt`.

Ben edits before build. Machine-readable edge list: `rybak_rules.json`
(same nodes/synapses schema as `w2l_equivalent_draft.json`).

## 1. Neuron populations (Rybak 2006a, reduced unilateral hindlimb model, Fig 2A + Model description)

Each population = 20 Hodgkin-Huxley spiking neurons (2006a, "Modelling
neural populations"). Per side:

| population | class | role |
|---|---|---|
| RG-E, RG-F | excitatory, INaP pacemaker | rhythm generation |
| Inrg-E, Inrg-F | inhibitory | RG reciprocal inhibition (laminated) |
| PF-E, PF-F | excitatory, INaP | pattern formation |
| Inpf-E, Inpf-F | inhibitory | PF reciprocal inhibition (laminated) |
| Ia-E, Ia-F | inhibitory | Ia inhibitory interneurons ("third level") |
| R-E, R-F | inhibitory | Renshaw cells |
| Mn-E, Mn-F | motoneurons | two-compartment (Booth et al. 1997 style) |
| MLR | external drive | tonic supraspinal excitation |

## 2. Full connection table (Rybak 2006a **Table 2**, exact values)

Weights are dimensionless per-spike increments (w_ji multiply the synaptic
conductance step; drives multiply d_mi) — NOT µS. Citation for every row:
Rybak 2006a Table 2 (local `lit_rybak2006_fulltext.txt`, "Table 2. Weights
of synaptic connections in the network").

| # | source | target | sign | weight | note |
|---|---|---|---|---|---|
| 1 | MLR | RG-E | exc | 1 (d_rg-e) | d_rg-e=0.5 (E-dominant) or 0.45 (F-dominant), Fig 2B/C legends |
| 2 | MLR | RG-F | exc | 1 (d_rg-f) | d_rg-f=0.43 / 0.51 |
| 3 | RG-E | RG-E | exc | 0.0125 | SELF-recurrent excitation |
| 4 | RG-E | RG-F | exc | 0.0125 | RG↔RG mutual EXCITATION |
| 5 | RG-F | RG-E | exc | 0.0125 | mirror of #4 |
| 6 | RG-F | RG-F | exc | 0.0125 | mirror of #3 |
| 7 | RG-F | Inrg-E | exc | 0.45 | flexor drives the IN that inhibits RG-E |
| 8 | Inrg-E | RG-E | inh | −0.115 | |
| 9 | RG-E | Inrg-F | exc | 0.45 | |
| 10 | Inrg-F | RG-F | inh | −0.115 | |
| 11 | MLR | PF-E | exc | 1 (d_pf-e = 0.5) | tonic drive to PF EXISTS in 2006 |
| 12 | RG-E | PF-E | exc | 0.0075 | WEAK homonymous RG→PF |
| 13 | Inrg-E | PF-E | inh | −0.05 | SAME Inrg also inhibits PF (opposite-phase gate) |
| 14 | Inpf-E | PF-E | inh | −0.35 | PF reciprocal inhibition |
| 15 | MLR | PF-F | exc | 1 (d_pf-f = 0.5) | |
| 16 | RG-F | PF-F | exc | 0.0075 | |
| 17 | Inrg-F | PF-F | inh | −0.05 | |
| 18 | Inpf-F | PF-F | inh | −0.35 | |
| 19 | PF-F | Inpf-E | exc | 0.2 | |
| 20 | PF-E | Inpf-F | exc | 0.2 | |
| 21 | PF-E | Ia-E | exc | 0.4 | **PF drives the Ia IN** (phase-gated reciprocal inh) |
| 22 | Ia-F | Ia-E | inh | −0.1 | **IaIN↔IaIN mutual inhibition** |
| 23 | R-E | Ia-E | inh | −0.1 | RC disinhibition of the Ia pathway |
| 24 | PF-F | Ia-F | exc | 0.4 | |
| 25 | Ia-E | Ia-F | inh | −0.1 | (Table prints "Ia-F (−0.1)" in the Ia-F row — read as Ia-E↔Ia-F mutual per the text; paper-typo flagged, verify against the PDF table before wiring) |
| 26 | R-F | Ia-F | inh | −0.1 | |
| 27 | Mn-E | R-E | exc | 0.25 | MN collateral → Renshaw |
| 28 | R-F | R-E | inh | −0.1 | **RC↔RC mutual inhibition** |
| 29 | Mn-F | R-F | exc | 0.25 | |
| 30 | R-E | R-F | inh | −0.1 | |
| 31 | PF-E | Mn-E | exc | 0.5 | |
| 32 | Ia-F | Mn-E | inh | −0.6 | reciprocal inhibition — strongest weight in table |
| 33 | R-E | Mn-E | inh | −0.2 | recurrent inhibition |
| 34 | PF-F | Mn-F | exc | 0.5 | |
| 35 | Ia-E | Mn-F | inh | −0.6 | |
| 36 | R-F | Mn-F | inh | −0.2 | |

NO afferent synapses in the reduced 2006a circuit — afferents/stimulation
enter as external drive terms d_mi with weights w_d_mi (2006a eqn (10)
description, local text); the dedicated afferent circuitry is in 2006b
(not obtained, see header).

## 3. Dynamics notes

- **Full HH spiking, not activity-based**: fast Na, K, CaN, CaL, Ca-dependent
  K, leak + INaP (2006a eqns (1)–(8), Table 1). MNs are two-compartment.
  Our SNS net is non-spiking activity-based — a 1:1 port is not intended;
  the wiring is the deliverable.
- **INaP kinetics (2006a Table 1)**: m∞ = (1+exp(−(V+47.1)/3.1))⁻¹
  instantaneous; h∞ = (1+exp((V+59)/8))⁻¹; τ_h = τ_h,max / cosh((V+59)/16).
  The τ_h,max value lives in the paper's appendix parameter table which my
  PMC text extraction did not capture — NOT verified here. The cosh bell has
  a large floor either way.
- **THE tau_h TRAP (project-known)**: Deng-style persistent-Na half-centers
  need tau_h FIXED ~350 ms; a tau_h(V) schedule that COLLAPSES at
  depolarized V quenches them (verified 2026-09-13, `spinal\_tau_h_check.py`
  + DESIGN.md). Rybak-lineage τ_h(V) bells keep 150–400 ms (Shevtsova 2026
  params) / 320–640 ms (Shinohara 2025 eq 6) — they never approach the
  ~0.1 ms collapse of sns_toolbox's schedule. Our build already runs the RG
  as `NonSpikingNeuronWithPersistentSodiumChannel` with FIXED
  TAU["rg_nap_h"]=0.35 s (`params.py:34`, `build_network.py:370-384`) —
  keep it fixed for any Rybak port.
- **Burst mechanism**: INaP-dependent pacemaker properties + RG↔RG mutual
  and self excitation (0.0125) sustain bursts; drive ASYMMETRY sets phase
  dominance (d_rg-e 0.5 vs d_rg-f 0.43 → T_E > T_F; swap → flexor-dominated;
  2006a Fig 2B/C).
- **The two-level diagnostic (why this paper matters to us)**: perturbations
  at the PF/MN level produce NON-resetting deletions; only RG-level inputs
  reset the rhythm (2006a; restated in Rybak 2015 Fig 1A and its caption —
  "perturbations and afferent signals acting at the PF or motoneuron level
  cannot reset the rhythm"). This is the literature basis for placing our
  phase-reset afferents at the RG (LIT_CIRCUIT_AUDIT.md §6 anchor quote).
- **Flexor dominance (2015 review)**: asymmetric flexor-driven RG
  organization is supported by deletion patterns (Zhong 2012) and phase
  durations (Talpalar 2013) — the review's §"asymmetric, flexor-dominated
  architectures". Shevtsova 2026 and Shinohara 2025 both allocate the
  rhythm to the FLEXOR center; our DRIVE biases E (DESIGN.md note).

## 4. Commissural/V-class frame (Rybak 2015 review, Figs 6-7)

- V0D = inhibitory CINs; V0V = excitatory CINs acting via V2a; V3 =
  excitatory CINs; V1 (IaIN/Renshaw family) + V2b = ipsilateral inhibitory;
  dI6 = third inhibitory CIN class ("not considered" in the 2015 models).
  (2015 review, class-definition sections + Fig 6 caption.)
- **Shevtsova 2015 Model 1 (2015 review Fig 7A, caption verbatim)**: left
  and right RGs interact via **CINe-F (V3)** = mutual EXCITATION between
  left and right flexor centers, **CINi-F (V0D)** = mutual inhibition, and
  **CINe1-F (V0V)**. Model 2 swaps the V0V path to CINe-E.
- **Rybak 2024 eLife sensory set (local lit_rybak2025.txt)**: SF-E1 =
  hip-flexor spindle afferents → EXC ipsilateral F half-center + INH
  contralateral F half-center (promotes E→F transition, grows with speed);
  SF-E2 = extensor group-Ib force-dependent EXCITATION of the ipsilateral
  E half-center (weight support); supraspinal drive PRESYNAPTICALLY
  INHIBITS all ipsilateral somatosensory feedback (α ≥ 0.35 → pattern mostly
  central). Regimes: state-machine (slow) → flexor-driven → classical
  half-center (local text, model description + results sections).

## 5. Mapping onto our SNS classes

| Rybak 2006a element | our element (build_network.py) | status |
|---|---|---|
| RG-E/RG-F NaP pacemakers | RG_E/RG_F, toolbox NaP class, fixed τ_h 350 ms | EXISTS |
| MLR drive | DRIVE input → RG-E 1.7 / RG-F 1.4 µS | EXISTS (different gain scale) |
| Inrg-E/F laminated RG inhibition | RG-E→InE (4)→RG-F (4) + mirror, NO direct RG↔RG synapses (build_network.py:385-391) | EXISTS — ours is single-gain two-stage; 2006a uses 0.45 pre / −0.115 post |
| PF-E/PF-F + Inpf lamination | PF_E1/E2/F1/F2 + PF_IN_E/PF_IN_F (build_network.py:508-527) | EXISTS (4 phase windows vs their 2 centers) |
| RG→PF weak exc 0.0075 | rg_to_pf = 2.4 (strong; audit row #4 "PF-forcED" deviation) | EXISTS-DEVIANT |
| Ia-E/Ia-F (PF-driven IaIN) | IaIN population, PF→IaIN 0.5 gate + Ia→IaIN 0.6 → antagonist MN 0.4 + IaIN↔IaIN 0.5 (conditional `ia_in`; the mutual edge also needs `full_rules`, build_network.py:795-800) | EXISTS (default 0) |
| R-E/R-F Renshaw set incl. RC→IaIN, RC↔RC | RC_{act}: MN→RC 1.0, RC→MN 0.5, RC→IaIN 0.5, RC↔RC 0.5 (conditional `renshaw`) | EXISTS (default 0) |
| (none — no afferent INs in 2006a) | Ia/II/Ib per-muscle encoders, IBEXC/LBIN, HEEL/TOE | OURS BEYOND 2006a (Di Russo/Deng/W2L heritage) |

## 6. EXISTS vs NEW (cross-checked against the compiled net THIS session)

Ground truth: `_net_edges.py` run in the myo env (exit 0; 82 edges at
default gains) + a variant probe with every conditional gain ON (exit 0;
**170 edges, 80 neurons** — renshaw 0.5, ia_in 0.62, heel 0.91, toe 0.63,
ib_rge 0.68, v3 0.3/0.3, ib_e_central 0.5, aff_* 0.4/0.4/0.3/0.3,
**full_rules 1.0** — the 2026-09-23 first probe omitted full_rules and
undercounted; corrected the same day after supervisor review).
Command: `C:\Users\Ben Bolen\.conda\envs\myo\python.exe _net_edges.py` from
`Code\MuJoCo_SNS\spinal`, plus `%TEMP%\net_edges_allon_20260923.py`.

**Already in our compiled net** (edge class present): DRIVE→RG-E/F;
RG-E→InE→RG-F lamination; RG→PF; PF_IN lamination; PF→MN; PF→IaIN;
Ia→IaIN→antagonist MN; **IaIN↔antagonist-IaIN mutual inhibition (g 0.5 —
build_network.py:795-800, "rule 1 addendum ... Deng A6", conditional on
ia_in>0 AND full_rules>0; verified in the full_rules probe dump:
`INH IaIN_IN -> IaIN_IN x2 g=0.5` in the 4-muscle representative build —
supervisor's full 92-muscle build measures 218 such edges)**; MN→RC→MN,
RC→IaIN, RC↔RC; KINH; Ib reversal (IBEXC) + Ib→IBIN→MN disynaptic
autogenic inhibition and **IBIN↔IBIN mutual inhibition (g 0.5, same
full_rules block)** + II→IIX→MN / II→IIIN→antagonist-MN (Di Russo rules
2a/2b); commissural CIN_F (RG-F→CIN_F→contra RG-F) and CIN_E
(RG-E→CIN_E→contra InE) — the latter is Shevtsova's V3-E→InE1 chain.

**NEW edges Rybak 2006a would add (not in our net at ANY gain):**
1. **RG-E↔RG-F mutual + self EXCITATION** (0.0125 ×3) — Table 2 rows 3-6.
   Ours has NO recurrent excitation anywhere in the RG. Smallest-change
   candidate for RG robustness at low drive; untested here — inference about
   its effect: NOT IN PAPER - inference.
2. **Opposite-RG inhibition of PF via the shared Inrg** (Inrg-E→PF-E −0.05;
   Table 2 row 13/17): our PF cells receive own-RG excitation and PF_IN
   cross-inhibition but NO direct opposite-RG-derived inhibition. Their
   mechanism closes each PF gate while the opposite RG bursts.
3. ~~IaIN↔IaIN mutual inhibition~~ **RECLASSIFIED 2026-09-23 (supervisor
   block): EXISTS-conditional** — already implemented at g 0.5 under
   ia_in>0 + full_rules>0 (build_network.py:795-800, Deng A6 rule-1
   addendum). Rybak 2006a Table 2 rows 22/25 (−0.1, E↔F antagonist pairs)
   is a SECOND literature source and a weight precedent (0.1 vs our 0.5)
   — merge the provenance, do not build a duplicate edge.
4. **Tonic DRIVE→PF** (MLR→PF, w=1 × d_pf=0.5; Table 2 rows 11/15): our
   DRIVE→PF 0.05 was REMOVED 2026-09-16 per audit row #23 ("not in
   Deng/Shevtsova"). Rybak 2006a Table 2 DOES include it — audit row #23's
   literature basis is incomplete. Correction note only; no code change
   without Ben's go.

## 7. Which Rybak circuits address our known gaps

- **Crossed flexor-excitatory edge (frozen left leg)**: 2015 review Fig 7A
  Model 1 CINe-F (V3) = l-RG-F ↔ r-RG-F MUTUAL EXCITATION; Shevtsova 2026
  Table 1 V3-F→c-RG-F (+0.03); Rybak 2024 SF-E1's contralateral-F
  INHIBITION is the same circuit with the opposite sign intent. Our net has
  every crossed path on the E side (CIN_E) or inhibitory (CIN_F); a
  crossed F→F excitation is NEW. Note the working W2L walker already used
  contact→contra RG_F excitation G=6 (`w2l_equivalent_draft.json`) — the
  behavior exists in our AnimatLab lineage, not in the MuJoCo net.
- **Sensory phase reset into the RG**: Rybak 2015 Fig 1 (reset ONLY at RG
  level) + Rybak 2024 SF-E1/SF-E2 give the canonical placement. Our current
  net: HEEL/TOE/LBIN→RG-E (extensor side, stance) + AFF_E/AFF_F→RG-E/RG-F
  semi-closed loops (v11b, `build_network.py:490-491`, default 0). The
  FLEXOR-side swing trigger (SF-E1 analog) exists only as the default-off
  AFF_F port with a lumped runner-side signal — no per-muscle fan-in, no
  length threshold. The per-MUSCLE flexor-length→RG-F wiring (Shinohara
  eq 10) is the concrete spec — see shinohara_draft.md.
- **IaIN↔IaIN and Ib-IN mutual inhibition**: IaIN↔IaIN = Table 2 rows
  22/25 (−0.1, E↔F) — second source for an edge we ALREADY have
  (g 0.5 under ia_in+full_rules, Deng A6-cited; see §6 item 3). Ib-IN
  mutual inhibition is NOT in Rybak 2006a (no afferent INs modeled); its
  literature source is Di Russo 2023 rule 3 (Ib INs "reciprocally inhibit
  the antagonist Ib-IN", LIT_CIRCUIT_AUDIT.md §10.2) — and it is ALSO
  already in our net under full_rules (IBIN↔IBIN g 0.5, verified in the
  corrected probe dump). Cross-reference only; nothing to build for either.

## 8. What Ben must decide

1. IaIN↔IaIN: CONFIRM/MERGE the gain provenance — the edge class already
   ships at g 0.5 under ia_in>0 + full_rules>0 (build_network.py:795-800,
   Deng A6); Rybak 2006a Table 2 adds a second source at −0.1. Decide
   whether the 0.5 stays or takes the Rybak/Deng average; do NOT build a
   duplicate edge. (Supervisor block 2026-09-23 — reworded from "wire it".)
2. Trial RG-E↔RG-F/self recurrent excitation (0.0125-scale, conditional
   default 0 — must respect the byte-identity-at-zero contract).
3. Opposite-RG→PF gate via the existing InE/InF (one new edge class).
4. Whether to reinstate DRIVE→PF now that 2006a supports it (audit #23
   amendment).
5. Whether the crossed F↔F excitation (V3-F/CINe-F) goes in as a CIN-class
   edge or via contact (W2L style).

## Verified / not verified

- Verified (ran/read): 2006a full text downloaded + Table 1/2 values read
  from it; 2015 review local text (Fig 1/6/7 captions + class sections);
  Rybak 2024 local text (SF-E1/SF-E2 description); `_net_edges.py` +
  all-conditionals probe incl. full_rules (exit 0 both; 170-edge dump
  shows IaIN↔IaIN x2 g=0.5 and IBIN↔IBIN x2 g=0.5 in the representative
  build).
- Not verified: Rybak 2006b (paywalled — not read); Rybak 2006a appendix
  numeric parameter table (τ_h,max etc. — extraction lost it; check the PDF
  before wiring dynamics); the 2006a Table-2 Ia-F row typo reading (verify
  against the rendered table); any simulation behavior claims (nothing was
  built or run for this draft). The full-92-muscle IaIN↔IaIN count (218
  edges) is the SUPERVISOR's measurement, not mine — my probe verified the
  edge class in the 4-muscle representative build only.
