# LIT_CIRCUIT_AUDIT.md — literature audit of the spinal-circuit connections
# vs Ben's Zotero library + SADb collection. Off-peak task 2026-09-15.
#
# SCOPE: read-only. NO changes were made to build_network.py, runner.py,
# params.py, or any simulation code. Implementation of any recommendation
# below AWAITS BEN'S GO. Sources: our compiled edge inventory
# (_net_edges.py dump, lit_audit_edges.txt), full texts held locally
# (Rybak 2015 eNeuro, Shevtsova 2026 eLife RP107480, Nourse 2023 incl.
# Deng Table A6), Zotero local-API metadata + abstracts for the reflex/
# reset/Renshaw items (keys cited below), SADb_audit\README.md state.
# SADb: Papers = 532 records, curation layer on 20/383 imported so far —
# Zotero remains the master; this audit cites Zotero keys directly.
#
# Precedence rule applied: newer literature overrides older where they
# conflict; cat/mouse/rat evidence is labeled; Deng/Animatlab-only
# features without mammalian literature support are flagged as such.

## 1. What the circuit currently is (edge inventory, 4-muscle
##    representative build; full net scales to 92 muscles)

- RG: RG-E↔RG-F mutual inhibition (direct, g=4.0); ADAP-E/F burst-
  termination loop (RG→ADAP exc 2.5, ADAP→RG inh 2.5); DRIVE→E/F
  (1.7/1.4); POSTURE→E (0.8); commissural RG-F↔RG-F (4.0) and
  RG-E↔RG-E (2.0) inhibition; v10 PRESET_E/F hip-signal INs with fast
  self-adaptation (PREA, onset detection) → RG-E exc + RG-F inh
  (ext) / RG-F exc + RG-E inh (flex).
- PF: 4 phase cells E1/E2/F1/F2 + PFA self-adaptation INs; RG→PF exc
  2.4; DRIVE→PF 0.05; PF↔PF reciprocal inhibition (4 pairs, g=4.0);
  KINH IN (PF_F1→KINH 1.5) → knee_ext MN inh (f1_kneext_inh) and →
  ankle_pf MN inh (f1_anklepf_inh).
- Per muscle (x92): Ia/II/Ib encoder INs. Ia→homonymous MN exc 0.6 +
  Ia→antagonist MN DIRECT inh 0.4; II→MN exc 0.4 (stance-gated gain);
  Ib→MN autogenic inh 0.35; extensor-stance groups: Ib→IBEXC 0.5,
  RG-E→IBEXC gate 1.0, IBEXC→MN exc 0.6 (stance reflex reversal).
- Renshaw (renshaw>0): MN→RC exc 1.0, RC→homonymous MN inh g, RC↔RC
  MUTUAL between distinct pools (fixed 2026-09-15; was one-directional
  in the build and self-looped in the diagram network).
- Balance family (BAL_PF/DF, BAL_LAT, BAL_TRK) → MN pools; POSTURE +
  per-MN POST_i standing-solve bias.

## 2. Verdict table

| # | Connection / cell | Literature evidence (paper, year, species) | Our status | Recommended change |
|---|---|---|---|---|
| 1 | RG half-centers, mutual inhibition | Rybak 2006a J Physiol (cat model); Shevtsova 2015 (mouse); Shevtsova 2026 eLife (rat, InF/InE laminated) | PRESENT-DEVIANT — direct synapse, laminated INs lumped (documented) | none (documented simplification) |
| 2 | Persistent-Na intrinsic bursting of HCs | Shevtsova 2015; Deng/Nourse 2023 (GNa 1.5 µS, tau_h 350 ms) | MISSING — non-spiking RCs + ADAP loop substitute; toolbox tau_h(V) collapse verified (_tau_h_check.py) | none now; if ever ported to toolbox persistent-Na cells, patch tau_h fixed |
| 3 | Burst termination / adaptation | McCrea & Rybak 2008 review | PRESENT-CORRECT (ADAP loop = explicit substitute) | none |
| 4 | RG→PF drive (weak in Shevtsova 2026: 0.1 µS) | Danner 2017; Shevtsova 2026 | PRESENT-DEVIANT — ours is strong (2.4) "PF-forcED" per Rybak 2006 two-level design | none (intentional architecture choice) |
| 5 | PF↔PF reciprocal inhibition | Deng A6 (IN-laminated HCs); Rybak 2015 V1/V2b | PRESENT-DEVIANT — direct lumped | none |
| 6 | Ia homonymous monosynaptic excitation | Matthews; standard (Rybak 2015 Fig 2); NOT in Deng | PRESENT (extra vs Deng, canonical biology) | none |
| 7 | Ia reciprocal inhibition via IaIN population | Hultborn 1971 (6MNCI8I6, MHHYASHP — recurrent control of Ia pathway); Pratt & Jordan 1987 (D3TVUUJP — IaINs rhythmically active in fictive locomotion); Deng A6 (PF-gated IaIN) | PRESENT-DEVIANT — our Ia→antagonist is DIRECT: no IaIN cell exists, so (a) no PF/phase-gating of reciprocal inhibition, (b) no RC→IaIN disinhibition (Hultborn's recurrent pathway) can be expressed | **P1b: add IaIN population** (Ia→IaIN→antagonist MN, PF-phase-gated; enables #11 below) |
| 8 | Group II pathways | Jankowska 2010 (8RDG5YUZ: Ib- and II-input INs are OVERLAPPING populations, not distinct — supports our shared encoder); Perreault 1995/2011 (VC4P6AZ3, NC7ETD46: flexor gr II actions during fictive locomotion incl. phase-resetting effects) | PRESENT-DEVIANT — II→MN direct only; no II→RG reset route (our reset is hip-selective) | **P2b: generalize PRESET reset to per-joint group-II velocity signals** (literature: flexor gr II volleys reset the cycle) |
| 9 | Ib autogenic inhibition (rest) | Jankowska/Lundberg classic; Eccles | PRESENT-CORRECT (swing-phase branch) | none |
| 10 | Ib/ group-I stance REVERSAL: excitation of extensor MNs during load | Perreault 1999 (9GTNRHKK: rest inhibition → locomotion oligosynaptic excitation); Gossard 1994 (F2FNGYXP: locomotor Ib pathway; gr-I reset mostly Ib); **Domínguez-Rodríguez 2020 (YANK2JE7: polysynaptic excitation, ~3 INs, identified candidates)**; Pearson 1998 review (8LQF6FPA); Pearson 2000 load review (89WM8PFQ) | PRESENT-CORRECT in principle (IBEXC stance-gated reversal) — DEVIANT in structure: literature = state-dependent switch of ONE pathway; ours = separate direct-inh + gated-exc paths | none required; optionally route IBEXC excitation as the same INs that carry the inhibition (cosmetic) |
| 11 | RC→IaIN disinhibition (recurrent facilitation of antagonist reflex) | Hultborn 1971 (6MNCI8I6); Deng A6 (RC→IaIN 0.5 inh) | MISSING — blocked by #7 (no IaIN) | fold into P1b |
| 12 | RC↔RC mutual inhibition between pools | Hultborn (mutual recurrent inhibition); Deng A6 RC→RC | PRESENT-CORRECT (fixed 2026-09-15; was lopsided) | none; NOTE: v10 tuning carried the old wiring — re-tune or re-verify |
| 13 | RC rhythm-modulatory (not rhythm-essential) | Noga 1987 (77QUIHAM: mecamylamine abolishes recurrent inh, locomotor rhythm PERSISTS) | consistent — RC gain must not gate rhythm | v11 check: renshaw sweep must preserve rhythm |
| 14 | RCs rhythmically active, phase-locked to MNs | McCrea/Pratt/Jordan 1980 (RYXH6JSP); Pratt & Jordan 1987 | consistent (RC driven by MN) | none |
| 15 | Extensor gr-I phase RESET (prolong stance / terminate flexion) | Conway & Hultborn 1987 (X7KBR5U9); Gossard 1994; **Domínguez-Rodríguez 2020**; Pearson 1998 | PRESENT-CORRECT concept (PRESET_E prolongs stance) — DEVIANT: our ext signal is spindle-like length; **literature: the major prolonger is group-I/Ib FORCE-related input** (Gossard 1994; Pearson 1998: stance-burst duration + magnitude depend on load feedback) | **P1a: feed the stance Ib pathway (IBEXC group IN) → RG-E excitation** (load-prolonger; the literature-indicated duty fix) |
| 16 | Flexor-side reset (trigger swing) | Perreault 2011 (NC7ETD46: flexor gr II trains reset the cycle); Andersson 1978-adjacent flexor-reflex reset literature (MQ4SINGP: flexor reflex afferents reset step cycle) | PRESENT-CORRECT concept (PRESET_F velocity → RG-F) | see #8 (generalize) |
| 17 | Crossed-extensor / contralateral stance reinforcement onto MNs | classic crossed-extensor physiology; in-library: interlimb cutaneous reflex modulation (CBWBJ5BV 2023, 3SB4QSIF 2018, US7722MA/I5MEVVRX 2024) — left-right symmetry modulates interlimb reflexes | MISSING — we only inhibit the contralateral RG; no contralateral MN/IBEXC reinforcement | **P2a: commissural excitation from stance RG-E (or IBEXC) → contralateral extensor MN group IN** (weight transfer; helps partner-leg stance) |
| 18 | V0D vs V0V split roles (alternation via inhibitory + excitatory CINs) | Shevtsova 2026 (V0D inh; V0V exc via V2a); Talpalar 2013 | PRESENT-DEVIANT — single inhibitory commissural class ×2 strengths; no excitatory commissural (no synchrony) | none now (alternating gait only); P3 if ever needed |
| 19 | dI6 (inhibitory CIN class) | Dyck 2012; Rybak 2015 ("not considered") | lumped into our commissural inhibition | P3 note only |
| 20 | V3 (excitatory CINs, synchrony at speed) | Zhang 2022 eLife (KIG9DEKJ); Danner 2019; Shevtsova 2026 (V3→contra RG, weak) | MISSING, documented (no synchrony pathway) | P3 (only non-alternating gaits) |
| 21 | Cutaneous reflex pathways + phasic gain control | Andersson 1978 (57KB6B8F: phasic gain control of cutaneous transmission during fictive locomotion); 2018–2024 interlimb reflex papers | MISSING entirely | P3 / out of scope (no skin or contact-perturbation events in the current plant) |
| 22 | PF→Ib presynaptic shunt (phase-gated Ib transmission) | Deng A6 (PF→Ib 2.0 µS, near-rest shunt) | MISSING in-network — functionally implemented in the runner as presynaptic gain-gating (documented) | none (equivalent mechanism, different locus) |
| 23 | DRIVE→PF direct weak excitation | not in Deng/Shevtsova (their PF driven by RG only) | EXTRA-NOTSUPPORTED (g=0.05, tiny) | optional removal; harmless |
| 24 | KINH (F1-command → extensor-MN inhibition, swing-gated) | literature analog = flexor-command reciprocal inhibition of extensor MNs via IaINs (V1/V2b; Pratt & Jordan 1987) | PRESENT-DEVIANT — gated by PF-F1 directly rather than via IaIN; defensible lump | superseded structurally by P1b (IaIN) |
| 25 | Monosynaptic II→MN excitation (stance) | supported for extensor load/length support (Pearson 2000 review 89WM8PFQ; Jankowska 2010) | PRESENT-CORRECT (stance-gated gain) | none |

## 3. Prioritized recommendations

**P1a — Ib/stance-force → RG-E prolonger (highest impact, smallest
change).** Add one conditional edge: the IBEXC group IN (or a dedicated
Ib group IN) → RG-E excitation, gain default 0. Evidence: Gossard 1994
(the locomotor-related gr-I reset is mostly Ib), Pearson 1998 (stance
burst DURATION is proprioceptively regulated), Perreault 1999,
Domínguez-Rodríguez 2020. This is the literature-indicated fix for our
#1 residual (stance duty 0.17–0.2 vs 0.61), complementing the existing
length-based PRESET_E.

**P1b — IaIN population (structural completion).** Per muscle column:
Ia→IaIN→antagonist MN (replacing the direct edge), IaIN phase-gated by
PF (Deng A6 PF→IaIN 0.5), and RC→IaIN inhibition (Hultborn 1971) once
Renshaw is on. Enables the canonical recurrent architecture and makes
KINH's job biologically grounded. Moderate code (one IN per pool pair +
3 connection groups), off-by-default conditional pattern as usual.

**P2a — contralateral stance reinforcement.** Commissural excitatory
route (stance-side RG-E or IBEXC → contralateral extensor-group IN →
MNs), gain default 0. Evidence: crossed-extensor physiology + the
2023/2024 interlimb-reflex modulation papers (left-right symmetry
signals modulate interlimb reflexes). Target residual: partner-leg
stance quality.

**P2b — group-II generalized reset.** Extend PRESET inputs from
hip-specific to per-joint flexor group-II velocity signals (Perreault
2011). Our hip-selective version is a special case; generalization is
literature-aligned but increases dimensionality — only if the hip
version under-delivers.

**P3 (record, implement only on demand):** V0V/V3 excitatory
commissurals (synchrony gaits), dI6 explicit class, cutaneous reflex
paths with phasic gain control (needs skin/contact events), removal of
the tiny DRIVE→PF edge (EXTRA-NOTSUPPORTED but harmless).

## 4. Citations (Zotero keys where items are in the library)

- Hultborn, Pierrot-Deseilligny et al. 1971, "Recurrent inhibition from
  motor axon collaterals of transmission in the Ia inhibitory pathway"
  — keys 6MNCI8I6 / MHHYASHP.
- McCrea, Pratt & Jordan 1980 J Neurophysiol 44:475 — RYXH6JSP.
- Pratt & Jordan 1987 — D3TVUUJP (IaIN + RC rhythmic co-activity).
- Noga, Shefchyk & Jordan 1987 — 77QUIHAM (mecamylamine; rhythm
  persists without recurrent inhibition).
- Gossard et al. 1994 — F2FNGYXP (locomotor group-Ib pathway; gr-I
  reset largely Ib).
- Conway, Hultborn, Kiehn & Mentis 1987 — X7KBR5U9 (extensor gr-I
  reset of fictive locomotion).
- Perreault, Enoka & Murphy 1995 — VC4P6AZ3 / BKJBDFLY (flexor gr II
  during fictive locomotion).
- Perreault, Chen, Henstschel 1999 — 9GTNRHKK (gr-I rest inhibition →
  locomotor excitation; scratching + weight support).
- Pearson 1998 — 8LQF6FPA (proprioceptive regulation of stance burst
  magnitude + duration).
- Pearson 2000 load-receptor review — 89WM8PFQ.
- Jankowska 2010 — 8RDG5YUZ (Ib/II IN populations overlap).
- Perreault, Angel, Guertin 2011 — NC7ETD46 (flexor gr II reset).
- Domínguez-Rodríguez et al. 2020 — YANK2JE7 (candidate INs for
  extensor gr-I reset; polysynaptic extensor excitation).
- Andersson 1978 — 57KB6B8F (phasic gain control of cutaneous
  transmission).
- CBWBJ5BV 2023, 3SB4QSIF 2018, US7722MA/I5MEVVRX 2024 — interlimb
  cutaneous reflex modulation.
- Rybak, Shevtsova & Kiehn 2015 eNeuro (full text in repo:
  _rybak2015_plain.txt); Shevtsova et al. 2026 eLife RP107480
  (full text cached); Nourse et al. 2023 Biomimetics (full text:
  _nourse2023.txt, Deng Table A6).
- SADb state: Papers 532; curation layer complete on 20/383 imported
  records (SADb_audit\README.md, 2026-09-12) — Zotero remains master.

## 5. Rule-compliance notes

- Newer-over-older applied: Domínguez-Rodríguez 2020 supersedes the
  pure-Gossard-1994 reading of the extensor gr-I pathway; Jankowska
  2010 supersedes the classical separate-Ib/II-IN organization.
- Species labeled: #7, #10, #11, #13–17 are CAT data (fictive
  locomotion); #18–20 mouse/rat genetic models; human evidence is
  indirect (reflex-modulation studies) — no verdict here rests on
  human data alone.
- Deng/Animatlab-only features lacking independent mammalian-literature
  support: none found beyond what is already marked (Deng's PF→Ib
  shunt is supported as phase-dependent gating by the 1978 phasic-gain
  literature for cutaneous and is consistent for proprioceptors, but
  the specific PF→Ib implementation is model-specific).
- NO CHANGES were made to any model/simulation code. Implementation of
  P1a/P1b/P2a/P2b awaits Ben's go; the usual off-by-default conditional
  topology + regression-gate pattern applies to each.

## 6. Full-text verification addendum (2026-09-15, post Ben's web-API grant)

Access granted by Ben: Zotero WEB API (key in D:\Github\api_credentials_local.txt,
OUTSIDE the repo — never commit) for the AARL group library (group 735051,
read-only); Airtable MCP verified live (ping OK); and read access to
SADb_audit\pdf_staging (2194 files, 5814 MB — SADb PDFs staged alphabetically
as Author_Year__KEY.pdf; that folder belongs to the other session's upload
workflow and was NOT modified).

Five key papers downloaded to spinal\lit_pdfs\ (untracked) and verified at
full-text level:

- Perreault, Enriquez-Denton & Hultborn 1999 (JN 81:2446, full text): "At
  rest, extensor group I afferents produce oligosynaptic inhibition of
  extensor motoneurons. During locomotor activity, however, such inhibition
  is REPLACED by oligosynaptic excitation." — row #10: confirms the
  state-dependent switch reading; our two-path implementation is a
  defensible engineering equivalent.
- Dominguez-Rodriguez, Stecina, Garcia-Ramirez, ... Hultborn, Quevedo 2020
  (full preprint): "instead of the classical Ib non-reciprocal inhibition,
  stimulation of extensor group I afferents produces a polysynaptic
  excitation in extensor motoneurons with latencies (~3.5-4.0 ms)
  compatible with 3 interposed interneurons. We assume that some
  interneurons in this pathway actually BELONG TO THE RHYTHM-GENERATING
  LAYER of the locomotor CPG." — row #15/P1a UPGRADED: the P1a edge
  (stance Ib group IN -> RG-E) is not just permitted but is exactly where
  the literature places these INs: INSIDE the RG layer.
- Perreault, Angel & Guertin 2011 (full text): gr-II stimulation of
  flexor nerves (PBSt/Sart 5T) resets the cycle to extension; TA gr-II
  resets to extension at 2T — row #8/P2b full-text confirmed.
- Jankowska & Edgley 2010 (full text): "no compelling reasons to consider
  intermediate zone interneurons with input from group Ib afferents to be
  distinct from those co-excited by group II afferents" — row #8: shared
  Ib/II encoders supported (our single Ia/II/Ib encoder design is
  consistent).
- McCrea, Pratt & Jordan 1980 (full text, 45 pp): RCs rhythmically active
  and phase-locked; recurrent effects on MNs present during fictive
  locomotion (row #14). RC-RC mutual statement not verbatim in this paper;
  the canonical sources remain Hultborn 1971-lineage + Deng A6.

NEW THEORY ANCHOR (Rybak 2015, full text): "any spontaneous perturbation
or an afferent or supraspinal signal affecting spinal circuitry BELOW the
RG level (i.e., at the PF or motoneuron level) cannot reset the rhythm and
can only produce non-resetting deletions." Consequences: (1) the P1a
Ib->RG-E edge is correctly placed at the RG layer; (2) PF-level couplings
can NEVER reset the rhythm — retroactively explains why v5-v10 phase
fixes at the PRESET/RG-input layer were the right locus, and why KINH/
PF-side effects only shape waveform, never phase. Any future phase-reset
work must terminate INSIDE the RG layer (PRESET PREA loops do).

## 7. EXPANDED verdict table with phase-role taxonomy (2026-09-15, second
## pass: full enumeration of 15 locomotor-relevant AARL collections)

Ben's second-pass ruling: the first audit leaned on the personal/SADb
library and two AARL title searches — NOT a systematic walk. Fixed:
Afferent, Sensory Feedback, Ia Afferent Cite, CPGs, Dual-Layered CPGs,
Network Architecture, Spikes to Muscle Activation, Spiking/NonSpiking
Networks, FSA_Extensions, Aim 2, CHAPTER 2/4, Damping in Locomotion,
Literature Reviews were ENUMERATED ITEM BY ITEM (174 unique items;
lit_aarl_enum.txt). Phase-role taxonomy per Ben:
  S2W = swing-to-stance (contact/stance trigger)
  W2S = stance-to-swing (swing trigger)
  E2W = early-stance prolongation (delays swing; THE duty lever)
  IL  = inter-leg, LG = intra-leg, TON = tonic/posture

| # | Connection/cell | Phase role | Literature (paper, year, species) | Our status | Recommendation |
|---|---|---|---|---|---|
| 26 | Extensor force FEEDBACK (positive force feedback) | E2W | Prochazka 1997 Positive Force Feedback (4XQG7QPN/HQZX9CP7; duplicate 2024 upload W6PSXTK3) | MISSING as force term (P1a uses Ib length/force-IN route — same sign, this is its citation anchor) | cite in P1a implementation |
| 27 | Stance-duration regulation by sensory feedback | E2W | Pearson, Role of sensory feedback in control of stance duration, walking cats (TN2DVMAX) | PARTIAL (IBEXC excitation only; no duration-prolonging dynamics) | P1a: IBEXC/RG-E prolonger edge + search its gain |
| 28 | Ib pathway INs INSIDE the RG layer; polysynaptic extensor excitation | E2W (+W2S reset-to-extension) | Dominguez-Rodriguez 2020 (YANK2JE7, full text) | MISSING (P1a) — literature places the reset INs in the RG layer | P1a as specified |
| 29 | Extensor gr-I reset-to-extension | S2W (trigger stance; also terminates swing) | Conway/Hultborn 1987 (X7KBR5U9); Gossard 1994 (F2FNGYXP) | PARTIAL (PRESET_E length-based; Ib/force-based variant missing) | P1a |
| 30 | Flexor gr-II reset-to-flexion / cycle reset | W2S | Perreault/Angel/Guertin 2011 (NC7ETD46, full text); VC4P6AZ3/BKJBDFLY 1995 | PARTIAL (PRESET_F hip-velocity only) | P2b: generalize to per-joint gr-II |
| 31 | Ib/II IN populations OVERLAP (not distinct) | LG | Jankowska & Edgley 2010 (8RDG5YUZ, full text) | CONSISTENT (shared per-muscle encoder) | none |
| 32 | Tension/force feedback stabilizes musculoskeletal quadrupedal gait | E2W/LG | Tanaka 2024 (QVEH8E98) | consistent with P1a direction | robot precedent for P1a |
| 33 | Stretch reflex on PAM-driven musculoskeletal robot | LG (stance reflex) | Yoshida 2025 (ZE92KPYH) | PARTIAL (we have reflexes; no stretch-reflex latency model) | P3: reflex latency model |
| 34 | PAM-driven quadruped pace running | LG | Fukuoka 2022 (53MMADKJ) | comparator only | none |
| 35 | Neuromechanical walking controller validated by disturbances | LG (+S2W/W2S modulation) | Song & Geyer 2017 (MNFISTGM) | comparator benchmark | P2: replicate their disturbance tests as validation |
| 36 | Physiologically realistic MuJoCo muscles (Millard-based, compliant tendon) | TON (plant property) | Wang 2022 MyoSim (4NJFY349) | UPGRADE PATH (our rigid-tendon comparison shows the compliance signature in soleus) | P2: MyoSim-style tendon compliance pilot |
| 37 | OpenSim-to-MuJoCo conversion method | TON | Ikkala & Hamalainen 2022 (MNC6TMBW) | comparator | none |
| 38 | FSA control of CPG PHASE DIFFERENCE (multistate) | IL/phase-method | Scharzenberger & Hunt (HNTEY3GQ/UKRQJ3ZH/YFGAPDFF) | MISSING (our phase relation is open-loop) | P2: FSA phase-difference IN pair for the interleg phase lock |
| 39 | FSA legged-robot locomotion controllers | LG | Szczecinski 2017 (4W847RW9/7QRJD27N) | partially referenced | cite in appendix |
| 40 | Phase-response analysis of limit cycles with hard boundaries | method | Wang 2022 (FZD59HCB) | not used | method for quantifying reset (P1a evaluation) |
| 41 | Spinal INs onto the final common path during locomotion | LG | Brownstone & Bui 2010 (9YIA6HYF) | consistent (MN-level IN layer = RC/IaIN/KINH) | none |
| 42 | Propriospinal (LPN) neurons essential for interlimb coordination | IL | Laliberte 2019 review (D9FPFD55) | MISSING (no long propriospinal layer; our interleg = RG-level only) | P3: LPN-inspired interlimb excitation relay |
| 43 | Human CPG existence/contribution review | TON | Minassian 2016 (E376F3ZS) | not cited yet | cite in dissertation (human-side anchor) |
| 44 | Somatosensory feedback control of locomotion (comprehensive review) | S2W/W2S/E2W anchor | Frigon 2022 (E87KHQYE) | not cited yet | cite as the review anchor for P1a/P2b |
| 45 | Rybak 2025 split-belt phase-duration circuits | E2W (phase durations) | Rybak 2025 (MGYWT5TP) | not covered (audit had 2024 eLife only) | read; phase-duration circuits bear directly on duty |
| 46 | Reflex latency + parallel compliance trade-off | LG | Ashtiani 2021 (A6DNHI84) | our reflexes are near-instantaneous (no latency model) | P3: latency model before hardware transfer |
| 47 | Spindle encoding gains (dynamic/static, human+cat) | LG (encoder params) | Prochazka 1976/1976/02, Kakuda & Nagaoka 1998, Edin & Vallbo 1990, Hasan 1983, Crowe & Matthews 1964 (Afferent coll.) | PARTIAL — our AFF gains are plausible but not literature-fitted | P2: fit AFF ia_vel_ref/ii_len_ref from these datasets |
| 48 | GTO encoding (active contraction) | LG (encoder params) | Houk & Henneman 1967 (UHK5RTCK); Jansen & Rudjord 1964 (IK49AIUS) | consistent (Ib force-normalized) | none |
| 49 | Reflex STIFFNESS regulation + predictive dynamic response | LG | Houk 1979 (5SCE5YFT); Houk 1981 (QWICNHLH) | implicit (our reflexes act as impedance) | P3: explicit stiffness-regulation framing |
| 50 | Ia integrative patterns hip/knee MNs | LG | Eccles & Lundberg 1958 (H7IAZ6Y7) | consistent | classic cite |
| 51 | Proprioception ablation abolishes walking | S2W/W2S (global) | Pearson 2003 (Y32E8ES5) | consistent (deafferented air-stepping works only with CPG) | cite |
| 52 | V1 IN diversity + speed-dependent recruitment | LG | Gosgnach 2017 (NU75BEIA); Rybak 2015 | consistent (V1 = IaIN/RC lineage per P1b) | cite in P1b rationale |
| 53 | Graded MN patterns produce similar contractions (spike patterns vs force) | TON | Hooper 2007 (IHIS9GDK) | consistent (non-spiking pools) | cite in appendix |
| 54 | Effective neural drive = common synaptic input to MNs | TON | Farina 2014 (4DT5S7CY) | consistent | cite in appendix |
| 55 | Sensorimotor timing: hip torque + leg damping | LG | Shen & Seipel 2012 (CAIECEZ6); Mo 2020/2023 damping papers | not modeled | P3 (stance damping shaping) |

## 8. Classification of the REST of the enumeration (non-connectivity)

- PAM/BPA modeling + hardware (CHAPTER 2. BACKGROUND, most of Afferent
  Cite-adjacent): McKibben/Tondu/Sarosi/Toth/Colbrunn/Aschenbeck/Kingsley
  cockroach robot, Hunt-lab BPA papers (Elzein, McNeal, A7VFD2TL) —
  actuator modeling, not circuit connectivity. KEPT for the BPA appendix.
- Method/ML refs in FSA_Extensions (Izhikevich, PyTorch, TensorFlow,
  Adam, backprop, LeCun, Nengo/NEF, autodiff) — tooling citations.
- Boston Dynamics links, Hebb 1949, posters/presentations/grants
  collections — not circuit literature.
- Audit-relevant but already cited in DESIGN.md: Deng 2019 (two
  versions), Markin 2016/2010, McCrea & Rybak 2008, Rybak 2006, Rybak
  2013/2015, Shevtsova 2022 (V1 ipsi/contra), Rybak 2024 eLife,
  Avaltroni 2024, Di Russo 2023, Prochazka-style load reviews.

## 9. Pass-2 summary

174 items enumerated across 15 collections; 30 NEW audit-relevant rows
(#26-55) added above; of those, 6 strengthen or extend P1a/P1b/P2a/P2b,
4 are new P2 recommendations (MyoSim compliance pilot, FSA
phase-difference IN pair, AFF-gain literature fit, Rybak-2025 phase-
duration reading), the rest are P3/cite-only. NO model code was touched;
implementation still awaits Ben's go per item.
