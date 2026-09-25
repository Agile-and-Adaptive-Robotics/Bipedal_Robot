# Feedback wiring audit: what targets the half-centers vs the layer interneurons? (2026-09-24)

Ben's question: "what was connected directly to the HC neurons and what was
connected to that layer's interneurons? The literature is a little vague."

Method: mined the ORIGINAL Deng walker directly from git history —
`AddingStepSensor_CoMorrow_stw` @ **bb192d48494fd04a14022ddbdf9f79293acfe4da**
("Added two models"; parent 407fa01 "Connected RG, PF, and Sensorimotor
Layers"). Files extracted read-only with GitHub Desktop's git to
`%TEMP%\deng_orig\` (W2L + Biped_2xCPG aproj at that commit; links CSV +
summary there too). Polarity = each SynapseType's EquilibriumPotential vs
the −60/−70 mV resting potentials (eq > −55 mV ⇒ excitatory).

## 1. Original Deng walker (W2L @ bb192d48) — the precise answer

**Nothing sensory touches an RG half-center. Nothing touches a PF
half-center either.** Every afferent lands at the reflex/MN level, and the
only inputs to the half-centers are intra-layer:

RG layer (L RG ext / L RG flx + RG-INs):
- HC↔HC: weak mutual EXCITATION g=0.1 ("RG to RG Excite", eq −40)
- HC→own IN EXC 2.749; IN→antagonist HC INH 2.749 ("RG Inhibit", eq −70)
- that is ALL. (The 2026 file adds only a kickoff TonicCurrent stimulus.)

PF layer (L Hip PF ext/flx + PF-INs):
- RG→PF EXC 0.1 ("RG to PF Excite")
- PF HC→own PF-IN EXC 2.749; PF-IN→antagonist PF HC INH 2.749
- that is ALL.

Afferent routing (hip joint only in this snapshot):
- Ia: StretchReceptor → adapter → Ia relay 1 (Hip ext/flx Ia) EXC 0.59 →
  Ia relay 2 ("Ia 2") EXC 0.59 → **antagonist MN, INH g=2** ("Ia to MN
  Inhibit", eq −100). Reciprocal inhibition AT THE MN.
- Ib: muscle force → Ib adapter → Ib cell → **own MN, EXCITATORY**
  (ext Ib→MN ext 0.59; flx Ib→MN flx 1.0 — the depolarizing autogenic
  choice AnimatLab-era; NOT the classic sign)
- Phase gating: **PF HC → Ia relay 2 EXC 0.5** ("PF to Ia Excite") —
  the PF layer gates the afferent RELAY neurons, not the HCs.
- Renshaw: MN→own RE EXC 0.5; RE→own MN INH 0.5; RE↔RE INH 0.5; and
  **RE → Ia relay 2 INH** (0.5 / 0.4232) — recurrent
  disinhibition of the antagonist via the relay.

So the original doctrine: **the CPG is deafferented; feedback sculpts
OUTPUT (MNs) and gates the afferent relays themselves.** The "subnetwork
layer" organization (per-joint Ia/Ib/PF/MN/RE chains) carries all of it.

## 2. Current s3k stack (build_network.py) — where it differs

- heel/toe → heel_IN → InE/InF lamination → RG (via INs; R5 route) —
  same doctrine as the original (via INs).
- stance-Ib → LBIN → RG-E (via IN).
- Ia → MN mono + IaIN → antagonist (via IN) — same as original but
  adding the mono corner.
- **II and Ia "central" knobs (ii_f_central, ii_e_central, ia_f_central)
  connect DIRECTLY to RG-F/RG-E HCs** — the clearest deviation from the
  original (which had zero afferent→HC) and from the literature pattern
  (aff → IN → center). These are the knobs to revisit.
- v5 phase reset (hip ext/flex signals → PRESET_E/F INs → RG HCs) —
  via dedicated INs, consistent with the Rybak phase-reset literature.
- PF structure: default = 4 phase-window cells/side (PF_E1/E2/F1/F2,
  single-neuron cells with adaptation) + PF_IN_E/F cross-inhibition;
  the **joint_pf>0 variant (built 09-18) already provides antagonistic
  per-joint HC pairs (HIP-E/F, KNEE-E/F, ANK-E/F)** — s3k production
  runs the phase-cell PF, not the joint variant.
- PF→MN: **every one of the 92 actuators has its own MN** and receives
  weighted PF input (W_PF_MN group weights, primary + 0.5×secondary);
  "six MNs with the same signal" was the EDITOR TEMPLATE's
  simplification, now fixed (templates redrawn per-muscle today).

## 3. Literature position (per our replication drafts; M3 gate applies)

- Shevtsova 2026 laminar core: NO afferents at all (core is
  deafferented; _crossref note in shevtsova_rules.json) — matches the
  original W2L doctrine.
- Deng A6 reflex layer (circuit_dengstyle / draw_circuit): Ia→IaIN→MN,
  Ib→IBIN→MN, II via collateral INs — afferents target INs and MNs, not
  HCs. Matches.
- Rybak 2006a: afferent access is to PF and MN levels; RG phase-reset
  pathways exist (hip afferents) but run through interneurons — matches
  the "via INs" doctrine, and is the citation base for PRESET_E/F.
- Shinohara 2025: per-muscle afferents INTO the centers — the exception
  that justifies per-muscle afferent fans if we go that way.

## Bottom line

The original answers the vague literature concretely: **HCs take no
afferents — layer INs and MNs take everything; PF gates the afferent
relays; Renshaw closes loops through the relays.** Our current stack
agrees except the direct II/Ia→HC central knobs. The joint-layer
antagonistic PF variant Ben asks for exists behind G["joint_pf"]; the
per-muscle MN diagram he asked for is in the editor templates as of
today (walker_v9…s3k redrawn: 233 nodes / 420 edges for s3k).
