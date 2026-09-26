"""All tunable parameters for the gait2392 spinal network, in one place.

Units follow sns_toolbox conventions: neuron inputs are nA, potentials mV,
time constants = membrane_capacitance / membrane_conductance in seconds
(the runner steps at dt = 0.005 s, so keep capacitances >= ~0.01 for stability).

Afferent/gain structure follows the references in Ben's Zotero
"Sensory Afferent Database" collection (see DESIGN.md):
  - Ia (spindle primary): dynamic, ~ dL/dt, monosynaptic excitation to
    homonymous MN + reciprocal inhibition of the antagonist (Bergmann et al.
    2007 review; Ross & Rybak mouse spindle-speed work 2018 for the
    speed-dependence).
  - II (spindle secondary): static length signal, excitation to homonymous MN,
    strongest in stance (load/length support).
  - Ib (GTO): force signal; autogenic inhibition by default, but during
    stance the extensor Ib signal is routed into group "load-sharing"
    excitatory interneurons (reflex reversal; Pearson 2000 / load-receptor
    review 2000).
Phase gating of afferent gains is implemented presynaptically in the runner
(gains multiplied by stance/swing gate read from RG potentials) - a stand-in
for presynaptic inhibition from the RG, documented in DESIGN.md.
"""
from __future__ import annotations

# ------------------------------------------------------------- global scales
E_HI = 5.0          # mV, "full activation" membrane potential for MNs
DT = 0.002          # s, matches the MJCF timestep (2 ms for stability)

# ------------------------------------------------------------- neuron time constants
TAU = dict(
    mn=0.03,        # motoneuron pool membrane (30 ms)
    afferent=0.02,  # afferent encoders (20 ms filter on the raw signal)
    rg=0.05,        # rhythm-generator half-center membrane (50 ms)
    rg_nap_h=0.35,  # RG persistent-Na h-gate time constant (FIXED,
                     # Deng semantics; ~3.2x this = cycle period: 0.35 s
                     # -> ~1.2 s. The PERIOD knob, measured 2026-09-16)
    pf=0.08,        # pattern-formation cells
    pf_adapt=0.5,   # (legacy: PFA removed 2026-09-16; kept for record)
    ib_exc=0.05,    # group Ib load-sharing interneuron
    descend=0.10,   # descending drive smoothing
    preset=0.04,    # (legacy: PRESET removed 2026-09-16; HEEL/TOE INs
                     # still use this value)
    preset_adapt=0.08,  # (legacy: PREA removed 2026-09-16)
    rg_adapt=1.9,   # (legacy: ADAP retired 2026-09-16; inert, kept so
                     # old jsons/studies still load without KeyError)
)

# ------------------------------------------------------------- RG NaP neurons
# Persistent-Na half-center parameters (Deng 2022 / Shinohara 2025 /
# Rybak 2024-25 conditional bursters), rescaled to OUR 0..5 mV operating
# range (e_ion near the plateau instead of ENa +50; slope_m/e_m set the
# activation threshold ~2 mV). tau_max_h is the burst-termination /
# PERIOD knob (fixed-tau semantics; tau_h ~0.35 s -> cycle ~1.2 s,
# ~2x tau_h scaling, measured _nap_fixed_test.py 2026-09-16).
NAP = dict(
    g_ion=12.0,        # NaP conductance (uS-scale on our range)
    e_ion=8.0,         # NaP reversal (plateau ceiling, mV)
    k_m=1.0, slope_m=0.8, e_m=2.0,    # fast activation gate
    k_h=1.0, slope_h=-2.0, e_h=3.5,   # slow inactivation gate
    tau_max_h=0.35,    # s, FIXED h time constant (Deng semantics)
)

# ------------------------------------------------------------- synapse conductances
# (max_conductance values; inhibitory synapses get reversal_potential=-E_HI)
G = dict(
    # half-center rhythm: mutual inhibition + slow self-adaptation
    rg_mutual_inh=4.0,       # RG-E <-> RG-F
    rg_adapt_inh=2.5,        # RG -> ADAP -> RG (slow negative feedback)
    rg_to_pf=2.4,            # half-center -> its two PF groups (2x: slow
                              # DRIVE keeps amplitudes saturated)
    pf_recip_inh=4.0,        # conflicting PF groups (E2<->F1, F2<->E1)
    pf_to_mn=2.0,            # PF cell -> MN (scaled per-muscle by W_PF_MN)
    posture_to_mn=1.0,       # POSTURE -> MN (scaled by W_POSTURE)
    ia_to_mn=0.6,            # Ia excitation, homonymous (nA-gain knob below)
    ia_to_antagonist=0.4,    # Ia reciprocal inhibition
    ii_to_mn=0.4,            # II excitation, homonymous
    ib_to_mn_inh=0.35,       # Ib autogenic inhibition (swing)
    ib_group_exc=0.5,        # Ib afferent -> group load-sharing IN
    ib_exc_to_mn=0.6,        # group IN -> extensor MNs (stance)
    bal_lat_to_abd=0.8,      # BAL_LAT -> stance-side hip abductor MNs
    bal_trunk=1.0,           # IMU trunk controller -> ercspn / obliques
    descend_to_rg_e=1.7,     # DRIVE -> RG-E (2x with slow ADAP: amplitude
                              # saturated at low DRIVE; E-duty knob - NOTE
                              # raising this SHORTENS the E burst via
                              # stronger adaptation-triggering)
    descend_to_rg_f=1.4,     # DRIVE -> RG-F
    posture_to_rg_e=0.8,     # POSTURE -> RG-E (standing stance bias)
    drive_to_pf=0.05,        # DRIVE -> PF cells (was 0.4: the tonic term
                             # kept all PF cells partially on through their
                             # OFF phase - constant co-contraction, no knee
                             # swing; 0.2 lets the windows close)
    posture_gain=1.0,
    # ---- v5 sensory phase-reset (2026-09-12 night): hip afferent signals
    # resetting the RG phase (the #1 lever from the v4 diagnosis: E-duty
    # 0.27 vs 0.61, cadence 1.18 vs 0.81 Hz, cycle-to-cycle jitter).
    # HIP_EXT_SIG (hip-extensor group length, stance-gated) EXCITES RG-E and
    # INHIBITS RG-F (prolongs stance); HIP_FLEX_SIG (hip-flexor group
    # positive shortening velocity) EXCITES RG-F and INHIBITS RG-E
    # (triggers swing). Wired through per-side PRESET_E/F interneurons.
    # DEFAULTS 0.0 = network behavior-identical to v4b (regression-gated).
    phase_reset_e=0.0,       # HIP_EXT_SIG gain (nA signal -> conductance)
    phase_reset_f=0.0,       # HIP_FLEX_SIG gain (nA signal -> conductance)
    # v5b/phase-3 swing-knee quad suppression (lever #2 from the v4
    # diagnosis: swing knee stays extension-dominant because F1 knee-flex
    # gains saturate against quad tone): F1 -> KINH inhibitory interneuron
    # -> knee_ext MN pools, phase-gated by F1 itself (F1 IS the swing
    # window). Default 0.0 = absent (conditional topology, v4-identical).
    f1_kneext_inh=0.0,       # KINH -> knee_ext MN inhibition
    # v6b: same swing-gated suppression applied to the ankle plantar-
    # flexor pools (reuses the F1-driven KINH interneuron). Motivated by
    # the air-stepping ankle test 2026-09-14: boosting swing DF drive
    # (x3-x5) does NOT produce dorsiflexion (9 PF muscles ~10 kN vs 3 DF
    # ~1.6 kN), but suppressing PF tone in swing should - the mechanism
    # that fixed the swing knee.
    f1_anklepf_inh=0.0,      # KINH -> ankle_pf MN inhibition
    # Renshaw recurrent inhibition (Ben's go, 2026-09-13; Deng Table A6:
    # MN->RC 0.5 exc, RC->MN 0.5 inh, RC<->RC 0.5 inh). Gain applies to
    # RC->homonymous-MN and RC<->RC; topology built only when > 0.
    renshaw=0.0,
    # Ankle standing-tone trim during walking (2026-09-14): the soleus/
    # tib_post posture overrides held a ~-45 deg PF ankle set-point
    # through gait (the tiptoe offset in the v9 overlay). This factor
    # multiplies the POST bias of ankle_pf-group muscles as drive rises:
    # 1.0 = v9 behavior (identical), 0 = no standing PF tone while
    # walking (physiologic: soleus tonic EMG drops with locomotor drive).
    ankle_post_walk_trim=1.0,
    # ---- v11 mechanosensory stance feedback (Ben's go 2026-09-15;
    # implements audit P1a + P1b, LIT_CIRCUIT_AUDIT.md rows 7/11/27/28/29).
    # ALL DEFAULTS 0 = topology absent = v10-identical (regression-gated).
    # P1a: heel/toe contact mechanosensors + stance-Ib prolonger.
    #   HEEL contact -> RG-E exc + RG-F inh (S2W trigger, Conway 1987)
    #   TOE contact -> RG-E exc (late-stance prolongation)
    #   LBIN (RG-layer stance-Ib group IN, per Dominguez 2020 "INs belong
    #   to the rhythm-generating layer") -> RG-E exc (duty prolonger;
    #   Gossard 1994/Pearson 1998 stance-duration regulation)
    heel_rge=0.0,            # heel contact -> RG-E exc / RG-F inh
    toe_rge=0.0,             # loaded toe -> RG-E exc
    ib_rge=0.0,              # stance-Ib group IN (LBIN) -> RG-E exc
    # ---- 2026-09-20 contact-EVENT (onset/offset) shaping of the same
    # heel/toe ports (runner-side; DESIGN.md stage-3 latch hypothesis:
    # a tonic-driven half-center will not oscillate under load - the
    # working reference RG is CONTACT-driven, AnimatLab lesson). 0 = off
    # (bit-identical); > 0 adds a decaying transient to HEEL_c/TOE_c at
    # the LOADING edge (heel strike, +) and UNLOADING edge (toe-off, -)
    # of each foot: brief E trigger at strike, brief E release + F
    # disinhibition at lift, instead of a tonic contact level that holds
    # the half-center. Units: nA of peak port current per unit gain.
    contact_onset=0.0,
    # ---- 2026-09-21 crossed swing trigger (t54 diagnosis: the left
    # stance leg never releases - its own unloading edge never fires
    # because the foot never unloads; meanwhile the v10 pf_gain 0.41
    # leaves PF drive below MN threshold for an unloaded leg). When
    # side A's foot LOADS (heel strike), the CONTRALATERAL side B gets
    # a decaying NEGATIVE kick on its HEEL port: heel_in_B -> rg_e exc
    # / rg_f inh, so -kick inhibits B's extensor half-center and
    # disinhibits its flexor = reset-to-swing (Aoi/Di Russo phase-
    # resetting family; the loading edge of the moving foot is the
    # reliable event). 0 = off (bit-identical).
    contra_swing=0.0,
    # ---- 2026-09-21 crossed KINH drive: contralateral heel-load IN ->
    # this side's KINH -> knee_ext/ankle_pf MN suppression (opposite
    # heel strike forces THIS leg's stance->swing transition at the MN
    # level). Needs f1_kneext_inh > 0 (the KINH cell's own-MN synapses
    # carry the f1 gains). 0 = edge absent (bit-identical build).
    contra_kinh=0.0,
    # ---- 2026-09-21 FULL LITERATURE CONNECTOME (G["full_rules"] > 0;
    # Ben: "wire it the way the working Deng model does"). Adds the
    # missing interneuron layer: II -> IIX (exc IN) -> same MN;
    # II -> IIIN (inh IN) -> antagonist MNs; Ib -> IBIN (inh IN) -> same
    # MN (direct II/Ib->MN edges REMOVED when on); Ib-IN <-> antagonist
    # Ib-IN and IaIN <-> antagonist IaIN mutual inhibition; heel/toe
    # routed through the laminated InE/InF layer instead of direct
    # half-center edges. 0 = legacy direct wiring (bit-identical).
    full_rules=0.0,
    # ---- 2026-09-21 per-side CONTACT-RESET PHASE MACHINE (the Di
    # Russo eq-7 analog; the build the s3c/s3d/s3e verdicts point to):
    # a runner-side phase variable per leg, advanced at 1/pm_T cycles/s,
    # RESET to 0 at that foot's own heel strike (loading onset), with a
    # gentle antiphase pull toward a half-cycle offset between legs.
    # The phase gates the MN OUTPUTS: inside the swing window
    # (phi 0.62-0.95, raised-cosine edges) extensor ctrl scales DOWN
    # by pm_gain and flexor ctrl scales UP by 0.6*pm_gain. This
    # guarantees every leg a swing window each cycle — the thing the
    # frozen stance leg never gets from the SNS half-centers alone.
    # pm_gain = 0 = OFF (bit-identical). If it works, port into SNS
    # topology (phase-oscillator neurons) per DESIGN.md 2026-09-21.
    pm_gain=0.0,
    pm_T=1.2,                # s per cycle (searched)
    # ---- 2026-09-21 pm v2 WEIGHT-SHIFT: (a) the swing window is
    # LOAD-GATED - the phase HOLDS at 0.55 until the CONTRALATERAL
    # foot carries >= 25% BW (only swing the left once the right
    # actually bears weight - real gait's lateral weight transfer);
    # (b) during stance prep (phi 0.42-0.62) the upcoming-swing side's
    # hip abductors scale DOWN and the upcoming-stance side's scale UP
    # by pm_ws (the lateral hip strategy that moves the CoP over the
    # next stance foot). 0 = off (v1-identical).
    pm_ws=0.0,
    # ---- 2026-09-21 pm v3 ADDITIVE flexor drive: the v2 flexor boost
    # was MULTIPLICATIVE on ctrl, and the swing-side flexor ctrl is ~0
    # (MNs below threshold - the original passive-flail diagnosis), so
    # multiplying ~0 stayed ~0: the "swing drive" never drove. pm_add
    # ADDS pm_add * w to flexor ctrl during the swing window (a real
    # flexor burst, forceful enough to break the loaded jam). The
    # multiplicative extensor cut (pm_gain) stays - it works on large
    # values. 0 = legacy v2 (bit-identical).
    pm_add=0.0,
    # ---- 2026-09-21 pm v4 AFFERENT DISFACILITATION: during a side's
    # swing window, scale that side's load-afferent INPUT channels
    # (HEEL_c, TOE_c, LOAD_c, AFF_E, AFF_F) by (1 - pm_aff*w). The
    # s3c..s3h verdict: swing commands (suppression + bursts) lose to
    # the leg's own load afferents re-latching RG-E/InE/MNs - a real
    # swing leg is UNLOADED, so its afferents should be silent; enforce
    # sensory consistency with the commanded swing. 0 = off.
    pm_aff=0.0,
    # P1b: IaIN population replaces the direct Ia->antagonist edge when
    # > 0 (Deng A6: Ia->IaIN->MN with PF_F1 phase gate; RC->IaIN inh
    # = recurrent disinhibition, Hultborn 1971).
    ia_in=0.0,
    # ---- 2026-09-16 interleg gain knobs (Shinohara c1/V3 commissurals).
    # Measured: at full strength (c1 1.0 / V3 0.5 of rg_mutual_inh) the
    # NaP network bilateral E-latches on ground (E-duty 1.0, knees pinned
    # -4..+11) - both commissurals are searched in stages 2-3 instead.
    c1_gain=1.0,              # RG_F -> c1 -> contra RG_F (antiphase lock)
    v3_gain=0.0,              # RG_E -> V3 -> contra InE (extensor sync;
                              # 0 = pathway absent - the latch culprit)
    # ---- 2026-09-16 (Ben, from Deng 2022 / Shinohara 2025 figures):
    # weak MUTUAL EXCITATION between the RG half-centers (Deng 2022 G_W;
    # raises the inhibited neuron's equilibrium = escape mode, and
    # neuromodulation of it raises frequency). Direct RG-E<->RG-F edges.
    rg_weak_exc=0.0,
    # ---- 2026-09-18 (T1/T3 joint-layer PF experiment, fsa_jointlayers.py):
    # replace the 4 phase-window PF cells (E1/E2/F1/F2) with 3 JOINT-layer
    # PF half-center pairs (HIP-E/F, KNEE-E/F, ANK-E/F; 6 HCs) whose
    # PF->MN weights come from the structured joint-layer fit
    # (joint_pf_weights.json). 0 = absent topology (bit-identical build);
    # 1 = layer HCs built. Fitted held-out centered VAF 0.747/0.713 vs
    # 0.934/0.925 unconstrained 6-synergy (T1 > merged-layer variants).
    joint_pf=0.0,
    # ---- 2026-09-24 evening PER-PF-LAYER CONTACT VARIANT (Ben: "build
    # it", coexists with the per-joint layering; source = his block-editor
    # drawing Neuromechanical_Models\Mujoco_SNS_models\
    # Circuit_rules_CONNECTOME_md__connectome.json, 97n/109e). Semantic
    # ruling: heel contact = stance-phase reset of the IPSILATERAL leg,
    # applied AT THE PF LAYER (not only the RG); toe contact ONLY
    # inhibits dorsiflexion; flexion afferent INs reinforce their own
    # joint's F half-center. Requires joint_pf > 0 (the micro-layers
    # these edges target exist only in that build) EXCEPT
    # heel_in_f_exc, which needs only the RG lamination. ALL DEFAULTS 0
    # = edges/neurons absent = bit-identical build (410/376/1186).
    heel_pf_layer=0.0,       # heel IN -> PF_IN_E exc (PF-layer reset;
                             # Ben's drawing g 0.5 to each micro-layer's
                             # E-lamination IN — ours is the one shared
                             # PF_IN_E per side, a documented lumping)
    toe_df_inh=0.0,          # toe IN -> TOEDF IN exc -> ANK-F inh
                             # (dorsiflexion inhibition ONLY; drawing:
                             # toe 5 -> IN-PF_dorsiflexion_inhibit ->
                             # inhib HC-PF-Dorsiflexion)
    heel_in_f_exc=0.0,       # heel IN -> InF EXCITATORY (drawing g 0.5;
                             # the full_rules branch has heel -> InF INH
                             # — the two coexist, this one adds the exc
                             # variant Ben drew)
    ia_pf_f=0.0,             # flexor-group IaIN -> same-joint PF-*-F exc
                             # (drawing: IN-IaIN -> HC-PF-F g 0.5;
                             # requires ia_in > 0 + joint_pf)
    ii_pf_f=0.0,             # flexor-group II exc IN (IIX) -> same-joint
                             # PF-*-F exc (drawing: IN-IIe -> HC-PF-F
                             # g 0.5; requires full_rules + joint_pf)
    # per-muscle afferent -> central feedback (Deng 2022 / Shinohara 2025
    # wiring): extensor muscles' Ib afferents project EXCITATORY to the
    # ipsilateral E-centers (PF_E1/E2, RG_E, InE) - force feedback
    # shifts the V-nullcline (Shinohara 4.2); flexor muscles' Ia+II
    # afferents project to the F-centers (PF_F1/F2, RG_F, InF) - length
    # feedback (same escape logic). Foot mechanosensors share the
    # extensor pathway (Ben). DEFAULTS 0 = absent.
    ib_e_central=0.0,        # extensor Ib -> PF_E / RG_E / InE (exc)
    ia_f_central=0.0,        # flexor Ia -> PF_F / RG_F / InF (exc)
    ii_f_central=0.0,        # flexor II -> PF_F / RG_F / InF (exc)
    ii_e_central=0.0,        # extensor II -> PF_E / RG_E / InE (exc,
                              # same-group; Ben 2026-09-16)
    # Rybak 2025 SF-E1 contralateral half: hip-flexor stretch afferents
    # also INHIBIT the CONTRALATERAL F half-center (promotes the E->F
    # transition + interleg coordination; grows with speed). Per-muscle
    # flexor Ia -> contralateral RG_F inhibition.
    ia_f_contra_f=0.0,
    # Rybak 2025 SF-E2 / audit P2a: the E-side commissural (V3) also
    # reinforces the CONTRALATERAL extensor MN groups (crossed-extensor
    # weight support): V3 IN -> contralateral IBEXC group INs.
    v3_to_ibexc=0.0,
    # ---- v11b semi-closed sensory loops (Shevtsova central principle):
    # afferent relay INs create three-layer positive feedback: muscle ->
    # afferent -> PF -> RG -> MN -> muscle. AFF_E receives extensor-side
    # afferent drive (force); AFF_F receives flexor-side (velocity/length).
    # ALL DEFAULTS 0 = absent.
    aff_e_rg=0.0,            # AFF_E -> RG-E exc (extensor afferent -> RG)
    aff_f_rg=0.0,            # AFF_F -> RG-F exc (flexor afferent -> RG)
    aff_e_pf=0.0,            # AFF_E -> PF-E exc (extensor afferent -> PF)
    aff_f_pf=0.0,            # AFF_F -> PF-F exc (flexor afferent -> PF)
    # ---- 2026-09-23 goal2 STANDING-BALANCE STAGE (Ben's request; modeled
    # on SCONE Tutorial 3a "Balance" = autogenic length reflexes + a
    # vestibular torso-point PD; see reports_20260923/goal2_balance_stage.md).
    # ALL DEFAULTS 0 = VEST cells absent + II loop unchanged = today's
    # behavior (byte-identity-at-0 contract).
    # vest_ext: VEST_{r,l} vestibular-analog cells (runner feeds them the
    #   rectified tilt-deviation/rate current) -> ipsilateral EXTENSOR-
    #   group MN excitation (antigravity tone; Ben: "vestibular analog =
    #   pelvis-tilt sensors driving extensor tone"). Directional ankle
    #   strategy stays with the existing BAL_PF/BAL_DF cells.
    # vest_flex_inh: same cells -> flexor-group MN inhibition (LVST
    #   reciprocal flexor inhibition, standard physiology - AWAITING
    #   Ben's connectome-spec confirmation).
    # vest_prop: stance-gated boost on the II length loop (SCONE T3a KL
    #   length-feedback analog; "proprioceptive balance = stance-gated
    #   Ia/II length-load loops"). Runner-side presynaptic gain, >0 only.
    vest_ext=0.0,
    vest_flex_inh=0.0,
    vest_prop=0.0,
    # ---- 2026-09-25 goal4 VARIANT SELECTOR (default 0 = stock build,
    # byte-identical; regression gate 410/376/1186): > 0 routes
    # build_network.build() to the 6-synergy variant builder
    # build_network_syn6.py (one RG per side driving six synergy PF
    # layers, Eq-18 W->conductance mapping). Env AARL_NET=syn6 selects
    # the same variant without this key.
    syn6=0.0,
    # goal4 syn6 sub-key (default 0 = OFF): Shevtsova brainstem
    # gamma/alpha cells folded onto the single DRIVE port. Measured
    # E-latch in the full runner air run when ON at fixed 0.5 gains
    # (reports_20260925/goal4_build_syn6.md) - the runner's
    # descend_to_rg_e/f edges already carry the descending-drive role.
    syn6_brainstem=0.0,
)

# Extra gains for the phase-reset pathways (not searched by default; the
# two search dimensions are params.G phase_reset_e / phase_reset_f).
PHASE_RESET = dict(
    inh=1.0,             # inhibitory branch scale relative to excitatory
    stance_gate=(0.3, 0.7),  # ext-signal gate = a + b * stance (II-style)
    adapt_g=1.5,         # v10 PRESET fast-adaptation loop gain (onset
                         # detection: PRESET -> PREA -> PRESET, fast tau)
)

# ------------------------------------------------------------- afferent encoding gains
# Current (nA) = baseline + gain * normalized_signal.
# Ia and Ib are PURE-SIGNAL encoders (no baseline: Ia fires on stretch
# velocity only, Ib on force only) - a resting tone here would drive
# constant reciprocal inhibition and crush the antagonist MN pools.
AFF = dict(
    i0_ii=1.0,          # nA baseline tone, II only (static stretch receptor)
    ia_vel_ref=0.6,     # m/s of tendon velocity -> full-scale Ia
    ii_len_ref=0.04,    # m of length deviation -> full-scale II
    ib_force_ref=0.5,   # fraction of Fmax -> full-scale Ib
    ia_gain=1.6,        # nA per unit velocity signal (speed-modulated, see MOD)
    ii_gain=1.2,
    ib_gain=1.5,
)

# Speed-dependent reflex modulation: gains_scale(t) = 1 + MOD * drive(t).
# This is the "reflex pathways change speed" knob (Ben; Bunz/Ijspeert/Schmitt
# 2026; Ross/Rybak 2018 spindle-gain speed effects).
MOD = dict(ia=0.8, ii=0.6, ib=0.5)

# ------------------------------------------------------------- PF -> MN weights
# Per phase group: {functional_group: weight}. MN weight = primary-group value
# + 0.5 * secondary-group value (biarticular muscles get both).
# Tuning note (2026-09-10, deafferented air run): knee EXTensors receive
# E1+E2+F2 (sum 1.15) vs flexors F1(+E1 0.15) -> co-contration standoff,
# joints pinned at +-2 deg in air. Sharpened: flexors carry F1 strongly,
# F2 only preps extension late in swing, E1 drops its knee-flex/ankle-DF
# co-contraction (stance DF drive fights the E2 push-off).
W_PF_MN: dict[str, dict[str, float]] = {
    "E1": dict(hip_ext=0.45, knee_ext=0.10, hip_abd=0.05, ankle_df=0.10,
               knee_flex=0.05, trunk_ext=0.20),
    "E2": dict(ankle_pf=0.25, hip_ext=0.50, knee_ext=0.15, hip_abd=0.05,
               trunk_ext=0.20),
    "F1": dict(hip_flex=0.75, knee_flex=1.80, ankle_df=0.55, hip_add=0.0),
    "F2": dict(knee_ext=0.0, ankle_df=0.55, hip_flex=0.15, trunk_flex=0.05),
}

# Posture (standing) tonic drive: {functional_group: weight}. Ankle_df gets
# a co-contraction term so the ankle is bidirectional impedance (a one-sided
# PF-held hinge tips backward freely). hip_add/hip_flex/knee_flex tones
# added 2026-09-11: without them the frontal plane was a free pendulum
# (legs splayed +-25 deg of adduction in air).
W_POSTURE = dict(hip_ext=0.22, knee_ext=0.62, ankle_pf=0.30, ankle_df=0.15,
                 hip_abd=0.05, trunk_ext=0.30,
                 hip_add=0.05, hip_flex=0.10, knee_flex=0.08)

# Per-muscle posture overrides (base name, no side suffix): quiet standing
# holds the ankle with the MONOarticular plantarflexors - gastrocnemii also
# flex the knee and will collapse it if co-active at standing levels.
POSTURE_OVERRIDE = {
    "soleus": 0.55, "tib_post": 0.35,
    "med_gas": 0.08, "lat_gas": 0.08,
    "flex_dig": 0.15, "flex_hal": 0.15,
    "per_brev": 0.15, "per_long": 0.15, "per_tert": 0.15,
}

# PF cells whose windows are phase-shifted: (tau multiplier, adapt multiplier).
# E1 fast/short, E2 slower onset + HARD adaptation (its tail was feeding the
# quads/gastrocs through swing - knee never swung), F1 mid, F2 slow/long
PF_SHAPE = dict(E1=(0.6, 0.7), E2=(0.9, 1.7), F1=(0.9, 0.9), F2=(1.6, 1.4))

# ------------------------------------------------------------- balance (standing)
# Supraspinal surrogate: ankle-strategy feedback from pelvis COM offset.
BAL = dict(
    x_ref=None,         # m; None -> runner anchors at the ankle axis
    kx=150.0,           # nA per m of x error
    kv=25.0,            # nA per (m/s) of x velocity
    max_current=6.0,    # nA clamp on the balance injection
    ff=1.0,             # nA initial forward (PF) shove, decays after release
    ff_tau=0.8,         # s, decay constant of the feedforward shove
    fade_with_drive=True,  # balance FB fades out as DRIVE rises (walking)
    ky_lat=400.0,       # nA per m of lateral CoG error (abductor strategy)
    kv_lat=120.0,       # lateral damping term coefficient
    # IMU/vestibular surrogate for trunk pitch (Ben 2026-09-10): PD on the
    # torso up-vector lean (+ = backward limbo lean), drives abdominals
    # when leaning back, erector spinae when leaning forward
    kp_trk=250.0,       # nA per rad of trunk lean
    kd_trk=60.0,        # nA per rad/s of trunk lean rate
    max_trk=10.0,       # nA clamp
    trk_ref=0.06,       # rad, slight forward-lean reference (gait posture)
)

# ------------------------------------------------------------- schedule (seconds)
SCHEDULE = dict(
    stand1=(0.0, 3.0),      # quiet standing, balance FB active
    ramp_up=(3.0, 5.0),     # DRIVE 0 -> 1
    walk=(5.0, 15.0),       # walking
    ramp_down=(15.0, 17.0), # DRIVE 1 -> 0
    stand2=(17.0, 22.0),    # back to quiet standing
)

# In air (no ground contact) the network still rhythms; contact realism comes
# from the load receptor gains. Set True to zero all afferents (pure CPG test).
AFFERENTS_ENABLED = True
