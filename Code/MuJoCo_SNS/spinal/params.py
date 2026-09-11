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
    rg=0.05,        # rhythm-generator half-center cells (50 ms)
    rg_adapt=1.9,   # slow adaptation interneuron (burst termination; slow
                     # for slow air-stepping rhythm - frequency knob that
                     # does NOT shrink burst amplitudes, unlike low DRIVE)
    pf=0.08,        # pattern-formation cells
    pf_adapt=0.5,   # PF burst self-adaptation
    ib_exc=0.05,    # group Ib load-sharing interneuron
    descend=0.10,   # descending drive smoothing
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
    "E1": dict(hip_ext=0.45, knee_ext=0.10, hip_abd=0.35, ankle_df=0.10,
               knee_flex=0.05, trunk_ext=0.20),
    "E2": dict(ankle_pf=0.35, hip_ext=0.50, knee_ext=0.15, hip_abd=0.30,
               trunk_ext=0.20),
    "F1": dict(hip_flex=0.45, knee_flex=1.80, ankle_df=0.45, hip_add=0.10),
    "F2": dict(knee_ext=0.0, ankle_df=0.45, hip_flex=0.05, trunk_flex=0.05),
}

# Posture (standing) tonic drive: {functional_group: weight}. Ankle_df gets
# a co-contraction term so the ankle is bidirectional impedance (a one-sided
# PF-held hinge tips backward freely).
W_POSTURE = dict(hip_ext=0.22, knee_ext=0.62, ankle_pf=0.30, ankle_df=0.15,
                 hip_abd=0.35, trunk_ext=0.30)

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
