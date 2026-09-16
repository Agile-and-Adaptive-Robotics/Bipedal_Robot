"""Build the spinal network for the Gait2392 MuJoCo model with SNS-Toolbox.

Architecture (one instance per side, plus shared descending cells):

  Descending (network inputs): DRIVE (MLR surrogate: RG frequency + PF gain),
  POSTURE (standing tonic co-activation), BAL_PF / BAL_DF (ankle-strategy
  balance feedback, positive/negative split).

  Rhythm generation (RG), per side: RG-E and RG-F half-center neurons with
  mutual inhibition; each excites its own slow adaptation interneuron
  (ADAP-E/F, tau ~0.9 s) which feeds back inhibition - a relaxation
  oscillator whose frequency rises with DRIVE. Cross-side mutual inhibition
  of the F cells (and weaker between E cells) gives left-right alternation.
  Sensory phase-reset (v5, gains default 0): the runner encodes per-side
  hip-extension (stance-prolonging) and hip-flexion-velocity (swing-
  triggering) signals into input ports HIP_EXT_SIG / HIP_FLEX_SIG; they
  reach the RG through PRESET_E / PRESET_F interneurons (ext: E up / F down;
  flex: F up / E down).

  Pattern formation (PF), per side: PF-E1/E2 driven by RG-E, PF-F1/F2 driven
  by RG-F. The four cells share the same drive but differ in membrane and
  adaptation time constants (PF_SHAPE), giving staggered burst windows within
  each half-cycle - the "similar but phase-shifted" PF family. Conflicting
  groups (E2<->F1, F1<->E1 boundary) inhibit each other.

  Motoneurons + proprioception, per muscle: one MN pool and three afferent
  encoder neurons (Ia spindle-primary ~ tendon velocity; II spindle-secondary
  ~ length; Ib GTO ~ force). Ia excites the homonymous MN and inhibits the
  antagonist MN pools (reciprocal pathway). II excites homonymous MN. Ib
  inhibits the homonymous MN; for stance extensor groups Ib also drives a
  group-level load-sharing interneuron (IB-EXC) that is gated by RG-E and
  excites that group's MNs - stance reflex reversal / intermuscular load
  sharing. Phase gating of afferent gains is applied presynaptically by the
  runner (see params.AFF / MOD and DESIGN.md).

  PF and POSTURE reach each MN through per-functional-group weight tables
  (params.W_PF_MN, params.W_POSTURE): MN weight = primary group value + 0.5 x
  secondary group value. These tables are the back-solving target of
  fit_synapses.py.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from muscle_map import (GROUPS, EXTENSOR_STANCE_GROUPS, MuscleInfo,
                        classify)
from params import (AFF, DT, E_HI, G, MOD, NAP, PF_SHAPE, PHASE_RESET,
                    POSTURE_OVERRIDE, TAU, W_PF_MN, W_POSTURE)

from sns_toolbox import backends as _backends
from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.neurons import NonSpikingNeuronWithPersistentSodiumChannel
from sns_toolbox.networks import Network


class SNS_NumpyFixedTau(_backends.SNS_Numpy):
    """SNS_Numpy with a CONSTANT h-gate time constant (Deng/Animatlab
    semantics; audit row #2's prescribed patch; measured 2026-09-16:
    the stock voltage-dependent tau_h(V) yields 10x-too-fast rhythms).
    Verbatim forward() copy of sns_toolbox 1.5.2's SNS_Numpy.forward;
    the single changed line is marked. If sns_toolbox is ever upgraded,
    re-diff against the new forward body."""

    def forward(self, x=None):
        self.V_last = np.copy(self.V)
        if x is None:
            i_app = 0
        else:
            i_app = np.matmul(self.input_connectivity, x)
        g_syn = np.maximum(0, np.minimum(
            self.g_max_non * ((self.V_last - self.e_lo) /
                              (self.e_hi - self.e_lo)), self.g_max_non))
        if self.spiking:
            self.theta_last = np.copy(self.theta)
            self.g_spike = self.g_spike * (1 - self.time_factor_synapse)
            g_syn += self.g_spike
        i_syn = np.sum(g_syn * self.del_e, axis=1) - \
            self.V_last * np.sum(g_syn, axis=1)
        if self.electrical:
            i_syn += (np.sum(self.g_electrical * self.V_last, axis=1) -
                      self.V_last * np.sum(self.g_electrical, axis=1))
        if self.electrical_rectified:
            mask = np.subtract.outer(self.V_last, self.V_last).T > 0
            masked_g = mask * self.g_rectified
            diag_masked = masked_g + masked_g.T - \
                np.diag(masked_g.diagonal())
            i_syn += np.sum(diag_masked * self.V_last, axis=1) - \
                self.V_last * np.sum(diag_masked, axis=1)
        if self.gated:
            a_inf = 1 / (1 + self.k_a * np.exp(
                self.slope_a * (self.e_a - self.V_last)))
            b_inf = 1 / (1 + self.k_b * np.exp(
                self.slope_b * (self.e_b - self.V_last)))
            c_inf = 1 / (1 + self.k_c * np.exp(
                self.slope_c * (self.e_c - self.V_last)))
            # *** THE PATCH: constant tau_b (was voltage-dependent) ***
            tau_b = self.tau_max_b
            tau_c = self.tau_max_c * c_inf * np.sqrt(
                self.k_c * np.exp(self.slope_c * (self.e_c - self.V_last)))
            self.b_gate_last = np.copy(self.b_gate)
            self.c_gate_last = np.copy(self.c_gate)
            self.b_gate = self.b_gate_last + self.dt * (
                (b_inf - self.b_gate_last) / tau_b)
            self.c_gate = self.c_gate_last + self.dt * (
                (c_inf - self.c_gate_last) / tau_c)
            i_ion = self.g_ion * (a_inf ** self.pow_a) * \
                (self.b_gate ** self.pow_b) * \
                (self.c_gate ** self.pow_c) * (self.e_ion - self.V_last)
            i_gated = np.sum(i_ion, axis=0)
            self.V = self.V_last + self.time_factor_membrane * (
                -self.g_m * (self.V_last - self.V_rest) + self.i_b +
                i_syn + i_app + i_gated)
        else:
            self.V = self.V_last + self.time_factor_membrane * (
                -self.g_m * (self.V_last - self.V_rest) + self.i_b +
                i_syn + i_app)
        if self.spiking:
            self.theta = self.theta_last + self.time_factor_threshold * (
                self.theta_leak * (self.theta_0 - self.theta_last) +
                self.m * (self.V_last - self.V_rest))
            self.spikes = np.sign(np.minimum(0, self.theta - self.V))
            if self.delay:
                self.spike_buffer = np.roll(self.spike_buffer, 1, axis=0)
                self.spike_buffer[0, :] = self.spikes
                self.delayed_spikes[self.spike_rows, self.spike_cols] = \
                    self.spike_buffer[self.buffer_steps, self.buffer_nrns]
                self.g_spike += np.minimum(
                    (-self.delayed_spikes * self.g_increment),
                    (-self.delayed_spikes) *
                    (self.g_max_spike - self.g_spike))
            else:
                self.g_spike += np.minimum(
                    (-self.spikes * self.g_increment),
                    (-self.spikes) * (self.g_max_spike - self.g_spike))
            self.V = ((self.V - self.V_reset) * (self.spikes + 1)) + \
                self.V_reset
            self.theta = np.maximum(
                self.theta_increment,
                self.theta_floor - self.theta) * (-self.spikes) + self.theta
        self.outputs = np.matmul(self.output_voltage_connectivity, self.V)
        if self.spiking:
            self.outputs += np.matmul(self.output_spike_connectivity,
                                      -self.spikes)
        return self.outputs

# Synapse transfer: I = g * clip((V_pre - e_lo)/(e_hi - e_lo), 0, 1) * (V_post - E_rev)
SYN_E_LO, SYN_E_HI = 0.0, E_HI   # presynaptic saturation at E_HI mV
E_REV_EXC = 8.0                  # mV, excitatory reversal (mild, keeps MNs in range)
E_REV_INH = -E_HI                # mV, inhibitory reversal

PF_PHASES = ("E1", "E2", "F1", "F2")

# functional group -> antagonist group(s)
ANTAGONIST = {
    "hip_ext": ("hip_flex",), "hip_flex": ("hip_ext",),
    "hip_abd": ("hip_add",), "hip_add": ("hip_abd",),
    "knee_ext": ("knee_flex",), "knee_flex": ("knee_ext",),
    "ankle_pf": ("ankle_df",), "ankle_df": ("ankle_pf",),
    "trunk_ext": ("trunk_flex",), "trunk_flex": ("trunk_ext",),
}


def _neu(tau: float) -> NonSpikingNeuron:
    """Non-spiking RC neuron with membrane time constant tau (s), G = 1."""
    return NonSpikingNeuron(membrane_capacitance=float(tau),
                            membrane_conductance=1.0,
                            resting_potential=0.0, bias=0.0)


def _syn(g: float, exc: bool) -> NonSpikingSynapse:
    return NonSpikingSynapse(
        max_conductance=float(g),
        reversal_potential=E_REV_EXC if exc else E_REV_INH,
        e_lo=SYN_E_LO, e_hi=SYN_E_HI)


def _group_weight(muscle: MuscleInfo, table: dict[str, float],
                  override: dict[str, float] | None = None) -> float:
    """Primary-group weight + half of any secondary-group weight.

    A per-muscle override (e.g. soleus vs gastroc during standing posture)
    replaces the primary-group weight but still adds the secondary.
    """
    if override and muscle.base in override:
        w = override[muscle.base]
    else:
        w = table.get(muscle.groups[0], 0.0)
    for sec in muscle.groups[1:]:
        w += 0.5 * table.get(sec, 0.0)
    return w


@dataclass
class SpinalNetwork:
    muscles: dict[str, MuscleInfo]                  # actuator name -> info
    sides: tuple[str, ...]
    interleg: bool = True                           # cross-side RG coupling
    net: Network = field(init=False)
    compiled: object = field(init=False, default=None)
    idx: dict[str, int] = field(default_factory=dict)
    inputs: list[str] = field(default_factory=list)  # ordered input-port names
    mn_names: dict[str, str] = field(default_factory=dict)
    aff_names: dict[str, dict[str, str]] = field(default_factory=dict)
    ib_exc_groups: dict[str, tuple[str, ...]] = field(default_factory=dict)
    phase_reset: bool = field(init=False, default=False)
    f1_kneext_inh: bool = field(init=False, default=False)

    # ------------------------------------------------------------------ build
    def __post_init__(self):
        self.net = Network(name="gait2392 spinal")
        n = self.net
        # v5 phase-reset topology is built ONLY when a gain is nonzero:
        # PRESET/PREA phase-reset pathway REMOVED 2026-09-16 (Ben; the
        # per-muscle afferent->central projections replace it). The
        # phase_reset_e/f gains are now INERT (kept so old jsons load).
        # same pattern for the phase-3 swing-knee quad suppression
        self.f1_kneext_inh = bool(G["f1_kneext_inh"] > 0.0
                                  or G["f1_anklepf_inh"] > 0.0)
        # Renshaw recurrent inhibition (Deng Table A6)
        self.renshaw = bool(G["renshaw"] > 0.0)
        # v11 mechanosensory stance feedback (heel/toe contact + stance-Ib
        # prolonger, audit P1a) and IaIN population (audit P1b);
        # ib_e_central also builds the mechanosensor neurons (they ride
        # the same extensor central pathway, Ben 2026-09-16)
        self.stance_fb = bool(G["heel_rge"] > 0.0 or G["toe_rge"] > 0.0
                              or G["ib_rge"] > 0.0
                              or G["ib_e_central"] > 0.0)
        self.ia_in = bool(G["ia_in"] > 0.0)
        self.aff_loops = bool(G["aff_e_rg"] > 0.0 or G["aff_f_rg"] > 0.0
                              or G["aff_e_pf"] > 0.0 or G["aff_f_pf"] > 0.0)

        # ---- descending / balance cells (shared, one each) ----
        for name, tau in (("DRIVE", TAU["descend"]), ("POSTURE", TAU["descend"]),
                          ("BAL_PF", TAU["descend"]), ("BAL_DF", TAU["descend"]),
                          ("BAL_TRK_EXT", TAU["descend"]),
                          ("BAL_TRK_FLX", TAU["descend"]),
                          ("BAL_LAT_R", TAU["descend"]),
                          ("BAL_LAT_L", TAU["descend"])):
            n.add_neuron(_neu(tau), name=name)
            self.idx[name] = len(self.idx)
            n.add_input(name)
            self.inputs.append(name)

        # ---- per-side circuitry ----
        for side in self.sides:
            self._build_rg(n, side)
            self._build_pf(n, side)

        # ---- muscles: create all neurons first, then wire (reciprocal
        # Ia inhibition references antagonist MNs across muscles) ----
        self._order = {a: i for i, a in enumerate(self.muscles)}
        for act, mi in self.muscles.items():
            self._add_muscle_neurons(n, act, mi)
        for act, mi in self.muscles.items():
            self._wire_muscle(n, act, mi)
        self._wire_balance(n)

        # ---- cross-side coordination, laminated through commissural INs
        # (Shinohara 2025 wiring): RG_F -> exc -> c1 IN -> INHIBIT ->
        # contralateral RG_F (flexor antiphase lock); RG_E -> exc -> V3
        # IN -> EXCITE -> contralateral RG_E (extensor synchronization).
        # interleg=False removes it entirely: independent left/right
        # rhythm generators (deafferented air-stepping preparation).
        if self.interleg:
            for a, b in (("r", "l"), ("l", "r")):
                cf, ce = f"CIN_F_{a}", f"CIN_E_{a}"
                self._add(cf, TAU["rg"], n)
                self._add(ce, TAU["rg"], n)
                n.add_connection(_syn(G["rg_mutual_inh"], exc=True),
                                 f"RG_F_{a}", cf)
                n.add_connection(_syn(G["rg_mutual_inh"], exc=False),
                                 cf, f"RG_F_{b}")
                n.add_connection(_syn(0.5 * G["rg_mutual_inh"], exc=True),
                                 f"RG_E_{a}", ce)
                n.add_connection(_syn(0.5 * G["rg_mutual_inh"], exc=True),
                                 ce, f"RG_E_{b}")
                # audit P2a / Rybak 2025 crossed-extensor reinforcement:
                # the V3 commissural also excites the CONTRALATERAL
                # extensor MN group INs (weight support on the partner
                # leg). Conditional on the gain.
                if G["v3_to_ibexc"] > 0.0:
                    for grp in self.ib_exc_groups.get(b, ()):
                        gname = f"IBEXC_{grp}_{b}"
                        if gname in self.idx:
                            n.add_connection(
                                _syn(G["v3_to_ibexc"], exc=True), ce, gname)

    # ------------------------------------------------------------------ parts
    def _add(self, name: str, tau: float, n: Network):
        n.add_neuron(_neu(tau), name=name)
        self.idx[name] = len(self.idx)

    def _build_rg(self, n: Network, side: str):
        rg_e, rg_f = f"RG_E_{side}", f"RG_F_{side}"
        ine, inf = f"InE_{side}", f"InF_{side}"
        # Deng 2022 / Shinohara 2025 / Rybak 2024-25 conditional bursters:
        # RG half-centers are PERSISTENT-SODIUM neurons (intrinsic burst
        # termination via the slow h-gate, FIXED tau per Animatlab
        # semantics; see SNS_NumpyFixedTau above + params.NAP). The ADAP
        # circuit substitute is RETIRED (2026-09-16).
        nap = NonSpikingNeuronWithPersistentSodiumChannel(
            membrane_capacitance=TAU["rg"], membrane_conductance=1.0,
            resting_potential=0.0, bias=0.0,
            g_ion=np.array([NAP["g_ion"]]),
            e_ion=np.array([NAP["e_ion"]]),
            k_m=np.array([NAP["k_m"]]),
            slope_m=np.array([NAP["slope_m"]]),
            e_m=np.array([NAP["e_m"]]),
            k_h=np.array([NAP["k_h"]]),
            slope_h=np.array([NAP["slope_h"]]),
            e_h=np.array([NAP["e_h"]]),
            tau_max_h=np.array([TAU["rg_nap_h"]]))
        n.add_neuron(nap, name=rg_e)
        self.idx[rg_e] = len(self.idx)
        n.add_neuron(nap, name=rg_f)
        self.idx[rg_f] = len(self.idx)
        for name in (ine, inf):
            self._add(name, TAU["rg"], n)

        # IN-laminated mutual inhibition (Shevtsova/Deng A6 architecture):
        # RG-E excites InE, InE inhibits RG-F; RG-F excites InF, InF
        # inhibits RG-E. NO direct inhibitory RG<->RG synapses.
        n.add_connection(_syn(G["rg_mutual_inh"], exc=True), rg_e, ine)
        n.add_connection(_syn(G["rg_mutual_inh"], exc=False), ine, rg_f)
        n.add_connection(_syn(G["rg_mutual_inh"], exc=True), rg_f, inf)
        n.add_connection(_syn(G["rg_mutual_inh"], exc=False), inf, rg_e)

        # weak MUTUAL EXCITATION between the half-centers (Deng 2022 G_W;
        # hypothesized counterpart of the V1/V2b-ablation synchronous
        # rhythm). Raises the inhibited neuron's equilibrium -> escape
        # mode; neuromodulating this conductance raises frequency
        # (Deng 2022 Fig 3). Direct RG-E<->RG-F edges; conditional on the
        # gain (zero-g synapses still perturb summation order - v5 lesson).
        if G["rg_weak_exc"] > 0.0:
            n.add_connection(_syn(G["rg_weak_exc"], exc=True), rg_e, rg_f)
            n.add_connection(_syn(G["rg_weak_exc"], exc=True), rg_f, rg_e)
        # descending drive raises frequency; stance-biased split sets duty
        n.add_connection(_syn(G["descend_to_rg_e"], exc=True), "DRIVE", rg_e)
        n.add_connection(_syn(G["descend_to_rg_f"], exc=True), "DRIVE", rg_f)
        # posture tonic bias keeps the E cell (load-bearing side) ready
        n.add_connection(_syn(G["posture_to_rg_e"], exc=True), "POSTURE", rg_e)

        # ---- v5/v10 PRESET/PREA hip phase-reset pathway REMOVED
        # (2026-09-16, Ben circled it in the panel figure): it was the
        # hip-selective special case of sensory phase reset; the general
        # literature mechanism is now wired directly - per-muscle
        # afferent -> central projections (ia_f_central / ii_f_central /
        # ib_e_central / ia_f_contra_f; Rybak 2025 SF-E1/SF-E2, Perreault
        # 2011 gr-II reset, Gossard 1994 Ib). Resets terminate INSIDE the
        # RG layer per Rybak 2015, as before.

        # v11 mechanosensory stance feedback (audit P1a; conditional
        # on the gains, defaults 0 = absent = v10-identical):
        #   HEEL_c -> HEEL_IN -> RG-E exc + RG-F inh (S2W trigger,
        #   Conway/Hultborn 1987 reset-to-extension) + the EXTENSOR
        #   central pathway (PF_E1/E2 + InE) when ib_e_central > 0
        #   (Ben 2026-09-16: foot mechanosensory feedback follows the
        #   same route as extensor Ib feedback)
        #   TOE_c  -> TOE_IN  -> RG-E exc (late-stance prolongation) +
        #   the same extensor central pathway
        #   LBIN   -> RG-E exc (stance-Ib group IN per Dominguez 2020:
        #   "INs belong to the rhythm-generating layer"; Gossard 1994 /
        #   Pearson 1998 stance-duration regulation). LBIN receives the
        #   stance-group IBEXC outputs (RG-E-gated Ib) in _wire_muscle
        #   and the contact-load current from the runner.
        if self.stance_fb or self.aff_loops:
            heel_in = f"HEEL_{side}"
            toe_in = f"TOE_{side}"
            lbin = f"LBIN_{side}"
            self._add(heel_in, TAU["preset"], n)
            self._add(toe_in, TAU["preset"], n)
            self._add(lbin, TAU["ib_exc"], n)
            n.add_input(heel_in)
            self.inputs.append("HEEL_c_" + side)
            n.add_input(toe_in)
            self.inputs.append("TOE_c_" + side)
            n.add_input(lbin)
            self.inputs.append("LOAD_c_" + side)
            n.add_connection(_syn(G["heel_rge"], exc=True), heel_in, rg_e)
            n.add_connection(_syn(G["heel_rge"] * PHASE_RESET.get("inh", 1.0),
                                  exc=False), heel_in, rg_f)
            n.add_connection(_syn(G["toe_rge"], exc=True), toe_in, rg_e)
            n.add_connection(_syn(G["ib_rge"], exc=True), lbin, rg_e)
            # mechanosensors ride the extensor central pathway (same
            # gain as extensor Ib -> central, per Ben's figure reading)
            if G["ib_e_central"] > 0.0:
                for src in (heel_in, toe_in):
                    n.add_connection(_syn(G["ib_e_central"], exc=True),
                                     src, f"PF_E1_{side}")
                    n.add_connection(_syn(G["ib_e_central"], exc=True),
                                     src, f"PF_E2_{side}")
                    n.add_connection(_syn(G["ib_e_central"], exc=True),
                                     src, ine)

        # ---- v11b semi-closed sensory loops (Shevtsova central principle):
        # the active phase's afferents excite that phase's PF and RG cells
        # through relay INs, creating a three-layer positive feedback loop:
        #   muscle -> afferent -> PF -> RG -> MN -> muscle
        # AFF_E receives extensor-side afferent drive (force/length);
        # AFF_F receives flexor-side afferent drive (velocity/length).
        # Both are input ports fed by the runner from the afferent
        # encoder signals. Gains default 0 = absent (conditional topology).
        if self.aff_loops:
            aff_e = f"AFF_E_{side}"
            aff_f = f"AFF_F_{side}"
            self._add(aff_e, TAU["afferent"], n)
            self._add(aff_f, TAU["afferent"], n)
            n.add_input(aff_e)
            self.inputs.append("AFF_E_" + side)
            n.add_input(aff_f)
            self.inputs.append("AFF_F_" + side)
            n.add_connection(_syn(G["aff_e_rg"], exc=True), aff_e, rg_e)
            n.add_connection(_syn(G["aff_f_rg"], exc=True), aff_f, rg_f)

    def _build_pf(self, n: Network, side: str):
        rg_of = {"E1": f"RG_E_{side}", "E2": f"RG_E_{side}",
                 "F1": f"RG_F_{side}", "F2": f"RG_F_{side}"}
        for phase in PF_PHASES:
            tau_m, _ = PF_SHAPE[phase]
            pf = f"PF_{phase}_{side}"
            self._add(pf, TAU["pf"] * tau_m, n)

            n.add_connection(_syn(G["rg_to_pf"], exc=True), rg_of[phase], pf)
            # (DRIVE->PF weak tonic edge REMOVED 2026-09-16: audit #23
            # EXTRA-NOTSUPPORTED - Deng/Shevtsova drive PF from RG only)

        # IN-laminated cross-reciprocal inhibition (Shevtsova/Deng A6):
        # PF-E cells excite PF_IN_E, PF_IN_E inhibits PF-F cells;
        # PF-F cells excite PF_IN_F, PF_IN_F inhibits PF-E cells.
        # NO direct PF↔PF synapses.
        pf_in_e = f"PF_IN_E_{side}"
        pf_in_f = f"PF_IN_F_{side}"
        self._add(pf_in_e, TAU["pf"], n)
        self._add(pf_in_f, TAU["pf"], n)
        for ph in ("E1", "E2"):
            n.add_connection(_syn(G["pf_recip_inh"], exc=True),
                             f"PF_{ph}_{side}", pf_in_e)
        for ph in ("F1", "F2"):
            n.add_connection(_syn(G["pf_recip_inh"], exc=True),
                             f"PF_{ph}_{side}", pf_in_f)
        for ph in ("F1", "F2"):
            n.add_connection(_syn(G["pf_recip_inh"], exc=False),
                             pf_in_e, f"PF_{ph}_{side}")
        for ph in ("E1", "E2"):
            n.add_connection(_syn(G["pf_recip_inh"], exc=False),
                             pf_in_f, f"PF_{ph}_{side}")

        # v11b Shevtsova semi-closed loops: AFF relay → PF cells
        # (afferent → PF excitation, completing the three-layer loop)
        if self.aff_loops:
            aff_e = f"AFF_E_{side}"
            aff_f = f"AFF_F_{side}"
            for ph in ("E1", "E2"):
                n.add_connection(_syn(G["aff_e_pf"], exc=True),
                                 aff_e, f"PF_{ph}_{side}")
            for ph in ("F1", "F2"):
                n.add_connection(_syn(G["aff_f_pf"], exc=True),
                                 aff_f, f"PF_{ph}_{side}")

    def _add_muscle_neurons(self, n: Network, act: str, mi: MuscleInfo):
        mn, ia, ii, ib = (f"MN_{act}", f"Ia_{act}", f"II_{act}", f"Ib_{act}")
        self.mn_names[act] = mn
        self.aff_names[act] = {"Ia": ia, "II": ii, "Ib": ib}
        for name, tau in ((mn, TAU["mn"] * (1.0 + 0.5 * mi.biarticular)),
                          (ia, TAU["afferent"]), (ii, TAU["afferent"]),
                          (ib, 2.0 * TAU["afferent"])):
            self._add(name, tau, n)
        # Renshaw cell per pool (Deng Table A6: recurrent inhibition +
        # RC<->RC mutual inhibition). Wired in _wire_muscle.
        if self.renshaw:
            self._add(f"RC_{act}", TAU["mn"], n)
        # MNs take one external current: the solved/standing posture bias
        n.add_input(mn)
        self.inputs.append("POST_" + act)
        for port, name in (("Ia", ia), ("II", ii), ("Ib", ib)):
            n.add_input(name)
            self.inputs.append(port + "_" + act)

    def _wire_muscle(self, n: Network, act: str, mi: MuscleInfo):
        mn, ia, ii, ib = (f"MN_{act}", f"Ia_{act}", f"II_{act}", f"Ib_{act}")

        # ---- PF -> MN (per phase groups covering this muscle) + POSTURE ----
        for phase in PF_PHASES:
            w = _group_weight(mi, W_PF_MN[phase])
            if w > 0.0:
                for side in self.sides:
                    if mi.side == side:
                        n.add_connection(_syn(G["pf_to_mn"] * w, exc=True),
                                         f"PF_{phase}_{side}", mn)
        w_post = _group_weight(mi, W_POSTURE, POSTURE_OVERRIDE)
        if w_post > 0.0:
            n.add_connection(_syn(G["posture_to_mn"] * w_post, exc=True),
                             "POSTURE", mn)

        # ---- proprioceptive pathways ----
        n.add_connection(_syn(G["ia_to_mn"], exc=True), ia, mn)
        n.add_connection(_syn(G["ii_to_mn"], exc=True), ii, mn)
        n.add_connection(_syn(G["ib_to_mn_inh"], exc=False), ib, mn)

        # ---- per-muscle afferent -> CENTRAL feedback (Deng 2022 /
        # Shinohara 2025 wiring, Ben 2026-09-16): force/length feedback
        # shifts the RG V-nullcline and triggers phase transitions
        # (Shinohara sec 4.2: force feedback excites the extensor center,
        # length feedback the flexor center). Extensor-family muscles'
        # Ib afferents project EXCITATORY to the ipsilateral E-centers
        # (PF_E1/E2, RG_E, InE); flexor-family muscles' Ia and II
        # afferents project to the F-centers (PF_F1/F2, RG_F, InF).
        # Foot mechanosensors share the extensor pathway (see _build_rg).
        # Conditional on the gains (default 0 = topology absent).
        grp = mi.groups[0]
        if G["ib_e_central"] > 0.0 and grp in EXTENSOR_STANCE_GROUPS:
            for tgt in (f"PF_E1_{mi.side}", f"PF_E2_{mi.side}",
                        f"RG_E_{mi.side}", f"InE_{mi.side}"):
                n.add_connection(_syn(G["ib_e_central"], exc=True),
                                 ib, tgt)
        if grp in ("hip_flex", "knee_flex", "ankle_df", "hip_add",
                   "trunk_flex"):
            for tgt in (f"PF_F1_{mi.side}", f"PF_F2_{mi.side}",
                        f"RG_F_{mi.side}", f"InF_{mi.side}"):
                if G["ia_f_central"] > 0.0:
                    n.add_connection(_syn(G["ia_f_central"], exc=True),
                                     ia, tgt)
                if G["ii_f_central"] > 0.0:
                    n.add_connection(_syn(G["ii_f_central"], exc=True),
                                     ii, tgt)
            # Rybak 2025 SF-E1 contralateral half: hip-flexor stretch
            # afferents also inhibit the CONTRALATERAL F half-center
            # (promotes the E->F transition + interleg coordination)
            if G["ia_f_contra_f"] > 0.0 and self.interleg:
                contra = "l" if mi.side == "r" else "r"
                n.add_connection(_syn(G["ia_f_contra_f"], exc=False),
                                 ia, f"RG_F_{contra}")
        # extensor II same-group excitation (Ben 2026-09-16: type II is
        # same-group excitatory, like the flexor-side Ia/II routing)
        if G["ii_e_central"] > 0.0 and grp in EXTENSOR_STANCE_GROUPS:
            for tgt in (f"PF_E1_{mi.side}", f"PF_E2_{mi.side}",
                        f"RG_E_{mi.side}", f"InE_{mi.side}"):
                n.add_connection(_syn(G["ii_e_central"], exc=True),
                                 ii, tgt)

        # ---- Renshaw recurrent inhibition (Deng Table A6: MN->RC 0.5
        # exc; RC->MN 0.5 inh; RC<->RC 0.5 inh), gain-scaled by
        # G["renshaw"]; topology present only when that gain > 0.
        if self.renshaw:
            rc = f"RC_{act}"
            g_r = G["renshaw"]
            n.add_connection(_syn(1.0, exc=True), mn, rc)
            n.add_connection(_syn(g_r, exc=False), rc, mn)
            # RC<->RC mutual inhibition between DIFFERENT pools (Deng A6;
            # Hultborn's cat data): each ordered pair wired once - both
            # directions emerge from the two muscles' wiring loops. No
            # self-synapse (an autapse is not standard RC connectivity).
            for act2, mi2 in self.muscles.items():
                if mi2.side == mi.side and act2 != act \
                        and f"RC_{act2}" in self.idx:
                    n.add_connection(_syn(g_r, exc=False), rc,
                                     f"RC_{act2}")

        # Ia reciprocal inhibition of antagonist MN pools.
        # v11 (audit P1b): when G["ia_in"] > 0 the reciprocal pathway is
        # routed through an IaIN population (Ia -> IaIN -> antagonist MN)
        # with PF_F1 phase-gating (Deng A6) and RC->IaIN recurrent
        # disinhibition (Hultborn 1971) when Renshaw is on; the direct
        # edge is used only when ia_in == 0 (v10 behavior).
        if self.ia_in:
            iain = f"IaIN_{act}"
            if iain not in self.idx:
                self._add(iain, TAU["afferent"], n)
                # afferent drive: the Ia afferent EXCITES its IaIN (the
                # afferent leg of "Ia -> IaIN -> antagonist MN"; same
                # conductance as the homonymous Ia->MN arc, Deng A6)
                n.add_connection(_syn(G["ia_to_mn"], exc=True), ia, iain)
                # phase gate: PF_F1 excites the IaIN (Deng A6 PF->IaIN 0.5)
                n.add_connection(_syn(0.5, exc=True), f"PF_F1_{mi.side}",
                                 iain)
                # recurrent disinhibition (RC -> IaIN inh)
                if self.renshaw and f"RC_{act}" in self.idx:
                    n.add_connection(_syn(G["renshaw"], exc=False),
                                     f"RC_{act}", iain)
        for ant in ANTAGONIST.get(mi.groups[0], ()):
            for act2, mi2 in self.muscles.items():
                if mi2.side == mi.side and mi2.groups[0] == ant:
                    if self.ia_in:
                        n.add_connection(
                            _syn(G["ia_to_antagonist"], exc=False),
                            f"IaIN_{act}", f"MN_{act2}")
                    else:
                        n.add_connection(
                            _syn(G["ia_to_antagonist"], exc=False),
                            ia, f"MN_{act2}")

        # ---- stance load sharing (extensor groups only) ----
        if mi.groups[0] in EXTENSOR_STANCE_GROUPS:
            grp = f"IBEXC_{mi.groups[0]}"
            gname = f"{grp}_{mi.side}"
            if gname not in self.idx:
                self._add(gname, TAU["ib_exc"], n)
                # gated by the stance (E) half-center
                n.add_connection(_syn(1.0, exc=True), f"RG_E_{mi.side}", gname)
            n.add_connection(_syn(G["ib_group_exc"], exc=True), ib, gname)
            n.add_connection(_syn(G["ib_exc_to_mn"], exc=True), gname, mn)
            # v11 P1a: the stance-group Ib IN also drives the RG-layer
            # load IN (LBIN) which excites RG-E — load prolongs stance
            if self.stance_fb:
                lb = f"LBIN_{mi.side}"
                if lb not in self.idx:
                    self._add(lb, TAU["ib_exc"], n)
                n.add_connection(_syn(0.5, exc=True), gname, lb)
            if mi.groups[0] not in self.ib_exc_groups.get(mi.side, ()):
                self.ib_exc_groups[mi.side] = self.ib_exc_groups.get(mi.side, ()) + (mi.groups[0],)

        # ---- phase-3 swing-knee quad suppression (gain default 0 =
        # topology absent): F1 (the swing flexor window) -> KINH
        # inhibitory interneuron -> every primary knee_ext MN of the side.
        # The suppression is phase-gated BY F1 itself - it can only act
        # during the swing window, releasing the quads so the knee can
        # flex (v4 diagnosis lever #2: swing knee extension-dominant).
        if self.f1_kneext_inh and mi.groups[0] == "knee_ext":
            kname = f"KINH_{mi.side}"
            if kname not in self.idx:
                self._add(kname, TAU["ib_exc"], n)
                n.add_connection(_syn(1.5, exc=True), f"PF_F1_{mi.side}",
                                 kname)
            n.add_connection(_syn(G["f1_kneext_inh"], exc=False), kname, mn)
        # v6b: same swing-gated suppression onto ankle PF pools (reuses
        # the KINH IN; separate gain)
        if self.f1_kneext_inh and mi.groups[0] == "ankle_pf":
            kname = f"KINH_{mi.side}"
            if kname not in self.idx:
                self._add(kname, TAU["ib_exc"], n)
                n.add_connection(_syn(1.5, exc=True), f"PF_F1_{mi.side}",
                                 kname)
            n.add_connection(_syn(G["f1_anklepf_inh"], exc=False), kname, mn)

    def _wire_balance(self, n: Network):
        """Balance inputs reach ankle + hip MNs (ankle + hip strategy).

        BAL_PF active = body swaying backward -> plantarflexion push + hip
        flexion pull the COM forward; BAL_DF = the opposite.
        BAL_LAT_R / BAL_LAT_L: frontal-plane strategy - lateral COM error
        drives the STANCE-side hip abductors (glut_med/min), pulling the
        CoG toward the stance foot (Trendelenburg mechanics; opensim +z =
        mujoco -y per Ben).
        Called once after all MNs exist.
        """
        for act, mi in self.muscles.items():
            g = mi.groups[0]
            if g in ("ankle_pf", "hip_flex"):
                n.add_connection(_syn(0.5 * G["posture_to_mn"] * 0.5, exc=True),
                                 "BAL_PF", f"MN_{act}")
            elif g in ("ankle_df", "hip_ext"):
                n.add_connection(_syn(0.5 * G["posture_to_mn"] * 0.5, exc=True),
                                 "BAL_DF", f"MN_{act}")
            elif g == "hip_abd":
                n.add_connection(_syn(G["bal_lat_to_abd"], exc=True),
                                 f"BAL_LAT_{mi.side.upper()}", f"MN_{act}")
            elif g == "trunk_ext":
                # IMU/vestibular surrogate -> erector spinae (Ben's trunk
                # stability plan; the trunk was doing the limbo without it)
                n.add_connection(_syn(G["bal_trunk"], exc=True),
                                 "BAL_TRK_EXT", f"MN_{act}")
            elif g == "trunk_flex":
                n.add_connection(_syn(G["bal_trunk"], exc=True),
                                 "BAL_TRK_FLX", f"MN_{act}")

    # ------------------------------------------------------------------ run
    def compile(self, dt: float = DT):
        self.compiled = self.net.compile(dt=dt, backend="numpy")
        # swap in the FIXED-tau_h stepper for the RG persistent-Na
        # h-gates (Deng semantics; see SNS_NumpyFixedTau docstring)
        self.compiled.__class__ = SNS_NumpyFixedTau
        assert self.compiled.V.shape[0] == len(self.idx), \
            "neuron index bookkeeping mismatch"
        return self.compiled

    def step(self, input_currents: np.ndarray) -> np.ndarray:
        """Advance one dt; input vector must match self.inputs order.

        Returns the full membrane-potential vector (indexed by self.idx).
        """
        self.compiled.forward(list(input_currents))
        return self.compiled.V

    # ------------------------------------------------------------------ helpers
    def input_index(self, port: str) -> int:
        return self.inputs.index(port)

    def make_inputs(self) -> np.ndarray:
        return np.zeros(len(self.inputs))


def build(model_actuators: list[str], dt: float = DT,
          interleg: bool = True) -> SpinalNetwork:
    """Classify actuators and build + compile the spinal network.

    interleg=False removes all cross-side RG coupling (independent
    half-centers per leg - the deafferented air-stepping preparation).
    """
    muscles: dict[str, MuscleInfo] = {}
    for act in model_actuators:
        mi = classify(act)
        if mi is None:
            raise ValueError(f"actuator {act!r} not in muscle_map")
        muscles[act] = mi
    sides = tuple(sorted({mi.side for mi in muscles.values()}))
    net = SpinalNetwork(muscles=muscles, sides=sides, interleg=interleg)
    net.compile(dt=dt)
    return net
