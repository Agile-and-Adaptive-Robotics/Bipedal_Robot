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

import json
from dataclasses import dataclass, field
from pathlib import Path

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

# ---- 2026-09-18 joint-layer PF (G["joint_pf"] > 0; default 0 = the
# phase-cell PF is built instead).  Three joint/functional half-center
# pairs replace the four phase-window cells; PF->MN weights come from
# the structured T1 fit (fsa_jointlayers.py -> joint_pf_weights.json,
# held-out centered VAF 0.747/0.713 vs 0.934/0.925 unconstrained).
JPF_HCS = ("HIP-E", "HIP-F", "KNEE-E", "KNEE-F", "ANK-E", "ANK-F")
JPF_GROUP2HC = {"hip_ext": ("HIP-E",), "hip_abd": ("HIP-E",),
                "hip_add": ("HIP-E", "HIP-F"), "hip_flex": ("HIP-F",),
                "knee_ext": ("KNEE-E",), "knee_flex": ("KNEE-F",),
                "ankle_pf": ("ANK-E",), "ankle_df": ("ANK-F",)}
JPF_CROSS = {"rect_fem": ("HIP-F",), "semimem": ("HIP-E",),
             "semiten": ("HIP-E",), "bifemsh": ("HIP-E",),
             "gas_med": ("KNEE-F",), "gas_lat": ("KNEE-F",),
             "grac": ("HIP-F",), "sart": ("KNEE-F",)}
_JPF_W = None


def _jpf_weights():
    """Muscle-level PF weights per joint-layer HC (fitted where the SO
    had signal, anatomical 1.0 fallback otherwise)."""
    global _JPF_W
    if _JPF_W is None:
        try:
            data = json.loads(
                (Path(__file__).parent / "joint_pf_weights.json")
                .read_text(encoding="utf-8"))
            _JPF_W = data["w"]
        except (FileNotFoundError, KeyError, ValueError):
            _JPF_W = {}
    return _JPF_W

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
        # the same extensor central pathway, Ben 2026-09-16). The
        # 2026-09-24 per-PF-layer contact variant keys also build the
        # HEEL/TOE INs (heel_pf_layer / toe_df_inh / heel_in_f_exc).
        self.stance_fb = bool(G["heel_rge"] > 0.0 or G["toe_rge"] > 0.0
                              or G["ib_rge"] > 0.0
                              or G["ib_e_central"] > 0.0
                              or G["heel_pf_layer"] > 0.0
                              or G["toe_df_inh"] > 0.0
                              or G["heel_in_f_exc"] > 0.0)
        self.ia_in = bool(G["ia_in"] > 0.0)
        self.aff_loops = bool(G["aff_e_rg"] > 0.0 or G["aff_f_rg"] > 0.0
                              or G["aff_e_pf"] > 0.0 or G["aff_f_pf"] > 0.0)
        # 2026-09-23 goal2 standing-balance stage (Ben's request; SCONE
        # Tutorial-3a analog): VEST_{r,l} vestibular-analog cells, built
        # ONLY when a vest gain > 0 (zero-gain synapses alone change BLAS
        # summation order - the v5 lesson; defaults 0 = network identical)
        self.vest = bool(G["vest_ext"] > 0.0 or G["vest_flex_inh"] > 0.0)
        # 2026-09-18 joint-layer PF (T1/T3 experiment; 0 = phase cells)
        self.joint_pf = bool(G.get("joint_pf", 0.0) > 0.0)
        if self.joint_pf:
            _jpf_weights()

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

        # ---- goal2 vestibular-analog cells (one per side; conditional
        # topology, defaults 0 = absent). The runner feeds each a current
        # rectified from the pelvis-tilt deviation + rate (otolith/canal
        # analog; SCONE T3a vestibular BodyPointReflex: torso-point PD,
        # 0.1 s delay ~ this 0.1 s membrane tau). Same cell family as the
        # BAL_* descending cells (brainstem surrogate, not spinal IN).
        if self.vest:
            for side in self.sides:
                vname = f"VEST_{side}"
                n.add_neuron(_neu(TAU["descend"]), name=vname)
                self.idx[vname] = len(self.idx)
                n.add_input(vname)
                self.inputs.append("VEST_c_" + side)

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
        # IN -> EXCITE -> contralateral InE (the RG-E interneuron, per
        # Ben's figure reading - NOT RG_E itself: V3->contra-RG_E direct
        # excitation forms a bilateral E<->E positive-feedback loop that
        # latches both NaP plateaus; the InE target supports the partner
        # leg's extensor side indirectly by suppressing its flexor).
        # interleg=False removes it entirely: independent left/right
        # rhythm generators (deafferented air-stepping preparation).
        if self.interleg:
            for a, b in (("r", "l"), ("l", "r")):
                # Names are indexed by SOURCE side: these are four distinct
                # directional cells (C1_r->l, C1_l->r, V3_r->l, V3_l->r),
                # never one shared bidirectional relay.
                cf, ce = f"CIN_F_{a}", f"CIN_E_{a}"
                self._add(cf, TAU["rg"], n)
                self._add(ce, TAU["rg"], n)
                n.add_connection(_syn(G["c1_gain"] * G["rg_mutual_inh"],
                                      exc=True), f"RG_F_{a}", cf)
                n.add_connection(_syn(G["c1_gain"] * G["rg_mutual_inh"],
                                      exc=False), cf, f"RG_F_{b}")
                if G["v3_gain"] > 0.0:
                    n.add_connection(
                        _syn(G["v3_gain"] * G["rg_mutual_inh"], exc=True),
                        f"RG_E_{a}", ce)
                    n.add_connection(
                        _syn(G["v3_gain"] * G["rg_mutual_inh"], exc=True),
                        ce, f"InE_{b}")
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
            # 2026-09-21 FULL RULES: reach the RG through the laminated
            # IN layer (heel -> InE excitation strengthens E's
            # suppression of F; heel -> InF inhibition releases RG-E),
            # not direct heel-IN -> half-center edges
            if G["full_rules"] > 0.0:
                n.add_connection(_syn(G["heel_rge"], exc=True),
                                 heel_in, ine)
                n.add_connection(_syn(G["heel_rge"], exc=False),
                                 heel_in, inf)
                n.add_connection(_syn(G["toe_rge"], exc=True),
                                 toe_in, ine)
            else:
                n.add_connection(_syn(G["heel_rge"], exc=True),
                                 heel_in, rg_e)
                n.add_connection(_syn(G["heel_rge"] * PHASE_RESET.get(
                    "inh", 1.0), exc=False), heel_in, rg_f)
                n.add_connection(_syn(G["toe_rge"], exc=True),
                                 toe_in, rg_e)
            n.add_connection(_syn(G["ib_rge"], exc=True), lbin, rg_e)
            # 2026-09-24 per-PF-layer contact variant (Ben's drawing,
            # heel_in_f_exc): heel IN -> InF EXCITATORY. The full_rules
            # branch above wires heel -> InF INHIBITORY; Ben's updated
            # drawing shows both InE and InF excited (g 0.5) - this
            # gain-gated edge adds the excitatory variant, coexisting.
            if G["heel_in_f_exc"] > 0.0:
                n.add_connection(_syn(G["heel_in_f_exc"], exc=True),
                                 heel_in, inf)
            # NOTE: the heel/toe -> PF_E/InE central-pathway extensions
            # are wired in _build_pf (the PF cells are created there,
            # AFTER this function runs - wiring them here was the
            # stage-2 crash of 2026-09-16 18:58, 'Population not found
            # by name PF_E1_l')

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
        if self.joint_pf:
            self._build_pf_layers(n, side)
            return
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

        # mechanosensors ride the extensor central pathway (same gain as
        # extensor Ib -> central, per Ben's figure reading 2026-09-16):
        # HEEL/TOE -> PF_E1/E2 + InE. Wired HERE because the PF cells are
        # created in this function (the 2026-09-16 stage-2 crash was this
        # block living in _build_rg before the PF populations existed).
        if (self.stance_fb or self.aff_loops) and G["ib_e_central"] > 0.0:
            for src in (f"HEEL_{side}", f"TOE_{side}"):
                n.add_connection(_syn(G["ib_e_central"], exc=True),
                                 src, f"PF_E1_{side}")
                n.add_connection(_syn(G["ib_e_central"], exc=True),
                                 src, f"PF_E2_{side}")
                n.add_connection(_syn(G["ib_e_central"], exc=True),
                                 src, f"InE_{side}")

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

    def _build_pf_layers(self, n: Network, side: str):
        """Joint-layer PF (G["joint_pf"] > 0, 2026-09-18): three
        joint/functional half-center pairs (HIP/KNEE/ANK x E/F) instead
        of four phase windows.  Same RG drive, IN lamination, and
        central-pathway attachments as the phase cells; MN weights come
        from the structured T1 fit (joint_pf_weights.json)."""
        for hc in JPF_HCS:
            tau_m, _ = PF_SHAPE["E1" if hc.endswith("E") else "F1"]
            pf = f"PF_{hc}_{side}"
            self._add(pf, TAU["pf"] * tau_m, n)
            n.add_connection(_syn(G["rg_to_pf"], exc=True),
                             f"RG_{'E' if hc.endswith('E') else 'F'}_{side}",
                             pf)
        pf_in_e = f"PF_IN_E_{side}"
        pf_in_f = f"PF_IN_F_{side}"
        self._add(pf_in_e, TAU["pf"], n)
        self._add(pf_in_f, TAU["pf"], n)
        for hc in JPF_HCS:
            n.add_connection(_syn(G["pf_recip_inh"], exc=True),
                             f"PF_{hc}_{side}",
                             pf_in_e if hc.endswith("E") else pf_in_f)
        for hc in JPF_HCS:
            n.add_connection(_syn(G["pf_recip_inh"], exc=False),
                             pf_in_e if hc.endswith("F")
                             else pf_in_f, f"PF_{hc}_{side}")
        # mechanosensors + afferent loops ride the same E/F families
        if (self.stance_fb or self.aff_loops) and G["ib_e_central"] > 0.0:
            for src in (f"HEEL_{side}", f"TOE_{side}"):
                for hc in ("HIP-E", "KNEE-E", "ANK-E"):
                    n.add_connection(_syn(G["ib_e_central"], exc=True),
                                     src, f"PF_{hc}_{side}")
                n.add_connection(_syn(G["ib_e_central"], exc=True),
                                 src, f"InE_{side}")
        if self.aff_loops:
            aff_e = f"AFF_E_{side}"
            aff_f = f"AFF_F_{side}"
            for hc in JPF_HCS:
                if hc.endswith("E"):
                    n.add_connection(_syn(G["aff_e_pf"], exc=True),
                                     aff_e, f"PF_{hc}_{side}")
                else:
                    n.add_connection(_syn(G["aff_f_pf"], exc=True),
                                     aff_f, f"PF_{hc}_{side}")
        # ---- 2026-09-24 evening PER-PF-LAYER CONTACT VARIANT (Ben:
        # "build it", coexists with the per-joint layering; source = his
        # block-editor drawing, 97n/109e). Heel = stance-phase reset of
        # the IPSILATERAL leg applied AT the PF layer; toe = dorsiflexion
        # inhibition ONLY. Defaults 0 = edges/neuron absent (bit-identical;
        # the TOEDF neuron is built ONLY when toe_df_inh > 0 - the v5
        # lesson: zero-g synapses alone change BLAS summation order).
        # heel IN -> PF_IN_E exc (drawing: g 0.5 into each micro-layer's
        # E-lamination IN; the runner keeps ONE shared PF_IN_E per side -
        # documented lumping, widen to per-joint INs if Ben wants).
        if self.stance_fb and G["heel_pf_layer"] > 0.0 and \
                f"HEEL_{side}" in self.idx:
            n.add_connection(_syn(G["heel_pf_layer"], exc=True),
                             f"HEEL_{side}", pf_in_e)
        # toe IN -> TOEDF IN -> ANK-F inhibition (drawing: toe g 5 ->
        # IN-PF_dorsiflexion_inhibit -> inhib HC-PF-Dorsiflexion; the
        # runner's dorsiflexor HC is PF_ANK-F). Both path legs scale
        # with the single toe_df_inh knob.
        if self.stance_fb and G["toe_df_inh"] > 0.0 and \
                f"TOE_{side}" in self.idx:
            toedf = f"TOEDF_{side}"
            if toedf not in self.idx:
                self._add(toedf, TAU["preset"], n)
            n.add_connection(_syn(G["toe_df_inh"], exc=True),
                             f"TOE_{side}", toedf)
            n.add_connection(_syn(G["toe_df_inh"], exc=False),
                             toedf, f"PF_ANK-F_{side}")

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
        if self.joint_pf:
            # 2026-09-18: joint-layer routing; weight = fitted (T1) or
            # anatomical-mapping 1.0 fallback
            base = act.rsplit("_", 1)[0]
            jw = _jpf_weights()
            grp_hcs = JPF_GROUP2HC.get(mi.groups[0], ())
            cross = JPF_CROSS.get(base, ())
            for hc in JPF_HCS:
                w = jw.get(hc, {}).get(base)
                if w is None:
                    w = 1.0 if (hc in grp_hcs or hc in cross) else 0.0
                if w > 0.0:
                    n.add_connection(
                        _syn(G["pf_to_mn"] * w, exc=True),
                        f"PF_{hc}_{mi.side}", mn)
        else:
            for phase in PF_PHASES:
                w = _group_weight(mi, W_PF_MN[phase])
                if w > 0.0:
                    for side in self.sides:
                        if mi.side == side:
                            n.add_connection(
                                _syn(G["pf_to_mn"] * w, exc=True),
                                f"PF_{phase}_{side}", mn)
        w_post = _group_weight(mi, W_POSTURE, POSTURE_OVERRIDE)
        if w_post > 0.0:
            n.add_connection(_syn(G["posture_to_mn"] * w_post, exc=True),
                             "POSTURE", mn)

        # ---- proprioceptive pathways ----
        # 2026-09-21 FULL LITERATURE RULES (G["full_rules"] > 0; Ben:
        # "how does the working Deng model connect things"). Deng A6 /
        # Di Russo rules 1-3: ONLY Ia homonymous is monosynaptic;
        # II and Ib reach their MNs through dedicated interneurons, and
        # the inhibitory INs mutually inhibit their antagonists.
        fr = G["full_rules"] > 0.0
        n.add_connection(_syn(G["ia_to_mn"], exc=True), ia, mn)
        if fr:
            # II -> (IIX exc IN) -> same MN   [rule 2a, disynaptic exc]
            iix = f"IIX_{act}"
            if iix not in self.idx:
                self._add(iix, TAU["afferent"], n)
                n.add_connection(_syn(1.0, exc=True), ii, iix)
            n.add_connection(_syn(G["ii_to_mn"], exc=True), iix, mn)
            # Ib -> (IBIN inh IN) -> same MN [rule 3, disynaptic autogenic]
            ibin = f"IBIN_{act}"
            if ibin not in self.idx:
                self._add(ibin, TAU["afferent"], n)
                n.add_connection(_syn(1.0, exc=True), ib, ibin)
            n.add_connection(_syn(G["ib_to_mn_inh"], exc=False),
                             ibin, mn)
            # Ib-IN <-> antagonist Ib-IN mutual inhibition [rule 3]
            # (IIIN -> antagonist MN wired in the antagonist loop below)
            for act2, mi2 in self.muscles.items():
                if mi2.side == mi.side and act2 != act and \
                        mi2.groups[0] in ANTAGONIST.get(mi.groups[0], ()) \
                        and f"IBIN_{act2}" in self.idx:
                    n.add_connection(_syn(0.5, exc=False), ibin,
                                     f"IBIN_{act2}")
        else:
            n.add_connection(_syn(G["ii_to_mn"], exc=True), ii, mn)
            n.add_connection(_syn(G["ib_to_mn_inh"], exc=False),
                             ib, mn)

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
        # joint-layer mode: route central projections to this muscle's
        # own joint half-center (PF_HIP-E_r etc.) instead of the absent
        # phase cells (the 09-18 build wired heel/toe but missed this
        # section - found 2026-09-21, "Population not found PF_E1_r")
        if self.joint_pf:
            tgt_e = ([f"PF_{h}_{mi.side}" for h in
                      JPF_GROUP2HC.get(grp, ()) if h.endswith("E")]
                     or [f"PF_{h}_{mi.side}" for h in
                         ("HIP-E", "KNEE-E", "ANK-E")])
            tgt_f = ([f"PF_{h}_{mi.side}" for h in
                      JPF_GROUP2HC.get(grp, ()) if h.endswith("F")]
                     or [f"PF_{h}_{mi.side}" for h in
                         ("HIP-F", "KNEE-F", "ANK-F")])
        else:
            tgt_e = [f"PF_E1_{mi.side}", f"PF_E2_{mi.side}"]
            tgt_f = [f"PF_F1_{mi.side}", f"PF_F2_{mi.side}"]
        if G["ib_e_central"] > 0.0 and grp in EXTENSOR_STANCE_GROUPS:
            for tgt in (*tgt_e, f"RG_E_{mi.side}", f"InE_{mi.side}"):
                n.add_connection(_syn(G["ib_e_central"], exc=True),
                                 ib, tgt)
        if grp in ("hip_flex", "knee_flex", "ankle_df", "hip_add",
                   "trunk_flex"):
            for tgt in (*tgt_f, f"RG_F_{mi.side}", f"InF_{mi.side}"):
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
            for tgt in (*tgt_e, f"RG_E_{mi.side}", f"InE_{mi.side}"):
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
                # phase gate: PF_F1 excites the IaIN (Deng A6 PF->IaIN
                # 0.5); joint-layer mode gates by the muscle's own F HC
                gate_f1 = (tgt_f[0] if self.joint_pf
                           else f"PF_F1_{mi.side}")
                n.add_connection(_syn(0.5, exc=True), gate_f1, iain)
                # recurrent disinhibition (RC -> IaIN inh)
                if self.renshaw and f"RC_{act}" in self.idx:
                    n.add_connection(_syn(G["renshaw"], exc=False),
                                     f"RC_{act}", iain)
        # II -> (IIIN inh IN) -> antagonist MN [rule 2b]: created once,
        # edges to the antagonist MNs in the loop below
        if fr:
            iiin = f"IIIN_{act}"
            if iiin not in self.idx:
                self._add(iiin, TAU["afferent"], n)
                n.add_connection(_syn(1.0, exc=True), ii, iiin)
        for ant in ANTAGONIST.get(mi.groups[0], ()):
            for act2, mi2 in self.muscles.items():
                if mi2.side == mi.side and mi2.groups[0] == ant:
                    if self.ia_in:
                        n.add_connection(
                            _syn(G["ia_to_antagonist"], exc=False),
                            f"IaIN_{act}", f"MN_{act2}")
                        # rule 1 addendum: IaIN <-> antagonist IaIN
                        # mutual inhibition (Deng A6)
                        if fr and f"IaIN_{act2}" in self.idx:
                            n.add_connection(_syn(0.5, exc=False),
                                             f"IaIN_{act}",
                                             f"IaIN_{act2}")
                    else:
                        n.add_connection(
                            _syn(G["ia_to_antagonist"], exc=False),
                            ia, f"MN_{act2}")
                    if fr and f"IIIN_{act}" in self.idx:
                        n.add_connection(
                            _syn(G["ia_to_antagonist"], exc=False),
                            f"IIIN_{act}", f"MN_{act2}")

        # 2026-09-24 per-PF-layer flexion-afferent variant (Ben's
        # drawing: IN-IaIN -> HC-PF-F g 0.5 and IN-IIe -> HC-PF-F
        # g 0.5): FLEXOR-group afferent INs reinforce their OWN joint's
        # F half-center (flexion afference supports the flexion layer).
        # Requires joint_pf (the target HC exists only there) plus the
        # IN population that each gain rides on; default 0 = absent.
        if self.joint_pf and mi.groups[0] in ("hip_flex", "knee_flex",
                                              "ankle_df"):
            f_hc = (f"PF_{JPF_GROUP2HC[mi.groups[0]][0]}_{mi.side}")
            if G["ia_pf_f"] > 0.0 and self.ia_in \
                    and f"IaIN_{act}" in self.idx:
                n.add_connection(_syn(G["ia_pf_f"], exc=True),
                                 f"IaIN_{act}", f_hc)
            if G["ii_pf_f"] > 0.0 and fr and f"IIX_{act}" in self.idx:
                n.add_connection(_syn(G["ii_pf_f"], exc=True),
                                 f"IIX_{act}", f_hc)

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
            # swing-gate source: the F1 phase cell, or in joint-layer
            # mode the knee flexor half-center (swing gate for the knee)
            gate = (f"PF_KNEE-F_{mi.side}" if self.joint_pf
                    else f"PF_F1_{mi.side}")
            if kname not in self.idx:
                self._add(kname, TAU["ib_exc"], n)
                n.add_connection(_syn(1.5, exc=True), gate, kname)
                # 2026-09-21 crossed KINH drive (gain G["contra_kinh"],
                # default 0 = edge absent): the CONTRALATERAL heel-load
                # IN also excites this side's KINH - opposite heel
                # strike suppresses THIS side's extensor MNs = forced
                # stance->swing transition. s3c/s3d verdict: RG-level
                # nudges (contra_swing) do not release a loaded jammed
                # leg; the suppression must reach the MN pools.
                if G["contra_kinh"] > 0.0:
                    other = "l" if mi.side == "r" else "r"
                    if f"HEEL_{other}" in self.idx:
                        n.add_connection(
                            _syn(G["contra_kinh"], exc=True),
                            f"HEEL_{other}", kname)
            n.add_connection(_syn(G["f1_kneext_inh"], exc=False), kname, mn)
        # v6b: same swing-gated suppression onto ankle PF pools (reuses
        # the KINH IN; separate gain)
        if self.f1_kneext_inh and mi.groups[0] == "ankle_pf":
            kname = f"KINH_{mi.side}"
            gate = (f"PF_ANK-F_{mi.side}" if self.joint_pf
                    else f"PF_F1_{mi.side}")
            if kname not in self.idx:
                self._add(kname, TAU["ib_exc"], n)
                n.add_connection(_syn(1.5, exc=True), gate, kname)
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

        # ---- goal2 vestibular-analog tone (2026-09-23; Ben's
        # standing-balance-stage request). Traces to (a) Ben's request
        # text ("vestibular analog = pelvis-tilt / COM sensors driving
        # extensor tone + ankle strategy") and (b) the SCONE Tutorial-3a
        # balance controller (vestibular BodyPointReflex from torso to
        # all major muscles; verified local copy in Documents\SCONE\
        # Tutorials\controllers\ControllerReflexBalance.scone lines
        # 25-43). Direct VEST->MN edges follow the established BAL_*
        # balance-cell pattern above (brainstem surrogates wire direct,
        # not through spinal INs). Each edge gated by ITS OWN gain > 0
        # (a 0-gain synapse still changes summation order).
        if self.vest:
            for act, mi in self.muscles.items():
                vsrc = f"VEST_{mi.side}"
                g = mi.groups[0]
                if G["vest_ext"] > 0.0 and g in ("knee_ext", "ankle_pf",
                                                 "hip_ext", "trunk_ext"):
                    # antigravity extensor tone (vestibulospinal)
                    n.add_connection(_syn(G["vest_ext"], exc=True),
                                     vsrc, f"MN_{act}")
                elif G["vest_flex_inh"] > 0.0 and g in (
                        "hip_flex", "knee_flex", "ankle_df", "trunk_flex"):
                    # reciprocal flexor inhibition (LVST) - awaiting
                    # Ben's connectome-spec confirmation
                    n.add_connection(_syn(G["vest_flex_inh"], exc=False),
                                     vsrc, f"MN_{act}")

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
