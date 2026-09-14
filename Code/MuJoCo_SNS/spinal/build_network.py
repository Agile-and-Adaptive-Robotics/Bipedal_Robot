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
from params import (AFF, DT, E_HI, G, MOD, PF_SHAPE, PHASE_RESET,
                    POSTURE_OVERRIDE, TAU, W_PF_MN, W_POSTURE)

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.networks import Network

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
        # at gains 0 the compiled network must be BYTE-IDENTICAL to v4b
        # (same neurons, same matrix size -> bit-identical trajectory;
        # zero-conductance synapses alone still change the BLAS summation
        # order and chaos amplifies the 1e-16 drift - reggate_v5_0 lesson)
        self.phase_reset = bool(G["phase_reset_e"] > 0.0
                                or G["phase_reset_f"] > 0.0)
        # same pattern for the phase-3 swing-knee quad suppression
        self.f1_kneext_inh = bool(G["f1_kneext_inh"] > 0.0
                                  or G["f1_anklepf_inh"] > 0.0)
        # Renshaw recurrent inhibition (Deng Table A6)
        self.renshaw = bool(G["renshaw"] > 0.0)

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

        # ---- cross-side coordination (F cells strongly, E cells weakly).
        # interleg=False removes it entirely: independent left/right
        # rhythm generators (deafferented air-stepping preparation).
        if self.interleg:
            for a, b in (("r", "l"), ("l", "r")):
                n.add_connection(_syn(G["rg_mutual_inh"], exc=False),
                                 f"RG_F_{a}", f"RG_F_{b}")
                n.add_connection(_syn(0.5 * G["rg_mutual_inh"], exc=False),
                                 f"RG_E_{a}", f"RG_E_{b}")

    # ------------------------------------------------------------------ parts
    def _add(self, name: str, tau: float, n: Network):
        n.add_neuron(_neu(tau), name=name)
        self.idx[name] = len(self.idx)

    def _build_rg(self, n: Network, side: str):
        rg_e, rg_f = f"RG_E_{side}", f"RG_F_{side}"
        ad_e, ad_f = f"ADAP_E_{side}", f"ADAP_F_{side}"
        for name, tau in ((rg_e, TAU["rg"]), (rg_f, TAU["rg"]),
                          (ad_e, TAU["rg_adapt"]), (ad_f, TAU["rg_adapt"])):
            self._add(name, tau, n)

        # half-center mutual inhibition
        n.add_connection(_syn(G["rg_mutual_inh"], exc=False), rg_e, rg_f)
        n.add_connection(_syn(G["rg_mutual_inh"], exc=False), rg_f, rg_e)
        # slow self-adaptation: burst termination (sets cycle period)
        n.add_connection(_syn(G["rg_adapt_inh"], exc=True), rg_e, ad_e)
        n.add_connection(_syn(G["rg_adapt_inh"], exc=True), rg_f, ad_f)
        n.add_connection(_syn(G["rg_adapt_inh"], exc=False), ad_e, rg_e)
        n.add_connection(_syn(G["rg_adapt_inh"], exc=False), ad_f, rg_f)
        # descending drive raises frequency; stance-biased split sets duty
        n.add_connection(_syn(G["descend_to_rg_e"], exc=True), "DRIVE", rg_e)
        n.add_connection(_syn(G["descend_to_rg_f"], exc=True), "DRIVE", rg_f)
        # posture tonic bias keeps the E cell (load-bearing side) ready
        n.add_connection(_syn(G["posture_to_rg_e"], exc=True), "POSTURE", rg_e)

        # ---- v5 sensory phase-reset (Rybak-style hip afferents switching
        # the half-centers): per-side input ports HIP_EXT_SIG / HIP_FLEX_SIG
        # drive PRESET_E / PRESET_F interneurons; PRESET_E excites RG-E and
        # inhibits RG-F (prolongs stance), PRESET_F excites RG-F and
        # inhibits RG-E (triggers swing). Topology is CONDITIONAL on the
        # gains (self.phase_reset): at 0/0 nothing is added and the
        # compiled network is byte-identical to v4b (regression-gated by
        # runner --fitted --best --eval).
        if self.phase_reset:
            pre_e, pre_f = f"PRESET_E_{side}", f"PRESET_F_{side}"
            self._add(pre_e, TAU["preset"], n)
            self._add(pre_f, TAU["preset"], n)
            n.add_input(pre_e)
            self.inputs.append("HIP_EXT_SIG_" + side)
            n.add_input(pre_f)
            self.inputs.append("HIP_FLEX_SIG_" + side)
            g_e, g_f = G["phase_reset_e"], G["phase_reset_f"]
            inh = PHASE_RESET.get("inh", 1.0)
            n.add_connection(_syn(g_e, exc=True), pre_e, rg_e)
            n.add_connection(_syn(g_e * inh, exc=False), pre_e, rg_f)
            n.add_connection(_syn(g_f, exc=True), pre_f, rg_f)
            n.add_connection(_syn(g_f * inh, exc=False), pre_f, rg_e)

    def _build_pf(self, n: Network, side: str):
        rg_of = {"E1": f"RG_E_{side}", "E2": f"RG_E_{side}",
                 "F1": f"RG_F_{side}", "F2": f"RG_F_{side}"}
        for phase in PF_PHASES:
            tau_m, tau_a = PF_SHAPE[phase]
            pf = f"PF_{phase}_{side}"
            pfa = f"PFA_{phase}_{side}"
            self._add(pf, TAU["pf"] * tau_m, n)
            self._add(pfa, TAU["pf_adapt"] * tau_a, n)

            n.add_connection(_syn(G["rg_to_pf"], exc=True), rg_of[phase], pf)
            n.add_connection(_syn(1.5, exc=True), pf, pfa)
            n.add_connection(_syn(1.5, exc=False), pfa, pf)
            n.add_connection(_syn(G["drive_to_pf"], exc=True), "DRIVE", pf)

        # conflicting phase windows (per side)
        for a, b in (("E2", "F1"), ("F2", "E1"), ("E1", "F1"), ("E2", "F2")):
            n.add_connection(_syn(G["pf_recip_inh"], exc=False),
                             f"PF_{a}_{side}", f"PF_{b}_{side}")
            n.add_connection(_syn(G["pf_recip_inh"], exc=False),
                             f"PF_{b}_{side}", f"PF_{a}_{side}")

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

        # ---- Renshaw recurrent inhibition (Deng Table A6: MN->RC 0.5
        # exc; RC->MN 0.5 inh; RC<->RC 0.5 inh), gain-scaled by
        # G["renshaw"]; topology present only when that gain > 0.
        if self.renshaw:
            rc = f"RC_{act}"
            g_r = G["renshaw"]
            n.add_connection(_syn(1.0, exc=True), mn, rc)
            n.add_connection(_syn(g_r, exc=False), rc, mn)
            # RC<->RC: wire each unordered same-side pair once (adding
            # both directions here would double the conductance)
            for act2, mi2 in self.muscles.items():
                if mi2.side == mi.side and act2 != act \
                        and f"RC_{act2}" in self.idx \
                        and self._order[act2] > self._order[act]:
                    n.add_connection(_syn(g_r, exc=False), rc,
                                     f"RC_{act2}")

        # Ia reciprocal inhibition of antagonist MN pools
        for ant in ANTAGONIST.get(mi.groups[0], ()):
            for act2, mi2 in self.muscles.items():
                if mi2.side == mi.side and mi2.groups[0] == ant:
                    n.add_connection(_syn(G["ia_to_antagonist"], exc=False),
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
