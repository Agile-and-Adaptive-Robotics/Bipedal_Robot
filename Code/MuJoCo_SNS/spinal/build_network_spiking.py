"""SPIKING MIRROR of build_network.py (Ben's SPIKING_MIRROR_PLAN.md,
2026-09-18; built 2026-10-02 goal-1 spiking campaign).

Hybrid architecture (the plan's ruling): the MN layer stays NON-spiking
and is driven through spike->graded synapses (the SNS-Toolbox hybrid
pattern); afferent encoders, RG/PF cells, INs, and commissurals become
SPIKING (SpikingNeuron LIF).  Topology is mirrored verbatim from
build_network.py - same population names, same conditional-gain checks,
same input ports - so _fix_check-style topology gates can be run against
both builds (tools/topology_mirror_check.py in
reports_spiking_20261002).

Cell-class mapping (SPIKING_MIRROR_PLAN.md table):
  DRIVE/POSTURE/BAL_*/VEST (brainstem surrogates)  -> NonSpikingNeuron
      (analog commands; graded synapses onto spiking cells carry them)
  RG-E/F  -> adapting LIF half-centers (threshold-increment SFA; the
      toolbox NaP class is non-spiking-only, so the plan's documented
      fallback "adapting-LIF pairs with mutual inhibition" is used;
      tau_theta reuses TAU["rg_nap_h"] = the analog period knob)
  PF cells (phase or joint-layer) -> adapting LIF (mild SFA)
  INs (InE/InF, PF_IN, CIN, IaIN, IIX, IBIN, IIIN, RC, KINH, LBIN,
      IBEXC, TOEDF, AFF_E/F) -> plain LIF relays
  Ia/II/Ib afferent encoders + HEEL/TOE mechanosensor INs -> LIF with
      real mV ranges; the runner's nA port currents map to firing rate
      (threshold placed so the documented AFF.i0_ii = 1 nA baseline
      tone sits just subthreshold - see SPIKE["aff_v_thr"])
  MN pools -> NonSpikingNeuron, UNCHANGED 0..5 mV frame, activation
      a = clip(V/E_HI, 0, 1) unchanged; spiking sources reach them via
      spike synapses whose g_increment is CALIBRATED (calibrate_spiking
      .py -> spiking_calibration.json) so a steady presynaptic train
      reproduces the non-spiking steady-state current (plan rule 2).

Real mV ranges (plan rule 1): spiking cells rest -70 mV, threshold
-50 mV, reset -60 mV, membrane_conductance 0.1 uS (Rin = 10 MOhm so
the existing nA-scale drive currents reach threshold).  Afferent
encoders use threshold -58 / reset -64 (the I->rate map chosen from
the existing afferent gain ranges, plan table row 1).

Sub-stepping (plan rule 3): the compiled network steps at NET_DT =
plant_dt / N_SUB (default 0.5 ms inside the 2 ms plant step);
SpikingNetworkSpiking.step() runs the sub-steps at a constant input
vector, exactly like the non-spiking single step.

Runner compatibility: the runner consumes net.compiled.V[net.idx[...]]
for the RG stance gates and neuro logging.  Spiking membranes are not
analog levels, so after the circuit is wired this class appends
readout tap neurons RD_<cell>_<side> (non-spiking sinks calibrated so a
50 Hz train reads full-scale E_HI) for RG_E/RG_F and every PF cell,
and REPOINTS net.idx[<cell>] at the readout while net.ridx[<cell>]
keeps the raw spiking cell (needed by spike-counting gates).  Readouts
are pure sinks - they cannot feed back into the circuit.

The non-spiking pipeline is untouched; this module is selected ONLY via
AARL_NET=spiking in build_network.build().
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

from sns_toolbox.connections import NonSpikingSynapse, SpikingSynapse
from sns_toolbox.neurons import NonSpikingNeuron, SpikingNeuron
from sns_toolbox.networks import Network

# Synapse transfer: I = g * clip((V_pre - e_lo)/(e_hi - e_lo), 0, 1) * (E_rev - V_post)
SYN_E_LO, SYN_E_HI = 0.0, E_HI   # presynaptic saturation window of ANALOG sources
E_REV_EXC = 8.0                  # mV, excitatory reversal in the MN 0..5 mV frame
E_REV_INH = -E_HI                # mV, inhibitory reversal in the MN frame
# Reversals in the spiking (real-mV) frame, plan rule 2:
S_REV_EXC = 0.0                  # mV above rest
S_REV_INH = -70.0                # mV = rest (shunting)

PF_PHASES = ("E1", "E2", "F1", "F2")

# ---- spiking cell parameters (real mV; plan rule 1) ----------------------
# MEASURED TRAPS (2026-10-02, probe_lif.py / probe_hybrid.py):
#  - SNS_Numpy.forward computes dV = (dt/(Cm/g_m)) * (-g_m*(V-Vr) + I),
#    which is only dimensionally correct at g_m = 1 uS (any other
#    conductance silently rescales every input current).  All spiking
#    cells therefore use membrane_conductance = 1.0 and get their input
#    resistance from a TONIC BIAS that parks the equilibrium just below
#    threshold (the physiological background-synaptic-current stand-in).
#  - threshold_time_constant is consumed in SECONDS by the same forward
#    (dt / tau_theta with dt in s) - pass seconds, not the docstring's ms.
SPIKE = dict(
    v_rest=-70.0,       # mV
    v_thr=-50.0,        # mV, main threshold
    v_reset=-60.0,      # mV
    bias=16.0,          # nA tonic background current: equilibrium
                        # V_eq = v_rest + bias = -54 mV, i.e. 4 mV below
                        # threshold - the analog cells' near-threshold
                        # operating point
    aff_v_thr=-52.0,    # afferent encoder threshold: 2 mV above V_eq so
    aff_v_reset=-58.0,  # the runner's 0..3 nA port currents map
                        # monotonically to ~0..25 Hz (the plan's chosen
                        # I->rate map); the AFF.i0_ii 1 nA baseline tone
                        # stays just subthreshold
    rc_v_thr=-52.5,     # Renshaw threshold (MN drive arrives attenuated
    rc_v_reset=-58.0,   # through the graded rate-encoder synapse)
    rg_thr_inc=1.0,     # mV per spike, RG burst-termination SFA (the
                        # analog tau_h role: steady lift ~ inc*f*tau_theta
                        # ~ 1.0*25*0.35 ~ 9 mV = burst termination)
    rg_thr_floor=-48.0, # mV, RG adaptation ceiling-ish floor
    pf_thr_inc=0.3,     # mV per spike, mild PF adaptation
    pf_tau_theta=0.15,  # s
    pf_thr_floor=-49.0, # mV
)
# synaptic time constants (s) of spike synapses; inhibition decays slower so
# lamination INs hold the opposite half-center down across interspike
# intervals (tau_inh was chosen so a 40 Hz train holds >80% mean conductance)
TAU_SYN = dict(exc=0.005, inh=0.020, readout=0.020)
F_FULL = 50.0         # Hz mapped to full-scale E_HI on the RD_ readout taps

# ---- calibration factors (calibrate_spiking.py -> spiking_calibration.json;
# defaults below are the analytic mean-conductance mapping at the nominal
# sat_ref/f_ref - the calibration script refines them by rig bisection) ----
_CAL_PATH = Path(__file__).parent / "spiking_calibration.json"
CAL_DEFAULTS = dict(
    f_ref=40.0,          # Hz, design within-burst rate of spiking RG/PF
    sat_ref=0.9567,      # active-phase mean saturation of analog sources
    k_s2s_exc=0.10,      # g_inc = k * g_ns   (spiking -> spiking, exc;
    k_s2s_inh=0.10,      # RATE-calibrated - see calibrate_spiking.py)
    k_sn_exc=5.41,       # (spiking -> MN, exc; mean-conductance match)
    k_sn_inh=1.23,       # (spiking -> MN, inh)
    k_ns=0.35,           # graded (analog -> spiking) DRIVE rate-map gain
    g_inc_readout=1.709, # readout tap increment (50 Hz -> V = E_HI)
    rg_thr_inc=4.0,      # RG SFA increment (tune_rg_pair.py sweep winner)
)
try:  # calibration file wins when present (written by calibrate_spiking.py)
    _cal_disk = json.loads(_CAL_PATH.read_text(encoding="utf-8"))
    CAL = {**CAL_DEFAULTS, **{k: v for k, v in _cal_disk.items()
                              if k in CAL_DEFAULTS}}
except (FileNotFoundError, ValueError):
    CAL = dict(CAL_DEFAULTS)

# ---- 2026-09-18 joint-layer PF mirror (same tables as build_network.py) --
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
    """Muscle-level PF weights per joint-layer HC (same loader as the
    non-spiking build)."""
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

ANTAGONIST = {
    "hip_ext": ("hip_flex",), "hip_flex": ("hip_ext",),
    "hip_abd": ("hip_add",), "hip_add": ("hip_abd",),
    "knee_ext": ("knee_flex",), "knee_flex": ("knee_ext",),
    "ankle_pf": ("ankle_df",), "ankle_df": ("ankle_pf",),
    "trunk_ext": ("trunk_flex",), "trunk_flex": ("trunk_ext",),
}


# --------------------------------------------------------------- cell makers
def _neu(tau: float) -> NonSpikingNeuron:
    """Non-spiking RC neuron, ANALOG 0..5 mV frame (identical to the
    non-spiking build): DRIVE/POSTURE/BAL/VEST cells, MN pools, readouts."""
    return NonSpikingNeuron(membrane_capacitance=float(tau),
                            membrane_conductance=1.0,
                            resting_potential=0.0, bias=0.0)


def _spk(tau: float, kind: str = "in") -> SpikingNeuron:
    """Spiking LIF in the real-mV frame. kind: 'rg' (adapting half-center),
    'pf' (mildly adapting PF cell), 'in' (plain relay), 'aff' (afferent
    encoder / mechanosensor thresholds), 'rc' (Renshaw thresholds).
    Membrane tau = the analog build's value for that class (plan rule 3);
    the toolbox's g_m=1 convention holds (see SPIKE comment)."""
    s = SPIKE
    if kind == "aff":
        thr, reset, floor = s["aff_v_thr"], s["aff_v_reset"], s["aff_v_thr"]
        inc, tau_theta = 0.0, 0.10
    elif kind == "rc":
        thr, reset, floor = s["rc_v_thr"], s["rc_v_reset"], s["rc_v_thr"]
        inc, tau_theta = 0.0, 0.10
    elif kind == "rg":
        thr, reset, floor = s["v_thr"], s["v_reset"], s["rg_thr_floor"]
        inc = float(CAL.get("rg_thr_inc", s["rg_thr_inc"]))
        tau_theta = float(TAU["rg_nap_h"])
    elif kind == "pf":
        thr, reset, floor = s["v_thr"], s["v_reset"], s["pf_thr_floor"]
        inc, tau_theta = s["pf_thr_inc"], s["pf_tau_theta"]
    else:
        thr, reset, floor = s["v_thr"], s["v_reset"], s["v_thr"]
        inc, tau_theta = 0.0, 0.10
    return SpikingNeuron(
        threshold_time_constant=tau_theta,           # SECONDS (measured)
        threshold_initial_value=thr,
        threshold_proportionality_constant=0.0,
        threshold_leak_rate=1.0,
        threshold_increment=inc,
        threshold_floor=floor,
        reset_potential=reset,
        membrane_capacitance=float(tau),
        membrane_conductance=1.0,
        resting_potential=s["v_rest"], bias=s["bias"])


# ------------------------------------------------------------ synapse makers
def _syn(g: float, exc: bool) -> NonSpikingSynapse:
    """Graded analog->analog synapse (unchanged from the non-spiking
    build): POSTURE/BAL/VEST -> MN."""
    return NonSpikingSynapse(
        max_conductance=float(g),
        reversal_potential=E_REV_EXC if exc else E_REV_INH,
        e_lo=SYN_E_LO, e_hi=SYN_E_HI)


def _ns2s(g: float, exc: bool) -> NonSpikingSynapse:
    """Graded ANALOG-source -> SPIKING-cell synapse (e.g. DRIVE -> RG,
    MN -> RC where the graded synapse acts as the rate encoder, plan
    table row 5).  Same analog saturation window on the source; reversal
    in the real-mV frame; conductance scaled by CAL['k_ns'] so the
    spiking cell's depolarization matches the analog operating-point
    mapping (calibrated, not hand-tuned)."""
    return NonSpikingSynapse(
        max_conductance=float(g) * CAL["k_ns"],
        reversal_potential=S_REV_EXC if exc else S_REV_INH,
        e_lo=SYN_E_LO, e_hi=SYN_E_HI)


def _ssyn(g: float, exc: bool, onto_analog: bool = False) -> SpikingSynapse:
    """Spike synapse from a spiking source.  onto_analog=True targets the
    MN 0..5 mV frame (hybrid pattern; reversal E_REV_*), else the
    real-mV spiking frame (reversal 0 / -70).  g_increment = k * g
    (mean-conductance matching: a steady f_ref train reproduces the
    analog synapse's saturating conductance - see calibrate_spiking.py);
    g_max = 8 x g_inc so a few spikes can sum without clipping."""
    k = (CAL["k_sn_exc"] if exc else CAL["k_sn_inh"]) if onto_analog \
        else (CAL["k_s2s_exc"] if exc else CAL["k_s2s_inh"])
    g_inc = float(g) * k
    g_max = max(8.0 * g_inc, 1e-6)
    return SpikingSynapse(
        max_conductance=g_max,
        reversal_potential=(E_REV_EXC if exc else E_REV_INH) if onto_analog
        else (S_REV_EXC if exc else S_REV_INH),
        time_constant=TAU_SYN["exc"] if exc else TAU_SYN["inh"],
        transmission_delay=0,
        conductance_increment=g_inc)


def _group_weight(muscle: MuscleInfo, table: dict[str, float],
                  override: dict[str, float] | None = None) -> float:
    """Primary-group weight + half of any secondary-group weight
    (identical to the non-spiking build)."""
    if override and muscle.base in override:
        w = override[muscle.base]
    else:
        w = table.get(muscle.groups[0], 0.0)
    for sec in muscle.groups[1:]:
        w += 0.5 * table.get(sec, 0.0)
    return w


@dataclass
class SpinalNetworkSpiking:
    """Mirror of build_network.SpinalNetwork with the hybrid cell mapping.
    Attributes consumed by runner.py are name-compatible (muscles, sides,
    idx, inputs, mn_names, compiled, step, make_inputs, input_index, vest,
    stance_fb)."""
    muscles: dict[str, MuscleInfo]
    sides: tuple[str, ...]
    interleg: bool = True
    net: Network = field(init=False)
    compiled: object = field(init=False, default=None)
    idx: dict[str, int] = field(default_factory=dict)
    ridx: dict[str, int] = field(default_factory=dict)   # raw spiking cells
    inputs: list[str] = field(default_factory=list)
    mn_names: dict[str, str] = field(default_factory=dict)
    aff_names: dict[str, dict[str, str]] = field(default_factory=dict)
    ib_exc_groups: dict[str, tuple[str, ...]] = field(default_factory=dict)
    n_sub: int = 4              # sub-steps per plant step (2 ms / 0.5 ms)
    spike_counts: np.ndarray = field(init=False, default=None)

    # ------------------------------------------------------------------ build
    def __post_init__(self):
        self.net = Network(name="gait2392 spinal (spiking mirror)")
        n = self.net
        # conditional-topology flags: IDENTICAL gain checks to the
        # non-spiking build (same G knobs decide what exists)
        self.f1_kneext_inh = bool(G["f1_kneext_inh"] > 0.0
                                  or G["f1_anklepf_inh"] > 0.0)
        self.renshaw = bool(G["renshaw"] > 0.0)
        self.stance_fb = bool(G["heel_rge"] > 0.0 or G["toe_rge"] > 0.0
                              or G["ib_rge"] > 0.0
                              or G["ib_e_central"] > 0.0
                              or G["heel_pf_layer"] > 0.0
                              or G["toe_df_inh"] > 0.0
                              or G["heel_in_f_exc"] > 0.0)
        self.ia_in = bool(G["ia_in"] > 0.0)
        self.aff_loops = bool(G["aff_e_rg"] > 0.0 or G["aff_f_rg"] > 0.0
                              or G["aff_e_pf"] > 0.0 or G["aff_f_pf"] > 0.0)
        self.vest = bool(G["vest_ext"] > 0.0 or G["vest_flex_inh"] > 0.0)
        self.joint_pf = bool(G.get("joint_pf", 0.0) > 0.0)
        if self.joint_pf:
            _jpf_weights()

        # ---- descending / balance cells (shared, ANALOG - unchanged) ----
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

        # ---- goal2 vestibular-analog cells (analog brainstem; unchanged)
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

        # ---- muscles: create all neurons first, then wire ----
        self._order = {a: i for i, a in enumerate(self.muscles)}
        for act, mi in self.muscles.items():
            self._add_muscle_neurons(n, act, mi)
        for act, mi in self.muscles.items():
            self._wire_muscle(n, act, mi)
        self._complete_in_mutual(n)
        self._wire_balance(n)

        # ---- cross-side commissurals (mirror of the non-spiking block) --
        if self.interleg:
            for a, b in (("r", "l"), ("l", "r")):
                cf, ce = f"CIN_F_{a}", f"CIN_E_{a}"
                self._add(cf, TAU["rg"], n)
                self._add(ce, TAU["rg"], n)
                n.add_connection(_ssyn(G["c1_gain"] * G["rg_mutual_inh"],
                                       exc=True), f"RG_F_{a}", cf)
                n.add_connection(_ssyn(G["c1_gain"] * G["rg_mutual_inh"],
                                       exc=False), cf, f"RG_F_{b}")
                if G["v3_gain"] > 0.0:
                    n.add_connection(
                        _ssyn(G["v3_gain"] * G["rg_mutual_inh"], exc=True),
                        f"RG_E_{a}", ce)
                    n.add_connection(
                        _ssyn(G["v3_gain"] * G["rg_mutual_inh"], exc=True),
                        ce, f"InE_{b}")
                if G["v3_to_ibexc"] > 0.0:
                    for grp in self.ib_exc_groups.get(b, ()):
                        gname = f"IBEXC_{grp}_{b}"
                        if gname in self.idx:
                            n.add_connection(
                                _ssyn(G["v3_to_ibexc"], exc=True), ce, gname)

        # ---- runner-facing analog readout taps (added LAST; pure sinks) -
        # RD_<cell> non-spiking neurons low-pass the spike trains; the
        # runner's stance gates + neuro logging read net.idx[...], which is
        # repointed at the readout.  Raw cells stay reachable via ridx.
        taps: list[str] = []
        for side in self.sides:
            taps += [f"RG_E_{side}", f"RG_F_{side}"]
            if self.joint_pf:
                taps += [f"PF_{hc}_{side}" for hc in JPF_HCS]
            else:
                taps += [f"PF_{ph}_{side}" for ph in PF_PHASES]
        for cell in taps:
            rd = f"RD_{cell}"
            n.add_neuron(_neu(0.05), name=rd)
            self.idx[rd] = len(self.idx)
            g_inc = float(CAL["g_inc_readout"])
            n.add_connection(
                SpikingSynapse(max_conductance=4.0 * g_inc,
                               reversal_potential=E_REV_EXC,
                               time_constant=TAU_SYN["readout"],
                               transmission_delay=0,
                               conductance_increment=g_inc),
                cell, rd)
            self.ridx[cell] = self.idx[cell]   # raw spiking cell
            self.idx[cell] = self.idx[rd]      # runner-facing analog level

    # ------------------------------------------------------------------ parts
    def _complete_in_mutual(self, n: Network, force: bool = False):
        """2026-10-03 audit-fix MIRROR (WIRING_RULINGS_20261003.md, s3k
        finding F3, propagated to the spiking twin the same day so the
        gate-1 edge-mirror contract holds against the FIXED non-spiking
        builder): complete the IaIN<->IaIN and IBIN<->IBIN MUTUAL
        inhibition (Deng Table A6 / Rybak 2006 Table 2: both directions).
        This builder has the same lazy-creation trap the audit measured in
        build_network.py (only the later->earlier direction of each pair
        existed; 218 of 436 directed edges per family under full_rules).
        Adds ONLY the missing directed edges as spike synapses at the same
        0.5 conductance-increment scale as the existing edges. Default
        build never enters (full_rules 0, ia_in 0: no IaIN/IBIN cells at
        all), so the stage-1 mirror counts (422 pops / 1186 edges) are
        unaffected."""
        if not force and not (G["full_rules"] > 0.0):
            return
        have = {(n.populations[c["source"]]["name"],
                 n.populations[c["destination"]]["name"])
                for c in n.connections}
        for act, mi in self.muscles.items():
            for ant in ANTAGONIST.get(mi.groups[0], ()):
                for act2, mi2 in self.muscles.items():
                    if mi2.side != mi.side or mi2.groups[0] != ant:
                        continue
                    for fam in ("IaIN", "IBIN"):
                        if f"{fam}_{act2}" not in self.idx:
                            continue
                        e = (f"{fam}_{act}", f"{fam}_{act2}")
                        if e not in have:
                            n.add_connection(_ssyn(0.5, exc=False), *e)

    def _add(self, name: str, tau: float, n: Network, kind: str = "in"):
        n.add_neuron(_spk(tau, kind), name=name)
        self.idx[name] = len(self.idx)

    def _build_rg(self, n: Network, side: str):
        rg_e, rg_f = f"RG_E_{side}", f"RG_F_{side}"
        ine, inf = f"InE_{side}", f"InF_{side}"
        # RG half-centers: adapting-LIF conditional bursters (the plan's
        # documented fallback - 1.5.2 ships no spiking bursting class; the
        # NaP class is non-spiking-only, verified 2026-10-02).  Burst
        # termination = threshold adaptation; tau_theta reuses the analog
        # rg_nap_h period knob, so curriculum rg_nap_h values transfer.
        self._add(rg_e, TAU["rg"], n, kind="rg")
        self._add(rg_f, TAU["rg"], n, kind="rg")
        for name in (ine, inf):
            self._add(name, TAU["rg"], n)

        # IN-laminated mutual inhibition (same topology; spike synapses)
        n.add_connection(_ssyn(G["rg_mutual_inh"], exc=True), rg_e, ine)
        n.add_connection(_ssyn(G["rg_mutual_inh"], exc=False), ine, rg_f)
        n.add_connection(_ssyn(G["rg_mutual_inh"], exc=True), rg_f, inf)
        n.add_connection(_ssyn(G["rg_mutual_inh"], exc=False), inf, rg_e)

        # weak MUTUAL EXCITATION between the half-centers (conditional)
        if G["rg_weak_exc"] > 0.0:
            n.add_connection(_ssyn(G["rg_weak_exc"], exc=True), rg_e, rg_f)
            n.add_connection(_ssyn(G["rg_weak_exc"], exc=True), rg_f, rg_e)
        # descending drive: ANALOG DRIVE cell -> graded synapse onto the
        # spiking RG (rate-map; calibrated k_ns)
        n.add_connection(_ns2s(G["descend_to_rg_e"], exc=True), "DRIVE", rg_e)
        n.add_connection(_ns2s(G["descend_to_rg_f"], exc=True), "DRIVE", rg_f)
        n.add_connection(_ns2s(G["posture_to_rg_e"], exc=True), "POSTURE", rg_e)

        # v11 mechanosensory stance feedback (conditional; same structure)
        if self.stance_fb or self.aff_loops:
            heel_in = f"HEEL_{side}"
            toe_in = f"TOE_{side}"
            lbin = f"LBIN_{side}"
            self._add(heel_in, TAU["preset"], n, kind="aff")
            self._add(toe_in, TAU["preset"], n, kind="aff")
            self._add(lbin, TAU["ib_exc"], n)
            n.add_input(heel_in)
            self.inputs.append("HEEL_c_" + side)
            n.add_input(toe_in)
            self.inputs.append("TOE_c_" + side)
            n.add_input(lbin)
            self.inputs.append("LOAD_c_" + side)
            if G["full_rules"] > 0.0:
                n.add_connection(_ssyn(G["heel_rge"], exc=True),
                                 heel_in, ine)
                n.add_connection(_ssyn(G["heel_rge"], exc=False),
                                 heel_in, inf)
                n.add_connection(_ssyn(G["toe_rge"], exc=True),
                                 toe_in, ine)
            else:
                n.add_connection(_ssyn(G["heel_rge"], exc=True),
                                 heel_in, rg_e)
                n.add_connection(_ssyn(G["heel_rge"] * PHASE_RESET.get(
                    "inh", 1.0), exc=False), heel_in, rg_f)
                n.add_connection(_ssyn(G["toe_rge"], exc=True),
                                 toe_in, rg_e)
            n.add_connection(_ssyn(G["ib_rge"], exc=True), lbin, rg_e)
            if G["heel_in_f_exc"] > 0.0:
                n.add_connection(_ssyn(G["heel_in_f_exc"], exc=True),
                                 heel_in, inf)

        # v11b semi-closed sensory loops (conditional; AFF relays spiking)
        if self.aff_loops:
            aff_e = f"AFF_E_{side}"
            aff_f = f"AFF_F_{side}"
            self._add(aff_e, TAU["afferent"], n, kind="aff")
            self._add(aff_f, TAU["afferent"], n, kind="aff")
            n.add_input(aff_e)
            self.inputs.append("AFF_E_" + side)
            n.add_input(aff_f)
            self.inputs.append("AFF_F_" + side)
            n.add_connection(_ssyn(G["aff_e_rg"], exc=True), aff_e, rg_e)
            n.add_connection(_ssyn(G["aff_f_rg"], exc=True), aff_f, rg_f)

    def _build_pf(self, n: Network, side: str):
        if self.joint_pf:
            self._build_pf_layers(n, side)
            return
        rg_of = {"E1": f"RG_E_{side}", "E2": f"RG_E_{side}",
                 "F1": f"RG_F_{side}", "F2": f"RG_F_{side}"}
        for phase in PF_PHASES:
            tau_m, _ = PF_SHAPE[phase]
            pf = f"PF_{phase}_{side}"
            self._add(pf, TAU["pf"] * tau_m, n, kind="pf")
            n.add_connection(_ssyn(G["rg_to_pf"], exc=True),
                             rg_of[phase], pf)

        pf_in_e = f"PF_IN_E_{side}"
        pf_in_f = f"PF_IN_F_{side}"
        self._add(pf_in_e, TAU["pf"], n)
        self._add(pf_in_f, TAU["pf"], n)
        for ph in ("E1", "E2"):
            n.add_connection(_ssyn(G["pf_recip_inh"], exc=True),
                             f"PF_{ph}_{side}", pf_in_e)
        for ph in ("F1", "F2"):
            n.add_connection(_ssyn(G["pf_recip_inh"], exc=True),
                             f"PF_{ph}_{side}", pf_in_f)
        for ph in ("F1", "F2"):
            n.add_connection(_ssyn(G["pf_recip_inh"], exc=False),
                             pf_in_e, f"PF_{ph}_{side}")
        for ph in ("E1", "E2"):
            n.add_connection(_ssyn(G["pf_recip_inh"], exc=False),
                             pf_in_f, f"PF_{ph}_{side}")

        # mechanosensors ride the extensor central pathway (conditional)
        if (self.stance_fb or self.aff_loops) and G["ib_e_central"] > 0.0:
            for src in (f"HEEL_{side}", f"TOE_{side}"):
                n.add_connection(_ssyn(G["ib_e_central"], exc=True),
                                 src, f"PF_E1_{side}")
                n.add_connection(_ssyn(G["ib_e_central"], exc=True),
                                 src, f"PF_E2_{side}")
                n.add_connection(_ssyn(G["ib_e_central"], exc=True),
                                 src, f"InE_{side}")

        if self.aff_loops:
            aff_e = f"AFF_E_{side}"
            aff_f = f"AFF_F_{side}"
            for ph in ("E1", "E2"):
                n.add_connection(_ssyn(G["aff_e_pf"], exc=True),
                                 aff_e, f"PF_{ph}_{side}")
            for ph in ("F1", "F2"):
                n.add_connection(_ssyn(G["aff_f_pf"], exc=True),
                                 aff_f, f"PF_{ph}_{side}")

    def _build_pf_layers(self, n: Network, side: str):
        """Joint-layer PF mirror (G['joint_pf'] > 0): identical structure."""
        for hc in JPF_HCS:
            tau_m, _ = PF_SHAPE["E1" if hc.endswith("E") else "F1"]
            pf = f"PF_{hc}_{side}"
            self._add(pf, TAU["pf"] * tau_m, n, kind="pf")
            n.add_connection(_ssyn(G["rg_to_pf"], exc=True),
                             f"RG_{'E' if hc.endswith('E') else 'F'}_{side}",
                             pf)
        pf_in_e = f"PF_IN_E_{side}"
        pf_in_f = f"PF_IN_F_{side}"
        self._add(pf_in_e, TAU["pf"], n)
        self._add(pf_in_f, TAU["pf"], n)
        for hc in JPF_HCS:
            n.add_connection(_ssyn(G["pf_recip_inh"], exc=True),
                             f"PF_{hc}_{side}",
                             pf_in_e if hc.endswith("E") else pf_in_f)
        for hc in JPF_HCS:
            n.add_connection(_ssyn(G["pf_recip_inh"], exc=False),
                             pf_in_e if hc.endswith("F")
                             else pf_in_f, f"PF_{hc}_{side}")
        if (self.stance_fb or self.aff_loops) and G["ib_e_central"] > 0.0:
            for src in (f"HEEL_{side}", f"TOE_{side}"):
                for hc in ("HIP-E", "KNEE-E", "ANK-E"):
                    n.add_connection(_ssyn(G["ib_e_central"], exc=True),
                                     src, f"PF_{hc}_{side}")
                n.add_connection(_ssyn(G["ib_e_central"], exc=True),
                                 src, f"InE_{side}")
        if self.aff_loops:
            aff_e = f"AFF_E_{side}"
            aff_f = f"AFF_F_{side}"
            for hc in JPF_HCS:
                if hc.endswith("E"):
                    n.add_connection(_ssyn(G["aff_e_pf"], exc=True),
                                     aff_e, f"PF_{hc}_{side}")
                else:
                    n.add_connection(_ssyn(G["aff_f_pf"], exc=True),
                                     aff_f, f"PF_{hc}_{side}")
        # per-PF-layer contact variant (conditional; mirror)
        if self.stance_fb and G["heel_pf_layer"] > 0.0 and \
                f"HEEL_{side}" in self.idx:
            n.add_connection(_ssyn(G["heel_pf_layer"], exc=True),
                             f"HEEL_{side}", pf_in_e)
        if self.stance_fb and G["toe_df_inh"] > 0.0 and \
                f"TOE_{side}" in self.idx:
            toedf = f"TOEDF_{side}"
            if toedf not in self.idx:
                self._add(toedf, TAU["preset"], n)
            n.add_connection(_ssyn(G["toe_df_inh"], exc=True),
                             f"TOE_{side}", toedf)
            n.add_connection(_ssyn(G["toe_df_inh"], exc=False),
                             toedf, f"PF_ANK-F_{side}")

    def _add_muscle_neurons(self, n: Network, act: str, mi: MuscleInfo):
        mn, ia, ii, ib = (f"MN_{act}", f"Ia_{act}", f"II_{act}", f"Ib_{act}")
        self.mn_names[act] = mn
        self.aff_names[act] = {"Ia": ia, "II": ii, "Ib": ib}
        # MN: NON-SPIKING, unchanged analog frame (the hybrid ruling)
        n.add_neuron(_neu(TAU["mn"] * (1.0 + 0.5 * mi.biarticular)),
                     name=mn)
        self.idx[mn] = len(self.idx)
        # afferent encoders: SPIKING (plan table row 1)
        for name, tau in ((ia, TAU["afferent"]), (ii, TAU["afferent"]),
                          (ib, 2.0 * TAU["afferent"])):
            self._add(name, tau, n, kind="aff")
        if self.renshaw:
            self._add(f"RC_{act}", TAU["mn"], n, kind="rc")
        n.add_input(mn)
        self.inputs.append("POST_" + act)
        for port, name in (("Ia", ia), ("II", ii), ("Ib", ib)):
            n.add_input(name)
            self.inputs.append(port + "_" + act)

    def _wire_muscle(self, n: Network, act: str, mi: MuscleInfo):
        mn, ia, ii, ib = (f"MN_{act}", f"Ia_{act}", f"II_{act}", f"Ib_{act}")

        # ---- PF -> MN (spike synapses onto the analog MN; calibrated) --
        if self.joint_pf:
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
                        _ssyn(G["pf_to_mn"] * w, exc=True, onto_analog=True),
                        f"PF_{hc}_{mi.side}", mn)
        else:
            for phase in PF_PHASES:
                w = _group_weight(mi, W_PF_MN[phase])
                if w > 0.0:
                    for side in self.sides:
                        if mi.side == side:
                            n.add_connection(
                                _ssyn(G["pf_to_mn"] * w, exc=True,
                                      onto_analog=True),
                                f"PF_{phase}_{side}", mn)
        # POSTURE -> MN stays graded analog->analog (unchanged)
        w_post = _group_weight(mi, W_POSTURE, POSTURE_OVERRIDE)
        if w_post > 0.0:
            n.add_connection(_syn(G["posture_to_mn"] * w_post, exc=True),
                             "POSTURE", mn)

        # ---- proprioceptive pathways (full_rules lamination mirror) ----
        fr = G["full_rules"] > 0.0
        n.add_connection(_ssyn(G["ia_to_mn"], exc=True, onto_analog=True),
                         ia, mn)
        if fr:
            iix = f"IIX_{act}"
            if iix not in self.idx:
                self._add(iix, TAU["afferent"], n)
                n.add_connection(_ssyn(1.0, exc=True), ii, iix)
            n.add_connection(_ssyn(G["ii_to_mn"], exc=True, onto_analog=True),
                             iix, mn)
            ibin = f"IBIN_{act}"
            if ibin not in self.idx:
                self._add(ibin, TAU["afferent"], n)
                n.add_connection(_ssyn(1.0, exc=True), ib, ibin)
            n.add_connection(_ssyn(G["ib_to_mn_inh"], exc=False,
                                   onto_analog=True), ibin, mn)
            for act2, mi2 in self.muscles.items():
                if mi2.side == mi.side and act2 != act and \
                        mi2.groups[0] in ANTAGONIST.get(mi.groups[0], ()) \
                        and f"IBIN_{act2}" in self.idx:
                    n.add_connection(_ssyn(0.5, exc=False), ibin,
                                     f"IBIN_{act2}")
        else:
            n.add_connection(_ssyn(G["ii_to_mn"], exc=True, onto_analog=True),
                             ii, mn)
            n.add_connection(_ssyn(G["ib_to_mn_inh"], exc=False,
                                   onto_analog=True), ib, mn)

        # ---- per-muscle afferent -> CENTRAL feedback (mirror) ----------
        grp = mi.groups[0]
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
                n.add_connection(_ssyn(G["ib_e_central"], exc=True),
                                 ib, tgt)
        if grp in ("hip_flex", "knee_flex", "ankle_df", "hip_add",
                   "trunk_flex"):
            for tgt in (*tgt_f, f"RG_F_{mi.side}", f"InF_{mi.side}"):
                if G["ia_f_central"] > 0.0:
                    n.add_connection(_ssyn(G["ia_f_central"], exc=True),
                                     ia, tgt)
                if G["ii_f_central"] > 0.0:
                    n.add_connection(_ssyn(G["ii_f_central"], exc=True),
                                     ii, tgt)
            if G["ia_f_contra_f"] > 0.0 and self.interleg:
                contra = "l" if mi.side == "r" else "r"
                n.add_connection(_ssyn(G["ia_f_contra_f"], exc=False),
                                 ia, f"RG_F_{contra}")
        if G["ii_e_central"] > 0.0 and grp in EXTENSOR_STANCE_GROUPS:
            for tgt in (*tgt_e, f"RG_E_{mi.side}", f"InE_{mi.side}"):
                n.add_connection(_ssyn(G["ii_e_central"], exc=True),
                                 ii, tgt)

        # ---- Renshaw recurrent inhibition (mirror; graded MN->RC edge =
        # the plan's "rate encoder", RC->MN spike inhibition) -----------
        if self.renshaw:
            rc = f"RC_{act}"
            g_r = G["renshaw"]
            n.add_connection(_ns2s(1.0, exc=True), mn, rc)
            n.add_connection(_ssyn(g_r, exc=False, onto_analog=True), rc, mn)
            for act2, mi2 in self.muscles.items():
                if mi2.side == mi.side and act2 != act \
                        and f"RC_{act2}" in self.idx:
                    n.add_connection(_ssyn(g_r, exc=False), rc,
                                     f"RC_{act2}")

        # Ia reciprocal inhibition through IaIN (mirror)
        if self.ia_in:
            iain = f"IaIN_{act}"
            if iain not in self.idx:
                self._add(iain, TAU["afferent"], n)
                n.add_connection(_ssyn(G["ia_to_mn"], exc=True), ia, iain)
                gate_f1 = (tgt_f[0] if self.joint_pf
                           else f"PF_F1_{mi.side}")
                n.add_connection(_ssyn(0.5, exc=True), gate_f1, iain)
                if self.renshaw and f"RC_{act}" in self.idx:
                    n.add_connection(_ssyn(G["renshaw"], exc=False),
                                     f"RC_{act}", iain)
        if fr:
            iiin = f"IIIN_{act}"
            if iiin not in self.idx:
                self._add(iiin, TAU["afferent"], n)
                n.add_connection(_ssyn(1.0, exc=True), ii, iiin)
        for ant in ANTAGONIST.get(mi.groups[0], ()):
            for act2, mi2 in self.muscles.items():
                if mi2.side == mi.side and mi2.groups[0] == ant:
                    if self.ia_in:
                        n.add_connection(
                            _ssyn(G["ia_to_antagonist"], exc=False,
                                  onto_analog=True),
                            f"IaIN_{act}", f"MN_{act2}")
                        if fr and f"IaIN_{act2}" in self.idx:
                            n.add_connection(_ssyn(0.5, exc=False),
                                             f"IaIN_{act}",
                                             f"IaIN_{act2}")
                    else:
                        n.add_connection(
                            _ssyn(G["ia_to_antagonist"], exc=False,
                                  onto_analog=True),
                            ia, f"MN_{act2}")
                    if fr and f"IIIN_{act}" in self.idx:
                        n.add_connection(
                            _ssyn(G["ia_to_antagonist"], exc=False,
                                  onto_analog=True),
                            f"IIIN_{act}", f"MN_{act2}")

        # per-PF-layer flexion-afferent variant (conditional; mirror)
        if self.joint_pf and mi.groups[0] in ("hip_flex", "knee_flex",
                                              "ankle_df"):
            f_hc = (f"PF_{JPF_GROUP2HC[mi.groups[0]][0]}_{mi.side}")
            if G["ia_pf_f"] > 0.0 and self.ia_in \
                    and f"IaIN_{act}" in self.idx:
                n.add_connection(_ssyn(G["ia_pf_f"], exc=True),
                                 f"IaIN_{act}", f_hc)
            if G["ii_pf_f"] > 0.0 and fr and f"IIX_{act}" in self.idx:
                n.add_connection(_ssyn(G["ii_pf_f"], exc=True),
                                 f"IIX_{act}", f_hc)

        # ---- stance load sharing (extensor groups only; mirror) -------
        if mi.groups[0] in EXTENSOR_STANCE_GROUPS:
            grp = f"IBEXC_{mi.groups[0]}"
            gname = f"{grp}_{mi.side}"
            if gname not in self.idx:
                self._add(gname, TAU["ib_exc"], n)
                n.add_connection(_ssyn(1.0, exc=True),
                                 f"RG_E_{mi.side}", gname)
            n.add_connection(_ssyn(G["ib_group_exc"], exc=True), ib, gname)
            n.add_connection(_ssyn(G["ib_exc_to_mn"], exc=True,
                                   onto_analog=True), gname, mn)
            if self.stance_fb:
                lb = f"LBIN_{mi.side}"
                if lb not in self.idx:
                    self._add(lb, TAU["ib_exc"], n)
                n.add_connection(_ssyn(0.5, exc=True), gname, lb)
            if mi.groups[0] not in self.ib_exc_groups.get(mi.side, ()):
                self.ib_exc_groups[mi.side] = \
                    self.ib_exc_groups.get(mi.side, ()) + (mi.groups[0],)

        # ---- swing-knee/ankle PF suppression via KINH (mirror) --------
        if self.f1_kneext_inh and mi.groups[0] == "knee_ext":
            kname = f"KINH_{mi.side}"
            gate = (f"PF_KNEE-F_{mi.side}" if self.joint_pf
                    else f"PF_F1_{mi.side}")
            if kname not in self.idx:
                self._add(kname, TAU["ib_exc"], n)
                n.add_connection(_ssyn(1.5, exc=True), gate, kname)
                if G["contra_kinh"] > 0.0:
                    other = "l" if mi.side == "r" else "r"
                    if f"HEEL_{other}" in self.idx:
                        n.add_connection(
                            _ssyn(G["contra_kinh"], exc=True),
                            f"HEEL_{other}", kname)
            n.add_connection(_ssyn(G["f1_kneext_inh"], exc=False,
                                   onto_analog=True), kname, mn)
        if self.f1_kneext_inh and mi.groups[0] == "ankle_pf":
            kname = f"KINH_{mi.side}"
            gate = (f"PF_ANK-F_{mi.side}" if self.joint_pf
                    else f"PF_F1_{mi.side}")
            if kname not in self.idx:
                self._add(kname, TAU["ib_exc"], n)
                n.add_connection(_ssyn(1.5, exc=True), gate, kname)
            n.add_connection(_ssyn(G["f1_anklepf_inh"], exc=False,
                                   onto_analog=True), kname, mn)

    def _wire_balance(self, n: Network):
        """Balance inputs reach ankle + hip MNs (ALL graded analog->analog;
        unchanged from the non-spiking build - the BAL/VEST cells stay
        analog)."""
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
                n.add_connection(_syn(G["bal_trunk"], exc=True),
                                 "BAL_TRK_EXT", f"MN_{act}")
            elif g == "trunk_flex":
                n.add_connection(_syn(G["bal_trunk"], exc=True),
                                 "BAL_TRK_FLX", f"MN_{act}")

        if self.vest:
            for act, mi in self.muscles.items():
                vsrc = f"VEST_{mi.side}"
                g = mi.groups[0]
                if G["vest_ext"] > 0.0 and g in ("knee_ext", "ankle_pf",
                                                 "hip_ext", "trunk_ext"):
                    n.add_connection(_syn(G["vest_ext"], exc=True),
                                     vsrc, f"MN_{act}")
                elif G["vest_flex_inh"] > 0.0 and g in (
                        "hip_flex", "knee_flex", "ankle_df", "trunk_flex"):
                    n.add_connection(_syn(G["vest_flex_inh"], exc=False),
                                     vsrc, f"MN_{act}")

    # ------------------------------------------------------------------ run
    def compile(self, dt: float = DT):
        """Compile at the NETWORK dt (0.5 ms), not the plant dt."""
        self.compiled = self.net.compile(dt=dt, backend="numpy")
        # clamp the float_max theta the compiler gives non-spiking neurons
        # (silences an overflow warning in the theta update; spikes for
        # those neurons remain impossible: 1e6 >> any reachable V)
        _ns_mask = np.array(
            [not self.net.populations[i]["type"].params.get("spiking", False)
             for i in range(len(self.net.populations))])
        for attr in ("theta_0", "theta", "theta_last"):
            arr = getattr(self.compiled, attr, None)
            if arr is not None:
                arr[_ns_mask] = 1.0e6
        assert self.compiled.V.shape[0] == len(self.idx), \
            "neuron index bookkeeping mismatch"
        self.spike_counts = np.zeros(len(self.idx))
        return self.compiled

    def step(self, input_currents: np.ndarray) -> np.ndarray:
        """Advance ONE plant step (DT) by running n_sub network sub-steps
        at a constant input vector (plan rule 3: sub-stepping inside the
        2 ms plant step).  Returns the membrane-potential vector; RG/PF
        entries are the READOUT analog levels (see module docstring)."""
        u = list(input_currents)
        for _ in range(self.n_sub):
            self.compiled.forward(u)
            self.spike_counts += (self.compiled.spikes == -1)
        return self.compiled.V

    # ------------------------------------------------------------------ helpers
    def input_index(self, port: str) -> int:
        return self.inputs.index(port)

    def make_inputs(self) -> np.ndarray:
        return np.zeros(len(self.inputs))


def build(model_actuators: list[str], dt: float = DT,
          interleg: bool = True) -> SpinalNetworkSpiking:
    """Classify actuators and build + compile the SPIKING mirror network.
    Same signature/semantics as build_network.build; selected via
    AARL_NET=spiking."""
    muscles: dict[str, MuscleInfo] = {}
    for act in model_actuators:
        mi = classify(act)
        if mi is None:
            raise ValueError(f"actuator {act!r} not in muscle_map")
        muscles[act] = mi
    sides = tuple(sorted({mi.side for mi in muscles.values()}))
    net = SpinalNetworkSpiking(muscles=muscles, sides=sides,
                               interleg=interleg)
    n_sub = max(1, int(round(dt / 0.0005)))
    net.n_sub = n_sub
    net.compile(dt=dt / n_sub)
    return net
