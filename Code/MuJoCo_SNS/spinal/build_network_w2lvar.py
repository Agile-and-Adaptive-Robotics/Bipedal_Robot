"""VARIANT 1 (2026-09-25 campaign, goal 4): the Walker_2_Layer_CPG layout
on the gait2392 body — w2lvar.

Per side:
  * ONE rhythm generator: persistent-Na half-centers RG-E/RG-F + laminated
    InE/InF — a verbatim copy of build_network.SpinalNetwork._build_rg
    (same params.NAP set, same FIXED tau_h via SNS_NumpyFixedTau).
  * TWO pattern-formation pairs instead of the parent's 4 phase cells /
    6 joint half-centers:
      - HIP pair   PF_HIP-E / PF_HIP-F   (drives hip_ext/hip_abd/hip_add /
        hip_flex)
      - KNEE+ANKLE synergy pair  PF_KNEE-E / PF_KNEE-F  — ONE cell drives
        BOTH knee and ankle muscle groups (the W2L biarticular-gas pattern):
        KNEE-E -> knee_ext + ankle_pf, KNEE-F -> knee_flex + ankle_df
        (w2l_cpg/README.md muscle-map table + parent JPF_GROUP2HC).
    Each pair carries its own IN lamination at Ben's drawing value 2.749
    (E cell -> IN -> inhibit partner F cell, both directions).
  * Per-muscle MNs stay per-muscle (all 92 actuators), with the
    Shinohara-style AUTogenic per-muscle afferent expansion: each muscle's
    own Ia/II/Ib encoder projects centrally (flexor Ia/II -> the F centers,
    extensor Ib -> the E centers; Ben's reading rule: "proprioception feeds
    back to ITSELF, autogenic per muscle — never group-broadcast").

Dressing per Ben's three rule files:
  * Circuit_rules (master, 90n/99e): heel = stance reset AT the PF layer
    (heel IN -> ipsi InE exc 0.5 + CONTRA InF exc 0.5 + the ipsilateral
    PF-layer E-lamination INs exc 0.5); toe = dorsiflexion inhibition ONLY
    (toe IN -> TOEDF exc 5 -> inh the dorsiflexion-driving HC at 2.749);
    Ib load unchanged (per-muscle Ib -> RG-E/InE/PF-E at 0.5); Renshaw;
    Ia homonymous + reciprocal via IaIN; II exc/inh disynaptic; Ib
    autogenic inhibition + stance reversal.
  * ben_shevtsova_20260924.json: the inter-RG wiring is the Shevtsova
    commissural set (replaces the parent's C1/V3 block):
      RG-F -> V2a -> V0V -> contralateral INI -> inh contra RG-F;
      RG-F -> V0D -> inh contra RG-F;
      RG-E -> V3-E -> small exc of contra RG-E + exc of contra InE;
      brainstem alpha (DRIVE) inhibits V0V/V0D (speed-dependent crossing).
  * ben_shinohara_20260924.json: per-muscle autogenic afferent expansion;
    duplicated drawing gains = the Ia vs II subpopulations (we reuse the
    stack's proven ia_to_mn / ii_to_mn calibration instead).

Selection: default-inert. The DEFAULT build is untouched; this module is
entered only when env AARL_NET=w2lvar (see the switch at the top of
build_network.build). Nothing here reads or writes any other module's
state except params.G gains listed in VARIANT_G (applied only in variant
runs — a variant process is dedicated; documented in the goal-4 report).

New-file guarantee: importing this module has no side effects; building
the DEFAULT network never executes any line of it.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from muscle_map import EXTENSOR_STANCE_GROUPS, MuscleInfo, classify
from params import (AFF, DT, E_HI, G, MOD, NAP, PF_SHAPE, PHASE_RESET,
                    POSTURE_OVERRIDE, TAU, W_PF_MN, W_POSTURE)

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.neurons import NonSpikingNeuronWithPersistentSodiumChannel
from sns_toolbox.networks import Network

import build_network as BN
from build_network import (ANTAGONIST, SNS_NumpyFixedTau, _group_weight,
                           _jpf_weights, _neu, _syn)

HERE = Path(__file__).parent

# ---------------------------------------------------------------------------
# Variant gains (module-local; NOT in params.G — see report §deviations).
# Relative magnitudes follow Ben's rule files; absolute conductances are
# calibrated on the stack's 0..5 mV toolbox scale per the w2l_cpg README
# convention ("gain mapping is a calibration, not a unit conversion").
GAINS = dict(
    # ---- PF layer ----
    # RG -> PF drive. Ben's drawing lists 0.1, but that is the SynAmp-scale
    # value of the contact-driven AnimatLab reference where PF-E ALSO gets
    # 0.5 of live load afference; on our scale the identical template edge
    # was calibrated to G["rg_to_pf"]=2.4 (w2l_cpg README pf_drive row) and
    # a literal 0.1 leaves PF cells sub-threshold (measured, gate-c probe).
    w2l_rg_to_pf=2.4,
    w2l_pf_recip=2.749,    # per-pair IN lamination (drawing literal)
    w2l_pf_to_mn=2.0,      # = parent G["pf_to_mn"] (stack-proven)
    # ---- contact rules (master rules file) ----
    w2l_heel_pf=0.5,       # heel IN -> the two ipsilateral PF-layer E-INs
    w2l_toe_df=5.0,        # toe IN -> TOEDF (drawing: toe g 5)
    w2l_toedf_inh=2.749,   # TOEDF -> inh dorsiflexion HC (drawing 2.749)
    # ---- per-muscle autogenic central projections (Shinohara eq10/11;
    # drawing rows SN-Ia/SN-II -> F trio and Ib grp -> E trio). These are
    # the drawing's AGGREGATE one-node conductances — the builder divides
    # by the per-side family count when it duplicates the edge per
    # muscle (autogenic expansion; measured co-latch fix 2026-09-25).
    w2l_ib_central=0.5,    # extensor Ib -> RG-E + InE + the two E PF cells
    w2l_ia_central=0.53,   # flexor Ia -> RG-F + InF + the two F PF cells
    w2l_ii_central=0.53,   # flexor II -> the same F trio
    # ---- Shevtsova commissural set (ben_shevtsova_20260924.json).
    # JSON relative gains x 4.0 (the w2l_cpg calibration: relative 1.0 ->
    # G["rg_mutual_inh"]=4.0), EXCEPT the antiphase/sync balance: at the
    # raw json ratios the V3-E sync leg (V3-E -> contra InE) overpowered
    # the crossed F inhibition and both legs stepped IN phase
    # (r(RG_E_r,RG_E_l)=+0.86, smoke2 2026-09-25); the crossed-F edges
    # were raised so the set locks ANTIPHASE (all four cell classes kept).
    w2l_comm=1.0,          # master switch for the whole crossed set
    w2l_v2a=4.0,           # RG-F -> V2a        (json 1.0)
    w2l_v0v=4.0,           # V2a -> V0V         (json 1.0)
    w2l_ini=2.4,           # V0V -> contra INI  (json 0.6)
    w2l_ini_inh=2.0,       # INI -> contra RG-F (json 0.075, raised for
                           #  antiphase lock)
    w2l_v0d=2.8,           # RG-F -> V0D        (json 0.7)
    w2l_v0d_inh=1.2,       # V0D -> contra RG-F (json 0.07, raised)
    w2l_v3e=1.4,           # RG-E -> V3-E       (json 0.35)
    w2l_v3e_contra=0.04,   # V3-E -> contra RG-E small exc (json 0.02)
    w2l_v3e_inE=1.0,       # V3-E -> contra InE exc (json 1.0, reduced:
                           #  the sync leg)
    w2l_alpha=2.0,         # DRIVE -> inh V0V/V0D (json alpha 0.5)
)

# params.G overlay applied (only) in variant runs. These are the
# runner-contract gains: the runner scales the HEEL_c/TOE_c port currents
# by G["heel_rge"]/G["toe_rge"] (runner.py:1387-1388), selects the neuro
# watch columns by G["joint_pf"] (runner.py:1145), and the RG builder
# semantics; values are starting points for the variant curriculum.
VARIANT_G = dict(
    joint_pf=1.0,        # runner watch names (PF_HIP-E_r / PF_KNEE-E_r /
                         # PF_KNEE-F_r / PF_ANK-F_r) + parent flag parity
    heel_rge=0.5,        # port-current scale; edges are wired variant-side
    toe_rge=0.5,
    ib_rge=0.0,          # LBIN port exists (runner contract) but inert;
                         # the drawing's Ib->RG-E is per-muscle
                         # (w2l_ib_central), not the group LBIN
    ia_in=0.5,           # IaIN population present (rules: reciprocal via
                         # IaIN, not a direct Ia->antagonist edge)
    renshaw=0.5,         # per-muscle RC (rules rc motif)
    full_rules=0.0,      # (inert here: the variant wires its own motifs)
    f1_kneext_inh=0.0,   # KINH swing suppression available, default off
    f1_anklepf_inh=0.0,
)

# ---------------------------------------------------------------------------
# W2L synergy PF mapping.
W2L_HCS = ("HIP-E", "HIP-F", "KNEE-E", "KNEE-F")

# functional group -> synergy half-centers (anatomical fallback when the
# fitted table has no row; trunk deliberately ABSENT = parent joint-layer
# parity: trunk rides POSTURE/BAL only, the W2L walker has no trunk).
W2L_GROUP2HC = {
    "hip_ext": ("HIP-E",), "hip_abd": ("HIP-E",),
    "hip_add": ("HIP-E", "HIP-F"), "hip_flex": ("HIP-F",),
    "knee_ext": ("KNEE-E",), "knee_flex": ("KNEE-F",),
    "ankle_pf": ("KNEE-E",), "ankle_df": ("KNEE-F",),
}

# biarticular crosses carried over from the parent's fitted JPF_CROSS
# (all targets exist in the 4-HC set; gas = the biarticular-gas pattern).
W2L_CROSS = {"rect_fem": ("HIP-F",), "semimem": ("HIP-E",),
             "semiten": ("HIP-E",), "bifemsh": ("HIP-E",),
             "gas_med": ("KNEE-F",), "gas_lat": ("KNEE-F",),
             "grac": ("HIP-F",), "sart": ("KNEE-F",)}

# fitted-table rows (parent joint_pf_weights.json) that feed each synergy
# cell: the merged KNEE+ANKLE cell delivers the KNEE-E row to knee muscles
# and the ANK-E row to ankle muscles (each muscle keeps its fitted weight).
W2L_CELL_ROWS = {"HIP-E": ("HIP-E",), "HIP-F": ("HIP-F",),
                 "KNEE-E": ("KNEE-E", "ANK-E"),
                 "KNEE-F": ("KNEE-F", "ANK-F")}


# ---------------------------------------------------------------------------
@dataclass
class W2LVarnet(BN.SpinalNetwork):
    """The w2lvar network: parent dataclass fields, variant topology."""

    # ---------------------------------------------------------------- build
    def __post_init__(self):
        self.net = Network(name="gait2392 spinal w2lvar")
        n = self.net
        # variant flags (runner contract: net.vest / net.stance_fb /
        # net.aff_loops are read by the runner to decide which ports to
        # feed). The variant ALWAYS has heel/toe mechanosensors (the
        # master rules' stance reset); no AFF relays, no VEST cells.
        self.stance_fb = True
        self.aff_loops = False
        self.vest = False
        self.joint_pf = True          # parity flag (runner reads G, not this)
        self.phase_reset = False
        self.f1_kneext_inh = bool(G["f1_kneext_inh"] > 0.0
                                  or G["f1_anklepf_inh"] > 0.0)
        self.renshaw = True           # rules rc motif, gain G["renshaw"]
        self.ia_in = True             # rules: reciprocal via IaIN
        self.ib_exc_groups = {}
        self._aliases = {}            # runner-watch name -> real neuron idx

        # Per-muscle decomposition of the drawing's SYMBOLIC afferent
        # nodes: SN-Ia/SN-II ("flexor ii/ia") and "extensor Ib" are ONE
        # aggregate node each per side in ben_shinohara_20260924.json /
        # ben_rules_20260924.json, so their g 0.5 / 0.53 edges are
        # AGGREGATE conductances. The autogenic expansion duplicates the
        # edge per muscle, so each per-muscle edge carries the aggregate
        # gain divided by the family's per-side muscle count (measured
        # co-latch 2026-09-25: broadcasting the aggregate from every
        # muscle pinned BOTH RG half-centers at ~3.0 mV — the w2l
        # co-latch signature).
        self._nE = {}                 # side -> #stance-extensor muscles
        self._nF = {}                 # side -> #flexor-family muscles
        for s in self.sides:
            self._nE[s] = sum(1 for mi in self.muscles.values()
                              if mi.side == s
                              and mi.groups[0] in EXTENSOR_STANCE_GROUPS)
            self._nF[s] = sum(1 for mi in self.muscles.values()
                              if mi.side == s
                              and mi.groups[0] in ("hip_flex", "knee_flex",
                                                   "ankle_df", "hip_add",
                                                   "trunk_flex"))

        # ---- descending / balance cells (shared; the runner drives all
        # of these ports unconditionally) — identical to the parent ----
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

        # ---- per-side circuitry (RG + the W2L PF pairs) ----
        for side in self.sides:
            self._build_rg_w2l(n, side)
            self._build_pf_w2l(n, side)

        # ---- muscles: per-muscle neurons first, then wire ----
        self._order = {a: i for i, a in enumerate(self.muscles)}
        for act, mi in self.muscles.items():
            self._add_muscle_neurons(n, act, mi)
        for act, mi in self.muscles.items():
            self._wire_muscle_w2l(n, act, mi)
        self._wire_balance(n)          # inherited verbatim (BAL_* -> MNs)

        # ---- Shevtsova commissural set (both sides exist now) ----
        if self.interleg and GAINS["w2l_comm"] > 0.0:
            self._build_commissural(n)
            # master-rules heel crossed leg: heel IN -> CONTRA InF EXC
            # (suppresses the partner's flexor center = crossed stance
            # support; drawing rows heel -> IN-InF_2813, g 0.5)
            for side in self.sides:
                other = "l" if side == "r" else "r"
                if f"HEEL_{side}" in self.idx and G["heel_rge"] > 0.0:
                    n.add_connection(_syn(G["heel_rge"], exc=True),
                                     f"HEEL_{side}", f"InF_{other}")

        # ---- runner-watch aliases (LOGGING only; no neurons, no
        # synapses). runner.py logs PF_KNEE-F_r/l and PF_ANK-F_r/l by
        # name; our dorsiflexion drive rides the merged KNEE-F cell, so
        # those names point at the same real neuron. Captured before
        # compile() so the neuron-count assert stays exact.
        self._n_true = len(self.idx)
        for side in self.sides:
            self._aliases[f"PF_ANK-F_{side}"] = self.idx[f"PF_KNEE-F_{side}"]
            self._aliases[f"PF_ANK-E_{side}"] = self.idx[f"PF_KNEE-E_{side}"]
        self.idx.update(self._aliases)

    # ------------------------------------------------------------------ RG
    def _build_rg_w2l(self, n: Network, side: str):
        """The parent _build_rg pattern, verbatim mechanics: NaP
        half-centers (params.NAP, FIXED tau_h) + InE/InF lamination +
        DRIVE/POSTURE edges — plus the variant's HEEL/TOE/LBIN mechanosensor
        INs and input ports (the runner's stance_fb contract)."""
        rg_e, rg_f = f"RG_E_{side}", f"RG_F_{side}"
        ine, inf = f"InE_{side}", f"InF_{side}"
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

        # IN-laminated mutual inhibition (Shevtsova/Deng A6): NO direct
        # RG<->RG inhibitory synapses.
        n.add_connection(_syn(G["rg_mutual_inh"], exc=True), rg_e, ine)
        n.add_connection(_syn(G["rg_mutual_inh"], exc=False), ine, rg_f)
        n.add_connection(_syn(G["rg_mutual_inh"], exc=True), rg_f, inf)
        n.add_connection(_syn(G["rg_mutual_inh"], exc=False), inf, rg_e)

        # weak mutual excitation (Deng G_W escape assist), conditional
        if G["rg_weak_exc"] > 0.0:
            n.add_connection(_syn(G["rg_weak_exc"], exc=True), rg_e, rg_f)
            n.add_connection(_syn(G["rg_weak_exc"], exc=True), rg_f, rg_e)
        # descending drive + posture stance bias (brainstem gamma/alpha
        # analog, shevtsova json rows "brainstem gamma/alpha -> RG")
        n.add_connection(_syn(G["descend_to_rg_e"], exc=True), "DRIVE", rg_e)
        n.add_connection(_syn(G["descend_to_rg_f"], exc=True), "DRIVE", rg_f)
        n.add_connection(_syn(G["posture_to_rg_e"], exc=True), "POSTURE", rg_e)

        # ---- mechanosensor INs + ports (runner stance_fb contract:
        # runner.py:1377-1389 feeds HEEL_c/TOE_c/LOAD_c scaled by
        # G["heel_rge"]/G["toe_rge"]/G["ib_rge"]) ----
        heel_in, toe_in, lbin = f"HEEL_{side}", f"TOE_{side}", f"LBIN_{side}"
        self._add(heel_in, TAU["preset"], n)
        self._add(toe_in, TAU["preset"], n)
        self._add(lbin, TAU["ib_exc"], n)
        n.add_input(heel_in)
        self.inputs.append("HEEL_c_" + side)
        n.add_input(toe_in)
        self.inputs.append("TOE_c_" + side)
        n.add_input(lbin)
        self.inputs.append("LOAD_c_" + side)
        # master-rules heel = stance reset:
        #   ipsi InE EXC (strengthens E's suppression of ipsi F)
        n.add_connection(_syn(G["heel_rge"], exc=True), heel_in, ine)
        #   (the CONTRA InF EXC edge is wired in __post_init__ after both
        #   sides exist — the partner's InF is not built yet here)
        #   LBIN -> RG-E kept from the parent contract (gain 0 = inert
        #   edge; the drawing's Ib->RG-E is per-muscle, w2l_ib_central)
        n.add_connection(_syn(G["ib_rge"], exc=True), lbin, rg_e)
        # NOTE: the heel -> PF-layer E-IN edges live in _build_pf_w2l
        # (the PF INs are created there — same ordering constraint that
        # caused the parent's 2026-09-16 stage-2 crash).

    # ------------------------------------------------------------------ PF
    def _build_pf_w2l(self, n: Network, side: str):
        """TWO PF pairs: HIP (E/F) and the KNEE+ANKLE synergy (E/F), each
        with its own IN lamination at the drawing's 2.749."""
        rg_of = {"HIP-E": f"RG_E_{side}", "KNEE-E": f"RG_E_{side}",
                 "HIP-F": f"RG_F_{side}", "KNEE-F": f"RG_F_{side}"}
        for hc in W2L_HCS:
            tau_m, _ = PF_SHAPE["E1" if hc.endswith("E") else "F1"]
            pf = f"PF_{hc}_{side}"
            self._add(pf, TAU["pf"] * tau_m, n)
            n.add_connection(_syn(GAINS["w2l_rg_to_pf"], exc=True),
                             rg_of[hc], pf)

        # per-pair IN lamination (drawing: HC-PF-E_x -> IN-PF_x ->
        # inh HC-PF-F_x, g 2.749 both directions; per-layer INs, not a
        # shared pair-wide IN)
        for pair, (e_hc, f_hc) in (("HIP", ("HIP-E", "HIP-F")),
                                   ("KNEE", ("KNEE-E", "KNEE-F"))):
            in_e = f"PF_IN_{pair}-E_{side}"
            in_f = f"PF_IN_{pair}-F_{side}"
            self._add(in_e, TAU["pf"], n)
            self._add(in_f, TAU["pf"], n)
            n.add_connection(_syn(GAINS["w2l_pf_recip"], exc=True),
                             f"PF_{e_hc}_{side}", in_e)
            n.add_connection(_syn(GAINS["w2l_pf_recip"], exc=False),
                             in_e, f"PF_{f_hc}_{side}")
            n.add_connection(_syn(GAINS["w2l_pf_recip"], exc=True),
                             f"PF_{f_hc}_{side}", in_f)
            n.add_connection(_syn(GAINS["w2l_pf_recip"], exc=False),
                             in_f, f"PF_{e_hc}_{side}")

        # ---- heel = stance reset AT the PF layer (master rules): the
        # heel IN excites the two E-side lamination INs (drawing rows
        # heel -> IN-PF_2810/2812/2814, g 0.5; two pairs -> two INs).
        if f"HEEL_{side}" in self.idx and GAINS["w2l_heel_pf"] > 0.0:
            for pair in ("HIP", "KNEE"):
                n.add_connection(_syn(GAINS["w2l_heel_pf"], exc=True),
                                 f"HEEL_{side}", f"PF_IN_{pair}-E_{side}")

        # ---- toe = dorsiflexion inhibition ONLY (drawing: toe IN ->
        # IN-PF_dorsiflexion_inhibit g 5 -> inh HC-PF-Dorsiflexion g
        # 2.749). The dorsiflexion drive rides the merged KNEE-F cell in
        # this variant, so the inhibition lands there (documented
        # synergy side-effect: knee-flex drive dips with it).
        if f"TOE_{side}" in self.idx and GAINS["w2l_toe_df"] > 0.0:
            toedf = f"TOEDF_{side}"
            if toedf not in self.idx:
                self._add(toedf, TAU["preset"], n)
            n.add_connection(_syn(GAINS["w2l_toe_df"], exc=True),
                             f"TOE_{side}", toedf)
            n.add_connection(_syn(GAINS["w2l_toedf_inh"], exc=False),
                             toedf, f"PF_KNEE-F_{side}")

    # ------------------------------------------------------------- commissural
    def _build_commissural(self, n: Network):
        """Shevtsova set (ben_shevtsova_20260924.json), replacing the
        parent's C1/V3 block. Per side a: V2a_a, V0V_a, V0D_a, V3E_a;
        per side b (target): INI_b. Edges (a->b both directions):
          RG_F_a -> V2a_a -> V0V_a -> INI_b -> inh RG_F_b   (json 1/1/0.6/0.075)
          RG_F_a -> V0D_a -> inh RG_F_b                     (json 0.7/0.07)
          RG_E_a -> V3E_a -> exc RG_E_b (small) + exc InE_b (json .35/.02/1)
          DRIVE  -> inh V0V_s + inh V0D_s  (brainstem alpha, json 0.5)
        """
        gs = GAINS
        for a, b in (("r", "l"), ("l", "r")):
            v2a, v0v, v0d, v3e = (f"V2a_{a}", f"V0V_{a}", f"V0D_{a}",
                                  f"V3E_{a}")
            ini = f"INI_{b}"
            for name, tau in ((v2a, TAU["pf"]), (v0v, TAU["pf"]),
                              (v0d, TAU["pf"]), (v3e, TAU["pf"]),
                              (ini, TAU["rg"])):
                if name not in self.idx:
                    self._add(name, tau, n)
            # V2a -> V0V -> INI -> contra RG-F
            n.add_connection(_syn(gs["w2l_v2a"], exc=True),
                             f"RG_F_{a}", v2a)
            n.add_connection(_syn(gs["w2l_v0v"], exc=True), v2a, v0v)
            n.add_connection(_syn(gs["w2l_ini"], exc=True), v0v, ini)
            n.add_connection(_syn(gs["w2l_ini_inh"], exc=False), ini,
                             f"RG_F_{b}")
            # V0D direct crossed F inhibition
            n.add_connection(_syn(gs["w2l_v0d"], exc=True), f"RG_F_{a}", v0d)
            n.add_connection(_syn(gs["w2l_v0d_inh"], exc=False), v0d,
                             f"RG_F_{b}")
            # V3-E crossed extensor support
            n.add_connection(_syn(gs["w2l_v3e"], exc=True), f"RG_E_{a}", v3e)
            n.add_connection(_syn(gs["w2l_v3e_contra"], exc=True), v3e,
                             f"RG_E_{b}")
            n.add_connection(_syn(gs["w2l_v3e_inE"], exc=True), v3e,
                             f"InE_{b}")
            # brainstem alpha inhibits V0V/V0D (speed-dependent crossing;
            # guarded by its own gain, conditional-topology rule)
            if gs["w2l_alpha"] > 0.0:
                n.add_connection(_syn(gs["w2l_alpha"], exc=False), "DRIVE",
                                 v0v)
                n.add_connection(_syn(gs["w2l_alpha"], exc=False), "DRIVE",
                                 v0d)

    # ------------------------------------------------------------- per muscle
    def _add_muscle_neurons(self, n: Network, act: str, mi: MuscleInfo):
        """Identical to the parent, except RCs are ALWAYS built (rules rc
        motif; gain = G["renshaw"], 0 leaves the RC silent-but-present —
        variant topology, not the default build, so the bit-exact-at-0
        contract does not apply here)."""
        mn, ia, ii, ib = (f"MN_{act}", f"Ia_{act}", f"II_{act}", f"Ib_{act}")
        self.mn_names[act] = mn
        self.aff_names[act] = {"Ia": ia, "II": ii, "Ib": ib}
        for name, tau in ((mn, TAU["mn"] * (1.0 + 0.5 * mi.biarticular)),
                          (ia, TAU["afferent"]), (ii, TAU["afferent"]),
                          (ib, 2.0 * TAU["afferent"])):
            self._add(name, tau, n)
        self._add(f"RC_{act}", TAU["mn"], n)
        n.add_input(mn)
        self.inputs.append("POST_" + act)
        for port, name in (("Ia", ia), ("II", ii), ("Ib", ib)):
            n.add_input(name)
            self.inputs.append(port + "_" + act)

    def _w2l_pf_weight(self, mi: MuscleInfo, base: str, cell: str,
                       jw: dict) -> float:
        """PF->MN weight of one synergy cell onto one muscle: the fitted
        row (parent joint_pf_weights.json) summed over the cell's source
        rows, else the anatomical fallback (1.0 when the muscle's group
        maps to the cell, plus the biarticular cross list)."""
        w = 0.0
        for row in W2L_CELL_ROWS[cell]:
            v = jw.get(row, {}).get(base)
            if v is not None:
                w += v
        if w <= 0.0:
            w = 1.0 if (cell in W2L_GROUP2HC.get(mi.groups[0], ())
                        or cell in W2L_CROSS.get(base, ())) else 0.0
        return w

    def _wire_muscle_w2l(self, n: Network, act: str, mi: MuscleInfo):
        """Per-muscle wiring: W2L synergy PF->MN routing + the rule-file
        per-muscle motifs (Ia homo + IaIN reciprocal, II disynaptic exc +
        IIIN antagonist inh, Ib autogenic inh + stance-reversal exc,
        Renshaw, Shinohara autogenic central projections)."""
        mn, ia, ii, ib = (f"MN_{act}", f"Ia_{act}", f"II_{act}", f"Ib_{act}")

        # ---- PF -> MN (synergy cells; fitted weights where available) --
        base = act.rsplit("_", 1)[0]
        jw = _jpf_weights()
        for cell in W2L_HCS:
            w = self._w2l_pf_weight(mi, base, cell, jw)
            if w > 0.0:
                n.add_connection(_syn(GAINS["w2l_pf_to_mn"] * w, exc=True),
                                 f"PF_{cell}_{mi.side}", mn)
        # POSTURE tonic (parent table + per-muscle overrides)
        w_post = _group_weight(mi, W_POSTURE, POSTURE_OVERRIDE)
        if w_post > 0.0:
            n.add_connection(_syn(G["posture_to_mn"] * w_post, exc=True),
                             "POSTURE", mn)

        # ---- proprioceptive autogenic pathways (rules motifs) ----
        # Ia homonymous monosynaptic excitation (rules ia_homo)
        n.add_connection(_syn(G["ia_to_mn"], exc=True), ia, mn)
        # Ia -> IaIN -> antagonist MN (rules ia_recip; Deng A6 phase gate
        # from the muscle's own F synergy cell; RC->IaIN disinhibition)
        iain = f"IaIN_{act}"
        self._add(iain, TAU["afferent"], n)
        n.add_connection(_syn(G["ia_to_mn"], exc=True), ia, iain)
        gate_f = (f"PF_HIP-F_{mi.side}" if mi.groups[0] in
                  ("hip_ext", "hip_abd", "hip_add", "hip_flex")
                  else f"PF_KNEE-F_{mi.side}")
        n.add_connection(_syn(0.5, exc=True), gate_f, iain)
        if self.renshaw:
            n.add_connection(_syn(G["renshaw"], exc=False),
                             f"RC_{act}", iain)
        # II -> IIX (exc IN) -> same MN (rules ii_exc)
        iix = f"IIX_{act}"
        self._add(iix, TAU["afferent"], n)
        n.add_connection(_syn(1.0, exc=True), ii, iix)
        n.add_connection(_syn(G["ii_to_mn"], exc=True), iix, mn)
        # II -> IIIN (inh IN) -> antagonist MNs (rules ii_inh)
        iiin = f"IIIN_{act}"
        self._add(iiin, TAU["afferent"], n)
        n.add_connection(_syn(1.0, exc=True), ii, iiin)
        # Ib -> IBIN (inh IN) -> same MN (rules ib_auto)
        ibin = f"IBIN_{act}"
        self._add(ibin, TAU["afferent"], n)
        n.add_connection(_syn(1.0, exc=True), ib, ibin)
        n.add_connection(_syn(G["ib_to_mn_inh"], exc=False), ibin, mn)

        # ---- Shinohara AUTogenic central projections (per muscle, to the
        # IPSILATERAL centers; reading rule: never group-broadcast). The
        # gains are the drawing's AGGREGATE node values divided by the
        # per-side family count (see __post_init__) so the sum over the
        # family reproduces the drawing's one-node conductance. ----
        grp = mi.groups[0]
        tgt_e = (f"PF_HIP-E_{mi.side}", f"PF_KNEE-E_{mi.side}",
                 f"RG_E_{mi.side}", f"InE_{mi.side}")
        tgt_f = (f"PF_HIP-F_{mi.side}", f"PF_KNEE-F_{mi.side}",
                 f"RG_F_{mi.side}", f"InF_{mi.side}")
        if GAINS["w2l_ib_central"] > 0.0 and grp in EXTENSOR_STANCE_GROUPS:
            g_ib = GAINS["w2l_ib_central"] / max(self._nE[mi.side], 1)
            for tgt in tgt_e:
                n.add_connection(_syn(g_ib, exc=True), ib, tgt)
        if grp in ("hip_flex", "knee_flex", "ankle_df", "hip_add",
                   "trunk_flex"):
            if GAINS["w2l_ia_central"] > 0.0:
                g_ia = GAINS["w2l_ia_central"] / max(self._nF[mi.side], 1)
                for tgt in tgt_f:
                    n.add_connection(_syn(g_ia, exc=True), ia, tgt)
            if GAINS["w2l_ii_central"] > 0.0:
                g_ii = GAINS["w2l_ii_central"] / max(self._nF[mi.side], 1)
                for tgt in tgt_f:
                    n.add_connection(_syn(g_ii, exc=True), ii, tgt)

        # ---- Renshaw recurrent inhibition (rules rc motif; gain
        # G["renshaw"] — 0 leaves the cell silent) ----
        if self.renshaw:
            rc = f"RC_{act}"
            g_r = G["renshaw"]
            n.add_connection(_syn(1.0, exc=True), mn, rc)
            n.add_connection(_syn(g_r, exc=False), rc, mn)
            for act2, mi2 in self.muscles.items():
                if mi2.side == mi.side and act2 != act \
                        and f"RC_{act2}" in self.idx:
                    n.add_connection(_syn(g_r, exc=False), rc,
                                     f"RC_{act2}")

        # ---- antagonist routing: IaIN + IIIN -> antagonist MNs; mutual
        # IaIN<->IaIN and IBIN<->IBIN inhibition (rules) ----
        for ant in ANTAGONIST.get(mi.groups[0], ()):
            for act2, mi2 in self.muscles.items():
                if mi2.side == mi.side and mi2.groups[0] == ant:
                    n.add_connection(
                        _syn(G["ia_to_antagonist"], exc=False),
                        iain, f"MN_{act2}")
                    n.add_connection(
                        _syn(G["ia_to_antagonist"], exc=False),
                        iiin, f"MN_{act2}")
                    if f"IaIN_{act2}" in self.idx:
                        n.add_connection(_syn(0.5, exc=False), iain,
                                         f"IaIN_{act2}")
                    if f"IBIN_{act2}" in self.idx:
                        n.add_connection(_syn(0.5, exc=False), ibin,
                                         f"IBIN_{act2}")

        # ---- stance load sharing (rules ib_rev; parent IBEXC pattern) --
        if grp in EXTENSOR_STANCE_GROUPS:
            gname = f"IBEXC_{grp}_{mi.side}"
            if gname not in self.idx:
                self._add(gname, TAU["ib_exc"], n)
                n.add_connection(_syn(1.0, exc=True),
                                 f"RG_E_{mi.side}", gname)
            n.add_connection(_syn(G["ib_group_exc"], exc=True), ib, gname)
            n.add_connection(_syn(G["ib_exc_to_mn"], exc=True), gname, mn)
            if grp not in self.ib_exc_groups.get(mi.side, ()):
                self.ib_exc_groups[mi.side] = \
                    self.ib_exc_groups.get(mi.side, ()) + (grp,)

        # ---- KINH swing suppression (parent lever, gain-gated; the
        # swing gate is the merged KNEE-F synergy cell for both knee_ext
        # and ankle_pf pools) ----
        if self.f1_kneext_inh and grp in ("knee_ext", "ankle_pf"):
            kname = f"KINH_{mi.side}"
            if kname not in self.idx:
                self._add(kname, TAU["ib_exc"], n)
                n.add_connection(_syn(1.5, exc=True),
                                 f"PF_KNEE-F_{mi.side}", kname)
            g = (G["f1_kneext_inh"] if grp == "knee_ext"
                 else G["f1_anklepf_inh"])
            if g > 0.0:
                n.add_connection(_syn(g, exc=False), kname, mn)

    # ------------------------------------------------------------------ run
    def compile(self, dt: float = DT):
        self.compiled = self.net.compile(dt=dt, backend="numpy")
        # fixed-tau_h stepper for the RG persistent-Na h-gates (Deng
        # semantics — the parent patch, single source of truth)
        self.compiled.__class__ = SNS_NumpyFixedTau
        # alias-aware assert (idx carries runner-watch aliases; the true
        # neuron count was captured in __post_init__ before aliasing)
        assert self.compiled.V.shape[0] == self._n_true, \
            "neuron index bookkeeping mismatch"
        return self.compiled


# ---------------------------------------------------------------------------
def build(model_actuators: list[str], dt: float = DT,
          interleg: bool = True) -> W2LVarnet:
    """Build + compile the w2lvar network. Applies the VARIANT_G overlay
    to params.G (documented; the variant process is dedicated — the
    runner scales HEEL_c/TOE_c port currents by G["heel_rge"]/["toe_rge"]
    and picks the neuro watch columns by G["joint_pf"], so the overlay
    must stay in force for the run)."""
    for k, v in VARIANT_G.items():
        G[k] = float(v)
    muscles: dict[str, MuscleInfo] = {}
    for act in model_actuators:
        mi = classify(act)
        if mi is None:
            raise ValueError(f"actuator {act!r} not in muscle_map")
        muscles[act] = mi
    sides = tuple(sorted({mi.side for mi in muscles.values()}))
    net = W2LVarnet(muscles=muscles, sides=sides, interleg=interleg)
    net.compile(dt=dt)
    return net
