"""Spinal circuit as IDIOMATIC SNS-Toolbox subnetworks (Tutorial 4
pattern: each layer is a Network subclass; the parent composes them with
add_network + parent-level connections and renders with
sns_toolbox.renderer).

Mirrors build_network.py exactly (same neurons, conductances from
params.G) for the REPRESENTATIVE 4-muscle build (knee ext/fl x 2 sides).
LAMINATED architecture (RG-E->InE->RG-F, PF-E->PF_IN_E->PF-F per
Shevtsova/Deng) with PERSISTENT-Na RG half-centers (Deng 2022 / fixed
tau_h semantics; the compiled dynamics require
build_network.SNS_NumpyFixedTau - this module is for RENDERING, use
build_network.py for simulation). This is the source of the
toolbox-rendered diagram; extend the muscle list for the full 92.
"""
from __future__ import annotations

import numpy as np

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import (NonSpikingNeuron,
                                 NonSpikingNeuronWithPersistentSodiumChannel)
from sns_toolbox.networks import Network

import params as P

E_HI = P.E_HI
C_EXT, C_FLEX = "cornflowerblue", "orangered"
C_GREEN, C_RED = "lightgreen", "lightpink"
C_AFF, C_IN = "gold", "plum"


def _neu(tau: float) -> NonSpikingNeuron:
    return NonSpikingNeuron(membrane_capacitance=float(tau),
                            membrane_conductance=1.0,
                            resting_potential=0.0, bias=0.0)


def _syn(g: float, exc: bool) -> NonSpikingSynapse:
    return NonSpikingSynapse(max_conductance=float(g),
                             reversal_potential=E_HI if exc else -E_HI,
                             e_lo=0.0, e_hi=E_HI)


class RhythmGeneratorNetwork(Network):
    """RG-E/F persistent-Na half-centers + InE/InF laminated mutual
    inhibition + weak mutual excitation G_W, one side (Deng 2022 /
    Shinohara 2025; NO direct inhibitory HC<->HC synapses; ADAP and
    PRESET/PREA retired 2026-09-16). NOTE: renders as plain circles -
    the renderer does not draw channel parameters; the NaP h-gate lives
    in the compiled dynamics (build_network.SNS_NumpyFixedTau)."""

    def __init__(self, side: str = "r"):
        super().__init__(name=f"Rhythm Generator {side.upper()}")
        G = P.G
        n = self
        nap = NonSpikingNeuronWithPersistentSodiumChannel(
            membrane_capacitance=P.TAU["rg"], membrane_conductance=1.0,
            resting_potential=0.0, bias=0.0,
            g_ion=np.array([P.NAP["g_ion"]]),
            e_ion=np.array([P.NAP["e_ion"]]),
            k_m=np.array([P.NAP["k_m"]]),
            slope_m=np.array([P.NAP["slope_m"]]),
            e_m=np.array([P.NAP["e_m"]]),
            k_h=np.array([P.NAP["k_h"]]),
            slope_h=np.array([P.NAP["slope_h"]]),
            e_h=np.array([P.NAP["e_h"]]),
            tau_max_h=np.array([P.TAU["rg_nap_h"]]))
        n.add_population(nap, shape=[1], name=f"RG-E_{side}", color=C_EXT)
        n.add_population(nap, shape=[1], name=f"RG-F_{side}", color=C_FLEX)
        n.add_population(_neu(P.TAU["rg"]), shape=[1], name=f"InE_{side}",
                         color=C_IN)
        n.add_population(_neu(P.TAU["rg"]), shape=[1], name=f"InF_{side}",
                         color=C_IN)
        # IN-laminated mutual inhibition: RG-E excites InE, InE inhibits
        # RG-F; RG-F excites InF, InF inhibits RG-E (never direct)
        n.add_connection(_syn(G["rg_mutual_inh"], True), f"RG-E_{side}",
                         f"InE_{side}")
        n.add_connection(_syn(G["rg_mutual_inh"], False), f"InE_{side}",
                         f"RG-F_{side}")
        n.add_connection(_syn(G["rg_mutual_inh"], True), f"RG-F_{side}",
                         f"InF_{side}")
        n.add_connection(_syn(G["rg_mutual_inh"], False), f"InF_{side}",
                         f"RG-E_{side}")
        # weak mutual excitation G_W (Deng 2022)
        if G["rg_weak_exc"] > 0.0:
            n.add_connection(_syn(G["rg_weak_exc"], True), f"RG-E_{side}",
                             f"RG-F_{side}")
            n.add_connection(_syn(G["rg_weak_exc"], True), f"RG-F_{side}",
                             f"RG-E_{side}")
        n.add_input(f"RG-E_{side}", name=f"DRIVE->E_{side}")
        n.add_input(f"RG-F_{side}", name=f"DRIVE->F_{side}")
        n.add_input(f"RG-E_{side}", name=f"POST->E_{side}")
        n.add_output(f"RG-E_{side}", name=f"RG-E out {side}")


class PatternFormationNetwork(Network):
    """4 phase-window cells + PF_IN_E/F laminated cross-reciprocal + KINH,
    one side (NO direct PF<->PF synapses; PFA self-adaptation loops
    REMOVED 2026-09-16 per Shevtsova/Deng - not in their PF layer)."""

    def __init__(self, side: str = "r"):
        super().__init__(name=f"Pattern Formation {side.upper()}")
        G = P.G
        n = self
        for ph, col in (("E1", C_EXT), ("E2", C_EXT),
                        ("F1", C_FLEX), ("F2", C_FLEX)):
            tau_m, _tau_a = P.PF_SHAPE[ph]
            n.add_population(_neu(P.TAU["pf"] * tau_m), shape=[1],
                             name=f"PF_{ph}_{side}", color=col)
        # IN-laminated cross-reciprocal: PF-E cells excite PF_IN_E, which
        # inhibits the PF-F cells; PF-F -> PF_IN_F -> PF-E (Shevtsova/Deng)
        n.add_population(_neu(P.TAU["pf"]), shape=[1], name=f"PF_IN_E_{side}",
                         color=C_IN)
        n.add_population(_neu(P.TAU["pf"]), shape=[1], name=f"PF_IN_F_{side}",
                         color=C_IN)
        for ph in ("E1", "E2"):
            n.add_connection(_syn(G["pf_recip_inh"], True),
                             f"PF_{ph}_{side}", f"PF_IN_E_{side}")
        for ph in ("F1", "F2"):
            n.add_connection(_syn(G["pf_recip_inh"], True),
                             f"PF_{ph}_{side}", f"PF_IN_F_{side}")
        for ph in ("F1", "F2"):
            n.add_connection(_syn(G["pf_recip_inh"], False),
                             f"PF_IN_E_{side}", f"PF_{ph}_{side}")
        for ph in ("E1", "E2"):
            n.add_connection(_syn(G["pf_recip_inh"], False),
                             f"PF_IN_F_{side}", f"PF_{ph}_{side}")
        n.add_population(_neu(P.TAU["ib_exc"]), shape=[1], name=f"KINH_{side}",
                         color=C_RED)
        n.add_connection(_syn(1.5, True), f"PF_F1_{side}", f"KINH_{side}")
        n.add_input(f"PF_E1_{side}", name=f"RG-E->PF_{side}")
        n.add_input(f"PF_F1_{side}", name=f"RG-F->PF_{side}")
        n.add_output(f"KINH_{side}", name=f"KINH out {side}")


class MotorColumnNetwork(Network):
    """One knee muscle column per muscle (MN + RC + Ia/II/Ib + IBEXC),
    extensor + flexor, one side."""

    def __init__(self, side: str = "r", ext: str = "knee_ext",
                 flx: str = "knee_flex"):
        super().__init__(name=f"Motor Circuit ({ext}/{flx}) {side.upper()}")
        G = P.G
        n = self
        self.pools = {}
        for label, grp in (("ext", ext), ("flx", flx)):
            m = f"{grp}_{side}"
            self.pools[label] = m
            n.add_population(_neu(P.TAU["mn"]), shape=[1], name=f"MN_{m}", color=C_GREEN)
            n.add_population(_neu(P.TAU["mn"]), shape=[1], name=f"RC_{m}",
                             color="lightgray")
            n.add_population(_neu(P.TAU["afferent"]), shape=[1], name=f"Ia_{m}",
                             color=C_AFF)
            n.add_population(_neu(P.TAU["afferent"]), shape=[1], name=f"II_{m}",
                             color=C_AFF)
            n.add_population(_neu(2 * P.TAU["afferent"]), shape=[1], name=f"Ib_{m}",
                             color=C_AFF)
            # Renshaw (Deng A6): MN->RC exc 1.0; RC->MN inh gain. NO
            # self-synapse - mutual inhibition is BETWEEN different RCs
            n.add_connection(_syn(1.0, True), f"MN_{m}", f"RC_{m}")
            n.add_connection(_syn(G["renshaw"], False), f"RC_{m}",
                             f"MN_{m}")
            # afferents -> homonymous MN
            n.add_connection(_syn(G["ia_to_mn"], True), f"Ia_{m}",
                             f"MN_{m}")
            n.add_connection(_syn(G["ii_to_mn"], True), f"II_{m}",
                             f"MN_{m}")
            n.add_connection(_syn(G["ib_to_mn_inh"], False), f"Ib_{m}",
                             f"MN_{m}")
            n.add_input(f"Ia_{m}", name=f"Ia in {m}")
            n.add_input(f"II_{m}", name=f"II in {m}")
            n.add_input(f"Ib_{m}", name=f"Ib in {m}")
            n.add_input(f"MN_{m}", name=f"POST {m}")
            n.add_output(f"MN_{m}", name=f"act {m}")
        # Renshaw cross-talk: MUTUAL between the two columns (Deng A6
        # RC->RC; Hultborn's mutual Renshaw inhibition)
        n.add_connection(_syn(G["renshaw"], False),
                         f"RC_{self.pools['ext']}",
                         f"RC_{self.pools['flx']}")
        n.add_connection(_syn(G["renshaw"], False),
                         f"RC_{self.pools['flx']}",
                         f"RC_{self.pools['ext']}")
        # Ia reciprocal inhibition of the antagonist: via a per-pool IaIN
        # (PF_F1 phase gate + RC->IaIN disinhibition) when G["ia_in"] > 0;
        # direct edge only at ia_in == 0 (v10 behavior)
        if G["ia_in"] > 0.0:
            for label, ant in (("ext", "flx"), ("flx", "ext")):
                m = self.pools[label]
                n.add_population(_neu(P.TAU["afferent"]), shape=[1],
                                 name=f"IaIN_{m}", color=C_IN)
                # afferent leg: Ia excites its IaIN (same gain as the
                # homonymous Ia->MN arc, mirroring build_network.py)
                n.add_connection(_syn(G["ia_to_mn"], True),
                                 f"Ia_{m}", f"IaIN_{m}")
                n.add_input(f"IaIN_{m}", name=f"PF-F1 gate {m}")
                if G["renshaw"] > 0.0:
                    n.add_connection(_syn(G["renshaw"], False),
                                     f"RC_{m}", f"IaIN_{m}")
                n.add_connection(_syn(G["ia_to_antagonist"], False),
                                 f"IaIN_{m}",
                                 f"MN_{self.pools[ant]}")
        else:
            # Ia reciprocal (direct/lumped; Deng routes via IaIN)
            n.add_connection(_syn(G["ia_to_antagonist"], False),
                             f"Ia_{self.pools['ext']}",
                             f"MN_{self.pools['flx']}")
            n.add_connection(_syn(G["ia_to_antagonist"], False),
                             f"Ia_{self.pools['flx']}",
                             f"MN_{self.pools['ext']}")
        # Ib load sharing (extensor stance group only)
        m = self.pools["ext"]
        n.add_population(_neu(P.TAU["ib_exc"]), shape=[1], name=f"IBEXC_{m}",
                         color=C_GREEN)
        n.add_connection(_syn(G["ib_group_exc"], True), f"Ib_{m}",
                         f"IBEXC_{m}")
        n.add_connection(_syn(G["ib_exc_to_mn"], True), f"IBEXC_{m}",
                         f"MN_{m}")
        n.add_input(f"IBEXC_{m}", name=f"RG-E gate {m}")
