"""Spinal circuit as IDIOMATIC SNS-Toolbox subnetworks (Tutorial 4
pattern: each layer is a Network subclass; the parent composes them with
add_network + parent-level connections and renders with
sns_toolbox.renderer).

Mirrors build_network.py exactly (same neurons, conductances from
params.G) for the REPRESENTATIVE 4-muscle build (knee ext/fl x 2 sides,
v5 PRESET + v6 KINH + Renshaw gains > 0). This is the source of the
toolbox-rendered diagram; extend the muscle list for the full 92.
"""
from __future__ import annotations

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.networks import Network

import params as P

E_HI = P.E_HI
C_ORANGE, C_BLUE = "orange", "cornflowerblue"
C_GREEN, C_RED = "lightgreen", "lightpink"
C_AFF = "gold"


def _neu(tau: float) -> NonSpikingNeuron:
    return NonSpikingNeuron(membrane_capacitance=float(tau),
                            membrane_conductance=1.0,
                            resting_potential=0.0, bias=0.0)


def _syn(g: float, exc: bool) -> NonSpikingSynapse:
    return NonSpikingSynapse(max_conductance=float(g),
                             reversal_potential=E_HI if exc else -E_HI,
                             e_lo=0.0, e_hi=E_HI)


class RhythmGeneratorNetwork(Network):
    """RG-E/F half-center + ADAP + PRESET INs, one side."""

    def __init__(self, side: str = "r"):
        super().__init__(name=f"Rhythm Generator {side.upper()}")
        G = P.G
        n = self
        n.add_population(_neu(P.TAU["rg"]), shape=[1], name=f"RG-E_{side}", color=C_ORANGE)
        n.add_population(_neu(P.TAU["rg"]), shape=[1], name=f"RG-F_{side}", color=C_BLUE)
        n.add_population(_neu(P.TAU["rg_adapt"]), shape=[1], name=f"ADAP-E_{side}")
        n.add_population(_neu(P.TAU["rg_adapt"]), shape=[1], name=f"ADAP-F_{side}")
        n.add_population(_neu(P.TAU["preset"]), shape=[1], name=f"PRESET-E_{side}",
                         color=C_RED)
        n.add_population(_neu(P.TAU["preset"]), shape=[1], name=f"PRESET-F_{side}",
                         color=C_RED)
        n.add_connection(_syn(G["rg_mutual_inh"], False), f"RG-E_{side}",
                         f"RG-F_{side}")
        n.add_connection(_syn(G["rg_mutual_inh"], False), f"RG-F_{side}",
                         f"RG-E_{side}")
        n.add_connection(_syn(G["rg_adapt_inh"], True), f"RG-E_{side}",
                         f"ADAP-E_{side}")
        n.add_connection(_syn(G["rg_adapt_inh"], False), f"ADAP-E_{side}",
                         f"RG-E_{side}")
        n.add_connection(_syn(G["rg_adapt_inh"], True), f"RG-F_{side}",
                         f"ADAP-F_{side}")
        n.add_connection(_syn(G["rg_adapt_inh"], False), f"ADAP-F_{side}",
                         f"RG-F_{side}")
        n.add_input(f"RG-E_{side}", name=f"DRIVE->E_{side}")
        n.add_input(f"RG-F_{side}", name=f"DRIVE->F_{side}")
        n.add_input(f"RG-E_{side}", name=f"POST->E_{side}")
        n.add_input(f"PRESET-E_{side}", name=f"HIP_EXT_{side}")
        n.add_input(f"PRESET-F_{side}", name=f"HIP_FLEX_{side}")
        n.add_output(f"RG-E_{side}", name=f"RG-E out {side}")


class PatternFormationNetwork(Network):
    """4 phase-window cells + PFA adaptation INs + KINH, one side."""

    def __init__(self, side: str = "r"):
        super().__init__(name=f"Pattern Formation {side.upper()}")
        G = P.G
        n = self
        for ph, col in (("E1", C_ORANGE), ("E2", C_ORANGE),
                        ("F1", C_BLUE), ("F2", C_BLUE)):
            tau_m, tau_a = P.PF_SHAPE[ph]
            n.add_population(_neu(P.TAU["pf"] * tau_m), shape=[1],
                             name=f"PF_{ph}_{side}", color=col)
            n.add_population(_neu(P.TAU["pf_adapt"] * tau_a), shape=[1],
                             name=f"PFA_{ph}_{side}")
            n.add_connection(_syn(1.5, True), f"PF_{ph}_{side}",
                             f"PFA_{ph}_{side}")
            n.add_connection(_syn(1.5, False), f"PFA_{ph}_{side}",
                             f"PF_{ph}_{side}")
        for a, b in (("E2", "F1"), ("F2", "E1"), ("E1", "F1"), ("E2", "F2")):
            n.add_connection(_syn(G["pf_recip_inh"], False),
                             f"PF_{a}_{side}", f"PF_{b}_{side}")
            n.add_connection(_syn(G["pf_recip_inh"], False),
                             f"PF_{b}_{side}", f"PF_{a}_{side}")
        n.add_population(_neu(P.TAU["ib_exc"]), shape=[1], name=f"KINH_{side}",
                         color=C_RED)
        n.add_connection(_syn(1.5, True), f"PF_F1_{side}", f"KINH_{side}")
        n.add_input(f"PF_E1_{side}", name=f"RG-E->PF_{side}")
        n.add_input(f"PF_F1_{side}", name=f"RG-F->PF_{side}")
        n.add_input(f"PF_E1_{side}", name=f"DRIVE->PF_{side}")
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
            # Renshaw (Deng A6): MN->RC exc 1.0; RC->MN inh gain; RC<->RC
            n.add_connection(_syn(1.0, True), f"MN_{m}", f"RC_{m}")
            n.add_connection(_syn(G["renshaw"], False), f"RC_{m}",
                             f"MN_{m}")
            n.add_connection(_syn(G["renshaw"], False), f"RC_{m}",
                             f"RC_{m}")
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
        # Renshaw cross-talk between the two columns (each pair once:
        # wire ext->flx; flx->ext is the mirror synapse, same effect
        # when both present, so add only one to match build_network's
        # one-synapse-per-pair)
        n.add_connection(_syn(G["renshaw"], False),
                         f"RC_{self.pools['ext']}",
                         f"RC_{self.pools['flx']}")
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
