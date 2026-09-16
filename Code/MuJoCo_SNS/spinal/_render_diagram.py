"""Compose the spinal circuit from spinal_layers.py subnetworks (Tutorial
4 pattern) and render with the OFFICIAL sns_toolbox.renderer (graphviz)
-> the sns_diagram.png style, structure-driven.

Representative circuit: right-side RG + PF + knee motor column, plus the
left RG for the commissural coupling. Cross-layer edges are parent-level
connections between flattened populations (names carry side suffixes).
"""
import io
import os
import sys
from pathlib import Path

# graphviz executables install into the conda env's Library\bin but that
# dir is only on PATH when the env is activated - add it for this process
_env_bin = Path(sys.executable).parent / "Library" / "bin"
if _env_bin.is_dir():
    os.environ["PATH"] = str(_env_bin) + os.pathsep + os.environ["PATH"]

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.networks import Network
from sns_toolbox.renderer import render

import params as P
from spinal_layers import (MotorColumnNetwork, PatternFormationNetwork,
                           RhythmGeneratorNetwork, _neu, _syn, E_HI)

HERE = Path(__file__).parent


def main():
    G = P.G
    net = Network(name="gait2392 spinal (representative)")
    net.add_network(RhythmGeneratorNetwork("r"), color="orange")
    net.add_network(RhythmGeneratorNetwork("l"), color="orange")
    net.add_network(PatternFormationNetwork("r"), color="skyblue")
    net.add_network(MotorColumnNetwork("r"), color="lightgreen")

    def syn(g, exc):
        return _syn(g, exc)

    # descending inputs -> RG (via the RG layer's named input ports)
    net.add_population(_neu(P.TAU["descend"]), shape=[1], name="DRIVE",
                       color="khaki")
    net.add_population(_neu(P.TAU["descend"]), shape=[1], name="POSTURE",
                       color="khaki")
    net.add_input("DRIVE", name="DRIVE (MLR)")
    net.add_input("POSTURE", name="POSTURE")
    net.add_connection(syn(G["descend_to_rg_e"], True), "DRIVE", "RG-E_r")
    net.add_connection(syn(G["descend_to_rg_f"], True), "DRIVE", "RG-F_r")
    net.add_connection(syn(G["posture_to_rg_e"], True), "POSTURE",
                       "RG-E_r")
    # hip-signal ports -> PRESET INs (v5)
    net.add_input("PRESET-E_r", name="HIP_EXT_SIG_r")
    net.add_input("PRESET-F_r", name="HIP_FLEX_SIG_r")
    # commissural (r<->l)
    net.add_connection(syn(G["rg_mutual_inh"], False), "RG-F_r", "RG-F_l")
    net.add_connection(syn(G["rg_mutual_inh"], False), "RG-F_l", "RG-F_r")
    net.add_connection(syn(0.5 * G["rg_mutual_inh"], False), "RG-E_r",
                       "RG-E_l")
    net.add_connection(syn(0.5 * G["rg_mutual_inh"], False), "RG-E_l",
                       "RG-E_r")
    # RG -> PF (both halves; parent-level, flattened names)
    for ph, g in (("E1", "RG-E_r"), ("E2", "RG-E_r"), ("F1", "RG-F_r"),
                  ("F2", "RG-F_r")):
        net.add_connection(syn(G["rg_to_pf"], True), g, f"PF_{ph}_r")
    net.add_connection(syn(G["drive_to_pf"], True), "DRIVE", "PF_E1_r")
    # PF -> MN (knee weights from the live table)
    ext, flx = "knee_ext_r", "knee_flex_r"
    for ph, pool, key in (("E1", ext, "knee_ext"), ("F1", flx, "knee_flex"),
                          ("F2", ext, "knee_ext")):
        w = P.W_PF_MN[ph].get(key, 0.0)
        if w > 0:
            net.add_connection(syn(G["pf_to_mn"] * w, True), f"PF_{ph}_r",
                               f"MN_{pool}")
    # KINH -> extensor MN (v6)
    net.add_connection(syn(G["f1_kneext_inh"], False), "KINH_r",
                       f"MN_{ext}")
    # IBEXC stance gate <- RG-E (parent edge into the motor layer port)
    net.add_connection(syn(1.0, True), "RG-E_r", f"IBEXC_{ext}")
    # POSTURE -> MNs (tonic; the per-muscle POST_i bias in the full net)
    net.add_connection(syn(0.3, True), "POSTURE", f"MN_{ext}")
    net.add_connection(syn(0.3, True), "POSTURE", f"MN_{flx}")
    # afferent input ports already exist inside the motor layer
    # (Ia/II/Ib per muscle); the runner formats muscle L/Ldot/F into them.

    out = HERE / "figures" / "sns_diagram_spinal"
    render(net, view=False, save=True, filename=str(out), img_format="png")
    print(f"rendered {out}.png")


if __name__ == "__main__":
    main()
