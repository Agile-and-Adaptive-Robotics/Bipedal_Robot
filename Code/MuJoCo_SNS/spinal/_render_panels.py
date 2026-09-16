"""Render each spinal subnetwork (RG / PF / motor column) separately with
the OFFICIAL sns_toolbox.renderer and compose the panels into one figure -
the tidy layered look of the toolbox docs, still 100% code-driven.

Usage: python _render_panels.py
"""
import io
import os
import sys
from pathlib import Path

_env_bin = Path(sys.executable).parent / "Library" / "bin"
if _env_bin.is_dir():
    os.environ["PATH"] = str(_env_bin) + os.pathsep + os.environ["PATH"]

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np
from PIL import Image, ImageDraw

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.networks import Network
from sns_toolbox.renderer import render

import params as P
from spinal_layers import (MotorColumnNetwork, PatternFormationNetwork,
                           RhythmGeneratorNetwork, _neu, _syn)

HERE = Path(__file__).parent
OUT = HERE / "figures"
G = P.G
# render the TUNED configuration (v7 winner) so edge labels are honest
G["f1_kneext_inh"] = 0.421
G["f1_anklepf_inh"] = 0.579
G["renshaw"] = 0.5


def add_drive(net, targets):
    """Descending-drive input node feeding the given layer inputs."""
    net.add_population(_neu(P.TAU["descend"]), shape=[1], name="DRIVE",
                       color="khaki")
    net.add_input("DRIVE", name="DRIVE (MLR)")
    for tgt, g in targets:
        net.add_connection(_syn(g, True), "DRIVE", tgt)


def render_layer(net, name):
    p = OUT / name
    render(net, view=False, save=True, filename=str(p), img_format="png")
    print(f"rendered {p}.png")
    return Path(str(p) + ".png")


def main():
    # ---- panel 1: RG layer (right) with descending + hip inputs ----
    rg = RhythmGeneratorNetwork("r")
    add_drive(rg, [("RG-E_r", G["descend_to_rg_e"]),
                   ("RG-F_r", G["descend_to_rg_f"])])
    rg.add_input("RG-E_r", name="POSTURE")
    p_rg = render_layer(rg, "sns_layer_rg")

    # ---- panel 2: PF layer (right) ----
    pf = PatternFormationNetwork("r")
    add_drive(pf, [("PF_E1_r", G["drive_to_pf"])])
    # RG drive shown as input ports at Deng's coupling strength
    pf.add_input("PF_E1_r", name=f"RG-E (g={G['rg_to_pf']:.2f})")
    pf.add_input("PF_F1_r", name=f"RG-F (g={G['rg_to_pf']:.2f})")
    p_pf = render_layer(pf, "sns_layer_pf")

    # ---- panel 3: knee motor column (right) ----
    mo = MotorColumnNetwork("r")
    mo.add_input("IBEXC_knee_ext_r", name=f"RG-E gate (1.0)")
    mo.add_population(_neu(P.TAU["descend"]), shape=[1], name="POSTURE",
                      color="khaki")
    mo.add_input("POSTURE", name="POSTURE + POST_i")
    mo.add_connection(_syn(0.3, True), "POSTURE", "MN_knee_ext_r")
    mo.add_connection(_syn(0.3, True), "POSTURE", "MN_knee_flex_r")
    # PF -> MN at the tuned weights
    for ph, pool, key in (("E1", "knee_ext_r", "knee_ext"),
                          ("E2", "knee_ext_r", "knee_ext"),
                          ("F1", "knee_flex_r", "knee_flex"),
                          ("F2", "knee_flex_r", "knee_flex")):
        w = P.W_PF_MN[ph].get(key, 0.0)
        if w > 0:
            mo.add_input(f"MN_{pool}", name=f"PF-{ph} (w={w:.3f})")
    mo.add_input("MN_knee_ext_r", name=f"KINH (g={G['f1_kneext_inh']:.2f})")
    p_mo = render_layer(mo, "sns_layer_motor")

    # ---- compose horizontally with a title strip ----
    imgs = [Image.open(p) for p in (p_rg, p_pf, p_mo)]
    scale = 1400 / max(im.height for im in imgs)
    imgs = [im.resize((max(1, int(im.width * scale)),
                       max(1, int(im.height * scale))), Image.LANCZOS)
            for im in imgs]
    pad = 24
    w = sum(im.width for im in imgs) + pad * 4
    h = max(im.height for im in imgs) + pad * 2 + 60
    canvas = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(canvas)
    labels = ["RHYTHM GENERATOR (right)", "PATTERN FORMATION (right)",
              "MOTOR CIRCUIT - knee column (right)"]
    x = pad
    for im, lab in zip(imgs, labels):
        canvas.paste(im, (x, pad + 40))
        d.text((x, pad + 6), lab, fill=(60, 60, 60))
        x += im.width + pad
    d.text((pad, h - 26),
           "gait2392 spinal SNS - subnetwork panels, rendered from the "
           "compiled Network objects (sns_toolbox.renderer). Left-right "
           "commissurals and the mirrored left layers omitted for "
           "clarity; conductances live from params.py.",
           fill=(120, 120, 120))
    out = OUT / "sns_diagram_panels.png"
    canvas.save(out)
    print(f"composed {out}")


if __name__ == "__main__":
    main()
