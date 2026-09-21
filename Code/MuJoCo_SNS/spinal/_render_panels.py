"""Render each spinal subnetwork (RG / PF / motor column) separately with
the OFFICIAL sns_toolbox.renderer and compose the panels into one figure -
the tidy layered look of the toolbox docs, still 100% code-driven.

Usage: python _render_panels.py
"""
import io
import os
import sys
import textwrap
from pathlib import Path

_env_bin = Path(sys.executable).parent / "Library" / "bin"
if _env_bin.is_dir():
    os.environ["PATH"] = str(_env_bin) + os.pathsep + os.environ["PATH"]

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.networks import Network
from sns_toolbox.renderer import render

import params as P
from spinal_layers import (MotorColumnNetwork, PatternFormationNetwork,
                           RhythmGeneratorNetwork, _neu, _syn)

HERE = Path(__file__).parent
OUT = HERE / "figures"
G = P.G


def tuned_gains():
    """Return gains and an honest source note for this render.

    A stage whose score is the -100 no-countable-cycles sentinel is not a
    winner and must not be presented as tuned. When stage 3 is invalid, use
    the successful stage-1 rhythm gains plus nonzero representative feedback
    gains so the conditional topology remains visible.
    """
    import json
    g = dict(f1_kneext_inh=0.6, f1_anklepf_inh=0.6,
             phase_reset_e=0.0, phase_reset_f=0.0, ia_in=0.6,
             heel_rge=0.6, toe_rge=0.4, ib_rge=0.6, renshaw=0.5,
             rg_weak_exc=0.4)
    note = "representative nonzero NaP-architecture gains"
    try:
        s1 = json.loads(
            (HERE / "curriculum_stage1.json").read_text("utf-8"))
        if float(s1.get("score", -100.0)) > -100.0:
            p = s1.get("params", {})
            if "desc_e" in p:
                g["descend_to_rg_e"] = float(p["desc_e"])
            if "desc_f" in p:
                g["descend_to_rg_f"] = float(p["desc_f"])
            if "rg_to_pf" in p:
                g["rg_to_pf"] = float(p["rg_to_pf"])
            if "rg_nap_h" in p:
                P.TAU["rg_nap_h"] = float(p["rg_nap_h"])
            note = "curriculum stage-1 winner + representative feedback gains"
            print("stage-1 rhythm winner applied: curriculum_stage1.json")
    except FileNotFoundError:
        print("no curriculum_stage1.json - using representative defaults")
    try:
        s3 = json.loads(
            (HERE / "curriculum_stage3.json").read_text("utf-8"))
        if float(s3.get("score", -100.0)) > -100.0:
            for k in ("phase_reset_e", "phase_reset_f", "heel_rge",
                      "toe_rge", "ib_rge", "ia_in"):
                if k in s3.get("params", {}):
                    g[k] = float(s3["params"][k])
            note = f"curriculum stage-3 winner (score {s3.get('score')})"
            print(f"stage-3 winner applied: curriculum_stage3.json "
                  f"(score {s3.get('score')})")
        else:
            note += f"; failed stage-3 score {s3.get('score')} ignored"
            print(f"ignored curriculum_stage3.json sentinel result "
                  f"(score {s3.get('score')})")
    except FileNotFoundError:
        print("no curriculum_stage3.json - laminated defaults only")
    return g, note


APPLIED, CONFIG_NOTE = tuned_gains()
G.update(APPLIED)
print("representative gains: " +
      " ".join(f"{k}={APPLIED[k]:.3f}" for k in sorted(APPLIED)))


def add_drive(net, targets):
    """Descending-drive input node feeding the given layer inputs."""
    net.add_population(_neu(P.TAU["descend"]), shape=[1], name="DRIVE",
                       color="khaki")
    net.add_input("DRIVE", name="DRIVE (MLR)")
    for tgt, g in targets:
        net.add_connection(_syn(g, True), "DRIVE", tgt)


def add_posture(net, targets):
    """Tonic POSTURE neuron used by the compiled network."""
    net.add_population(_neu(P.TAU["descend"]), shape=[1], name="POSTURE",
                       color="khaki")
    net.add_input("POSTURE", name="POSTURE (tonic)")
    for tgt, g in targets:
        net.add_connection(_syn(g, True), "POSTURE", tgt)


def render_layer(net, name):
    p = OUT / name
    render(net, view=False, save=True, filename=str(p), img_format="png")
    print(f"rendered {p}.png")
    return Path(str(p) + ".png")


def main():
    # ---- panel 1: RG layer (right) with descending + hip inputs ----
    rg = RhythmGeneratorNetwork("r")
    # The reusable RG class exposes direct composition ports.  This panel
    # instead draws the actual compiled DRIVE and POSTURE populations, so
    # suppress the otherwise duplicate direct-port invhouses.
    rg.inputs.clear()
    add_drive(rg, [("RG-E_r", G["descend_to_rg_e"]),
                   ("RG-F_r", G["descend_to_rg_f"])])
    add_posture(rg, [("RG-E_r", G["posture_to_rg_e"])])
    p_rg = render_layer(rg, "sns_layer_rg")

    # ---- panel 2: PF layer (right) ----
    pf = PatternFormationNetwork("r")
    # PF is driven by RG only. The old tonic DRIVE->PF edge was removed
    # from build_network.py on 2026-09-16 and must not appear here.
    for port in pf.inputs:
        if port["name"].startswith("RG-E"):
            port["name"] = f"RG-E (g={G['rg_to_pf']:.2f})"
        elif port["name"].startswith("RG-F"):
            port["name"] = f"RG-F (g={G['rg_to_pf']:.2f})"
    p_pf = render_layer(pf, "sns_layer_pf")

    # ---- panel 3: knee motor column (right) ----
    mo = MotorColumnNetwork("r")
    mo.add_input("IBEXC_knee_ext_r", name=f"RG-E gate (1.0)")
    mo.add_population(_neu(P.TAU["descend"]), shape=[1], name="POSTURE",
                      color="khaki")
    mo.add_input("POSTURE", name="POSTURE (tonic)")
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

    # ---- page-scale composition: RG and PF above, wider motor panel below ----
    # A single horizontal strip made the verified diagrams unreadably small
    # when fitted to a dissertation page.  Preserve the official renderer
    # output, but arrange it as two rows with explicit subpanel labels.
    rg, pf, mo = [Image.open(p).convert("RGB") for p in (p_rg, p_pf, p_mo)]
    content_w, pad, gap = 2700, 64, 48
    top_scale = (content_w - gap) / (rg.width + pf.width)
    rg = rg.resize((int(rg.width * top_scale), int(rg.height * top_scale)),
                   Image.LANCZOS)
    pf = pf.resize((int(pf.width * top_scale), int(pf.height * top_scale)),
                   Image.LANCZOS)
    motor_scale = content_w / mo.width
    mo = mo.resize((content_w, int(mo.height * motor_scale)), Image.LANCZOS)

    from matplotlib import font_manager
    font_path = font_manager.findfont("DejaVu Sans")
    title_font = ImageFont.truetype(font_path, 38)
    label_font = ImageFont.truetype(font_path, 29)
    note_font = ImageFont.truetype(font_path, 18)
    title_h, label_h, row_gap, footer_h = 58, 44, 56, 126
    top_h = max(rg.height, pf.height)
    w = content_w + 2 * pad
    h = (pad + title_h + label_h + top_h + row_gap + label_h +
         mo.height + footer_h + pad)
    canvas = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(canvas)
    d.text((pad, pad), "gait2392 spinal SNS — toolbox-native subnetworks",
           font=title_font, fill=(25, 25, 25))

    y_top_label = pad + title_h
    d.text((pad, y_top_label), "A   Rhythm generator (right)",
           font=label_font, fill=(55, 55, 55))
    x_pf = pad + rg.width + gap
    d.text((x_pf, y_top_label), "B   Pattern formation (right)",
           font=label_font, fill=(55, 55, 55))
    y_top = y_top_label + label_h
    canvas.paste(rg, (pad, y_top))
    canvas.paste(pf, (x_pf, y_top))

    y_motor_label = y_top + top_h + row_gap
    d.text((pad, y_motor_label), "C   Motor/reflex circuit — representative right-knee columns",
           font=label_font, fill=(55, 55, 55))
    y_motor = y_motor_label + label_h
    canvas.paste(mo, (pad, y_motor))

    note = (
        "Rendered with sns_toolbox.renderer from the Network objects. "
        "Laminated architecture: RG-E→InE→RG-F and "
        "PF-E→PF_IN_E→PF-F (no direct half-center inhibitory synapses). "
        "Mirrored left layers and left–right commissurals are omitted for "
        "clarity. Conductances: " + CONFIG_NOTE + ".")
    note_y = y_motor + mo.height + 28
    for line in textwrap.wrap(note, width=165):
        d.text((pad, note_y), line, font=note_font, fill=(100, 100, 100))
        note_y += 24
    out = OUT / "sns_diagram_panels.png"
    canvas.save(out)
    canvas.save(OUT / "sns_diagram_panels.pdf", resolution=300.0)
    print(f"composed {out} and {OUT / 'sns_diagram_panels.pdf'}")


if __name__ == "__main__":
    main()
