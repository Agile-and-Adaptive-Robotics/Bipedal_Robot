"""Reader-facing spinal-circuit schematic grounded in the literature.

Panel A follows the bilateral layered organization used by Rybak,
Shevtsova, Shinohara, and Deng: supraspinal drive -> RG -> PF -> MN ->
muscle, with c1/V3 commissural relays centered at the midline.
Panel B follows the mechanism-at-a-time presentation used by Deng and
Di Russo: one complete antagonist knee motif with no dangling interface
paths.  The dense ``circuit_dengstyle`` figure remains the audit view.

All displayed neural edge classes are asserted against a freshly compiled
representative SNS network.  Plant-interface paths are outside SNS-Toolbox's
connection list, so they have a separate explicit completeness contract.
"""
from __future__ import annotations

import argparse

from matplotlib.patches import FancyBboxPatch

import build_network as bn
import draw_circuit as dc
import params as P


def _compiled_representative():
    tg = dc._tuned_gains()
    gain_keys = (
        "phase_reset_e", "phase_reset_f", "f1_kneext_inh", "renshaw",
        "ia_in", "heel_rge", "toe_rge", "ib_rge", "rg_weak_exc",
        "ib_e_central", "ia_f_central", "ii_f_central", "ii_e_central",
        "ia_f_contra_f", "v3_to_ibexc", "c1_gain", "v3_gain",
    )
    saved = {k: P.G[k] for k in gain_keys}
    try:
        for k in gain_keys:
            P.G[k] = tg[k]
        net = bn.build(
            ["vas_lat_r", "semimem_r", "vas_lat_l", "semimem_l"],
            interleg=True,
        )
    finally:
        P.G.update(saved)
    return net, tg


def _panel_frame(cv, x0, x1, y0, y1, title):
    cv.ax.add_patch(FancyBboxPatch(
        (x0, y0), x1 - x0, y1 - y0,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        fc="white", ec="0.72", lw=1.0, zorder=0.1,
    ))
    cv.ax.text(x0 + 0.18, y1 - 0.18, title, ha="left", va="top",
               fontsize=10.5, fontweight="bold", color="0.18", zorder=8)


def make_figure(fmts):
    net, tg = _compiled_representative()
    compiled = set(dc._edge_groups(net))
    drawn = set()
    plant_drawn = set()

    def neural(cv, p, q, key, **kwargs):
        assert key in compiled, f"displayed neural edge absent from net: {key}"
        dc.syn(cv, p, q, key[2] == "exc", **kwargs)
        drawn.add(key)

    def plant(cv, p, q, key, **kwargs):
        dc.syn(cv, p, q, True, **kwargs)
        plant_drawn.add(key)

    cv = dc.Canvas(16.2, 10.4)
    _panel_frame(cv, 0.35, 10.25, 0.48, 10.0,
                 "A   Bilateral spinal architecture")
    _panel_frame(cv, 10.48, 15.87, 0.48, 10.0,
                 "B   Representative right-knee sensorimotor motif")

    # ------------------------------------------------------------------ A
    dc.layer_band(cv, 0.55, 10.05, 8.82, 9.58, dc.OI_SKY,
                  "SUPRASPINAL", fs=7.3)
    dc.layer_band(cv, 0.55, 10.05, 6.62, 8.68, dc.OI_ORANGE,
                  "RHYTHM GENERATOR", fs=7.3)
    dc.layer_band(cv, 0.55, 10.05, 4.55, 6.48, dc.OI_SKY,
                  "PATTERN FORMATION", fs=7.3)
    dc.layer_band(cv, 0.55, 10.05, 2.70, 4.41, dc.OI_GREEN,
                  "MOTOR POOLS", fs=7.3)
    dc.layer_band(cv, 0.55, 10.05, 0.72, 2.56, dc.OI_VERM,
                  "MUSCULOSKELETAL PLANT + SENSORY FEEDBACK", fs=7.3)
    cv.ax.plot([5.30, 5.30], [0.78, 8.63], color="0.62", lw=0.9,
               ls=(0, (3, 3)), zorder=1)
    cv.ax.text(2.75, 8.56, "LEFT", ha="center", va="top", fontsize=7.2,
               color="0.35", fontweight="bold")
    cv.ax.text(7.85, 8.56, "RIGHT", ha="center", va="top", fontsize=7.2,
               color="0.35", fontweight="bold")

    drive = dc.input_box(cv, 4.28, 9.17, 1.55, 0.42, "DRIVE",
                         "MLR-like command")
    posture = dc.input_box(cv, 6.35, 9.17, 1.55, 0.42, "POSTURE / BALANCE",
                           "tonic + feedback")

    sides = {}
    for side, cx in (("L", 2.75), ("R", 7.85)):
        ext_x, flx_x = cx - 0.78, cx + 0.78
        rge = dc.neuron(cv, ext_x, 7.83, "RG-E", "NaP",
                        fc=dc._tint(dc.W_E), r=0.34, fs=7.0)
        rgf = dc.neuron(cv, flx_x, 7.83, "RG-F", "NaP",
                        fc=dc._tint(dc.W_F), r=0.34, fs=7.0)
        ine = dc.neuron(cv, cx - 0.30, 6.98, "InE",
                        fc=dc._tint(dc.OI_PURPLE), r=0.19, fs=5.7)
        inf = dc.neuron(cv, cx + 0.30, 6.98, "InF",
                        fc=dc._tint(dc.OI_PURPLE), r=0.19, fs=5.7)
        neural(cv, rge, ine, ("RG-E", "RG-IN", "exc"),
               color=dc.W_E, lw=1.1)
        neural(cv, ine, rgf, ("RG-IN", "RG-F", "inh"),
               color=dc.W_INH, lw=1.1, rad=-0.10)
        neural(cv, rgf, inf, ("RG-F", "RG-IN", "exc"),
               color=dc.W_F, lw=1.1)
        neural(cv, inf, rge, ("RG-IN", "RG-E", "inh"),
               color=dc.W_INH, lw=1.1, rad=0.10)
        pfe = dc.input_box(cv, ext_x, 5.60, 1.12, 0.52,
                           "PF-E1 | PF-E2", "stance phases")
        pff = dc.input_box(cv, flx_x, 5.60, 1.12, 0.52,
                           "PF-F1 | PF-F2", "swing phases")
        pine = dc.neuron(cv, cx - 0.27, 4.90, "IN-E",
                         fc=dc._tint(dc.OI_PURPLE), r=0.16, fs=5.2)
        pinf = dc.neuron(cv, cx + 0.27, 4.90, "IN-F",
                         fc=dc._tint(dc.OI_PURPLE), r=0.16, fs=5.2)
        neural(cv, rge, pfe, ("RG-E", "PF", "exc"),
               color=dc.W_E, lw=1.2)
        neural(cv, rgf, pff, ("RG-F", "PF", "exc"),
               color=dc.W_F, lw=1.2)
        neural(cv, pfe, pine, ("PF", "PF-IN", "exc"),
               color=dc.W_E, lw=0.9)
        neural(cv, pine, pff, ("PF-IN", "PF", "inh"),
               color=dc.W_INH, lw=0.9, rad=-0.14)
        neural(cv, pff, pinf, ("PF", "PF-IN", "exc"),
               color=dc.W_F, lw=0.9)
        neural(cv, pinf, pfe, ("PF-IN", "PF", "inh"),
               color=dc.W_INH, lw=0.9, rad=0.14)

        mne = dc.neuron(cv, ext_x, 3.55, "MN-E\npools", fc="white",
                        r=0.34, fs=6.4)
        mnf = dc.neuron(cv, flx_x, 3.55, "MN-F\npools", fc="white",
                        r=0.34, fs=6.4)
        neural(cv, pfe, mne, ("PF", "MN", "exc"),
               color=dc.W_E, lw=1.25)
        neural(cv, pff, mnf, ("PF", "MN", "exc"),
               color=dc.W_F, lw=1.25)
        me = dc.muscle(cv, ext_x, 1.86, w=0.82, h=0.28)
        mf = dc.muscle(cv, flx_x, 1.86, w=0.82, h=0.28)
        cv.ax.text(ext_x, 1.54, "extensor\nmuscles", ha="center",
                   va="top", fontsize=5.6)
        cv.ax.text(flx_x, 1.54, "flexor\nmuscles", ha="center",
                   va="top", fontsize=5.6)
        plant(cv, mne, me, (f"MN-{side}-E", f"MUSCLE-{side}-E"),
              color="0.30", lw=1.2, label="activation", lfs=5.2,
              loff=(0.34, 0.0), lbox=False)
        plant(cv, mnf, mf, (f"MN-{side}-F", f"MUSCLE-{side}-F"),
              color="0.30", lw=1.2, label="activation", lfs=5.2,
              loff=(0.34, 0.0), lbox=False)
        sens = dc.input_box(cv, cx, 0.93, 1.62, 0.34, "Ia / II / Ib + foot contact",
                            "phase-dependent feedback")
        plant(cv, me, sens, (f"MUSCLE-{side}-E", f"SENSORY-{side}"),
              color=dc.W_E, lw=0.8, rad=-0.12)
        plant(cv, mf, sens, (f"MUSCLE-{side}-F", f"SENSORY-{side}"),
              color=dc.W_F, lw=0.8, rad=0.12)
        # Clean aggregate feedback arrows; detailed destinations are in B.
        dc.syn(cv, sens, pfe, True, color=dc.W_E, lw=0.85, dashed=True,
               rad=0.22)
        dc.syn(cv, sens, pff, True, color=dc.W_F, lw=0.85, dashed=True,
               rad=-0.22)
        dc.syn(cv, sens, rge, True, color=dc.W_E, lw=0.70, dashed=True,
               rad=0.34)
        dc.syn(cv, sens, rgf, True, color=dc.W_F, lw=0.70, dashed=True,
               rad=-0.34)
        sides[side] = dict(rge=rge, rgf=rgf, ine=ine, inf=inf)

    # Four DISTINCT directional commissural relays.  A single shared V3 or
    # C1 glyph is wrong: the left-to-right and right-to-left paths are
    # separate interneurons in both the compiled network and the literature
    # architecture.
    relays = {
        "LR": {
            "c1": dc.neuron(cv, 4.43, 7.18, "C1\nL->R",
                             fc=dc._tint(dc.OI_PURPLE), r=0.20, fs=4.6),
            "v3": dc.neuron(cv, 4.43, 8.20, "V3\nL->R",
                             fc=dc._tint(dc.OI_GREEN), r=0.20, fs=4.6),
        },
        "RL": {
            "c1": dc.neuron(cv, 6.17, 7.18, "C1\nR->L",
                             fc=dc._tint(dc.OI_PURPLE), r=0.20, fs=4.6),
            "v3": dc.neuron(cv, 6.17, 8.20, "V3\nR->L",
                             fc=dc._tint(dc.OI_GREEN), r=0.20, fs=4.6),
        },
    }
    for direction, src, dst in (
            ("LR", sides["L"], sides["R"]),
            ("RL", sides["R"], sides["L"])):
        c1 = relays[direction]["c1"]
        v3 = relays[direction]["v3"]
        neural(cv, src["rgf"], c1, ("RG-F", "CIN", "exc"),
               color=dc.W_F, lw=0.9, rad=0.10)
        neural(cv, c1, dst["rgf"], ("CIN", "RG-F", "inh"),
               color=dc.W_INH, lw=1.0, rad=-0.10)
        neural(cv, src["rge"], v3, ("RG-E", "CIN", "exc"),
               color=dc.W_E, lw=0.9, rad=-0.10)
        neural(cv, v3, dst["ine"], ("CIN", "RG-IN", "exc"),
               color=dc.W_EXC, lw=1.0, rad=0.10)

    for side in sides.values():
        neural(cv, drive, side["rge"], ("DRIVE", "RG-E", "exc"),
               color=dc.W_DESC, lw=0.9, rad=0.10)
        neural(cv, drive, side["rgf"], ("DRIVE", "RG-F", "exc"),
               color=dc.W_DESC, lw=0.9, rad=-0.10)
        neural(cv, posture, side["rge"], ("POSTURE", "RG-E", "exc"),
               color=dc.W_DESC, lw=0.8, dashed=True, rad=-0.10)

    # ------------------------------------------------------------------ B
    cv.ax.text(13.18, 9.45,
               "one complete antagonist pair; pathways repeated for all pools",
               ha="center", va="center", fontsize=6.0, color="0.35",
               style="italic")
    pfe = dc.input_box(cv, 11.55, 8.78, 1.05, 0.42, "PF-E", "E1/E2")
    pff = dc.input_box(cv, 14.72, 8.78, 1.05, 0.42, "PF-F", "F1/F2")
    mne = dc.neuron(cv, 11.55, 7.25, "MN\nknee-ext", r=0.35, fs=6.2)
    mnf = dc.neuron(cv, 14.72, 7.25, "MN\nknee-flx", r=0.35, fs=6.2)
    neural(cv, pfe, mne, ("PF", "MN", "exc"), color=dc.W_E, lw=1.25)
    neural(cv, pff, mnf, ("PF", "MN", "exc"), color=dc.W_F, lw=1.25)
    kinh = dc.neuron(cv, 13.05, 8.15, "KINH\nconditional", r=0.20,
                     fs=4.7)
    neural(cv, pff, kinh, ("PF", "KINH", "exc"), color=dc.W_F, lw=0.9)
    neural(cv, kinh, mne, ("KINH", "MN", "inh"),
           color=dc.W_INH, lw=1.1, rad=0.12)

    rce = dc.neuron(cv, 12.35, 7.22, "RC", fc="#dddddd", r=0.17, fs=5.4)
    rcf = dc.neuron(cv, 13.92, 7.22, "RC", fc="#dddddd", r=0.17, fs=5.4)
    neural(cv, mne, rce, ("MN", "RC", "exc"), color=dc.W_EXC, lw=0.9)
    neural(cv, rce, mne, ("RC", "MN", "inh"),
           color=dc.W_INH, lw=0.9, rad=0.30)
    neural(cv, mnf, rcf, ("MN", "RC", "exc"), color=dc.W_EXC, lw=0.9)
    neural(cv, rcf, mnf, ("RC", "MN", "inh"),
           color=dc.W_INH, lw=0.9, rad=-0.30)
    neural(cv, rce, rcf, ("RC", "RC", "inh"),
           color=dc.W_INH, lw=0.75, rad=-0.22)
    neural(cv, rcf, rce, ("RC", "RC", "inh"),
           color=dc.W_INH, lw=0.75, rad=-0.22)

    iae = dc.neuron(cv, 10.98, 4.72, "Ia", fc=dc._tint(dc.OI_ORANGE),
                    r=0.16, fs=5.2)
    iie = dc.neuron(cv, 11.55, 4.72, "II", fc=dc._tint(dc.OI_ORANGE),
                    r=0.16, fs=5.2)
    ibe = dc.neuron(cv, 12.12, 4.72, "Ib", fc=dc._tint(dc.OI_ORANGE),
                    r=0.16, fs=5.2)
    iaf = dc.neuron(cv, 14.15, 4.72, "Ia", fc=dc._tint(dc.OI_ORANGE),
                    r=0.16, fs=5.2)
    iif = dc.neuron(cv, 14.72, 4.72, "II", fc=dc._tint(dc.OI_ORANGE),
                    r=0.16, fs=5.2)
    ibf = dc.neuron(cv, 15.29, 4.72, "Ib", fc=dc._tint(dc.OI_ORANGE),
                    r=0.16, fs=5.2)
    for ia, ii, ib, mn in ((iae, iie, ibe, mne), (iaf, iif, ibf, mnf)):
        neural(cv, ia, mn, ("Ia", "MN", "exc"), color=dc.W_EXC,
               lw=0.9, rad=-0.10)
        neural(cv, ii, mn, ("II", "MN", "exc"), color=dc.W_EXC, lw=0.9)
        neural(cv, ib, mn, ("Ib", "MN", "inh"), color=dc.W_INH,
               lw=0.9, rad=0.10)

    iaine = dc.neuron(cv, 12.55, 5.72, "IaIN-E",
                      fc=dc._tint(dc.OI_PURPLE), r=0.20, fs=5.1)
    iainf = dc.neuron(cv, 13.72, 5.72, "IaIN-F",
                      fc=dc._tint(dc.OI_PURPLE), r=0.20, fs=5.1)
    neural(cv, iae, iaine, ("Ia", "IaIN", "exc"), color=dc.W_EXC,
           lw=0.9)
    neural(cv, iaf, iainf, ("Ia", "IaIN", "exc"), color=dc.W_EXC,
           lw=0.9)
    neural(cv, iaine, mnf, ("IaIN", "MN", "inh"),
           color=dc.W_INH, lw=1.1, rad=-0.12)
    neural(cv, iainf, mne, ("IaIN", "MN", "inh"),
           color=dc.W_INH, lw=1.1, rad=0.12)
    neural(cv, rce, iaine, ("RC", "IaIN", "inh"),
           color=dc.W_INH, lw=0.8, rad=-0.10)
    neural(cv, rcf, iainf, ("RC", "IaIN", "inh"),
           color=dc.W_INH, lw=0.8, rad=0.10)
    neural(cv, pff, iaine, ("PF", "IaIN", "exc"),
           color=dc.W_F, lw=0.75, dashed=True, rad=0.14)
    neural(cv, pff, iainf, ("PF", "IaIN", "exc"),
           color=dc.W_F, lw=0.75, dashed=True, rad=-0.14)

    ibexc = dc.neuron(cv, 10.90, 5.65, "IB-EXC", fc=dc._tint(dc.OI_GREEN),
                      r=0.20, fs=5.0)
    rge_gate = dc.input_box(cv, 10.93, 8.20, 0.70, 0.34, "RG-E",
                            "stance gate")
    neural(cv, ibe, ibexc, ("Ib", "IBEXC", "exc"), color=dc.W_EXC,
           lw=0.9)
    neural(cv, ibexc, mne, ("IBEXC", "MN", "exc"), color=dc.W_EXC,
           lw=0.9, rad=-0.14)
    neural(cv, rge_gate, ibexc, ("RG-E", "IBEXC", "exc"),
           color=dc.W_E, lw=0.8, dashed=True, rad=0.10)
    cv.ax.text(10.82, 6.18, "stance-extensor\nload sharing", ha="center",
               va="center", fontsize=5.0, color="0.34", style="italic")

    me = dc.muscle(cv, 11.55, 1.55, w=0.90, h=0.30)
    mf = dc.muscle(cv, 14.72, 1.55, w=0.90, h=0.30)
    cv.ax.text(11.55, 1.16, "knee extensor", ha="center", fontsize=5.8)
    cv.ax.text(14.72, 1.16, "knee flexor", ha="center", fontsize=5.8)
    plant(cv, mne, me, ("MN-E", "MUSCLE-E"), color="0.28", lw=1.15,
          label="a = clip(V/5 mV, 0, 1)", lfs=5.1, lpos=0.72,
          loff=(0.58, 0.0))
    plant(cv, mnf, mf, ("MN-F", "MUSCLE-F"), color="0.28", lw=1.15,
          label="a = clip(V/5 mV, 0, 1)", lfs=5.1, lpos=0.72,
          loff=(-0.58, 0.0))
    for side, mus, affs, col in (("E", me, (iae, iie, ibe), dc.W_E),
                                  ("F", mf, (iaf, iif, ibf), dc.W_F)):
        for kind, aff in zip(("Ia", "II", "Ib"), affs):
            plant(cv, mus, aff, (f"MUSCLE-{side}", f"{kind}-{side}"),
                  color=col, lw=0.78, rad=0.06 if kind == "Ia" else
                  (-0.06 if kind == "Ib" else 0.0))
    cv.ax.text(13.14, 3.45,
               "muscle velocity / length / force\nencoded as Ia / II / Ib currents",
               ha="center", va="center", fontsize=5.5, color="0.32",
               style="italic",
               bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none",
                         alpha=0.88))

    # ----------------------------------------------------------- contracts
    required_neural = {
        ("RG-E", "RG-IN", "exc"), ("RG-IN", "RG-F", "inh"),
        ("RG-F", "RG-IN", "exc"), ("RG-IN", "RG-E", "inh"),
        ("RG-E", "PF", "exc"), ("RG-F", "PF", "exc"),
        ("PF", "PF-IN", "exc"), ("PF-IN", "PF", "inh"),
        ("PF", "MN", "exc"), ("RG-F", "CIN", "exc"),
        ("CIN", "RG-F", "inh"), ("RG-E", "CIN", "exc"),
        ("CIN", "RG-IN", "exc"), ("MN", "RC", "exc"),
        ("RC", "MN", "inh"), ("RC", "RC", "inh"),
        ("RC", "IaIN", "inh"),
        ("Ia", "MN", "exc"), ("II", "MN", "exc"),
        ("Ib", "MN", "inh"), ("Ia", "IaIN", "exc"),
        ("IaIN", "MN", "inh"), ("Ib", "IBEXC", "exc"),
        ("IBEXC", "MN", "exc"), ("RG-E", "IBEXC", "exc"),
        ("PF", "KINH", "exc"),
        ("KINH", "MN", "inh"),
    }
    assert required_neural <= drawn, (
        f"reader schematic omitted required mechanism classes: "
        f"{sorted(required_neural - drawn)}")
    assert drawn <= compiled, (
        f"reader schematic contains noncompiled neural edges: "
        f"{sorted(drawn - compiled)}")

    required_plant = {
        (f"MN-{s}-{p}", f"MUSCLE-{s}-{p}")
        for s in ("L", "R") for p in ("E", "F")
    } | {
        (f"MUSCLE-{s}-{p}", f"SENSORY-{s}")
        for s in ("L", "R") for p in ("E", "F")
    } | {
        ("MN-E", "MUSCLE-E"), ("MN-F", "MUSCLE-F")
    } | {
        (f"MUSCLE-{s}", f"{a}-{s}")
        for s in ("E", "F") for a in ("Ia", "II", "Ib")
    }
    assert plant_drawn == required_plant, (
        f"reader schematic plant-interface mismatch: missing="
        f"{sorted(required_plant - plant_drawn)}, extra="
        f"{sorted(plant_drawn - required_plant)}")

    # Compact legend matching the glyphs actually used in this figure.
    y = 0.24
    cv.ax.add_patch(dc.Circle((0.75, y), 0.14, fc="white", ec="black",
                              lw=1.2))
    cv.ax.text(0.98, y, "neuron / aggregate pool", fontsize=6.2,
               va="center")
    cv.ax.add_patch(dc.FancyArrowPatch((2.62, y), (3.10, y),
                                       arrowstyle="-", lw=1.2,
                                       color=dc.W_EXC))
    cv.ax.add_patch(dc.Polygon([(3.34, y + 0.075), (3.34, y - 0.075),
                                (3.14, y)], closed=True, fc="white",
                               ec=dc.W_EXC, lw=1.1))
    cv.ax.text(3.45, y, "excitatory", fontsize=6.2, va="center")
    cv.ax.add_patch(dc.FancyArrowPatch((4.55, y), (5.05, y),
                                       arrowstyle="-", lw=1.2,
                                       color=dc.W_INH))
    cv.ax.add_patch(dc.Circle((5.15, y), 0.07, fc="black", ec="black"))
    cv.ax.text(5.30, y, "inhibitory", fontsize=6.2, va="center")
    cv.ax.add_patch(dc.Ellipse((6.62, y), 0.52, 0.20,
                               fc=dc._tint(dc.OI_SKY), ec="black", lw=1.0))
    cv.ax.text(7.00, y, "muscle", fontsize=6.2, va="center")
    cv.ax.plot([7.88, 8.38], [y, y], color="0.35", lw=1.1)
    cv.ax.text(8.48, y, "motor activation", fontsize=6.2, va="center")
    cv.ax.plot([9.78, 10.28], [y, y], color=dc.W_E, lw=0.9,
               ls=(0, (3, 2)))
    cv.ax.text(10.38, y, "aggregate sensory feedback", fontsize=6.2,
               va="center")
    cv.ax.text(15.62, 0.22,
               "all solid neural pathways checked against the compiled SNS network",
               ha="right", va="center", fontsize=5.4, color="0.38",
               style="italic")
    return cv.save("circuit_literature", fmts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fmt", default="pdf,svg,png")
    args = ap.parse_args()
    fmts = [x.strip() for x in args.fmt.split(",") if x.strip()]
    dc.OUT.mkdir(exist_ok=True)
    outs = make_figure(fmts)
    print(f"reader schematic: {outs}")


if __name__ == "__main__":
    main()
