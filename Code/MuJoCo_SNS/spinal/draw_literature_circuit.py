"""Reader-facing spinal-circuit schematic grounded in the literature.

Panel A follows the bilateral layered organization used by Rybak,
Shevtsova, Shinohara, and Deng: supraspinal drive -> RG -> PF -> MN ->
muscle, with c1/V3 commissural relays centered at the midline.
Panel B follows the mechanism-at-a-time presentation used by Deng and
Di Russo: one complete antagonist knee motif with no dangling interface
paths.  The dense ``circuit_dengstyle`` figure remains the audit view.

2026-09-17 bilateral repair (Ben's annotated review):
  * the RIGHT E/F columns are mirrored so both RG-F half-centers face
    the midline (extensor columns on the outside, flexor inside);
  * every box edge leaves/arrives at an explicit NAMED PORT on the
    rectangle boundary (RG input top-center, MN output bottom-center,
    PF-IN excitation bottom-inner, cross-inhibition + aggregate
    feedback on the inner side edge) -- no glyph-to-glyph shrink
    guessing, which is what made PF connections look unattached;
  * aggregate sensory feedback is drawn as one explicit bus per side
    (sensory box -> vertical bus -> PF junction -> RG junction -> four
    short labeled branches) instead of four long sweeping arcs;
  * V3 -> contralateral InE runs through a waypoint below the In row so
    it cannot graze the InF / RG-F glyphs after the mirror.

Contracts: every displayed neural edge CLASS is asserted against a
freshly compiled representative SNS network; every per-side PATH is
asserted in an instance-level contract (side, source glyph, target
glyph, sign); and every wire's recorded geometric start/end point is
asserted to lie on its registered glyph boundary (boxes: exact port;
circles: within a ring) with nonzero length and no self-edge.
Plant-interface paths have the separate explicit completeness contract.
"""
from __future__ import annotations

import argparse

import numpy as np
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


def make_figure(fmts, tag=""):
    net, tg = _compiled_representative()
    compiled = set(dc._edge_groups(net))
    drawn, plant_drawn, agg_drawn = set(), set(), set()
    inst = set()
    geom, refs = [], []

    cv = dc.Canvas(16.2, 10.4)
    _panel_frame(cv, 0.35, 10.25, 0.48, 10.0,
                 "A   Bilateral spinal architecture")
    _panel_frame(cv, 10.48, 15.87, 0.48, 10.0,
                 "B   Representative right-knee sensorimotor motif")

    # ------------------------------------------------ glyph + wire registry
    GLYPH = {}

    def neuron(name, x, y, *a, **k):
        g = dc.neuron(cv, x, y, *a, **k)
        GLYPH[name] = ("circ", g)
        return g

    def box(name, x, y, w, h, *a, **k):
        g = dc.input_box(cv, x, y, w, h, *a, **k)
        GLYPH[name] = ("box", (x, y, w, h))
        return g

    def muscle(name, x, y, **k):
        g = dc.muscle(cv, x, y, **k)
        GLYPH[name] = ("circ", (x, y, k.get("h", 0.26) / 2.0))
        return g

    def bport(name, fx, fy):
        """Named port: (fx, fy) in half-extent fractions, |fx|==1 or
        |fy|==1 so the point lies exactly on the rectangle boundary."""
        cx, cy, w, h = GLYPH[name][1]
        onx, ony = abs(abs(fx) - 1) < 1e-9, abs(abs(fy) - 1) < 1e-9
        assert (onx and abs(fy) <= 1 + 1e-9) or \
               (ony and abs(fx) <= 1 + 1e-9), (name, fx, fy)
        return (cx + fx * w / 2, cy + fy * h / 2, 0.0)

    def wp(name, x, y):
        GLYPH[name] = ("pt", (x, y))
        cv.ax.add_patch(dc.Circle((x, y), 0.045, fc="0.55", ec="none",
                                  zorder=2.5))
        return (x, y, 0.0)

    def pt_ref(name, x, y):
        GLYPH[name] = ("pt", (x, y))
        return (x, y, 0.0)

    def neural(p, q, key, src=None, dst=None, **kw):
        assert key in compiled, f"displayed neural edge absent from net: {key}"
        dc.syn(cv, p, q, key[2] == "exc", record=geom, **kw)
        drawn.add(key)
        if src is not None:
            inst.add((src, dst, key[2]))
        refs.append((geom[-1], src, dst, key))

    def plant(p, q, key, src=None, dst=None, **kw):
        dc.syn(cv, p, q, True, record=geom, **kw)
        plant_drawn.add(key)
        refs.append((geom[-1], src, dst, key))

    def agg(p, q, src, dst, color, lw, rad):
        """One aggregate sensory-feedback branch (separate contract)."""
        dc.syn(cv, p, q, True, color=color, lw=lw, rad=rad, dashed=True,
               record=geom)
        agg_drawn.add((src, dst))
        refs.append((geom[-1], src, dst, ("AGG", src, dst, "exc")))

    def wire_line(x0, y0, x1, y1, src, dst, **kw):
        cv.ax.plot([x0, x1], [y0, y1], zorder=1.6, **kw)
        refs.append((dict(start=np.array([x0, y0]), end=np.array([x1, y1]),
                          exc=True, rad=0.0, dashed=True), src, dst,
                     ("AGG", src, dst, "line")))

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

    drive = box("DRIVE", 4.28, 9.17, 1.55, 0.42, "DRIVE", "MLR-like command")
    posture = box("POSTURE", 6.35, 9.17, 1.55, 0.42, "POSTURE / BALANCE",
                  "tonic + feedback")

    sides = {}
    for side, cx in (("L", 2.75), ("R", 7.85)):
        # Mirror the whole right functional column: flexor half-centers
        # face the midline on BOTH sides (defect 1 of the 2026-09-17
        # annotated review).
        med = 1.0 if side == "L" else -1.0
        ext_x, flx_x = cx - med * 0.78, cx + med * 0.78
        rge = neuron(f"RG-E_{side}", ext_x, 7.83, "RG-E", "NaP",
                     fc=dc._tint(dc.W_E), r=0.34, fs=7.0)
        rgf = neuron(f"RG-F_{side}", flx_x, 7.83, "RG-F", "NaP",
                     fc=dc._tint(dc.W_F), r=0.34, fs=7.0)
        ine = neuron(f"InE_{side}", cx - med * 0.30, 6.98, "InE",
                     fc=dc._tint(dc.OI_PURPLE), r=0.19, fs=5.7)
        inf = neuron(f"InF_{side}", cx + med * 0.30, 6.98, "InF",
                     fc=dc._tint(dc.OI_PURPLE), r=0.19, fs=5.7)
        neural(rge, ine, ("RG-E", "RG-IN", "exc"),
               src=f"RG-E_{side}", dst=f"InE_{side}", color=dc.W_E, lw=1.1)
        neural(ine, rgf, ("RG-IN", "RG-F", "inh"),
               src=f"InE_{side}", dst=f"RG-F_{side}",
               color=dc.W_INH, lw=1.1, rad=-0.05 * med)
        neural(rgf, inf, ("RG-F", "RG-IN", "exc"),
               src=f"RG-F_{side}", dst=f"InF_{side}", color=dc.W_F, lw=1.1)
        neural(inf, rge, ("RG-IN", "RG-E", "inh"),
               src=f"InF_{side}", dst=f"RG-E_{side}",
               color=dc.W_INH, lw=1.1, rad=-0.05 * med)
        pfe = box(f"PF-E_{side}", ext_x, 5.60, 1.12, 0.52,
                  "PF-E1 | PF-E2", "stance phases")
        pff = box(f"PF-F_{side}", flx_x, 5.60, 1.12, 0.52,
                  "PF-F1 | PF-F2", "swing phases")
        pine = neuron(f"IN-E_{side}", cx - med * 0.27, 4.90, "IN-E",
                      fc=dc._tint(dc.OI_PURPLE), r=0.16, fs=5.2)
        pinf = neuron(f"IN-F_{side}", cx + med * 0.27, 4.90, "IN-F",
                      fc=dc._tint(dc.OI_PURPLE), r=0.16, fs=5.2)
        # RG input on the box TOP-CENTER port; MN output leaves the
        # BOTTOM-CENTER port (defect 2: no more floating stems).
        neural(rge, bport(f"PF-E_{side}", 0, 1), ("RG-E", "PF", "exc"),
               src=f"RG-E_{side}", dst=f"PF-E_{side}",
               color=dc.W_E, lw=1.2)
        neural(rgf, bport(f"PF-F_{side}", 0, 1), ("RG-F", "PF", "exc"),
               src=f"RG-F_{side}", dst=f"PF-F_{side}",
               color=dc.W_F, lw=1.2)
        # PF -> PF-IN excitation leaves the BOTTOM-INNER port; the
        # triangle lands on the interneuron perimeter (defect 3).
        neural(bport(f"PF-E_{side}", med * 0.75, -1), pine,
               ("PF", "PF-IN", "exc"), src=f"PF-E_{side}",
               dst=f"IN-E_{side}", color=dc.W_E, lw=0.9)
        neural(bport(f"PF-F_{side}", -med * 0.75, -1), pinf,
               ("PF", "PF-IN", "exc"), src=f"PF-F_{side}",
               dst=f"IN-F_{side}", color=dc.W_F, lw=0.9)
        # cross-inhibition arrives on the OPPOSITE box's inner side edge,
        # routed as a visibly separate curve from the excitation stems.
        neural(pine, bport(f"PF-F_{side}", -med, 0.45),
               ("PF-IN", "PF", "inh"), src=f"IN-E_{side}",
               dst=f"PF-F_{side}", color=dc.W_INH, lw=0.9, rad=-0.14 * med)
        neural(pinf, bport(f"PF-E_{side}", med, 0.45),
               ("PF-IN", "PF", "inh"), src=f"IN-F_{side}",
               dst=f"PF-E_{side}", color=dc.W_INH, lw=0.9, rad=0.14 * med)

        mne = neuron(f"MN-E_{side}", ext_x, 3.55, "MN-E\npools", fc="white",
                     r=0.34, fs=6.4)
        mnf = neuron(f"MN-F_{side}", flx_x, 3.55, "MN-F\npools", fc="white",
                     r=0.34, fs=6.4)
        neural(bport(f"PF-E_{side}", 0, -1), mne, ("PF", "MN", "exc"),
               src=f"PF-E_{side}", dst=f"MN-E_{side}",
               color=dc.W_E, lw=1.25)
        neural(bport(f"PF-F_{side}", 0, -1), mnf, ("PF", "MN", "exc"),
               src=f"PF-F_{side}", dst=f"MN-F_{side}",
               color=dc.W_F, lw=1.25)
        me = muscle(f"MUS-E_{side}", ext_x, 1.86, w=0.82, h=0.28)
        mf = muscle(f"MUS-F_{side}", flx_x, 1.86, w=0.82, h=0.28)
        cv.ax.text(ext_x, 1.54, "extensor\nmuscles", ha="center",
                   va="top", fontsize=5.6)
        cv.ax.text(flx_x, 1.54, "flexor\nmuscles", ha="center",
                   va="top", fontsize=5.6)
        plant(mne, me, (f"MN-{side}-E", f"MUSCLE-{side}-E"),
              src=f"MN-E_{side}", dst=f"MUS-E_{side}",
              color="0.30", lw=1.2, label="activation", lfs=5.2,
              loff=(0.34 * med, 0.0), lbox=False)
        plant(mnf, mf, (f"MN-{side}-F", f"MUSCLE-{side}-F"),
              src=f"MN-F_{side}", dst=f"MUS-F_{side}",
              color="0.30", lw=1.2, label="activation", lfs=5.2,
              loff=(0.34 * med, 0.0), lbox=False)
        sens = box(f"SENS_{side}", cx, 0.93, 1.62, 0.34,
                   "Ia / II / Ib + foot contact", "phase-dependent feedback")
        # encoder paths terminate exactly on the sensory-box top edge
        # (defect 4); feedback leaves one bus instead of sweeping arcs.
        plant(me, bport(f"SENS_{side}", -med * 0.62, 1),
              (f"MUSCLE-{side}-E", f"SENSORY-{side}"),
              src=f"MUS-E_{side}", dst=f"SENS_{side}",
              color=dc.W_E, lw=0.8, rad=-0.12 * med)
        plant(mf, bport(f"SENS_{side}", med * 0.62, 1),
              (f"MUSCLE-{side}-F", f"SENSORY-{side}"),
              src=f"MUS-F_{side}", dst=f"SENS_{side}",
              color=dc.W_F, lw=0.8, rad=0.12 * med)
        # explicit feedback bus: sensory box top-center -> vertical bus
        # through the limb midline -> PF junction -> RG junction, with
        # four short branches to named ports (defect 4).
        jpf = pt_ref(f"JPF_{side}", cx, 5.20)
        jrg = pt_ref(f"JRG_{side}", cx, 7.40)
        btop = pt_ref(f"BUSS_{side}", cx, 1.10)
        busend = pt_ref(f"BUSE_{side}", cx, 7.40)
        wire_line(btop[0], btop[1], busend[0], busend[1],
                  f"BUSS_{side}", f"BUSE_{side}",
                  color="0.45", lw=0.9, ls=(0, (3, 2)))
        for jy in (jpf[1], jrg[1]):
            cv.ax.add_patch(dc.Circle((cx, jy), 0.05, fc="0.45",
                                      ec="none", zorder=2.5))
        cv.ax.text(cx + 0.12, 4.28, "aggregate\nsensory\nfeedback",
                   fontsize=4.6, ha="left", va="center", color="0.40",
                   style="italic", zorder=6)
        agg(jpf, bport(f"PF-E_{side}", med, -0.45),
            f"JPF_{side}", f"PF-E_{side}", dc.W_E, 0.85, rad=-0.20 * med)
        agg(jpf, bport(f"PF-F_{side}", -med, -0.45),
            f"JPF_{side}", f"PF-F_{side}", dc.W_F, 0.85, rad=0.20 * med)
        # RG branches arrive at the circle BOTTOM so the excitatory
        # triangle cannot collide with the InE/InF inhibitory dots on
        # the lower-left / lower-right sectors.
        agg(jrg, (ext_x, 7.83 - 0.34, 0.0), f"JRG_{side}", f"RG-E_{side}",
            dc.W_E, 0.70, rad=-0.25 * med)
        agg(jrg, (flx_x, 7.83 - 0.34, 0.0), f"JRG_{side}", f"RG-F_{side}",
            dc.W_F, 0.70, rad=0.25 * med)
        sides[side] = dict(rge=rge, rgf=rgf, ine=ine, inf=inf, med=med)

    # Four DISTINCT directional commissural relays.  A single shared V3 or
    # C1 glyph is wrong: the left-to-right and right-to-left paths are
    # separate interneurons in both the compiled network and the literature
    # architecture.  V3 -> contra-InE runs through a waypoint below the In
    # row (the mirrored InF glyph sits on the direct chord).
    wp_e = {"LR": wp("WP-IN_R", 8.15, 6.38), "RL": wp("WP-IN_L", 2.45, 6.38)}
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
    GLYPH["C1_LR"] = ("circ", relays["LR"]["c1"])
    GLYPH["V3_LR"] = ("circ", relays["LR"]["v3"])
    GLYPH["C1_RL"] = ("circ", relays["RL"]["c1"])
    GLYPH["V3_RL"] = ("circ", relays["RL"]["v3"])
    for direction, src, dst in (
            ("LR", sides["L"], sides["R"]),
            ("RL", sides["R"], sides["L"])):
        c1, v3 = relays[direction]["c1"], relays[direction]["v3"]
        md = dst["med"]
        neural(src["rgf"], c1, ("RG-F", "CIN", "exc"),
               src=f"RG-F_{direction[0]}", dst=f"C1_{direction}",
               color=dc.W_F, lw=0.9, rad=0.10)
        neural(c1, dst["rgf"], ("CIN", "RG-F", "inh"),
               src=f"C1_{direction}", dst=f"RG-F_{direction[1]}",
               color=dc.W_INH, lw=1.0, rad=-0.05 * md)
        neural(src["rge"], v3, ("RG-E", "CIN", "exc"),
               src=f"RG-E_{direction[0]}", dst=f"V3_{direction}",
               color=dc.W_E, lw=0.9, rad=-0.22 * src["med"])
        neural(v3, wp_e[direction], ("CIN", "RG-IN", "exc"),
               src=f"V3_{direction}", dst=f"WP-IN_{direction[1]}",
               color=dc.W_EXC, lw=1.0, rad=-0.10 * md)
        neural(wp_e[direction], dst["ine"], ("CIN", "RG-IN", "exc"),
               src=f"WP-IN_{direction[1]}", dst=f"InE_{direction[1]}",
               color=dc.W_EXC, lw=1.0)

    for side in sides.values():
        med = side["med"]
        neural(bport("DRIVE", 0, -1), side["rge"], ("DRIVE", "RG-E", "exc"),
               src="DRIVE", dst=f"RG-E_{'L' if med > 0 else 'R'}",
               color=dc.W_DESC, lw=0.9, rad=0.12 * med)
        neural(bport("DRIVE", 0, -1), side["rgf"], ("DRIVE", "RG-F", "exc"),
               src="DRIVE", dst=f"RG-F_{'L' if med > 0 else 'R'}",
               color=dc.W_DESC, lw=0.9, rad=-0.12 * med)
        neural(bport("POSTURE", 0, -1), side["rge"],
               ("POSTURE", "RG-E", "exc"), src="POSTURE",
               dst=f"RG-E_{'L' if med > 0 else 'R'}",
               color=dc.W_DESC, lw=0.8, dashed=True, rad=-0.10 * med)

    # ------------------------------------------------------------------ B
    cv.ax.text(13.18, 9.45,
               "one complete antagonist pair; pathways repeated for all pools",
               ha="center", va="center", fontsize=6.0, color="0.35",
               style="italic")
    pfe = box("PF-E_B", 11.55, 8.78, 1.05, 0.42, "PF-E", "E1/E2")
    pff = box("PF-F_B", 14.72, 8.78, 1.05, 0.42, "PF-F", "F1/F2")
    mne = neuron("MN-E_B", 11.55, 7.25, "MN\nknee-ext", r=0.35, fs=6.2)
    mnf = neuron("MN-F_B", 14.72, 7.25, "MN\nknee-flx", r=0.35, fs=6.2)
    neural(bport("PF-E_B", 0, -1), mne, ("PF", "MN", "exc"),
           src="PF-E_B", dst="MN-E_B", color=dc.W_E, lw=1.25)
    neural(bport("PF-F_B", 0, -1), mnf, ("PF", "MN", "exc"),
           src="PF-F_B", dst="MN-F_B", color=dc.W_F, lw=1.25)
    kinh = neuron("KINH_B", 13.05, 8.15, "KINH\nconditional", r=0.20,
                  fs=4.7)
    neural(bport("PF-F_B", -0.65, -1), kinh, ("PF", "KINH", "exc"),
           src="PF-F_B", dst="KINH_B", color=dc.W_F, lw=0.9)
    neural(kinh, mne, ("KINH", "MN", "inh"), src="KINH_B", dst="MN-E_B",
           color=dc.W_INH, lw=1.1, rad=0.12)

    rce = neuron("RC_E", 12.35, 7.22, "RC", fc="#dddddd", r=0.17, fs=5.4)
    rcf = neuron("RC_F", 13.92, 7.22, "RC", fc="#dddddd", r=0.17, fs=5.4)
    neural(mne, rce, ("MN", "RC", "exc"), src="MN-E_B", dst="RC_E",
           color=dc.W_EXC, lw=0.9)
    neural(rce, mne, ("RC", "MN", "inh"), src="RC_E", dst="MN-E_B",
           color=dc.W_INH, lw=0.9, rad=0.30)
    neural(mnf, rcf, ("MN", "RC", "exc"), src="MN-F_B", dst="RC_F",
           color=dc.W_EXC, lw=0.9)
    neural(rcf, mnf, ("RC", "MN", "inh"), src="RC_F", dst="MN-F_B",
           color=dc.W_INH, lw=0.9, rad=-0.30)
    neural(rce, rcf, ("RC", "RC", "inh"), src="RC_E", dst="RC_F",
           color=dc.W_INH, lw=0.75, rad=-0.22)
    neural(rcf, rce, ("RC", "RC", "inh"), src="RC_F", dst="RC_E",
           color=dc.W_INH, lw=0.75, rad=-0.22)

    iae = neuron("Ia_E", 10.98, 4.72, "Ia", fc=dc._tint(dc.OI_ORANGE),
                 r=0.16, fs=5.2)
    iie = neuron("II_E", 11.55, 4.72, "II", fc=dc._tint(dc.OI_ORANGE),
                 r=0.16, fs=5.2)
    ibe = neuron("Ib_E", 12.12, 4.72, "Ib", fc=dc._tint(dc.OI_ORANGE),
                 r=0.16, fs=5.2)
    iaf = neuron("Ia_F", 14.15, 4.72, "Ia", fc=dc._tint(dc.OI_ORANGE),
                 r=0.16, fs=5.2)
    iif = neuron("II_F", 14.72, 4.72, "II", fc=dc._tint(dc.OI_ORANGE),
                 r=0.16, fs=5.2)
    ibf = neuron("Ib_F", 15.29, 4.72, "Ib", fc=dc._tint(dc.OI_ORANGE),
                 r=0.16, fs=5.2)
    for ia, ii, ib, mn, s in ((iae, iie, ibe, mne, "E"),
                              (iaf, iif, ibf, mnf, "F")):
        neural(ia, mn, ("Ia", "MN", "exc"), src=f"Ia_{s}", dst=f"MN-{s}_B",
               color=dc.W_EXC, lw=0.9, rad=-0.10)
        neural(ii, mn, ("II", "MN", "exc"), src=f"II_{s}", dst=f"MN-{s}_B",
               color=dc.W_EXC, lw=0.9)
        neural(ib, mn, ("Ib", "MN", "inh"), src=f"Ib_{s}", dst=f"MN-{s}_B",
               color=dc.W_INH, lw=0.9, rad=0.10)

    iaine = neuron("IaIN-E_B", 12.55, 5.72, "IaIN-E",
                   fc=dc._tint(dc.OI_PURPLE), r=0.20, fs=5.1)
    iainf = neuron("IaIN-F_B", 13.72, 5.72, "IaIN-F",
                   fc=dc._tint(dc.OI_PURPLE), r=0.20, fs=5.1)
    neural(iae, iaine, ("Ia", "IaIN", "exc"), src="Ia_E", dst="IaIN-E_B",
           color=dc.W_EXC, lw=0.9)
    neural(iaf, iainf, ("Ia", "IaIN", "exc"), src="Ia_F", dst="IaIN-F_B",
           color=dc.W_EXC, lw=0.9)
    neural(iaine, mnf, ("IaIN", "MN", "inh"), src="IaIN-E_B", dst="MN-F_B",
           color=dc.W_INH, lw=1.1, rad=-0.12)
    neural(iainf, mne, ("IaIN", "MN", "inh"), src="IaIN-F_B", dst="MN-E_B",
           color=dc.W_INH, lw=1.1, rad=0.12)
    neural(rce, iaine, ("RC", "IaIN", "inh"), src="RC_E", dst="IaIN-E_B",
           color=dc.W_INH, lw=0.8, rad=-0.10)
    neural(rcf, iainf, ("RC", "IaIN", "inh"), src="RC_F", dst="IaIN-F_B",
           color=dc.W_INH, lw=0.8, rad=0.10)
    neural(bport("PF-F_B", -0.85, -1), iaine, ("PF", "IaIN", "exc"),
           src="PF-F_B", dst="IaIN-E_B",
           color=dc.W_F, lw=0.75, dashed=True, rad=0.14)
    neural(bport("PF-F_B", -0.35, -1), iainf, ("PF", "IaIN", "exc"),
           src="PF-F_B", dst="IaIN-F_B",
           color=dc.W_F, lw=0.75, dashed=True, rad=-0.14)

    ibexc = neuron("IBEXC_B", 10.90, 5.65, "IB-EXC",
                   fc=dc._tint(dc.OI_GREEN), r=0.20, fs=5.0)
    rge_gate = box("RGGATE_B", 10.93, 8.20, 0.70, 0.34, "RG-E",
                   "stance gate")
    neural(ibe, ibexc, ("Ib", "IBEXC", "exc"), src="Ib_E", dst="IBEXC_B",
           color=dc.W_EXC, lw=0.9)
    neural(ibexc, mne, ("IBEXC", "MN", "exc"), src="IBEXC_B", dst="MN-E_B",
           color=dc.W_EXC, lw=0.9, rad=-0.14)
    neural(bport("RGGATE_B", 0, -1), ibexc, ("RG-E", "IBEXC", "exc"),
           src="RGGATE_B", dst="IBEXC_B",
           color=dc.W_E, lw=0.8, dashed=True, rad=0.10)
    cv.ax.text(10.82, 6.18, "stance-extensor\nload sharing", ha="center",
               va="center", fontsize=5.0, color="0.34", style="italic")

    me = muscle("MUS-E_B", 11.55, 1.55, w=0.90, h=0.30)
    mf = muscle("MUS-F_B", 14.72, 1.55, w=0.90, h=0.30)
    cv.ax.text(11.55, 1.16, "knee extensor", ha="center", fontsize=5.8)
    cv.ax.text(14.72, 1.16, "knee flexor", ha="center", fontsize=5.8)
    plant(mne, me, ("MN-E", "MUSCLE-E"), src="MN-E_B", dst="MUS-E_B",
          color="0.28", lw=1.15,
          label="a = clip(V/5 mV, 0, 1)", lfs=5.1, lpos=0.72,
          loff=(0.58, 0.0))
    plant(mnf, mf, ("MN-F", "MUSCLE-F"), src="MN-F_B", dst="MUS-F_B",
          color="0.28", lw=1.15,
          label="a = clip(V/5 mV, 0, 1)", lfs=5.1, lpos=0.72,
          loff=(-0.58, 0.0))
    for side, mus, affs, col in (("E", me, (iae, iie, ibe), dc.W_E),
                                 ("F", mf, (iaf, iif, ibf), dc.W_F)):
        for kind, aff in zip(("Ia", "II", "Ib"), affs):
            plant(mus, aff, (f"MUSCLE-{side}", f"{kind}-{side}"),
                  src=f"MUS-{side}_B", dst=f"{kind}_{side}",
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

    # Instance-level contract: exact source glyph + exact target glyph +
    # sign, for every per-side pathway and both commissural directions
    # (report step 4).
    required_inst = set()
    for sd in ("L", "R"):
        required_inst |= {
            (f"RG-E_{sd}", f"InE_{sd}", "exc"),
            (f"InE_{sd}", f"RG-F_{sd}", "inh"),
            (f"RG-F_{sd}", f"InF_{sd}", "exc"),
            (f"InF_{sd}", f"RG-E_{sd}", "inh"),
            (f"RG-E_{sd}", f"PF-E_{sd}", "exc"),
            (f"RG-F_{sd}", f"PF-F_{sd}", "exc"),
            (f"PF-E_{sd}", f"IN-E_{sd}", "exc"),
            (f"IN-E_{sd}", f"PF-F_{sd}", "inh"),
            (f"PF-F_{sd}", f"IN-F_{sd}", "exc"),
            (f"IN-F_{sd}", f"PF-E_{sd}", "inh"),
            (f"PF-E_{sd}", f"MN-E_{sd}", "exc"),
            (f"PF-F_{sd}", f"MN-F_{sd}", "exc"),
        }
    required_inst |= {
        ("RG-F_L", "C1_LR", "exc"), ("C1_LR", "RG-F_R", "inh"),
        ("RG-E_L", "V3_LR", "exc"), ("V3_LR", "WP-IN_R", "exc"),
        ("WP-IN_R", "InE_R", "exc"),
        ("RG-F_R", "C1_RL", "exc"), ("C1_RL", "RG-F_L", "inh"),
        ("RG-E_R", "V3_RL", "exc"), ("V3_RL", "WP-IN_L", "exc"),
        ("WP-IN_L", "InE_L", "exc"),
    }
    missing_inst = required_inst - inst
    assert not missing_inst, \
        f"instance contract missing paths: {sorted(missing_inst)}"
    for e in inst:
        assert e[0] != e[1], f"self-edge in instance contract: {e}"

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

    required_agg = {
        (f"JPF_{s}", f"PF-E_{s}") for s in ("L", "R")
    } | {(f"JPF_{s}", f"PF-F_{s}") for s in ("L", "R")} | {
        (f"JRG_{s}", f"RG-E_{s}") for s in ("L", "R")
    } | {(f"JRG_{s}", f"RG-F_{s}") for s in ("L", "R")}
    assert agg_drawn == required_agg, (
        f"aggregate feedback contract mismatch: missing="
        f"{sorted(required_agg - agg_drawn)}, extra="
        f"{sorted(agg_drawn - required_agg)}")

    # Geometry contract (report step 5): every wire's recorded start/end
    # lies on its registered glyph boundary; no degenerate paths.
    for rec, sg, dg, key in refs:
        s, e = rec["start"], rec["end"]
        assert float(np.linalg.norm(e - s)) > 0.12, \
            f"degenerate path {key} {sg}->{dg}"
        for pterm, gname, is_end in ((s, sg, False), (e, dg, True)):
            kind, par = GLYPH[gname]
            if kind == "circ":
                x, y, r = par
                dd = float(np.linalg.norm(pterm - (x, y)))
                lo, hi = (r - 0.03, r + 0.34) if is_end \
                    else (r - 0.03, r + 0.12)
                assert lo <= dd <= hi, \
                    f"{key} endpoint not on {g}: d={dd:.3f}"
            elif kind == "box":
                bx, by, bw, bh = par
                dx = abs(float(pterm[0]) - bx) - bw / 2
                dy = abs(float(pterm[1]) - by) - bh / 2
                if dx > 0 or dy > 0:
                    dist = float(np.hypot(max(dx, 0.0), max(dy, 0.0)))
                else:
                    dist = float(min(-dx, -dy))   # depth inside the box
                # ends may sit up to one synapse-marker setback off the
                # boundary (exc base 0.02, inhibitory dot 0.10)
                lim = 0.16 if is_end else 0.04
                assert dist <= lim, \
                    f"{key} endpoint not on box {par}: ({dx:.3f},{dy:.3f})"
            else:
                assert float(np.linalg.norm(
                    pterm - np.asarray(par))) < 0.04, \
                    f"{key} endpoint not on waypoint {par}"

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
               "solid edges checked against the compiled SNS network;\n"
               "every wire boundary-checked against its glyph",
               ha="right", va="center", fontsize=5.4, color="0.38",
               style="italic")
    name = "circuit_literature" + (f"_{tag}" if tag else "")
    return cv.save(name, fmts)


def inst_missing_guard(required, have):
    """Return a set-like that passes only when required <= have (keeps the
    failure message readable)."""
    missing = required - have
    assert not missing, f"instance contract missing paths: {sorted(missing)}"
    return required


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fmt", default="pdf,svg,png")
    ap.add_argument("--tag", default="",
                    help="render to circuit_literature_<tag>.* (review "
                         "renders; the promoted name stays untouched)")
    args = ap.parse_args()
    fmts = [x.strip() for x in args.fmt.split(",") if x.strip()]
    dc.OUT.mkdir(exist_ok=True)
    outs = make_figure(fmts, tag=args.tag)
    print(f"reader schematic: {outs}")


if __name__ == "__main__":
    main()
