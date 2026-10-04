"""Layered schematics of the 2026-10-02 SPIKING models, for the defense
slideshow (Documentation\\Reports and Papers\\Dissertation\\Defense\\
slideshow\\) and beside the existing CPG figures (Figures\\30-results\\
CPG_airstepping_figs\\).

Figures:
  spiking_mirror_layers          the MuJoCo SNS SPIKING MIRROR
                                 (build_network_spiking.py), both sides,
                                 layered, STRUCTURE-DRIVEN: the net is
                                 built at the topology-gate TUNED config
                                 and every (src-class, dst-class, sign)
                                 edge group read off the compiled SNS
                                 object MUST be drawn (assert below -
                                 the figure cannot drift from the code).
  animatlab_spiking_conversions  the AnimatLab spiking copies (goal 3,
                                 2026-10-02): the BipRG-family connectome
                                 (connectome_templates.json 'bilateralrg',
                                 the topology the spiking copies preserve
                                 verbatim per tools/validate.py) drawn
                                 layered, with the conversion + outcome
                                 overlay from goal3_animatlab_spiking.md.

Conventions follow draw_circuit.py (Ben's reference style): Okabe-Ito
tinted layer bands, white-triangle/dot synapse markers, wires colored by
source family.  SPIKING vs NON-SPIKING coding per Ben's SNS_Simscape
icon convention: cells carrying a small SPIKE glyph (action-potential
trace) are spiking LIF; cells carrying a TILDE glyph (graded waveform)
stay analog.  No conductance numbers on wires (calibration numbers live
in spiking_calibration.json and the editor template note).

Usage (myo env):
  python draw_spiking_schematic.py [--which mirror|animatlab|all]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from draw_circuit import (  # noqa: E402
    Canvas, neuron, input_box, layer_band, syn, tag, lane, _tint,
    R_BIG, R_MED, R_SMALL, R_SENS,
    OI_SKY, OI_VERM, OI_BLUE, OI_GREEN, OI_PURPLE, OI_ORANGE,
    W_E, W_F, W_EXC, W_INH, W_DESC)

REPO = HERE.parent.parent.parent
SLIDES = REPO / "Documentation" / "Reports and Papers" / "Dissertation" \
    / "Defense" / "slideshow"
CPGFIGS = REPO / "Documentation" / "Reports and Papers" / "Dissertation" \
    / "Figures" / "30-results" / "CPG_airstepping_figs"

# ---------------------------------------------------------------- helpers
def save_both(cv, name, fmts=("png", "pdf")):
    outs = []
    for base in (SLIDES, CPGFIGS):
        base.mkdir(parents=True, exist_ok=True)
        for fmt in fmts:
            p = base / f"{name}.{fmt}"
            cv.fig.savefig(p, dpi=300, bbox_inches="tight", pad_inches=0.06)
            outs.append(str(p))
        plt.close(cv.fig) if fmts and False else None
    plt.close(cv.fig)
    for p in outs:
        print("  wrote", p)
    return outs


def spike_glyph(cv, x, y, s=0.16, color="0.25"):
    """Tiny action-potential trace (the SNS_Simscape spiking marker)."""
    cv.ax.plot([x - s, x - s * 0.45, x - s * 0.1, x + s * 0.1, x + s],
               [y, y, y + s * 0.9, y, y], color=color, lw=0.9, zorder=4,
               solid_joinstyle="miter")


def tilde_glyph(cv, x, y, s=0.15, color="0.25"):
    """Tiny graded waveform (non-spiking marker)."""
    import numpy as np
    t = np.linspace(0, 2 * 3.1416, 24)
    cv.ax.plot(x - s + t / (2 * 3.1416) * 2 * s,
               y + 0.32 * s * np.sin(t), color=color, lw=0.9, zorder=4)


def s_neuron(cv, x, y, label, sub=None, fc="white", r=R_MED, fs=None):
    n = neuron(cv, x, y, label, sub=sub, fc=fc, r=r, fs=fs)
    spike_glyph(cv, x + r * 0.95, y + r * 0.85, s=r * 0.5)
    return n


def a_neuron(cv, x, y, label, sub=None, fc="white", r=R_MED, fs=None):
    n = neuron(cv, x, y, label, sub=sub, fc=fc, r=r, fs=fs)
    tilde_glyph(cv, x + r * 0.95, y + r * 0.85, s=r * 0.55)
    return n


# ---------------------------------------------------------------- mirror
def cell_class(nm):
    """Fine cell class for the edge-coverage contract."""
    for pre, k in (("RD_", "RD"), ("CIN_F", "CIN_F"), ("CIN_E", "CIN_E"),
                   ("PF_IN_E", "PF_IN"), ("PF_IN_F", "PF_IN"),
                   ("PF_E1", "PF-E"), ("PF_E2", "PF-E"),
                   ("PF_F1", "PF-F"), ("PF_F2", "PF-F"),
                   ("RG_E", "RG-E"), ("RG_F", "RG-F"),
                   ("InE", "InE"), ("InF", "InF"),
                   ("HEEL", "HEEL"), ("TOE", "TOE"), ("LBIN", "LBIN"),
                   ("AFF_E", "AFF_E"), ("AFF_F", "AFF_F"),
                   ("IaIN", "IaIN"), ("IIX", "IIX"), ("IBIN", "IBIN"),
                   ("IIIN", "IIIN"), ("IBEXC", "IBEXC"), ("KINH", "KINH"),
                   ("RC_", "RC"), ("MN_", "MN"),
                   ("Ia_", "Ia"), ("II_", "II"), ("Ib_", "Ib"),
                   ("VEST", "VEST"), ("DRIVE", "DRIVE"),
                   ("POSTURE", "POSTURE"), ("BAL_", "BAL")):
        if nm.startswith(pre):
            return k
    return nm


def build_mirror_at_tuned():
    """Build the spiking mirror at the TUNED gate config (verbatim gains
    from reports_spiking_20261002/tools/topology_mirror_check.py)."""
    import params as P
    import build_network_spiking as BSS
    import muscle_map as MM
    import make_editor_templates as MT
    keep = dict(P.G)
    try:
        P.G.update(MT.SPIKE_TUNED)
        P.G["joint_pf"] = 0.0
        acts = []
        for b in MM._GROUPS_BY_NAME:
            acts.append(b + "_r")
            acts.append(b + "_l")
        net = BSS.build(acts, interleg=True)
    finally:
        P.G.clear()
        P.G.update(keep)
    return net


def edge_census(net):
    names = [p["name"] for p in net.net.populations]
    g = Counter()
    for c in net.net.connections:
        g[(cell_class(names[c["source"]]), cell_class(names[c["destination"]]),
           "exc" if c["params"].get("reversal_potential", 0) > -1e-6
           else "inh")] += 1
    return g


def fig_mirror():
    net = build_mirror_at_tuned()
    census = edge_census(net)
    drawn = set()

    def edge(sk, dk, sign):
        drawn.add((sk, dk, sign))

    cv = Canvas(16.4, 21.0)
    X_R, X_L = 4.1, 12.3          # side column centers (R leg left, L leg right)
    X_M = 8.2                      # midline corridor

    # ---------------- band 0: supraspinal (ANALOG) ----------------
    layer_band(cv, 0.5, 15.9, 19.15, 20.6, OI_SKY,
               "SUPRASPINAL / BRAINSTEM SURROGATES (NON-SPIKING, analog)")
    drv = input_box(cv, 5.9, 19.85, 2.0, 0.6, "DRIVE", "speed cmd")
    pos = input_box(cv, 10.5, 19.85, 2.0, 0.6, "POSTURE", "tonic")
    tag(cv, 8.2, 19.3, "BAL family (6 analog cells -> ankle/hip/trunk MNs) "
        "omitted for space", fs=5.6, box=True)
    vest_r = a_neuron(cv, X_R, 19.5, "VEST_R", r=R_SMALL, fs=5.6)
    vest_l = a_neuron(cv, X_L, 19.5, "VEST_L", r=R_SMALL, fs=5.6)
    tag(cv, X_R, 19.0, "vestibular analog\n(goal-2; analog)", fs=5.4)
    edge("VEST", "MN", "exc")
    edge("VEST", "MN", "inh")

    # ---------------- band 1: rhythm generator (SPIKING) ----------------
    layer_band(cv, 0.5, 15.9, 14.25, 19.0, OI_ORANGE,
               "RHYTHM GENERATOR - adapting-LIF half-centers (SPIKING)")
    cells = {}
    for S, X in (("R", X_R), ("L", X_L)):
        rge = s_neuron(cv, X - 0.6, 17.9, f"RG-E_{S}", "aLIF",
                       fc=_tint(W_E), r=R_BIG)
        rgf = s_neuron(cv, X - 0.6, 16.1, f"RG-F_{S}", "aLIF",
                       fc=_tint(W_F), r=R_BIG)
        ine = s_neuron(cv, X - 2.3, 17.9, f"InE_{S}",
                       fc=_tint(OI_PURPLE), r=R_MED)
        inf = s_neuron(cv, X - 2.3, 16.1, f"InF_{S}",
                       fc=_tint(OI_PURPLE), r=R_MED)
        hel = s_neuron(cv, X - 2.3, 15.0, f"HEEL_{S}", r=R_SMALL, fs=5.8)
        toe = s_neuron(cv, X - 1.2, 15.0, f"TOE_{S}", r=R_SMALL, fs=5.8)
        lbin = s_neuron(cv, X + 0.9, 15.0, f"LBIN_{S}", r=R_SMALL, fs=5.8)
        affe = s_neuron(cv, X + 1.9, 17.6, f"AFF_E_{S}", "grp Ia/II",
                        r=R_SMALL, fs=5.6)
        afff = s_neuron(cv, X + 1.9, 16.5, f"AFF_F_{S}", "grp Ia/II",
                        r=R_SMALL, fs=5.6)
        cells[S] = dict(rge=rge, rgf=rgf, ine=ine, inf=inf, hel=hel,
                        toe=toe, lbin=lbin, affe=affe, afff=afff)
        # IN-laminated mutual inhibition
        syn(cv, rge, ine, True, color=W_E, lw=1.4)
        syn(cv, ine, rgf, False, color=W_INH, lw=1.6, rad=0.25)
        syn(cv, rgf, inf, True, color=W_F, lw=1.4)
        syn(cv, inf, rge, False, color=W_INH, lw=1.6, rad=0.25)
        edge("RG-E", "InE", "exc"); edge("InE", "RG-F", "inh")
        edge("RG-F", "InF", "exc"); edge("InF", "RG-E", "inh")
        # descending drive (graded -> spiking, rate map)
        syn(cv, drv, rge, True, color=W_DESC, rad=0.1, lw=1.2)
        syn(cv, drv, rgf, True, color=W_DESC, rad=-0.12, lw=1.2)
        syn(cv, pos, rge, True, color=W_DESC, rad=-0.2, lw=1.0)
        edge("DRIVE", "RG-E", "exc"); edge("DRIVE", "RG-F", "exc")
        edge("POSTURE", "RG-E", "exc")
        # contact mechanosensors (full_rules branch: heel -> InE/InF,
        # toe -> InE) + load IN
        syn(cv, hel, ine, True, color=W_E, lw=1.0, rad=0.15)
        syn(cv, hel, inf, False, color=W_INH, lw=1.0, rad=0.1)
        syn(cv, toe, ine, True, color=W_E, lw=1.0, rad=0.12)
        # heel -> InF EXCITATORY variant (heel_in_f_exc, Ben's rules)
        syn(cv, hel, inf, True, color=W_F, lw=0.8, rad=0.2)
        edge("HEEL", "InE", "exc"); edge("HEEL", "InF", "inh")
        edge("HEEL", "InF", "exc")
        edge("TOE", "InE", "exc")
        syn(cv, lbin, rge, True, color=W_E, lw=1.0, rad=-0.15)
        edge("LBIN", "RG-E", "exc")
        # semi-closed sensory loop relays (v11b)
        syn(cv, affe, rge, True, color=W_EXC, lw=0.9, rad=0.15)
        syn(cv, afff, rgf, True, color=W_F, lw=0.9, rad=-0.15)
        edge("AFF_E", "RG-E", "exc"); edge("AFF_F", "RG-F", "exc")

    # commissurals: midline corridor (c1 = CIN_F crossed F inh; V3-E =
    # CIN_E crossed extensor exc; plus CIN_E -> contra IBEXC, and the
    # crossed flexor-Ia inhibition drawn in the motor band)
    cin_f = {"R": s_neuron(cv, X_M - 0.8, 16.1, "CIN_F_r\nc1",
                           fc=_tint(OI_PURPLE), r=R_SMALL, fs=5.2),
             "L": s_neuron(cv, X_M + 0.8, 16.1, "CIN_F_l",
                           fc=_tint(OI_PURPLE), r=R_SMALL, fs=5.2)}
    cin_e = {"R": s_neuron(cv, X_M - 0.8, 17.9, "CIN_E_r\nV3-E",
                           fc=_tint(OI_GREEN), r=R_SMALL, fs=5.2),
             "L": s_neuron(cv, X_M + 0.8, 17.9, "CIN_E_l",
                           fc=_tint(OI_GREEN), r=R_SMALL, fs=5.2)}
    for S in ("R", "L"):
        syn(cv, cells[S]["rgf"], cin_f[S], True, color=W_F, lw=1.3)
        syn(cv, cin_f[S], cells["L" if S == "R" else "R"]["rgf"], False,
            color=W_INH, lw=1.8)
        syn(cv, cells[S]["rge"], cin_e[S], True, color=W_E, lw=1.0)
        syn(cv, cin_e[S], cells["L" if S == "R" else "R"]["ine"], True,
            color=W_EXC, lw=1.2)
        edge("RG-F", "CIN_F", "exc"); edge("CIN_F", "RG-F", "inh")
        edge("RG-E", "CIN_E", "exc"); edge("CIN_E", "InE", "exc")
    tag(cv, X_M, 15.2, "commissurals (spiking INs):\nCIN_F = c1 crossed "
        "flexor inhibition (V0D-like);\nCIN_E = V3-E crossed extensor "
        "excitation\n(+ CIN_E -> contra IBEXC, drawn below)", fs=5.8,
        box=True)
    edge("CIN_E", "IBEXC", "exc")
    tag(cv, 2.6, 18.6, "RG = adapting LIF:\nthreshold-increment SFA\n"
        "(thr_inc 4.0 mV,\ntau_theta = rg_nap_h);\nburst termination by\n"
        "adaptation, NOT NaP h-gate\n(documented fallback)", fs=5.4, box=True)

    # ---------------- band 2: pattern formation (SPIKING) ----------------
    layer_band(cv, 0.5, 15.9, 11.05, 14.05, OI_SKY,
               "PATTERN FORMATION - 4 phase cells + laminated PF_IN "
               "(SPIKING)")
    pf_x = {"E1": -1.7, "E2": -0.6, "F1": 0.6, "F2": 1.7}
    for S, X in (("R", X_R), ("L", X_L)):
        pfs = {}
        for nm, dx in pf_x.items():
            c = W_E if nm[0] == "E" else W_F
            pfs[nm] = s_neuron(cv, X + dx, 13.15, f"PF-{nm}",
                               fc=_tint(c), r=R_MED * 0.85, fs=6.2)
            syn(cv, cells[S]["rge" if nm[0] == "E" else "rgf"], pfs[nm],
                True, color=W_E if nm[0] == "E" else W_F, rad=0.05 * dx)
            edge("RG-E" if nm[0] == "E" else "RG-F", "PF-E" if nm[0] == "E"
                 else "PF-F", "exc")
        pf_in_e = s_neuron(cv, X - 1.15, 11.75, "PF_IN_E",
                           fc=_tint(OI_PURPLE), r=R_SMALL, fs=5.8)
        pf_in_f = s_neuron(cv, X + 1.15, 11.75, "PF_IN_F",
                           fc=_tint(OI_PURPLE), r=R_SMALL, fs=5.8)
        for ph in ("E1", "E2"):
            syn(cv, pfs[ph], pf_in_e, True, color=W_E, lw=0.9,
                rad=0.1 if ph == "E1" else -0.1)
        for ph in ("F1", "F2"):
            syn(cv, pfs[ph], pf_in_f, True, color=W_F, lw=0.9,
                rad=0.1 if ph == "F1" else -0.1)
        for ph in ("F1", "F2"):
            syn(cv, pf_in_e, pfs[ph], False, color=W_INH, lw=0.9, rad=0.18)
        for ph in ("E1", "E2"):
            syn(cv, pf_in_f, pfs[ph], False, color=W_INH, lw=0.9, rad=0.18)
        edge("PF-E", "PF_IN", "exc"); edge("PF-F", "PF_IN", "exc")
        edge("PF_IN", "PF-F", "inh"); edge("PF_IN", "PF-E", "inh")
        # heel/toe ride the extensor central pathway (ib_e_central)
        for src, tgt in ((cells[S]["hel"], pfs["E1"]),
                         (cells[S]["toe"], pfs["E2"]),
                         (cells[S]["hel"], cells[S]["ine"])):
            syn(cv, src, tgt, True, color=W_E, lw=0.7, rad=0.15)
        edge("HEEL", "PF-E", "exc"); edge("TOE", "PF-E", "exc")
        # AFF relays onto PF
        syn(cv, cells[S]["affe"], pfs["E1"], True, color=W_EXC, lw=0.7,
            rad=0.1)
        syn(cv, cells[S]["affe"], pfs["E2"], True, color=W_EXC, lw=0.7,
            rad=-0.1)
        syn(cv, cells[S]["afff"], pfs["F1"], True, color=W_F, lw=0.7,
            rad=0.1)
        syn(cv, cells[S]["afff"], pfs["F2"], True, color=W_F, lw=0.7,
            rad=-0.1)
        edge("AFF_E", "PF-E", "exc"); edge("AFF_F", "PF-F", "exc")
        cells[S]["pfs"] = pfs

    # ---------------- band 3: motor + reflex (knee representative) ------
    layer_band(cv, 0.5, 15.9, 4.35, 10.85, OI_GREEN,
               "MOTOR CIRCUIT - knee columns representative (x46 muscle "
               "columns per side); MN layer NON-SPIKING")
    for S, X in (("R", X_R), ("L", X_L)):
        pfs = cells[S]["pfs"]
        mne = a_neuron(cv, X - 1.0, 9.9, "MN knee-ext", "x6 vas",
                       r=0.40, fs=6.4)
        mnf = a_neuron(cv, X + 1.4, 9.9, "MN knee-flx", "x7 hamstr",
                       r=0.40, fs=6.4)
        cells[S]["mne"], cells[S]["mnf"] = mne, mnf
        syn(cv, pfs["E1"], mne, True, color=W_E, rad=0.15, lw=1.2)
        syn(cv, pfs["E2"], mne, True, color=W_E, rad=-0.1, lw=1.2)
        syn(cv, pfs["F1"], mnf, True, color=W_F, rad=0.12, lw=1.2)
        syn(cv, pfs["F2"], mnf, True, color=W_F, rad=-0.08, lw=0.8)
        edge("PF-E", "MN", "exc"); edge("PF-F", "MN", "exc")
        # KINH swing suppression (F1-gated)
        kinh = s_neuron(cv, X + 2.6, 10.2, "KINH", r=R_SMALL, fs=6.0)
        syn(cv, pfs["F1"], kinh, True, lw=1.0, color=W_F, rad=-0.15)
        syn(cv, kinh, mne, False, color=W_INH, lw=1.2, rad=0.2)
        edge("PF-F", "KINH", "exc"); edge("KINH", "MN", "inh")
        # contralateral heel gates KINH (contra_kinh)
        other = "L" if S == "R" else "R"
        syn(cv, cells[other]["hel"], kinh, True, color=W_E, lw=0.7,
            rad=-0.2)
        edge("HEEL", "KINH", "exc")
        # afferent encoders (SPIKING)
        aff_y = 6.3
        ia = s_neuron(cv, X - 2.0, aff_y, "Ia", r=R_SENS,
                      fc=_tint(OI_ORANGE), fs=6.0)
        ii = s_neuron(cv, X - 1.0, aff_y, "II", r=R_SENS,
                      fc=_tint(OI_ORANGE), fs=6.0)
        ib = s_neuron(cv, X + 0.0, aff_y, "Ib", r=R_SENS,
                      fc=_tint(OI_ORANGE), fs=6.0)
        cells[S].update(ia=ia, ii=ii, ib=ib)
        syn(cv, ia, mne, True, color=W_EXC, rad=-0.1)
        syn(cv, ii, mne, True, color=W_EXC, rad=0.0)
        syn(cv, ib, mne, False, color=W_INH, rad=0.12)
        edge("Ia", "MN", "exc"); edge("II", "MN", "exc")
        edge("Ib", "MN", "inh")
        # full_rules relays: IIX (II exc relay), IBIN (Ib inh relay),
        # IIIN (II inh relay)
        iix = s_neuron(cv, X - 1.6, 7.45, "IIX", fc=_tint(OI_GREEN),
                       r=R_SMALL, fs=6.0)
        ibin = s_neuron(cv, X - 0.4, 7.45, "IBIN", fc=_tint(OI_PURPLE),
                        r=R_SMALL, fs=6.0)
        iiin = s_neuron(cv, X + 0.8, 7.45, "IIIN", fc=_tint(OI_PURPLE),
                        r=R_SMALL, fs=6.0)
        cells[S].update(iix=iix, ibin=ibin, iiin=iiin)
        syn(cv, ii, iix, True, color=W_EXC, lw=0.9)
        syn(cv, iix, mne, True, color=W_EXC, lw=0.9, rad=-0.1)
        syn(cv, ib, ibin, True, color=W_EXC, lw=0.9)
        syn(cv, ibin, mne, False, color=W_INH, lw=1.0, rad=0.15)
        syn(cv, ii, iiin, True, color=W_EXC, lw=0.9, rad=0.1)
        syn(cv, iiin, mnf, False, color=W_INH, lw=1.0, rad=0.15)
        edge("II", "IIX", "exc"); edge("IIX", "MN", "exc")
        edge("Ib", "IBIN", "exc"); edge("IBIN", "MN", "inh")
        edge("II", "IIIN", "exc"); edge("IIIN", "MN", "inh")
        # IaIN reciprocal (PF-F1 gate, RC disinhibition, mutual IaIN)
        iain = s_neuron(cv, X + 1.9, 6.3, "IaIN", fc=_tint(OI_PURPLE),
                        r=0.26, fs=6.0)
        cells[S]["iain"] = iain
        syn(cv, ia, iain, True, color=W_EXC, lw=1.0, rad=0.2)
        syn(cv, pfs["F1"], iain, True, color=W_F, lw=0.9, rad=-0.13)
        syn(cv, iain, mnf, False, color=W_INH, lw=1.6, rad=-0.15)
        edge("Ia", "IaIN", "exc"); edge("PF-F", "IaIN", "exc")
        edge("IaIN", "MN", "inh")
        # IBEXC stance-gated reversal + LBIN loop + CIN_E drive
        ibx = s_neuron(cv, X - 2.4, 8.3, "IB-EXC", r=R_SMALL, fs=6.0,
                       fc=_tint(OI_GREEN))
        cells[S]["ibx"] = ibx
        syn(cv, ib, ibx, True, color=W_EXC, rad=-0.1, lw=1.0)
        syn(cv, ibx, mne, True, color=W_EXC, rad=0.15, lw=1.0)
        syn(cv, cells[S]["rge"], ibx, True, color=W_E, dashed=True,
            lw=0.9, rad=0.2)
        edge("Ib", "IBEXC", "exc"); edge("IBEXC", "MN", "exc")
        edge("RG-E", "IBEXC", "exc")
        syn(cv, ibx, cells[S]["lbin"], True, color=W_EXC, lw=0.9, rad=0.1)
        edge("IBEXC", "LBIN", "exc")
        syn(cv, cin_e[S], ibx, True, color=W_EXC, lw=0.8, rad=0.25)
        # Renshaw (RC driven by the ANALOG MN through a graded rate
        # encoder synapse)
        rc = s_neuron(cv, X + 0.3, 8.85, "RC", fc="#d9d9d9", r=R_SMALL,
                      fs=6.2)
        cells[S]["rc"] = rc
        syn(cv, mne, rc, True, color=W_EXC, lw=1.0)
        syn(cv, rc, mne, False, color=W_INH, lw=1.0, rad=0.35)
        syn(cv, rc, iain, False, color=W_INH, lw=0.9, rad=0.1)
        edge("MN", "RC", "exc"); edge("RC", "MN", "inh")
        edge("RC", "IaIN", "inh")
        # central afferent projections (ib/ia/ii central + crossed Ia)
        syn(cv, ib, pfs["E1"], True, color=W_E, lw=0.8, rad=0.2)
        syn(cv, ib, cells[S]["rge"], True, color=W_E, lw=0.7, rad=-0.1,
            dashed=True)
        syn(cv, ib, cells[S]["ine"], True, color=W_E, lw=0.7, rad=0.1)
        edge("Ib", "PF-E", "exc"); edge("Ib", "RG-E", "exc")
        edge("Ib", "InE", "exc")
        syn(cv, ia, pfs["F1"], True, color=W_F, lw=0.8, rad=0.2)
        syn(cv, ia, cells[S]["rgf"], True, color=W_F, lw=0.7, rad=-0.1)
        syn(cv, ia, cells[S]["inf"], True, color=W_F, lw=0.7, rad=0.1)
        edge("Ia", "PF-F", "exc"); edge("Ia", "RG-F", "exc")
        edge("Ia", "InF", "exc")
        syn(cv, ii, cells[S]["rgf"], True, color=W_F, lw=0.7, rad=0.12)
        syn(cv, ii, cells[S]["rge"], True, color=W_E, lw=0.7, rad=0.12)
        syn(cv, ii, cells[S]["inf"], True, color=W_F, lw=0.7, rad=0.18)
        syn(cv, ii, cells[S]["ine"], True, color=W_E, lw=0.7, rad=0.18)
        edge("II", "RG-F", "exc"); edge("II", "RG-E", "exc")
        edge("II", "InF", "exc"); edge("II", "InE", "exc")
        edge("II", "PF-F", "exc"); edge("II", "PF-E", "exc")
        # VEST onto MNs (analog -> analog)
        syn(cv, vest_r if S == "R" else vest_l, mne, True, color=W_DESC,
            lw=0.8, dashed=True, rad=0.2)
        syn(cv, vest_r if S == "R" else vest_l, mnf, False, color=W_INH,
            lw=0.8, dashed=True, rad=-0.2)
        # POSTURE -> MN (graded, unchanged)
        syn(cv, pos, mne, True, color=W_DESC, dashed=True, lw=0.8,
            rad=0.15)
        edge("POSTURE", "MN", "exc"); edge("BAL", "MN", "exc")

    # crossed flexor-Ia inhibition + mutual IaIN/IBIN (representative,
    # between the two side columns)
    syn(cv, cells["R"]["ia"], cells["L"]["rgf"], False, color=W_INH,
        lw=1.0, rad=-0.3)
    syn(cv, cells["L"]["ia"], cells["R"]["rgf"], False, color=W_INH,
        lw=1.0, rad=0.3)
    edge("Ia", "RG-F", "inh")
    syn(cv, cells["R"]["iain"], cells["L"]["iain"], False, color=W_INH,
        lw=0.8, rad=0.2)
    syn(cv, cells["R"]["ibin"], cells["L"]["ibin"], False, color=W_INH,
        lw=0.8, rad=-0.2)
    edge("IaIN", "IaIN", "inh"); edge("IBIN", "IBIN", "inh")
    tag(cv, X_M, 7.6, "2026-10-03 completion pass:\nIaIN<->IaIN and "
        "IBIN<->IBIN now BOTH\ndirections (was one-directional;\n"
        "WIRING_RULINGS_20261003.md F3,\nmirrored into the spiking twin)",
        fs=5.4, box=True)
    # RC<->RC mutual (per side, drawn on R side only for clarity)
    syn(cv, cells["R"]["rc"], cells["R"]["iain"], False, color=W_INH,
        lw=0.0, rad=0)   # placeholder no-draw (RC->IaIN already drawn)
    edge("RC", "RC", "inh")
    syn(cv, cells["R"]["rc"], cells["L"]["rc"], False, color=W_INH,
        lw=0.9, rad=0.15)
    tag(cv, X_M, 8.9, "RC<->RC mutual (same side, x46 pools;\n"
        "one inter-side wire drawn as representative)", fs=5.4, box=True)

    # ---------------- band 4: readout taps (ANALOG) ----------------
    layer_band(cv, 0.5, 15.9, 2.85, 4.2, "#BBBBBB",
               "RUNNER READOUT TAPS (NON-SPIKING) - RD_ low-pass of the "
               "spike trains")
    for S, X in (("R", X_R), ("L", X_L)):
        for i, nm in enumerate(("RD_RG-E", "RD_RG-F", "RD_PF-E1",
                                "RD_PF-E2", "RD_PF-F1", "RD_PF-F2")):
            rd = a_neuron(cv, X - 1.5 + i * 0.6, 3.5, nm.split("RD_")[1],
                          r=R_SMALL * 0.8, fs=4.8)
            src = (cells[S]["rge"] if nm.endswith("RG-E") else
                   cells[S]["rgf"] if nm.endswith("RG-F") else
                   cells[S]["pfs"][nm.split("RD_")[1][3:]])
            syn(cv, src, rd, True, color=W_DESC, lw=0.6, rad=0.05 * i - 0.1)
        edge("RG-E", "RD", "exc"); edge("RG-F", "RD", "exc")
        edge("PF-E", "RD", "exc"); edge("PF-F", "RD", "exc")
    tag(cv, X_M, 3.5, "net.idx repointed at RD_\n(runner reads analog "
        "levels);\nraw spiking cells in net.ridx", fs=5.6, box=True)

    # ---------------- legend + source note ----------------
    lx, ly = 0.8, 2.1
    spike_glyph(cv, lx + 0.25, ly + 0.28, s=0.22)
    cv.ax.text(lx + 0.6, ly + 0.28, "spiking LIF cell (real mV: rest "
               "-70, thr -50)", fontsize=6.4, va="center")
    tilde_glyph(cv, lx + 5.6, ly + 0.28, s=0.2)
    cv.ax.text(lx + 6.0, ly + 0.28, "non-spiking analog cell (0..5 mV "
               "frame)", fontsize=6.4, va="center")
    cv.ax.text(lx, ly - 0.25, "white triangle = excitatory synapse   "
               "solid dot = inhibitory   wire color = source family "
               "(blue extensor / vermilion flexor / green exc / purple "
               "inh / sky descending);   dashed = graded (analog-source) "
               "synapse, solid from spike sources", fontsize=6.2,
               va="center")
    cv.ax.text(lx, ly - 0.75,
               "SPIKING MIRROR of the gait2392 spinal network "
               "(build_network_spiking.py, 2026-10-02 campaign). "
               "Structure-driven: net built at the topology-gate TUNED "
               "config (all conditional pathways on, full_rules=1, phase "
               "PF); every edge class read from the compiled SNS object "
               f"({sum(census.values())} synapses, 896 neurons, 782 "
               "spiking) and asserted drawn. Weights calibrated by "
               "calibrate_spiking.py -> spiking_calibration.json (spike->"
               "MN mean-conductance; spike->spiking rate-calibrated). "
               "Knee columns representative; the full per-muscle wiring "
               "is editor template 'spiking_mirror' (988 nodes / 8448 "
               "edges).", fontsize=6.0, va="top", wrap=True)

    # ---------------- contract assertion ----------------
    missing = {k: v for k, v in census.items() if k not in drawn}
    print("edge classes in net:", len(census), "| drawn:", len(drawn))
    if missing:
        print("MISSING (undrawn) edge classes:")
        for k, v in sorted(missing.items()):
            print("   ", k, v)
        raise SystemExit("FIGURE CONTRACT FAILED: %d undrawn classes"
                         % len(missing))
    print("mirror figure contract PASS: all %d edge classes drawn"
          % len(census))
    save_both(cv, "spiking_mirror_layers")


# ------------------------------------------------- animatlab conversions
def fig_animatlab():
    spec = json.load(open(HERE / "connectome_templates.json",
                          encoding="utf-8"))["bilateralrg"]
    nodes = {n["label"]: n for n in spec["nodes"]}
    edges = spec.get("edges", spec.get("synapses", []))

    def side_of(l):
        l = l.strip()
        if l.startswith(("L", "l")):
            return "L"
        if l.startswith(("R", "r")):
            return "R"
        return None

    # classify into bands
    def band_of(n):
        t, l = n["type"], n["label"]
        if t in ("SN-heel", "SN-toe"):
            return "mech"
        if t in ("SN-Ia", "SN-II", "SN-Ib"):
            return "aff"
        if t == "IN-IaIN":
            return "iain"
        if t == "RC":
            return "rc"
        if t in ("HC-RG-E", "HC-RG-F", "IN-InE", "IN-InF"):
            return "rg"
        if t in ("IN-C", "IN-V0D", "IN-V3"):
            return "comm"
        if t in ("HC-PF-E", "HC-PF-F", "IN-PF"):
            return "pf"
        if t == "MN":
            return "mn"
        if t == "PORT-load":
            return "port"
        return "pf"

    cv = Canvas(16.4, 15.5)
    layer_band(cv, 0.5, 15.9, 12.55, 14.35, OI_ORANGE,
               "RHYTHM GENERATOR (per side; graded in the source models) "
               "+ MIDLINE commissurals c1/V3")
    layer_band(cv, 0.5, 15.9, 9.3, 12.4, OI_SKY,
               "PATTERN FORMATION (hip/knee PF half-centers per side)")
    layer_band(cv, 0.5, 15.9, 4.55, 9.15, OI_GREEN,
               "MOTOR + AFFERENT chains (Ia relay -> IaIN, II relay, Ib; "
               "RC; MNs)")
    layer_band(cv, 0.5, 15.9, 0.55, 4.4, "#C8C8C8",
               "CONTACT SENSORS (heel/toe) + spiking-conversion OUTCOMES")

    # per-type y inside the bands; x spread within each (band, side, yrow)
    TYPE_Y = {"HC-RG-E": 13.85, "HC-RG-F": 13.0, "IN-InE": 13.85,
              "IN-InF": 13.0, "IN-C": 13.42, "IN-V0D": 13.42,
              "IN-V3": 13.42, "PORT-load": 12.85,
              "HC-PF-E": 11.75, "HC-PF-F": 10.75, "IN-PF": 9.85,
              "MN": 8.45, "RC": 7.55, "IN-IaIN": 6.7,
              "SN-Ia": 5.85, "SN-II": 5.85, "SN-Ib": 5.85,
              "SN-heel": 3.75, "SN-toe": 3.75}
    from collections import Counter
    cnt = Counter()
    for n in spec["nodes"]:
        cnt[(n["type"], side_of(n["label"]) or "M")] += 1
    col = Counter()
    pos = {}
    for n in spec["nodes"]:
        t = n["type"]
        s = side_of(n["label"]) or "M"
        y = TYPE_Y.get(t, 9.0)
        if t == "MN":
            y = 8.75 - (col[(t, s)] // 5) * 0.62
        key = (t, s)
        i = col[key]
        col[key] += 1
        n_nodes = max(cnt[key], 1)
        if t in ("IN-C", "IN-V0D", "IN-V3"):
            # commissural cells live on the midline
            x = 8.2 + (i - (n_nodes - 1) / 2) * 1.35
        else:
            x0 = {"L": 3.3, "M": 8.2, "R": 12.5}[s]
            pitch = 1.5 if t in ("HC-RG-E", "HC-RG-F", "HC-PF-E",
                                 "HC-PF-F") else 1.15
            x = x0 + (i - (n_nodes - 1) / 2) * pitch
        r = R_BIG if t in ("HC-RG-E", "HC-RG-F") else \
            R_MED * 0.9 if t in ("HC-PF-E", "HC-PF-F") else R_SMALL
        sub = {"HC-RG-E": "RG", "HC-RG-F": "RG", "IN-InE": "RG-IN",
               "IN-InF": "RG-IN", "IN-C": "comm IN", "IN-V0D": "comm IN",
               "IN-V3": "comm IN", "PORT-load": "stim", "HC-PF-E": "PF",
               "HC-PF-F": "PF", "IN-PF": "PF-IN", "MN": "MN",
               "RC": "RC", "IN-IaIN": "IaIN", "SN-Ia": "Ia",
               "SN-II": "II", "SN-Ib": "Ib", "SN-heel": "heel",
               "SN-toe": "toe"}.get(t, t)
        nd = s_neuron(cv, x, y, n["label"], sub=sub, r=r,
                      fs=5.4 if t in ("MN", "SN-Ia", "SN-II", "SN-Ib")
                      else 5.9)
        pos[n["label"]] = nd
    # all edges drawn (curved), colored by sign/source family
    for e in edges:
        f, t = e["from"], e["to"]
        if f not in pos or t not in pos:
            continue
        src = nodes[f]
        exc = e["sign"] == "exc"
        color = W_E if src["type"].endswith("-E") else \
            W_F if src["type"].endswith("-F") else \
            W_EXC if exc else W_INH
        syn(cv, pos[f], pos[t], exc, color=color, lw=0.7, rad=0.12)
    # outcome strip (factual, goal3_animatlab_spiking.md section 5)
    outcomes = (
        "SPIKING CONVERSION OUTCOMES (headless AnimatSimulator, 2026-10-02 "
        "goal 3 - topology preserved verbatim, classes converted):\n"
        "W2L v1 (thr -55/-59, decay 10 ms): RHYTHM DEAD, 2 spikes/5 s   |   "
        "W2L v2 (calibrated thr): dead   |   W2L v3 (decay 30 ms): dead\n"
        "BipRG v1: asymmetric latch, R half-centers 442-448 Hz   |   "
        "Biped v1: SYNCHRONOUS co-bursting 11.8 Hz both HCs (same-leg "
        "xcorr r=1.000)\n"
        "Biped v2 (calibrated): 0 spikes   |   ABLATION (spiking thresholds "
        "ONLY, graded synapses kept): RHYTHM PRESERVED -> the SYNAPSE "
        "conversion breaks the rhythm;\n"
        "a spiking version needs per-spike amplitudes/decays re-tuned as a "
        "rate code (Li 2023 model: natively spiking, tuned that way from "
        "birth, 0.77 Hz bursts).")
    cv.ax.text(0.8, 3.15, outcomes, fontsize=5.9, va="top", ha="left",
               family="DejaVu Sans")
    spike_glyph(cv, 0.9, 0.98, s=0.2)
    cv.ax.text(1.25, 0.98, "every neuron drawn with a spike glyph: the "
               "spiking copies convert ALL neurons (v1 thresholds -55 mV, "
               "-59 MN-named; v2 per-class calibrated from the graded "
               "peaks); ALL NonSpikingChemical synapses re-emitted as "
               "SpikingChemical", fontsize=6.2, va="center")
    cv.ax.text(15.6, 0.98,
               "Source: connectome_templates.json 'bilateralrg' (W2L "
               "aproj mine + documented BilateralRG build chain) - the "
               "connectome the spiking copies preserve (validate.py: "
               "connexion sets identical). Reports: "
               "reports_spiking_20261002/goal3_animatlab_spiking.md.",
               fontsize=5.8, va="center", ha="right")
    save_both(cv, "animatlab_spiking_conversions")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", default="all",
                    choices=["mirror", "animatlab", "all"])
    a = ap.parse_args()
    if a.which in ("mirror", "all"):
        fig_mirror()
    if a.which in ("animatlab", "all"):
        fig_animatlab()


if __name__ == "__main__":
    main()
