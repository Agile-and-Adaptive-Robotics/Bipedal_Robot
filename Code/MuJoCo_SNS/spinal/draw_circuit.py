"""Publication-grade schematics of the gait2392 spinal circuit.

Style follows Ben's reference pictures (Rybak / SNS-hierarchy style,
2026-09-12):
  * tinted layer bands with plain labels (no boxes everywhere)
  * populations drawn as clusters of small circles, real counts pulled
    LIVE from muscle_map.py (never hardcoded)
  * size hierarchy: half-centers big, pattern cells medium, sensors and
    local interneurons small; small loop cells for adaptation
  * strict rows, mirrored left/right around a dashed midline in the
    "full" panel, commissural wires gathered in one pale central corridor
  * wires colored by source (Okabe-Ito CVD-safe); synapse markers keep
    Ben's conventions: white triangle with its BASE flat against the
    excited cell (snsfig.m 2026-09-09 inversion), solid black dot =
    inhibitory. Triangle bases are aligned with the wire's ARRIVAL
    TANGENT so curved wires stay attached.
  * clean wires: no conductance numbers on schematics -- numbers live in
    the weights figure

Figures (--which):
  core     one side, hierarchy: DRIVE/POSTURE -> persistent-Na RG
           half-centers -> laminated PF cells -> MN pool clusters -> muscles; Ia/II/Ib
           sensor triplet (x46 each) + IB-EXC load-sharing cells.
  full     both sides drawn fully (no ghost half): mirrored rows,
           midline commissural corridor, BAL family, trunk pools.
  deng     Deng/Nourse-Fig-6A-style layered schematic, edges read back
           from the COMPILED laminated network (assert contract).
  weights  W_PF_MN phase x group heatmap + W_POSTURE column (numbers).

Numbers are never hardcoded from memory: --source composites the live
config exactly like `runner --fitted --best`; pool counts come from
muscle_map.py. Every figure carries a gray source note naming the files.

Usage: python draw_circuit.py [--which core|full|weights|all]
                             [--source params|fitted|best]
                             [--fmt pdf,svg,png]
Outputs into figures/.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import (Circle, Ellipse, FancyArrowPatch, FancyBboxPatch,
                                Polygon)

import params
import muscle_map

HERE = Path(__file__).parent
OUT = HERE / "figures"

# ---------------------------------------------------------------- source
def effective_tables(source: str) -> dict:
    """Composite the live configuration exactly as runner --fitted --best
    does. Returns G, W (W_PF_MN), WPOST, TAU, note."""
    G = dict(params.G)
    W = {ph: dict(tbl) for ph, tbl in params.W_PF_MN.items()}
    WPOST = dict(params.W_POSTURE)
    TAU = dict(params.TAU)
    note = "params.py defaults"
    fit_keys = None
    if source in ("fitted", "best"):
        fit = json.loads((HERE / "fitted_walk_params.json").read_text("utf-8"))
        for ph, tbl in fit["W_PF_MN"].items():
            for g, w in tbl.items():
                W[ph][g] = float(w)
        for g, w in fit["W_POSTURE"].items():
            WPOST[g] = float(w)
        # pf_gain scales ONLY fitted-file entries (params-default trunk
        # keys are not part of the back-solved table; matching runner)
        fit_keys = ({ph: set(tbl) for ph, tbl in fit["W_PF_MN"].items()},
                    set(fit["W_POSTURE"]))
        note = "fitted_walk_params.json (IK/NNLS back-solve refit)"
    if source == "best":
        best = json.loads((HERE / "best_walk_params.json").read_text("utf-8"))
        gain = float(best.get("pf_gain", 1.0))
        if gain != 1.0:
            for ph in W:
                for g in fit_keys[0].get(ph, ()):
                    W[ph][g] *= gain
            for g in fit_keys[1]:
                WPOST[g] *= gain
        p = best["params"]
        W["E2"]["ankle_pf"] = p["e2_pf"]
        W["F1"]["ankle_df"] = p["f1_df"]
        W["F1"]["knee_flex"] = p["f1_kf"]
        WPOST["knee_ext"] = p["post_kneext"]
        if "post_hipext" in p:
            WPOST["hip_ext"] = p["post_hipext"]
        G["descend_to_rg_e"] = p["desc_e"]
        G["rg_to_pf"] = p["rg_to_pf"]
        TAU["rg_adapt"] = p["rg_adapt"]
        note = (f"best_walk_params.json (study {best.get('study', '?')}, "
                f"pf_gain {gain:.2f}) on fitted refit")
    return dict(G=G, W=W, WPOST=WPOST, TAU=TAU, note=note)


# ------------------------------------------------- live pool counts
def live_counts() -> dict:
    """Per-group primary/secondary muscle lists from muscle_map.py."""
    prim, sec = {}, {}
    for g in muscle_map.GROUPS:
        prim[g] = sorted(b for b, gs in muscle_map._GROUPS_BY_NAME.items()
                         if gs[0] == g)
        sec[g] = sorted(b for b, gs in muscle_map._GROUPS_BY_NAME.items()
                        if g in gs[1:])
    return dict(prim=prim, sec=sec)


CNT = live_counts()
N_PER_SIDE = sum(len(v) for v in CNT["prim"].values())          # 46
N_MN_TOTAL = 2 * N_PER_SIDE                                      # 92
# structural cell counts per side (build_network.py layout, NaP
# architecture 2026-09-16: RG = 2 NaP half-centers + InE/InF + CIN-E/F;
# ADAP and PRESET/PREA retired; PFA removed)
N_RG, N_ADAP, N_PF, N_PFA, N_IBEXC = 6, 0, 6, 0, 5
N_SHARED = 8                       # DRIVE,POSTURE,BAL_PF,BAL_DF,BAL_TRKx2,BAL_LATx2
N_SIDE = N_PER_SIDE * 4 + N_RG + N_ADAP + N_PF + N_PFA + N_IBEXC   # 201
N_TOTAL = 2 * N_SIDE + N_SHARED                                      # 410
N_INPUTS = 4 * N_MN_TOTAL + N_SHARED                                 # 376

ABBREV = {"hip_ext": "HIP-E", "hip_flex": "HIP-F", "hip_abd": "HIP-AB",
          "hip_add": "HIP-AD", "knee_ext": "KNEE-E", "knee_flex": "KNEE-F",
          "ankle_pf": "ANK-PF", "ankle_df": "ANK-DF",
          "trunk_ext": "TRK-E", "trunk_flex": "TRK-F"}
# display order: sagittal chain first, then frontal, then trunk
GROUP_ORDER = ("hip_ext", "hip_flex", "knee_ext", "knee_flex",
               "ankle_pf", "ankle_df", "hip_abd", "hip_add",
               "trunk_ext", "trunk_flex")

# Okabe-Ito CVD-safe accents (tints only; nodes stay white/light)
OI_SKY, OI_VERM, OI_BLUE = "#56B4E9", "#D55E00", "#0072B2"
OI_GREEN, OI_PURPLE, OI_ORANGE = "#009E73", "#CC79A7", "#E69F00"
TINT = 0.25
# wire colors: by source family / sign
W_E, W_F = OI_BLUE, OI_VERM          # extensor / flexor (Rybak convention)
W_EXC, W_INH, W_DESC = OI_GREEN, OI_PURPLE, OI_SKY   # other exc / inh / descending

W_MIN = 0.03        # draw PF->MN edges down to this weight


def _tint(hex_color: str) -> str:
    r = int(hex_color[1:3], 16) / 255
    g = int(hex_color[3:5], 16) / 255
    b = int(hex_color[5:7], 16) / 255
    f = 1 - TINT
    return f"#{int((f * r + TINT) * 255):02x}{int((f * g + TINT) * 255):02x}" \
           f"{int((f * b + TINT) * 255):02x}"


# ---------------------------------------------------------------- glyphs
R_MED = 0.30                    # PF cells
R_SMALL = 0.18                  # ADAP/PFA/IB-EXC
R_BIG = 0.44                    # RG half-centers
R_POOL = 0.14                   # circles inside an MN pool cluster
R_SENS = 0.16                   # Ia/II/Ib
POOL_GAP = 0.34                 # pitch between pool circles
LBL_BBOX = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85)


class Canvas:
    """Circuit canvas: full-bleed coordinates, axis off."""

    def __init__(self, w, h):
        self.fig, self.ax = plt.subplots(figsize=(w, h))
        self.ax.set_xlim(0, w)
        self.ax.set_ylim(0, h)
        self.ax.axis("off")
        self.fig.subplots_adjust(left=0.004, right=0.996, top=0.996,
                                 bottom=0.004)

    def save(self, name, fmts):
        outs = []
        for fmt in fmts:
            p = OUT / f"{name}.{fmt}"
            tmp = p.with_name(f".{p.stem}.{os.getpid()}.tmp{p.suffix}")
            self.fig.savefig(tmp, dpi=300, bbox_inches="tight", pad_inches=0.04)
            last_error = None
            for _ in range(20):
                try:
                    os.replace(tmp, p)
                    last_error = None
                    break
                except OSError as exc:
                    last_error = exc
                    time.sleep(0.15)
            if last_error is not None:
                raise OSError(
                    f"could not replace locked figure {p}; completed render "
                    f"is preserved at {tmp}") from last_error
            outs.append(p.name)
        plt.close(self.fig)
        return outs


def neuron(cv, x, y, label, sub=None, fc="white", r=R_MED, ghost=False,
           fs=None, lw=1.5):
    cv.ax.add_patch(Circle((x, y), r, fc="0.97" if ghost else fc,
                           ec="0.55" if ghost else "black",
                           lw=lw, zorder=3))
    if ghost:
        return (x, y, r)
    fs = fs if fs else min(11.0, max(6.0, 9.5 * r / R_MED))
    cv.ax.text(x, y, label, ha="center", va="center", fontsize=fs, zorder=4)
    if sub:
        cv.ax.text(x, y - r - 0.17, sub, ha="center", va="center",
                   fontsize=fs - 0.8, style="italic", color="0.35", zorder=4)
    return (x, y, r)


def input_box(cv, x, y, w, h, label, sub=None):
    cv.ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                   boxstyle="round,pad=0.02",
                                   fc="0.94", ec="black", lw=1.2, zorder=3))
    cv.ax.text(x, y + (0.09 if sub else 0), label, ha="center",
               va="center", fontsize=7.8, zorder=4)
    if sub:
        cv.ax.text(x, y - 0.13, sub, ha="center", va="center",
                   fontsize=6.4, style="italic", color="0.35", zorder=4)
    return (x, y, max(w, h) * 0.55)


def muscle(cv, x, y, w=0.72, h=0.26):
    cv.ax.add_patch(Ellipse((x, y), w, h, fc=_tint(OI_SKY), ec="black",
                            lw=1.0, zorder=3))
    return (x, y, h / 2)


def cluster_hw(group, cap=4):
    k = min(len(CNT["prim"][group]), cap)
    return (k - 1) / 2 * POOL_GAP + R_POOL


def pool_cluster(cv, x, y, group, cap=4):
    """MN pool cluster: small circles = real motoneuron pools (primary
    members). Biarticular secondary members are summarized in the count
    text "(+n)". Returns anchor dict."""
    n_prim = len(CNT["prim"][group])
    n_sec = len(CNT["sec"][group])
    k = min(n_prim, cap)
    xs = [x + (i - (k - 1) / 2) * POOL_GAP for i in range(k)]
    for xc in xs:
        cv.ax.add_patch(Circle((xc, y), R_POOL, fc=_tint(OI_SKY),
                               ec="black", lw=0.9, zorder=3))
    cnt_txt = f"x{n_prim}" + (f"\n(+{n_sec})" if n_sec else "")
    cv.ax.text(x, y - 0.48, f"{ABBREV[group]}\n{cnt_txt}", ha="center",
               va="center", fontsize=6.4, zorder=4)
    hw = cluster_hw(group, cap)
    return dict(c=(x, y), top=(x, y + R_POOL + 0.05),
                left=(x - hw, y), right=(x + hw, y),
                bot=(x, y - R_POOL - 0.05), hw=hw)


def layout_pools(cv, x0, y, groups, cap=4, gap=0.24, min_pitch=1.05):
    """Width-aware single-row pool layout with a minimum pitch so
    two-line labels never collide; returns {group: anchors}."""
    pools, x = {}, x0 + cluster_hw(groups[0], cap)
    for i, g in enumerate(groups):
        pools[g] = pool_cluster(cv, x, y, g, cap=cap)
        if i + 1 < len(groups):
            x += max(cluster_hw(g, cap) + gap + cluster_hw(groups[i + 1], cap),
                     min_pitch)
    return pools


def layer_band(cv, x0, x1, y0, y1, hex_color, label, lab_dx=0.18,
               lab_dy=0.22, fs=8.0):
    cv.ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                   boxstyle="round,pad=0.02,rounding_size=0.18",
                                   fc=_tint(hex_color), ec="none", zorder=0.5))
    if label:
        cv.ax.text(x0 + lab_dx, y1 - lab_dy, label, fontsize=fs, ha="left",
                   va="top", color="0.30", style="italic", zorder=1)


def syn(cv, p, q, exc: bool, label=None, lpos=0.5, loff=(0.1, 0.1), lw=1.4,
        rad=0.0, ghost=False, dashed=False, shrink_a=None, shrink_b=None,
        color=None, lbox=True, lfs=6.4, record=None):
    """Edge p->q. Excitatory: WHITE TRIANGLE WITH ITS BASE FLAT AGAINST THE
    TARGET CELL, apex pointing back along the wire (Ben's snsfig.m
    convention, 2026-09-09). Inhibitory: solid black dot at the target.
    The base/dot are aligned with the wire's arrival tangent (arc3
    control point geometry), so curved wires stay attached. p/q may be
    (x, y) or the (x, y, r) tuples returned by glyph helpers.
    record: optional list; the wire's true geometric start/end points
    (after shrink) are appended as dicts for boundary-contract checks."""
    col = "0.55" if ghost else (color if color else "0.15")
    ra = p[2] if len(p) > 2 and shrink_a is None else (shrink_a if shrink_a
                                                       is not None else 0.0)
    rb = q[2] if len(q) > 2 and shrink_b is None else (shrink_b if shrink_b
                                                       is not None else 0.0)
    pa, qa = np.asarray(p[:2], float), np.asarray(q[:2], float)
    d = qa - pa
    L = max(np.linalg.norm(d), 1e-9)
    u = d / L
    perp = np.array([-u[1], u[0]])
    a = pa + u * (ra + 0.02)
    # arrival tangent of the arc3 bezier at q (matplotlib bulges the arc
    # to the CW side of the chord for positive rad -- verified against
    # the over-the-top arcs note in the previous revision)
    ctrl = 0.5 * (pa + qa) - perp * (rad * L)
    ut = qa - ctrl
    ut = ut / max(np.linalg.norm(ut), 1e-9)
    pt = np.array([-ut[1], ut[0]])
    base = dot = None
    if exc:
        base = qa - ut * (rb + 0.02)
        tri = min(0.20, max(0.09, 0.45 * (L - ra - rb)))
        apex = base - ut * tri
        cv.ax.add_patch(FancyArrowPatch(
            tuple(a), tuple(apex), arrowstyle="-", lw=lw, color=col, zorder=2,
            connectionstyle=f"arc3,rad={rad}",
            linestyle=(0, (3, 2)) if dashed else "solid"))
        if not ghost:
            cv.ax.add_patch(Polygon(
                [tuple(base + pt * 0.095), tuple(base - pt * 0.095),
                 tuple(apex)], closed=True, fc="white", ec=col, lw=1.2,
                zorder=3))
    else:
        dot = qa - ut * (rb + 0.10)
        cv.ax.add_patch(FancyArrowPatch(
            tuple(a), tuple(dot), arrowstyle="-", lw=lw, color=col, zorder=2,
            connectionstyle=f"arc3,rad={rad}",
            linestyle=(0, (3, 2)) if dashed else "solid"))
        if not ghost:
            cv.ax.add_patch(Circle(tuple(dot), 0.085, fc="black", ec="black",
                                   zorder=3))
    if record is not None and not ghost:
        record.append(dict(start=a.copy(), end=(base if exc else dot).copy(),
                           exc=bool(exc), rad=float(rad),
                           dashed=bool(dashed)))
    if label and not ghost:
        mid = pa + (qa - pa) * lpos
        if rad:
            mid = mid + perp * rad * L * 0.5
        mid = mid + np.asarray(loff, float)
        cv.ax.text(*mid, label, fontsize=lfs, ha="center", va="center",
                   style="italic", color="0.25", zorder=5,
                   bbox=LBL_BBOX if lbox else None)


def lane(cv, x, y_top, y_bot, color=W_DESC):
    """A quiet dashed margin lane (vertical) for descending/bias wires."""
    cv.ax.plot([x, x], [y_bot, y_top], ls=(0, (3, 2)), color=color, lw=1.0,
               zorder=1.5)


def tag(cv, x, y, text, fs=6.6, color="0.30", ha="center", box=False):
    cv.ax.text(x, y, text, fontsize=fs, ha=ha, va="center", style="italic",
               color=color, zorder=6 if box else 4,
               bbox=LBL_BBOX if box else None)


def legend_row(cv, y, x0=0.6):
    ax = cv.ax
    ax.add_patch(Circle((x0, y), 0.16, fc="white", ec="black", lw=1.3))
    ax.text(x0 + 0.26, y, "neuron", fontsize=7.4, va="center")
    x = x0 + 1.35
    ax.add_patch(FancyArrowPatch((x, y), (x + 0.55, y), arrowstyle="-",
                                 lw=1.3, color=W_EXC))
    ax.add_patch(Polygon([(x + 0.85, y + 0.085), (x + 0.85, y - 0.085),
                          (x + 0.62, y)], closed=True, fc="white",
                         ec=W_EXC, lw=1.2))
    ax.text(x + 1.0, y, "excitatory", fontsize=7.4, va="center")
    x = x0 + 3.6
    ax.add_patch(FancyArrowPatch((x, y), (x + 0.62, y), arrowstyle="-",
                                 lw=1.3, color=W_INH))
    ax.add_patch(Circle((x + 0.72, y), 0.085, fc="black"))
    ax.text(x + 0.95, y, "inhibitory", fontsize=7.4, va="center")
    x = x0 + 6.1
    ax.add_patch(Circle((x, y), R_POOL, fc=_tint(OI_SKY), ec="black", lw=0.9))
    ax.add_patch(Circle((x + 0.34, y), R_POOL, fc=_tint(OI_SKY), ec="black",
                        lw=0.9))
    ax.text(x + 0.75, y, "MN pool cluster", fontsize=7.4, va="center")
    x = x0 + 9.0
    ax.add_patch(Ellipse((x + 0.2, y), 0.55, 0.22, fc=_tint(OI_SKY),
                         ec="black", lw=1.0))
    ax.text(x + 0.75, y, "muscle", fontsize=7.4, va="center")


def wire_color_key(cv, x, y):
    ax = cv.ax
    items = [(W_E, "extensor-half drive"), (W_F, "flexor-half drive"),
             (W_EXC, "other excitatory"), (W_INH, "inhibitory"),
             (W_DESC, "descending / balance")]
    for i, (c, name) in enumerate(items):
        xx = x + i * 2.35
        ax.add_patch(FancyArrowPatch((xx, y), (xx + 0.5, y), arrowstyle="-",
                                     lw=1.6, color=c))
        ax.text(xx + 0.6, y, name, fontsize=6.6, va="center", color="0.25")


def key_box(cv, x, y, width=50):
    """Group key: abbreviation -> real member muscles (live), wrapped."""
    lines = ["MN pool key (live from muscle_map.py)"]
    for g in GROUP_ORDER:
        head = f"{ABBREV[g]} ({len(CNT['prim'][g])}): "
        words = (", ".join(CNT["prim"][g])
                 + (f" (+{', '.join(CNT['sec'][g])} sec)"
                    if CNT["sec"][g] else "")).split(" ")
        line = head
        for wd in words:
            if len(line) + len(wd) > width:
                lines.append(line)
                line = " " * (len(head) + 2) + wd
            else:
                line += wd + " "
        lines.append(line)
    txt = "\n".join(lines)
    cv.ax.text(x, y, txt, fontsize=5.5, va="top", ha="left", color="0.25",
               family="monospace", zorder=4, linespacing=1.15,
               bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.75",
                         lw=0.7))


def source_note(cv, x, y, t, extra=True):
    tot = (f"network size (live): {N_PER_SIDE} MN + {N_PER_SIDE}x(3) sensory "
           f"per side, +{N_RG + N_ADAP + N_PF + N_PFA} rhythm/pattern cells "
           f"+ {N_IBEXC} load-sharing INs per side; both sides + "
           f"{N_SHARED} shared = {N_TOTAL} neurons, {N_INPUTS} inputs"
           ) if extra else ""
    cv.ax.text(x, y, t + ("   |   " + tot if tot else ""), fontsize=6.2,
               ha="left", va="center", color="0.40", style="italic")


# ------------------------------------------------- V-class annotations
# Correspondence of our circuit elements to the genetically identified
# locomotor-CPG interneuron classes, grounded in Shevtsova et al. 2026
# (eLife RP107480, Fig. 2 schematic: InF/InE, V2a, V0V, V0D, V3-E/F, Ini,
# InE1) and Rybak et al. 2015 (eNeuro review: V0D/V0V/V0C, V1, V2a, V2b,
# V3, dI6 roles). OUR circuit has NO V3 analog: every cross-side
# connection is inhibitory.
VCLASS_REF = ("Shevtsova et al. 2026 eLife RP107480 Fig. 2; "
              "Rybak et al. 2015 eNeuro")


def vclass_strip(cv):
    """Top strip: class-by-class mapping, one column per V-class."""
    y0 = 14.05
    cv.ax.plot([0.55, 18.45], [y0 + 0.52, y0 + 0.52], color="0.75", lw=0.7)
    cv.ax.text(0.55, y0 + 0.74, "V-class correspondence of this circuit "
               "(--vclasses; cf. " + VCLASS_REF + ")",
               fontsize=7.0, ha="left", va="center", color="0.25",
               style="italic")
    cols = [
        (0.55, "V0$_D$ (V0$_c$) / dI6  —  CIN, inhibitory",
         "cross-side RG inhibition (F$\\leftrightarrow$F strong,\n"
         "E$\\leftrightarrow$E weak): left–right alternation\n"
         "= Shevtsova '26 RG-F$\\to$V0D$\\dashrightarrow$c-RG-F (−0.07)"),
        (5.05, "V1 / V2b  —  ipsilateral, inhibitory",
         "RG-E$\\leftrightarrow$RG-F half-center + PF reciprocal\n"
         "inhibition + Ia reciprocal inhibition onto MNs\n"
         "(flexor–extensor alternation; Renshaw-like)"),
        (10.05, "V2$_a$  —  ipsilateral, excitatory",
         "RG$\\to$PF and PF$\\to$MN excitatory relay\n"
         "= Shevtsova '26 RG-F$\\to$V2a$\\to$V0V chain\n"
         "(also drives MN pools in the full models)"),
        (14.55, "V3  —  CIN, excitatory",
         "ABSENT in this circuit: every cross-side\n"
         "connection is inhibitory (their model uses\n"
         "V3-E/F for left–right synchrony at speed)"),
    ]
    for x, head, body in cols:
        cv.ax.text(x, y0 + 0.30, head, fontsize=7.2, ha="left", va="center",
                   color="black", family="sans-serif", weight="bold")
        cv.ax.text(x, y0 - 0.18, body, fontsize=6.3, ha="left", va="top",
                   color="0.30")
    note_txt = ("we have NO V3 analog — all cross-side connections are "
                "inhibitory (alternation only, no synchrony pathway)")
    cv.ax.text(18.45, y0 + 0.74, note_txt, fontsize=6.6, ha="right",
               va="center", color="0.45", style="italic")


# ------------------------------------------------------- shared blocks
def rg_block(cv, x, y, suffix="", big=True, mirror=False):
    """Half-centers + adaptation loop cells. Returns (rg_e, rg_f)."""
    r = R_BIG if big else R_MED
    dx = -1.55 if not mirror else 1.55
    rge = neuron(cv, x + dx, y + 0.62, f"RG-E{suffix}", fc=_tint(W_E), r=r)
    rgf = neuron(cv, x + dx, y - 0.78, f"RG-F{suffix}", fc=_tint(W_F), r=r)
    syn(cv, rge, rgf, False, rad=0.25, color=W_INH)
    syn(cv, rgf, rge, False, rad=0.25, color=W_INH)
    sgn = 1 if not mirror else -1
    tag(cv, x + dx + sgn * 1.25, y - 0.08, "mutual\ninhibition", fs=6.0,
        box=True)
    adap_x = x + dx - sgn * 1.75
    ade = neuron(cv, adap_x, y + 0.62, "ADAP-E", r=R_SMALL, fs=6.0)
    adf = neuron(cv, adap_x, y - 0.78, "ADAP-F", r=R_SMALL, fs=6.0)
    syn(cv, rge, ade, True, rad=-0.35, lw=1.0, color=W_E)
    syn(cv, ade, rge, False, rad=-0.35, lw=1.0, color=W_INH)
    syn(cv, rgf, adf, True, rad=0.35, lw=1.0, color=W_F)
    syn(cv, adf, rgf, False, rad=0.35, lw=1.0, color=W_INH)
    tag(cv, adap_x, y - 0.08, "adaptation\n(fatigue)", fs=5.8, box=True)
    return rge, rgf


def pf_block(cv, x, y, suffix="", mirror=False):
    """Four PF phase cells + PFA loop cells + reciprocal inhibition.
    Returns dict of PF cells."""
    names = (("E1", 1.35, W_E), ("E2", 0.45, W_E),
             ("F1", -0.5, W_F), ("F2", -1.35, W_F))
    out = {}
    for nm, dy, c in names:
        out[nm] = neuron(cv, x, y + dy, f"PF-{nm}{suffix}", fc=_tint(c))
    syn(cv, out["E2"], out["F1"], False, color=W_INH, rad=0.12)
    syn(cv, out["F1"], out["E2"], False, color=W_INH, rad=0.12)
    syn(cv, out["F2"], out["E2"], False, color=W_INH, rad=-0.35)
    syn(cv, out["E2"], out["F2"], False, color=W_INH, rad=-0.35)
    syn(cv, out["E1"], out["F1"], False, color=W_INH, rad=-0.3)
    sgn = 1 if not mirror else -1
    for nm in out:
        pfa = neuron(cv, x + sgn * 1.05, out[nm][1], "", r=0.13,
                     fc="white", lw=0.9)
        cv.ax.text(pfa[0] + sgn * 0.42, pfa[1], f"PFA-{nm}",
                   fontsize=5.4, ha="center", va="center", color="0.35",
                   zorder=6, bbox=LBL_BBOX)
        syn(cv, out[nm], pfa, True, lw=0.8,
            color=W_E if nm[0] == "E" else W_F)
        syn(cv, pfa, out[nm], False, lw=0.8, color=W_INH, rad=0.25)
    return out


def sensor_triplet(cv, x, y, suffix="", sub=True):
    """Ia/II/Ib population circles (x46 each per side)."""
    ia = neuron(cv, x - 1.5, y, f"Ia{suffix}", f"x{N_PER_SIDE}", r=R_SENS,
                fc=_tint(OI_ORANGE), fs=6.6)
    ii = neuron(cv, x, y, f"II{suffix}", f"x{N_PER_SIDE}", r=R_SENS,
                fc=_tint(OI_ORANGE), fs=6.6)
    ib = neuron(cv, x + 1.5, y, f"Ib{suffix}", f"x{N_PER_SIDE}", r=R_SENS,
                fc=_tint(OI_ORANGE), fs=6.6)
    if sub:
        tag(cv, x, y - 0.72, f"one Ia, II, Ib per pool "
            f"({N_PER_SIDE} of each per side)", fs=6.0, box=True)
    return dict(ia=ia, ii=ii, ib=ib)


def ibexc_cluster(cv, x, y, short=False):
    """The five load-sharing interneurons (one per extensor-stance group).
    Stance-gated by the RG-E signal (noted in tag; wire omitted for
    clarity)."""
    xs = [x + (i - 2) * 0.42 for i in range(N_IBEXC)]
    for xc in xs:
        cv.ax.add_patch(Circle((xc, y), R_SMALL, fc=_tint(OI_GREEN),
                               ec="black", lw=0.9, zorder=3))
    cv.ax.text(x, y, "IB-EXC", ha="center", va="center", fontsize=5.4,
               zorder=4)
    if short:
        tag(cv, x, y - 0.36, "x5: one per extensor-stance group",
            fs=5.8, box=True)
    else:
        tag(cv, x, y - 0.46, f"{N_IBEXC} cells: one per extensor-stance "
            "group;\nRG-E stance-gated Ib reversal (load sharing)", fs=6.0,
            box=True)
    return (x, y, R_SMALL)


# ------------------------------------------------------------ fig: core
def fig_core(t, fmts):
    W = t["W"]
    cv = Canvas(15.6, 14.0)

    # ---------------- inputs (top-left) ----------------
    drv = input_box(cv, 1.35, 12.75, 1.9, 0.62, "DRIVE", "MLR surrogate")
    pos = input_box(cv, 4.05, 12.75, 1.9, 0.62, "POSTURE", "tonic")
    post_i = input_box(cv, 13.9, 12.75, 1.9, 0.62, "POST$_i$",
                       "standing bias")

    # POST_i margin lane down the right side (outside all bands)
    lane(cv, 14.35, 12.55, 4.30)

    # ---------------- rhythm generator ----------------
    layer_band(cv, 0.5, 8.6, 9.6, 12.0, OI_ORANGE,
               "RHYTHM GENERATOR (RG)")
    rge, rgf = rg_block(cv, 4.3, 10.85)
    syn(cv, drv, rge, True, color=W_DESC, rad=0.1)
    syn(cv, drv, rgf, True, color=W_DESC, rad=-0.15)
    syn(cv, pos, rge, True, color=W_DESC, rad=-0.2)
    tag(cv, 7.3, 11.55, "frequency rises with DRIVE;\nstance-biased drive +\n"
        "adaptation set the duty", box=True)

    # ---------------- pattern formation ----------------
    layer_band(cv, 0.5, 8.6, 6.05, 9.45, OI_SKY,
               "PATTERN FORMATION (PF)  --  phase-window cells")
    pfs = pf_block(cv, 4.3, 7.75)
    syn(cv, rge, pfs["E1"], True, color=W_E, rad=-0.08)
    syn(cv, rge, pfs["E2"], True, color=W_E, rad=-0.12)
    syn(cv, rgf, pfs["F1"], True, color=W_F, rad=-0.08)
    syn(cv, rgf, pfs["F2"], True, color=W_F, rad=-0.12)
    syn(cv, drv, pfs["E1"], True, color=W_DESC, rad=0.3, dashed=True)
    tag(cv, 7.3, 8.05, "PF cells = phase windows;\n"
        "(tau, adaptation) stagger\ntheir bursts", box=True)

    # ---------------- interneurons + sensors (flanks) ----------------
    ib = ibexc_cluster(cv, 10.0, 5.35)
    sens = sensor_triplet(cv, 12.6, 5.5)

    # ---------------- motoneuron pools ----------------
    layer_band(cv, 0.5, 15.1, 1.6, 3.95, OI_GREEN,
               "MOTONEURON POOLS  --  one cluster per functional group, "
               "real pool counts")
    pools = layout_pools(cv, 0.7, 3.0, GROUP_ORDER, cap=4, gap=0.26)
    for ph, node in pfs.items():
        for g, w in W[ph].items():
            if g in pools and w >= W_MIN:
                syn(cv, node, pools[g]["top"], True,
                    lw=0.7 + min(w, 0.25) * 6, shrink_b=0.05,
                    color=W_E if ph[0] == "E" else W_F)
    # POST_i lane arrows into the right-end pools + tag
    syn(cv, (14.35, 4.30), pools["ankle_df"]["top"], True, color=W_DESC,
        dashed=True, rad=-0.25, shrink_a=0.0)
    syn(cv, (14.35, 4.30), pools["trunk_ext"]["top"], True, color=W_DESC,
        dashed=True, rad=0.1, shrink_a=0.0)
    tag(cv, 12.1, 4.55, "POST$_i$: solved standing bias\ninto every pool "
        "(dashed)", ha="center", box=True)
    # reflex pathways (representative; every pool has Ia/II/Ib)
    syn(cv, sens["ia"], pools["knee_ext"]["top"], True, color=W_EXC,
        rad=-0.25, label="Ia homonymous", lpos=0.45, loff=(0.0, 0.55))
    syn(cv, sens["ia"], (pools["knee_flex"]["c"][0] + 0.30,
                         pools["knee_flex"]["c"][1] + 0.28), False,
        color=W_INH, rad=0.2, shrink_b=0.0,
        label="Ia reciprocal", lpos=0.55, loff=(0.3, 0.5))
    syn(cv, sens["ii"], pools["hip_ext"]["top"], True, color=W_EXC,
        rad=0.25, label="II (stance-gated)", lpos=0.45, loff=(0.0, 0.6))
    syn(cv, sens["ib"], (pools["ankle_pf"]["c"][0] + 0.32,
                         pools["ankle_pf"]["c"][1] + 0.26), False,
        color=W_INH, rad=-0.2, shrink_b=0.0,
        label="Ib autogenic", lpos=0.5, loff=(0.35, 0.55))
    syn(cv, sens["ib"], ib, True, color=W_EXC, rad=-0.15, lw=1.0)
    syn(cv, ib, pools["hip_ext"]["top"], True, color=W_EXC, rad=0.2, lw=1.0,
        label="Ib reversal $\\rightarrow$ extensors", lpos=0.55,
        loff=(-0.1, 0.5))
    # muscles row
    mx = 0.7 + cluster_hw(GROUP_ORDER[0], 4)
    for i, g in enumerate(GROUP_ORDER):
        muscle(cv, mx, 1.1)
        syn(cv, pools[g]["bot"], (mx, 1.24), True, color="0.4", lw=0.8,
            shrink_a=0.04, shrink_b=0.10)
        if i + 1 < len(GROUP_ORDER):
            mx += (cluster_hw(g, 4) + 0.26
                   + cluster_hw(GROUP_ORDER[i + 1], 4))

    # ---------------- legend / key / source ----------------
    legend_row(cv, 0.68)
    wire_color_key(cv, 0.7, 0.25)
    key_box(cv, 9.6, 8.7)
    source_note(cv, 0.55, 13.42, f"source: {t['note']}")
    return cv.save("circuit_core", fmts)


# ------------------------------------------------------------ fig: full
ROW1 = ("hip_ext", "hip_flex", "knee_ext", "knee_flex", "ankle_pf")
ROW2 = ("ankle_df", "hip_abd", "hip_add", "trunk_ext", "trunk_flex")


def fig_full(t, fmts, vc=False):
    W = t["W"]
    cv = Canvas(19.0, 14.95 if vc else 13.4)
    MID = 9.5

    # midline + commissural corridor (RG row height only)
    cv.ax.plot([MID, MID], [1.4, 11.9], ls=(0, (4, 3)), color="0.6", lw=1.0,
               zorder=1)
    layer_band(cv, MID - 0.55, MID + 0.55, 8.6, 11.4, OI_PURPLE, "")
    cv.ax.text(MID, 11.65, "interleg commissural inhibition",
               fontsize=7.0, ha="center", va="center", color="0.35",
               style="italic", zorder=1)
    tag(cv, MID, 9.3, "F$\\leftrightarrow$F strong\nE$\\leftrightarrow$E weak",
        fs=5.8, box=True)

    # shared inputs (top)
    drv = input_box(cv, MID - 1.9, 12.9, 1.6, 0.6, "DRIVE", "MLR surrogate")
    pos = input_box(cv, MID + 1.9, 12.9, 1.6, 0.6, "POSTURE", "+ POST$_i$")

    sides = {"r": dict(x0=10.4, x1=18.6, mirror=False, sfx="$_r$",
                       lane=18.8),
             "l": dict(x0=0.4, x1=8.6, mirror=True, sfx="$_l$", lane=0.2)}
    rg, pfs, sens, ib, pools = {}, {}, {}, {}, {}
    for sd, cfg in sides.items():
        cx = (cfg["x0"] + cfg["x1"]) / 2
        outer = cfg["x1"] if not cfg["mirror"] else cfg["x0"]
        # balance bus lane: OUTSIDE the side bands
        lane(cv, cfg["lane"], 12.6, 4.3)
        b1 = input_box(cv, outer - (1.05 if not cfg["mirror"] else -1.05),
                       12.75, 1.6, 0.55, "BAL$_{PF/DF}$", "sagittal COM PD")
        b2 = input_box(cv, outer - (1.05 if not cfg["mirror"] else -1.05),
                       12.1, 1.6, 0.55, "BAL$_{LAT}$/TRK", "abduct/trunk")
        # RG
        layer_band(cv, cfg["x0"], cfg["x1"], 9.35, 11.75, OI_ORANGE,
                   f"RHYTHM GENERATOR {sd.upper()}")
        rge, rgf = rg_block(cv, cx, 10.6, suffix=cfg["sfx"],
                            mirror=cfg["mirror"])
        rg[sd] = (rge, rgf)
        syn(cv, drv, rge, True, color=W_DESC, rad=0.15)
        syn(cv, drv, rgf, True, color=W_DESC, rad=-0.1)
        syn(cv, pos, rge, True, color=W_DESC, rad=-0.25)
        # PF
        layer_band(cv, cfg["x0"], cfg["x1"], 5.85, 9.25, OI_SKY,
                   f"PATTERN FORMATION {sd.upper()}")
        pf = pf_block(cv, cx - (0.9 if not cfg["mirror"] else -0.9), 7.55,
                      suffix=cfg["sfx"], mirror=cfg["mirror"])
        pfs[sd] = pf
        syn(cv, rge, pf["E1"], True, color=W_E)
        syn(cv, rge, pf["E2"], True, color=W_E, rad=-0.08)
        syn(cv, rgf, pf["F1"], True, color=W_F)
        syn(cv, rgf, pf["F2"], True, color=W_F, rad=-0.08)
        # interneurons + sensors: outer flank of the side
        layer_band(cv, cfg["x0"], cfg["x1"], 4.3, 5.75, OI_GREEN,
                   f"IN + SENSORS {sd.upper()}")
        x_in = outer - (2.2 if not cfg["mirror"] else -2.2)
        ib[sd] = ibexc_cluster(cv, x_in, 5.5, short=True)
        sens[sd] = sensor_triplet(cv, x_in - (0.9 if not cfg["mirror"]
                                              else -0.9), 4.85,
                                  suffix=cfg["sfx"], sub=False)
        tag(cv, x_in - (0.9 if not cfg["mirror"] else -0.9), 4.5,
            f"Ia/II/Ib x{N_PER_SIDE} each", fs=5.8, box=True)
        # MN pools: two rows of five per side
        layer_band(cv, cfg["x0"], cfg["x1"], 1.35, 4.15, OI_SKY,
                   f"MOTONEURON POOLS {sd.upper()}")
        pools[sd] = {}
        pools[sd].update(layout_pools(cv, cfg["x0"] + 0.2, 3.5, ROW1,
                                      cap=4, gap=0.24))
        pools[sd].update(layout_pools(cv, cfg["x0"] + 0.2, 2.1, ROW2,
                                      cap=4, gap=0.24))
        for ph, node in pf.items():
            for g, w in W[ph].items():
                if g in pools[sd] and w >= W_MIN:
                    row = 0 if g in ROW1 else 1
                    syn(cv, node, pools[sd][g]["top"], True,
                        lw=0.7 + min(w, 0.25) * 6, shrink_b=0.05,
                        rad=0.0 if row == 0 else -0.12,
                        color=W_E if ph[0] == "E" else W_F)
        # representative reflex wiring
        syn(cv, sens[sd]["ia"], pools[sd]["knee_ext"]["top"], True,
            color=W_EXC, lw=0.9, rad=-0.15)
        syn(cv, sens[sd]["ia"], pools[sd]["knee_flex"]["top"], False,
            color=W_INH, lw=0.9, rad=0.12)
        syn(cv, sens[sd]["ib"], pools[sd]["ankle_pf"]["top"], False,
            color=W_INH, lw=0.9, rad=-0.12)
        syn(cv, ib[sd], pools[sd]["hip_ext"]["top"], True, color=W_EXC,
            lw=0.9, rad=0.15)
        # balance bus: boxes feed the outer lane; lane feeds pools
        syn(cv, b1, (cfg["lane"], 12.55), True, color=W_DESC, dashed=True,
            shrink_a=0.0, rad=0.0)
        syn(cv, b2, (cfg["lane"], 12.4), True, color=W_DESC, dashed=True,
            shrink_a=0.0, rad=0.0)
        syn(cv, (cfg["lane"], 12.3), pools[sd]["ankle_pf"]["top"], True,
            color=W_DESC, dashed=True, rad=-0.3, shrink_a=0.0)
        syn(cv, (cfg["lane"], 4.3), pools[sd]["hip_abd"]["right"], True,
            color=W_DESC, dashed=True, rad=0.1, shrink_a=0.0)
        syn(cv, (cfg["lane"], 4.3), pools[sd]["trunk_ext"]["right"], True,
            color=W_DESC, dashed=True, rad=0.1, shrink_a=0.0)

    # commissural wires through the corridor
    syn(cv, rg["r"][1], rg["l"][1], False, color=W_INH, lw=2.0, rad=-0.12)
    syn(cv, rg["r"][0], rg["l"][0], False, color=W_INH, lw=1.1, rad=0.08)

    # V-class annotations (--vclasses): tags adjacent to the annotated
    # elements + the class-by-class mapping strip at the top
    if vc:
        vclass_strip(cv)
        tag(cv, MID, 8.78, "V0$_D$/dI6 analog", fs=6.4, box=True)
        tag(cv, MID, 8.30, "no V3 analog:\nall cross-side inhibitory",
            fs=5.9, box=True)
        tag(cv, 16.55, 10.52, "V1/V2b analog", fs=6.6, box=True)
        tag(cv, 16.55, 8.55, "V2$_a$ analog\nRG$\\to$PF · PF$\\to$MN exc.",
            fs=6.2, box=True)
        tag(cv, 11.85, 6.45, "V1/V2b analog\n(PF recip. inhib.)", fs=6.2,
            box=True)

    # legend / key / source
    legend_row(cv, 0.55)
    wire_color_key(cv, 9.0, 0.02)
    tag(cv, 15.3, 1.1, "pool abbreviations & member muscles:\n"
        "see the circuit_core key", fs=6.0, box=True)
    source_note(cv, 0.55, 13.15, f"source: {t['note']}")
    return cv.save("circuit_full", fmts)


# ------------------------------------------------------- fig: deng-style
def map_circle(cv, x, y, text, sub=None, r=0.52):
    """Red conversion-map circle: continuous quantity <-> neuron signal
    (Deng/Nourse Fig 6 red circles; here with the ACTUAL formula)."""
    cv.ax.add_patch(plt.Circle((x, y), r, fc="#fde3e0", ec=OI_VERM,
                               lw=1.4, zorder=3))
    cv.ax.text(x, y + (0.13 if sub else 0), text, ha="center", va="center",
               fontsize=6.2, zorder=4, family="monospace", linespacing=1.2)
    if sub:
        cv.ax.text(x, y - 0.24, sub, ha="center", va="center",
                   fontsize=5.4, zorder=4, style="italic", color="0.35")
    return (x, y, r)


def ghost(cv, x, y, label, sub=None, r=0.30):
    """Dashed ghost cell: exists in Deng et al./Nourse 2023, lumped or
    absent in ours."""
    cv.ax.add_patch(plt.Circle((x, y), r, fc="white", ec="0.55",
                               lw=1.1, ls=(0, (3, 2)), zorder=3))
    cv.ax.text(x, y, label, ha="center", va="center", fontsize=6.0,
               color="0.45", zorder=4)
    if sub:
        cv.ax.text(x, y - r - 0.15, sub, ha="center", va="center",
                   fontsize=5.4, color="0.45", style="italic", zorder=4)
    return (x, y, r)


def _edge_groups(net):
    """(src_kind, dst_kind, sign) -> count, from the COMPILED network
    object (structure-driven figure contract: every group must be drawn)."""
    from collections import defaultdict
    names = [p["name"] for p in net.net.populations]

    def kind(nm):
        for pre, k in (("PRESET_E", "PRESET_E"), ("PRESET_F", "PRESET_F"),
                       ("PREA", "PREA"),
                       ("CIN_F", "CIN"), ("CIN_E", "CIN"),
                       ("PF_IN_E", "PF-IN"), ("PF_IN_F", "PF-IN"),
                       ("KINH", "KINH"), ("IBEXC", "IBEXC"),
                       ("LBIN", "LBIN"), ("HEEL", "HEEL"), ("TOE", "TOE"),
                       ("AFF_E", "AFF"), ("AFF_F", "AFF"),
                       ("ADAP", "ADAP"), ("RG_E", "RG-E"), ("RG_F", "RG-F"),
                       ("InE", "RG-IN"), ("InF", "RG-IN"),
                       ("IaIN", "IaIN"), ("RC_", "RC"),
                       ("PF_", "PF"), ("PFA", "PFA"), ("MN_", "MN"),
                       ("Ia_", "Ia"), ("II_", "II"), ("Ib_", "Ib")):
            if nm.startswith(pre):
                return k
        return nm

    g = defaultdict(int)
    for c in net.net.connections:
        g[(kind(names[c["source"]]), kind(names[c["destination"]]),
           "exc" if c["params"].get("reversal_potential", 0) > -1e-6
           else "inh")] += 1
    return dict(g)


DENG_A6 = (
    "Deng et al./Nourse 2023 Table A6 (rat hindlimb SNS) - reference values:\n"
    "  HC->IN 2.75 uS exc | IN->HC 2.75 uS inh   (RG+PF mutual inhibition IS LAMINATED\n"
    "  through dedicated INs; no direct HC<->HC synapses, no mutual excitation)\n"
    "  RGHC->PFHC 0.10 uS exc | PF->MN 1.5-4.9 uS exc (per joint)\n"
    "  PF->Ia 0.50 uS exc (phase gate) | IaIN->MN 2.0 uS inh | IaIN<->IaIN 0.5 inh\n"
    "  MN->RC 0.50 uS exc | RC->MN 0.50 uS inh | RC<->RC 0.5 inh | RC->IaIN 0.5 inh\n"
    "  IbIN->MN 0.59 uS EXCITATORY (positive force feedback) | PF->Ib 2.0 uS shunt\n"
    "  MN->muscle: act = 1/(1+exp(s(x0-V)))+y0, s=0.153, x0=-70 mV (Fig 6B)\n"
    "  feedback: MuJoCo muscle TENSION formatted as Ia and Ib input (no II)"
)


def _tuned_gains() -> dict:
    """Representative gains for the Deng schematic.

    Use the current curriculum winners only when their score is above the
    -100 no-countable-cycles sentinel. Otherwise use the successful stage-1
    rhythm gains plus explicitly labelled, nonzero representative feedback
    gains so every conditional pathway can be audited in one diagram.
    """
    g = dict(phase_reset_e=0.0, phase_reset_f=0.0, f1_kneext_inh=0.59,
             f1_anklepf_inh=0.6, renshaw=0.5, ia_in=0.6, heel_rge=0.6,
             toe_rge=0.4, ib_rge=0.6, desc_e=1.7, desc_f=1.4, rg_to_pf=2.4,
             # The directional V3/C1 commissural relays provide the
             # bilateral coordination shown in the literature figure; do
             # not add a separate ipsilateral RG-E<->RG-F excitatory pair.
             rg_weak_exc=0.0, ib_e_central=0.5, ia_f_central=0.5,
             ii_f_central=0.3, ii_e_central=0.3, ia_f_contra_f=0.4,
             v3_to_ibexc=0.4, c1_gain=0.6, v3_gain=0.25,
             src="representative defaults (NaP architecture)")
    try:
        s1 = json.loads(
            (HERE / "curriculum_stage1.json").read_text("utf-8"))
        if float(s1.get("score", -100.0)) > -100.0:
            for k in ("desc_e", "desc_f", "rg_to_pf"):
                if k in s1.get("params", {}):
                    g[k] = float(s1["params"][k])
            g["src"] = ("curriculum stage-1 winner + representative "
                        "nonzero feedback gains")
    except FileNotFoundError:
        pass
    try:
        s3 = json.loads(
            (HERE / "curriculum_stage3.json").read_text("utf-8"))
        if float(s3.get("score", -100.0)) > -100.0:
            for k in ("phase_reset_e", "phase_reset_f", "heel_rge",
                      "toe_rge", "ib_rge", "ia_in", "desc_e", "desc_f",
                      "rg_to_pf", "ib_e_central", "ia_f_central",
                      "ii_f_central", "ii_e_central", "ia_f_contra_f",
                      "v3_to_ibexc", "c1_gain", "v3_gain"):
                if k in s3.get("params", {}):
                    g[k] = float(s3["params"][k])
            g["src"] = (f"curriculum_stage3.json winner (score "
                        f"{s3.get('score')})")
        else:
            g["src"] += (f"; failed stage-3 score {s3.get('score')} ignored")
    except FileNotFoundError:
        pass
    return g


def fig_deng(t, fmts):
    """Deng/Nourse-style layered schematic, STRUCTURE-DRIVEN: builds the
    representative network (knee extensor + flexor column per side, with
    the LAMINATED architecture live: InE/InF RG relays, PF_IN_E/F PF
    relays, Renshaw + IaIN + heel/toe/LBIN stance feedback at the tuned
    gains), reads every connection from the compiled SNS object, and
    asserts each edge group is drawn."""
    import build_network as bn
    import params as _P
    tg = _tuned_gains()
    gain_keys = ("phase_reset_e", "phase_reset_f", "f1_kneext_inh",
                 "renshaw", "ia_in", "heel_rge", "toe_rge", "ib_rge",
                 "rg_weak_exc", "ib_e_central", "ia_f_central",
                 "ii_f_central", "ii_e_central", "ia_f_contra_f",
                 "v3_to_ibexc", "c1_gain", "v3_gain")
    saved_gains = {k: _P.G[k] for k in gain_keys}
    for k in gain_keys:
        _P.G[k] = tg[k]
    ACTS = ["vas_lat_r", "semimem_r", "vas_lat_l", "semimem_l"]
    try:
        net = bn.build(ACTS, interleg=True)
    finally:
        _P.G.update(saved_gains)
    print(f"deng representative gains: {tg['src']}")
    Gc = _edge_groups(net)
    drawn = set()

    def edge(sk, dk, sign):
        drawn.add((sk, dk, sign))

    cv = Canvas(13.8, 18.6)

    # ---------------- band 0: supraspinal ----------------
    layer_band(cv, 0.5, 13.3, 16.6, 18.2, OI_SKY, "SUPRASPINAL (inputs)")
    drv = input_box(cv, 2.6, 17.4, 2.0, 0.6, "DRIVE", "MLR speed cmd")
    pos = input_box(cv, 6.0, 17.4, 2.0, 0.6, "POSTURE", "tonic")
    tag(cv, 2.6, 16.85, "BAL family (COM PD, IMU trunk, lateral ->\n"
        "ankle/hip/abd/trunk MNs in the full net) omitted here",
        fs=5.6, box=True)

    # ---------------- band 1: rhythm generator ----------------
    layer_band(cv, 0.5, 13.3, 12.7, 16.45, OI_ORANGE,
               "RHYTHM GENERATOR (per side; right shown)")
    rge = neuron(cv, 3.6, 15.35, "RG-E", "NaP", fc=_tint(W_E), r=R_BIG)
    rgf = neuron(cv, 3.6, 13.55, "RG-F", "NaP", fc=_tint(W_F), r=R_BIG)
    # IN-laminated mutual inhibition (Shevtsova/Deng A6): RG-E excites
    # InE, InE inhibits RG-F; RG-F excites InF, InF inhibits RG-E.
    # NO direct inhibitory RG<->RG synapses (removed 2026-09-15).
    ine = neuron(cv, 1.75, 15.35, "InE", fc=_tint(OI_PURPLE), r=R_MED)
    inf = neuron(cv, 1.75, 13.55, "InF", fc=_tint(OI_PURPLE), r=R_MED)
    syn(cv, rge, ine, True, color=W_E, lw=1.4)
    syn(cv, ine, rgf, False, color=W_INH, lw=1.6, rad=0.25)
    syn(cv, rgf, inf, True, color=W_F, lw=1.4)
    syn(cv, inf, rge, False, color=W_INH, lw=1.6, rad=0.25)
    edge("RG-E", "RG-IN", "exc"); edge("RG-IN", "RG-F", "inh")
    edge("RG-F", "RG-IN", "exc"); edge("RG-IN", "RG-E", "inh")
    tag(cv, 5.15, 14.45, "IN-laminated mutual inhibition\n(Shevtsova/Deng "
        "A6: RG-E->InE->RG-F,\nRG-F->InF->RG-E; no direct\nHC<->HC "
        "synapses)", fs=6.2, box=True)
    # weak mutual excitation G_W between the half-centers (Deng 2022)
    if tg["rg_weak_exc"] > 0.0:
        syn(cv, rge, rgf, True, color=W_E, lw=1.0, rad=-0.45)
        syn(cv, rgf, rge, True, color=W_F, lw=1.0, rad=-0.45)
        edge("RG-E", "RG-F", "exc"); edge("RG-F", "RG-E", "exc")
    tag(cv, 7.3, 15.35, "RG = persistent-Na\ncond. bursters (Deng 2022,\n"
        "Shinohara 2025, Rybak 2024):\nintrinsic h-gate burst\ntermination, "
        "FIXED tau_h\n(= period knob) + weak\nmutual excitation G_W\n"
        "(escape mode; speed\nneuromodulation).\nADAP retired 2026-09-16",
        fs=5.2, box=True)
    syn(cv, drv, rge, True, color=W_DESC, rad=0.1,
        label=f"{tg['desc_e']:.1f}", lfs=5.4)
    syn(cv, drv, rgf, True, color=W_DESC, rad=-0.12,
        label=f"{tg['desc_f']:.1f}", lfs=5.4)
    syn(cv, pos, rge, True, color=W_DESC, rad=-0.2,
        label=f"{params.G['posture_to_rg_e']:.1f}", lfs=5.4)
    edge("DRIVE", "RG-E", "exc"); edge("DRIVE", "RG-F", "exc")
    edge("POSTURE", "RG-E", "exc")
    # v11 mechanosensory stance feedback (audit P1a; conditional on the
    # tuned gains): heel contact -> RG-E exc + RG-F inh (S2W trigger),
    # loaded toe -> RG-E exc (late-stance prolongation) + both ride the
    # extensor central pathway (PF_E + InE) when ib_e_central > 0
    if tg["heel_rge"] > 0.0 or tg["toe_rge"] > 0.0:
        hel = neuron(cv, 1.75, 12.95, "HEEL", r=R_SMALL, fs=5.4)
        toe = neuron(cv, 3.3, 12.95, "TOE", r=R_SMALL, fs=5.4)
        if tg["heel_rge"] > 0.0:
            syn(cv, hel, rge, True, color=W_E, lw=1.1, rad=0.18)
            syn(cv, hel, rgf, False, color=W_INH, lw=1.0, rad=0.12)
            edge("HEEL", "RG-E", "exc"); edge("HEEL", "RG-F", "inh")
        if tg["toe_rge"] > 0.0:
            syn(cv, toe, rge, True, color=W_E, lw=1.1, rad=0.15)
            edge("TOE", "RG-E", "exc")
        tag(cv, 5.35, 13.0, "v11 contact mechanosensors\n(S2W trigger + "
            "late-stance\nprolongation; same extensor\ncentral pathway as "
            "Ib)", fs=5.2, box=True)
    # commissural stubs (right margin), laminated through CIN INs
    # (Shinohara 2025 wiring): RG-F -> c1 IN -> INHIBIT -> contra RG-F;
    # RG-E -> V3 IN -> EXCITE -> contra RG-E (+ contra extensor MN groups
    # via IBEXC, audit P2a, when v3_to_ibexc > 0)
    cin_f = neuron(cv, 12.3, 13.55, "c1", fc=_tint(OI_PURPLE),
                   r=R_SMALL, fs=5.4)
    cin_e = neuron(cv, 12.3, 15.35, "V3", fc=_tint(OI_GREEN),
                   r=R_SMALL, fs=5.4)
    syn(cv, rgf, cin_f, True, color=W_F, lw=1.6)
    syn(cv, cin_f, (13.15, 13.55), False, color=W_INH, lw=2.0, shrink_b=0.0)
    syn(cv, rge, cin_e, True, color=W_E, lw=1.1)
    syn(cv, cin_e, (13.15, 15.35), True, color=W_EXC, lw=1.4, shrink_b=0.0)
    tag(cv, 12.35, 12.95, f"c1_r -> INHIBIT RG-F_l\n(g="
        f"{tg['c1_gain'] * params.G['rg_mutual_inh']:.1f}); V3_r ->\n"
        f"EXCITE InE_l (g={tg['v3_gain'] * params.G['rg_mutual_inh']:.1f})"
        "\n+ contra extensor MNs "
        "(P2a)", fs=5.4, box=True)
    edge("RG-F", "CIN", "exc"); edge("CIN", "RG-F", "inh")
    edge("RG-E", "CIN", "exc"); edge("CIN", "RG-IN", "exc")

    # ---------------- band 2: pattern formation ----------------
    layer_band(cv, 0.5, 13.3, 9.5, 12.55, OI_SKY,
               "PATTERN FORMATION (4 phase-window cells, laminated "
               "cross-reciprocal via PF_IN)")
    pfx = {"E1": 2.3, "E2": 4.4, "F1": 6.9, "F2": 9.0}
    pfs = {}
    for nm, x in pfx.items():
        c = W_E if nm[0] == "E" else W_F
        pfs[nm] = neuron(cv, x, 11.7, f"PF-{nm}", fc=_tint(c))
    syn(cv, rge, pfs["E1"], True, color=W_E, rad=-0.06)
    syn(cv, rge, pfs["E2"], True, color=W_E, rad=0.1)
    syn(cv, rgf, pfs["F1"], True, color=W_F, rad=-0.1)
    syn(cv, rgf, pfs["F2"], True, color=W_F, rad=0.06)
    edge("RG-E", "PF", "exc"); edge("RG-F", "PF", "exc")
    # foot mechanosensors ride the extensor central pathway (PF_E + InE)
    # when ib_e_central > 0 (Ben 2026-09-16: same route as extensor Ib)
    if (tg["heel_rge"] > 0.0 or tg["toe_rge"] > 0.0) \
            and tg["ib_e_central"] > 0.0:
        syn(cv, hel, pfs["E1"], True, color=W_E, lw=0.7, rad=0.15)
        syn(cv, hel, ine, True, color=W_E, lw=0.7, rad=0.1)
        syn(cv, toe, pfs["E2"], True, color=W_E, lw=0.7, rad=0.12)
        syn(cv, toe, ine, True, color=W_E, lw=0.7, rad=-0.1)
        edge("HEEL", "PF", "exc"); edge("HEEL", "RG-IN", "exc")
        edge("TOE", "PF", "exc"); edge("TOE", "RG-IN", "exc")
    # (DRIVE->PF weak tonic edge REMOVED 2026-09-16, audit #23:
    # EXTRA-NOTSUPPORTED - PF is driven by RG only, per Deng/Shevtsova)
    # direct PF<->PF inhibition REMOVED (laminated 2026-09-15)
    pf_in_e = neuron(cv, 5.65, 12.42, "PF_IN_E", fc=_tint(OI_PURPLE),
                     r=R_SMALL, fs=5.6)
    pf_in_f = neuron(cv, 5.65, 10.85, "PF_IN_F", fc=_tint(OI_PURPLE),
                     r=R_SMALL, fs=5.6)
    syn(cv, pfs["E1"], pf_in_e, True, color=W_E, lw=0.9, rad=0.1)
    syn(cv, pfs["E2"], pf_in_e, True, color=W_E, lw=0.9, rad=-0.1)
    syn(cv, pfs["F1"], pf_in_f, True, color=W_F, lw=0.9, rad=0.1)
    syn(cv, pfs["F2"], pf_in_f, True, color=W_F, lw=0.9, rad=-0.1)
    syn(cv, pf_in_e, pfs["F1"], False, color=W_INH, lw=0.9, rad=-0.22)
    syn(cv, pf_in_e, pfs["F2"], False, color=W_INH, lw=0.9, rad=0.18)
    syn(cv, pf_in_f, pfs["E1"], False, color=W_INH, lw=0.9, rad=0.22)
    syn(cv, pf_in_f, pfs["E2"], False, color=W_INH, lw=0.9, rad=-0.18)
    edge("PF", "PF-IN", "exc"); edge("PF-IN", "PF", "inh")
    tag(cv, 4.45, 9.78, "IN-laminated PF cross-reciprocal\n(Shevtsova/Deng "
        "A6; no direct PF<->PF)", fs=5.8, box=True)
    kinh = neuron(cv, 11.6, 10.5, "KINH\ncond.", r=R_SMALL, fs=5.2)
    syn(cv, pfs["F1"], kinh, True, lw=1.0, color=W_F, rad=-0.15)
    edge("PF", "KINH", "exc")
    tag(cv, 11.6, 9.7, "v6 swing quad-\nsuppression IN", fs=5.6, box=True)

    # ---------------- band 3: motor circuit (knee representative) ----------------
    layer_band(cv, 0.5, 13.3, 4.9, 9.35, OI_GREEN,
               "MOTOR CIRCUIT - knee column (representative; x46 muscle"
               " columns per side)")
    mne = neuron(cv, 3.0, 8.5, "MN\nknee-ext", fc="white", r=0.42, fs=6.6)
    mnf = neuron(cv, 9.6, 8.5, "MN\nknee-flx", fc="white", r=0.42, fs=6.6)
    syn(cv, pfs["E1"], mne, True, color=W_E, rad=0.15, lw=1.2)
    syn(cv, pfs["E2"], mne, True, color=W_E, rad=-0.1, lw=1.2)
    syn(cv, pfs["F1"], mnf, True, color=W_F, rad=0.12, lw=1.2)
    syn(cv, pfs["F2"], mnf, True, color=W_F, rad=-0.08, lw=0.8)
    edge("PF", "MN", "exc")
    tag(cv, 4.9, 9.0, "PF->MN x W_PF_MN(group)\n(F2 = late-swing prep)",
        fs=5.8, box=True)
    syn(cv, pos, (1.0, 17.4), True, color=W_DESC, dashed=True,
        shrink_b=0.0, lw=0.8)
    lane(cv, 1.0, 17.4, 8.5, color=W_DESC)
    syn(cv, (1.0, 8.6), mne, True, color=W_DESC, dashed=True, shrink_a=0.0)
    syn(cv, (1.0, 8.6), mnf, True, color=W_DESC, dashed=True, shrink_a=0.0,
        rad=-0.1)
    edge("POSTURE", "MN", "exc")
    tag(cv, 1.0, 8.95, "POSTURE (+ per-MN\nPOST_i standing bias)", fs=5.6,
        box=True)
    syn(cv, kinh, mne, False, color=W_INH, lw=1.4, rad=0.18)
    edge("KINH", "MN", "inh")
    # afferents (extensor column fully labeled; flexor mirrored)
    aff_y = 6.55
    ia_e = neuron(cv, 2.0, aff_y, "Ia", r=R_SENS, fc=_tint(OI_ORANGE))
    ii_e = neuron(cv, 3.0, aff_y, "II", r=R_SENS, fc=_tint(OI_ORANGE))
    ib_e = neuron(cv, 4.0, aff_y, "Ib", r=R_SENS, fc=_tint(OI_ORANGE))
    ia_f = neuron(cv, 8.6, aff_y, "Ia", r=R_SENS, fc=_tint(OI_ORANGE))
    ii_f = neuron(cv, 9.6, aff_y, "II", r=R_SENS, fc=_tint(OI_ORANGE))
    ib_f = neuron(cv, 10.6, aff_y, "Ib", r=R_SENS, fc=_tint(OI_ORANGE))
    tag(cv, 2.55, 7.1, "per-muscle afferent\nencoders (x92)", fs=5.6,
        box=True)
    syn(cv, ia_e, mne, True, color=W_EXC, rad=-0.1)
    syn(cv, ii_e, mne, True, color=W_EXC, rad=0.0)
    syn(cv, ib_e, mne, False, color=W_INH, rad=0.12)
    edge("Ia", "MN", "exc"); edge("II", "MN", "exc"); edge("Ib", "MN", "inh")
    syn(cv, ia_f, mnf, True, color=W_EXC, rad=-0.1)
    syn(cv, ii_f, mnf, True, color=W_EXC)
    syn(cv, ib_f, mnf, False, color=W_INH, rad=0.12)
    tag(cv, 2.0, 5.9, "Ia homonymous exc\nII length exc\nIb autogenic inh",
        fs=5.6, box=True)
    # Renshaw recurrent inhibition (Deng A6; IMPLEMENTED 2026-09-14):
    # per-pool RC, MN->RC exc, RC->MN inh, mutual RC<->RC between
    # distinct pools (no autapse)
    if tg["renshaw"] > 0.0:
        rc_e = neuron(cv, 4.55, 8.5, "RC", fc="#d9d9d9", r=R_SMALL,
                      fs=6.0)
        rc_f = neuron(cv, 8.05, 8.5, "RC", fc="#d9d9d9", r=R_SMALL,
                      fs=6.0)
        syn(cv, mne, rc_e, True, color=W_EXC, lw=1.0)
        syn(cv, rc_e, mne, False, color=W_INH, lw=1.0, rad=0.35)
        syn(cv, mnf, rc_f, True, color=W_EXC, lw=1.0)
        syn(cv, rc_f, mnf, False, color=W_INH, lw=1.0, rad=-0.35)
        syn(cv, rc_e, rc_f, False, color=W_INH, lw=0.9, rad=-0.3)
        syn(cv, rc_f, rc_e, False, color=W_INH, lw=0.9, rad=-0.3)
        edge("MN", "RC", "exc"); edge("RC", "MN", "inh")
        edge("RC", "RC", "inh")
        tag(cv, 6.3, 8.62, "Renshaw: MN->RC 1.0, RC->MN "
            f"g={tg['renshaw']:.2f}\nmutual RC<->RC (distinct pools)",
            fs=5.4, box=True)
    # Ia reciprocal inhibition of the antagonist: via per-pool IaIN
    # (PF-F1 phase gate + RC->IaIN disinhibition) when ia_in > 0;
    # the direct edge exists only at ia_in == 0 (v10 behavior)
    if tg["ia_in"] > 0.0:
        iain_e = neuron(cv, 5.9, 6.55, "IaIN", fc=_tint(OI_PURPLE),
                        r=0.26, fs=5.8)
        iain_f = neuron(cv, 7.5, 6.55, "IaIN", fc=_tint(OI_PURPLE),
                        r=0.26, fs=5.8)
        syn(cv, pfs["F1"], iain_e, True, color=W_F, lw=0.9,
            rad=0.13)
        syn(cv, pfs["F1"], iain_f, True, color=W_F, lw=0.9,
            rad=-0.13)
        edge("PF", "IaIN", "exc")
        tag(cv, 6.7, 7.28, "PF-F1 phase gate (0.5)", fs=5.2, box=True)
        # afferent leg (fixed 2026-09-16): Ia excites its IaIN
        syn(cv, ia_e, iain_e, True, color=W_EXC, lw=1.0, rad=0.2)
        syn(cv, ia_f, iain_f, True, color=W_EXC, lw=1.0, rad=0.2)
        edge("Ia", "IaIN", "exc")
        syn(cv, iain_e, mnf, False, color=W_INH, lw=1.6, rad=-0.15)
        syn(cv, iain_f, mne, False, color=W_INH, lw=1.6, rad=-0.15)
        edge("IaIN", "MN", "inh")
        if tg["renshaw"] > 0.0:
            syn(cv, rc_e, iain_e, False, color=W_INH, lw=1.0, rad=0.1)
            syn(cv, rc_f, iain_f, False, color=W_INH, lw=1.0, rad=0.1)
            edge("RC", "IaIN", "inh")
        tag(cv, 6.7, 5.92, "Ia reciprocal inhibition via IaIN: Ia->IaIN\n"
            f"(g={params.G['ia_to_mn']:.2f}), PF-F1-gated; IaIN->MN "
            f"g={params.G['ia_to_antagonist']:.2f};\nRC->IaIN = recurrent "
            "disinhibition", fs=5.4, box=True)
    else:
        # reciprocal: Ia_ext -> MN_flex (direct/lumped)
        syn(cv, ia_e, mnf, False, color=W_INH, lw=1.6, rad=-0.22)
        syn(cv, ia_f, mne, False, color=W_INH, lw=1.6, rad=-0.22)
        edge("Ia", "MN", "inh")
        tag(cv, 6.3, 7.6, "Ia reciprocal inhibition (DIRECT -\n"
            "ia_in == 0 configuration)", fs=5.8, box=True)
    # Ib load sharing
    ibx = neuron(cv, 5.2, 7.6, "IB-EXC", r=R_SMALL, fs=6.0,
                 fc=_tint(OI_GREEN))
    syn(cv, ib_e, ibx, True, color=W_EXC, rad=-0.1, lw=1.0)
    edge("Ib", "IBEXC", "exc")
    syn(cv, ibx, mne, True, color=W_EXC, rad=0.15, lw=1.0)
    edge("IBEXC", "MN", "exc")
    # RG-E stance gate: route down the far-left lane (below the POSTURE
    # lane) so the long wire does not slash across the PF band
    lane(cv, 0.72, 15.1, 8.15, color=W_E)
    syn(cv, rge, (0.72, 15.1), True, color=W_E, dashed=True,
        shrink_b=0.0, lw=0.8)
    syn(cv, (0.72, 8.3), ibx, True, color=W_E, dashed=True, shrink_a=0.0,
        rad=0.0, lw=0.8)
    tag(cv, 5.2, 8.25, "RG-E stance gate g=1.0\n(loaded extensor groups "
        "only; x5 groups)", fs=5.6, box=True)
    edge("RG-E", "IBEXC", "exc")
    tag(cv, 4.95, 7.0, "stance load sharing\n(reflex reversal)", fs=5.6,
        box=True)
    # v11 P1a: stance-Ib group IN living in the RG layer (Dominguez 2020):
    # IBEXC -> LBIN -> RG-E (load prolongs stance; the duty lever)
    if tg["ib_rge"] > 0.0:
        lbin = neuron(cv, 6.35, 7.6, "LBIN", fc=_tint(OI_GREEN), r=R_SMALL,
                      fs=5.8)
        syn(cv, ibx, lbin, True, color=W_EXC, lw=0.9, rad=-0.1)
        edge("IBEXC", "LBIN", "exc")
        lane(cv, 6.35, 15.35, 7.78, color=W_E)
        syn(cv, (6.35, 15.35), rge, True, color=W_E, dashed=True,
            shrink_a=0.0, shrink_b=0.0, lw=0.9)
        syn(cv, lbin, (6.35, 7.78), True, color=W_E, dashed=True,
            shrink_a=0.0, shrink_b=0.0, lw=0.9)
        edge("LBIN", "RG-E", "exc")
        tag(cv, 7.65, 7.9, "LBIN stance-Ib IN\n(RG layer, Dominguez 2020;\n"
            "load prolongs stance)", fs=5.2, box=True)

    # ---- per-muscle afferent -> CENTRAL feedback (Deng 2022 figure /
    # Shinohara 2025 / Rybak 2025 SF-E1/SF-E2, Ben 2026-09-16): thin
    # arrows from the afferent row up into the pattern/rhythm layers
    if tg["ib_e_central"] > 0.0:
        syn(cv, ib_e, pfs["E1"], True, color=W_E, lw=0.8, rad=0.2)
        edge("Ib", "PF", "exc")
        lane(cv, 0.95, 15.35, 6.7, color=W_E)
        syn(cv, ib_e, (0.95, 6.7), True, color=W_E, dashed=True, lw=0.7,
            shrink_a=0.0, shrink_b=0.0)
        syn(cv, (0.95, 15.35), rge, True, color=W_E, dashed=True, lw=0.7,
            shrink_a=0.0, shrink_b=0.0)
        edge("Ib", "RG-E", "exc")
        syn(cv, ib_e, ine, True, color=W_E, lw=0.7, rad=0.15)
        edge("Ib", "RG-IN", "exc")
    if tg["ia_f_central"] > 0.0:
        syn(cv, ia_f, pfs["F1"], True, color=W_F, lw=0.8, rad=0.2)
        edge("Ia", "PF", "exc")
        syn(cv, ia_f, rgf, True, color=W_F, lw=0.7, rad=-0.1)
        edge("Ia", "RG-F", "exc")
        syn(cv, ia_f, inf, True, color=W_F, lw=0.7, rad=0.12)
        edge("Ia", "RG-IN", "exc")
    if tg["ii_f_central"] > 0.0:
        syn(cv, ii_f, pfs["F1"], True, color=W_F, lw=0.7, rad=0.25)
        edge("II", "PF", "exc")
        syn(cv, ii_f, rgf, True, color=W_F, lw=0.7, rad=0.1)
        edge("II", "RG-F", "exc")
    if tg["ii_e_central"] > 0.0:
        syn(cv, ii_e, rge, True, color=W_E, lw=0.7, rad=-0.15)
        edge("II", "RG-E", "exc")
        syn(cv, ii_e, ine, True, color=W_E, lw=0.7, rad=-0.1)
        edge("II", "RG-IN", "exc")
    if tg["ia_f_contra_f"] > 0.0:
        syn(cv, ia_f, (13.15, 13.55), False, color=W_INH, lw=1.0,
            shrink_b=0.0, rad=-0.1)
        edge("Ia", "RG-F", "inh")
    tag(cv, 6.9, 5.55, "per-muscle afferent -> central (same-group exc):\n"
        "extensor Ib(+II) -> PF_E / RG_E / InE (Rybak SF-E2,\nforce "
        "feedback prolongs stance); flexor Ia+II -> PF_F / RG_F /\nInF "
        "(SF-E1, length feedback triggers swing); flexor Ia\nalso inhibits "
        "the CONTRALATERAL RG-F (SF-E1, interleg)", fs=5.2, box=True)

    # All other hip/ankle/trunk muscle columns repeat this complete motif.
    # Do not draw partial MN glyphs without their plant/reflex connections.
    # V3 -> contralateral extensor MN groups (audit P2a)
    if tg["v3_to_ibexc"] > 0.0:
        syn(cv, cin_e, ibx, True, color=W_EXC, lw=0.8, rad=-0.15,
            dashed=True)
        edge("CIN", "IBEXC", "exc")
    tag(cv, 12.0, 8.55, "same complete MN--muscle--afferent motif\n"
        "repeats for the other 44 hip, ankle,\nfrontal-plane, and trunk "
        "muscle pools", fs=5.2, box=True)

    # ---------------- band 4: neuromuscular interface ----------------
    layer_band(cv, 0.5, 13.3, 1.95, 4.75, OI_VERM,
               "NEUROMUSCULAR INTERFACE (conversion maps: red circles)")
    m1 = muscle(cv, 3.0, 2.65, w=0.9, h=0.3)
    m2 = muscle(cv, 9.6, 2.65, w=0.9, h=0.3)
    cv.ax.text(3.0, 2.16, "knee-ext muscle\n(MuJoCo Hill F-L-V,\nrigid tendon,\n"
               "Fmax = OpenSim)", ha="center", va="center", fontsize=5.6,
               zorder=4)
    cv.ax.text(9.6, 2.16, "knee-flx muscle\n(same MuJoCo model)",
               ha="center", va="center",
               fontsize=5.6, zorder=4)
    plant_drawn = set()

    def plant_syn(p, q, src, dst, **kwargs):
        syn(cv, p, q, True, **kwargs)
        plant_drawn.add((src, dst))

    # MN -> activation map -> muscle, shown completely for BOTH antagonists.
    act_e = map_circle(cv, 5.15, 3.35, "a=clip(V/5mV,\n    0,1)",
                       "MN -> activation", r=0.45)
    act_f = map_circle(cv, 11.85, 3.35, "a=clip(V/5mV,\n    0,1)",
                       "MN -> activation", r=0.45)
    plant_syn(mne, act_e, "MN-E", "ACT-E", color="0.35", lw=1.2,
              rad=-0.12)
    plant_syn(act_e, m1, "ACT-E", "MUSCLE-E", color="0.35", lw=1.2,
              rad=-0.08)
    plant_syn(mnf, act_f, "MN-F", "ACT-F", color="0.35", lw=1.2,
              rad=0.12)
    plant_syn(act_f, m2, "ACT-F", "MUSCLE-F", color="0.35", lw=1.2,
              rad=0.08)

    # Muscle Ldot/L/F -> Ia/II/Ib maps, complete on both antagonist sides.
    enc_e = [map_circle(cv, x, 4.15, txt, sub, r=0.40)
             for x, txt, sub in ((1.25, "clip(Ldot/0.6)", "-> Ia"),
                                 (2.35, "(L-Lmid)/Lhalf", "-> II"),
                                 (3.45, "clip(F/Fmax)", "-> Ib"))]
    enc_f = [map_circle(cv, x, 4.15, txt, sub, r=0.40)
             for x, txt, sub in ((7.35, "clip(Ldot/0.6)", "-> Ia"),
                                 (8.45, "(L-Lmid)/Lhalf", "-> II"),
                                 (9.55, "clip(F/Fmax)", "-> Ib"))]
    for side, mus, maps, affs, col in (
            ("E", m1, enc_e, (ia_e, ii_e, ib_e), W_E),
            ("F", m2, enc_f, (ia_f, ii_f, ib_f), W_F)):
        for kind, mp, aff in zip(("IA", "II", "IB"), maps, affs):
            plant_syn(mus, mp, f"MUSCLE-{side}", f"{kind}-MAP-{side}",
                      color=col, lw=0.9, rad=0.04)
            plant_syn(mp, aff, f"{kind}-MAP-{side}", f"{kind}-{side}",
                      color=col, lw=0.9, rad=-0.04)
    tag(cv, 6.25, 4.45, "pure-signal encoders; gains are phase- and\n"
        "speed-gated presynaptically in runner.py", fs=5.4, box=True)
    tag(cv, 6.25, 3.25, "MuJoCo supplies muscle length, velocity, and force;\n"
        "the complete interface repeats for all 92 muscles", fs=5.4,
        box=True)

    required_plant = {
        ("MN-E", "ACT-E"), ("ACT-E", "MUSCLE-E"),
        ("MN-F", "ACT-F"), ("ACT-F", "MUSCLE-F"),
    }
    for side in ("E", "F"):
        for kind in ("IA", "II", "IB"):
            required_plant.add((f"MUSCLE-{side}", f"{kind}-MAP-{side}"))
            required_plant.add((f"{kind}-MAP-{side}", f"{kind}-{side}"))
    assert plant_drawn == required_plant, (
        "incomplete neuromuscular interface drawing: missing="
        f"{sorted(required_plant - plant_drawn)}, extra="
        f"{sorted(plant_drawn - required_plant)}")

    # ---------------- legend + Deng reference strip ----------------
    ax = cv.ax
    y = 1.55
    ax.add_patch(Circle((0.9, y), 0.16, fc="white", ec="black", lw=1.3))
    ax.text(1.15, y, "neuron", fontsize=7.0, va="center")
    ax.add_patch(FancyArrowPatch((2.2, y), (2.75, y), arrowstyle="-",
                                 lw=1.3, color=W_EXC))
    ax.add_patch(Polygon([(3.05, y + 0.085), (3.05, y - 0.085), (2.82, y)],
                         closed=True, fc="white", ec=W_EXC, lw=1.2))
    ax.text(3.15, y, "excitatory", fontsize=7.0, va="center")
    ax.add_patch(FancyArrowPatch((4.6, y), (5.22, y), arrowstyle="-",
                                 lw=1.3, color=W_INH))
    ax.add_patch(Circle((5.32, y), 0.085, fc="black"))
    ax.text(5.45, y, "inhibitory", fontsize=7.0, va="center")
    ax.add_patch(Circle((6.9, y), 0.2, fc="#fde3e0", ec=OI_VERM, lw=1.2))
    ax.text(7.2, y, "conversion map\n(continuous <-> neural)", fontsize=6.4,
            va="center")
    ax.text(9.9, y, "all drawn cells are IMPLEMENTED in build_network.py\n"
            "(edges read back from the compiled network)", fontsize=6.4,
            va="center", style="italic", color="0.35")
    ax.text(0.6, 0.35, DENG_A6, fontsize=5.6, va="bottom", ha="left",
            family="monospace", color="0.30", linespacing=1.25,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.75",
                      lw=0.7))
    source_note(cv, 0.55, 18.47, f"source: {tg['src']} - edges drawn from "
                 "the COMPILED network object (build_network.py, LAMINATED "
                 "architecture), representative knee column x46/side")

    # ---------------- structure-driven contract ----------------
    compiled = set(Gc)
    missing = compiled - drawn
    extra = drawn - compiled
    assert not missing, f"undrawn edge groups: {sorted(missing)}"
    assert not extra, f"drawn edge groups absent from compiled net: {sorted(extra)}"
    return cv.save("circuit_dengstyle", fmts)


# --------------------------------------------------------- fig: weights
def fig_weights(t, fmts):
    W, WPOST = t["W"], t["WPOST"]
    groups = [g for g in GROUP_ORDER
              if any(W[ph].get(g, 0.0) > 0 or g in WPOST
                     for ph in ("E1", "E2", "F1", "F2"))]
    cols = ["E1", "E2", "F1", "F2", "POST"]
    M = np.zeros((len(groups), len(cols)))
    for i, g in enumerate(groups):
        for j, ph in enumerate(("E1", "E2", "F1", "F2")):
            M[i, j] = W[ph].get(g, 0.0)
        M[i, 4] = WPOST.get(g, 0.0)

    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    fig.subplots_adjust(left=0.17, right=0.98, top=0.97, bottom=0.17)
    im = ax.imshow(M, cmap="Oranges", vmin=0, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, fontsize=8)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([f"{ABBREV[g]} ({len(CNT['prim'][g])})" for g in groups],
                       fontsize=8)
    ax.axvline(3.5, color="black", lw=1.2)
    for i in range(len(groups)):
        for j in range(len(cols)):
            v = M[i, j]
            ax.text(j, i, f"{v:.2f}" if v > 0 else "–", ha="center",
                    va="center", fontsize=7,
                    color="white" if v > 0.55 * M.max() else "0.15")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(groups), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.4)
    ax.tick_params(which="both", length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("PF$\\to$MN weight (pf_to_mn multiplier)", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    ax.set_xlabel("PF phase cells (E1 early stance, E2 push-off, F1 early "
                  "swing, F2 late swing); POST = tonic posture column",
                  fontsize=7, style="italic", color="0.3", labelpad=8)
    fig.text(0.17, 0.035, f"source: {t['note']}; pools: "
             f"{N_MN_TOTAL} MN total, counts live from muscle_map.py",
             fontsize=6.3, color="0.4", style="italic")
    outs = []
    for fmt in fmts:
        p = OUT / f"circuit_weights.{fmt}"
        fig.savefig(p, dpi=300, bbox_inches="tight", pad_inches=0.05)
        outs.append(p.name)
    plt.close(fig)
    return outs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", default="all",
                    choices=["core", "full", "deng", "weights", "all"])
    ap.add_argument("--source", default="best",
                    choices=["params", "fitted", "best"])
    ap.add_argument("--fmt", default="pdf,svg,png")
    ap.add_argument("--vclasses", action="store_true",
                    help="annotate circuit_full with the V-class "
                         "correspondence (Shevtsova 2026 / Rybak 2015)")
    args = ap.parse_args()
    fmts = [f.strip() for f in args.fmt.split(",") if f.strip()]
    OUT.mkdir(exist_ok=True)
    t = effective_tables(args.source)
    print(f"source: {t['note']}")
    print(f"network: {N_TOTAL} neurons ({N_SIDE}/side + {N_SHARED} shared), "
          f"{N_INPUTS} inputs, pools live from muscle_map.py")
    builders = dict(core=fig_core, full=fig_full, deng=fig_deng,
                    weights=fig_weights)
    if args.vclasses:
        builders["full"] = lambda tt, ff: fig_full(tt, ff, vc=True)
    which = list(builders) if args.which == "all" else [args.which]
    for nm in which:
        outs = builders[nm](t, fmts)
        print(f"  figures/{outs}")


if __name__ == "__main__":
    main()
