"""Publication-grade schematics of the gait2392 spinal circuit.

Style follows Ben's reference pictures (Rybak / SNS-hierarchy style,
2026-09-12):
  * tinted layer bands with plain labels (no boxes everywhere)
  * populations drawn as clusters of small circles, real counts pulled
    LIVE from muscle_map.py (never hardcoded)
  * size hierarchy: half-centers big, pattern cells medium, sensors and
    local interneurons small; small loop cells for adaptation
  * strict rows in the dense panel, mirrored left/right around a dashed
    midline, commissural wires gathered in one pale central corridor
  * wires colored by source (Okabe-Ito CVD-safe); synapse markers keep
    Ben's conventions: white triangle with its BASE flat against the
    excited cell (snsfig.m 2026-09-09 inversion), solid black dot =
    inhibitory
  * clean wires: no conductance numbers on schematics -- numbers live in
    the weights figure

Figures (--which):
  core     one side, hierarchy: DRIVE/POSTURE -> RG half-centers (+ADAP)
           -> PF cells (+PFA) -> MN pool clusters -> muscles; Ia/II/Ib
           sensor triplet (x46 each) + IB-EXC load-sharing cells.
  full     both sides drawn fully (no ghost half): mirrored rows,
           midline commissural corridor, BAL family, trunk pools.
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
    if source in ("fitted", "best"):
        fit = json.loads((HERE / "fitted_walk_params.json").read_text("utf-8"))
        for ph, tbl in fit["W_PF_MN"].items():
            for g, w in tbl.items():
                W[ph][g] = float(w)
        for g, w in fit["W_POSTURE"].items():
            WPOST[g] = float(w)
        note = "fitted_walk_params.json (IK/NNLS back-solve refit)"
    if source == "best":
        best = json.loads((HERE / "best_walk_params.json").read_text("utf-8"))
        gain = float(best.get("pf_gain", 1.0))
        if gain != 1.0:
            for tbl in W.values():
                for g in tbl:
                    tbl[g] *= gain
            for g in WPOST:
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
# structural cell counts per side (build_network.py layout)
N_RG, N_ADAP, N_PF, N_PFA, N_IBEXC = 2, 2, 4, 4, 5
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
W_E, W_F = OI_VERM, OI_BLUE          # extensor / flexor half-center families
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
            self.fig.savefig(p, dpi=300, bbox_inches="tight", pad_inches=0.04)
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
    fs = fs if fs else max(6.0, 10.0 * r / R_MED)
    cv.ax.text(x, y, label, ha="center", va="center", fontsize=fs, zorder=4)
    if sub:
        cv.ax.text(x, y - r - 0.16, sub, ha="center", va="center",
                   fontsize=fs - 1.0, style="italic", color="0.35", zorder=4)
    return (x, y, r)


def input_box(cv, x, y, w, h, label, sub=None):
    cv.ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                   boxstyle="round,pad=0.02",
                                   fc="0.94", ec="black", lw=1.2, zorder=3))
    cv.ax.text(x, y + (0.09 if sub else 0), label, ha="center",
               va="center", fontsize=7.8, zorder=4)
    if sub:
        cv.ax.text(x, y - 0.13, sub, ha="center", va="center",
                   fontsize=6.6, style="italic", color="0.35", zorder=4)
    return (x, y, min(w, h) / 2)


def muscle(cv, x, y, w=0.72, h=0.26, label=None):
    cv.ax.add_patch(Ellipse((x, y), w, h, fc=_tint(OI_SKY), ec="black",
                            lw=1.0, zorder=3))
    if label:
        cv.ax.text(x, y, label, ha="center", va="center", fontsize=6.0,
                   zorder=4)
    return (x, y, h / 2)


def pool_cluster(cv, x, y, group, cap=5):
    """MN pool cluster: small circles = real motoneuron pools (primary
    members), + a small pale ellipse per biarticular secondary member.
    Returns anchor dict (top/left/right/bottom points)."""
    n_prim = len(CNT["prim"][group])
    n_sec = len(CNT["sec"][group])
    k = min(n_prim, cap)
    xs = [x + (i - (k - 1) / 2) * (2 * R_POOL + 0.055) for i in range(k)]
    for xc in xs:
        cv.ax.add_patch(Circle((xc, y), R_POOL, fc=_tint(OI_SKY),
                               ec="black", lw=0.9, zorder=3))
    # secondary (biarticular) pools: one shared pale ellipse to the right
    if n_sec:
        cv.ax.add_patch(Ellipse((xs[-1] + 0.42, y), 0.34, 0.20,
                                fc="white", ec="0.45", lw=0.8, zorder=3))
    lbl = ABBREV[group]
    cnt_txt = f"x{n_prim}" + (f" (+{n_sec})" if n_sec else "")
    cv.ax.text(x, y - 0.42, f"{lbl}  {cnt_txt}", ha="center", va="center",
               fontsize=6.8, zorder=4)
    half_w = (k - 1) / 2 * (2 * R_POOL + 0.055) + R_POOL + (0.5 if n_sec else 0)
    return dict(c=(x, y), top=(x, y + R_POOL + 0.06),
                left=(x - half_w, y), right=(x + half_w, y),
                bot=(x, y - R_POOL - 0.06), hw=half_w)


def layer_band(cv, x0, x1, y0, y1, hex_color, label, lab_x=None,
               lab_above=True, fs=8.0):
    cv.ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                   boxstyle="round,pad=0.02,rounding_size=0.18",
                                   fc=_tint(hex_color), ec="none", zorder=0.5))
    lx = lab_x if lab_x is not None else x0 + 0.15
    ly = y1 + 0.14 if lab_above else y1 - 0.20
    cv.ax.text(lx, ly, label, fontsize=fs, ha="left", va="center",
               color="0.30", style="italic", zorder=1)


def syn(cv, p, q, exc: bool, label=None, lpos=0.5, loff=(0.1, 0.1), lw=1.4,
        rad=0.0, ghost=False, dashed=False, shrink_a=None, shrink_b=None,
        color=None, lbox=True, lfs=6.4):
    """Edge p->q. Excitatory: WHITE TRIANGLE WITH ITS BASE FLAT AGAINST THE
    TARGET CELL, apex pointing back along the wire (Ben's snsfig.m
    convention, 2026-09-09). Inhibitory: solid black dot at the target.
    p/q may be (x, y) or the (x, y, r) tuple returned by the glyph
    helpers; shrink_a/shrink_b override the source/target radii."""
    col = "0.55" if ghost else (color if color else ("0.15" if exc else "0.15"))
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
    if exc:
        # base plane just off the target edge; triangle extends back along
        # the wire; the wire stops at the apex
        base = qa - u * (rb + 0.02)
        tri = min(0.20, max(0.08, 0.45 * (L - ra - rb)))
        apex = base - u * tri
        cv.ax.add_patch(FancyArrowPatch(
            tuple(a), tuple(apex), arrowstyle="-", lw=lw, color=col, zorder=2,
            connectionstyle=f"arc3,rad={rad}",
            linestyle=(0, (3, 2)) if dashed else "solid"))
        if not ghost:
            cv.ax.add_patch(Polygon(
                [tuple(base + perp * 0.095), tuple(base - perp * 0.095),
                 tuple(apex)], closed=True, fc="white", ec=col, lw=1.2,
                zorder=3))
    else:
        dot = qa - u * (rb + 0.09)
        cv.ax.add_patch(FancyArrowPatch(
            tuple(a), tuple(dot), arrowstyle="-", lw=lw, color=col, zorder=2,
            connectionstyle=f"arc3,rad={rad}",
            linestyle=(0, (3, 2)) if dashed else "solid"))
        if not ghost:
            cv.ax.add_patch(Circle(tuple(dot), 0.085, fc="black", ec="black",
                                   zorder=3))
    if label and not ghost:
        mid = pa + (qa - pa) * lpos
        if rad:
            mid = mid + perp * rad * L * 0.5
        mid = mid + np.asarray(loff, float)
        cv.ax.text(*mid, label, fontsize=lfs, ha="center", va="center",
                   style="italic", color="0.25", zorder=5,
                   bbox=LBL_BBOX if lbox else None)


def tag(cv, x, y, text, fs=6.6, color="0.30", ha="center"):
    cv.ax.text(x, y, text, fontsize=fs, ha=ha, va="center", style="italic",
               color=color, zorder=4)


def legend_row(cv, y, x0=0.6):
    ax = cv.ax
    ax.add_patch(Circle((x0, y), 0.16, fc="white", ec="black", lw=1.3))
    ax.text(x0 + 0.26, y, "neuron", fontsize=7.4, va="center")
    x = x0 + 1.35
    # excitatory: wire ending in a triangle whose base faces the target
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
    """Small colored-wire key: source/sign colors."""
    ax = cv.ax
    items = [(W_E, "extensor-half drive"), (W_F, "flexor-half drive"),
             (W_EXC, "other excitatory"), (W_INH, "inhibitory"),
             (W_DESC, "descending / balance")]
    for i, (c, name) in enumerate(items):
        xx = x + i * 2.35
        ax.add_patch(FancyArrowPatch((xx, y), (xx + 0.5, y), arrowstyle="-",
                                     lw=1.6, color=c))
        ax.text(xx + 0.6, y, name, fontsize=6.6, va="center", color="0.25")


def key_box(cv, x, y):
    """Group key: abbreviation -> real member muscles (live)."""
    lines = []
    for g in GROUP_ORDER:
        prim = ", ".join(CNT["prim"][g])
        extra = f"  (+{', '.join(CNT['sec'][g])} secondary)" if CNT["sec"][g] else ""
        lines.append(f"{ABBREV[g]} ({len(CNT['prim'][g])}): {prim}{extra}")
    txt = "MN pool key (live from muscle_map.py)\n" + "\n".join(lines)
    cv.ax.text(x, y, txt, fontsize=5.9, va="top", ha="left", color="0.25",
               family="monospace", zorder=4,
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


# ------------------------------------------------------- shared blocks
def rg_block(cv, x, y, suffix="", big=True, mirror=False):
    """Half-centers + adaptation loop cells. Returns (rg_e, rg_f)."""
    r = R_BIG if big else R_MED
    dx = -1.55 if not mirror else 1.55
    rge = neuron(cv, x + dx, y + 0.62, f"RG-E{suffix}", fc=_tint(OI_VERM), r=r)
    rgf = neuron(cv, x + dx, y - 0.78, f"RG-F{suffix}", fc=_tint(OI_BLUE), r=r)
    syn(cv, rge, rgf, False, rad=0.25, color=W_INH)
    syn(cv, rgf, rge, False, rad=0.25, color=W_INH)
    tag(cv, x + dx + 0.55, y - 0.08, "mutual\ninhibition", fs=6.0)
    # adaptation loop cells (fatigue/adaptation onto each half-center)
    sgn = 1 if not mirror else -1
    ade = neuron(cv, x + dx - sgn * 1.45, y + 0.62, "ADAP-E", r=R_SMALL, fs=6.0)
    adf = neuron(cv, x + dx - sgn * 1.45, y - 0.78, "ADAP-F", r=R_SMALL, fs=6.0)
    syn(cv, rge, ade, True, rad=-0.3, lw=1.0, color=W_E)
    syn(cv, ade, rge, False, rad=-0.3, lw=1.0, color=W_INH)
    syn(cv, rgf, adf, True, rad=0.3, lw=1.0, color=W_F)
    syn(cv, adf, rgf, False, rad=0.3, lw=1.0, color=W_INH)
    return rge, rgf


def pf_block(cv, x, y, suffix="", mirror=False):
    """Four PF phase cells + PFA loop cells + reciprocal inhibition.
    Returns dict of PF cells."""
    names = (("E1", 1.5, OI_VERM), ("E2", 0.5, OI_VERM),
             ("F1", -0.6, OI_BLUE), ("F2", -1.6, OI_BLUE))
    out = {}
    for nm, dy, c in names:
        out[nm] = neuron(cv, x, y + dy, f"PF-{nm}{suffix}", fc=_tint(c))
    syn(cv, out["E2"], out["F1"], False, color=W_INH, rad=0.12)
    syn(cv, out["F1"], out["E2"], False, color=W_INH, rad=0.12)
    syn(cv, out["F2"], out["E2"], False, color=W_INH, rad=-0.35)
    syn(cv, out["E2"], out["F2"], False, color=W_INH, rad=-0.35)
    syn(cv, out["E1"], out["F1"], False, color=W_INH, rad=-0.3)
    # PFA adaptation loop cells under each PF cell (small)
    sgn = 1 if not mirror else -1
    for nm in out:
        pfa = neuron(cv, x + sgn * 1.05, out[nm][1], f"", r=0.13,
                     fc="white", lw=0.9)
        cv.ax.text(pfa[0], pfa[1] - 0.32, f"PFA-{nm}", fontsize=5.6,
                   ha="center", color="0.35", zorder=4)
        syn(cv, out[nm], pfa, True, lw=0.8,
            color=W_E if nm[0] == "E" else W_F)
        syn(cv, pfa, out[nm], False, lw=0.8, color=W_INH, rad=0.25)
    return out


def sensor_triplet(cv, x, y, suffix=""):
    """Ia/II/Ib population circles (x46 each per side)."""
    ia = neuron(cv, x - 1.5, y, f"Ia{suffix}", f"x{N_PER_SIDE}", r=R_SENS,
                fc=_tint(OI_ORANGE), fs=6.6)
    ii = neuron(cv, x, y, f"II{suffix}", f"x{N_PER_SIDE}", r=R_SENS,
                fc=_tint(OI_ORANGE), fs=6.6)
    ib = neuron(cv, x + 1.5, y, f"Ib{suffix}", f"x{N_PER_SIDE}", r=R_SENS,
                fc=_tint(OI_ORANGE), fs=6.6)
    tag(cv, x, y - 0.62, "muscle sensors: one Ia, II, Ib per pool "
        f"({N_PER_SIDE} each per side)", fs=6.0)
    return dict(ia=ia, ii=ii, ib=ib)


def ibexc_cluster(cv, x, y, suffix=""):
    """The five load-sharing interneurons (one per extensor-stance group)."""
    xs = [x + (i - 2) * 0.42 for i in range(N_IBEXC)]
    for xc in xs:
        cv.ax.add_patch(Circle((xc, y), R_SMALL, fc=_tint(OI_GREEN),
                               ec="black", lw=0.9, zorder=3))
    cv.ax.text(x, y, "IB-EXC", ha="center", va="center", fontsize=5.4,
               zorder=4)
    tag(cv, x, y - 0.42, f"{N_IBEXC} cells: one per extensor-stance group\n"
        "(stance-gated Ib reversal, load sharing)", fs=6.0)
    return (x, y, R_SMALL)


# ------------------------------------------------------------ fig: core
def fig_core(t, fmts):
    W = t["W"]
    cv = Canvas(14.0, 13.6)

    # ---------------- inputs ----------------
    drv = input_box(cv, 1.25, 12.9, 1.7, 0.6, "DRIVE", "MLR surrogate")
    pos = input_box(cv, 3.85, 12.9, 1.7, 0.6, "POSTURE", "tonic")
    post_i = input_box(cv, 6.45, 12.9, 1.75, 0.6, "POST$_i$",
                       "solved standing bias")

    # ---------------- rhythm generator ----------------
    layer_band(cv, 0.5, 8.3, 9.35, 11.95, OI_ORANGE,
               "RHYTHM GENERATOR (RG)  --  half-centers set timing & duty")
    rge, rgf = rg_block(cv, 4.0, 10.65)
    syn(cv, drv, rge, True, color=W_DESC, rad=0.1)
    syn(cv, drv, rgf, True, color=W_DESC, rad=-0.15)
    syn(cv, pos, rge, True, color=W_DESC, rad=-0.2)
    tag(cv, 7.6, 11.5, "frequency rises with DRIVE;\nstance-biased drive +\n"
        "adaptation set the duty")

    # ---------------- pattern formation ----------------
    layer_band(cv, 0.5, 8.3, 5.7, 9.1, OI_SKY,
               "PATTERN FORMATION (PF)  --  four phase-window cells per side")
    pfs = pf_block(cv, 2.6, 7.4)
    syn(cv, rge, pfs["E1"], True, color=W_E, rad=-0.1)
    syn(cv, rge, pfs["E2"], True, color=W_E, rad=-0.2)
    syn(cv, rgf, pfs["F1"], True, color=W_F, rad=-0.1)
    syn(cv, rgf, pfs["F2"], True, color=W_F, rad=-0.2)
    syn(cv, drv, pfs["E1"], True, color=W_DESC, rad=0.3, dashed=True)
    tag(cv, 6.9, 7.4, "PF cells = phase windows;\n(tau, adaptation) stagger\n"
        "their bursts; PFA loops\nshape each cell")
    # interneurons row
    ib = ibexc_cluster(cv, 6.9, 6.35, "")
    syn(cv, rge, ib, True, color=W_E, rad=0.25, dashed=True,
        label="stance gate", lpos=0.45, loff=(-0.3, 0.25))

    # ---------------- sensors ----------------
    sens = sensor_triplet(cv, 3.4, 4.55)

    # ---------------- motoneuron pools ----------------
    layer_band(cv, 0.5, 13.5, 1.7, 3.5, OI_GREEN,
               "MOTONEURON POOLS  --  one cluster per functional group "
               "(real pool counts)")
    pools, x = {}, 1.55
    for g in GROUP_ORDER:
        pools[g] = pool_cluster(cv, x, 2.75, g)
        x += 1.30 if len(CNT["prim"][g]) <= 5 else 1.42
    # PF -> MN fan (color by phase family, thickness by weight, no numbers)
    for ph, node in pfs.items():
        for g, w in W[ph].items():
            if g in pools and w >= W_MIN:
                syn(cv, node, pools[g]["top"], True,
                    lw=0.7 + min(w, 0.25) * 6, shrink_b=0.06,
                    color=W_E if ph[0] == "E" else W_F)
    # POST_i dashed standing bias into every pool (one drawn + tag)
    syn(cv, post_i, pools["HIP-E"]["bot"], True, color=W_DESC, dashed=True,
        rad=0.18, shrink_b=0.04)
    syn(cv, post_i, pools["ANK-DF"] if "ANK-DF" in pools else
        pools[list(pools)[-1]]["bot"], True, color=W_DESC, dashed=True,
        rad=-0.1, shrink_b=0.04)
    tag(cv, 9.6, 4.1, "dashed: POST$_i$ standing bias into every pool",
        ha="left")
    # reflex pathways (representative; every pool has Ia/II/Ib)
    syn(cv, sens["ia"], pools["KNEE-E"]["top"], True, color=W_EXC, rad=0.15,
        label="Ia homonymous", lpos=0.4, loff=(-0.45, 0.15))
    syn(cv, sens["ia"], pools["KNEE-F"]["top"], False, color=W_INH, rad=-0.2,
        label="Ia reciprocal", lpos=0.4, loff=(-0.5, -0.1))
    syn(cv, sens["ii"], pools["HIP-E"]["top"], True, color=W_EXC, rad=0.25,
        label="II (stance-gated)", lpos=0.5, loff=(0.5, 0.2))
    syn(cv, sens["ib"], pools["ANK-PF"]["top"], False, color=W_INH, rad=0.12,
        label="Ib autogenic", lpos=0.35, loff=(0.55, 0.0))
    syn(cv, sens["ib"], ib, True, color=W_EXC, rad=-0.25, lw=1.0)
    syn(cv, ib, pools["HIP-E"]["top"], True, color=W_EXC, rad=0.3, lw=1.0,
        label="Ib reversal -> extensors", lpos=0.6, loff=(0.7, 0.15))
    # muscles row
    mx = 1.55
    for g in GROUP_ORDER:
        muscle(cv, mx, 1.15)
        syn(cv, pools[g]["bot"], (mx, 1.30), True, color="0.4", lw=0.8,
            shrink_a=0.04, shrink_b=0.10)
        mx += 1.30 if len(CNT["prim"][g]) <= 5 else 1.42

    # ---------------- legend / key / source ----------------
    legend_row(cv, 0.55)
    wire_color_key(cv, 0.7, 0.02)
    key_box(cv, 0.6, 12.75)
    source_note(cv, 0.55, 13.35, f"source: {t['note']}")
    return cv.save("circuit_core", fmts)


# ------------------------------------------------------------ fig: full
def fig_full(t, fmts):
    W = t["W"]
    cv = Canvas(19.0, 13.4)
    MID = 9.5

    # midline + commissural corridor
    cv.ax.plot([MID, MID], [1.4, 11.7], ls=(0, (4, 3)), color="0.6", lw=1.0,
               zorder=1)
    layer_band(cv, MID - 0.55, MID + 0.55, 6.0, 11.7, OI_PURPLE, "", )
    cv.ax.text(MID, 11.95, "interleg commissural inhibition",
               fontsize=7.0, ha="center", va="center", color="0.35",
               style="italic", zorder=1)
    tag(cv, MID, 6.3, "F$\\leftrightarrow$F strong\nE$\\leftrightarrow$E weak",
        fs=6.2)

    # ---------------- shared inputs (top) ----------------
    drv = input_box(cv, MID - 1.9, 12.9, 1.6, 0.6, "DRIVE", "MLR surrogate")
    pos = input_box(cv, MID + 1.9, 12.9, 1.6, 0.6, "POSTURE", "+ POST$_i$")

    # ---------------- per-side stacks ----------------
    sides = {"r": dict(x0=10.4, x1=18.6, mirror=False, sfx="$_r$"),
             "l": dict(x0=0.4, x1=8.6, mirror=True, sfx="$_l$")}
    rg, pfs, sens, ib, pools = {}, {}, {}, {}, {}
    for sd, cfg in sides.items():
        cx = (cfg["x0"] + cfg["x1"]) / 2          # side center
        # balance inputs at the outer edge
        bx = cfg["x1"] if not cfg["mirror"] else cfg["x0"]
        b1 = input_box(cv, bx - (0.95 if not cfg["mirror"] else -0.95),
                       12.05, 1.55, 0.55, "BAL$_{PF/DF}$", "sagittal COM")
        b2 = input_box(cv, bx - (0.95 if not cfg["mirror"] else -0.95),
                       11.25, 1.55, 0.55, "BAL$_{LAT}$/TRK",
                       "abduct / trunk")
        # RG
        layer_band(cv, cfg["x0"], cfg["x1"], 8.55, 10.45, OI_ORANGE,
                   f"RHYTHM GENERATOR {sd.upper()}")
        rge, rgf = rg_block(cv, cx, 9.5, suffix=cfg["sfx"],
                            mirror=cfg["mirror"])
        rg[sd] = (rge, rgf)
        syn(cv, drv, rge, True, color=W_DESC, rad=0.15)
        syn(cv, drv, rgf, True, color=W_DESC, rad=-0.2)
        syn(cv, pos, rge, True, color=W_DESC, rad=-0.25)
        # PF
        layer_band(cv, cfg["x0"], cfg["x1"], 6.35, 8.35, OI_SKY,
                   f"PATTERN FORMATION {sd.upper()}")
        pf = pf_block(cv, cx - (0.9 if not cfg["mirror"] else -0.9), 7.35,
                      suffix=cfg["sfx"], mirror=cfg["mirror"])
        pfs[sd] = pf
        syn(cv, rge, pf["E1"], True, color=W_E)
        syn(cv, rge, pf["E2"], True, color=W_E, rad=-0.15)
        syn(cv, rgf, pf["F1"], True, color=W_F)
        syn(cv, rgf, pf["F2"], True, color=W_F, rad=-0.15)
        # interneurons + sensors
        layer_band(cv, cfg["x0"], cfg["x1"], 4.55, 6.15, OI_GREEN,
                   f"LOAD-SHARING IN + SENSORS {sd.upper()}")
        ib[sd] = ibexc_cluster(cv, cx - (1.7 if not cfg["mirror"] else -1.7),
                               5.55, suffix=cfg["sfx"])
        sens[sd] = sensor_triplet(
            cv, cx + (1.7 if not cfg["mirror"] else -1.7), 5.55,
            suffix=cfg["sfx"])
        syn(cv, rge, ib[sd], True, color=W_E, dashed=True, rad=0.2, lw=1.0)
        syn(cv, sens[sd]["ib"], ib[sd], True, color=W_EXC, lw=1.0)
        # MN pools
        layer_band(cv, cfg["x0"], cfg["x1"], 1.75, 4.35, OI_SKY,
                   f"MOTONEURON POOLS {sd.upper()} (real counts)")
        pools[sd], x = {}, cfg["x0"] + 0.85
        for g in GROUP_ORDER:
            pools[sd][g] = pool_cluster(cv, x, 3.35, g, cap=4)
            x += 0.98 if len(CNT["prim"][g]) <= 4 else 1.06
        for ph, node in pf.items():
            for g, w in W[ph].items():
                if g in pools[sd] and w >= W_MIN:
                    syn(cv, node, pools[sd][g]["top"], True,
                        lw=0.7 + min(w, 0.25) * 6, shrink_b=0.05,
                        color=W_E if ph[0] == "E" else W_F)
        # representative reflex + balance wiring
        syn(cv, sens[sd]["ia"], pools[sd]["KNEE-E"]["top"], True,
            color=W_EXC, lw=0.9, rad=0.1)
        syn(cv, sens[sd]["ia"], pools[sd]["KNEE-F"]["top"], False,
            color=W_INH, lw=0.9, rad=-0.15)
        syn(cv, sens[sd]["ib"], pools[sd]["ANK-PF"]["top"], False,
            color=W_INH, lw=0.9, rad=0.1)
        syn(cv, ib[sd], pools[sd]["HIP-E"]["top"], True, color=W_EXC,
            lw=0.9, rad=0.25)
        syn(cv, b1, pools[sd]["ANK-PF"]["top"], True, color=W_DESC,
            dashed=True, rad=0.3, lw=1.0)
        syn(cv, b2, pools[sd]["HIP-AB"]["top"], True, color=W_DESC,
            dashed=True, rad=0.35, lw=1.0)
        syn(cv, b2, pools[sd]["TRK-E"]["top"], True, color=W_DESC,
            dashed=True, rad=0.25, lw=1.0)
        syn(cv, pos, pools[sd]["TRK-F"]["top"], True, color=W_DESC,
            dashed=True, rad=-0.3, lw=0.9)

    # commissural wires through the corridor
    syn(cv, rg["r"][1], rg["l"][1], False, color=W_INH, lw=2.0, rad=-0.18)
    syn(cv, rg["r"][0], rg["l"][0], False, color=W_INH, lw=1.1, rad=0.15)

    # muscles row
    for sd, cfg in sides.items():
        mx = cfg["x0"] + 0.85
        for g in GROUP_ORDER:
            muscle(cv, mx, 1.25, w=0.5, h=0.2)
            syn(cv, pools[sd][g]["bot"], (mx, 1.39), True, color="0.4",
                lw=0.7, shrink_a=0.03, shrink_b=0.09)
            mx += 0.98 if len(CNT["prim"][g]) <= 4 else 1.06

    # legend / key / source
    legend_row(cv, 0.55)
    wire_color_key(cv, 9.0, 0.02)
    key_box(cv, 15.0, 12.75)
    source_note(cv, 0.55, 13.15, f"source: {t['note']}")
    return cv.save("circuit_full", fmts)


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
                    choices=["core", "full", "weights", "all"])
    ap.add_argument("--source", default="best",
                    choices=["params", "fitted", "best"])
    ap.add_argument("--fmt", default="pdf,svg,png")
    args = ap.parse_args()
    fmts = [f.strip() for f in args.fmt.split(",") if f.strip()]
    OUT.mkdir(exist_ok=True)
    t = effective_tables(args.source)
    print(f"source: {t['note']}")
    print(f"network: {N_TOTAL} neurons ({N_SIDE}/side + {N_SHARED} shared), "
          f"{N_INPUTS} inputs, pools live from muscle_map.py")
    builders = dict(core=fig_core, full=fig_full, weights=fig_weights)
    which = list(builders) if args.which == "all" else [args.which]
    for nm in which:
        outs = builders[nm](t, fmts)
        print(f"  figures/{outs}")


if __name__ == "__main__":
    main()
