"""Prune-campaign figures (2026-09-30): (1) leave-one-out delta bars per
variant/mode; (2) combo summary panel. Reads prune_results_*.jsonl (local
copies) + prune_analysis.md numbers; writes campaigns/20260930/figs/.
Reuses the scoring rules of analyze_prune.py.
"""
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
OUT = HERE / "campaigns" / "20260930" / "figs"
OUT.mkdir(parents=True, exist_ok=True)

CUTS = ["noaff", "interleg", "contact", "ia", "ii", "ib", "renshaw",
        "rgweak"]
PUSHES = ["+x", "-x", "+y", "-y"]
VARIANTS = ["s3k", "w2lvar", "syn6"]
COLORS = {"s3k": "#0072B2", "w2lvar": "#E69F00", "syn6": "#009E73"}


def _load(v):
    cells = {}
    fn = HERE / f"prune_results_{v}.jsonl"
    if not fn.exists():
        return cells
    for line in open(fn, encoding="utf-8", errors="replace"):
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            r = json.loads(line)
        except Exception:
            continue
        cells[(r["config"], r["mode"], r.get("push", ""))] = r
    return cells


def _walk_score(m):
    if not m or m.get("nan"):
        return -400.0
    if not (m.get("kine_components") or {}):
        return -320.0
    s = max(float(m.get("kine_score", -25.0)), -315.0)
    if m.get("kz", 1.0) < 0.62:
        s -= 20.0
    if m.get("tilt_max", 0.0) > 40.0:
        s -= 10.0
    return s


def main():
    # ---- Figure 1: WALK deltas (leave-one-out) per variant
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    for ax, v in zip(axes, VARIANTS):
        cells = _load(v)
        ref = _walk_score((cells.get(("full", "WALK", "")) or
                           {}).get("metrics", {}))
        xs, ys, cols = [], [], []
        for i, c in enumerate(CUTS):
            r = cells.get((c, "WALK", ""))
            if not r:
                continue
            d = _walk_score(r.get("metrics", {})) - ref
            xs.append(c)
            ys.append(d)
            cols.append("#D55E00" if d >= -3 else COLORS[v])
        if v == "s3k":
            r = cells.get(("combo", "WALK", ""))
            if r:
                xs.append("COMBO")
                ys.append(_walk_score(r.get("metrics", {})) - ref)
                cols.append("#CC79A7")
        ax.bar(range(len(xs)), ys, color=cols)
        ax.axhline(0, color="k", lw=0.8)
        ax.axhline(-3, color="gray", ls="--", lw=0.8)
        ax.set_xticks(range(len(xs)))
        ax.set_xticklabels(xs, rotation=45, ha="right", fontsize=8)
        ax.set_title(f"{v}\n(full walk = "
                     f"{_walk_score((cells.get(('full', 'WALK', '')) or {}).get('metrics', {})):.1f})",
                     fontsize=10)
        ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("walk kine  delta  (cut - full)\n"
                       "bars at/above dashed line = component prunable")
    fig.suptitle("Leave-one-out ablation: ground-walk cost of each "
                 "component (higher = less needed)", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "prune_walk_deltas.png", dpi=160)
    print("wrote prune_walk_deltas.png")

    # ---- Figure 2: AIR rhythm scores + stand + push recovery panel
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, v in zip(axes, VARIANTS):
        cells = _load(v)
        xs, ys, cols = [], [], []
        ref_air = (cells.get(("full", "AIR", "")) or {}) \
            .get("air", {}).get("air_score")
        for i, c in enumerate(CUTS):
            r = cells.get((c, "AIR", ""))
            if not r or "air" not in r:
                continue
            xs.append(c)
            ys.append(r["air"].get("air_score") or -200)
            cols.append(COLORS[v])
        ax.bar(range(len(xs)), ys, color=cols)
        ax.set_xticks(range(len(xs)))
        ax.set_xticklabels(xs, rotation=45, ha="right", fontsize=8)
        ax.set_title(f"{v} AIR (full = {ref_air:.1f})", fontsize=10)
        ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("air-stepping objective\n(rhythm gate: <3 bursts "
                       "collapses)")
    fig.suptitle("Air-stepping rhythm survives every cut? (flat bars = "
                 "rhythm died)", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "prune_air_scores.png", dpi=160)
    print("wrote prune_air_scores.png")

    # ---- Figure 3: push-recovery sway, all variants, FULL vs COMBO
    fig, ax = plt.subplots(figsize=(8, 4.5))
    w = 0.2
    for vi, v in enumerate(VARIANTS):
        cells = _load(v)
        ys = []
        for ax_ in PUSHES:
            m = (cells.get(("full", "PUSH", ax_)) or {}).get("metrics", {})
            ys.append(m.get("bal_push_sway", np.nan))
        ax.bar(np.arange(4) + (vi - 1) * w, ys, width=w,
               color=COLORS[v], label=v)
    cells = _load("s3k")
    ys = [(cells.get(("combo", "PUSH", ax_)) or {}).get("metrics", {})
          .get("bal_push_sway", np.nan) for ax_ in PUSHES]
    ax.bar(np.arange(4) + 2 * w, ys, width=w, color="#CC79A7",
           label="s3k COMBO (pruned)")
    ax.set_xticks(range(4))
    ax.set_xticklabels([f"push {a} (40 N)" for a in PUSHES])
    ax.set_ylabel("post-push COM sway radius [m]")
    ax.set_title("Push recovery: sway radius after the pulse "
                 "(lower = tighter recovery)")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "prune_push_recovery.png", dpi=160)
    print("wrote prune_push_recovery.png")


if __name__ == "__main__":
    main()
