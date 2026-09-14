"""Refit W_PF_MN + W_POSTURE from the IK/NNLS back-solved human pattern.

Pipeline (Ben 2026-09-11 staged plan, step after bsolve_ik.py):
 1. Load the back-solved activations (bsolve_out.npz) and bin them by gait
    phase into per-functional-group profiles (already stored per group).
 2. Extract the network's own PF-cell phase windows (E1/E2/F1/F2 potentials
    vs RG_E phase) from the last recorded run (spinal_run.npz, walk window).
 3. Per functional group, per leg: NNLS for non-negative weights over the 4
    PF windows so their mix reproduces the human group activation profile.
    Right and left fits are averaged into one table (the network shares it).
 4. W_POSTURE: tonic set to half the cycle-mean of each group's human
    activation (the PF table carries the full phase profile; the old values
    were as large as the human PEAK - the direct cause of the hip-extensor
    saturation / pelvis-limbo).
 5. Write fitted_walk_params.json; runner --fitted and optuna_walk.py v3
    consume it.

Usage: python fit_pf.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import nnls

HERE = Path(__file__).parent


def pf_windows_from_run(npz, side="r", n_bins=20):
    """PF-cell phase windows over the RG_E cycle, from a recorded run."""
    neuro, names, t = npz["neuro"], list(npz["neuro_names"]), npz["t"]
    col = {n: i for i, n in enumerate(names)}
    s = side.upper()
    rg = neuro[:, col[f"RG_E_{side}"]]
    # walk window only
    walk = (t > 5.0) & (t < 15.0)
    rgw = np.where(walk, rg, 0.0)
    on = rgw > 0.5 * rgw.max()
    rises = np.flatnonzero(np.diff(on.astype(int)) == 1)
    falls = np.flatnonzero(np.diff(on.astype(int)) == -1)
    if len(rises) < 3:
        print("WARNING: no clean RG cycles in spinal_run.npz - "
              "using analytic gaussian windows")
        ph = np.linspace(0, 1, n_bins)
        cen = dict(E1=0.08, E2=0.38, F1=0.66, F2=0.86)
        wdt = dict(E1=0.07, E2=0.12, F1=0.09, F2=0.12)
        return {p: np.exp(-0.5 * ((ph - c) / wdt[p]) ** 2)
                for p, c in cen.items()}
    cyc = np.diff(t[rises])
    print(f"PF windows from {len(rises)} RG_E cycles, period "
          f"{cyc.mean():.3f} s")
    bins = [[] for _ in range(n_bins)]
    for k in range(len(t) - 1):
        if not walk[k]:
            continue
        prev = t[rises][t[rises] <= t[k]]
        if len(prev) < 2:
            continue
        t0 = prev[-1]
        nxt = t[rises][t[rises] > t[k]]
        if len(nxt) == 0:
            continue
        ph = (t[k] - t0) / (nxt[0] - t0)
        if 0 <= ph < 1:
            bins[min(int(ph * n_bins), n_bins - 1)].append(k)
    out = {}
    for p in ("E1", "E2", "F1", "F2"):
        v = neuro[:, col[f"PF_{p}_{side}"]]
        prof = np.array([np.max(v[b]) if len(b) else 0.0 for b in bins])
        if prof.max() <= 0:
            prof = np.ones(n_bins) * 0.01
        out[p] = np.clip(prof / prof.max(), 0.0, 1.0)
    return out


def fit_leg(grp_prof, win, n_bins=20):
    """Per-group NNLS weights over PF windows; returns {group: {phase: w}}."""
    P = np.column_stack([win[p] for p in ("E1", "E2", "F1", "F2")])  # [bins,4]
    tab = {}
    for g, target in grp_prof.items():
        w, _ = nnls(P, target)
        resid = np.linalg.norm(P @ w - target) / max(np.linalg.norm(target), 1e-9)
        tab[g] = (dict(zip(("E1", "E2", "F1", "F2"), np.round(w, 4))),
                  float(resid), float(target.max()))
    return tab


def main():
    d = np.load(HERE / "bsolve_out.npz", allow_pickle=True)
    run = np.load(HERE / "spinal_run.npz", allow_pickle=True)
    if not d["grp_prof_names"].shape[0]:
        raise SystemExit("bsolve_out.npz has no group profiles - "
                         "no gait cycle detected; run bsolve_ik.py first")

    grp_r = {n: d["grp_prof"][i] for i, n in enumerate(d["grp_prof_names"])}
    grp_l = {n: d["grp_prof_l"][i] for i, n in enumerate(d["grp_prof_l_names"])}

    win_r = pf_windows_from_run(run, "r")
    # the runner only logs right-side PF cells; the network builds the same
    # window shapes per side, so the right windows stand in for the left fit
    win_l = win_r

    fit_r = fit_leg(grp_r, win_r)
    fit_l = fit_leg(grp_l, win_l)

    # merge (network table is shared across sides): average where both
    W_PF_MN = {p: {} for p in ("E1", "E2", "F1", "F2")}
    report = []
    groups = sorted(set(fit_r) | set(fit_l))
    for g in groups:
        ws = []
        for fit in (fit_r, fit_l):
            if g in fit:
                ws.append(fit[g][0])
        avg = {p: float(np.mean([w[p] for w in ws])) for p in
               ("E1", "E2", "F1", "F2")}
        resid = np.mean([fit[g][1] for fit in (fit_r, fit_l) if g in fit])
        peak = max([fit[g][2] for fit in (fit_r, fit_l) if g in fit])
        for p in ("E1", "E2", "F1", "F2"):
            W_PF_MN[p][g] = round(avg[p], 4)
        report.append((g, avg, resid, peak))

    # posture tonic: half the human cycle-mean (PF table carries the peak)
    W_POSTURE = {}
    for g in groups:
        profs = [x[g] for x in (grp_r, grp_l) if g in x]
        mean_act = float(np.mean([p.mean() for p in profs]))
        W_POSTURE[g] = round(0.5 * mean_act, 4)

    print(f"\n{'group':10s} {'E1':>6s} {'E2':>6s} {'F1':>6s} {'F2':>6s} "
          f"{'resid':>6s} {'peak':>5s} {'posture':>8s}")
    for g, avg, resid, peak in report:
        print(f"{g:10s} {avg['E1']:6.3f} {avg['E2']:6.3f} {avg['F1']:6.3f} "
              f"{avg['F2']:6.3f} {resid:6.3f} {peak:5.2f} "
              f"{W_POSTURE[g]:8.3f}")

    # saturation audit: hip_ext old stack vs human peak
    import params
    old = (params.W_PF_MN["E1"]["hip_ext"] + params.W_PF_MN["E2"]["hip_ext"]
           + params.W_POSTURE["hip_ext"])
    new = (W_PF_MN["E1"]["hip_ext"] + W_PF_MN["E2"]["hip_ext"]
           + W_POSTURE["hip_ext"])
    print(f"\nhip_ext stance drive stack: old {old:.3f} -> fitted {new:.3f} "
          f"(human peak {max(fit_r['hip_ext'][2], fit_l['hip_ext'][2]):.2f})")

    out = dict(W_PF_MN=W_PF_MN, W_POSTURE=W_POSTURE, pf_gain=1.0,
               meta=dict(
        source="bsolve_ik.py NNLS back-solve of subject01 IK + GRF",
        duty_r=float(d["duty_r"]), fit="per-group NNLS over PF windows",
        note="runner --fitted / optuna_walk v3 consume this; pf_gain is "
             "the global table gain (1.0 = as-fitted amplitudes)"))
    with open(HERE / "fitted_walk_params.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print("saved fitted_walk_params.json")


if __name__ == "__main__":
    main()
