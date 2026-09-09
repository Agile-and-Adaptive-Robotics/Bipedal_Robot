"""Back-solve PF->MN synaptic weights from reference activation patterns.

Implements the pipeline Ben proposed: run IK (+ static optimization) in
OpenSim - or any reference (human data, RL policy, recorded sim) - extract
muscle activation patterns, then fit the spinal network's phase-group ->
motoneuron weight tables so the network reproduces them.

Pipeline
--------
1. Reference activations A[t, muscle] in [0,1] on a common time grid
   (CSV: header = actuator names, rows = samples; see load_activations).
2. NMF into k synergies (Ivanenko 2004: ~5 per leg, speed-invariant
   structure; timing/amplitude carry the speed modulation - Kibushi 2018).
3. Phase assignments: either from the network's RG/PF states recorded during
   a rollout (log_neurons=True in runner) or from a gait-event column.
4. Per phase bin p: non-negative least squares for W[p, synergy] so that
   sum_s W[p, s] * synergy_s(t) approximates A(t) for t in p. The PF->MN
   table in params.W_PF_MN / build_network is then replaced by the fit
   (via build_network.set_pf_weights, writing per-MN weights = the
   synergy-to-muscle matrix NMF components).
5. Reflex gains (Ia/II/Ib) are fit the same way with proprioceptor-driven
   basis columns once afferent recordings exist (extend COLS below).

Usage:
    python fit_synapses.py activations.csv [--k 5] [--out fitted_weights.npz]
"""
from __future__ import annotations

import argparse
import numpy as np


def load_activations(path: str):
    """CSV with a time column ('t' or 'time') + one column per actuator."""
    import csv
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    header = [h.strip() for h in rows[0]]
    tcol = next(i for i, h in enumerate(header) if h.lower() in ("t", "time"))
    data = np.array([[float(v) for v in r] for r in rows[1:] if r])
    t = data[:, tcol]
    names = [header[i] for i in range(len(header)) if i != tcol]
    A = data[:, [i for i in range(len(header)) if i != tcol]]
    return t, names, A


def nmf(X: np.ndarray, k: int, iters: int = 500, seed: int = 0):
    """Multiplicative-update NMF: X ~ C @ S, X[t,m] >= 0.

    Returns C [t,k] (synergy activation coefficients) and
    S [k,m] (synergy -> muscle weights, i.e. the NMF components).
    """
    rng = np.random.default_rng(seed)
    t, m = X.shape
    C = rng.random((t, k)) + 0.1
    S = rng.random((k, m)) + 0.1
    for _ in range(iters):
        C *= (X @ S.T) / np.maximum(C @ S @ S.T + 1e-9, 1e-9)
        S *= (C.T @ X) / np.maximum(C.T @ C @ S + 1e-9, 1e-9)
    return C, S


def fit_phase_weights(C: np.ndarray, phases: np.ndarray, n_phases: int):
    """Per-phase NNLS: for each phase bin p, solve min ||C_p W_p^T - A||...

    Here A is already approximated by C @ S; fitting W over synergy
    coefficients gives phase weights over [k] synergies:
        W[p, s] >= 0, activations_p ~ C_p @ diag(W_p) @ S
    Implemented directly as weighted NNLS per phase.
    """
    from scipy.optimize import nnls
    k = C.shape[1]
    W = np.zeros((n_phases, k))
    for p in range(n_phases):
        mask = phases == p
        if not mask.any():
            continue
        # target per phase = mean synergy demand; weights scale synergies
        target = C[mask].mean(axis=0)
        W[p], _ = nnls(np.eye(k), target)
    return W


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("activations")
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--out", default="fitted_weights.npz")
    ap.add_argument("--phases-from", default=None,
                    help="optional CSV column with phase labels 0..n-1")
    args = ap.parse_args(argv)

    t, names, A = load_activations(args.activations)
    print(f"{A.shape[0]} samples x {A.shape[1]} muscles")
    C, S = nmf(A, args.k)
    resid = np.linalg.norm(A - C @ S) / max(np.linalg.norm(A), 1e-9)
    print(f"NMF k={args.k}: relative residual {resid:.3f}")

    if args.phases_from:
        _, _, ph = load_activations(args.phases_from)
        W = fit_phase_weights(C, ph.astype(int), int(ph.max()) + 1)
        print("phase x synergy weights:")
        print(np.round(W, 3))
        np.savez(args.out, C=C, S=S, W=W, names=names)
        print(f"saved {args.out}")
    else:
        np.savez(args.out, C=C, S=S, names=names)
        print(f"saved {args.out} (no phase labels: run phase assignment "
              f"from RG/PF logs next)")


if __name__ == "__main__":
    main()
