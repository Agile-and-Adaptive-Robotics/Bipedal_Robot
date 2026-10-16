"""Held-out and multi-seed robustness audit for the PF/synergy count."""
import json
import warnings

import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import NMF
from sklearn.exceptions import ConvergenceWarning

import fsa_backsolve as fsa


RANKS = range(3, 9)
SEEDS = range(10)


def fit(train, test, rank, seed):
    model = NMF(n_components=rank, init="random", random_state=seed,
                max_iter=10000, tol=1e-6)
    model.fit(train)
    coeff = model.transform(test)
    return coeff @ model.components_, model.components_


def cosine_rows(a, b):
    aa = a / np.maximum(np.linalg.norm(a, axis=1, keepdims=True), 1e-12)
    bb = b / np.maximum(np.linalg.norm(b, axis=1, keepdims=True), 1e-12)
    return aa @ bb.T


def audit_side(x, side):
    splits = ((np.arange(len(x)) % 2 == 0),
              (np.arange(len(x)) % 2 == 1))
    rows = []
    for rank in RANKS:
        values = []
        for test_mask in splits:
            train, test = x[~test_mask], x[test_mask]
            for seed in SEEDS:
                pred, _ = fit(train, test, rank, seed)
                values.append(fsa.centered_vaf(test, pred))
        rows.append((rank, np.mean(values), np.std(values), np.min(values),
                     np.max(values)))
    print(f"\nside {side} held-out interleaved-frame VAF")
    print("rank   mean     sd     min     max")
    for row in rows:
        print(f"{row[0]:4d} {row[1]:7.3f} {row[2]:6.3f} "
              f"{row[3]:7.3f} {row[4]:7.3f}")

    # Full-data six-component spatial stability after optimal matching.
    components = []
    for seed in SEEDS:
        _, h = fit(x, x, 6, seed)
        components.append(h)
    ref = components[0]
    similarities = []
    for h in components[1:]:
        sim = cosine_rows(ref, h)
        r, c = linear_sum_assignment(-sim)
        similarities.extend(sim[r, c])
    print(f"side {side} K=6 spatial component cosine after matching: "
          f"mean={np.mean(similarities):.3f}, min={np.min(similarities):.3f}")
    return rows, similarities


def main():
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    data = np.load(fsa.HERE / "bsolve_out.npz", allow_pickle=True)
    time = np.asarray(data["t"], dtype=float)
    payload = {}
    for side in ("r", "l"):
        _, activation, _ = fsa.side_data(data, side)
        x = np.clip(fsa.smooth_matrix(time, activation), 0.0, 1.0)
        rows, similarities = audit_side(x, side)
        payload[side] = {
            "heldout": [dict(rank=r, mean=mean, sd=sd, min=lo, max=hi)
                        for r, mean, sd, lo, hi in rows],
            "k6_spatial_cosine_mean": float(np.mean(similarities)),
            "k6_spatial_cosine_min": float(np.min(similarities)),
        }
    path = fsa.OUT / "fsa_rank_robustness.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
