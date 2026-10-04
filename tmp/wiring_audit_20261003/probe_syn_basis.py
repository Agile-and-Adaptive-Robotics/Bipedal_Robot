"""Audit probe 1: synergy_basis.npz vs fsa_results/fsa_backsolve.npz.

Checks (syn6 rule verification):
  A. W_{side} shape == (n_muscles, 6), names end with _side, all classifiable.
  B. W == gain.T verbatim from the fsa npz (synergy_model.stage_basis layout).
  C. max W < 1.6 (Eq-18 validity bound dE/R = 8/5).
  D. rank(W) per side (rank-6 claim).
  E. pruned runner muscles absent from basis names.
  F. stored-H replay R2/VAF vs fsa target (basis of record consistency).
  G. fsa npz phase keys {side}_pf_phase_mean / {side}_phase_grid exist;
     stance_frac per channel; family assignment the builder computes;
     S5 asymmetry claim (0.58 r vs 0.24 l).
  H. Eq-18 sample: hand-compute g for a few W entries vs analytical_conductance.
"""
import io
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, str(SPINAL))

PRUNE = {"quad_fem_r", "quad_fem_l", "gem_r", "gem_l", "peri_r", "peri_l"}

basis = np.load(SPINAL / "synergy_basis.npz", allow_pickle=True)
fsa = np.load(SPINAL / "fsa_results" / "fsa_backsolve.npz", allow_pickle=True)

print("== A/B: shapes, names, verbatim W ==")
for side in ("r", "l"):
    W = np.asarray(basis[f"W_{side}"], dtype=float)
    names = [str(x) for x in basis[f"muscle_names_{side}"]]
    gain = np.asarray(fsa[f"{side}_gain"], dtype=float)
    fnames = [str(x) for x in fsa[f"{side}_names"]]
    print(f"side {side}: W {W.shape}, {len(names)} names, "
          f"gain {gain.shape}; all end _{side}: "
          f"{all(n.endswith('_' + side) for n in names)}; "
          f"W == gain.T verbatim: {np.array_equal(W, gain.T)}; "
          f"names == fsa names: {names == fnames}")
    import muscle_map as mm
    unclassifiable = [n for n in names if mm.classify(n) is None]
    print(f"  unclassifiable names: {unclassifiable}")
    pruned_here = sorted(set(names) & PRUNE)
    print(f"  pruned muscles present in basis: {pruned_here}")

print("\n== C/D: Eq-18 validity bound + rank ==")
for side in ("r", "l"):
    W = np.asarray(basis[f"W_{side}"], dtype=float)
    print(f"side {side}: W.min={W.min():.6f} W.max={W.max():.6f} "
          f"(bound 1.6), n_nonzero={int((W > 0).sum())}, "
          f"rank={np.linalg.matrix_rank(W)}")

print("\n== F: stored-H replay vs target ==")


def r2c(t, p):
    sse = float(np.sum((t - p) ** 2))
    sst = float(np.sum((t - t.mean()) ** 2))
    return 1.0 - sse / max(sst, 1e-12)


def vafu(t, p):
    return 1.0 - float(np.sum((t - p) ** 2)) / max(float(np.sum(t ** 2)), 1e-12)


for side in ("r", "l"):
    W = np.asarray(basis[f"W_{side}"], dtype=float)
    H = np.asarray(basis[f"H_{side}"], dtype=float)
    tgt = np.asarray(basis[f"target_{side}"], dtype=float)
    rec = H.T @ W.T
    print(f"side {side}: stored-H replay R2={r2c(tgt, rec):.4f} "
          f"VAF={vafu(tgt, rec):.4f}")

print("\n== G: fsa phase keys + family assignment (builder logic) ==")
print("fsa keys:", sorted(str(k) for k in fsa.files))
phase = {}
for side in ("r", "l"):
    try:
        pm = np.asarray(fsa[f"{side}_pf_phase_mean"], dtype=float)
        grid = np.asarray(fsa[f"{side}_phase_grid"], dtype=float)
    except KeyError as e:
        print(f"side {side}: MISSING KEY {e}")
        continue
    st = pm[grid < 50].mean(axis=0)
    sw = pm[grid >= 50].mean(axis=0)
    sf = st / np.maximum(st + sw, 1e-12)
    peak = [float(grid[int(np.argmax(pm[:, k]))]) for k in range(pm.shape[1])]
    phase[side] = (sf, peak)
    print(f"side {side}: stance_frac={np.round(sf, 3).tolist()}")
    print(f"          peak_phase ={np.round(peak, 1).tolist()}")
if all(phase.get(s) for s in ("r", "l")):
    mean_sf = 0.5 * (phase["r"][0] + phase["l"][0])
    fam = ["E" if mean_sf[k] > 0.5 else "F" for k in range(6)]
    print(f"SYMMETRIZED stance_frac={np.round(mean_sf, 3).tolist()}")
    print(f"FAMILIES = {fam}")
    print(f"S5 stance frac r={phase['r'][0][4]:.3f} l={phase['l'][0][4]:.3f} "
          f"(claim: 0.58 r vs 0.24 l)")

print("\n== H: Eq-18 hand check vs analytical_conductance ==")
from fsa_backsolve import analytical_conductance  # noqa: E402
W = np.asarray(basis["W_r"], dtype=float)
flat = W[W > 0]
sample = flat[np.linspace(0, flat.size - 1, 5).astype(int)]
for k in sample:
    g_hand = k * 5.0 * 1.0 / (8.0 - k * 5.0)
    g_fn, valid = analytical_conductance(np.array([k]))
    print(f"k={k:.6f}  g_hand={g_hand:.6f}  g_fn={g_fn[0]:.6f}  "
          f"valid={bool(valid[0])}")
g_all, valid_all = analytical_conductance(W)
print(f"all W entries valid for Eq-18: {bool(valid_all.all())}; "
      f"invalid count={int((~valid_all).sum())}")
gs = g_all[W > 0]
print(f"Eq-18 conductance stats over W>0 (r): n={gs.size} "
      f"min={gs.min():.4f} max={gs.max():.4f} mean={gs.mean():.4f}")

print("\n== extra: df-mass channel per side (F family) ==")
import muscle_map as mm
for side in ("r", "l"):
    names = [str(x) for x in basis[f"muscle_names_{side}"]]
    W = np.asarray(basis[f"W_{side}"], dtype=float)
    df_mass = np.zeros(6)
    for i, act in enumerate(names):
        info = mm.classify(act)
        if info and info.groups[0] == "ankle_df":
            df_mass += W[i]
    print(f"side {side}: ankle_df W mass per channel="
          f"{np.round(df_mass, 3).tolist()}")
