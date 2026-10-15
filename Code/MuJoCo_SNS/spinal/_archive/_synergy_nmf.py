"""NMF synergy analysis on the backsolved muscle activations.
Determines the natural number of synergy groups and their muscle
compositions, to inform the PF layer architecture.
"""
import io
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

OSIM_DIR = Path(r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim"
                r"\Gait2392_Robotbody\ResultsBSolve")


def read_sto(path):
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    end = next(i for i, ln in enumerate(lines) if ln.strip() == "endheader")
    names = lines[end + 1].split()
    rows = []
    for ln in lines[end + 2:]:
        s = ln.split()
        if s:
            rows.append([float(x) for x in s])
    data = np.array(rows)
    return data[:, 0], names, data[:, 1:]


t, act_names, so_act = read_sto(OSIM_DIR / "zz_bsolve_StaticOptimization_activation.sto")
# first name is 'time', vals exclude it
muscle_names = act_names[1:]
print(f"data: {so_act.shape[0]} frames x {so_act.shape[1]} muscles")
print(f"time range: {t[0]:.2f} to {t[-1]:.2f} s")

# filter: remove near-silent muscles (max activation < 0.05)
active_mask = so_act.max(axis=0) > 0.05
active_names = [muscle_names[i] for i in range(len(muscle_names))
                if active_mask[i]]
X = so_act[:, active_mask]
print(f"active muscles (peak > 0.05): {X.shape[1]} of {len(act_names)}")

# NMF
from sklearn.decomposition import NMF

print(f"\n{'n':>3s} {'recon_err':>10s} {'explained':>10s}")
results = {}
for n_comp in range(2, 9):
    nmf = NMF(n_components=n_comp, init="nndsvda", max_iter=500,
              random_state=42)
    W = nmf.fit_transform(X)   # [T, n_comp] activation of each synergy
    H = nmf.components_         # [n_comp, n_muscles] muscle weights
    Xr = W @ H
    err = float(np.sqrt(np.mean((X - Xr) ** 2)))
    tot_var = float(np.var(X))
    exp = 1.0 - err ** 2 / max(tot_var, 1e-12)
    results[n_comp] = (W, H, err, exp)
    print(f"{n_comp:3d} {err:10.4f} {exp:10.3f}")

# detailed composition at the best few n values
for n_comp in (3, 4, 5, 6):
    W, H, err, exp = results[n_comp]
    print(f"\n===== {n_comp} synergies (explained var {exp:.3f}) =====")
    for c in range(n_comp):
        weights = H[c]
        # sort muscles by weight, show top contributors
        order = np.argsort(weights)[::-1]
        top = [(active_names[i], weights[i]) for i in order
               if weights[i] > 0.15 * weights[order[0]]]
        names_str = ", ".join(f"{n}({w:.2f})" for n, w in top[:10])
        duty = float(np.mean(W[:, c] > 0.3 * max(W[:, c].max(), 1e-9)))
        peak_t = t[int(np.argmax(W[:, c]))]
        print(f"  synergy {c}: duty {duty:.2f}, peak t={peak_t:.2f}s")
        print(f"    muscles: {names_str}")
