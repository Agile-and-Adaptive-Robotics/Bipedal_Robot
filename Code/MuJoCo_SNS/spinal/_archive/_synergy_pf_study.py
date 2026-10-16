"""PF-layer design study (Ben 2026-09-16): NMF synergies on the OpenSim
SO backsolved activations, cycle-phase-folded so each synergy reads as a
candidate PATTERN-FORMATION layer (a PF layer = a set of muscles co-driven
during a phase window). Also prints the CURRENT 4-cell PF mapping and the
hip+knee biarticular membership (primary full + secondary 0.5 weight)."""
import io
import json
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

HERE = Path(__file__).parent
OSIM_DIR = Path(r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim"
                r"\Gait2392_Robotbody\ResultsBSolve")
CYCLE = 1.23      # s, subject01_walk1 cycle (GRF onsets, kine_ref)
NBINS = 20


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
muscle_names = act_names[1:]
print(f"data: {so_act.shape[0]} frames x {so_act.shape[1]} muscles, "
      f"t {t[0]:.2f}-{t[-1]:.2f} s")

active_mask = so_act.max(axis=0) > 0.05
active_names = [muscle_names[i] for i in range(len(muscle_names))
                if active_mask[i]]
X = so_act[:, active_mask]
print(f"active muscles (peak > 0.05): {X.shape[1]}")

from sklearn.decomposition import NMF

print(f"\n{'n':>3s} {'recon_err':>10s} {'explained':>10s} "
      f"{'BIC-ish':>9s}")
results = {}
T, M = X.shape
for n_comp in range(2, 9):
    nmf = NMF(n_components=n_comp, init="nndsvda", max_iter=500,
              random_state=42)
    W = nmf.fit_transform(X)
    H = nmf.components_
    Xr = W @ H
    err2 = float(np.mean((X - Xr) ** 2))
    exp = 1.0 - err2 / max(float(np.var(X)), 1e-12)
    k = n_comp * (T + M)
    bic = T * M * np.log(max(err2, 1e-12)) + k * np.log(T * M)
    results[n_comp] = (W, H, exp)
    print(f"{n_comp:3d} {np.sqrt(err2):10.4f} {exp:10.3f} {bic:9.0f}")


def phase_profile(w):
    """Cycle-fold one synergy's coefficient; return NBINS mean profile."""
    ph = np.mod(t, CYCLE) / CYCLE
    prof = np.zeros(NBINS)
    for b in range(NBINS):
        m = (ph >= b / NBINS) & (ph < (b + 1) / NBINS)
        prof[b] = w[m].mean() if m.any() else 0.0
    return prof / max(prof.max(), 1e-9)


def window(prof, thr=0.5):
    on = np.where(prof >= thr)[0]
    if len(on) == 0:
        return "-"
    return f"{on[0] * 100 // NBINS}-{(on[-1] + 1) * 100 // NBINS}%cycle"


for n_comp in (3, 4, 5, 6):
    W, H, exp = results[n_comp]
    print(f"\n===== {n_comp} synergies (explained var {exp:.3f}) =====")
    for c in range(n_comp):
        weights = H[c]
        order = np.argsort(weights)[::-1]
        top = [(active_names[i], weights[i]) for i in order
               if weights[i] > 0.2 * weights[order[0]]]
        prof = phase_profile(W[:, c])
        bar = "".join("#" if v >= 0.75 else ("+" if v >= 0.5 else ".")
                      for v in prof)
        names_str = ", ".join(f"{n}" for n, _ in top[:9])
        print(f"  syn{c}: peak@{t[int(np.argmax(W[:, c]))]:.2f}s "
              f"win>{window(prof)}  [{bar}] 5-95%cycle")
        print(f"       {names_str}")

# ---- current 4-cell PF mapping (fitted refit weights), biarticular focus
fit = json.loads((HERE / "fitted_walk_params.json").read_text("utf-8"))
WPF = fit["W_PF_MN"]
WPOST = fit["W_POSTURE"]
import muscle_map
print("\n===== CURRENT 4-cell PF columns (fitted refit; primary + 0.5x "
      "secondary) =====")
grp_of = {}
for base, groups in muscle_map._GROUPS_BY_NAME.items():
    grp_of[base] = groups
print(f"{'muscle':22s} {'groups':16s} {'E1':>6s} {'E2':>6s} {'F1':>6s} "
      f"{'F2':>6s} {'POST':>6s}")
for base in sorted(grp_of):
    prim, sec = grp_of[base][0], grp_of[base][1:]
    row = []
    for ph in ("E1", "E2", "F1", "F2"):
        w = WPF[ph].get(prim, 0.0) + 0.5 * sum(WPF[ph].get(s, 0.0)
                                               for s in sec)
        row.append(w)
    post = WPOST.get(prim, 0.0) + 0.5 * sum(WPOST.get(s, 0.0) for s in sec)
    if max(row) > 0.02 or post > 0.02:
        bi = "+".join(g[:7] for g in grp_of[base])
        print(f"{base:22s} {bi:16s} " +
              " ".join(f"{w:6.3f}" for w in row) + f" {post:6.3f}")
