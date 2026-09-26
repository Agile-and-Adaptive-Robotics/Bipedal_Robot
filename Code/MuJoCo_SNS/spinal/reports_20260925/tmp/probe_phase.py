import numpy as np
from pathlib import Path
HERE = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
f = np.load(HERE / "fsa_results" / "fsa_backsolve.npz", allow_pickle=True)
for sd in ("r", "l"):
    pm = np.asarray(f[f"{sd}_pf_phase_mean"], float)
    grid = np.asarray(f[f"{sd}_phase_grid"], float)
    print(f"side {sd}: grid {grid[0]:.0f}..{grid[-1]:.0f} n={grid.size}")
    for k in range(6):
        col = pm[:, k]
        print(f"  S{k+1}: argmax={grid[int(np.argmax(col))]:.1f}%  "
              f"val@0/25/50/75/100="
              f"{[round(float(col[int(np.argmin(np.abs(grid-g)))]), 3) for g in (0,25,50,75,100)]}")
    # symmetrized stance fraction
for sd in ("r", "l"):
    pm = np.asarray(f[f"{sd}_pf_phase_mean"], float)
    grid = np.asarray(f[f"{sd}_phase_grid"], float)
    sf = pm[grid < 50].mean(axis=0) / np.maximum(
        pm[grid < 50].mean(axis=0) + pm[grid >= 50].mean(axis=0), 1e-12)
    print(f"side {sd} stance_frac: {np.round(sf, 3).tolist()}")
sf = {}
for sd in ("r", "l"):
    pm = np.asarray(f[f"{sd}_pf_phase_mean"], float)
    grid = np.asarray(f[f"{sd}_phase_grid"], float)
    sf[sd] = pm[grid < 50].mean(axis=0) / np.maximum(
        pm[grid < 50].mean(axis=0) + pm[grid >= 50].mean(axis=0), 1e-12)
mean_sf = (sf["r"] + sf["l"]) / 2
print("mean stance_frac:", np.round(mean_sf, 3).tolist(),
      "-> families:", ["E" if v > 0.5 else "F" for v in mean_sf])

