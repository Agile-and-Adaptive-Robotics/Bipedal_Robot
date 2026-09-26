"""Inspect synergy_basis.npz + fsa_results + muscle_map classify (syn6 prep)."""
import numpy as np
import json
from pathlib import Path

HERE = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")

print("== synergy_basis.npz ==")
d = np.load(HERE / "synergy_basis.npz", allow_pickle=True)
for k in d.files:
    a = d[k]
    print(f"  {k}: shape={a.shape} dtype={a.dtype}")
names_r = [str(x) for x in d["muscle_names_r"]] if "muscle_names_r" in d.files else None
for cand in ("muscle_names_r", "names_r", "muscles_r", "names"):
    if cand in d.files:
        print(f"  names key={cand}: {list(d[cand])[:8]} ...")
# print W_r stats
for wk in ("W_r", "W_l"):
    if wk in d.files:
        W = d[wk]
        print(f"  {wk}: min={W.min():.4f} max={W.max():.4f} mean={W.mean():.4f}")
        print(f"  {wk} col max: {np.round(W.max(axis=0), 3)}")
        print(f"  {wk} col min: {np.round(W.min(axis=0), 3)}")
# name keys?
for k in d.files:
    if "name" in k.lower():
        arr = d[k]
        print(f"  {k} -> first 6: {[str(x) for x in arr.ravel()[:6]]}")
        print(f"  {k} -> len {arr.size}")

print()
print("== fsa_results ==")
fsa = HERE / "fsa_results"
for p in sorted(fsa.glob("*")):
    print("  ", p.name, p.stat().st_size if p.is_file() else "<dir>")
if (fsa / "fsa_backsolve.npz").exists():
    f = np.load(fsa / "fsa_backsolve.npz", allow_pickle=True)
    for k in f.files:
        a = f[k]
        if a.dtype.kind in "US" or a.size <= 20:
            print(f"  fsa {k}: {a.shape} {a.dtype} {[str(x) for x in np.ravel(a)[:8]]}")
        else:
            print(f"  fsa {k}: shape={a.shape} dtype={a.dtype} min={np.nanmin(a):.4g} max={np.nanmax(a):.4g} nan%={100*np.mean(np.isnan(a.astype(float))):.1f}")
