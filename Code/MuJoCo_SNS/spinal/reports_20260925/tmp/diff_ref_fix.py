import json
import numpy as np
a = json.load(open(r"reports_20260925/tmp/ref_before_fix.json"))
b = json.load(open(r"reports_20260925/tmp/ref_after_fix.json"))
print("scalar diffs (should be LEFT-side only):")
for k in sorted(set(a) | set(b)):
    if isinstance(a.get(k), dict) or isinstance(b.get(k), dict):
        continue
    va, vb = a.get(k), b.get(k)
    if va != vb:
        print(f"  {k}: {va} -> {vb}")
print("cycle diffs:")
for side in ("r", "l"):
    for j in ("hip", "knee", "ankle"):
        xa, xb = np.array(a[side][j]), np.array(b[side][j])
        d = np.abs(xa - xb)
        print(f"  ref['{side}'][{j}]: max|diff| {d.max():.2f} deg, "
              f"first-40% mean|diff| {d[:40].mean():.2f}, "
              f"last-60% mean|diff| {d[40:].mean():.2f}")
lead = np.array(b["l"]["knee"])[:39]
print(f"left knee first-39-gridpt spread after fix: "
      f"{lead.max()-lead.min():.2f} deg (was ~0 = frozen)")
