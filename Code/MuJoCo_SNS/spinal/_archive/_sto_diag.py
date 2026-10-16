"""Diagnose the TendonForce.sto column layout."""
import io
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
p = Path(r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody"
         r"\ResultsBSolve\zz_bsolve_MuscleAnalysis_TendonForce.sto")
lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
end = next(i for i, ln in enumerate(lines) if ln.strip() == "endheader")
names = lines[end + 1].split()
print("n columns:", len(names))
print("first 8 names:", names[:8])
print("time col?", names[0])
rows = []
for ln in lines[end + 2:]:
    s = ln.split()
    if s:
        rows.append([float(x) for x in s])
data = np.array(rows)
print("data shape:", data.shape)
peaks = data.max(axis=0)
nz = [(names[i], peaks[i]) for i in range(1, len(names)) if peaks[i] > 1.0]
print(f"columns with peak > 1 N: {len(nz)}")
for nm, pk in nz[:15]:
    print(f"  {nm:16s} peak {pk:.0f}")
for target in ("soleus_r", "vas_lat_r", "rect_fem_r"):
    if target in names:
        i = names.index(target)
        print(f"{target}: header idx {i}, vals idx {i-1}, peak "
              f"{peaks[i]:.1f}, peak(vals idx-1) {peaks[i-1]:.1f}")
