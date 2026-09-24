"""Twin check: are ext/flx chart columns in the modern W2L run identical traces?"""
import os
import numpy as np

FOLDER = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG"

def load(name):
    path = os.path.join(FOLDER, name)
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split("\t")
    data = np.loadtxt(path, skiprows=1)
    datacols = data[:, 2:] if data.shape[1] > 2 else data
    nz = np.where(np.any(datacols != 0.0, axis=1))[0]
    npop = data.shape[0] - 1 - (nz[-1] if nz.size else -1)
    if npop > 0:
        data = data[: nz[-1] + 1]
    return header, data

for name in ["Rhythm Generator.txt", "L Hip PF.txt", "L Knee PF.txt"]:
    h, d = load(name)
    print(f"=== {name} ===")
    cols = [c for c in h if c not in ("TimeSlice", "Time")]
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            a, b = d[:, h.index(cols[i])], d[:, h.index(cols[j])]
            print(f"  {cols[i]} vs {cols[j]}: max|diff|={np.max(np.abs(a-b)):.2e} V "
                  f"corr={np.corrcoef(a,b)[0,1]:+.6f}")
