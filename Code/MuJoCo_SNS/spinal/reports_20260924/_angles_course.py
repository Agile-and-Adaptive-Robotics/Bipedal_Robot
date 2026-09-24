"""Joint-angle time course check: repeated stepping vs one-time settling.
Usage: _angles_course.py <folder-with-L Angles.txt> [<folder2> ...]
"""
import os
import sys
import numpy as np

THR_NMV = -0.055

def load(folder, name):
    path = os.path.join(folder, name)
    if not os.path.exists(path):
        return None, None
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split("\t")
    data = np.loadtxt(path, skiprows=1)
    datacols = data[:, 2:] if data.shape[1] > 2 else data
    nz = np.where(np.any(datacols != 0.0, axis=1))[0]
    npop = data.shape[0] - 1 - (nz[-1] if nz.size else -1)
    if npop > 0:
        data = data[: nz[-1] + 1]
    return header, data

for folder in sys.argv[1:]:
    print(f"===== {folder} =====")
    header, data = load(folder, "L Angles.txt")
    if header is None:
        print("  no L Angles.txt")
        continue
    t = data[:, 1]
    for j, col in enumerate(header[2:], start=2):
        v = data[:, j]
        # count positive-going peaks: crossings of (min + 0.5*range) upward
        half = v.min() + 0.5 * (v.max() - v.min())
        cross = np.where((v[:-1] < half) & (v[1:] >= half))[0]
        print(f"  {col:20s} min={np.degrees(v.min()):8.2f} max={np.degrees(v.max()):8.2f} deg "
              f"half-range upcrossings={len(cross)}")
        if len(cross):
            periods = np.diff(t[cross])
            if len(periods):
                print(f"      upcrossing times (deg>half-range): "
                      f"{np.round(t[cross], 2).tolist()}")
                print(f"      intervals median={np.median(periods):.3f} s" if len(periods) else "")
    # short time course every 0.5 s for hip
    ih = header.index("hip_L") if "hip_L" in header else 2
    ik = header.index("knee_L") if "knee_L" in header else 3
    print("  t(s)   hip_L(deg)  knee_L(deg)")
    for tt in np.arange(0, t[-1] + 0.001, 0.5):
        k = np.searchsorted(t, tt)
        k = min(k, len(t) - 1)
        print(f"  {t[k]:5.2f}  {np.degrees(data[k, ih]):9.2f}  {np.degrees(data[k, ik]):9.2f}")
