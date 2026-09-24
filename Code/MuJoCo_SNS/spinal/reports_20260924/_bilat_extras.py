"""BilateralRG extras: L RG flx vs R RG flx antiphase + joint-angle movement check."""
import os
import sys
import numpy as np

FOLDER = sys.argv[1]
THR = -0.055

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

def onsets(t, v, thr):
    return t[np.where((v[:-1] < thr) & (v[1:] >= thr))[0] + 1]

h, d = load("Rhythm Generator.txt")
t = d[:, 1]
il = h.index("L RG flx")
ir = h.index("R RG flx")
onl = onsets(t, d[:, il], THR)
onr = onsets(t, d[:, ir], THR)
print(f"L RG flx onsets n={len(onl)}: {np.round(onl, 3).tolist()}")
print(f"R RG flx onsets n={len(onr)}: {np.round(onr, 3).tolist()}")
if len(onl) >= 2 and len(onr) >= 2:
    per = np.median(np.diff(onl))
    lags = []
    for o in onl:
        k = np.argmin(np.abs(onr - o))
        lags.append(((onr[k] - o) / per + 0.5) % 1.0 - 0.5)
    lags = np.array(lags)
    print(f"period L={per:.4f} s R={np.median(np.diff(onr)):.4f} s")
    print(f"L->R onset lag median={np.median(lags):+.3f} cycle "
          f"({np.median(lags)*per*1000:+.0f} ms)  (|0.5| = antiphase)")

# joint movement
p = os.path.join(FOLDER, "L Angles.txt")
if os.path.exists(p):
    h2, d2 = load("L Angles.txt")
    print(f"\nL Angles.txt columns: {h2[2:]}")
    for j, col in enumerate(h2[2:], start=2):
        v = d2[:, j]
        print(f"  {col:20s} min={v.min():9.4f} max={v.max():9.4f} range={v.max()-v.min():.4f}")
