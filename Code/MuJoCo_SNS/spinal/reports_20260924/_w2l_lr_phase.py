"""Modern W2L: L vs R phase from R/L Hip PF flexor traces (sub-threshold rhythm,
so use full-trace correlation sweep + low-threshold crossing times)."""
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

hL, dL = load("L Hip PF.txt")
hR, dR = load("R Hip PF.txt")
t = dL[:, 1]
vl = dL[:, hL.index("L Hip PF flx")]
vr = dR[:, hR.index("R Hip PF flx")]
m = t >= 0.5  # post-transient
a = vl[m] - vl[m].mean()
b = vr[m] - vr[m].mean()
n = len(a)
corr = np.correlate(a, b, mode="full")
lags = np.arange(-n + 1, n)
dt = float(np.median(np.diff(t)))
kmin = np.argmin(corr)   # most negative = antiphase lag
kmax = np.argmax(corr)
norm = np.linalg.norm(a) * np.linalg.norm(b)
per = 1 / 2.222
print(f"L/R Hip PF flx post-transient: corr most-negative at lag {lags[kmin]*dt*1000:+.0f} ms "
      f"(r={corr[kmin]/norm:+.3f}) = {lags[kmin]*dt/per:+.3f} cycle")
print(f"corr most-positive at lag {lags[kmax]*dt*1000:+.0f} ms (r={corr[kmax]/norm:+.3f})")

# crossing times at a mid-swing level (-58.5 mV = midway between -61 and -56)
THR2 = -0.0585
onl = t[m][np.where((vl[m][:-1] < THR2) & (vl[m][1:] >= THR2))[0] + 1]
onr = t[m][np.where((vr[m][:-1] < THR2) & (vr[m][1:] >= THR2))[0] + 1]
print(f"mid-swing (-58.5 mV) upcrossings: L n={len(onl)}, R n={len(onr)}")
if len(onl) > 2 and len(onr) > 2:
    lags2 = []
    for o in onl:
        kk = np.argmin(np.abs(onr - o))
        lags2.append(((onr[kk] - o) / per + 0.5) % 1.0 - 0.5)
    lags2 = np.array(lags2)
    print(f"L->R mid-swing lag median={np.median(lags2):+.3f} cycle "
          f"({np.median(lags2)*per*1000:+.0f} ms), n={len(lags2)} (|0.5|=antiphase)")
    print(f"  median interval L={np.median(np.diff(onl)):.4f} s, R={np.median(np.diff(onr)):.4f} s")
