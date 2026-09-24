"""Deep dive: modern W2L run - are the neural traces oscillating sub-threshold and
does their frequency match the joint oscillation?

- post-transient (t>=0.5 s) max of each neural trace
- dominant frequency (FFT periodogram, 0.5-10 Hz band) of neural traces and joint angles
- cross-correlation (per-frequency match) between neural trace and knee/hip angle
"""
import os
import numpy as np

FOLDER = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG"
T_TRANS = 0.5  # s, drop initialization transient

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

def dom_freq(t, v, fmin=0.5, fmax=10.0):
    v = v - v.mean()
    dt = float(np.median(np.diff(t)))
    n = len(v)
    freqs = np.fft.rfftfreq(n, dt)
    amp = np.abs(np.fft.rfft(v))
    band = (freqs >= fmin) & (freqs <= fmax)
    k = np.argmax(amp[band])
    f = freqs[band][k]
    # significance: peak amplitude vs rms of band
    rms = np.sqrt(np.mean(amp[band] ** 2))
    return f, amp[band][k] / (rms if rms > 0 else 1)

neural = {}
for name, cols in [
    ("Rhythm Generator.txt", ["L RG ext", "L RG flx"]),
    ("L Hip PF.txt", ["L Hip PF ext", "L Hip PF flx"]),
    ("L Knee PF.txt", ["L Knee PF flx", "L Knee PF ext"]),
    ("R Hip PF.txt", ["R Hip PF flx"]),
    ("R Knee PF.txt", ["R Knee PF flx"]),
]:
    h, d = load(name)
    t = d[:, 1]
    for c in cols:
        j = h.index(c)
        neural[c] = (t, d[:, j])

ha, da = load("L Angles.txt")
ta = da[:, 1]
angles = {c: (ta, da[:, ha.index(c)]) for c in ("hip_L", "knee_L", "ankle_L")}

print(f"post-transient window t>={T_TRANS} s")
print(f"{'trace':18s} {'post max (mV)':>13s} {'post min (mV)':>13s} {'swing (mV)':>10s} {'dom f (Hz)':>10s} {'peak/rms':>8s}")
for c, (t, v) in {**neural, **angles}.items():
    m = t >= T_TRANS
    vv = v[m]
    if c in angles:
        print(f"{c:18s} {np.degrees(vv.max()):13.2f} {np.degrees(vv.min()):13.2f} "
              f"{np.degrees(vv.max()-vv.min()):10.2f}", end="")
        unit = "deg"
    else:
        print(f"{c:18s} {vv.max()*1000:13.2f} {vv.min()*1000:13.2f} "
              f"{(vv.max()-vv.min())*1000:10.2f}", end="")
        unit = "mV"
    f, s = dom_freq(t[m], vv)
    print(f" {f:10.3f} {s:8.1f}  ({unit})")

# cross-correlation: knee angle vs L Knee PF flx; hip angle vs L RG ext (drive side)
print("\ncross-correlation peak (drive->angle interpretation is loose; frequency match is the point):")
for ac, nc in [("knee_L", "L Knee PF flx"), ("hip_L", "L RG ext"), ("hip_L", "L Hip PF ext")]:
    (t1, v1), (t2, v2) = angles[ac], neural[nc]
    m1 = t1 >= T_TRANS
    a = v1[m1] - v1[m1].mean()
    # neural sampled on same grid; angles chart may differ in rows - align by time grid
    m2 = np.searchsorted(t2, t1[m1])
    if m2[-1] >= len(t2):
        m2 = m2[:-1]; a = a[:-1]
    b = v2[m2] - v2[m2].mean()
    n = len(a)
    corr = np.correlate(a, b, mode="full")
    lags = np.arange(-n + 1, n)
    dt = float(np.median(np.diff(t1)))
    k = np.argmax(np.abs(corr))
    print(f"  {nc} vs {ac}: peak |corr| at lag {lags[k]*dt*1000:+.1f} ms "
          f"(r={corr[k]/ (np.linalg.norm(a)*np.linalg.norm(b)):+.3f})")
