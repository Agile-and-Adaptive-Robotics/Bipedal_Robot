"""Analyze the 2023 W2L asim DataTool_7.txt - part 2: trailing zeros, body motion,
L/R antiphase (burst-onset lag + fixed correlation sweep)."""
import numpy as np

PATH = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\DataTool_7.txt"
with open(PATH, "r", encoding="utf-8", errors="replace") as f:
    header = f.readline().strip().split("\t")
data = np.loadtxt(PATH, skiprows=1)
t = data[:, 1]

# trailing zero-fill check
neuro_cols = list(range(2, 8))  # the six MembraneVoltage columns
allzero_rows = np.where(np.all(data[:, neuro_cols] == 0.0, axis=1))[0]
print(f"rows where ALL six neuron cols == 0 exactly: {allzero_rows.size}")
if allzero_rows.size:
    print(f"  first at row {allzero_rows[0]} (t={t[allzero_rows[0]]:.3f} s), last at row "
          f"{allzero_rows[-1]} (t={t[allzero_rows[-1]]:.3f} s)")

h, d = data[:, 8], data[:, 9]
print("\nheight/distance sampled every 1 s (before any zero-fill tail):")
for tt in np.arange(0, 10.01, 1.0):
    k = np.searchsorted(t, tt)
    k = min(k, len(t) - 1)
    print(f"  t={t[k]:5.2f}  height={h[k]:7.4f} m  distance={d[k]:8.4f} m")

def spikes(v, thr=0.0, tt=t):
    return tt[np.where((v[:-1] < thr) & (v[1:] >= thr))[0] + 1]

def bursts(st, gap=0.05):
    if st.size == 0:
        return []
    ids = np.where(np.diff(st) > gap)[0]
    return np.split(st, ids + 1)

sr = spikes(data[:, 4])   # R_hip middle
sl = spikes(data[:, 6])   # L_hip middle
br, bl = bursts(sr), bursts(sl)
onr = np.array([b[0] for b in br])
onl = np.array([b[0] for b in bl])
print(f"\nR_hip middle bursts={len(br)} onsets={np.round(onr, 2).tolist()}")
print(f"L_hip middle bursts={len(bl)} onsets={np.round(onl, 2).tolist()}")
per = np.median(np.diff(onr))
print(f"median period (R onsets) = {per:.3f} s ({1/per:.3f} Hz)")

# antiphase via onset lags: for each L onset, nearest R onset, wrapped to +/-0.5 cycle
lags = []
for o in onl:
    k = np.argmin(np.abs(onr - o))
    lagf = ((onr[k] - o) / per + 0.5) % 1.0 - 0.5
    lags.append(lagf)
lags = np.array(lags)
print(f"L->R onset lag: median={np.median(lags):+.3f} cycle ({np.median(lags)*per*1000:+.0f} ms), "
      f"n={len(lags)} (|0.5| = antiphase, 0 = in-phase)")

# fixed correlation sweep
dt = float(np.median(np.diff(t)))
maxlag = int(round(2.0 / dt))
a = data[:, 4][maxlag:len(t) - maxlag]
r_best, r_worst = (-2.0, None), (2.0, None)
for lag in range(-maxlag, maxlag + 1):
    b = data[:, 6][maxlag + lag:len(t) - maxlag + lag]
    c = np.corrcoef(a, b)[0, 1]
    if c > r_best[0]:
        r_best = (c, lag * dt)
    if c < r_worst[0]:
        r_worst = (c, lag * dt)
print(f"full-trace corr: peak r={r_best[0]:+.3f} @ {r_best[1]*1000:+.0f} ms; "
      f"min r={r_worst[0]:+.3f} @ {r_worst[1]*1000:+.0f} ms")
