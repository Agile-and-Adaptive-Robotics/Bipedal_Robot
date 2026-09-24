"""Analyze the 2023 W2L asim (walk new new tester added 2 axis) DataTool_7.txt chart.

Columns: contact neurons + hip 'middle' neurons (MembraneVoltage, V), height/distance (m).
These are SPIKING neurons (Equil -50 mV class per AGENTS.md), so rhythm is detected as
burst sequences of spike threshold crossings, plus graded analysis as fallback.
"""
import os
import numpy as np

PATH = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\DataTool_7.txt"

with open(PATH, "r", encoding="utf-8", errors="replace") as f:
    header = f.readline().strip().split("\t")
data = np.loadtxt(PATH, skiprows=1)
t = data[:, 1]
print(f"rows={data.shape[0]}  time {t[0]:.3f}..{t[-1]:.3f} s  dt={np.median(np.diff(t))*1000:.1f} ms")

sig = {}
for j, col in enumerate(header):
    if col in ("TimeSlice", "Time"):
        continue
    v = data[:, j]
    print(f"{col:24s} min={v.min():9.4f}  max={v.max():9.4f}  mean={v.mean():9.4f}")
    sig[col] = v

# Spike detection: spiking neurons in AnimatLab fire at threshold ~ 0 mV on the
# membrane trace (rest -60 mV, peak ~ +40..50 mV). Count UPWARD crossings of 0 V.
def spikes(v, thr=0.0):
    return t[np.where((v[:-1] < thr) & (v[1:] >= thr))[0] + 1]

def bursts(st, gap=0.05):
    """Group spike times into bursts split when ISI > gap."""
    if st.size == 0:
        return []
    ids = np.where(np.diff(st) > gap)[0]
    return np.split(st, ids + 1)

print("\n=== spike/burst summary (threshold 0 V) ===")
for col in header[2:]:
    if col in ("TimeSlice", "Time", "height", "distance"):
        continue
    st = spikes(sig[col])
    bs = bursts(st)
    if len(bs) >= 3:
        onsets = np.array([b[0] for b in bs])
        per = np.diff(onsets)
        print(f"{col:24s} spikes={len(st):5d} bursts={len(bs):3d} "
              f"burst onset span {onsets[0]:.2f}..{onsets[-1]:.2f} s "
              f"period median={np.median(per):.3f} s mean={per.mean():.3f}+/-{per.std(ddof=1):.3f} s "
              f"({(len(bs)-1)/(onsets[-1]-onsets[0]):.2f} Hz) "
              f"burst dur median={np.median([b[-1]-b[0] for b in bs])*1000:.0f} ms")
    elif len(bs) > 0:
        print(f"{col:24s} spikes={len(st):5d} bursts={len(bs):3d}  (too few for period)")
    else:
        print(f"{col:24s} spikes=0")

# L/R antiphase on the hip middle neurons: cross-correlation of the full traces
print("\n=== L/R hip middle antiphase (full-trace correlation sweep) ===")
rl, ll = sig.get("R_hip middle"), sig.get("L_hip middle")
if rl is not None and ll is not None:
    dt = float(np.median(np.diff(t)))
    per_est = t[-1] / max(1, len(bursts(spikes(rl))))
    maxlag = int(round(min(2.0, per_est) / dt))
    n = len(t) - maxlag - 1
    best = (None, -2.0, None)
    for lag in range(-maxlag, maxlag + 1):
        a = rl[maxlag:maxlag + n]
        b = ll[maxlag + lag:maxlag + lag + n]
        c = np.corrcoef(a, b)[0, 1]
        if c > best[1]:
            best = (lag * dt, c, lag)
    print(f"  corr peak lag={best[0]*1000:+.1f} ms (r={best[1]:+.3f}); "
          f"anticorrelation minimum for reference:")
    worst = (None, 2.0)
    for lag in range(-maxlag, maxlag + 1):
        a = rl[maxlag:maxlag + n]
        b = ll[maxlag + lag:maxlag + lag + n]
        c = np.corrcoef(a, b)[0, 1]
        if c < worst[1]:
            worst = (lag * dt, c)
    print(f"  corr min lag={worst[0]*1000:+.1f} ms (r={worst[1]:+.3f})")

# Body motion: height/distance range (air-walking = legs swing under a hanging body)
h, d = sig["height"], sig["distance"]
print(f"\nheight {h.min():.4f}..{h.max():.4f} m (range {h.max()-h.min()*1:.4f} m)")
print(f"distance {d.min():.4f}..{d.max():.4f} m (drift {d.max()-d.min():.4f} m)")
