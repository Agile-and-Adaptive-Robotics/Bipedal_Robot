"""Analyze W2L-family chart .txt files (Rhythm Generator, L/R Hip PF, L/R Knee PF)
for flexor half-center bursts.

Usage: _w2l_analyze_charts.py [folder]   (default: Walker_2_Layer_CPG modern folder)

Threshold: -55 mV per the 09-18 notes = -0.055 in chart units (charts store VOLTS;
neurons rest at -60 mV per the asim's RestingPot). Count FLEXOR HC upward crossings
of the threshold (extensor HCs are tonic and not counted). Trailing zero-fill rows
are popped.
"""
import os
import sys
import numpy as np

FOLDER = sys.argv[1] if len(sys.argv) > 1 else \
    r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG"
THR = -0.055

CHARTS = [
    "Rhythm Generator.txt",   # L RG only (no R columns in this chart)
    "L Hip PF.txt",
    "R Hip PF.txt",
    "L Knee PF.txt",
    "R Knee PF.txt",
]

def load_chart(name):
    path = os.path.join(FOLDER, name)
    if not os.path.exists(path):
        return None, None, None
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split("\t")
    data = np.loadtxt(path, skiprows=1)
    datacols = data[:, 2:] if data.shape[1] > 2 else data
    nonzero = np.where(np.any(datacols != 0.0, axis=1))[0]
    n_pop = data.shape[0] - 1 - (nonzero[-1] if nonzero.size else -1)
    if n_pop > 0:
        data = data[: nonzero[-1] + 1]
    return header, data, n_pop

def upward_crossings(t, v, thr):
    return t[np.where((v[:-1] < thr) & (v[1:] >= thr))[0] + 1]

def burst_stats(t, v, thr):
    on = upward_crossings(t, v, thr)
    down = t[np.where((v[:-1] >= thr) & (v[1:] < thr))[0] + 1]
    per = np.diff(on)
    duty = (v >= thr).astype(float).mean() if v.size else np.nan
    return on, down, per, duty

flexor_cols, signals = [], {}
for name in CHARTS:
    header, data, n_pop = load_chart(name)
    if header is None:
        print(f"=== {name}: MISSING ===")
        continue
    t = data[:, 1]
    print(f"=== {name} ===")
    print(f"  rows={data.shape[0]} (popped {n_pop} trailing zero-fill) "
          f"time {t[0]:.3f}..{t[-1]:.3f} s, dt={np.median(np.diff(t))*1000:.1f} ms")
    for j, col in enumerate(header):
        if col in ("TimeSlice", "Time"):
            continue
        v = data[:, j]
        print(f"  {col:22s} min={v.min():8.4f}  max={v.max():8.4f}  mean={v.mean():8.4f}  (chart units V)")
        signals[col] = (t, v)
        if "flx" in col.lower() and not col.lower().endswith("in"):
            flexor_cols.append((name, col))
    for _, col in [c for c in flexor_cols if c[0] == name]:
        t, v = signals[col]
        on, down, per, duty = burst_stats(t, v, THR)
        n = len(on)
        if n >= 2:
            print(f"  --> FLEXOR {col}: onsets={n}  onsets@[{on[0]:.3f}..{on[-1]:.3f}] s")
            print(f"      period median={np.median(per):.4f} s  mean={per.mean():.4f}+/-{per.std(ddof=1):.4f}"
                  f"  min={per.min():.4f}  max={per.max():.4f}")
            print(f"      mean period over span = {(on[-1]-on[0])/(n-1):.4f} s ({(n-1)/(on[-1]-on[0]):.3f} Hz)")
        elif n == 1:
            print(f"  --> FLEXOR {col}: onsets=1 @{on[0]:.3f} s (no period)")
        else:
            print(f"  --> FLEXOR {col}: onsets=0 (no -55 mV upward crossings)")
        print(f"      duty(above {THR} V)={duty:.3f}  downward crossings={len(down)}")
    print()

print("=== L/R antiphase (onset-lag method) ===")
pairs = [("L Hip PF flx", "R Hip PF flx"), ("L Knee PF flx", "R Knee PF flx")]
for lc, rc in pairs:
    if lc not in signals or rc not in signals:
        print(f"  {lc} vs {rc}: chart missing one side - skipped")
        continue
    tl, vl = signals[lc]
    tr, vr = signals[rc]
    onl, _, perl, _ = burst_stats(tl, vl, THR)
    onr, _, perr, _ = burst_stats(tr, vr, THR)
    if len(onl) < 2 or len(onr) < 2:
        print(f"  {lc} vs {rc}: too few onsets (L={len(onl)}, R={len(onr)}) - no rhythm to phase")
        continue
    per = np.median(perl)
    lags = []
    for o in onl:
        d = onr - o
        k = np.argmin(np.abs(d))
        lags.append(d[k])
    lags = np.array(lags)
    fracw = ((lags / per) + 0.5) % 1.0 - 0.5
    print(f"  {lc} vs {rc}: median period L={per:.4f} s, R={np.median(perr):.4f} s")
    print(f"    L->R lag median={np.median(lags)*1000:+.1f} ms = {np.median(fracw):+.3f} cycle "
          f"(n={len(lags)}; |0.5| = antiphase)")

if "L RG flx" in signals and "L Hip PF flx" in signals:
    tg, vg = signals["L RG flx"]
    th, vh = signals["L Hip PF flx"]
    ong, _, perg, _ = burst_stats(tg, vg, THR)
    onh, _, perh, _ = burst_stats(th, vh, THR)
    print("=== Left-side consistency: L RG flx vs L Hip PF flx ===")
    print(f"  onsets RG={len(ong)}, HipPF={len(onh)}; median period RG={np.median(perg):.4f} s, "
          f"HipPF={np.median(perh):.4f} s  (RG chart logs LEFT side only; no R RG columns)")
