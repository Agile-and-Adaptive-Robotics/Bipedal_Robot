"""Rhythm-metric analyzer for AnimatLab chart .txt byproducts.

Usage:  python analyze.py <asim> [--mode graded|spiking] [--json out.json]

Reads the asim to resolve each chart's columns, then analyzes the tab-separated
chart .txt files written next to the asim by AnimatSimulator.

Metrics:
  - per MembraneVoltage column: mean/min/max mV (graded) or spike count/rate (spiking)
  - bursts: graded = up-crossings of -60 mV (merged with 50 ms gap); spiking =
    envelope of -40 mV crossings in 10 ms bins, bursts where rate > 20 Hz
  - period = median inter-onset interval; duty = burst length / period
  - antiphase lag: cross-correlation of the two given envelopes
  - JointRotation columns: min/max/range in degrees
Chart .txt values are VOLTS (mV = *1000); rotations are radians.
"""
import io
import sys
import re
import json
import os
import math
import xml.etree.ElementTree as ET

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


def txt(e, t):
    x = e.find(t)
    return x.text if x is not None else None


def read_chart(path):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    # strip trailing all-zero-fill rows (skill trap: they add fake crossings/peaks)
    while rows and all(x == 0.0 for x in rows[-1][2:]):
        rows.pop()
    return header, rows


def upcrossings(t, v, thr):
    """Times of rising crossings of thr."""
    out = []
    above = v[0] > thr
    for i in range(1, len(v)):
        a = v[i] > thr
        if a and not above:
            # linear interp crossing time
            out.append(t[i - 1] + (thr - v[i - 1]) * (t[i] - t[i - 1]) / (v[i] - v[i - 1]))
        above = a
    return out


def merge_events(times, gap):
    out = []
    for x in times:
        if out and x - out[-1] < gap:
            continue
        out.append(x)
    return out


def bursts_from_onsets(onsets, t_end, min_len=0.02):
    """Burst = onset to next onset (or t_end)."""
    out = []
    for i, s in enumerate(onsets):
        e = onsets[i + 1] if i + 1 < len(onsets) else t_end
        if e - s >= min_len:
            out.append((s, e))
    return out


def spike_times(t, v, thr=-0.040, refract=0.001):
    cr = upcrossings(t, v, thr)
    return merge_events(cr, refract)


def envelope(t, events, width=0.01):
    """Events/sec in fixed bins; returns bin centers + rate."""
    if not t:
        return [], []
    t0, t1 = t[0], t[-1]
    nb = max(1, int((t1 - t0) / width))
    rate = [0.0] * nb
    for x in events:
        b = int((x - t0) / width)
        if 0 <= b < nb:
            rate[b] += 1.0 / width
    centers = [t0 + (i + 0.5) * width for i in range(nb)]
    return centers, rate


def xcorr_lag(t1, y1, t2, y2, max_lag):
    """Lag (s) maximizing cross-correlation of two equal-grid envelopes (y2 shifted)."""
    n = len(y1)
    if n == 0 or len(y2) != n:
        return None, None
    m1 = sum(y1) / n
    m2 = sum(y2) / n
    dt = t1[1] - t1[0] if n > 1 else 1.0
    best_lag, best_r = 0, -2
    for k in range(-max_lag, max_lag + 1):
        sxy = sxx = syy = 0.0
        for i in range(n):
            j = i + k
            if 0 <= j < n:
                sxy += (y1[i] - m1) * (y2[j] - m2)
                sxx += (y1[i] - m1) ** 2
                syy += (y2[j] - m2) ** 2
        if sxx > 0 and syy > 0:
            r = sxy / math.sqrt(sxx * syy)
            if r > best_r:
                best_r, best_lag = r, k * dt
    return best_lag, best_r


def analyze_col(colname, t, v, mode):
    d = {"name": colname}
    mv = [x * 1000 for x in v]
    d["mean_mV"] = round(sum(mv) / len(mv), 2)
    d["min_mV"] = round(min(mv), 2)
    d["max_mV"] = round(max(mv), 2)
    if mode == "graded":
        onsets = merge_events(upcrossings(t, v, -0.060), 0.05)
        d["crossings_-60mV"] = len(onsets)
    else:
        st = spike_times(t, v)
        d["spikes"] = len(st)
        dur = t[-1] - t[0] if t else 1
        d["rate_Hz"] = round(len(st) / dur, 1) if dur > 0 else None
    return d


def main():
    args = sys.argv[1:]
    asim = args[0]
    mode = "graded"
    jsonout = None
    if "--mode" in args:
        mode = args[args.index("--mode") + 1]
    if "--json" in args:
        jsonout = args[args.index("--json") + 1]
    root = ET.parse(asim).getroot()
    nid2name = {txt(n, "ID"): txt(n, "Name") for n in root.iter("Neuron")}
    result = {"asim": asim, "mode": mode, "charts": {}}
    folder = os.path.dirname(os.path.abspath(asim))
    for ch in root.iter("DataChart"):
        cname = txt(ch, "Name")
        fname = txt(ch, "OutputFilename")
        path = os.path.join(folder, fname)
        if not os.path.exists(path):
            result["charts"][cname] = {"file": fname, "missing": True}
            continue
        header, rows = read_chart(path)
        if not rows:
            result["charts"][cname] = {"file": fname, "empty": True}
            continue
        entry = {"file": fname, "rows": len(rows), "t_end": rows[-1][1], "cols": []}
        t = [r[1] for r in rows]
        for ci, col in enumerate(header[2:], start=2):
            v = [r[ci] for r in rows if len(r) > ci]
            tt = [r[1] for r in rows if len(r) > ci]
            dt = txt(ch, "DataType")
            # find the DataColumn to learn its DataType
            ddt = None
            for c in ch.iter("DataColumn"):
                if txt(c, "ColumnName") == col:
                    ddt = txt(c, "DataType")
                    break
            if ddt == "MembraneVoltage":
                entry["cols"].append(analyze_col(col, tt, v, mode))
            elif ddt in ("JointRotation", "JointAngle", "ActualPosition"):
                if (col or "").lower().endswith("_deg"):
                    deg = v  # phase1 instrumentation already scales to degrees
                else:
                    deg = [x * 180 / math.pi for x in v]
                entry["cols"].append({
                    "name": col, "type": "angle_deg",
                    "min": round(min(deg), 1), "max": round(max(deg), 1),
                    "range": round(max(deg) - min(deg), 1)})
            elif ddt and "Contact" in (ddt or ""):
                entry["cols"].append({
                    "name": col, "type": ddt,
                    "min": min(v), "max": max(v), "mean": round(sum(v) / len(v), 3)})
            else:
                entry["cols"].append({"name": col, "type": ddt or "?", "mean": round(sum(v) / len(v), 5)})
        result["charts"][cname] = entry

    # RG/PF/MN burst metrics across all charts
    result["rg_bursts"] = {}
    envelopes = {}
    tgrid = None
    for cname, entry in result["charts"].items():
        if entry.get("missing") or entry.get("empty"):
            continue
        header, rows = read_chart(os.path.join(folder, entry["file"]))
        if not rows:
            continue
        t = [r[1] for r in rows]
        for i, nm in enumerate(header[2:], start=2):
            if not re.search(r"(RG|PF_|CPG|MN_)", nm) or "IN" in nm:
                continue
            v = [r[i] for r in rows if len(r) > i]
            if mode == "graded":
                onsets = merge_events(upcrossings(t, v, -0.060), 0.05)
            else:
                st = spike_times(t, v)
                if not st:
                    result["rg_bursts"][nm] = {"bursts": 0}
                    continue
                onsets = [st[0]]
                for a, b in zip(st, st[1:]):
                    if b - a > 0.05:
                        onsets.append(b)
            if len(onsets) >= 2:
                iv = sorted([b - a for a, b in zip(onsets, onsets[1:])])
                med = iv[len(iv) // 2]
                result["rg_bursts"][nm] = {
                    "bursts": len(onsets), "period_s": round(med, 4),
                    "freq_Hz": round(1 / med, 3) if med > 0 else None}
            elif len(onsets) == 1:
                result["rg_bursts"][nm] = {"bursts": 1, "note": "single onset (latch?)"}
            else:
                result["rg_bursts"][nm] = {"bursts": 0}
            if mode == "graded":
                grid, env = envelope(t, upcrossings(t, v, -0.060), 0.02)
            else:
                grid, env = envelope(t, spike_times(t, v), 0.02)
            envelopes[nm] = env
            tgrid = grid
    # antiphase: prefer RG_L_stance vs RG_R_stance, else L/R flx pair
    pair = None
    keys = list(envelopes)
    for a in keys:
        for b in keys:
            if a >= b:
                continue
            same_kind = (("stance" in a and "stance" in b) or
                         ("flx" in a.lower() and "flx" in b.lower()))
            lr = ("RG_L" in a and "RG_R" in b) or ("RG_R" in a and "RG_L" in b)
            if lr and (same_kind or "inhabit" not in a + b):
                if same_kind:
                    pair = (a, b)
                    break
                pair = pair or (a, b)
        if pair and ("stance" in pair[0] and "stance" in pair[1] or
                     "flx" in pair[0].lower() and "flx" in pair[1].lower()):
            break
    if pair is None:
        flxnames = [n for n in keys if "flx" in n.lower()]
        for a in range(len(flxnames)):
            for b in range(a + 1, len(flxnames)):
                if flxnames[a][0] != flxnames[b][0]:
                    pair = (flxnames[a], flxnames[b])
    if pair and tgrid:
        a, b = pair
        n = min(len(envelopes[a]), len(envelopes[b]))
        ea, eb = envelopes[a][:n], envelopes[b][:n]
        tt = [0.02 * i for i in range(n)]
        lag, r = xcorr_lag(tt, ea, tt, eb, int(2.0 / 0.02))
        per = result["rg_bursts"].get(a, {}).get("period_s")
        cyc = None
        if lag is not None and per:
            cyc = (lag / per + 0.5) % 1.0 - 0.5  # fold into [-0.5, 0.5)
        result["antiphase"] = {"pair": [a, b], "lag_s": lag,
                               "r": round(r, 3) if r is not None else None,
                               "lag_cycle_folded": round(cyc, 3) if cyc is not None else None}
    print(json.dumps(result, indent=1))
    if jsonout:
        with open(jsonout, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=1)


if __name__ == "__main__":
    main()
