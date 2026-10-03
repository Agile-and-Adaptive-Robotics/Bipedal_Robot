"""Same-leg half-center relation for a spiking-copy run (persisted convention).

Convention (fixed so the number is reproducible):
  - chart file: "Rhythm Generator.txt" next to the asim
  - spikes: rising crossings of -40 mV with 1 ms refractory
  - envelopes: spike counts in 0.02 s bins (events/s)
  - lag scan: +/- 2.0 s in 0.02 s steps; Pearson r of ext-envelope vs SHIFTED
    flx-envelope; report the maximizing lag and r.
r > 0 near lag 0 (or lag = +/- one burst period) = co-bursting (in-phase);
r < 0 near lag = half the burst period = alternation.

Usage: python sameleg_xcorr.py <asim> <metrics.json>
Writes the `same_leg_xcorr` key into the given metrics JSON.
"""
import io
import sys
import json
import os
import math
import xml.etree.ElementTree as ET

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

BIN = 0.02
SCAN = 2.0


def txt(e, t):
    x = e.find(t)
    return x.text if x is not None else None


def read_chart(path):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    while rows and all(x == 0.0 for x in rows[-1][2:]):
        rows.pop()
    return header, rows


def spike_times(t, v, thr=-0.040, refract=0.001):
    out = []
    above = v[0] > thr
    for i in range(1, len(v)):
        a = v[i] > thr
        if a and not above:
            ti = t[i - 1] + (thr - v[i - 1]) * (t[i] - t[i - 1]) / (v[i] - v[i - 1])
            if not out or ti - out[-1] >= refract:
                out.append(ti)
        above = a
    return out


def envelope(t, events, width=BIN):
    t0, t1 = t[0], t[-1]
    n = max(1, int((t1 - t0) / width))
    rate = [0.0] * n
    for x in events:
        b = int((x - t0) / width)
        if 0 <= b < n:
            rate[b] += 1.0 / width
    return rate


def xcorr(ea, eb):
    n = len(ea)
    m1 = sum(ea) / n
    m2 = sum(eb) / n
    best_r, best_lag = -2.0, 0
    for k in range(-int(SCAN / BIN), int(SCAN / BIN) + 1):
        sxy = sxx = syy = 0.0
        for i in range(n):
            j = i + k
            if 0 <= j < n:
                sxy += (ea[i] - m1) * (eb[j] - m2)
                sxx += (ea[i] - m1) ** 2
                syy += (eb[j] - m2) ** 2
        if sxx > 0 and syy > 0:
            r = sxy / math.sqrt(sxx * syy)
            if r > best_r:
                best_r, best_lag = r, k * BIN
    return best_r, best_lag


def main():
    asim, metrics_path = sys.argv[1], sys.argv[2]
    root = ET.parse(asim).getroot()
    folder = os.path.dirname(os.path.abspath(asim))
    chart = next(ch for ch in root.iter("DataChart")
                 if (txt(ch, "OutputFilename") or "").startswith("Rhythm Generator"))
    header, rows = read_chart(os.path.join(folder, txt(chart, "OutputFilename")))
    t = [r[1] for r in rows]
    res = {}
    legs = {}
    for i, h in enumerate(header[2:], start=2):
        if h.endswith(" RG ext"):
            legs.setdefault(h.split(" ")[0], {})["ext"] = i
        if h.endswith(" RG flx"):
            legs.setdefault(h.split(" ")[0], {})["flx"] = i
    for leg, idxs in sorted(legs.items()):
        se = spike_times(t, [r[idxs["ext"]] for r in rows if len(r) > idxs["ext"]])
        sf = spike_times(t, [r[idxs["flx"]] for r in rows if len(r) > idxs["flx"]])
        ee = envelope(t, se)
        ef = envelope(t, sf)
        n = min(len(ee), len(ef))
        r, lag = xcorr(ee[:n], ef[:n])
        res[f"{leg} RG ext-vs-flx"] = {
            "ext_spikes": len(se), "flx_spikes": len(sf),
            "r": round(r, 3), "lag_s": round(lag, 3),
            "convention": "0.02 s spike-count envelopes, -40 mV rising crossings, 1 ms refractory, +/-2 s scan",
        }
    out = {"same_leg_xcorr": res}
    print(json.dumps(out, indent=1))
    if os.path.exists(metrics_path):
        with open(metrics_path, encoding="utf-8") as f:
            m = json.load(f)
        m.update(out)
        with open(metrics_path, "w", encoding="utf-8") as f:
            json.dump(m, f, indent=1)
        print(f"updated {metrics_path}")


if __name__ == "__main__":
    main()
