"""Parse AnimatLab DataTool/chart .txt traces (tab-separated, header row).

Streams line-by-line. Data-cleaning rule (documented): TRAILING rows whose
data columns (everything except TimeSlice/Time) are all exactly 0.0 are
dropped as DataTool shutdown-flush artifacts; the count and first-dropped
timestamp are reported so a genuine mid-run collapse stays visible (only
the tail is stripped).

Burst rule: adaptive threshold at the midpoint of [p05,p95] over the window,
re-armed after dropping below p05 + 0.25*range; onsets, mean period, mean
on-duration, duty fraction. Fixed-level crossings (-50/-55 mV, 1 mV
hysteresis) as cross-checks on volt-scale series. Zero-lag Pearson r plus a
lagged normalized cross-correlation on a 1 ms grid (+/- 1.5 s). FFT dominant
frequency on the median-dt uniform grid (Hann window). Post-transient window
t >= T0 = 1.0 s in addition to full trace.

Usage: python parse_animatlab_traces.py <json_spec_file>
Prints one JSON object to stdout.
"""
import json
import sys

import numpy as np

T0 = 1.0  # post-transient window start [s]


def read_tsv(path):
    with open(path, "r") as f:
        header = f.readline().rstrip("\r\n").split("\t")
        idx = {n: i for i, n in enumerate(header)}
        data_idx = [i for n, i in idx.items() if n not in ("TimeSlice", "Time")]
        cols = {n: [] for n in header}
        ncol = len(header)
        nrows = 0
        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < ncol:
                continue
            try:
                vals = [float(parts[i]) for i in range(ncol)]
            except ValueError:
                continue
            nrows += 1
            for i, n in enumerate(header):
                cols[n].append(vals[i])
    # strip TRAILING all-zero data rows (shutdown flush artifact)
    dropped = 0
    first_dropped_t = None
    rows = nrows
    while rows > 0:
        if all(cols[header[i]][rows - 1] == 0.0 for i in data_idx):
            if first_dropped_t is None:
                first_dropped_t = cols["Time"][rows - 1] if "Time" in cols else None
            rows -= 1
            dropped += 1
        else:
            break
    out = {"n_rows_raw": nrows, "n_trailing_zero_rows_dropped": dropped,
           "first_dropped_time_s": first_dropped_t,
           "cols": {n: np.asarray(cols[n][:rows], dtype=float)
                    for n in header if cols[n]}}
    return out


def post(t, v, t0=T0):
    m = t >= t0
    return t[m], v[m]


def burst_scan(t, v):
    lo, hi = np.percentile(v, [5, 95])
    thr = 0.5 * (lo + hi)
    rearm = lo + 0.25 * (hi - lo)
    onsets, offsets = [], []
    armed = True
    onset_t = None
    prev = v[0]
    for ti, vi in zip(t, v):
        if armed and vi >= thr and prev < thr:
            onset_t = ti
            armed = False
        elif (not armed) and vi < rearm:
            offsets.append(ti)
            onsets.append(onset_t)
            armed = True
        prev = vi
    if not armed:
        offsets.append(t[-1])
        onsets.append(onset_t)
    onsets = np.asarray(onsets, dtype=float)
    durs = np.asarray(offsets, dtype=float) - onsets
    return {"min": round(float(v.min()), 5), "max": round(float(v.max()), 5),
            "p05": round(float(lo), 5), "p95": round(float(hi), 5),
            "thr": round(float(thr), 5),
            "bursts": int(len(onsets)),
            "mean_period_s": round(float(np.diff(onsets).mean()), 4) if len(onsets) >= 2 else None,
            "std_period_s": round(float(np.diff(onsets).std()), 4) if len(onsets) >= 2 else None,
            "mean_dur_s": round(float(durs.mean()), 4) if len(durs) else None,
            "duty_frac": round(float(durs.sum() / (t[-1] - t[0])), 4) if len(durs) else None}


def fixed_crossings(v, level, rearm_delta):
    thr = level
    rearm = level - rearm_delta
    n = 0
    armed = True
    prev = v[0]
    for vi in v:
        if armed and vi >= thr and prev < thr:
            n += 1
            armed = False
        elif not armed and vi < rearm:
            armed = True
        prev = vi
    return n


def fft_peak(t, v, fmin=0.1, fmax=15.0):
    dt = float(np.median(np.diff(t)))
    n = len(v)
    g = np.arange(n) * dt
    vi = np.interp(g, t, v)
    vi = vi - vi.mean()
    if np.allclose(vi, 0):
        return None
    sp = np.abs(np.fft.rfft(vi * np.hanning(n)))
    fr = np.fft.rfftfreq(n, dt)
    m = (fr >= fmin) & (fr <= fmax)
    if not m.any():
        return None
    k = int(np.argmax(sp[m]))
    return round(float(fr[m][k]), 4)


def pearson(a, b):
    if len(a) != len(b) or len(a) < 3:
        return None
    a = a - a.mean()
    b = b - b.mean()
    den = np.sqrt(float((a * a).sum()) * float((b * b).sum()))
    if den == 0:
        return None
    return round(float((a * b).sum() / den), 4)


def xcorr_lag(t, a, b, max_lag_s=1.5, ds=0.001):
    """Normalized cross-correlation c(tau) = corr(a(t), b(t+tau)) on a uniform
    ds grid. tau > 0 at the peak => b LAGS a by tau."""
    ta = np.arange(0.0, (t[-1] - t[0]), ds) + t[0]
    aa = np.interp(ta, t, a)
    bb = np.interp(ta, t, b)
    aa = aa - aa.mean()
    bb = bb - bb.mean()
    n = len(aa)
    K = int(max_lag_s / ds)
    K = min(K, n - 3)
    lags = np.arange(-K, K + 1)
    rs = np.full(len(lags), np.nan)
    saa = float((aa * aa).sum())
    for j, k in enumerate(lags):
        if k >= 0:
            x, y = aa[:n - k], bb[k:]
        else:
            x, y = aa[-k:], bb[:n + k]
        sbb = float((y * y).sum())
        den = np.sqrt(saa * sbb) * (len(x) / n)
        if den > 0:
            rs[j] = float((x * y).sum()) / (np.sqrt(float((x * x).sum()) * sbb))
    j0 = int(np.nanargmax(rs))
    return {"r_zero_lag": round(float(rs[lags == 0]), 4),
            "r_max": round(float(rs[j0]), 4),
            "lag_at_rmax_s": round(float(lags[j0] * ds), 4)}


def main():
    spec_path = sys.argv[1]
    with open(spec_path) as f:
        spec = json.load(f)
    out = {"T0_post_transient_s": T0, "models": []}
    for entry in spec:
        raw = read_tsv(entry["file"])
        data = raw["cols"]
        res = {"label": entry["label"], "file": entry["file"],
               "n_rows_raw": raw["n_rows_raw"],
               "n_trailing_zero_rows_dropped": raw["n_trailing_zero_rows_dropped"],
               "first_dropped_time_s": raw["first_dropped_time_s"],
               "series": list(data.keys())}
        t = data.get("Time")
        if t is not None and len(t):
            res["t_min_s"] = round(float(t.min()), 3)
            res["t_max_s"] = round(float(t.max()), 3)
            res["dt_ms"] = round(float(np.median(np.diff(t))) * 1e3, 4)
        an = entry.get("analyses", {})
        periods = {}
        for col in an.get("bursts", []):
            if col not in data or t is None:
                continue
            full = burst_scan(t, data[col])
            full["xings_-50mV_full"] = fixed_crossings(data[col], -0.050, 1e-3)
            postw = burst_scan(*post(t, data[col]))
            postw["xings_-55mV_post"] = fixed_crossings(post(t, data[col])[1], -0.055, 1e-3)
            res.setdefault("bursts", []).append({"col": col, "full": full, "post_T0": postw})
            if full["mean_period_s"]:
                periods[col] = full["mean_period_s"]
        for col in an.get("flexcheck", []):
            if col in data and t is not None:
                vv_full = data[col]
                vv_post = post(t, vv_full)[1]
                res.setdefault("flexcheck", []).append(
                    {"col": col,
                     "full_max_mV": round(float(vv_full.max()) * 1e3, 2),
                     "post_max_mV": round(float(vv_post.max()) * 1e3, 2),
                     "xings_-55mV_full": fixed_crossings(vv_full, -0.055, 1e-3),
                     "xings_-55mV_post": fixed_crossings(vv_post, -0.055, 1e-3)})
        for a, b in an.get("corr", []):
            if a in data and b in data and t is not None:
                va = post(t, data[a])[1]
                vb = post(t, data[b])[1]
                n = min(len(va), len(vb))
                xc = xcorr_lag(post(t, data[a])[0], va[:n], vb[:n])
                xc.update({"a": a, "b": b})
                per = periods.get(a) or periods.get(b)
                if per and xc["lag_at_rmax_s"] is not None:
                    xc["lag_cycles"] = round(xc["lag_at_rmax_s"] / per, 3)
                    half = per / 2.0
                res.setdefault("corr_post", []).append(xc)
        for col in an.get("freq", []):
            if col in data and t is not None:
                tt, vv = post(t, data[col])
                pk = fft_peak(tt, vv)
                res.setdefault("freq_hz_post", []).append(
                    {"col": col, "peak_hz": pk,
                     "period_from_peak_s": round(1.0 / pk, 4) if pk else None})
        for col in an.get("range", []):
            if col in data and t is not None:
                tt, vv = post(t, data[col])
                vf = data[col]
                res.setdefault("ranges", []).append(
                    {"col": col,
                     "full_min": round(float(vf.min()), 4),
                     "full_max": round(float(vf.max()), 4),
                     "post_min": round(float(vv.min()), 4),
                     "post_max": round(float(vv.max()), 4),
                     "post_mean": round(float(vv.mean()), 4),
                     "post_ptp": round(float(vv.ptp()), 4),
                     "post_ptp_deg_if_rad": round(float(np.degrees(vv.ptp())), 2)})
        for col in an.get("duty", []):
            if col in data and t is not None:
                tt, vv = post(t, data[col])
                lo, hi = float(vv.min()), float(vv.max())
                res.setdefault("duty_post", []).append(
                    {"col": col, "min": round(lo, 4), "max": round(hi, 4),
                     "frac_gt_0.5": round(float((vv > 0.5).mean()), 4)})
        for col in an.get("ends", []):
            if col in data and t is not None:
                tt, vv = post(t, data[col])
                res.setdefault("endpoints_post", []).append(
                    {"col": col, "first": round(float(vv[0]), 4),
                     "last": round(float(vv[-1]), 4),
                     "t_first": round(float(tt[0]), 3),
                     "t_last": round(float(tt[-1]), 3),
                     "abs_rate": round(abs(float(vv[-1]) - float(vv[0])) /
                                       float(tt[-1] - tt[0]), 4)})
        out["models"].append(res)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
