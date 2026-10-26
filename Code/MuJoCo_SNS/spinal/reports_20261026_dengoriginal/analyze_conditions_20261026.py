r"""Burst analysis for the Deng-2019-original runs (Order B1).

Chart file: deng_rg_mn.txt - col0 TimeSlice, col1 Time(s), then mV-scale
voltages (chart stores VOLTS; multiply by 1000). Columns in build order.
Burst onset = upward crossing of -60 mV (skill convention for NonSpiking),
burst end = next downward crossing. Periods from onsets. In-window stats for
the 2.0-2.5 s perturbation window. Pop trailing all-zero rows first.

Usage: python analyze_conditions_20261026.py <run_folder> [--json out.json]
"""
from __future__ import annotations
import json, os, sys
import numpy as np

COLS = ["LH_RG_E", "LH_RG_F", "RH_RG_E", "RH_RG_F",
        "LH_KneeMN_E", "LH_KneeMN_F", "RH_KneeMN_E", "RH_KneeMN_F",
        "LH_KAPF_E", "LH_KAPF_F", "LH_HipPF_E"]
THRESH_V = -0.060     # volts
WIN = (2.0, 2.5)

def bursts(t, v, thresh=THRESH_V):
    on = (v[:-1] < thresh) & (v[1:] >= thresh)
    off = (v[:-1] >= thresh) & (v[1:] < thresh)
    onsets = t[:-1][on]
    offsets = t[:-1][off]
    # pair them: drop an offset before the first onset
    if offsets.size and onsets.size and offsets[0] < onsets[0]:
        offsets = offsets[1:]
    n = min(onsets.size, offsets.size)
    return onsets[:n], offsets[:n]

def stats_for(t, v):
    onsets, offsets = bursts(t, v)
    d = {}
    d["n_bursts"] = int(onsets.size)
    if onsets.size >= 2:
        per = np.diff(onsets)
        d["period_mean_s"] = float(np.mean(per))
        d["period_std_s"] = float(np.std(per))
        n2 = min(offsets.size, per.size)
        d["duty_mean"] = float(np.mean((offsets[:n2] - onsets[:n2]) / per[:n2])) if n2 else None
        d["onsets"] = [float(x) for x in onsets]
    d["v_mean_mV"] = float(np.mean(v) * 1000)
    d["v_min_mV"] = float(np.min(v) * 1000)
    d["v_max_mV"] = float(np.max(v) * 1000)
    return d, onsets, offsets

def main() -> int:
    rd = sys.argv[1]
    f = os.path.join(rd, "deng_rg_mn.txt")
    raw = np.genfromtxt(f, delimiter="\t", names=None)
    if raw.ndim == 1:
        raw = raw[None, :]
    # drop trailing all-zero rows and any NaN rows
    nonzero = np.where(np.abs(np.nan_to_num(raw)).sum(axis=1) > 0)[0]
    raw = raw[:nonzero[-1] + 1] if nonzero.size else raw
    raw = raw[~np.isnan(raw).any(axis=1)]
    t = raw[:, 1]
    out = {"run": rd, "n_rows": int(raw.shape[0]), "t_end": float(t[-1])}
    m = (t >= WIN[0]) & (t <= WIN[1])
    for k, name in enumerate(COLS):
        v = raw[:, 2 + k]
        st, onsets, offsets = stats_for(t, v)
        st["bursts_in_window"] = int(((onsets >= WIN[0]) & (onsets <= WIN[1])).sum())
        st["v_in_window_mean_mV"] = float(np.mean(v[m]) * 1000) if m.any() else None
        out[name] = st
    # antiphase check L vs R RG flexors
    if "LH_RG_F" in out and "RH_RG_F" in out and out["LH_RG_F"]["n_bursts"] > 1 and out["RH_RG_F"]["n_bursts"] > 1:
        lo = np.array(out["LH_RG_F"]["onsets"]); ro = np.array(out["RH_RG_F"]["onsets"])
        lags = []
        for r in ro:
            i = np.argmin(np.abs(lo - r))
            lags.append(lo[i] - r)
        lags = np.array(lags)
        lags = np.where(lags > 0.5 * out["LH_RG_F"]["period_mean_s"], lags - out["LH_RG_F"]["period_mean_s"], lags)
        out["L_vs_R_RGflx_lag_mean_s"] = float(np.mean(lags))
    if "--json" in sys.argv:
        o = sys.argv[sys.argv.index("--json") + 1]
        json.dump(out, open(o, "w"), indent=2)
    # compact print
    for name in COLS:
        st = out[name]
        print("%-12s bursts=%-3d in-win=%-2d per=%.3f+-%.3f s  v=[%.1f..%.1f] mV" % (
            name, st["n_bursts"], st["bursts_in_window"],
            st.get("period_mean_s", float("nan")), st.get("period_std_s", float("nan")),
            st["v_min_mV"], st["v_max_mV"]))
    if "L_vs_R_RGflx_lag_mean_s" in out:
        print("L vs R RG-flexor lag: %.3f s" % out["L_vs_R_RGflx_lag_mean_s"])
    # within-leg E vs F alternation
    for side in ("LH", "RH"):
        a, b = out.get(side + "_RG_E"), out.get(side + "_RG_F")
        if a and b and a["n_bursts"] > 1 and b["n_bursts"] > 1:
            lo = np.array(a["onsets"]); ro = np.array(b["onsets"])
            lags = []
            for r in ro:
                i = np.argmin(np.abs(lo - r))
                lg = lo[i] - r
                per = a["period_mean_s"]
                if lg > 0.5 * per: lg -= per
                if lg < -0.5 * per: lg += per
                lags.append(lg)
            out[side + "_EvsF_lag_mean_s"] = float(np.mean(lags))
            print("%s RG E-vs-F lag: %.3f s (period %.3f)" % (side, out[side + "_EvsF_lag_mean_s"], a["period_mean_s"]))
    return 0

if __name__ == "__main__":
    sys.exit(main())
