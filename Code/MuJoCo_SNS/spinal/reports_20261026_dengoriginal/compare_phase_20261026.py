r"""Phase-offset comparison of every perturbation run vs the baseline.

For each recorded column: per-burst onset offset vs the NEAREST baseline
onset, wrapped to +-half the baseline period, reported as % of the baseline
period. Segments: in-window (2.0-2.5 s), post-window (>3.0 s steady).
Positive % = the perturbed run fires LATER than baseline (phase delay),
negative = earlier (phase advance). Paper benchmarks: RG test ~15% phase
DELAY; knee-PF test ~5% phase ADVANCE (2023 model) / no effect (2019 family).

Usage: python compare_phase_20261026.py
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
B0 = "run_B0_baseline"
RUNS = ["run_P1_rgflexor", "run_P2_kneepfflexor", "run_P3_authorsblock",
        "run_P1b_rgflexor_2p1", "run_P2b_kneepfflexor_2p1"]
KEY = ["LH_RG_E", "LH_RG_F", "RH_RG_E", "RH_RG_F",
       "LH_KneeMN_E", "LH_KAPF_E", "LH_KAPF_F", "LH_HipPF_E"]

def onsets(run, col):
    a = json.load(open(os.path.join(HERE, run, "analysis.json")))
    st = a[col]
    return (np.array(st["onsets"]) if st["n_bursts"] > 1 else np.array([])), st

def phase_offsets(o, b, per):
    if len(o) == 0 or len(b) == 0:
        return None
    offs = []
    for x in o:
        d = x - b[np.argmin(np.abs(b - x))]
        while d > 0.5 * per: d -= per
        while d < -0.5 * per: d += per
        offs.append(d / per * 100.0)
    return np.array(offs)

def main():
    print("%-24s %-11s %10s %10s %10s" % ("run", "column", "in-win %", "post %", "post sd"))
    for run in RUNS:
        for col in KEY:
            o, ost = onsets(run, col)
            b, bst = onsets(B0, col)
            if len(o) == 0 or len(b) == 0 or "period_mean_s" not in bst:
                print("%-24s %-11s %10s" % (run, col, "n/a"))
                continue
            per = bst["period_mean_s"]
            inw = phase_offsets(o[(o >= 2.0) & (o <= 2.5)], b, per)
            post = phase_offsets(o[o > 3.0], b, per)
            f = lambda v: ("%+7.1f%%" % np.mean(v)) if v is not None and len(v) else "   --"
            fs = lambda v: ("%5.1f" % np.std(v)) if v is not None and len(v) else " --"
            print("%-24s %-11s %10s %10s %10s" % (run, col, f(inw), f(post), fs(post)))

if __name__ == "__main__":
    main()
