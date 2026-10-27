r"""W2L robot WEAN campaign -- single-run evaluator (overnight campaign
2026-10-16, egg; stage "feedback optimizer (robot)").

Runs test_w2l_ground.run_phase("walk") IN-PROCESS with the given knobs and
computes the campaign metrics from the RAW log, using EXACTLY the bootstrap
stage's metric definitions (reports_20261016\w2l_walk_sweep\
_write_best_params.py): primary stride cycles = heel-force rising edges
>20 N merged at 0.45 s (double-tap immune), inter-leg phase = L primary
strike phase within consecutive R primary strikes, sustained speed =
linear fit of COM x over t >= 8 s, heel duty = frac > 5 N, combined duty
= heel OR toe > 5 N.

SURVIVAL (the ask's alternating-progression criterion, >= 6 cycles):
    finite AND no fall AND min(cycL,cycR) >= 6 AND inter-leg phase median
    in [0.35, 0.65] AND forward disp > 0.05 m.

OBJECTIVE (maximize; multi-term by design -- Ben's anti-lottery rule; the
per-trial diagnostics CSV + the wean ladder itself are the structural
checks):
    -100                                      if non-finite or fell
    5*min(cyc,12)                             cycle count (gate-like, capped)
    + 200*clip(slope, 0, 0.2)                 sustained progression speed
    + 20*(1 - clip(|phase-0.5|/0.5, 0, 1))    antiphase quality
    + 5*clip(mean_combined_duty/0.4, 0, 1)    loading toward campaign target
    - 0.5*max(tilt_max-16, 0)                 lean penalty past bootstrap's

Every run appends a full diagnostics row to wean_trials.csv.

CLI (standalone probe):
    D:\Anaconda\envs\myo\python.exe _w2l_wean_eval.py --id=NAME
        [--csv=path.csv] [--dscale=1.0] [--dur=30] [--g-heel=..] ...
(knob names = test_w2l_ground KV keys, dashes or underscores).
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
W2L = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
if W2L not in sys.path:
    sys.path.insert(0, W2L)

import numpy as np                     # noqa: E402

DEFAULT_CSV = os.path.join(HERE, "wean_trials.csv")
# bootstrap winner (reports_20261016\w2l_walk_sweep\w2l_walk_best_params.json)
BOOT_TE, BOOT_TF = 4.0, 5.0

METRIC_KEYS = ("finite", "fall", "fall_t", "cycL", "cycR", "rawL", "rawR",
               "perL", "perR", "phase_med", "phase_n", "rg_r", "fwd",
               "v_mean", "slope", "span", "dutyhL", "dutyhR", "dutycL",
               "dutycR", "tilt_max", "pelz_min", "harness", "survive",
               "obj", "wall_s")
KNOB_KEYS = ("dscale", "te", "tf", "g_heel", "g_toe_in", "g_toe_df",
             "g_toe_df_mn", "g_ib", "comm", "ctn_amp", "cref", "ibnA",
             "rig", "walker_yield", "sink", "legcoll", "dur_walk")

_argv = sys.argv[:]
sys.argv = ["test_w2l_ground.py"]      # keep the test's argv parser quiet
import test_w2l_ground as T            # noqa: E402
sys.argv = _argv

T.write_ground_xml()                   # idempotent; ensures w2l_ground.xml

# base config = the bootstrap winner's protocol/body knobs (fixed this
# stage; only feedback gains + drive scale vary)
BASE = dict(legcoll=0.0, sink=0.06, rig=0.25, walker_yield=1.5)


def _primary(f: np.ndarray, t: np.ndarray, thr: float = 20.0,
             merge_s: float = 0.45):
    """Rising edges of (f>thr); PRIMARY = merged when closer than merge_s
    (the winner's heel double-tap sits ~0.3-0.45 s after the strike)."""
    on = f > thr
    st = np.flatnonzero(on[1:] & ~on[:-1]) + 1
    prim = []
    for x in st:
        if not prim or (t[x] - t[prim[-1]]) > merge_s:
            prim.append(int(x))
    return prim, [int(x) for x in st]


def _phase(tw: np.ndarray, pL, pR):
    """L primary-strike phase within consecutive R primary strikes."""
    ph = []
    for a in range(len(pR) - 1):
        for b in pL:
            if tw[pR[a]] <= tw[b] <= tw[pR[a + 1]]:
                ph.append((tw[b] - tw[pR[a]]) /
                          (tw[pR[a + 1]] - tw[pR[a]]))
    return ph


def objective_from(m: dict) -> float:
    if (not m["finite"]) or m["fall"]:
        return -100.0
    cyc = min(m["cycL"], m["cycR"])
    v = m["slope"] if m["slope"] is not None else 0.0
    aph = abs(m["phase_med"] - 0.5) if m["phase_med"] is not None else 0.5
    duty = (m["dutycL"] + m["dutycR"]) / 2.0
    return (5.0 * min(cyc, 12)
            + 200.0 * min(max(v, 0.0), 0.2)
            + 20.0 * (1.0 - min(aph / 0.5, 1.0))
            + 5.0 * min(duty / 0.4, 1.0)
            - 0.5 * max(m["tilt_max"] - 16.0, 0.0))


def _j(x):
    if isinstance(x, float) and not np.isfinite(x):
        return None
    if isinstance(x, (float, np.floating)):
        return round(float(x), 5)
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    return x


def evaluate(knobs: dict, tag: str = "", csv_path: str | None = None,
             want_log: bool = False) -> dict:
    """One walk-phase run. knobs keys = test_w2l_ground KV names."""
    # Register supported conductances at the evaluation boundary, including
    # fresh-process verification (which does not import the search driver).
    if str(T.KV.get("net")) == "bilateralrg2":
        sys.path.insert(0, os.path.join(T.HERE, "..", "w2l_cpg"))
        import build_w2l_net as B
        for key in knobs:
            if key in B.W2L_GAINS:
                T.KV.setdefault(key, float(B.W2L_GAINS[key]))
    unknown = set(knobs) - set(T.KV) - {"dscale"}
    if unknown:
        raise ValueError(f"Unsupported evaluation knobs: {sorted(unknown)}")
    kv = {k: float(v) for k, v in knobs.items()
          if k in T.KV and k != "phase"}
    T.KV.update(kv)
    t0 = time.time()
    wlog = T.run_phase("walk", T.KV["dur_walk"], T.KV["rig"], None)
    wall = time.time() - t0

    t = wlog["t"]
    w0 = T.KV["hold"] + T.KV["lower"] + 0.5
    w = t >= w0
    tw = t[w]
    heel = wlog["heelF"][w]
    toe = wlog["toeF"][w]
    com = wlog["com"][w]
    pel = wlog["pel"][w]
    tilt = wlog["tilt"][w]
    rg = wlog["rgE"][w]
    span = float(tw[-1] - tw[0]) if len(tw) else 0.0

    pL, sL = _primary(heel[:, 0], tw)
    pR, sR = _primary(heel[:, 1], tw)
    perL = float(np.median(np.diff(tw[pL]))) if len(pL) >= 3 else None
    perR = float(np.median(np.diff(tw[pR]))) if len(pR) >= 3 else None
    ph = _phase(tw, pL, pR)
    phase_med = float(np.median(ph)) if ph else None
    rg_r = (float(np.corrcoef(rg[::5, 0], rg[::5, 1])[0, 1])
            if len(rg) > 10 else None)
    fwd = float(com[-1, 0] - com[0, 0])
    m2 = tw >= 8.0
    slope = (float(np.polyfit(tw[m2], com[m2, 0], 1)[0])
             if int(m2.sum()) >= 10 else None)
    m = dict(
        finite=bool(wlog["finite"]),
        fall=wlog["fall_t"] is not None,
        fall_t=(round(float(wlog["fall_t"]), 2)
                if wlog["fall_t"] is not None else None),
        cycL=len(pL), cycR=len(pR), rawL=len(sL), rawR=len(sR),
        perL=perL, perR=perR,
        phase_med=phase_med, phase_n=len(ph), rg_r=rg_r,
        fwd=fwd, v_mean=(fwd / span if span else None), slope=slope,
        span=span,
        dutyhL=float((heel[:, 0] > 5.0).mean()),
        dutyhR=float((heel[:, 1] > 5.0).mean()),
        dutycL=float(((heel[:, 0] > 5.0) | (toe[:, 0] > 5.0)).mean()),
        dutycR=float(((heel[:, 1] > 5.0) | (toe[:, 1] > 5.0)).mean()),
        tilt_max=float(tilt.max()), pelz_min=float(pel[:, 2].min()),
        harness=float(100.0 * wlog["rigFz"][w].mean() / 411.0),
        wall_s=round(wall, 1))
    cyc = min(m["cycL"], m["cycR"])
    m["survive"] = bool(m["finite"] and not m["fall"] and cyc >= 6
                        and phase_med is not None
                        and 0.35 <= phase_med <= 0.65 and fwd > 0.05)
    m["obj"] = objective_from(m)

    # ---- diagnostics CSV row (anti-lottery logging: every run recorded)
    rowd = {"tag": tag, "ts": time.strftime("%H:%M:%S")}
    rowd.update({k: _j(T.KV.get(k, knobs.get(k))) for k in KNOB_KEYS
                 if k != "dscale"})
    rowd["dscale"] = _j(knobs.get("dscale", T.KV.get("te", 0.0) / BOOT_TE))
    rowd.update({k: _j(m[k]) for k in METRIC_KEYS})
    cp = csv_path or DEFAULT_CSV
    new = not os.path.exists(cp)
    with open(cp, "a", newline="", encoding="utf-8") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(rowd.keys()))
        if new:
            wtr.writeheader()
        wtr.writerow(rowd)

    print(f"[{tag}] survive={m['survive']} obj={m['obj']:.1f} "
          f"cyc {m['cycL']}/{m['cycR']} phase "
          f"{m['phase_med'] if m['phase_med'] is None else round(m['phase_med'], 3)} "
          f"slope {m['slope'] if m['slope'] is None else round(m['slope'], 4)} "
          f"fwd {m['fwd']:+.3f} duty {m['dutycL']:.2f}/{m['dutycR']:.2f} "
          f"tilt {m['tilt_max']:.1f} ({m['wall_s']}s)", flush=True)
    if want_log:
        m["_log"] = wlog
    return m


def main() -> int:
    cid, csv_path, knobs = "probe", None, dict(BASE)
    dur = 30.0
    dscale = 1.0
    for a in sys.argv[1:]:
        if a.startswith("--id="):
            cid = a[5:]
        elif a.startswith("--csv="):
            csv_path = a[6:]
        elif a.startswith("--dur="):
            dur = float(a[6:])
        elif a.startswith("--dscale="):
            dscale = float(a[9:])
        elif a.startswith("--") and "=" in a:
            k, v = a[2:].split("=", 1)
            try:
                v = float(v)
            except ValueError:
                pass
            knobs[k.replace("-", "_")] = v
    knobs.update(dur_walk=dur, te=BOOT_TE * dscale, tf=BOOT_TF * dscale,
                 dscale=dscale)
    m = evaluate(knobs, tag=cid, csv_path=csv_path)
    print("RESULT " + json.dumps({k: _j(m[k]) for k in METRIC_KEYS}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
