r"""Analyze the 2026-10-20 tuned-winner AnimatLab drop run (30 s, headless).

Reads the chart .txt byproducts next to the asim copy (run_tuned\):
  Root Track.txt        RootX/Y/Z (world) + hip_L/hip_R JointRotation
  Rhythm Generator.txt  L/R RG ext IN, flx IN, ext, flx, L/R c1, L/R V3 (V)
  Contact.txt           toe/foot contact counts (L/R)
  L Angles.txt          ankle_L, hip_L, knee_L JointRotation (rad)

Metrics (conventions from the 09-16/18/24 session notes + the 2026-10-16
analyzer): volts->mV; trailing all-zero rows popped; rhythm health =
FLEXOR RG up-crossings of -55 mV; period = median inter-onset interval;
antiphase = median L-onset phase within consecutive R onsets.

Drop protocol verdicts:
  - STAY-UP: RootY (height) stays above a 0.55 m floor (half the initial
    0.993 m pelvis height) through t=30 s.
  - TOPPLE (Ben's GUI failure mode, "fell, circling on its side"):
    RootY collapsing + hip rotations drifting to a rail + contacts going
    continuous/zero while the RG still bursts.
  - PROGRESSION: RootX net displacement + linear-fit speed.

Usage: python analyze_run_20261020.py <run_folder> [--json out.json]
"""
from __future__ import annotations

import json
import os
import statistics as st
import sys


def load_chart(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None, None
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = [ln.rstrip("\n").rstrip("\r") for ln in f if ln.strip()]
    hdr = lines[0].split("\t")
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) != len(hdr):
            continue
        try:
            rows.append([float(v) for v in parts])
        except ValueError:
            continue

    def all_zero(r):
        return all(v == 0.0 for v in r[2:])

    while rows and all_zero(rows[-1]):
        rows.pop()
    return hdr, rows


def crossings(t, v, thr):
    out = []
    for i in range(1, len(v)):
        if v[i - 1] <= thr < v[i]:
            out.append(t[i])
    return out


def main():
    folder = sys.argv[1]
    json_out = None
    if "--json" in sys.argv:
        json_out = sys.argv[sys.argv.index("--json") + 1]
    res = {"folder": folder}

    # ---- Root Track: drop / topple / progression --------------------------
    hdr, rows = load_chart(os.path.join(folder, "Root Track.txt"))
    if rows is None:
        res["root"] = "MISSING/EMPTY"
    else:
        t = [r[1] for r in rows]
        X = [r[2] for r in rows]
        Y = [r[3] for r in rows]
        Z = [r[4] for r in rows]
        hL = [r[5] for r in rows]
        hR = [r[6] for r in rows]
        n = len(t)
        # linear fit speed on X
        mx, my = st.mean(t), st.mean(X)
        sxx = sum((a - mx) ** 2 for a in t)
        slope = sum((a - mx) * (b - my) for a, b in zip(t, X)) / max(sxx, 1e-12)
        res["root"] = dict(
            t_end=t[-1], n=n,
            RootY_start=round(Y[0], 4), RootY_min=round(min(Y), 4),
            RootY_min_t=round(t[Y.index(min(Y))], 3),
            RootY_end=round(Y[-1], 4),
            RootY_mean_last10s=round(st.mean(Y[int(n * 2 / 3):]), 4),
            RootX_start=round(X[0], 4), RootX_end=round(X[-1], 4),
            RootX_net=round(X[-1] - X[0], 4), RootX_slope_mps=round(slope, 5),
            RootZ_start=round(Z[0], 4), RootZ_end=round(Z[-1], 4),
            RootZ_absmax=round(max(abs(v) for v in Z), 4),
            hip_L_rad_minmax=(round(min(hL), 3), round(max(hL), 3)),
            hip_L_end=round(hL[-1], 3),
            hip_R_rad_minmax=(round(min(hR), 3), round(max(hR), 3)),
            hip_R_end=round(hR[-1], 3),
        )
        res["root"]["stay_up"] = bool(min(Y) > 0.55)

    # ---- Rhythm Generator -------------------------------------------------
    hdr, rows = load_chart(os.path.join(folder, "Rhythm Generator.txt"))
    if rows is None:
        res["rhythm_chart"] = "MISSING/EMPTY"
    else:
        t = [r[1] for r in rows]
        cols = {}
        for ci, name in enumerate(hdr[2:], start=2):
            v = [r[ci] for r in rows]
            mv = [x * 1000.0 for x in v]
            cols[name] = dict(mean_mV=round(st.mean(mv), 3),
                              min_mV=round(min(mv), 3),
                              max_mV=round(max(mv), 3))
        res["columns"] = cols
        thr = -0.055
        on = {}
        for name in ("L RG flx", "R RG flx", "L RG ext", "R RG ext"):
            ci = hdr.index(name)
            on[name] = crossings(t, [r[ci] for r in rows], thr)
        res["bursts_-55mV"] = {k: len(v) for k, v in on.items()}
        per = {}
        for k in ("L RG flx", "R RG flx"):
            if len(on[k]) >= 3:
                d = [b - a for a, b in zip(on[k], on[k][1:])]
                per[k] = round(st.median(d), 4)
        res["period_s_flx"] = per
        oL, oR = on.get("L RG flx", []), on.get("R RG flx", [])
        ph = []
        for a in range(len(oR) - 1):
            for b in oL:
                if oR[a] <= b <= oR[a + 1]:
                    ph.append((b - oR[a]) / (oR[a + 1] - oR[a]))
        if ph:
            res["antiphase_cycle"] = round(st.median(ph), 4)
            res["antiphase_n"] = len(ph)
        # late-window bursts (does the rhythm LAST?)
        late = [b for b in oL + oR if b > t[-1] - 10]
        res["bursts_last10s"] = len(late)

    # ---- Contact duty -----------------------------------------------------
    hdr, rows = load_chart(os.path.join(folder, "Contact.txt"))
    if rows is None:
        res["contact_chart"] = "MISSING/EMPTY"
    else:
        duty = {}
        for ci, name in enumerate(hdr[2:], start=2):
            v = [r[ci] for r in rows]
            n3 = len(v) // 3
            duty[name] = dict(
                all=round(sum(1 for x in v if x > 0) / max(len(v), 1), 4),
                last10s=round(sum(1 for x in v[-n3:] if x > 0) / max(n3, 1), 4))
        res["contact_duty"] = duty

    # ---- L Angles ---------------------------------------------------------
    hdr, rows = load_chart(os.path.join(folder, "L Angles.txt"))
    if rows is not None:
        ang = {}
        for ci, name in enumerate(hdr[2:], start=2):
            v = [r[ci] for r in rows]
            n3 = len(v) // 3
            osc = 0
            for i in range(1, len(v)):
                if (v[i - 1] < 0 <= v[i]) or (v[i - 1] > 0 >= v[i]):
                    osc += 1
            ang[name] = dict(min=round(min(v), 3), max=round(max(v), 3),
                             end=round(v[-1], 3),
                             zero_crossings=osc,
                             std_last10s=round(st.pstdev(v[-n3:]), 4))
        res["L_angles"] = ang

    print(json.dumps(res, indent=2))
    if json_out:
        with open(json_out, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
