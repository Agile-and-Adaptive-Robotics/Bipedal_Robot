r"""Analyze one AnimatLab BilateralRG headless run folder (2026-10-16 campaign).

Reads the chart .txt byproducts AnimatSimulator writes next to the asim:
  Rhythm Generator.txt  (cols: L/R RG ext IN, flx IN, ext, flx, L/R c1, L/R V3)
  Contact.txt           (cols: toe_L/R_contact, foot_L/R_contact counts)

Metrics (definitions from the 2026-09-16/18/24 session notes):
  - values are VOLTS in the file -> mV here
  - trailing all-zero chart rows are popped before analysis (the known
    end-of-file zero-fill trap, 2026-09-18 note)
  - rhythm health = FLEXOR RG up-crossings of -55 mV (the documented rule:
    extensor HCs ride tonic in the Rybak regime; count the flexors)
  - period = median inter-onset interval; antiphase = median phase of L
    onsets within consecutive R onsets (0.5 = perfect antiphase)
  - latch check = per-column mean/min/max (vmean.pl equivalent)

Usage:  python analyze_run_20261016.py <folder> [--json out.json]
"""
from __future__ import annotations

import json
import os
import statistics as st
import sys


def load_chart(path: str):
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
    if not rows:
        return hdr, None
    # pop trailing rows where every value column (2..) is exactly 0
    def all_zero(r):
        return all(v == 0.0 for v in r[2:])
    while rows and all_zero(rows[-1]):
        rows.pop()
    return hdr, rows


def crossings(t, v, thr):
    """up-crossing onset times of v over thr (volts)"""
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

    hdr, rows = load_chart(os.path.join(folder, "Rhythm Generator.txt"))
    if rows is None:
        res["rhythm_chart"] = "MISSING/EMPTY"
        print("Rhythm Generator.txt MISSING or empty -- run failed?")
    else:
        t = [r[1] for r in rows]
        cols = {}
        for ci, name in enumerate(hdr[2:], start=2):
            v = [r[ci] for r in rows]
            mv = [x * 1000.0 for x in v]
            cols[name] = dict(
                mean_mV=round(st.mean(mv), 3),
                min_mV=round(min(mv), 3),
                max_mV=round(max(mv), 3),
            )
        res["columns"] = cols

        thr = -0.055  # volts
        on = {}
        for ci, name in enumerate(hdr[2:], start=2):
            if name in ("L RG flx", "R RG flx", "L RG ext", "R RG ext"):
                on[name] = crossings(t, [r[ci] for r in rows], thr)
        bursts = {k: len(v) for k, v in on.items()}
        res["bursts_-55mV"] = bursts
        per = {}
        for k in ("L RG flx", "R RG flx"):
            if len(on[k]) >= 3:
                d = [b - a for a, b in zip(on[k], on[k][1:])]
                per[k] = round(st.median(d), 4)
        res["period_s_flx"] = per

        # antiphase: L onset phase within consecutive R onsets
        oL, oR = on.get("L RG flx", []), on.get("R RG flx", [])
        ph = []
        for a in range(len(oR) - 1):
            for b in oL:
                if oR[a] <= b <= oR[a + 1]:
                    ph.append((b - oR[a]) / (oR[a + 1] - oR[a]))
        if ph:
            res["antiphase_cycle"] = round(st.median(ph), 4)
            res["antiphase_n"] = len(ph)

    hdr2, rows2 = load_chart(os.path.join(folder, "Contact.txt"))
    if rows2 is None:
        res["contact_chart"] = "MISSING/EMPTY"
    else:
        duty = {}
        for ci, name in enumerate(hdr2[2:], start=2):
            v = [r[ci] for r in rows2]
            duty[name] = round(sum(1 for x in v if x > 0) / max(len(v), 1), 4)
        res["contact_duty"] = duty

    print(json.dumps(res, indent=2))
    if json_out:
        with open(json_out, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
