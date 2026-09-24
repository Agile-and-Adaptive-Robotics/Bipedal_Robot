"""Ong/SCONE predictive-walking reference builder (2026-09-24 night).

The results-speeds states.sto files carry everything the GRF-phased
route needs IN ONE FILE: gait2392-named joint angles (radians) AND
per-leg vertical GRF (Leg1_r.grf_y / Leg0_l.grf_y, body-weight
normalized, +up). This converts each SpecifiedSpeeds/SelfSelectedSpeeds
states file into a kine_ref-schema reference npz in gait_refs/,
matching the thumb-drive campaign convention. Deficit runs are
intentionally skipped for now (Ben: crouch library not helpful yet).
"""
import io
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import kine_ref as KR


def _read_sto_rad(path: Path):
    """OpenSim .sto parser that ACCEPTS radians (the shared _read_mot
    refuses inDegrees=no by design)."""
    lines = Path(path).read_text(encoding="utf-8",
                                  errors="ignore").splitlines()
    end = next(i for i, ln in enumerate(lines) if ln.strip() == "endheader")
    names = lines[end + 1].split()
    rows = []
    for ln in lines[end + 2:]:
        s = ln.strip()
        if not s:
            continue
        try:
            rows.append([float(v) for v in s.split()])
        except ValueError:
            break
    data = np.asarray(rows)
    return data[:, 0], names, data[:, 1:]

SPEEDS = Path(r"D:\temp\gait_lib_staging\downloads\results-speeds")
REFS_OUT = HERE / "gait_refs"
JOINTS = {"r": ("hip_flexion_r", "knee_angle_r", "ankle_angle_r"),
          "l": ("hip_flexion_l", "knee_angle_l", "ankle_angle_l")}
GRF = {"r": "Leg1_r.grf_y", "l": "Leg0_l.grf_y"}
LOAD_BW = 0.10          # loading threshold as fraction of body weight
BW_N = 750.0            # scale BW-normalized GRF into newtons so the
                        # kine_ref fixed LOAD_N=50 N onset detector works
RAD2DEG = 180.0 / np.pi


def build_ong_ref(states: Path, load_bw: float = LOAD_BW):
    t, names, vals = _read_sto_rad(states)
    col = {n: i for i, n in enumerate(names[1:])}
    # SpecifiedSpeeds use "Leg1_r.grf_y", SelfSelected "leg1_r.grf_y"
    low = {n.lower(): n for n in col}
    GRF = {"r": low.get("leg1_r.grf_y"), "l": low.get("leg0_l.grf_y")}
    for side in ("r", "l"):
        for j in JOINTS[side]:
            if j not in col:
                raise ValueError(f"missing column {j}")
        if GRF[side] is None or GRF[side] not in col:
            raise ValueError(f"missing column {GRF[side]}")
    vy = {s: vals[:, col[GRF[s]]] * BW_N for s in ("r", "l")}
    thr = load_bw * BW_N
    deg = {s: {j.split("_")[0]: vals[:, col[j]] * RAD2DEG
               for j in JOINTS[s]} for s in ("r", "l")}
    on = {s: KR._loading_onsets(t, vy[s]) for s in ("r", "l")}
    # single-cycle periodic: double the trace so each side sees 2 onsets
    if len(on["r"]) >= 2 and len(on["l"]) >= 2:
        use_t, use_vy, use_deg, use_on, periodic = t, vy, deg, on, False
    else:
        T = float(t[-1] - t[0] + np.mean(np.diff(t)))
        t2 = np.concatenate([t, t + T])
        vy2 = {s: np.concatenate([vy[s], vy[s]]) for s in ("r", "l")}
        deg2 = {s: {j: np.concatenate([deg[s][j], deg[s][j]])
                    for j in deg[s]} for s in ("r", "l")}
        on2 = {s: KR._loading_onsets(t2, vy2[s]) for s in ("r", "l")}
        use_t, use_vy, use_deg, use_on, periodic = t2, vy2, deg2, on2, True
    if len(use_on["r"]) < 2 or len(use_on["l"]) < 2:
        raise RuntimeError("no full gait cycle detected")

    ref = {"grf_style": "ong_states", "periodic": periodic}
    for side in ("r", "l"):
        t0, t1 = use_on[side][0], use_on[side][1]
        cyc = {j: np.interp(KR.GRID, (use_t - t0) / (t1 - t0) * 100.0,
                            use_deg[side][j]) for j in ("hip", "knee", "ankle")}
        ref[side] = cyc
        vyc = np.interp(KR.GRID, (use_t - t0) / (t1 - t0) * 100.0, use_vy[side])
        ref[f"duty_{side}"] = float(np.mean(vyc > thr))
        ref[f"T_{side}"] = float(t1 - t0)
        for j in ("hip", "knee", "ankle"):
            ref[f"{j}_range_{side}"] = float(np.ptp(cyc[j]))
            ref[f"mean_{j}_{side}"] = float(np.mean(cyc[j]))
        ref[f"knee_min_{side}"] = float(np.min(cyc["knee"]))
    t_r0, T_r = use_on["r"][0], ref["T_r"]
    after = use_on["l"][use_on["l"] >= t_r0]
    ref["lag_rl"] = float((after[0] - t_r0) / T_r) if after.size else 0.5
    i0 = int(np.searchsorted(use_t, use_on["r"][0]))
    i1 = int(np.searchsorted(use_t, use_on["r"][1]))
    ref["ds"] = float(np.mean((use_vy["r"][i0:i1] > thr) &
                              (use_vy["l"][i0:i1] > thr)))
    return ref


def main():
    REFS_OUT.mkdir(exist_ok=True)
    targets = []
    for d in sorted(SPEEDS.glob("SpecifiedSpeeds/*")):
        st = sorted(d.glob("*_states.sto")) or sorted(d.glob("states.sto"))
        if st:
            targets.append((f"ong_speed_{d.name}", st[0]))
    for d in sorted(SPEEDS.glob("SelfSelectedSpeeds/*")):
        st = sorted(d.glob("states.sto")) or sorted(d.glob("*_states.sto"))
        if st:
            targets.append((f"ong_selfsel_{d.name.split('-')[-1]}", st[0]))
    print(f"{len(targets)} Ong states targets")
    for name, st in targets:
        try:
            ref = build_ong_ref(st)
        except Exception as e:
            print(f"FAILED {name}: {e}")
            continue
        out = {}
        for side in ("r", "l"):
            for j in ("hip", "knee", "ankle"):
                out[f"{side}_{j}"] = ref[side][j]
        for k, v in ref.items():
            if k not in ("r", "l"):
                out[k] = v
        np.savez(REFS_OUT / f"{name}.npz", **out)
        print(f"OK {name}: T_r {ref['T_r']:.3f} s, duty "
              f"{ref['duty_r']:.2f}/{ref['duty_l']:.2f}, knee_min "
              f"{ref['knee_min_r']:.1f} deg, hip range "
              f"{ref['hip_range_r']:.1f} deg")


if __name__ == "__main__":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
    main()
