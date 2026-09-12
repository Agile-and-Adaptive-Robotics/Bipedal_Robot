"""Kinematics similarity vs the OpenSim IK benchmark (Ben 2026-09-12:
"fine-tune until the kinematics are similar to OpenSim").

Reference: subject01_walk1_ik.mot (degrees, flexion-negative knee - the
same convention the converted model preserved) phased by the measured
GRF: the right gait cycle runs from one right-foot loading onset
(vertical GRF > 50 N, rising) to the next. The sim side: RG-E_r burst
onsets over the walk window of a runner log, each cycle interpolated to
a common phase grid and averaged.

Comparison (all on the MEAN cycle, shape = mean-offset removed):
    rmse_hip / rmse_knee / rmse_arm  [deg]   shape RMSE vs reference
    knee_min    [deg]   peak flexion (ref ~ -60)
    *_range     [deg]   hip / ankle excursions of the mean cycle
    duty                fraction of the cycle RG-E (stance) is on
    kine_score          single number, HIGHER IS BETTER (0 = perfect):
        -(1.2 rmse_knee + 1.0 rmse_hip + 0.8 rmse_ankle)
        -0.15 |knee_min err| - 0.10 (|hip range err| + |ankle range err|)
        -3.0 |duty err|

Used by runner --eval (metrics["kine"]) and optuna_walk v4. No
matplotlib here - safe to import inside the sim loop.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

REPO = Path(__file__).parents[3]
IK_MOT = (REPO / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
          / "subject01_walk1_ik.mot")
GRF_MOT = (REPO / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
           / "subject01_walk1_grf.mot")
NPHASE = 100
GRID = np.linspace(0.0, 100.0, NPHASE)


def _read_mot(path: Path):
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
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


def load_reference():
    """Reference mean cycle (degrees) + ranges + duty from the IK + GRF."""
    t_ik, names, vals = _read_mot(IK_MOT)
    # vals excludes the time column -> build the map from names[1:]
    col = {n: i for i, n in enumerate(names[1:])}
    t_g, gn, gv = _read_mot(GRF_MOT)
    # gv excludes the time column -> name index in gn is off by one
    vy = gv[:, gn.index("ground_force_vy") - 1]
    on = (vy > 50.0).astype(int)
    onsets = t_g[np.flatnonzero(np.diff(on) == 1)]
    if len(onsets) < 2:
        raise RuntimeError("no full gait cycle in the GRF reference")
    t0, t1 = onsets[0], onsets[1]
    m = (t_ik >= t0) & (t_ik < t1)

    def cyc(name):
        return np.interp(GRID, (t_ik[m] - t0) / (t1 - t0) * 100.0,
                         vals[m, col[name]])

    ref = dict(hip=cyc("hip_flexion_r"), knee=cyc("knee_angle_r"),
               ankle=cyc("ankle_angle_r"))
    on_cyc = np.interp(GRID, (t_g - t0) / (t1 - t0) * 100.0, vy)
    ref["duty"] = float(np.mean(on_cyc > 50.0))
    ref["t0"], ref["t1"] = float(t0), float(t1)
    for j in ("hip", "knee", "ankle"):
        ref[f"{j}_range"] = float(np.ptp(ref[j]))
    ref["knee_min"] = float(np.min(ref["knee"]))
    return ref


def sim_cycles(t, q_deg, neuro, walk_start, rge_col=2):
    """Mean simulated right-leg cycle (degrees) over RG-E_r onsets in the
    walk window. q_deg columns: 3 hip_flexion_r, 4 knee_angle_r,
    5 ankle_angle_r (runner KEY_JOINTS order). Returns (mean cycle dict,
    n_cycles, duty) or (None, 0, nan)."""
    m = t >= walk_start
    if not m.any():
        return None, 0, float("nan")
    rge = neuro[m, rge_col]
    on = rge > 0.5 * max(np.max(rge), 1e-9)
    rises = np.flatnonzero(np.diff(on.astype(int)) == 1)
    falls = np.flatnonzero(np.diff(on.astype(int)) == -1)
    if len(rises) < 3:
        return None, 0, float("nan")
    tt = t[m]
    cycles = {j: [] for j in ("hip", "knee", "ankle")}
    duty_fracs = []
    for a, b in zip(rises[:-1], rises[1:]):
        dur = tt[b] - tt[a]
        if not (0.35 <= dur <= 2.5) or b >= len(tt):
            continue
        ph = (tt[a:b] - tt[a]) / dur * 100.0
        for j, c in zip(("hip", "knee", "ankle"), (3, 4, 5)):
            cycles[j].append(np.interp(GRID, ph, q_deg[m][a:b, c]))
        seg = on[a:b]
        duty_fracs.append(float(np.mean(seg)))
    if len(cycles["hip"]) < 2:
        return None, 0, float("nan")
    mean = {j: np.mean(cycles[j], axis=0) for j in cycles}
    return mean, len(cycles["hip"]), float(np.mean(duty_fracs))


def compare(t, q_deg, neuro, walk_start, ref=None):
    """Full comparison dict (incl. kine_score) or None if the run has no
    usable rhythm."""
    ref = ref or REF_CACHE
    if ref is None:
        ref = load_reference()
    mean, n_cyc, duty = sim_cycles(t, q_deg, neuro, walk_start)
    if mean is None:
        return None
    out = dict(n_cycles=n_cyc, duty=duty, duty_ref=ref["duty"],
               knee_min=float(np.min(mean["knee"])),
               knee_min_ref=ref["knee_min"])
    total = 0.0
    wts = dict(hip=1.0, knee=1.2, ankle=0.8)
    for j, w in wts.items():
        s = mean[j] - np.mean(mean[j])
        r = ref[j] - np.mean(ref[j])
        rmse = float(np.sqrt(np.mean((s - r) ** 2)))
        out[f"rmse_{j}"] = rmse
        out[f"range_{j}"] = float(np.ptp(mean[j]))
        out[f"range_{j}_ref"] = ref[f"{j}_range"]
        total += w * rmse
        if j in ("hip", "ankle"):
            total += 0.10 * abs(out[f"range_{j}"] - out[f"range_{j}_ref"])
    total += 0.15 * abs(out["knee_min"] - ref["knee_min"])
    total += 3.0 * abs(duty - ref["duty"])
    out["kine_score"] = float(-total)
    return out


REF_CACHE = None


def ref_cached():
    global REF_CACHE
    if REF_CACHE is None:
        REF_CACHE = load_reference()
    return REF_CACHE


if __name__ == "__main__":
    r = ref_cached()
    print(f"reference cycle {r['t0']:.2f}-{r['t1']:.2f} s, duty "
          f"{r['duty']:.2f}, knee_min {r['knee_min']:.1f} deg, ranges "
          f"hip {r['hip_range']:.1f} / ankle {r['ankle_range']:.1f} deg")
