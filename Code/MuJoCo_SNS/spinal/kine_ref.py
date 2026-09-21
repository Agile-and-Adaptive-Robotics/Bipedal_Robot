"""Kinematics similarity vs the OpenSim IK benchmark — v2 (2026-09-20).

Ben's verdict on v1 (after seeing the s3b "winner" figures): "The legs
are split, one always flexed at the hip, the other not alternating, the
left side fairly static, ankles dragging. Lower the walker so it makes
contact. If our gait cycles for BOTH legs better matched OpenSim (mean,
amplitude, phase, periodicity), I'll be happy."

v1 blind spots that produced that gait (all fixed here):
  - RIGHT LEG ONLY (columns 3/4/5) -> a frozen left leg cost nothing;
  - mean-offset REMOVED before comparing -> a hip riding at +38 deg
    flexion scored like a level hip;
  - duty was the NEURAL RG-E duty, not real foot contact -> the left
    foot (0% contact frames) scored duty 0.69;
  - no phase or periodicity terms at all.

v2: BOTH legs, each cycle-detected from that foot's own CONTACT loading
onsets (runner now logs per-side heel+toe normal force; neural RG-E
onsets are the fallback), compared to that leg's own OpenSim reference
cycle on MEAN, AMPLITUDE (shape RMSE + excursion), PHASE (per-joint
cycle lag + interleg antiphase lag), and PERIODICITY (period error +
cycle-to-cycle variability), with both-foot contact guards and a
double-support (drag) term.

Conventions: OpenSim / ISB throughout. Joint angles are the OpenSim
coordinates (hip_flexion +flexion, knee_angle NEGATIVE = flexion,
ankle_angle +dorsiflexion); global axes X anterior, Y up, Z right —
figures must state this and use these signs.

kine_score: HIGHER IS BETTER, 0 = perfect. Guards make one-legged or
dragging gaits score below every genuine two-legged walker.
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
LOAD_N = 50.0        # N: loading threshold (same as the ref analysis)
CONTACT_MIN_FRAC = 0.15   # per-foot contact fraction required


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


def _loading_onsets(t, force):
    on = (force > LOAD_N).astype(int)
    return t[np.flatnonzero(np.diff(on) == 1)]


def load_reference():
    """Reference mean cycles for BOTH legs + duties, interleg lag, period.

    The right cycle runs right-loading-onset to next right-loading-onset
    (vertical GRF > 50 N, rising); the left cycle likewise from the left
    GRF column ('1_ground_force_vy' in subject01_walk1_grf.mot)."""
    t_ik, names, vals = _read_mot(IK_MOT)
    col = {n: i for i, n in enumerate(names[1:])}
    t_g, gn, gv = _read_mot(GRF_MOT)
    vy_r = gv[:, gn.index("ground_force_vy") - 1]
    vy_l = gv[:, gn.index("1_ground_force_vy") - 1]
    on_r = _loading_onsets(t_g, vy_r)
    on_l = _loading_onsets(t_g, vy_l)
    if len(on_r) < 2 or len(on_l) < 2:
        raise RuntimeError("no full gait cycle in the GRF reference")

    ref = {}
    for side, onsets, vcol, joints in (
            ("r", on_r, vy_r, ("hip_flexion_r", "knee_angle_r",
                               "ankle_angle_r")),
            ("l", on_l, vy_l, ("hip_flexion_l", "knee_angle_l",
                               "ankle_angle_l"))):
        t0, t1 = onsets[0], onsets[1]
        m = (t_ik >= t0) & (t_ik < t1)
        cyc = {j.split("_")[0]: np.interp(
            GRID, (t_ik[m] - t0) / (t1 - t0) * 100.0, vals[m, col[j]])
            for j in joints}
        ref[side] = cyc
        vcyc = np.interp(GRID, (t_g - t0) / (t1 - t0) * 100.0, vcol)
        ref[f"duty_{side}"] = float(np.mean(vcyc > LOAD_N))
        ref[f"T_{side}"] = float(t1 - t0)
        for j in ("hip", "knee", "ankle"):
            ref[f"{j}_range_{side}"] = float(np.ptp(cyc[j]))
        ref[f"knee_min_{side}"] = float(np.min(cyc["knee"]))
        ref[f"mean_{j}_{side}"] = None  # filled below

    for side in ("r", "l"):
        for j in ("hip", "knee", "ankle"):
            ref[f"mean_{j}_{side}"] = float(np.mean(ref[side][j]))

    # interleg lag: first left loading onset after the right cycle start,
    # as a fraction of the right period (contralateral heel strike)
    t_r0, T_r = on_r[0], ref["T_r"]
    after = on_l[on_l >= t_r0]
    ref["lag_rl"] = float((after[0] - t_r0) / T_r) if after.size else 0.5

    # reference double-support fraction (both vertical GRFs > 50 N) over
    # the right cycle
    i0 = int(np.searchsorted(t_g, on_r[0]))
    i1 = int(np.searchsorted(t_g, on_r[1]))
    ds = float(np.mean((vy_r[i0:i1] > LOAD_N) & (vy_l[i0:i1] > LOAD_N)))
    ref["ds"] = ds
    return ref


def _cycles_from_onsets(t, sig, onsets, cols, q_deg):
    """Mean cycle dict + periods for the columns `cols` of q_deg, cycles
    cut at `onsets` (seconds). Mirrors the v1 acceptance windows."""
    cycles = {j: [] for j in ("hip", "knee", "ankle")}
    periods = []
    for a, b in zip(onsets[:-1], onsets[1:]):
        dur = b - a
        if not (0.35 <= dur <= 2.5):
            continue
        ia = int(np.searchsorted(t, a))
        ib = int(np.searchsorted(t, b))
        if ib - ia < 8:
            continue
        ph = (t[ia:ib] - a) / dur * 100.0
        for j, c in zip(("hip", "knee", "ankle"), cols):
            cycles[j].append(np.interp(GRID, ph, q_deg[ia:ib, c]))
        periods.append(dur)
    if len(cycles["hip"]) < 2:
        return None, [], []
    return ({j: np.mean(v, axis=0) for j, v in cycles.items()},
            periods, cycles)


def sim_side(t, q_deg, neuro, walk_start, side, contact):
    """One leg's mean cycle from ITS OWN contact loading onsets (falls
    back to its RG-E burst onsets when contact is absent — flagged so
    the guard can still fail it). Returns (mean, periods, duty_contact,
    contact_frac, used_contact, n_cycles, onsets)."""
    m = t >= walk_start
    if not m.any():
        return None, [], float("nan"), 0.0, False, 0, np.array([])
    tt = t[m]
    si = 0 if side == "r" else 1
    rge_col = 2 if side == "r" else 4
    cols = (3, 4, 5) if side == "r" else (8, 9, 10)
    force = contact[m, si] if contact is not None else np.zeros(len(tt))
    cfrac = float(np.mean(force > 20.0))
    on = _loading_onsets(tt, force)
    used_c = True
    if len(on) < 3:
        # fallback: that side's RG-E burst onsets (contact-free cycle
        # detection; the contact guard still fails this trial)
        rge = neuro[m, rge_col]
        thr = 0.5 * max(np.max(rge), 1e-9)
        bon = (rge > thr).astype(int)
        on = tt[np.flatnonzero(np.diff(bon) == 1)]
        used_c = False
    mean, periods, _ = _cycles_from_onsets(t[m], tt, on, cols, q_deg[m])
    duty_c = float(np.mean(force > LOAD_N))
    return mean, periods, duty_c, cfrac, used_c, \
        (len(periods) if periods else 0), on


def _phase_lag(sim, refc):
    """Circular lag (percent of cycle, signed shortest) of the DC-removed
    sim cycle vs the ref cycle, by SSE-minimizing roll."""
    s = sim - np.mean(sim)
    r = refc - np.mean(refc)
    lags = np.arange(NPHASE)
    sse = [np.sum((s - np.roll(r, L)) ** 2) for L in lags]
    L = int(lags[int(np.argmin(sse))])
    if L > NPHASE // 2:
        L -= NPHASE
    return float(L)


def compare(t, q_deg, neuro, walk_start, ref=None, contact=None):
    """Both-leg comparison dict (kine_score higher = better) or None if
    neither leg yields cycles."""
    ref = ref or REF_CACHE
    if ref is None:
        ref = load_reference()
    legs = {}
    for side in ("r", "l"):
        legs[side] = sim_side(t, q_deg, neuro, walk_start, side, contact)
    if legs["r"][0] is None and legs["l"][0] is None:
        return None

    out = dict(contact_frac_r=legs["r"][3], contact_frac_l=legs["l"][3])
    total = 0.0
    W_SHAPE = dict(hip=1.0, knee=1.2, ankle=0.8)
    W_RANGE = dict(hip=0.5, knee=0.3, ankle=0.5)
    W_PHASE = dict(hip=0.25, knee=0.20, ankle=0.15)
    COLS = {"r": (3, 4, 5), "l": (8, 9, 10)}
    mw = t >= walk_start

    for side in ("r", "l"):
        mean, periods, duty_c, cfrac, used_c, n_cyc, on = legs[side]
        # ---- contact guards (Ben: both feet must actually walk) ----
        if cfrac < 0.05:
            total += 20.0          # foot never loads -> fail this leg
        elif cfrac < CONTACT_MIN_FRAC:
            total += 12.0
        if mean is None:
            # no cycles even on fallback onsets (2026-09-20 fix): score
            # this leg from its WALK-WINDOW statistics so a FROZEN leg
            # cannot dodge its joint penalties - the s3c round-1 winner's
            # left leg froze (duty 1.0, no cycles) and thereby scored
            # BETTER than a leg that cycles but misses would have.
            total += 8.0
            for j, c in zip(("hip", "knee", "ankle"), COLS[side]):
                seg = q_deg[mw, c]
                total += 0.3 * abs(float(np.mean(seg))
                                   - ref[f"mean_{j}_{side}"])
                total += W_RANGE[j] * abs(
                    float(np.ptp(seg)) - ref[f"{j}_range_{side}"])
            out[f"n_cycles{f'_{side}'}"] = 0
            out[f"frozen_{side}"] = True
            continue
        p = f"_{side}"
        out[f"n_cycles{p}"] = n_cyc
        out[f"duty{p}"] = duty_c
        out[f"duty{p}_ref"] = ref[f"duty{p}"]
        out[f"knee_min{p}"] = float(np.min(mean["knee"]))
        out[f"knee_min{p}_ref"] = ref[f"knee_min{p}"]
        for j in ("hip", "knee", "ankle"):
            s, r = mean[j], ref[side][j]
            # SHAPE (amplitude of the modulation)
            s0 = s - np.mean(s)
            r0 = r - np.mean(r)
            rmse = float(np.sqrt(np.mean((s0 - r0) ** 2)))
            out[f"rmse_{j}{p}"] = rmse
            out[f"range_{j}{p}"] = float(np.ptp(s))
            out[f"range_{j}{p}_ref"] = ref[f"{j}_range{p}"]
            total += W_SHAPE[j] * rmse
            total += W_RANGE[j] * abs(out[f"range_{j}{p}"]
                                      - out[f"range_{j}{p}_ref"])
            # MEAN (Ben's "mean": the DC level the v1 score removed)
            out[f"mean_{j}{p}"] = float(np.mean(s))
            out[f"mean_{j}{p}_ref"] = ref[f"mean_{j}{p}"]
            total += 0.3 * abs(out[f"mean_{j}{p}"] - out[f"mean_{j}{p}_ref"])
            # PHASE within the cycle
            total += W_PHASE[j] * min(abs(_phase_lag(s, r)), 50.0) / 50.0 \
                * 10.0
        total += 0.15 * abs(out[f"knee_min{p}"] - out[f"knee_min{p}_ref"])
        # PERIODICITY: period error + cycle-to-cycle variability
        if periods:
            T = float(np.mean(periods))
            out[f"T{p}"] = T
            out[f"T{p}_ref"] = ref[f"T{p}"]
            total += 4.0 * abs(T - ref[f"T{p}"]) / ref[f"T{p}"]
            if len(periods) > 2 and T > 0:
                cv = float(np.std(periods) / T)
                out[f"period_cv{p}"] = cv
                total += 4.0 * min(cv, 0.5)
        # duty from REAL contact
        total += 1.5 * abs(duty_c - ref[f"duty{p}"])

    # INTERLEG phase (antiphase): lag of left loading onsets vs right,
    # compared to the reference contralateral heel-strike lag
    on_r, on_l = legs["r"][6], legs["l"][6]
    if len(on_r) >= 2 and len(on_l) >= 1:
        T_r = float(np.mean(on_r[1:] - on_r[:-1])) if len(on_r) > 2 else \
            ref["T_r"]
        after = on_l[on_l >= on_r[0]]
        if after.size and T_r > 0:
            lag = float((after[0] - on_r[0]) / T_r) % 1.0
            out["lag_rl"] = lag
            out["lag_rl_ref"] = ref["lag_rl"]
            total += 6.0 * min(abs(lag - ref["lag_rl"]), 0.5)

    # DOUBLE SUPPORT / drag: both feet loaded simultaneously
    if contact is not None:
        m = t >= walk_start
        both = (contact[m, 0] > LOAD_N) & (contact[m, 1] > LOAD_N)
        ds = float(np.mean(both))
        out["ds"] = ds
        out["ds_ref"] = ref["ds"]
        total += 4.0 * min(abs(ds - ref["ds"]), 0.6)

    out["kine_score"] = float(-total)
    # legacy keys (right leg) so older readers keep working
    if "n_cycles_r" in out:
        out["n_cycles"] = out["n_cycles_r"]
        out["duty"] = out["duty_r"]
        out["duty_ref"] = out["duty_r_ref"]
        out["knee_min"] = out["knee_min_r"]
        out["knee_min_ref"] = out["knee_min_r_ref"]
        for j in ("hip", "knee", "ankle"):
            for k in ("rmse", "range"):
                out[f"{k}_{j}"] = out[f"{k}_{j}_r"]
            out[f"range_{j}_ref"] = out[f"range_{j}_r_ref"]
    return out


REF_CACHE = None


def ref_cached():
    global REF_CACHE
    if REF_CACHE is None:
        REF_CACHE = load_reference()
    return REF_CACHE


if __name__ == "__main__":
    r = ref_cached()
    print(f"reference: T_r {r['T_r']:.3f} s T_l {r['T_l']:.3f} s, "
          f"duty r/l {r['duty_r']:.2f}/{r['duty_l']:.2f}, "
          f"interleg lag {r['lag_rl']:.2f}, double support {r['ds']:.2f}")
    for side in ("r", "l"):
        print(f"  {side}: knee_min {r[f'knee_min_{side}']:.1f} deg, "
              f"ranges " + " ".join(
                  f"{j} {r[f'{j}_range_{side}']:.1f}" for j in
                  ("hip", "knee", "ankle")))
