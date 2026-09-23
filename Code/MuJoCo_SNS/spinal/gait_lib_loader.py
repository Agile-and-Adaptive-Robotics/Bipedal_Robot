"""Generalized gait-library loader (2026-09-23, goal-5 integration).

Turns any (IK motion .mot, vertical-GRF .mot) pair into the SAME reference
dict schema kine_ref.load_reference() produces (per-leg hip/knee/ankle
100-point phase cycles cut at loading onsets, duty/T/ranges/knee_min/
means/lag_rl/ds), so gait-library datasets (Falisse predictsim_mtp staged
at D:\\temp\\gait_lib_staging, Arnold muscfib_walkrun once Ben downloads)
drop straight in as additional reference cycles.

GRF column auto-detection: subject01 style ("ground_force_vy" right,
"1_ground_force_vy" left) or Falisse style ("r_/l_ground_force_vy").

New file only - kine_ref.py and the core four are untouched.
"""
from __future__ import annotations

import io
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import kine_ref as KR

GRF_STYLES = {
    "subject01": ("ground_force_vy", "1_ground_force_vy"),
    "falisse": ("r_ground_force_vy", "l_ground_force_vy"),
}
JOINTS = {"r": ("hip_flexion_r", "knee_angle_r", "ankle_angle_r"),
          "l": ("hip_flexion_l", "knee_angle_l", "ankle_angle_l")}


def _read_mot(path: Path):
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    end = next(i for i, ln in enumerate(lines) if ln.strip() == "endheader")
    hdr = " ".join(lines[:end]).lower()
    if "indegrees=no" in hdr:
        raise ValueError(f"{path.name}: inDegrees=no - radians not expected")
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


def _detect_grf(names):
    for style, (r, l) in GRF_STYLES.items():
        if r in names and l in names:
            return style, r, l
    raise ValueError(f"no known GRF vy column pair in {names[:8]}...")


def load_reference_general(ik_path: Path, grf_path: Path,
                           load_n: float = KR.LOAD_N):
    """kine_ref.load_reference() schema from any mot pair."""
    t_ik, names, vals = _read_mot(Path(ik_path))
    col = {n: i for i, n in enumerate(names[1:])}
    t_g, gn, gv = _read_mot(Path(grf_path))
    style, cr, cl = _detect_grf(gn)
    vy_r = gv[:, gn.index(cr) - 1]
    vy_l = gv[:, gn.index(cl) - 1]
    on_r = KR._loading_onsets(t_g, vy_r)
    on_l = KR._loading_onsets(t_g, vy_l)
    if len(on_r) < 2 or len(on_l) < 2:
        # Single-cycle periodic datasets (predictive sims are periodic by
        # construction): double the trace with period = file duration and
        # re-detect - each side then sees its onset twice, one cycle apart.
        T = float(t_g[-1] - t_g[0] + np.mean(np.diff(t_g)))
        t_g = np.concatenate([t_g, t_g + T])
        vy_r = np.concatenate([vy_r, vy_r])
        vy_l = np.concatenate([vy_l, vy_l])
        t_ik = np.concatenate([t_ik, t_ik + T])
        vals = np.vstack([vals, vals])
        on_r = KR._loading_onsets(t_g, vy_r)
        on_l = KR._loading_onsets(t_g, vy_l)
        ref_periodic = True
    else:
        ref_periodic = False
    if len(on_r) < 2 or len(on_l) < 2:
        raise RuntimeError(f"no full gait cycle in {Path(grf_path).name}")

    ref = {"grf_style": style, "periodic": ref_periodic}
    for side, onsets, vcol in (("r", on_r, vy_r), ("l", on_l, vy_l)):
        t0, t1 = onsets[0], onsets[1]
        m = (t_ik >= t0) & (t_ik < t1)
        cyc = {j.split("_")[0]: np.interp(
            KR.GRID, (t_ik[m] - t0) / (t1 - t0) * 100.0, vals[m, col[j]])
            for j in JOINTS[side]}
        ref[side] = cyc
        vcyc = np.interp(KR.GRID, (t_g - t0) / (t1 - t0) * 100.0, vcol)
        ref[f"duty_{side}"] = float(np.mean(vcyc > load_n))
        ref[f"T_{side}"] = float(t1 - t0)
        for j in ("hip", "knee", "ankle"):
            ref[f"{j}_range_{side}"] = float(np.ptp(cyc[j]))
            ref[f"mean_{j}_{side}"] = float(np.mean(cyc[j]))
        ref[f"knee_min_{side}"] = float(np.min(cyc["knee"]))
    t_r0, T_r = on_r[0], ref["T_r"]
    after = on_l[on_l >= t_r0]
    ref["lag_rl"] = float((after[0] - t_r0) / T_r) if after.size else 0.5
    i0 = int(np.searchsorted(t_g, on_r[0]))
    i1 = int(np.searchsorted(t_g, on_r[1]))
    ref["ds"] = float(np.mean((vy_r[i0:i1] > load_n) & (vy_l[i0:i1] > load_n)))
    return ref


def _stats_line(tag, ref):
    return (f"{tag}: style={ref.get('grf_style','subject01')} "
            f"T_r={ref['T_r']:.3f}s duty_r={ref['duty_r']:.2f} "
            f"duty_l={ref['duty_l']:.2f} lag_rl={ref['lag_rl']:.2f} "
            f"ds={ref['ds']:.2f} | knee_min r/l "
            f"{ref['knee_min_r']:.1f}/{ref['knee_min_l']:.1f} deg "
            f"hip_range r/l {ref['hip_range_r']:.1f}/{ref['hip_range_l']:.1f}")


def main():
    out = []
    # 1) regression: subject01 through the general loader must equal
    #    kine_ref.load_reference() (same file, same algorithm).
    mine = load_reference_general(KR.IK_MOT, KR.GRF_MOT)
    theirs = KR.load_reference()
    bad = []
    for side in ("r", "l"):
        for j in ("hip", "knee", "ankle"):
            d = float(np.max(np.abs(mine[side][j] - theirs[side][j])))
            if d > 1e-9:
                bad.append(f"{side}.{j}:{d:.2e}")
    for k in ("duty_r", "duty_l", "T_r", "T_l", "lag_rl", "ds",
              "knee_min_r", "knee_min_l"):
        if abs(mine[k] - theirs[k]) > 1e-9:
            bad.append(f"{k}:{mine[k]} vs {theirs[k]}")
    verdict = ("REGRESSION PASS (max dev < 1e-9 vs kine_ref)" if not bad
               else "REGRESSION FAIL: " + ", ".join(bad))
    print(verdict)
    out.append("## subject01 regression\n" + verdict)
    print(_stats_line("subject01 (kine_ref of record)", theirs))

    # 2) Falisse predicted walking (Case_40) as a second reference source.
    FAL = Path(r"D:\temp\gait_lib_staging\falisse\predictsim_mtp-master"
               r"\Results\Case_40")
    fal = load_reference_general(FAL / "motion.mot", FAL / "GRF.mot")
    print(_stats_line("Falisse Case_40 (predicted)", fal))
    out.append("\n## Falisse Case_40 (predicted walking)\n"
               + _stats_line("", fal))

    # 3) synergy cross-check: do OUR six per-leg synergies explain the
    #    Falisse predicted activations (independent dataset, same 43
    #    muscles/leg gait2392 lineage)?
    try:
        zb = np.load(HERE / "synergy_basis.npz", allow_pickle=True)
        print("synergy_basis keys:", sorted(zb.files))
        t_f, fn, fv = _read_mot(FAL / "motion.mot")
        for side, key in (("r", "W_r"), ("l", "W_l")):
            if key not in zb.files:
                continue
            W = zb[key]                       # (43 muscles, 6)
            names_key = f"muscle_names_{side}" if f"muscle_names_{side}" in zb.files \
                else ("muscle_names" if "muscle_names" in zb.files else None)
            if names_key is None:
                continue
            mnames = [str(s) for s in zb[names_key]]
            cols = [fn.index(f"{m}/activation") - 1 for m in mnames]
            A = fv[:, cols].T                 # (43, time)
            from scipy.optimize import nnls
            C = np.zeros((6, A.shape[1]))
            for k in range(A.shape[1]):
                C[:, k] = nnls(W, A[:, k])[0]
            R = W @ C
            vaf = 1.0 - np.sum((A - R) ** 2) / np.sum((A - A.mean(
                axis=1, keepdims=True)) ** 2)
            resid = 1.0 - np.sum((A - R) ** 2) / np.sum(A ** 2)
            print(f"  our W_{side} vs Falisse activations: "
                  f"VAF {vaf:.3f}, R2(vs 0) {resid:.3f}, "
                  f"mean co-act {C.mean():.3f}")
            out.append(f"- our W_{side} explains Falisse {side}-leg "
                       f"activations: VAF {vaf:.3f}, R2 {resid:.3f}")
    except Exception as e:
        print(f"(synergy cross-check skipped: {e})")
        out.append(f"(synergy cross-check skipped: {e})")

    dst = HERE / "reports_20260923" / "goal5_loader_smoke.txt"
    dst.parent.mkdir(exist_ok=True)
    dst.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"saved {dst}")


if __name__ == "__main__":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
    main()
