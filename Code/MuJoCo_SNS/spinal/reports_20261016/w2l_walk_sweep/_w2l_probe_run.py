r"""Instrumented probe of test_w2l_ground.run_phase for the drive sweep
(campaign 2026-10-16). Runs the walk phase with given CLI-style args,
then analyzes the RAW log (not the test's own window logic) to explain
the contradictory c-series metrics and to get strike times / antiphase /
progression for finalists.

Usage:
    python _w2l_probe_run.py <outname> --te=5 --tf=3 --alpha-l=0.1 ...
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco")

import numpy as np                     # noqa: E402

OUT = sys.argv[1]
ARGV = ["test_w2l_ground.py"] + sys.argv[2:]
sys.argv = ARGV
import test_w2l_ground as T            # noqa: E402  (parses KV from argv)

# ---- optional runtime leg-leg collision exclusion (--legcoll=0) -------
# The W2L aproj is a PLANAR walker: both legs share the sagittal plane and
# the aproj has NO leg-leg contact. The MJCF transport kept everything at
# contype=1/conaffinity=1, so the near-parallel tibia capsules overlap at
# spawn and MuJoCo wedges them apart with ~3.4 kN (measured,
# _w2l_static_pose.py). Re-mask so L<->R pairs are off, everything else
# unchanged: L ct=1/ca=1, R ct=2/ca=2, Root+ground ct=3/ca=3.
LEGCOLL = T.KV.get("legcoll", 1.0)
if LEGCOLL == 0:
    import mujoco
    _orig = mujoco.MjModel.from_xml_path

    L_GEOMS = ("femur_L", "tibia_L", "foot_L", "toe_L",
               "toe_L_contact", "foot_L_contact")
    R_GEOMS = ("femur_R", "tibia_R", "foot_R", "toe_R",
               "toe_R_contact", "foot_R_contact")

    def _patched(path, *a, **k):
        m = _orig(path, *a, **k)
        for names, ct in ((L_GEOMS, 1), (R_GEOMS, 2),
                          (("Root", "ground"), 3)):
            for nm in names:
                gid = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, nm)
                if gid >= 0:
                    m.geom_contype[gid] = ct
                    m.geom_conaffinity[gid] = ct
        return m

    mujoco.MjModel.from_xml_path = _patched
    print("[probe] leg-leg collision EXCLUDED (runtime contype remask)")

DT = T.DT_PHY


def strikes_of(f, thr=20.0, ref_s=0.3):
    on = f > thr
    st = np.flatnonzero(on[1:] & ~on[:-1]) + 1
    keep = []
    for x in st:
        if not keep or (x - keep[-1]) * DT > ref_s:
            keep.append(x)
    return np.array(keep)


def main() -> int:
    print("KV:", {k: T.KV[k] for k in ("te", "tf", "te_l", "te_r", "tf_l",
                                       "tf_r", "alpha_l", "alpha_r", "rig",
                                       "dur_walk", "wsup", "cap", "acap")})
    wlog = T.run_phase("walk", T.KV["dur_walk"], T.KV["rig"], None)
    t = wlog["t"]
    print(f"n={len(t)} t[0]={t[0]:.3f} t[-1]={t[-1]:.3f} "
          f"monotonic={bool(np.all(np.diff(t) >= 0))}")
    print(f"finite flag={wlog['finite']} fall_t={wlog['fall_t']}")
    nfin = int((~np.isfinite(t)).sum())
    print(f"non-finite t entries: {nfin}")
    q = wlog["q"]
    print(f"non-finite q rows: {int((~np.isfinite(q).all(axis=1)).sum())}")
    heelF, toeF = wlog["heelF"], wlog["toeF"]
    com, pel, tilt = wlog["com"], wlog["pel"], wlog["tilt"]
    rg = wlog["rgE"]
    t0 = T.KV["hold"] + T.KV["lower"] + 0.5
    w = t >= t0
    print(f"window t>={t0}: {int(w.sum())} samples, "
          f"span {t[w][-1] - t[w][0]:.3f} s" if w.sum() else
          f"window EMPTY (max t {t.max():.3f})")
    for nm, arr in (("heelF", heelF), ("toeF", toeF)):
        for k, s in enumerate(("L", "R")):
            print(f"{nm} {s}: max {arr[:, k].max():8.1f} N | "
                  f"mean {arr[:, k].mean():7.2f} | "
                  f">20N frac {float((arr[:, k] > 20).mean()):.3f} | "
                  f">5N frac {float((arr[:, k] > 5).mean()):.3f}")
    for k, s in enumerate(("L", "R")):
        stf = strikes_of(heelF[:, k])
        print(f"heel strikes {s} (FULL log): {len(stf)} at "
              f"{[round(float(x*DT),2) for x in stf]}")
        stw = strikes_of(heelF[w][:, k]) if w.sum() else np.array([])
        print(f"heel strikes {s} (window): {len(stw)} at "
              f"{[round(float(t[w][x]),2) for x in stw]}")
    # antiphase at heel level: phase of L strikes within R cycle
    stL = strikes_of(heelF[w][:, 0]) if w.sum() else np.array([])
    stR = strikes_of(heelF[w][:, 1]) if w.sum() else np.array([])
    if len(stL) >= 3 and len(stR) >= 3:
        tt = t[w]
        pl = np.median(np.diff(tt[stL])) if len(stL) > 2 else float("nan")
        pr = np.median(np.diff(tt[stR])) if len(stR) > 2 else float("nan")
        # L strike phase relative to R cycle
        ph = []
        for x in stR[:-1]:
            nx = stR[np.searchsorted(stR, x) + 1] if \
                np.searchsorted(stR, x) + 1 < len(stR) else None
            if nx is None:
                continue
            for y in stL:
                if x <= y <= nx:
                    ph.append((tt[y] - tt[x]) / (tt[nx] - tt[x]))
        print(f"strike intervals L {pl:.3f} s R {pr:.3f} s | "
              f"L-in-R-cycle phases {[round(p,2) for p in ph]}")
    print(f"RG-E: L max {rg[:,0].max():.2f} mV R max {rg[:,1].max():.2f} mV"
          f" | corr(full,dec) "
          f"{float(np.corrcoef(rg[::5,0], rg[::5,1])[0,1]):+.3f}")
    for lo in range(0, int(t[-1]), 2):
        m = (t >= lo) & (t < lo + 2)
        if m.sum():
            print(f"  t {lo:2d}-{lo+2:2d}s: rgE L {rg[m,0].mean():6.2f}"
                  f" R {rg[m,1].mean():6.2f} mV | com x {com[m,0].mean():+.3f}"
                  f" z {com[m,2].mean():.3f} | heelF L {heelF[m,0].mean():6.1f}"
                  f" R {heelF[m,1].mean():6.1f} | rigFz "
                  f"{wlog['rigFz'][m].mean():6.1f} N")
    print(f"com x first/last (window): {com[w][0,0]:+.4f} -> "
          f"{com[w][-1,0]:+.4f} (fwd {com[w][-1,0]-com[w][0,0]:+.4f} m)"
          if w.sum() else "no window")
    np.savez(os.path.join(HERE, OUT + ".npz"), **wlog)
    print("saved", OUT + ".npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
