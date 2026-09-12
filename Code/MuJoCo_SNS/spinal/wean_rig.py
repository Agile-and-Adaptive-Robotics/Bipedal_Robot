"""Wean the support rig toward self-supported walking.

Runs the stand->walk->stand schedule at decreasing --rig-scale (every RIG
element: pelvis translations + rotations + lumbar/hip-rotation stabilizers;
stiffness *S, damping *sqrt(S)) with the balanced-pattern configuration
(--fitted --best). A stage PASSES when the run completes its full schedule
without NaN, COM height stays above 0.62 m, and pelvis tilt stays under
35 deg. Stops at the first failing stage and prints the ladder.

Usage: python wean_rig.py [scale2 scale3 ...]   (default 1.0 .6 .4 .25 .15)
Each stage's spinal_run.npz is preserved as wean_stage{i}_S{s}.npz.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np

import runner as R

HERE = Path(__file__).parent


def stage_metrics():
    d = np.load(HERE / "spinal_run.npz", allow_pickle=True)
    t, com, q = d["t"], d["com"], d["q"]       # q already in degrees
    n_done = int(np.sum(t > 0)) + 1
    n_done = min(n_done, len(t))
    full = t[n_done - 1] >= 0.95 * t[-1] and np.all(np.isfinite(q[:n_done]))
    kz = float(np.min(com[:n_done, 2]))
    tilt = float(np.max(q[:n_done, 0]))
    knee_min = float(np.min(q[:n_done, 4]))
    hip_amp = float(np.max(q[:n_done, 3]) - np.min(q[:n_done, 3]))
    dx = float(com[n_done - 1, 0] - com[int(np.searchsorted(t[:n_done], 5)), 0])
    return dict(full=bool(full), t_end=float(t[n_done - 1]), kz=kz,
                tilt=tilt, knee_min=knee_min, hip_amp=hip_amp, dx=dx)


def main(argv):
    stages = ([float(x) for x in argv] if argv
              else [1.0, 0.6, 0.4, 0.25, 0.15])
    rows = []
    for i, s in enumerate(stages):
        print(f"\n=== wean stage {i}: rig-scale {s} ===", flush=True)
        R.main(["--fitted", "--best", "--rig-scale", str(s)])
        m = stage_metrics()
        ok = m["full"] and m["kz"] > 0.62 and m["tilt"] < 35.0
        rows.append((s, m, ok))
        shutil.copy(HERE / "spinal_run.npz",
                    HERE / f"wean_stage{i}_S{s:g}.npz")
        print(f"stage {i} S={s}: {'PASS' if ok else 'FAIL'} "
              f"(t_end {m['t_end']:.1f} s, kz {m['kz']:.2f}, tilt "
              f"{m['tilt']:.1f} deg, knee {m['knee_min']:.1f} deg, "
              f"hip amp {m['hip_amp']:.1f} deg, dx {m['dx']:.2f} m)",
              flush=True)
        if not ok:
            break
    print("\n== weaning ladder ==")
    for s, m, ok in rows:
        print(f"  S={s:5g}: {'PASS' if ok else 'FAIL'} kz {m['kz']:.2f} "
              f"tilt {m['tilt']:5.1f} knee {m['knee_min']:6.1f} "
              f"hipamp {m['hip_amp']:5.1f} dx {m['dx']:+.2f}")


if __name__ == "__main__":
    main(sys.argv[1:])
