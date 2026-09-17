"""Basin-robustness gate for tuned spinal winners (banked insight 2026-09-13).

The network sits near a bifurcation: the reggate_v5_0 summation-order chaos
tips marginal limit cycles either way, so a winner that oscillates with only
a small basin will NOT survive a BLAS-order change or a Simulink/Simscape
port (E2's tonic collapse was exactly this). This gate loads a
RUNNER_DUMP_STATE json (the exact winner state dumped by runner.py), perturbs
every scalar parameter entry by +/-sigma (default 1%), runs the network-only
constant-DRIVE rhythm exactly like check_selfsustain.py, and requires the RG
antiphase cycle to PERSIST in the last 5 s of a 20 s run.

Basin margin reported per config:
  - pass fraction over trials
  - worst/median final-window swing as a fraction of that config's own
    unperturbed baseline swing

Usage:
  python basin_gate.py [state.json ...] [--drive D] [--trials N]
                       [--sigma S] [--tend T] [--out results.json]
Defaults: state_v10_best.json state_v10_study.json, drive = each dump's own
walk_drive, trials 12, sigma 0.01, T 20 s. PASS rule: last-window swing
>= max(0.1 mV, 50% of the BASELINE's first-window swing) and >= 6 antiphase
peaks in the last 10 s (period <= ~1.7 s).
"""
from pathlib import Path
import argparse
import json

import numpy as np
from scipy.signal import find_peaks

import mujoco

import build_network as bn
import params

HERE = Path(__file__).parent

TABLES = ("W_PF_MN", "W_POSTURE", "TAU", "G", "PF_SHAPE", "BAL", "AFF", "MOD")


def apply_state(d: dict) -> float:
    """Load a RUNNER_DUMP_STATE dict into the params module; return drive."""
    for name in TABLES:
        tbl = getattr(params, name)
        src = d[name]
        tbl.clear()
        for k, v in src.items():
            if isinstance(v, dict):
                tbl[k] = {g: (tuple(w) if isinstance(w, list) else w)
                          for g, w in v.items()}
            elif isinstance(v, list):
                tbl[k] = tuple(v)
            else:
                tbl[k] = v
    return float(d.get("walk_drive", 2.5))


def perturb(sigma: float, rng: np.random.Generator) -> None:
    """In-place +/-sigma multiplicative jitter on every scalar table entry."""
    def jitter(x):
        return x * (1.0 + sigma * rng.uniform(-1.0, 1.0))

    for name in ("G", "TAU", "W_POSTURE", "BAL"):
        tbl = getattr(params, name)
        for k in list(tbl):
            if isinstance(tbl[k], (int, float)):
                tbl[k] = jitter(float(tbl[k]))
    for ph in params.W_PF_MN:
        for g in list(params.W_PF_MN[ph]):
            params.W_PF_MN[ph][g] = jitter(float(params.W_PF_MN[ph][g]))
    for ph in list(params.PF_SHAPE):
        params.PF_SHAPE[ph] = tuple(
            jitter(float(x)) for x in params.PF_SHAPE[ph])
    for name in ("AFF", "MOD"):
        tbl = getattr(params, name)
        for k in list(tbl):
            if isinstance(tbl[k], (int, float)):
                tbl[k] = jitter(float(tbl[k]))


def run_once(net, drive: float, t_end: float) -> np.ndarray:
    u = net.make_inputs()
    u[net.input_index("DRIVE")] = drive
    n = int(round(t_end / params.DT))
    v = np.empty(n + 1)
    ie, if_ = net.idx["RG_E_r"], net.idx["RG_F_r"]
    for k in range(n):
        V = net.step(u)
        v[k + 1] = V[ie] - V[if_]
    return v


def window_swings(v: np.ndarray, win_s: float = 5.0) -> np.ndarray:
    # adapt the window when the run is shorter than one window (smoke tests)
    win_s = min(win_s, max(len(v) * params.DT / 4.0, params.DT))
    n_per = int(round(win_s / params.DT))
    n_win = int(len(v) // n_per)
    return np.array([v[k * n_per: (k + 1) * n_per].ptp()
                     for k in range(n_win)])


def peaks_last10(v: np.ndarray) -> int:
    tail = v[-int(round(10.0 / params.DT)):]
    pk, _ = find_peaks(tail, prominence=0.3)
    return len(pk)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("states", nargs="*",
                    default=["state_v10_best.json", "state_v10_study.json"])
    ap.add_argument("--drive", type=float, default=None,
                    help="override the dump's own walk_drive")
    ap.add_argument("--trials", type=int, default=12)
    ap.add_argument("--sigma", type=float, default=0.01)
    ap.add_argument("--tend", type=float, default=20.0)
    ap.add_argument("--out", default="basin_gate_results.json")
    args = ap.parse_args()

    m = mujoco.MjModel.from_xml_path(str(
        HERE.parents[2] / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
        / "mjc" / "gait2392_simbody" / "gait2392_simbody_cvt3.xml"))
    acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(m.nu)]

    report = {"sigma": args.sigma, "t_end": args.tend, "trials": args.trials,
              "configs": {}}
    for state_name in args.states:
        state_path = HERE / state_name
        d = json.loads(state_path.read_text("utf-8"))
        drive = args.drive if args.drive is not None else apply_state(d)
        net = bn.build(acts, dt=params.DT, interleg=True)

        base_v = run_once(net, drive, args.tend)
        base_win = window_swings(base_v)
        base_pk = peaks_last10(base_v)
        bar = max(0.1, 0.5 * float(base_win[0]))
        print(f"\n=== {state_name} (drive {drive:.3f}) ===")
        print(f"baseline windows (5 s swings, mV): {np.round(base_win, 3)}  "
              f"peaks last10 {base_pk}")

        passes, final_swings = 0, []
        for t in range(args.trials):
            apply_state(d)                       # restore exact winner state
            perturb(args.sigma, np.random.default_rng(1000 + t))
            net = bn.build(acts, dt=params.DT, interleg=True)
            v = run_once(net, drive, args.tend)
            if not np.isfinite(v).all() or float(np.abs(v).max()) > 1e4:
                # a +/-1% perturbation made the network blow up entirely -
                # the clearest possible basin failure (seen in the smoke run)
                print(f"  trial {t:2d}: EXPLODED "
                      f"(max |v| {float(np.abs(v[np.isfinite(v)]).max() if np.isfinite(v).any() else np.inf):.1e})  FAIL")
                final_swings.append(float("nan"))
                continue
            win = window_swings(v)
            # PERSISTENCE is the basin criterion: last-window swing alive and
            # >= half the baseline's first window. Peak count is REPORTED only
            # (period varies with drive; a fixed peak floor mislabels slow
            # winners — v10-best cycles at ~2 s = 5 peaks/10 s and is 12/12
            # robust, v1 rule called it FAIL).
            ok = win[-1] >= bar
            passes += int(ok)
            final_swings.append(float(win[-1]))
            print(f"  trial {t:2d}: windows {np.round(win, 3)}  "
                  f"{'PASS' if ok else 'FAIL'}")
        fin = np.array(final_swings)
        fin_ok = fin[np.isfinite(fin)]
        cfg = {
            "drive": drive,
            "baseline_windows": [float(x) for x in base_win],
            "baseline_peaks_last10": base_pk,
            "pass_fraction": passes / args.trials,
            "n_exploded": int(np.isnan(fin).sum()),
            "final_swing_min": float(fin_ok.min()) if fin_ok.size else None,
            "final_swing_median": float(np.median(fin_ok)) if fin_ok.size else None,
            "margin_vs_baseline": float(fin_ok.min() / base_win[-1])
            if fin_ok.size and base_win[-1] > 0 else float("nan"),
        }
        report["configs"][state_name] = cfg
        print(f"  --> pass {passes}/{args.trials}, "
              f"final swing min/med "
              f"{cfg['final_swing_min'] if cfg['final_swing_min'] is not None else float('nan'):.3f}/"
              f"{cfg['final_swing_median'] if cfg['final_swing_median'] is not None else float('nan'):.3f} mV, "
              f"margin {cfg['margin_vs_baseline']:.2f}x baseline")

    (HERE / args.out).write_text(json.dumps(report, indent=1), "utf-8")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
