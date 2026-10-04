"""GATE 3 (SPIKING_MIRROR_PLAN.md): basin_gate.py perturbation protocol
on the SPIKING RG.

Adaptation of spinal/basin_gate.py (the 2026-09-13 banked insight: the
network sits near a bifurcation, so a winner that oscillates with only a
small basin will not survive a BLAS-order change or a Simulink port).
The original gate loads a RUNNER_DUMP_STATE winner json; the spiking
mirror has NO tuned winner yet, so the perturbed config is the same
baseline the rhythm gate uses: effective_tables('best') + the mirror's
SPIKE/CAL scalars.  Protocol otherwise identical: +/-1% multiplicative
jitter on every scalar parameter, network-only constant-DRIVE 20 s run,
PASS = last-5 s swing >= max(0.1 mV, 50% of the config's own baseline
first-window swing).

The spiking addition: the mirror's OWN scalars (SPIKE thresholds /
adaptation / bias + CAL calibration factors) are perturbed too - they
are exactly the parameters a port to AnimatLab/Simulink would re-round.

Usage: D:\\Anaconda\\envs\\myo\\python.exe basin_gate_spiking.py
        [--trials N] [--sigma S] [--drive D]
"""
import argparse
import io
import os
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

SPINAL = Path(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
os.chdir(SPINAL)
sys.path.insert(0, str(SPINAL))
assert "AARL_NET" not in os.environ

import numpy as np

import mujoco

import params
import build_network_spiking as bs
from draw_circuit import effective_tables

T = effective_tables("best")
BASE_G = dict(T["G"])
BASE_TAU = dict(T["TAU"])
BASE_W = {ph: dict(tbl) for ph, tbl in T["W"].items()}
BASE_WPOST = dict(T["WPOST"])
BASE_SPIKE = dict(bs.SPIKE)
BASE_CAL = dict(bs.CAL)


def apply_base():
    params.G.clear(); params.G.update(BASE_G)
    params.TAU.clear(); params.TAU.update(BASE_TAU)
    for ph in params.W_PF_MN:
        params.W_PF_MN[ph].clear()
        params.W_PF_MN[ph].update(BASE_W[ph])
    params.W_POSTURE.clear(); params.W_POSTURE.update(BASE_WPOST)
    bs.SPIKE.clear(); bs.SPIKE.update(BASE_SPIKE)
    bs.CAL.clear(); bs.CAL.update(BASE_CAL)


def perturb(sigma: float, rng: np.random.Generator) -> None:
    def jitter(x):
        return x * (1.0 + sigma * rng.uniform(-1.0, 1.0))

    for name in ("G", "TAU", "W_POSTURE"):
        tbl = getattr(params, name)
        for k in list(tbl):
            if isinstance(tbl[k], (int, float)):
                tbl[k] = jitter(float(tbl[k]))
    for ph in params.W_PF_MN:
        for g in list(params.W_PF_MN[ph]):
            params.W_PF_MN[ph][g] = jitter(float(params.W_PF_MN[ph][g]))
    for k in list(bs.SPIKE):
        if isinstance(bs.SPIKE[k], (int, float)):
            bs.SPIKE[k] = jitter(float(bs.SPIKE[k]))
    for k in list(bs.CAL):
        if isinstance(bs.CAL[k], (int, float)):
            bs.CAL[k] = jitter(float(bs.CAL[k]))


def run_once(net, drive, t_end):
    u = net.make_inputs()
    u[net.input_index("DRIVE")] = drive
    n = int(round(t_end / params.DT))
    v = np.empty(n + 1)
    for k in range(n):
        V = net.step(u)
        v[k + 1] = V[net.idx["RG_E_r"]] - V[net.idx["RG_F_r"]]
    return v


def window_swings(v, win_s=5.0):
    n_per = int(round(win_s / params.DT))
    n_win = int(len(v) // n_per)
    return np.array([v[k * n_per:(k + 1) * n_per].ptp()
                     for k in range(n_win)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=12)
    ap.add_argument("--sigma", type=float, default=0.01)
    ap.add_argument("--drive", type=float, default=2.5)
    ap.add_argument("--tend", type=float, default=20.0)
    args = ap.parse_args()

    m = mujoco.MjModel.from_xml_path(str(
        SPINAL.parents[2] / "Solid_Models" / "OpenSim" /
        "Gait2392_Robotbody" / "mjc" / "gait2392_simbody" /
        "gait2392_simbody_cvt3.xml"))
    acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(m.nu)]

    apply_base()
    net = bs.build(acts, dt=params.DT, interleg=True)
    base_v = run_once(net, args.drive, args.tend)
    base_win = window_swings(base_v)
    bar = max(0.1, 0.5 * float(base_win[0]))
    print(f"baseline windows (5 s swings, mV): "
          f"{np.round(base_win, 3)}  bar={bar:.3f} mV")

    passes, finals = 0, []
    for tr in range(args.trials):
        apply_base()
        perturb(args.sigma, np.random.default_rng(2000 + tr))
        net = bs.build(acts, dt=params.DT, interleg=True)
        v = run_once(net, args.drive, args.tend)
        if not np.isfinite(v).all() or float(np.abs(v).max()) > 1e4:
            print(f"  trial {tr:2d}: EXPLODED  FAIL")
            finals.append(float("nan"))
            continue
        win = window_swings(v)
        ok = win[-1] >= bar
        passes += int(ok)
        finals.append(float(win[-1]))
        print(f"  trial {tr:2d}: windows {np.round(win, 3)}  "
              f"{'PASS' if ok else 'FAIL'}")
    fin = np.array(finals)
    fin_ok = fin[np.isfinite(fin)]
    print(f"--> pass {passes}/{args.trials}; final swing min/med "
          f"{fin_ok.min() if fin_ok.size else float('nan'):.3f}/"
          f"{np.median(fin_ok) if fin_ok.size else float('nan'):.3f} mV; "
          f"margin {fin_ok.min() / base_win[-1] if fin_ok.size and base_win[-1] > 0 else float('nan'):.2f}x")
    out = {"sigma": args.sigma, "trials": args.trials,
           "drive": args.drive,
           "baseline_windows": [float(x) for x in base_win],
           "pass_fraction": passes / args.trials,
           "final_swing_min": float(fin_ok.min()) if fin_ok.size else None,
           "final_swing_median":
               float(np.median(fin_ok)) if fin_ok.size else None}
    import json
    (SPINAL / "reports_spiking_20261002" / "basin_gate_spiking.json"
     ).write_text(json.dumps(out, indent=1), "utf-8")
    ok_all = passes == args.trials
    print("GATE3 BASIN (spiking):",
          "PASS" if ok_all else
          f"PARTIAL/FAIL ({passes}/{args.trials})")
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
