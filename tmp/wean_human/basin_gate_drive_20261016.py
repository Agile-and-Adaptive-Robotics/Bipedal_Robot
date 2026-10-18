"""Basin gate ADAPTED to the hand-drive winners (basin_gate.py pattern):
loads a wean round winner json via the campaign recipe (set_stage(4,
params) + knobs x scale), perturbs every scalar param +/-1% (basin_gate.
perturb), runs the network-only constant-drive rhythm 20 s feeding DRIVE
PLUS the GAMMA_c/ALPHA_c port currents exactly like the runner step loop
(runner.py:1366-1378), and requires the RG antiphase cycle to persist.

PASS rule (basin_gate.py): last-5s swing >= max(0.1 mV, half the
baseline's first-window swing); peaks reported.

Usage: python basin_gate_drive_20261016.py <winner.json> [--trials 12]
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import wean_tune_20261016 as WT
import basin_gate as BG
import params
import mujoco
import build_network as bn

ap = argparse.ArgumentParser()
ap.add_argument("winner")
ap.add_argument("--trials", type=int, default=12)
ap.add_argument("--sigma", type=float, default=0.01)
ap.add_argument("--tend", type=float, default=20.0)
ap.add_argument("--out", default=None)
args = ap.parse_args()

win = json.loads(open(args.winner, encoding="utf-8").read())
scale = float(win.get("scale", 1.0))
knobs = win["knobs"]

# apply the winner state into params (set_stage(4, base+fb) + knobs x
# scale + alpha_to_mn — fb injection required, see wean_tune apply_config)
fb = {k: float(win["params"].get(k, 0.0)) for k in WT.FB_KEYS}
WT.apply_config(knobs, scale, fb=fb)
P, R = WT.MODS["P"], WT.MODS["R"]
drive = WT.BASE_DRIVE

m = mujoco.MjModel.from_xml_path(str(
    os.path.join(os.path.dirname(HERE), "..", "..", "Solid_Models",
                 "OpenSim", "Gait2392_Robotbody", "mjc", "gait2392_simbody",
                 "gait2392_simbody_cvt3.xml")))
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
        for i in range(m.nu)]


def build_net():
    net = bn.build(acts, dt=params.DT, interleg=True)
    return net


def run_once(net):
    """Network-only rhythm with DRIVE + hand-drive ports (runner recipe)."""
    u = net.make_inputs()
    u[net.input_index("DRIVE")] = drive
    for s in getattr(net, "sides", ("l", "r")):
        for pre, gk in (("GAMMA_c", f"gamma_{s}"), ("ALPHA_c", f"alpha_{s}")):
            try:
                idx = net.input_index(f"{pre}_{s}")
            except Exception:
                continue
            u[idx] = drive * float(params.G[gk])
    n = int(round(args.tend / params.DT))
    v = np.empty(n + 1)
    ie, if_ = net.idx["RG_E_r"], net.idx["RG_F_r"]
    for k in range(n):
        V = net.step(u)
        v[k + 1] = V[ie] - V[if_]
    return v


net = build_net()
base_v = run_once(net)
base_win = BG.window_swings(base_v)
base_pk = BG.peaks_last10(base_v)
bar = max(0.1, 0.5 * float(base_win[0]))
print(f"=== {os.path.basename(args.winner)} scale {scale:g} drive {drive:.4f}")
print(f"baseline windows (5 s swings, mV): {np.round(base_win, 3)}  "
      f"peaks last10 {base_pk}")

# snapshot the applied state so each trial restores EXACTLY the winner
snapshot = {}
for name in ("G", "TAU", "W_POSTURE", "BAL", "PF_SHAPE", "W_PF_MN",
             "AFF", "MOD"):
    tbl = getattr(params, name)
    snapshot[name] = (json.dumps(tbl, default=lambda o: list(o)
                                 if isinstance(o, tuple) else o))


def restore():
    for name in ("G", "TAU", "W_POSTURE", "BAL", "PF_SHAPE", "W_PF_MN",
                 "AFF", "MOD"):
        tbl = getattr(params, name)
        tbl.clear()
        d = json.loads(snapshot[name])
        for k, v in d.items():
            if isinstance(v, list):
                v = tuple(v)
            elif isinstance(v, dict):
                v = {g: (tuple(w) if isinstance(w, list) else w)
                     for g, w in v.items()}
            tbl[k] = v


passes, final_swings = 0, []
for t in range(args.trials):
    restore()
    BG.perturb(args.sigma, np.random.default_rng(2000 + t))
    net = build_net()
    v = run_once(net)
    if not np.isfinite(v).all() or float(np.abs(v).max()) > 1e4:
        print(f"  trial {t:2d}: EXPLODED FAIL")
        passes += 0
        final_swings.append(float("nan"))
        continue
    win_sw = BG.window_swings(v)
    ok = win_sw[-1] >= bar
    passes += int(ok)
    final_swings.append(float(win_sw[-1]))
    print(f"  trial {t:2d}: windows {np.round(win_sw, 3)}  "
          f"{'PASS' if ok else 'FAIL'}")

fin = np.array(final_swings)
fin_ok = fin[np.isfinite(fin)]
report = dict(
    winner=os.path.abspath(args.winner), scale=scale, drive=drive,
    sigma=args.sigma, trials=args.trials,
    baseline_windows=[float(x) for x in base_win],
    baseline_peaks_last10=int(base_pk),
    pass_fraction=passes / args.trials,
    n_exploded=int(np.isnan(fin).sum()),
    final_swing_min=float(fin_ok.min()) if fin_ok.size else None,
    final_swing_median=float(np.median(fin_ok)) if fin_ok.size else None,
    margin_vs_baseline=float(fin_ok.min() / base_win[-1])
    if fin_ok.size and base_win[-1] > 0 else None,
)
out = args.out or os.path.join(
    HERE, "reports_20261016", "wean_human",
    f"basin_{os.path.basename(args.winner)}.json")
json.dump(report, open(out, "w", encoding="utf-8"), indent=1)
print(f"--> pass {passes}/{args.trials}  wrote {out}")
