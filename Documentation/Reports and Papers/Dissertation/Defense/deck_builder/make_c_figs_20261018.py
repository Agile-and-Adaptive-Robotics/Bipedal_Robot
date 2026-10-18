"""Regenerate batch-C deck figures with projection-legible fonts (QC round 1).

Reads the run artifacts (npz) directly and writes big-font PNGs into deck_assets/figs:
  c_bpa_ns_knees_big.png   / c_bpa_spk_knees_big.png  - BPA walker, both knees + RG-E
  c_w2lvar_joints_big.png  / c_syn6_joints_big.png    - right-leg hip/knee/ankle vs OpenSim ref

Knee columns in the BPA npz are identified from the run's own metrics JSON
(knee_rng_L/R); no column semantics are guessed.
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPINAL = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "deck_assets", "figs")
sys.path.insert(0, SPINAL)
from kine_ref import load_reference  # noqa: E402

SIM_C, REF_C = "#0072B2", "#E69F00"  # Tol blue / orange
TICK, LABEL, LEG = 17, 18, 15


def _style(ax):
    ax.tick_params(labelsize=TICK)
    ax.grid(True, alpha=0.3)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def bpa_knees(npz_name, out_name):
    d = np.load(os.path.join(SPINAL, npz_name), allow_pickle=True)
    q, t, rg = d["q"], d["t"], d["rgE"]
    m = json.loads(str(d["metrics"]))
    rng = q.max(axis=0) - q.min(axis=0)
    # column order is fixed by the campaign code (w2l_mujoco/_curriculum_bpa.py:59,230-232):
    # JNTS = (hip_L, knee_L, ankle_L, hip_R, knee_R, ankle_R); knee range asserted vs metrics
    kl, kr = 1, 4
    assert 55.0 < rng[kl] < 70.0 and 55.0 < rng[kr] < 70.0, (rng[kl], rng[kr])
    assert abs(rng[kl] - m["knee_rng_L"]) < 1.5 and abs(rng[kr] - m["knee_rng_R"]) < 1.5, (rng[kl], rng[kr])

    fig, axes = plt.subplots(3, 1, figsize=(4.8, 4.4), sharex=True)
    axes[0].plot(t, q[:, kl], color=SIM_C, lw=2.2)
    axes[0].set_ylabel("left knee\n(deg)", fontsize=LABEL)
    axes[1].plot(t, q[:, kr], color=SIM_C, lw=2.2)
    axes[1].set_ylabel("right knee\n(deg)", fontsize=LABEL)
    axes[2].plot(t, rg * 1e3, color=SIM_C, lw=2.2)
    axes[2].set_ylabel("RG-E\n(mV)", fontsize=LABEL)
    axes[2].set_xlabel("time (s)", fontsize=LABEL)
    for ax in axes:
        _style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, out_name), dpi=200)
    plt.close(fig)
    print("wrote", out_name, "| knee cols", kl, kr, "| ranges", round(rng[kl], 1), round(rng[kr], 1))


def joints_ref(npz_name, out_name, ref):
    d = np.load(os.path.join(SPINAL, npz_name), allow_pickle=True)
    q, t = d["q"], d["t"]
    kj = [str(x) for x in d["key_joints"]]
    cols = [("hip_flexion_r", "hip (deg)"), ("knee_angle_r", "knee (deg)"), ("ankle_angle_r", "ankle (deg)")]
    T_ref = ref.get("T_r", 1.2334)
    t0 = t - t[0]
    fig, axes = plt.subplots(3, 1, figsize=(4.8, 4.4), sharex=True)
    for ax, (cname, lab) in zip(axes, cols):
        jname = cname.split("_")[0]
        cyc = np.asarray(ref["r"][jname])
        n, dt = len(cyc), T_ref / len(cyc)
        ncyc = int(np.ceil(t0[-1] / T_ref)) + 1
        rt = np.concatenate([np.arange(n) * dt + c * T_ref for c in range(ncyc)])
        rv = np.concatenate([cyc] * ncyc)
        mask = rt <= t0[-1]
        ax.plot(t0, q[:, kj.index(cname)], color=SIM_C, lw=2.2, label="SNS model")
        ax.plot(rt[mask], rv[mask], color=REF_C, lw=1.6, ls="--", label="OpenSim ref")
        ax.set_ylabel(lab, fontsize=LABEL)
        _style(ax)
    axes[0].legend(fontsize=LEG, frameon=False, loc="upper right")
    axes[-1].set_xlabel("time (s)", fontsize=LABEL)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, out_name), dpi=200)
    plt.close(fig)
    print("wrote", out_name)


if __name__ == "__main__":
    bpa_knees("bpa_ns_best_20261003.npz", "c_bpa_ns_knees_big.png")
    bpa_knees("bpa_spk_best_20261003.npz", "c_bpa_spk_knees_big.png")
    ref = load_reference()
    joints_ref("spinal_run_w2lvar_s1_best_20261003.npz", "c_w2lvar_joints_big.png", ref)
    joints_ref("spinal_run_syn6_best_20261003.npz", "c_syn6_joints_big.png", ref)
