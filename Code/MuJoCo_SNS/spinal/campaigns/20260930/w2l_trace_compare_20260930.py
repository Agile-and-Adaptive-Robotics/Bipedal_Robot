"""Simulink-vs-numpy trace comparison for the SNS_W2L_CPG validation.

Loads the numpy full-net reference (w2l_numpy_ref.npz) and the Simulink
representative-core run (SNS_Simscape\\results\\SNS_W2L_CPG_run_20260930.mat,
written by dev\\build_sns_w2l_cpg_20260930.m) and reports, for each of the
four RG half-centers over the 2-12 s window:
  - peak membrane potential (mV) both sides,
  - pointwise RMSE (mV) - expected to be phase-dominated (AGENTS.md chaos
    note: pointwise agreement between integrators dies at ~0.4 s),
  - best cross-correlation LAG of Simulink vs numpy (ms; negative = Simulink
    leads) within +/-200 ms,
plus the burst-onset offsets on L RG ext. Also writes the overlay figure
sns_w2l_traces.png (2-8 s).
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

HERE = Path(__file__).parent
SLX = Path(r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Code\Matlab"
           r"\SNS_Simscape\results\SNS_W2L_CPG_run_20260930.mat")

ref = np.load(HERE / "w2l_numpy_ref.npz", allow_pickle=True)
names = list(ref["names"])
t_np = ref["t"]
V_np = ref["V"]

from scipy.io import loadmat  # noqa: E402

mat = loadmat(SLX)
t_sl = mat["t"].ravel()
tr = {"L_RG_ext": mat["vLE"].ravel(), "L_RG_flx": mat["vLF"].ravel(),
      "R_RG_ext": mat["vRE"].ravel(), "R_RG_flx": mat["vRF"].ravel()}

SKIP = 2.0
wnp = t_np >= SKIP
wsl = t_sl >= SKIP

report = {}
NAME = {"L_RG_ext": "L RG ext", "L_RG_flx": "L RG flx",
        "R_RG_ext": "R RG ext", "R_RG_flx": "R RG flx"}
for nm, vsl in tr.items():
    j = names.index(NAME[nm])
    vnp = V_np[wnp, j]
    vsl_w = vsl[wsl]
    n = min(vnp.size, vsl_w.size)      # 2 ms grids, +-1 boundary sample
    vnp, vsl_w = vnp[:n], vsl_w[:n]
    rmse = float(np.sqrt(np.mean((vsl_w - vnp) ** 2)))
    # cross-correlation lag within +/-100 samples (200 ms)
    a = vsl_w - vsl_w.mean()
    b = vnp - vnp.mean()
    lags = range(-100, 101)
    cc = [float(np.corrcoef(a[max(0, -l):a.size - max(0, l)],
                            b[max(0, l):b.size - max(0, -l)])[0, 1])
          for l in lags]
    best = int(np.nanargmax(cc))
    lag_ms = (best - 100) * 2.0
    report[nm] = dict(peak_sim_mV=float(vsl_w.max()),
                      peak_np_mV=float(vnp.max()), rmse_mV=rmse,
                      xcorr_lag_ms=lag_ms, xcorr_r=cc[best])
    print(f"{nm:10s} peak sim {vsl_w.max():7.3f} / np {vnp.max():7.3f} mV | "
          f"RMSE {rmse:6.3f} mV | best xcorr r={cc[best]:.3f} at "
          f"lag {lag_ms:+.0f} ms (sim vs np)")

# burst onsets on L RG ext (smoke convention) for the onset-offset table
for tag, sig, t_ in (("sim", tr["L_RG_ext"][wsl], t_sl[wsl]),
                     ("np", V_np[wnp, names.index("L RG ext")], t_np[wnp])):
    on = sig > 0.5 * sig.max()
    st = t_[np.flatnonzero(on[1:] & ~on[:-1]) + 1]
    print(f"L RG ext burst onsets ({tag}):",
          " ".join(f"{x:.2f}" for x in st))
    report[f"onsets_{tag}"] = [float(x) for x in st]

(HERE / "w2l_trace_compare.json").write_text(json.dumps(report, indent=2),
                                             encoding="utf-8")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(4, 1, figsize=(9, 7), sharex=True)
    for a, nm in zip(ax, ["L_RG_ext", "L_RG_flx", "R_RG_ext", "R_RG_flx"]):
        m = (t_np >= 2) & (t_np <= 8)
        j = names.index(NAME[nm])
        a.plot(t_np[m], V_np[m, j], lw=1.0, label="numpy full net (95n)")
        a.plot(t_sl[(t_sl >= 2) & (t_sl <= 8)], tr[nm][(t_sl >= 2) & (t_sl <= 8)],
               lw=1.0, ls="--", label="Simulink core (26 cells)")
        a.set_ylabel(nm, fontsize=8)
        a.grid(True, alpha=0.3)
    ax[0].legend(fontsize=7, loc="upper right")
    ax[-1].set_xlabel("t [s]")
    fig.suptitle("SNS_W2L_CPG (Simulink representative core) vs numpy full "
                 "W2L net - same drive protocol")
    fig.tight_layout()
    fig.savefig(HERE / "sns_w2l_traces.png", dpi=130)
    print(f"figure: {HERE / 'sns_w2l_traces.png'}")
except Exception as e:  # figure must never flip the verdict
    print(f"(figure skipped: {e})")
