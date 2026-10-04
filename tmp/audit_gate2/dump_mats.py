"""Audit: dump the worker's saved artifact numbers (before re-runs)."""
import numpy as np
from scipy.io import loadmat

REP = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_spiking_20261002"

print("=" * 30, "goal2_baselines.mat")
m = loadmat(REP + r"\goal2_baselines.mat", squeeze_me=True)
T = m["T"]
for row in T:
    print(" | ".join(str(x) for x in row))
k = m["k"]
print(f"knee demo: rise={k['rise']:.4f} settleMean={k['settleMean']:.3f} "
      f"settleMin={k['settleMin']:.2f} settleMax={k['settleMax']:.2f} "
      f"Ae={k['AeWin']:.3f} Af={k['AfWin']:.3f}")
cpg = m["cpg"]
print(f"cpg: thetaMin={cpg['thetaMin']:.2f} thetaMax={cpg['thetaMax']:.2f} nSw={cpg['nSw']}")
runs = m["runs"]
for r in runs:
    print(f"beer {r['label']}: sagAfter2={float(r['sagAfter2']):.3f} final={float(r['final']):.3f} "
          f"(timeseries A/A_t are MatlabOpaque; numbers in the T table row above)")

print("=" * 30, "goal2_knee_spiking.mat")
m = loadmat(REP + r"\goal2_knee_spiking.mat", squeeze_me=True)
S = m["S"]
for f in ["rise", "settleMean", "settleMin", "settleMax", "final", "AeWin", "AfWin",
          "rateIa", "rateIb", "VmeMean", "VmfMean", "altHz"]:
    print(f"  {f} = {S[f]}")

print("=" * 30, "goal2_beer_spiking.mat")
m = loadmat(REP + r"\goal2_beer_spiking.mat", squeeze_me=True)
print("  keys:", [k_ for k_ in m if not k_.startswith("__")])
for k_ in m:
    if k_.startswith("__"):
        continue
    v = m[k_]
    try:
        print(f"  {k_} = {v}")
    except Exception:
        print(f"  {k_}: (unprintable)")

print("=" * 30, "goal2_rgmn_spiking.mat")
m = loadmat(REP + r"\goal2_rgmn_spiking.mat", squeeze_me=True)
print("  keys:", [k_ for k_ in m if not k_.startswith("__")])
for k_ in m:
    if k_.startswith("__"):
        continue
    try:
        print(f"  {k_} = {m[k_]}")
    except Exception:
        print(f"  {k_}: (unprintable)")

print("=" * 30, "results/units_ref_spiking_sim.mat")
m = loadmat(r"D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\results\units_ref_spiking_sim.mat",
            squeeze_me=True)
tspk = np.atleast_1d(m["tspk"])
print(f"  sim spikes: {len(tspk)}, ISI mean {np.diff(tspk).mean()*1e3:.4f} ms, "
      f"vb end {np.asarray(m['vb_d']).ravel()[-1]:.5f} mV")
