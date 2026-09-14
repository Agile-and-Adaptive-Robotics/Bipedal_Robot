"""Diagnose the fixed-tau Deng RG trace: peak times, dwell states, NaNs."""
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

import deng_cpg_ode as D

out, hh = D.simulate("fixed", t_end=20000.0)
tt = np.arange(out.shape[0]) * 0.1
v = out[:, 0]
print("NaNs:", int(np.isnan(out).sum()), " h range [%.2f, %.2f]"
      % (np.nanmin(hh), np.nanmax(hh)))
pk, _ = find_peaks(v, prominence=1.0)
print("HC_ext peak times (s):", np.round(pk * 0.1 / 1000, 2))
pk2, _ = find_peaks(-v, prominence=1.0)
print("HC_ext valley times (s):", np.round(pk2 * 0.1 / 1000, 2))
# histogram of V values -> which states are occupied
hist, edges = np.histogram(v, bins=24)
for i in range(24):
    print(f"  V {edges[i]:7.1f}: {'#' * int(60 * hist[i] / hist.max())}")
