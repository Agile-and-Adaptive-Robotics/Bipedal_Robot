"""Compare the original (unscaled-replay) and retargeted backsolves."""
import io
import json
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

old = np.load("bsolve_out_unscaled_backup.npz", allow_pickle=True)
new = np.load("bsolve_out_retarget.npz", allow_pickle=True)

print("signs old:", json.dumps(json.loads(str(old["signs"]))))
print("signs new:", json.dumps(json.loads(str(new["signs"]))))

ro = np.asarray(old["res_norm"], float)
rn = np.asarray(new["res_norm"], float)
print(f"\nID/activation residual: old mean {ro.mean():.3f} -> "
      f"new mean {rn.mean():.3f}")

ao = np.asarray(old["acts"], float)
an = np.asarray(new["acts"], float)
names = [str(s) for s in old["act_names"]]

def peak_phase(a, t, hs):
    hs0, hs1 = hs[0], hs[1]
    ph = np.full(len(t), np.nan)
    m = (t >= hs0) & (t <= hs1)
    ph[m] = (t[m] - hs0) / (hs1 - hs0) * 100.0
    ok = np.isfinite(ph)
    if not ok.any() or a[ok].max() < 0.03:
        return float("nan")
    return float(ph[ok][int(np.argmax(a[ok]))])

hs = np.atleast_1d(np.asarray(old["hs_r"], float))
t = np.asarray(old["t"], float)
print("\nmuscle        mean-old  mean-new   peakphase old->new (right cycle)")
for mus in ("tib_ant_r", "soleus_r", "med_gas_r", "vas_lat_r",
            "rect_fem_r", "semimem_r", "iliacus_r", "glut_max1_r"):
    i = names.index(mus)
    po = peak_phase(ao[:, i], t, hs)
    pn = peak_phase(an[:, i], t, hs)
    print(f"{mus:<14s} {ao[:, i].mean():7.3f}  {an[:, i].mean():7.3f}   "
          f"{po:5.1f} -> {pn:5.1f}")

print(f"\ncorr(all acts, pooled): "
      f"{float(np.corrcoef(ao.ravel(), an.ravel())[0, 1]):.3f}")
