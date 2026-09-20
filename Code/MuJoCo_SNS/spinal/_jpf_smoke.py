"""Flag-on smoke for the joint-layer PF build (G["joint_pf"] = 1).

Builds, compiles (numpy), and steps the network-only circuit at 2 ms for
3 s with the same constant-DRIVE recipe as check_selfsustain.py; checks
finite dynamics and RG-E/RG-F alternation, and reports per-HC burst
amplitudes.  NOT a tuning result - a topology/plumbing check only.
"""
import io
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import build_network as bn
import params as P

P.G["joint_pf"] = 1.0
net = bn.build(["vas_lat_r", "semimem_r", "vas_lat_l", "semimem_l"],
               interleg=True)
names = [p["name"] for p in net.net.populations]
print(f"populations: {len(names)} (phase-cell build: 50)")
layer = [n for n in names if n.startswith("PF_HIP") or n.startswith("PF_KNEE")
         or n.startswith("PF_ANK")]
print(f"layer HCs: {sorted(layer)}")

be = net
DT = P.DT
u = net.make_inputs()
u[net.input_index("DRIVE")] = 2.5
steps = int(3.0 / DT)
hcs_avail = [f"PF_{hc}_r" for hc in bn.JPF_HCS
             if f"PF_{hc}_r" in net.idx]
rec = {nm: np.zeros(steps) for nm in
       ["RG_E_r", "RG_F_r"] + hcs_avail}
idx = {nm: net.idx[nm] for nm in rec}
for k in range(steps):
    V = net.step(u)
    for nm, i in idx.items():
        rec[nm][k] = V[i]
e = rec["RG_E_r"]
f = rec["RG_F_r"]
print(f"RG_E_r range {e.min():.2f}..{e.max():.2f} mV, "
      f"RG_F_r range {f.min():.2f}..{f.max():.2f} mV")
tail = slice(int(1.0 / DT), None)
em, fm = e[tail], f[tail]
alt = float(np.corrcoef(em, fm)[0, 1]) if em.std() > 1e-6 else 0.0
print(f"last-2s E/F antiphase correlation: {alt:.2f} "
      f"({'ALTERNATING' if alt < -0.5 else 'NOT alternating'})")
for nm in rec:
    if nm.startswith("PF_"):
        a = rec[nm]
        print(f"  {nm}: {a.min():.2f}..{a.max():.2f} mV "
              f"(last-1s swing {a[-int(1.0/DT):].std():.3f})")
assert np.isfinite(np.array(list(rec.values()))).all(), "NONFINITE dynamics"
print("SMOKE OK (finite; topology intact)")

