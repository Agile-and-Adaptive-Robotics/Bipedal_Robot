"""Dynamic audit of the phase-3 KINH pathway (F1 -> KINH -> knee_ext MN).

Builds two small networks (gain 0 = no KINH topology, gain 1.5 = KINH),
rhythms them at DRIVE=4 and compares a knee_ext MN pool (vas_lat) during
KINH-active windows. Expect: with KINH present, the knee_ext MN potential
is LOWER during PF-F1 (swing) bursts than without.
"""
from __future__ import annotations

import numpy as np

import build_network as bn
from params import DT, E_HI, G

AUDIT_ACTS = ["vas_lat_r", "semimem_r", "soleus_r", "psoas_r",
              "vas_lat_l", "semimem_l", "soleus_l", "psoas_l"]


def run(gain: float, dur: float = 16.0):
    G["f1_kneext_inh"] = gain
    net = bn.build(AUDIT_ACTS, dt=DT, interleg=True)
    u = net.make_inputs()
    u[net.input_index("DRIVE")] = 4.0
    u[net.input_index("POSTURE")] = 1.0
    n = int(dur / DT)
    names = ("PF_F1_r", "KINH_r", "MN_vas_lat_r")
    rec = np.zeros((n, len(names)))
    for k in range(n):
        v = net.step(u)
        rec[k] = (v[net.idx["PF_F1_r"]],
                  v[net.idx.get("KINH_r", -1)] if "KINH_r" in net.idx else 0.0,
                  v[net.idx["MN_vas_lat_r"]])
    tail = rec[n // 2:]
    f1 = tail[:, 0]
    on = f1 > 0.5 * max(f1.max(), 1e-9)
    swing = tail[on, 2]
    stance = tail[~on, 2]
    return dict(kinh_max=float(tail[:, 1].max()),
                mn_swing=float(swing.mean()), mn_stance=float(stance.mean()),
                mn_drop=float(stance.mean() - swing.mean()))


if __name__ == "__main__":
    a = run(0.0)
    print(f"gain 0.0 (no KINH): MN_vas_lat swing {a['mn_swing']:.3f} "
          f"stance {a['mn_stance']:.3f} drop {a['mn_drop']:.3f} "
          f"kinh_max {a['kinh_max']:.3f}")
    b = run(1.5)
    print(f"gain 1.5 (KINH)  : MN_vas_lat swing {b['mn_swing']:.3f} "
          f"stance {b['mn_stance']:.3f} drop {b['mn_drop']:.3f} "
          f"kinh_max {b['kinh_max']:.3f}")
    ok = b["kinh_max"] > 0.5 and b["mn_swing"] < a["mn_swing"] - 0.05
    print("KINH suppression:", "PASS" if ok else "FAIL/weak")
