"""Introspect the compiled spinal network: dump EVERY connection with
sign and conductance, grouped by (src_kind, dst_kind).

This is the ground truth the Deng-style figure renders from — the
figure builder asserts its drawn edge groups against this inventory.
"""
import io
import sys
from collections import defaultdict

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import build_network as bn
import params as P
from muscle_map import classify

# representative pair per side incl. both v5/v6 topologies (gains > 0)
P.G["phase_reset_e"] = 0.22
P.G["phase_reset_f"] = 0.22
P.G["f1_kneext_inh"] = 0.59
ACTS = ["vas_lat_r", "semimem_r", "vas_lat_l", "semimem_l"]


def kind(name: str, sides=None) -> str:
    if name.startswith(("DRIVE", "POSTURE", "BAL")):
        return name
    for pre, k in (("PRESET_E", "PRESET_E_IN"), ("PRESET_F", "PRESET_F_IN"),
                   ("KINH", "KINH_IN"), ("IBEXC", "IBEXC_IN"),
                   ("ADAP", "ADAP_IN"), ("RG_E", "RG-E"), ("RG_F", "RG-F"),
                   ("PF_", "PF"), ("PFA", "PFA_IN"), ("MN_", "MN"),
                   ("Ia_", "Ia-aff"), ("II_", "II-aff"), ("Ib_", "Ib-aff")):
        if name.startswith(pre):
            return k
    return name


def main():
    net = bn.build(ACTS, interleg=True)
    n = net.net
    conns = getattr(n, "connections", None)
    print(f"neurons {len(net.idx)}, inputs {len(net.inputs)}, "
          f"connections attr type: {type(conns)}")
    if conns is None:
        print("attrs:", [a for a in dir(n) if not a.startswith("_")])
        return
    groups = defaultdict(list)
    pop_names = [p["name"] for p in n.populations]
    for c in conns:
        src = pop_names[c["source"]]
        dst = pop_names[c["destination"]]
        syn = c["params"]          # {'max_conductance': g, ...}
        g = float(syn.get("max_conductance", np.nan))
        er = float(syn.get("reversal_potential", np.nan))
        groups[(kind(src), kind(dst),
                "exc" if er > -1e-6 else "inh")].append((src, dst, g))
    print(f"total connections: {sum(len(v) for v in groups.values())}")
    for (sk, dk, sign), lst in sorted(groups.items()):
        g0 = lst[0][2]
        print(f"{sign.upper():3s} {sk:12s} -> {dk:12s}  x{len(lst):3d}"
              f"  g={g0:.3g}")
    # neuron name examples per kind
    print("\nname examples:")
    seen = set()
    for name in net.idx:
        k = kind(name, None)
        if k not in seen:
            seen.add(k)
            print(f"  {k:12s} {name}")


if __name__ == "__main__":
    main()
