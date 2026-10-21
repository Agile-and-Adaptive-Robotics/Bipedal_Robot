r"""Session scratch (tmp\, disposable): census + knob checks for the
2026-10-20 RA/SA + presynaptic edits to w2l_cpg\build_w2l_net.py."""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_cpg")

import numpy as np              # noqa: E402
import build_w2l_net as BW      # noqa: E402


def census(**kw):
    net = BW.build(template_key="bilateralrg2", **kw)
    return len(net.idx), net.n_synapses, len(net.inputs), net


c0 = census()
assert c0[:3] == (96, 204, 10), c0[:3]
c_ra = census(ra=True)
c_sa = census(sa=True)
c_both = census(ra=True, sa=True)
assert c_ra[:3] == (100, 224, 14), c_ra[:3]
assert c_sa[:3] == (100, 224, 14), c_sa[:3]
assert c_both[:3] == (104, 244, 18), c_both[:3]
print("census default / ra / sa / both:", c0[:3], c_ra[:3], c_sa[:3],
      c_both[:3])

net = c_both[3]
sc = dict(net.synapse_counts)
assert sc.get("ra_contact") == 20 and sc.get("sa_contact") == 20, sc
want_ports = sorted(f"{k} {s} {d}" for k in ("RA", "SA")
                    for s in ("heel", "toe") for d in ("L", "R"))
got_ports = sorted(p for p in net.inputs if p.startswith(("RA ", "SA ")))
assert got_ports == want_ports, got_ports
print("RA/SA ports:", got_ports)
# wiring check (when the toolbox exposes the connection list): every
# RA/SA channel must target exactly the 5 heel contact_hc extensor
# targets of its side
if hasattr(net.net, "connections"):
    from collections import Counter
    per_src = Counter()
    tgts_of = {}
    for cn in net.net.connections:
        pass
    # sns_toolbox Network.connections is list[Connection]; names resolved
    # through net.net.populations order = self.idx order
    names = list(net.idx)
    for cn in net.net.connections:
        pre = getattr(cn, "pre", None)
        post = getattr(cn, "post", None)
        if pre is None or post is None:
            per_src = None
            break
        tgts_of.setdefault(pre, []).append(post)
    if tgts_of and "RA heel L" in tgts_of:
        want = {"L RG ext", "L Hip PF ext", "L Knee PF ext",
                "L Hip MN ext", "L Knee MN ext"}
        assert set(tgts_of["RA heel L"]) == want, tgts_of["RA heel L"]
        assert set(tgts_of["SA toe R"]) == {"R RG ext", "R Hip PF ext",
                                            "R Knee PF ext", "R Hip MN ext",
                                            "R Knee MN ext"}, \
            tgts_of["SA toe R"]
        print("RA/SA wiring targets = heel contact_hc family: OK")
u = net.make_inputs()
u[net.input_index("RA heel L")] = 4.0
u[net.input_index("SA toe R")] = 4.0
V = net.step(u)
assert np.isfinite(V).all()
print("one-step finite with RA/SA drive: OK (V max %.4f mV)" % float(V.max()))
# tau check
assert BW.TAU["ra"] == 0.02 and BW.TAU["sa"] == 0.20
print("TAU ra/sa:", BW.TAU["ra"], BW.TAU["sa"])

# ---- presyn: default exact-identity + override takes effect
n_def = c0[3]
n_ovr = BW.build(template_key="bilateralrg2",
                 gains=dict(presyn_lo=0.5, presyn_hi=4.5,
                            contact_hc_presyn_lo=1.0,
                            contact_hc_presyn_hi=4.0))
gm = n_def.compiled.g_max_non != 0
assert (n_def.compiled.e_lo[gm] == 0.0).all() and \
       (n_def.compiled.e_hi[gm] == 5.0).all()
print(f"default presyn window on all {int(gm.sum())} synapse slots "
      f"exactly 0.0/5.0: OK (bit-identity)")
lo_arr = np.round(n_ovr.compiled.e_lo, 9)
hi_arr = np.round(n_ovr.compiled.e_hi, 9)
gmo = n_ovr.compiled.g_max_non != 0
print("override e_lo values present:",
      sorted(set(lo_arr[gmo].tolist())),
      "| e_hi values present:", sorted(set(hi_arr[gmo].tolist())))
assert 0.5 in lo_arr[gmo] and 1.0 in lo_arr[gmo]
assert 4.5 in hi_arr[gmo] and 4.0 in hi_arr[gmo]
# contact_hc edges (10) get 1.0/4.0; the rest 0.5/4.5
n_lo = int((lo_arr[gmo] == 1.0).sum())
n_hi = int((hi_arr[gmo] == 4.0).sum())
assert n_lo == 10 and n_hi == 10, (n_lo, n_hi)
print("contact_hc family (10 synapses) carries the family override: OK")

# old template still builds (shared-code regression)
n_old = BW.build()          # default key 'bilateralrg' (95/208 legacy)
assert (len(n_old.idx), n_old.n_synapses) == (95, 208)
print("legacy 'bilateralrg' default build 95/208: OK")
print("ALL CHECKS PASS")
