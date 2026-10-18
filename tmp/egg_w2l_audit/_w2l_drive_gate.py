r"""Pre/post BIT-EXACT PROBE for the L/R alpha+gamma drive knobs
(campaign 2026-10-16, egg). Run BEFORE the edit and AFTER (defaults):
identical stdout = bit-identical build + dynamics at default-off.

Probes (each: census + 4 s deterministic stepping + sha256 of the sampled
watch-traces, 3 repeats to prove in-process determinism):
  A. w2l_cpg bilateralrg: smoke-style drive (tonic 2/3, kickoff, 1 Hz
     alternating heel trains at 20 nA)
  B. w2l_mujoco M4 split: te/tf 3/4 + antiphase kickoff
  C. w2l_mujoco M5 aff net: split drive + scripted heel/toe pulses +
     constant Ib 3 nA (exercises every afferent family in the hash)

POST-EDIT ONLY (guarded by hasattr on the build signature): alpha/gamma
sanity — alpha_l=0.6 lifts the L MNs but not R; te asymmetry changes the
L/R RG-E burst counts.
"""
from __future__ import annotations

import hashlib
import inspect
import io
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = Path(__file__).resolve().parent
SPINAL = HERE.parent.parent
sys.path.insert(0, str(SPINAL))
sys.path.insert(0, str(SPINAL / "w2l_cpg"))
sys.path.insert(0, str(SPINAL / "w2l_mujoco"))

NSTEP = 2000          # 4 s at dt 2 ms
SAMPLE = 10           # hash every 10th step


def trace_hash(net, watch, drive, n=NSTEP):
    acc = []
    wi = [net.idx[w] for w in watch]
    for k in range(n):
        V = net.step(drive(net, k))
        if k % SAMPLE == 0:
            acc.append(V[wi].copy())
    A = np.asarray(acc, dtype=np.float64)
    return (hashlib.sha256(A.tobytes()).hexdigest(),
            " ".join(f"{x:.12f}" for x in A[-1]))


def probe_cpg(alpha_l=0.0, alpha_r=0.0):
    import build_w2l_net as M
    kwargs = dict()
    if "alpha_l" in inspect.signature(M.build).parameters:
        kwargs = dict(alpha_l=alpha_l, alpha_r=alpha_r)
    net = M.build(**kwargs)
    u = net.make_inputs()
    iS1 = net.input_index("Stimulus_1")
    iS2 = net.input_index("Stimulus_2 (10nA)")
    itE = [net.input_index("TONIC " + w) for w in ("L RG ext", "R RG ext")]
    itF = [net.input_index("TONIC " + w) for w in ("L RG flx", "R RG flx")]
    iH = (net.input_index("L heel contact"),
          net.input_index("R heel contact"))
    ai = []
    if alpha_l:
        ai += [(net.input_index(p), alpha_l) for p in net.inputs
               if p.startswith("ALPHA L ")]
    if alpha_r:
        ai += [(net.input_index(p), alpha_r) for p in net.inputs
               if p.startswith("ALPHA R ")]

    def drive(net_, k):
        t = k * 0.002
        u[:] = 0.0
        u[iS1] = 10.0 if t < 0.01 else 0.0
        u[iS2] = 10.0 if t < 0.01 else 0.0
        for i in itE:
            u[i] = 2.0
        for i in itF:
            u[i] = 3.0
        ph = t % 1.0
        u[iH[0]] = 20.0 if ph < 0.30 else 0.0
        u[iH[1]] = 20.0 if 0.5 <= ph < 0.8 else 0.0
        for i, val in ai:
            u[i] = val
        return u
    watch = ["L RG ext", "L RG flx", "R RG ext", "R RG flx"]
    return net, watch, drive


def probe_split(alpha_l=0.0, alpha_r=0.0):
    import build_w2l_split_net as M
    kwargs = {}
    if "alpha_l" in inspect.signature(M.build).parameters:
        kwargs = dict(alpha_l=alpha_l, alpha_r=alpha_r)
    net = M.build(**kwargs)
    return finish_split(net, alpha_l, alpha_r)


def probe_aff(alpha_l=0.0, alpha_r=0.0, te_l=3.0, te_r=3.0):
    import build_w2l_aff_net as M
    kwargs = {}
    if "alpha_l" in inspect.signature(M.build).parameters:
        kwargs = dict(alpha_l=alpha_l, alpha_r=alpha_r)
    net = M.build(**kwargs)
    return finish_aff(net, alpha_l, alpha_r, te_l, te_r)


def finish_split(net, alpha_l, alpha_r):
    u = net.make_inputs()
    iS1 = net.input_index("Stimulus_1")
    iS2 = net.input_index("Stimulus_2")
    itE = [net.input_index("TONIC " + w) for w in ("L RG ext", "R RG ext")]
    itF = [net.input_index("TONIC " + w) for w in ("L RG flx", "R RG flx")]
    ai_l = [net.input_index(p) for p in net.inputs
            if p.startswith("ALPHA L ")]
    ai_r = [net.input_index(p) for p in net.inputs
            if p.startswith("ALPHA R ")]

    def drive(net_, k):
        t = k * 0.002
        u[:] = 0.0
        u[iS1] = 10.0 if t < 0.01 else 0.0
        u[iS2] = 10.0 if t < 0.01 else 0.0
        for i in itE:
            u[i] = 3.0
        for i in itF:
            u[i] = 4.0
        for i in ai_l:
            u[i] = alpha_l
        for i in ai_r:
            u[i] = alpha_r
        return u
    return net, ["L RG ext", "R RG ext"], drive


def finish_aff(net, alpha_l, alpha_r, te_l=3.0, te_r=3.0):
    net, watch, drive = finish_split(net, alpha_l, alpha_r)
    # rebind with per-side te + afferent scripting on top
    ihe = [net.input_index(f"PORT heel {s}") for s in ("L", "R")]
    ito = [net.input_index(f"PORT toe {s}") for s in ("L", "R")]
    iib = [net.input_index(f"PORT Ib {s}") for s in ("L", "R")]
    itE = [net.input_index("TONIC " + w) for w in ("L RG ext", "R RG ext")]
    u = net.make_inputs()
    base = drive

    def drive2(net_, k):
        base(net_, k)
        t = k * 0.002
        u[itE[0]] = te_l
        u[itE[1]] = te_r
        ph = t % 1.2
        u[ihe[0]] = 4.0 if (t > 1.5 and ph < 0.4) else 0.0
        u[ito[0]] = 4.0 if (t > 1.5 and 0.7 < ph < 1.0) else 0.0
        u[ihe[1]] = 4.0 if (t > 1.5 and 0.6 <= ph < 1.0) else 0.0
        u[ito[1]] = 4.0 if (t > 1.5 and 0.1 < ph < 0.4) else 0.0
        u[iib[0]] = 3.0
        u[iib[1]] = 2.0
        return u
    return net, watch + ["L RG flx", "R RG flx"], drive2


def run_probe(label, factory):
    results = []
    for rep in range(3):
        net, watch, drive = factory()   # FRESH net each repeat: the
        results.append(trace_hash(net, watch, drive))  # compiled net is
        # stateful (V carries over) - a carried-over repeat is meaningless
    (h0, last0), (h1, _), (h2, _) = results
    print(f"[{label}] neurons={len(net.idx)} synapses={net.n_synapses} "
          f"inputs={len(net.inputs)}")
    print(f"   hash={h0}  repeats-equal={h0 == h1 == h2}  last={last0}")
    return net


def main():
    print("== PRE/POST DEFAULT PROBES (must be identical pre-edit vs "
          "post-edit-default) ==")
    run_probe("A bilateralrg", lambda: probe_cpg())
    run_probe("B split", lambda: probe_split())
    run_probe("C aff", lambda: probe_aff())

    import build_w2l_net as M
    if "alpha_l" not in inspect.signature(M.build).parameters:
        print("\n(post-edit sanity skipped: alpha knobs not present yet)")
        return 0

    print("\n== POST-EDIT SANITY (knobs ON; NOT expected to match anything) ==")
    net, watch, drive = probe_cpg(alpha_l=0.6)
    h, last = trace_hash(net, watch, drive)
    mns_l = [net.idx[net.muscle_outputs[a]] for a in net.muscle_outputs
             if a.startswith("L")]
    mns_r = [net.idx[net.muscle_outputs[a]] for a in net.muscle_outputs
             if a.startswith("R")]
    V = None
    u = net.make_inputs()
    for k in range(500):
        V = net.step(drive(net, k))
    print(f"[A alpha_l=0.6] hash={h[:16]}... L MN mean "
          f"{np.mean(V[mns_l]):.3f} mV vs R MN mean "
          f"{np.mean(V[mns_r]):.3f} mV (want L > R), finite="
          f"{np.isfinite(V).all()}")

    net, watch, drive = probe_aff(alpha_l=0.0, alpha_r=0.5)
    h, last = trace_hash(net, watch, drive)
    print(f"[C alpha_r=0.5] hash={h[:16]}... last={last} finite-check done")
    V = net.step(drive(net, 100))
    print(f"   finite={np.isfinite(V).all()}")

    net, watch, drive = probe_aff(te_l=3.0, te_r=1.0)
    # count RG-E bursts per side over 8 s to show asymmetry takes effect
    wi = {s: net.idx[f"{s} RG ext"] for s in ("L", "R")}
    tr = {s: np.zeros(4000) for s in ("L", "R")}
    for k in range(4000):
        V = net.step(drive(net, k))
        for s in wi:
            tr[s][k] = V[wi[s]]

    def bursts(sig):
        on = sig > 0.5 * sig.max()
        return int((np.flatnonzero(on[1:] & ~on[:-1]) + 1).size)
    print(f"[C te_l=3.0 te_r=1.0] RG-E burst starts L={bursts(tr['L'])} "
          f"R={bursts(tr['R'])} (asymmetry knob reachable)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
