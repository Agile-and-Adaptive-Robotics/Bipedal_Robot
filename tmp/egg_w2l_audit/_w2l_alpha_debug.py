r"""Debug probe: verify the ALPHA ports exist and actually depolarize the
targeted side's MN pools (window-mean comparison, not single timepoint).
Run on egg post-edit."""
from __future__ import annotations

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


def run(net, nsteps, drive):
    Vm = np.zeros((nsteps, len(net.idx)))
    u = net.make_inputs()
    for k in range(nsteps):
        Vm[k] = net.step(drive(net, k, u))
    return Vm


def base_drive(net, te, tf, extra_alpha):
    iS1 = net.input_index("Stimulus_1")
    iS2 = net.input_index("Stimulus_2" if "Stimulus_2" in net.inputs
                          else "Stimulus_2 (10nA)")
    itE = [net.input_index("TONIC " + w) for w in ("L RG ext", "R RG ext")]
    itF = [net.input_index("TONIC " + w) for w in ("L RG flx", "R RG flx")]
    ai = extra_alpha

    def d(net_, k, u):
        t = k * 0.002
        u[:] = 0.0
        u[iS1] = 10.0 if t < 0.01 else 0.0
        u[iS2] = 10.0 if t < 0.01 else 0.0
        for i in itE:
            u[i] = te
        for i in itF:
            u[i] = tf
        for i, val in ai:
            u[i] = val
        return u
    return d


def mn_labels(net, side):
    # muscle label -> MN label; bilateralrg keys are MUSCLE node labels
    # ("L Hip ext"), aff/split keys are actuator names ("hip_L_ext") -
    # both side-tagged but in different positions
    return sorted(mn for mus, mn in net.muscle_outputs.items()
                  if mus.startswith(side) or f"_{side}_" in mus)


def report(tag, builder, alpha_l, alpha_r):
    net = builder(alpha_l=alpha_l, alpha_r=alpha_r)
    ai = ([(net.input_index(p), alpha_l)
           for p in net.inputs if p.startswith("ALPHA L ")] +
          [(net.input_index(p), alpha_r)
           for p in net.inputs if p.startswith("ALPHA R ")])
    print(f"[{tag}] alpha_l={alpha_l} alpha_r={alpha_r}: "
          f"neurons={len(net.idx)} synapses={net.n_synapses} "
          f"inputs={len(net.inputs)} alpha_ports={len(ai)}")
    if ai:
        print(f"   ports: {[p for p in net.inputs if p.startswith('ALPHA')]}")
    d = base_drive(net, 3.0, 4.0, ai)
    Vm = run(net, 4000, d)          # 8 s
    w = Vm[2000:]                   # window-mean over the last 4 s
    for side in ("L", "R"):
        labs = mn_labels(net, side)
        mean = float(np.mean([w[:, net.idx[m]] for m in labs]))
        print(f"   {side} MN window-mean (6 pools): {mean:.4f} mV")
    return net


print("signature check:",
      "alpha_l" in inspect.signature(__import__(
          "build_w2l_net").build).parameters)

import build_w2l_net as CPG
import build_w2l_aff_net as AFF

print("\n== A bilateralrg: off vs alpha_l=0.6 ==")
report("A-off", lambda **kw: CPG.build(**kw), 0.0, 0.0)
report("A-alphaL", lambda **kw: CPG.build(**kw), 0.6, 0.0)

print("\n== C aff net: off vs alpha_r=0.5 ==")
report("C-off", lambda **kw: AFF.build(**kw), 0.0, 0.0)
report("C-alphaR", lambda **kw: AFF.build(**kw), 0.0, 0.5)

print("\n== gamma asym reachability: aff net te_l=3 te_r=2 vs symmetric ==")
for te_l, te_r in ((3.0, 3.0), (3.0, 2.0), (3.0, 1.0)):
    net = AFF.build()
    ai = []
    d = base_drive(net, te_l, te_r, ai)
    # NOTE base_drive applies te to BOTH; override per side:
    itE = [net.input_index("TONIC " + w) for w in ("L RG ext", "R RG ext")]
    orig = d

    def d2(net_, k, u, orig=orig, itE=itE, te_l=te_l, te_r=te_r):
        orig(net_, k, u)
        u[itE[0]] = te_l
        u[itE[1]] = te_r
        return u
    Vm = run(net, 6000, d2)         # 12 s
    eL, eR = Vm[:, net.idx["L RG ext"]], Vm[:, net.idx["R RG ext"]]

    def bursts(sig):
        on = sig > 0.5 * sig.max()
        return np.flatnonzero(on[1:] & ~on[:-1]) + 1
    bL, bR = bursts(eL), bursts(eR)
    perL = np.diff(bL).mean() * 0.002 if len(bL) > 1 else float("nan")
    perR = np.diff(bR).mean() * 0.002 if len(bR) > 1 else float("nan")
    print(f"   te L/R {te_l}/{te_r}: L RG-E bursts {len(bL)} "
          f"(period {perL:.3f} s), R {len(bR)} (period {perR:.3f} s)")
print("\nDEBUG DONE")
