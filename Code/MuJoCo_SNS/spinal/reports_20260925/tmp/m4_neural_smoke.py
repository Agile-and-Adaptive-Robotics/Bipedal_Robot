# M4 neural-only smoke: does the split net oscillate on both sides (coupled
# and ablated) before it goes into the MuJoCo gate?
import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco")
import numpy as np
import build_w2l_split_net as S

DT = 0.002
def run(comm, te=3.0, tf=4.0, dur=21.0, tau=0.25):
    S.NAP["tau_max_h"] = tau
    net = S.build(comm=comm)
    u = net.make_inputs()
    iE = [net.input_index("TONIC " + w) for w in ("L RG ext", "R RG ext")]
    iF = [net.input_index("TONIC " + w) for w in ("L RG flx", "R RG flx")]
    iS1 = net.input_index("Stimulus_1"); iS2 = net.input_index("Stimulus_2")
    n = int(dur / DT)
    lE = np.zeros(n); rE = np.zeros(n)
    for k in range(n):
        t = k * DT
        for i in iE: u[i] = te
        for i in iF: u[i] = tf
        u[iS1] = 10.0 if t < 0.01 else 0.0
        u[iS2] = 10.0 if t < 0.01 else 0.0
        V = net.step(u)
        if not np.isfinite(V).all():
            print(f"comm={comm}: NON-FINITE at t={t:.3f}"); return
        lE[k] = V[net.idx["L RG ext"]]; rE[k] = V[net.idx["R RG ext"]]
    win = slice(int(2.0 / DT), None)
    def bursts(sig, dt, gap=0.4):
        on = sig > 0.5 * sig.max()
        st = np.flatnonzero(on[1:] & ~on[:-1]) + 1
        if len(st) == 0: return st
        keep = [st[0]]
        for s in st[1:]:
            if (s - keep[-1]) * dt > gap: keep.append(s)
        return np.array(keep)
    bL = bursts(lE[win], DT); bR = bursts(rE[win], DT)
    pL = float(np.diff(bL).mean() * DT) if len(bL) >= 3 else float("nan")
    pR = float(np.diff(bR).mean() * DT) if len(bR) >= 3 else float("nan")
    r = float(np.corrcoef(lE[win], rE[win])[0, 1])
    print(f"comm={comm}: L bursts={len(bL)} period {pL:.3f}s | "
          f"R bursts={len(bR)} period {pR:.3f}s | r(L,RG-E,R RG-E)={r:+.3f} | "
          f"max L {lE[win].max():.2f} R {rE[win].max():.2f} mV")

run(comm=1.0)
run(comm=0.0)
