"""M3 neural-only rhythm sweep: does the single-LH RG free-run endogenously
after the verbatim 10 nA / 10 ms kickoff (+ tonic drive), and at what period?
No MuJoCo. Sweeps tonic_e/tonic_f x tau_rg_nap_h.
"""
import io, os, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import importlib
import numpy as np

sys.path.insert(0, r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco")
sys.path.insert(0, r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
import build_w2l_orig_net as B
importlib.reload(B)

DT = 0.002
DUR = 24.0
SKIP = 3.0


def run(te, tf, tau, kick=10.0, kick_t=0.01):
    B.NAP["tau_max_h"] = tau
    net = B.build()
    u = net.make_inputs()
    iE = net.input_index("TONIC L RG ext")
    iF = net.input_index("TONIC L RG flx")
    iS = net.input_index("Stimulus_1")
    n = int(DUR / DT)
    e = np.zeros(n); f = np.zeros(n)
    for k in range(n):
        t = k * DT
        u[iE] = te; u[iF] = tf
        u[iS] = kick if t < kick_t else 0.0
        V = net.step(u)
        e[k] = V[net.idx["L RG ext"]]
        f[k] = V[net.idx["L RG flx"]]
    w = np.arange(n) * DT >= SKIP
    ew, fw = e[w], f[w]
    if ew.max() < 1e-6:
        return dict(flat=True, emax=ew.max(), fmax=fw.max())
    def bursts(sig):
        on = sig > 0.5 * sig.max()
        return np.flatnonzero(on[1:] & ~on[:-1]) + 1
    se, sf = bursts(ew), bursts(fw)
    r = float(np.corrcoef(ew, fw)[0, 1])
    per = float(np.diff(se).mean() * DT) if len(se) >= 3 else float("nan")
    return dict(flat=False, emax=round(ew.max(), 2), fmax=round(fw.max(), 2),
                bursts_e=len(se), bursts_f=len(sf), per=round(per, 3),
                r=round(r, 3), duty=round(float((ew > 0.5 * ew.max()).mean()), 2))


print("te  tf   tau   | eMax  fMax | burstsE burstsF | per(s) r     duty")
for tau in (0.25, 0.35):
    for te, tf in ((0.0, 0.0), (1.0, 2.0), (2.0, 3.0), (2.0, 4.0), (3.0, 4.0)):
        d = run(te, tf, tau)
        if d.get("flat"):
            print(f"{te:.1f} {tf:.1f} {tau:.2f} | FLAT eMax={d['emax']:.2g} fMax={d['fmax']:.2g}")
        else:
            print(f"{te:.1f} {tf:.1f} {tau:.2f} | {d['emax']:5.2f} {d['fmax']:5.2f} |"
                  f" {d['bursts_e']:7d} {d['bursts_f']:7d} | {d['per']:.3f} {d['r']:+.3f} {d['duty']:.2f}")
