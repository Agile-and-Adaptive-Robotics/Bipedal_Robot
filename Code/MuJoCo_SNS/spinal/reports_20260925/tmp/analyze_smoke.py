"""Analyze a w2lvar smoke npz: rhythm metrics from the neuro columns."""
import sys
import numpy as np

path = sys.argv[1] if len(sys.argv) > 1 else "w2lvar_smoke.npz"
d = np.load(path, allow_pickle=True)
neuro = d["neuro"]
names = [str(x) for x in d["neuro_names"]]
t = d["t"]
q = d["q"]
print("neuro columns:", names)
col = {n: i for i, n in enumerate(names)}

t0, t1 = 5.0, 15.0          # walk window
w = (t >= t0) & (t <= t1)
print(f"walk window {t0}-{t1} s: {w.sum()} samples")


def burst_count(x, frac=0.5):
    m = max(np.max(x), 1e-6)
    on = (x > frac * m).astype(int)
    return int(np.sum(np.diff(on) == 1)), float(np.mean(x > frac * m))


def corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    d = np.std(a) * np.std(b)
    return float(np.mean(a * b) / d) if d > 1e-12 else float("nan")


pairs = ["RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l",
         "PF_HIP-E_r", "PF_KNEE-E_r", "PF_KNEE-F_r",
         "PF_HIP-F_r" if "PF_HIP-F_r" in col else "PF_KNEE-F_r"]
for n in pairs:
    x = neuro[w, col[n]]
    bc, duty = burst_count(x)
    print(f"{n:14s} min {x.min():6.3f} max {x.max():6.3f} mean {x.mean():6.3f} "
          f"bursts {bc:3d} duty {duty:.2f}")

print()
print("correlations (walk window):")
for a, b in (("RG_E_r", "RG_F_r"), ("RG_E_l", "RG_F_l"),
             ("RG_E_r", "RG_E_l"), ("RG_F_r", "RG_F_l"),
             ("PF_HIP-E_r", "PF_KNEE-E_r"),
             ("PF_HIP-E_r", "PF_KNEE-F_r"),
             ("PF_KNEE-E_r", "PF_KNEE-F_r")):
    if a in col and b in col:
        print(f"  r({a:12s}, {b:12s}) = {corr(neuro[w, col[a]], neuro[w, col[b]]):+.3f}")

print()
print("joint ranges (walk window, deg):")
for j, i in (("hip_r", 3), ("knee_r", 4), ("ankle_r", 5),
             ("hip_l", 8), ("knee_l", 9), ("ankle_l", 10)):
    print(f"  {j:8s} {q[w, i].min():7.2f} .. {q[w, i].max():7.2f}")
