import numpy as np

def Ry(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

rng = np.random.default_rng(0)
# Claim 1: diag(a,b,a) invariant under rotation about y
worst = 0.0
for _ in range(1000):
    a, b = rng.uniform(1, 1e6, 2)
    K = np.diag([a, b, a])
    t = rng.uniform(-np.pi, np.pi)
    worst = max(worst, np.max(np.abs(Ry(t) @ K @ Ry(t).T - K)) / np.max(np.abs(K)))
print("max rel. deviation, y-symmetric array:", worst)

# Claim 2: diag(a,b,c) with a != c is NOT invariant (deviation nonzero for generic theta)
count = 0
for _ in range(1000):
    a, b, c = rng.uniform(1, 1e6, 3)
    if np.isclose(a, c):
        continue
    K = np.diag([a, b, c])
    t = rng.uniform(0.2, np.pi / 2 - 0.2)
    dev = np.max(np.abs(Ry(t) @ K @ Ry(t).T - K)) / np.max(np.abs(K))
    if dev < 1e-12:
        count += 1
print("cases with ~zero deviation despite a!=c:", count, "of 1000 (expected 0)")
