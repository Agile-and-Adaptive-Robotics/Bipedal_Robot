"""Independent numerical verification of the 30-torque.tex math changes (2026-10-16).
Run with: C:\\Users\\Ben Bolen\\.conda\\envs\\myo\\python.exe tmp\\verify_torque_math_20261016.py
"""
import numpy as np

rng = np.random.default_rng(20261016)
res = {}


def reldev(A, B):
    A, B = np.asarray(A, float), np.asarray(B, float)
    return np.linalg.norm(A - B) / np.linalg.norm(B)


# --- 1. eq:M_vec / eq:momentarm_def vector identity: ||d x f|| = F ||d x u_hat|| ---
worst = 0.0
for _ in range(20000):
    d = rng.normal(size=3)
    u = rng.normal(size=3)
    u /= np.linalg.norm(u)
    F = rng.uniform(1, 900)
    f = F * u
    worst = max(worst, reldev(np.linalg.norm(np.cross(d, f)), F * np.linalg.norm(np.cross(d, u))))
res['1. ||d x f|| = F*||d x u_hat|| (max rel dev, 20000 draws)'] = worst

# --- 2. Planar case: r = ||d x u_hat|| coincides with every moment-arm reading ---
z = np.array([0.0, 0.0, 1.0])
worst = 0.0
for _ in range(20000):
    d2 = rng.normal(size=2)
    u2 = rng.normal(size=2)
    u2 /= np.linalg.norm(u2)
    d3 = np.array([*d2, 0.0])
    u3 = np.array([*u2, 0.0])
    r_cross = np.linalg.norm(np.cross(d3, u3))            # eq:momentarm_def
    r_axiscomp = abs(np.dot(np.cross(d3, u3), z))         # |z.(d x u)| = M_z/F
    r_ptline = np.linalg.norm(d3 - np.dot(d3, u3) * u3)   # point-to-line dist (joint center)
    n = np.cross(z, u3)
    r_axline = abs(np.dot(d3, n)) / np.linalg.norm(n)     # axis-line to line-of-action dist
    worst = max(worst, abs(r_cross - r_axiscomp), abs(r_cross - r_ptline), abs(r_cross - r_axline))
res['2. planar: r=||dxu||=|z.(dxu)|=pt-line=ax-line (max abs diff)'] = worst

# 3-D case: ||d x u_hat|| vs the axis-component arm -- documents planarity caveat
cnt, worst = 0, 0.0
for _ in range(20000):
    d3 = rng.normal(size=3)
    u3 = rng.normal(size=3)
    u3 /= np.linalg.norm(u3)
    r_cross = np.linalg.norm(np.cross(d3, u3))
    r_axis = abs(np.dot(np.cross(d3, u3), z))
    if r_cross > 1e-6:
        cnt += 1
        worst = max(worst, abs(r_cross - r_axis) / r_cross)
res['3. 3D: fraction where ||dxu|| != |z.(dxu)|'] = cnt / 20000
res['3. 3D: max rel gap ||dxu|| vs axis arm'] = worst

# eq:M_scalar in planar geometry: |M_z| = |z.(d x f)| = r F ; sign is a sensing convention
worst, pos = 0.0, 0
for _ in range(20000):
    d2 = rng.normal(size=2)
    u2 = rng.normal(size=2)
    u2 /= np.linalg.norm(u2)
    d3 = np.array([*d2, 0.0])
    f3 = F * np.array([*u2, 0.0])
    Mz = np.dot(np.cross(d3, f3), z)
    rF = np.linalg.norm(f3) * np.linalg.norm(np.cross(d3, np.array([*u2, 0.0])))
    worst = max(worst, reldev(abs(Mz), rF))
    pos += (Mz > 0)
res['4. planar |M_z| = r F (max rel dev)'] = worst
res['4. M_z sign positive fraction (sensing convention)'] = pos / 20000


# --- 5. y-symmetry: R_y K R_y^T = K for all theta iff kx = kz (30-torque.tex:135) ---
def Ry(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def Rz(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


ths = np.linspace(0, np.pi, 2001)
maxsym = 0.0
for _ in range(200):                       # symmetric family diag(a,b,a)
    a, b = 10 ** rng.uniform(4, 6), 10 ** rng.uniform(4, 6)
    K = np.diag([a, b, a])
    for t in ths:
        maxsym = max(maxsym, reldev(Ry(t) @ K @ Ry(t).T, K))
res['5. symmetric diag(a,b,a): max rel dev of R_y K R_y^T (200x2001)'] = maxsym

inv = 0
for _ in range(1000):                      # random asymmetric diagonals
    K = np.diag(10 ** rng.uniform(4, 6, 3))
    dev = max(reldev(Ry(t) @ K @ Ry(t).T, K) for t in ths)
    if dev < 1e-12:
        inv += 1
res['5. asymmetric random diagonals: invariant cases out of 1000'] = inv

# frames "differ by exactly a rotation about y_br": R2 = Rz(tz) Ry(-ty), R1 = Rz(tz)
worst = 0.0
for _ in range(2000):
    p = rng.normal(size=3)
    tz = np.arctan2(p[1], p[0])
    ty = rng.uniform(-1.4, 1.4)
    R1, R2 = Rz(tz), Rz(tz) @ Ry(-ty)
    worst = max(worst, reldev(R1.T @ R2, Ry(-ty)), np.linalg.norm(R1[:, 1] - R2[:, 1]))
res['5. two-rot vs pitch-only frame: relative rot about shared y (max dev)'] = worst

# stiffness seen by the two conventions: K1 = Q K2 Q^T with Q = Ry(-ty)
worst = 0.0
for _ in range(2000):
    a, b = 10 ** rng.uniform(4, 6, 2)
    K = np.diag([a, b, a])
    Q = Ry(-rng.uniform(-1.4, 1.4))
    worst = max(worst, reldev(Q @ K @ Q.T, K))
res['5. y-symmetric K identical under convention change (max rel dev)'] = worst

# --- 6. eq:equivalentstiffness series sum ---
worst = 0.0
for _ in range(2000):
    k1, k2 = 10 ** rng.uniform(4, 6, 2)
    F = rng.uniform(10, 500)
    x_total = F / k1 + F / k2                      # same force through both springs
    k_eq = (1 / k1 + 1 / k2) ** -1
    worst = max(worst, reldev(F / k_eq, x_total))
res['6. series k_eq = (1/k1+1/k2)^-1 (max rel dev)'] = worst

# --- 7. deflection chain eqs bktforcevector..keff with K = diag(chi1, chi2, chi1) ---
worst_c = worst_d = worst_k = worst_s = 0.0
chi1, chi2 = 3.998e4, 1.473e4
K = np.diag([chi1, chi2, chi1])
Kinv = np.linalg.inv(K)
for _ in range(5000):
    u = rng.normal(size=3)
    u /= np.linalg.norm(u)
    Fm = rng.uniform(10, 800)
    f_br = Fm * u                                   # eq:bktforcevector
    eps = Kinv @ f_br                               # eq:bktdeflection
    comp = np.array([Fm * u[0] / chi1, Fm * u[1] / chi2, Fm * u[2] / chi1])
    worst_c = max(worst_c, reldev(eps, comp))       # eq:bktcomponentdeflection
    d_br = u @ eps                                  # eq:bktprojecteddeflection
    worst_d = max(worst_d, reldev(d_br, Fm * (u @ Kinv @ u)))   # eq:bktprojectedcompliance
    worst_k = max(worst_k, reldev(Fm / d_br, (u @ Kinv @ u) ** -1))  # eq:keff
    k_ten = 10 ** rng.uniform(4, 6)
    worst_s = max(worst_s, reldev(1 / (1 / ((u @ Kinv @ u) ** -1) + 1 / k_ten),
                                  1 / (u @ Kinv @ u + 1 / k_ten)))
res['7. eq:bktcomponentdeflection = K^-1 f_br (max rel dev)'] = worst_c
res['7. eq:bktprojectedcompliance delta_br = F u^T K^-1 u (max rel dev)'] = worst_d
res['7. eq:keff F/delta_br = (u^T K^-1 u)^-1 (max rel dev)'] = worst_k
res['7. eq:equivalentstiffness series form (max rel dev)'] = worst_s

# --- 8. eq:deltaL / eq:comp: zeta = 1 - eps*, loss ~ zeta^2 ---
eps_star = np.linspace(0, 1, 11)
loss = (1 - eps_star) ** 2
res['8. loss at eps*=0 (resting length)'] = loss[0]
res['8. loss at eps*=1 (full contraction)'] = loss[-1]
res['8. loss monotone decreasing toward full contraction'] = bool(np.all(np.diff(loss) < 0))

for k, v in res.items():
    print(f'{k:<62s} {v!r}')
