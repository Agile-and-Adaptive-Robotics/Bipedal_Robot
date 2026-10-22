"""Forensic read of one coact-stand npz: what does the 'fall' actually
look like? COM/pelvis height trace, joint angles at key times, contact
forces, and when each DoF first leaves its post-release value."""
import sys
import numpy as np

f = sys.argv[1]
d = np.load(f)
t, com, q, qf = d["t"], d["com"], d["q"], d["qfull"]
ct, ctrl = d["contact"], d["ctrl"]
kj = [str(x) for x in d["key_joints"]]
n = min(int(np.sum(t > 0)) + 1, len(t))
t, com, q, qf, ct = t[:n], com[:n], q[:n], qf[:n], ct[:n]
print(f"file {f}: {n} steps, t_end {t[-1]:.2f}")
print("pelvis_ty adr guess: qfull col where? printing com + first qpos cols")
for tt in (0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.8, 2.1, 3.0, 4.0, 6.0, 7.9):
    k = int(np.searchsorted(t, tt))
    if k >= n:
        continue
    print(f"t={tt:4.1f} com=({com[k,0]:+.2f},{com[k,1]:+.2f},{com[k,2]:+.3f}) "
          f"ct=({ct[k,0]:6.0f},{ct[k,1]:6.0f})N "
          + " ".join(f"{kj[i]}={q[k,i]:+.0f}" for i in
                     (0, 3, 4, 5, 8, 9, 10)))
# pelvis translational dofs live at qpos[0:3] (free joint first)
print("\npelvis free-joint qpos [x y z qw qx qy qz] sample:")
for tt in (0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 7.9):
    k = int(np.searchsorted(t, tt))
    if k >= n:
        continue
    print(f"t={tt:4.1f} qpos0..6 = {np.round(qf[k, :7], 3)}")
# first divergence per key joint from t=0.7 s values
k0 = int(np.searchsorted(t, 0.7))
print("\nfirst |dq|>3 deg after t=0.7 (name, t, deg):")
rows = []
for i in range(q.shape[1]):
    dd = np.abs(q[k0:, i] - q[k0, i])
    w = np.flatnonzero(dd > 3.0)
    if w.size:
        rows.append((float(t[k0 + w[0]]), kj[i],
                     float(q[k0 + w[0], i] - q[k0, i])))
rows.sort()
for r in rows:
    print(f"  {r[0]:5.2f} s  {r[1]:18s} {r[2]:+7.1f}")
print("\nctrl mean over time (all 92):", round(float(ctrl[k0:].mean()), 3))
