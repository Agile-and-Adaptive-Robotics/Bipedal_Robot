import io, sys
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
p = r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\DataTool_7.txt'
d = np.genfromtxt(p, skip_header=1)
t = d[:, 1]
V_THR = -0.055


def episodes(v, t, thr=V_THR, gap=0.15):
    on = v > thr
    # merge: on-periods separated by < gap belong to one episode
    raw = []
    for i in range(1, len(on)):
        if on[i] and not on[i - 1]:
            raw.append([t[i], t[i]])
        elif on[i] and on[i - 1]:
            raw[-1][1] = t[i]
    merged = []
    for s, e in raw:
        if merged and s - merged[-1][1] < gap:
            merged[-1][1] = e
        else:
            merged.append([s, e])
    return [tuple(m) for m in merged]


names = dict(R_foot=2, R_toe=3, R_hipmid=4, L_toe=5, L_hipmid=6, L_foot=7)
eps = {}
for name, c in names.items():
    eps[name] = [x for x in episodes(d[:, c], t) if x[1] - x[0] > 0.05]
    print(f"{name}: {len(eps[name])} episodes")
    print("   ", [(round(s, 3), round(e, 3)) for s, e in eps[name][:8]])

# heel-contact (stance) period + L/R phase from the foot channels
for name in ("L_foot", "R_foot"):
    onsets = [s for s, e in eps[name]]
    if len(onsets) > 2:
        iv = np.diff(onsets)
        print(f"{name} stance onsets: {len(onsets)}, median interval "
              f"{np.median(iv):.3f} s => period {np.median(iv):.3f} s "
              f"({1/np.median(iv):.2f} Hz), durations median "
              f"{np.median([e-s for s,e in eps[name]]):.3f} s")

# L/R antiphase: overlap of stance episodes
lf, rf = eps["L_foot"], eps["R_foot"]
both = sum(1 for s1, e1 in lf for s2, e2 in rf
           if min(e1, e2) - max(s1, s2) > 0.02)
print(f"stance episodes with >20 ms L/R overlap: {both} "
      f"(of {len(lf)} L / {len(rf)} R)")

# duty factor from the hip-middle channels as swing markers (long episodes)
for name in ("L_hipmid", "R_hipmid"):
    onsets = [s for s, e in eps[name]]
    if len(onsets) > 2:
        print(f"{name} episodes: {len(onsets)}, median interval "
              f"{np.median(np.diff(onsets)):.3f} s")
