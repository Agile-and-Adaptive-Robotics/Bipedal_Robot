import io, sys
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
p = r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\DataTool_7.txt'
d = np.genfromtxt(p, skip_header=1)
t = d[:, 1]
cols = dict(R_foot=2, R_toe=3, R_hipmid=4, L_toe=5, L_hipmid=6, L_foot=7,
            height=8, distance=9)
V_THR = -0.055  # spike threshold in chart volts


def bursts(v, t, thr=V_THR):
    """Return (start, end) times of episodes where v > thr."""
    on = v > thr
    starts, ends = [], []
    for i in range(1, len(on)):
        if on[i] and not on[i - 1]:
            starts.append(t[i])
        if (not on[i]) and on[i - 1]:
            ends.append(t[i])
    if on[0]:
        starts.insert(0, t[0])
    if on[-1]:
        ends.append(t[-1])
    return list(zip(starts, ends))


# active window (tail zeros => sim recording artefact; find last t with nonzero)
nz = np.any(np.abs(d[:, 2:]) > 1e-12, axis=1)
t_end = t[np.where(nz)[0][-1]]
print(f"active data window: 0 .. {t_end:.3f} s  (rows {len(t)}, dt={t[1]-t[0]:.4f})")
print(f"height: start {d[0, 8]:.4f} min {d[:np.where(nz)[0][-1], 8].min():.4f} "
      f"max {d[:np.where(nz)[0][-1], 8].max():.4f} end(active) "
      f"{d[np.where(nz)[0][-1], 8]:.4f}")
print(f"distance: start {d[0, 9]:.4f} end(active) {d[np.where(nz)[0][-1], 9]:.4f} m"
      f"  => mean speed {(d[np.where(nz)[0][-1], 9]-d[0, 9]) / max(t_end, 1e-9):.4f} m/s")

for name, c in cols.items():
    if name in ("height", "distance"):
        continue
    b = bursts(d[:, c], t)
    # keep bursts in the active window
    b = [(s, e) for s, e in b if s < t_end]
    durs = [e - s for s, e in b]
    onsets = [s for s, e in b]
    iv = np.diff(onsets)
    print(f"\n{name}: {len(b)} bursts (thr {V_THR} V)")
    print(f"   onsets (s): {[round(x, 3) for x in onsets[:12]]}")
    if len(iv):
        print(f"   onset intervals: median {np.median(iv):.4f} s "
              f"mean {iv.mean():.4f} s  => rate {1/np.median(iv):.3f} Hz")
    if durs:
        print(f"   burst durations: median {np.median(durs):.4f} s")
