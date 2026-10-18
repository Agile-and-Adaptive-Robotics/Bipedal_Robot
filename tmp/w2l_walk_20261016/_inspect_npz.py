import numpy as np
import sys

z = np.load(sys.argv[1])
t = z["t"]
print("len", len(t))
d = np.diff(t)
bad = np.flatnonzero(d < 0)
print("negative diffs:", len(bad), "first at idx", bad[:5] if len(bad) else "-")
for i in (0, 1000, 4490, 4500, 4600, 8000, 10000, 12000, 12250, 12500, 15000, 19999):
    print(f"t[{i}] = {t[i]:.4f}")
print("tail zeros from idx:", np.flatnonzero(t == 0.0).min() if (t == 0).any() else "none")
pel = z["pel"]
print("pel z: [0]", pel[0], "[4500]", pel[4500], "[10000]", pel[10000], "[19999]", pel[19999])
print("pel x: [0]", pel[0, 0], "[10000]", pel[10000, 0], "[19999]", pel[19999, 0])
rg = z["rgE"]
print("rgE L max", rg[:, 0].max(), "at idx", rg[:, 0].argmax())
print("rgE L [4500]", rg[4500, 0], "[10000]", rg[10000, 0])
print("rgE nonzero count L", int((np.abs(rg[:, 0]) > 1e-9).sum()),
      "R", int((np.abs(rg[:, 1]) > 1e-9).sum()))
heel = z["heelF"]
print("heelF L max", heel[:, 0].max(), "at t", heel[:, 0].argmax() * 0.001)
print("heelF L>20 count", int((heel[:, 0] > 20).sum()))
tilt = z["tilt"]
print("tilt max", tilt.max(), "at idx", tilt.argmax())
com = z["com"]
print("com x [0]", com[0, 0], "min", com[:, 0].min(), "max", com[:, 0].max())
print("com z [0]", com[0, 2], "min", com[:, 2].min(), "max", com[:, 2].max())
