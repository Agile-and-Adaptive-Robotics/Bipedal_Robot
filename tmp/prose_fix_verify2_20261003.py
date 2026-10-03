# Part 2 of verification: tauAbs census + torque field names (squeeze_me=True)
import numpy as np
from scipy.io import loadmat

res = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"

m = loadmat(res + r"\BiPulley_opensim_bifemsh_r_20260926_1702.mat", squeeze_me=True, struct_as_record=False)
st = m["st"]
print("st type:", type(st), "fields:", st._fieldnames if hasattr(st, "_fieldnames") else "-")
tau = np.asarray(st.tauAbs)
print("tauAbs shape", tau.shape, "total", tau.size, "NaN", int(np.isnan(tau).sum()),
      "finite", int(np.isfinite(tau).sum()), "max finite %.6f" % np.nanmax(tau))
print("needPair:", m.get("needPair"), " tauArmNominal:", m.get("tauArmNominal"),
      " tauTarget:", m.get("tauTarget"))
print("fBest:", m.get("fBest"), " exitflagP:", m.get("exitflagP"), " minPair:", m.get("minPair"),
      " pairName:", m.get("pairName"))

for f in ("Bifemsh_20mm_Result_pulley_20260926_1652.mat", "Bifemsh_20mm_Result.mat"):
    mm = loadmat(res + "\\" + f, squeeze_me=True, struct_as_record=False)
    keys = [k for k in mm if not k.startswith("__")]
    print("\n==", f, "top keys:", keys)
    def walk(prefix, v, depth=0):
        if depth > 2:
            return
        if hasattr(v, "_fieldnames"):
            for fn in v._fieldnames:
                arr = getattr(v, fn)
                if ("orque" in fn or "Tau" in fn) and not hasattr(arr, "_fieldnames"):
                    a = np.asarray(arr)
                    if a.dtype.kind in "fiu" and a.size:
                        print("   %s%s: shape %s min %.4f max %.4f NaN %d" %
                              (prefix, fn, a.shape, np.nanmin(a), np.nanmax(a), int(np.isnan(a).sum())))
                    else:
                        print("   %s%s: %s" % (prefix, fn, type(arr)))
                elif hasattr(arr, "_fieldnames"):
                    walk(prefix + fn + ".", arr, depth + 1)
    for k in keys:
        walk(k + ".", mm[k])
