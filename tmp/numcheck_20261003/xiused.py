import scipy.io as sio
import numpy as np
R = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"
for f in ["Vas_Pam_20mm_Result.mat", "Bifemsh_20mm_Result.mat",
          "Bifemsh_20mm_Result_pulley_20260926_1652.mat"]:
    d = sio.loadmat(R + "\\" + f, squeeze_me=True, struct_as_record=False)
    print(f, "->")
    print("  keys:", sorted([k for k in d if not k.startswith("__")]))
    if "XiUsed" in d:
        print("  XiUsed =", np.asarray(d["XiUsed"]).ravel())
    # also hunt nested
    def walk(o, p="root", dep=0):
        if dep > 3: return
        if hasattr(o, "_fieldnames"):
            for fn in o._fieldnames:
                v = getattr(o, fn)
                if fn in ("XiUsed", "removedText", "removedNext"):
                    print("   nested %s.%s = %s" % (p, fn, np.asarray(v).ravel() if hasattr(v, "ravel") else v))
                walk(v, p + "." + fn, dep + 1)
        elif isinstance(o, np.ndarray) and o.dtype == object:
            for i, v in enumerate(o.ravel()[:4]):
                walk(v, "%s[%d]" % (p, i), dep + 1)
    walk(d.get("predBest", None), "predBest")
    walk(d.get("routeInfo", None), "routeInfo")
