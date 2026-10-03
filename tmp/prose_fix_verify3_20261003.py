import numpy as np
from scipy.io import loadmat

res = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"
m = loadmat(res + r"\Vas_Pam_20mm_Result.mat", squeeze_me=True, struct_as_record=False)
ctx = m["ctx"]
print("targetName:", ctx.targetName if hasattr(ctx, "targetName") else "(none)")
h = np.asarray(ctx.humanTorqueAbs)
print("extensor humanTorqueAbs: min %.4f max %.4f" % (h.min(), h.max()))
print("fBest:", m["fBest"], "exitflagP:", m["exitflagP"])
mm = np.asarray(m["torqueMarginFraction"]) if hasattr(m, "torqueMarginFraction") else None
print("min margin fraction:", mm.min() if mm is not None else "(field absent)")
