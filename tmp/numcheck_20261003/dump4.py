import scipy.io as sio
import numpy as np
R = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"
d = sio.loadmat(R + r"\Bifemsh_20mm_Result_pulley_20260926_1652.mat",
                squeeze_me=True, struct_as_record=False)
xb = np.asarray(d["xBest"], dtype=float).ravel()
lb = np.asarray(d["ctx"].lb, dtype=float).ravel()
ub = np.asarray(d["ctx"].ub, dtype=float).ravel()
n = min(xb.size, lb.size)
print("xBest:", xb)
print("lb   :", lb)
print("ub   :", ub)
viol_lo = xb[:n] < lb[:n] - 1e-12
viol_hi = xb[:n] > ub[:n] + 1e-12
print("bound violations at xBest:", int(viol_lo.sum() + viol_hi.sum()))
pb = d["predBest"]
print("predBest.infeasibleCount =", getattr(pb, "infeasibleCount", None))
print("predBest.slackCount =", getattr(pb, "slackCount", None))
print("predBest.failReason =", getattr(pb, "failReason", None))
print("predBest.ok =", getattr(pb, "ok", None))
# collision info min clearance
ci = getattr(d.get("collisionInfo", None), "minClearance", None)
print("collisionInfo.minClearance =", ci if ci is None else np.asarray(ci).ravel()[:3])
mm = getattr(d.get("collisionInfo", None), "minRouteSeparation", None)
print("collisionInfo.minRouteSeparation =", mm if mm is None else np.asarray(mm).ravel()[:3])
rc = getattr(d.get("collisionInfo", None), "routeClearanceRequired", None)
print("collisionInfo.routeClearanceRequired =", rc if rc is None else np.asarray(rc).ravel()[:3])
