import scipy.io as sio
import numpy as np
R = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"

d = sio.loadmat(R + r"\Vas_Pam_20mm_Result_20260925.mat", squeeze_me=True, struct_as_record=False)
pb = d["predBest"]
for f in ["pathLength0", "pathLength", "maxPathLength", "idxMaxPath", "maxPathAngleD",
          "maxPathLength0", "idxMaxPath0", "maxPathAngleD0", "rest", "tendon", "KMAX"]:
    if hasattr(pb, f):
        v = np.asarray(getattr(pb, f), dtype=float)
        print("predBest.%s = %s (n=%d)" % (f, np.array2string(v.ravel()[:6], precision=5), v.size))
ang = np.asarray(d["ctx"].humanAngleD, dtype=float).ravel()
pl = np.asarray(pb.pathLength, dtype=float).ravel()
i0 = int(np.argmin(np.abs(ang - 0.0)))
print("pathLength at theta=0 (idx %d, angle %.3f): %.4f m" % (i0, ang[i0], pl[i0]))
pl0 = np.asarray(pb.pathLength0, dtype=float).ravel()
print("pathLength0 at theta=0: %.4f m" % pl0[i0])

print()
for f in ["BiPulley_opensim_bifemsh_r_20260926_1702.mat", "BiPulley_opensim_med_gas_r_20260926_1655.mat"]:
    dd = sio.loadmat(R + "\\" + f, squeeze_me=True, struct_as_record=False)
    sp = dd["spec"]
    print(f)
    print("  spec.name =", sp.name, "| spec.fmax =", sp.fmax, "| fmaxSource =", sp.fmaxSource,
          "| footprintWidth =", sp.footprintWidth)
    print("  tauTarget =", dd["tauTarget"], " needPair =", dd["needPair"], " minPair =", dd["minPair"],
          " pairName =", dd["pairName"], " f0 =", dd["f0"], " fBest =", dd["fBest"])
    ta = np.asarray(dd["st"].tauAbs, dtype=float)
    print("  st.tauAbs: shape", ta.shape, "nan", int(np.isnan(ta).sum()), "/", ta.size,
          "max finite %.3f" % np.nanmax(ta))
