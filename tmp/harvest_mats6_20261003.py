import scipy.io as sio
import numpy as np

f = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\BiPulley_opensim_bifemsh_r_20260926_1702.mat"
m = sio.loadmat(f, squeeze_me=True, struct_as_record=True)
st = m["st"]
for name in ("obj", "tauAbs"):
    v = st[name]
    print("st.%s type=%s" % (name, type(v).__name__))
    a = np.atleast_1d(np.asarray(v, dtype=object)).ravel()
    for i, x in enumerate(a):
        try:
            xx = np.asarray(x, dtype=float).ravel()
            print("  [%d] min=%.6g max=%.6g n=%d" % (i, xx.min(), xx.max(), xx.size))
        except Exception:
            print("  [%d] %r" % (i, x))
spec = m["spec"]
for name in ("name", "source", "notes", "fmaxSource", "pulley", "kin", "cross", "diameter", "frames"):
    if name in spec.dtype.names:
        print("spec.%s = %r" % (name, spec[name]))
print("predBest fields absent; T shape", np.asarray(m["T"]).shape)
