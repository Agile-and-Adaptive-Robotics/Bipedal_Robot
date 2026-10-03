import scipy.io as sio
import numpy as np

R = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"
d = sio.loadmat(R + r"\Vas_Pam_20mm_Result_20260925.mat", squeeze_me=True, struct_as_record=False)
ri = d["routeInfo"]
act = np.asarray(ri.active, dtype=float)          # expect (9, 100): points x frames
ea = np.asarray(ri.eliminatedAngleD, dtype=float).ravel()
ang = np.asarray(d["ctx"].phiD, dtype=float).ravel()
print("active shape:", act.shape)
print("removedText =", d["removedText"], " removedNext =", np.asarray(d["removedNext"]).ravel())
names = ["p1","p2","p3","p4","p5","p6","p7","p8","p9"]
for i in [2, 7]:  # p3, p8 (0-based)
    a = act[i] if act.shape[0] == 9 else act[:, i]
    n_active = int((a > 0).sum())
    print("%s: active frames = %d / %d ; eliminatedAngleD = %+.4f deg (rounds to %+.2f / 1dp %+.1f)" %
          (names[i], n_active, a.size, ea[i], ea[i], ea[i]))
    off_idx = np.where(a < 0.5)[0]
    if off_idx.size:
        print("   first off frame: 1-based idx %d at angle %+.4f deg" %
              (off_idx[0]+1, ang[off_idx[0]]))
        print("   last active frame: 1-based idx %d at angle %+.4f deg" %
              (off_idx[0], ang[off_idx[0]-1]))
# panel angle before p8 release = frame before first-off
a8 = act[7] if act.shape[0] == 9 else act[:, 7]
off8 = np.where(a8 < 0.5)[0]
print("p8 release-event panel (frame before release) angle = %+.4f deg" % ang[off8[0]-1])
# sanity: total frames, spacing
print("frames n=%d, angles %.3f..%.3f, spacing %.4f" % (ang.size, ang[0], ang[-1], ang[1]-ang[0]))
