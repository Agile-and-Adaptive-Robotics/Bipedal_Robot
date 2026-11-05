"""Validate the extracted cam geometry (READ-ONLY check of the dumped design of record).

Verifies, per pose:
  1. pivot (ICR) in femur frame = T_Pam(1:3,4)
  2. moment arm = perpendicular distance from pivot to the tibia-side line of
     action (last active tibial contact -> p9), compared with the evaluator's
     predBest.momentArm.
  3. route transforms reproduce the figure-frame geometry (femur frame).
"""
import json
import numpy as np

D = json.load(open(r"C:/Users/Ben/Documents/GitHub/Bipedal_Robot/tmp/cam_visuals_dump_20261104.json"))

route = np.array(D["route_raw"])          # 9 x 3 x N
active = np.array(D["route_active"])      # 9 x N
T_Pam = np.array(D["T_Pam"])              # 4 x 4 x N
T_ICR_t1 = np.array(D["T_ICR_t1"])
phiD = np.array(D["phiD"])
mA_eval = np.array(D["momentArm_m"])
N = D["n_poses"]


def rvt(T, p):
    """RowVecTrans equivalent: row-vector transform."""
    p = np.asarray(p, float)
    return p @ T[:3, :3].T + T[:3, 3]


def to_femur(p_t1, ii):
    return rvt(T_Pam[:, :, ii], rvt(T_ICR_t1[:, :, ii], p_t1))


print("phiD     mA_eval   mA_evalrep    mA_geom(ins)  wrap_r_mm  contacts")
errs = []
errs2 = []
T_Pam_inv = np.array(D["T_Pam"])
# compute inverse per pose
Tinv = np.linalg.inv(np.transpose(T_Pam[:3, :3, :], (2, 0, 1)))
for ii in range(0, N, 10):
    Tp = T_Pam[:, :, ii]
    Tpi = np.eye(4)
    Tpi[:3, :3] = Tp[:3, :3].T
    Tpi[:3, 3] = -Tp[:3, :3].T @ Tp[:3, 3]
    # evaluator replication: ICR frame, p6 line of action toward p5
    p5_icr = rvt(Tpi, route[4, :, ii])       # p5 femur -> ICR
    p6_icr = rvt(T_ICR_t1[:, :, ii], route[5, :, ii])  # p6 t1 -> ICR
    u = p5_icr - p6_icr
    u = u / np.linalg.norm(u)
    mA_rep = float(np.linalg.norm(p6_icr - np.dot(p6_icr, u) * u))
    # teaching insertion-line distance (femur frame)
    pivot = Tp[:3, 3]
    act = active[:, ii].astype(bool)
    js = [j for j in (8, 7, 6) if act[j - 1]]
    j = js[0] if js else 6
    p_last = to_femur(route[j - 1, :, ii], ii)
    p9 = to_femur(route[8, :, ii], ii)
    v = p9 - p_last
    v = v / np.linalg.norm(v)
    w = pivot - p_last
    mA_geom = float(np.linalg.norm(w - np.dot(w, v) * v))
    wrf = np.array(D["wrap_radius_full"])   # 4 x N: [p2 femur cyl; p6-8 tibia; p3-5 femur; p5-6 ICR]
    wr = wrf[1, ii]
    print(f"{phiD[ii]:+7.2f}  {mA_eval[ii]:.4f}   {mA_rep:.4f}      {mA_geom:.4f}      {wr:6.1f}    "
          + "".join("1" if a else "0" for a in active[5:9, ii].astype(int)))
    errs.append(mA_geom - mA_eval[ii])
    errs2.append(mA_rep - mA_eval[ii])

errs = np.array(errs)
errs2 = np.array(errs2)
print(f"\ninsertion-line: mean diff {errs.mean()*1000:+.1f} mm   std {errs.std()*1000:.1f} mm")
print(f"p6-line replica: mean diff {errs2.mean()*1000:+.3f} mm   max abs {np.abs(errs2).max()*1000:.3f} mm")
print(f"mA_eval range: {mA_eval.min()*1000:.1f} .. {mA_eval.max()*1000:.1f} mm  (record: 59.4-102.8)")

# pivot migration range
pivots = T_Pam[:3, 3, :].T
print(f"pivot x range {pivots[:,0].min()*1000:.1f}..{pivots[:,0].max()*1000:.1f} mm, "
      f"y range {pivots[:,1].min()*1000:.1f}..{pivots[:,1].max()*1000:.1f} mm (femur frame)")

# wrap radius structure
wrf = D["wrap_radius_full"]
print("wrap_radius_full shape:", np.array(wrf).shape)
# geo key numbers
g = D["geo"]
print("femurCylCenter", np.array(g["femurCylCenter"]), "femurCylRadius(m)", g["femurCylRadius"],
      "clear", g["femurCylClearRadius"])
print("tibiaLowerCenter", np.array(g["tibiaLowerCenter"]), "r", g["tibiaLowerRadius"],
      "clear", g["tibiaLowerClearRadius"])
print("tibiaUpperCenter", np.array(g["tibiaUpperCenter"]), "r", g["tibiaUpperRadius"],
      "clear", g["tibiaUpperClearRadius"])
print("seedP6", np.array(g["seedP6XY"]), "seedP7", np.array(g["seedP7XY"]), "seedP8", np.array(g["seedP8XY"]))
print("femurProfileCenter", np.array(g["femurProfileCenter"]))
print("Xi", D.get("Xi"))
print("eliminatedAngleD", D["eliminatedAngleD"])
