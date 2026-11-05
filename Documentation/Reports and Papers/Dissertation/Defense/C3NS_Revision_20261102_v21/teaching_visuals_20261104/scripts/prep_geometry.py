"""prep_geometry.py — derive Blender-ready geometry from the validated sources.

Inputs (all of record, read-only):
  data/cam_visuals_dump_20261104.json   (dumped from Vas_Pam_20mm_Result_cam_20261014_0416.mat)
  Geometry STLs from the MyoConverter conversion of gait2392_simbody.osim
  Gait2392 knee_r SimmSpline tables (gait2392_simbody.osim)

Outputs:
  data/visual1_geom.npz   cam-extensor per-pose geometry in the femur frame
  data/visual3_geom.npz   Gait2392 knee translation tables on a dense q grid
"""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parents[1]
D = json.load(open(HERE / "data" / "cam_visuals_dump_20261104.json"))

route = np.array(D["route_raw"])            # 9 x 3 x N  (rows 1-5 femur frame, 6-9 t1 frame)
active = np.array(D["route_active"])        # 9 x N bool
T_Pam = np.array(D["T_Pam"])                # 4 x 4 x N  (ICR -> femur)
T_ICR_t1 = np.array(D["T_ICR_t1"])          # 4 x 4 x N  (t1 -> ICR, pure translation)
phiD = np.array(D["phiD"])                  # deg, -120 .. +10
N = D["n_poses"]
geo = D["geo"]

# shaped cam drum of record (Ben, 2026-10-14 wrap-radius rule):
#   center = circumcenter of the contacted full-flexion contacts p6/p7/p8 (t1 frame)
#   radius = |p7x - p9x| at full extension = 96.9 mm (the p7x-p9x rule;
#            wrapRadius row 2 in the dumped mat = 96.9 mm constant)
seeds = np.array([geo["seedP6XY"], geo["seedP7XY"], geo["seedP8XY"]])  # 3x2 [m]
(x1, y1), (x2, y2), (x3, y3) = seeds
a = np.array([[x2 - x1, y2 - y1], [x3 - x1, y3 - y1]]) * 2.0
b = np.array([x2**2 - x1**2 + y2**2 - y1**2, x3**2 - x1**2 + y3**2 - y1**2])
cx, cy = np.linalg.solve(a, b)
r_loss = abs(float(route[6, 0, N - 1] - route[8, 0, N - 1]))   # 96.9 mm loss-model radius of record
r_drum = float(np.hypot(x1 - cx, y1 - cy))                     # clip-part surface arc through the contacts
assert abs(r_loss - 0.0969) < 2e-3, f"loss radius {r_loss:.4f} m != 96.9 mm of record"
print(f"cam clip arc: center t1 ({cx*1000:.2f}, {cy*1000:.2f}) mm, radius {r_drum*1000:.2f} mm; "
      f"loss-model drum radius of record {r_loss*1000:.2f} mm (p7x-p9x rule)")

# per-pose quantities in the femur frame
def rvt(T, p):
    return np.asarray(p, float) @ T[:3, :3].T + T[:3, 3]

route_f = np.zeros_like(route)                # (9, 3, N) femur frame
t1toICR = np.zeros((N, 3))                    # t1 origin offset from the ICR
for i in range(N):
    t1toICR[i] = -T_ICR_t1[:3, 3, i]
    for j in range(9):
        if j < 5:
            route_f[j, :, i] = route[j, :, i]
        else:
            route_f[j, :, i] = rvt(T_Pam[:, :, i], rvt(T_ICR_t1[:, :, i], route[j, :, i]))

pivot = T_Pam[:3, 3, :].T.copy()             # (N,3) migrating knee pivot, femur frame
tibia_rot = np.deg2rad(phiD).copy()          # tibia rotation about z, per pose

# moment arm of record: perpendicular from pivot to line of action through p6
# along the incoming femoral segment (validated: mean |diff| 0.18 mm vs evaluator)
mA = np.array(D["momentArm_m"])
loa_p = route_f[5, :, :].T.copy()            # (N,3) point on line = p6 (femur frame)
loa_d = np.zeros((N, 3))
for i in range(N):
    u = route_f[4, :, i] - route_f[5, :, i]  # pull direction p6 -> p5
    loa_d[i] = u / np.linalg.norm(u)
foot = np.zeros((N, 3))
for i in range(N):
    w = pivot[i] - loa_p[i]
    foot[i] = loa_p[i] + np.dot(w, loa_d[i]) * loa_d[i]

np.savez_compressed(
    HERE / "data" / "visual1_geom.npz",
    phiD=phiD, route_f=route_f, active=active, pivot=pivot, tibia_rot=tibia_rot,
    t1toICR=t1toICR,
    drum_center_t1=np.array([cx, cy, 0.0]), drum_radius=r_drum, loss_radius=r_loss,
    seeds_t1=np.column_stack([seeds, np.zeros(3)]),
    p_end_t1=route[8, :, 0].copy(),
    femur_pulley_center=np.array([geo["femurCylCenter"][0], geo["femurCylCenter"][1], 0.0]),
    femur_pulley_radius=float(geo["femurCylRadius"]),
    moment_arm=mA, loa_p=loa_p, loa_d=loa_d, foot=foot,
    femur_reg_off=np.array(D["bone_registration"]["femur_offset_used"] + [0.0]),
    tibia_reg_off=np.array(D["bone_registration"]["tibia_offset_used"] + [0.0]),
    p1_femur=route_f[:, 0, :].copy(),
)
print(f"visual1_geom.npz: {N} poses, mA {mA.min()*1000:.1f}..{mA.max()*1000:.1f} mm")

# ---- TibiaClip_002 hardware registration (actual CAD part) ---------------
# Part-frame wrap center from fixed-radius Hough on the real mesh (59.83 mm
# radius, matching the evaluator circumcircle of the adopted contacts).
# Rotation: the evaluator documents p8 as the "clip-edge" contact, so the
# registered arc ENDPOINT lands on p8's angular position; p6/p7 land strictly
# inside the surface (verified in scripts/registration checks, tmp history).
# z: the two halves straddle z 0..92.2 mm with a 0.5 mm cable slot at the
# mid-plane -> scene z = part z - 46.1.
C_P = np.array([119.1, 52.9])
DELTA_DEG = 145.95
Z_MID = 46.1
REPO = HERE.parents[5]
np.savez_compressed(
    HERE / "data" / "clip_registration.npz",
    C_P=C_P, delta_deg=DELTA_DEG, C_T1=np.array([cx, cy]), z_mid=Z_MID,
    stl_lateral=str(REPO / "Solid_Models/Biomimetics_2022-Knee_Test/STL Files/"
                    "TibiaClip_002_lateral.STL"),
    stl_medial=str(REPO / "Solid_Models/Biomimetics_2022-Knee_Test/STL Files/"
                    "TibiaClip_002_medial.STL"),
)
print(f"clip registration written (delta {DELTA_DEG} deg, z_mid {Z_MID} mm)")

# ---- visual 3: Gait2392 knee translation tables, dense grid --------------
x1 = np.array([-2.0944, -1.74533, -1.39626, -1.0472, -0.698132, -0.349066, -0.174533,
               0.197344, 0.337395, 0.490178, 1.52146, 2.0944])
y1 = np.array([-0.0032, 0.00179, 0.00411, 0.0041, 0.00212, -0.001, -0.0031,
               -0.005227, -0.005435, -0.005574, -0.005435, -0.00525])
x2 = np.array([-2.0944, -1.22173, -0.523599, -0.349066, -0.174533, 0.159149, 2.0944])
y2 = np.array([-0.4226, -0.4082, -0.399, -0.3976, -0.3966, -0.395264, -0.396])
f1 = CubicSpline(x1, y1, bc_type="natural")   # OpenSim SimmSpline = natural cubic spline
f2 = CubicSpline(x2, y2, bc_type="natural")
qg = np.linspace(-2.0944, 0.17453293, 881)
np.savez_compressed(
    HERE / "data" / "visual3_geom.npz",
    q_grid=qg, tx=f1(qg), ty=f2(qg),
)
print("visual3_geom.npz: dense knee tables written")
print("DONE_PREP")
