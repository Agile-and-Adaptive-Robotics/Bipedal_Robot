"""Verify eq:app_Tk prose claims against the archived per-frame Tk (kf(1)):
rotation block R_z(theta_k) sweeping -120..+10 deg, constant translation
[4.5, -395.8, 0] mm, bottom row [0 0 0 1]."""
import scipy.io as sio
import numpy as np

d = sio.loadmat(r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo\FlxPinBPASet.mat")
Tk = d["kf"].flat[0]["Tk"]  # 4x4x100
Ak = d["kf"].flat[0]["Ak"]  # knee angle per frame presumably
print("Tk shape:", Tk.shape, " Ak shape:", Ak.shape)

thetas = np.degrees(np.arctan2(Tk[1, 0, :], Tk[0, 0, :]))  # from R_z block
print("theta from R_z block: first=%.3f last=%.3f min=%.3f max=%.3f deg"
      % (thetas[0], thetas[-1], thetas.min(), thetas.max()))

# rotation block checks: R = [[c,-s,0],[s,c,0],[0,0,1]]
R = Tk[:3, :3, :]
for i in (0, 50, 99):
    Ri = R[:, :, i]
    orth = np.allclose(Ri @ Ri.T, np.eye(3), atol=1e-12)
    det = np.linalg.det(Ri)
    print(f"frame {i}: orth={orth} det={det:.12f} R={np.round(Ri,4).tolist()}")

# translation column
p = Tk[:3, 3, :] * 1000.0  # mm
print("translation col: mean =", np.round(p.mean(axis=1), 4),
      " ptp =", np.ptp(p, axis=1))

# bottom row
print("bottom row unique:", np.unique(Tk[3, :, :]))

# Ak vs theta
print("Ak first/last/min/max (deg): %.3f %.3f %.3f %.3f"
      % (Ak.flat[0], Ak.flat[-1], Ak.min(), Ak.max()))

# cos/sin consistency with theta for frame 0 and last
for i in (0, 99):
    c, s = np.cos(np.radians(thetas[i])), np.sin(np.radians(thetas[i]))
    print(f"frame {i}: theta={thetas[i]:.4f}  Tk00={Tk[0,0,i]:.6f} (cos={c:.6f})  "
          f"Tk10={Tk[1,0,i]:.6f} (sin={s:.6f})  Tk01={Tk[0,1,i]:.6f} (-sin={-s:.6f})")
