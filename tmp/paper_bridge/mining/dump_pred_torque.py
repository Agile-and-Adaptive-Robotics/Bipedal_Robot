"""Extract predBest torque spans (READ-ONLY)."""
import scipy.io as sio
import numpy as np

MESH = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"


def span(f, lab, fields):
    m = sio.loadmat(MESH + f, squeeze_me=True, struct_as_record=False)
    pb = m["predBest"]
    print(f"===== {lab} =====")
    for fld in fields:
        try:
            a = np.atleast_1d(np.asarray(getattr(pb, fld), dtype=float)).ravel()
            a = a[np.isfinite(a)]
            if a.size:
                print(f"  predBest.{fld}: n={a.size} min={a.min():.4f} max={a.max():.4f}")
        except Exception as e:
            print(f"  predBest.{fld}: ({e})")


span(r"\Bifemsh_20mm_Result_20260921_2037.mat", "FLX DIRECT 0921",
     ["Torque", "TorqueZ", "offAxisTorque", "momentArm"])
span(r"\Bifemsh_20mm_Result_pulley_20260926_1652.mat", "FLX PULLEY 0926",
     ["Torque", "TorqueZ", "TorqueInsZ", "offAxisTorque", "FinsMag", "ReactionMag"])
span(r"\Vas_Pam_20mm_Result_20260925.mat", "EXT DIRECT 0925",
     ["Torque", "TorqueZ", "offAxisTorque", "momentArm"])
