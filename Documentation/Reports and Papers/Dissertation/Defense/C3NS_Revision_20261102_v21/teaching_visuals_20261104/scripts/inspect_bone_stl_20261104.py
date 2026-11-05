"""Inspect OpenSim bone STLs (from mjc Geometry) vs the registered point clouds."""
import numpy as np
import struct
import json

ROOT = r"C:/Users/Ben/Documents/GitHub/Bipedal_Robot"


def read_stl(path):
    with open(path, "rb") as f:
        head = f.read(84)
        n = struct.unpack("<I", head[80:84])[0]
        data = f.read()
    if len(data) >= n * 50:
        arr = np.frombuffer(data[: n * 50], dtype=np.uint8).reshape(n, 50)
        tri = arr[:, 12:48].copy().view("<f4").reshape(n, 3, 3)
        return tri.astype(np.float64), n, "binary"
    # ascii fallback
    with open(path, "r", errors="ignore") as f:
        verts = []
        for line in f:
            p = line.split()
            if p and p[0] == "vertex":
                verts.append([float(v) for v in p[1:4]])
        return np.array(verts).reshape(-1, 3, 3), len(verts) // 3, "ascii"


for name in ("femur_r", "tibia_r"):
    p = rf"{ROOT}/Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/Geometry/{name}.stl"
    tri, n, fmt = read_stl(p)
    v = tri.reshape(-1, 3)
    print(f"{name}.stl: {fmt}, {n} tris, {len(v)} verts")
    print("  bounds x [%.4f %.4f] y [%.4f %.4f] z [%.4f %.4f]"
          % (v[:, 0].min(), v[:, 0].max(), v[:, 1].min(), v[:, 1].max(), v[:, 2].min(), v[:, 2].max()))

D = json.load(open(rf"{ROOT}/tmp/cam_visuals_dump_20261104.json"))
be = D["bone_extents"]
print("\ncloud femur bounds:", np.round(np.array(be["femur"]), 4).tolist())
print("cloud tibia bounds:", np.round(np.array(be["tibia"]), 4).tolist())
