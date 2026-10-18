"""STL geometry extraction v2: snapped-point loop chaining, clean circle fitting."""
import struct, math, json
from collections import defaultdict
import numpy as np

SRC = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\CAM\05_01_TI_R_008\05_01_TI_R_008.stl"

with open(SRC, "rb") as f:
    f.read(80)
    n = struct.unpack("<I", f.read(4))[0]
    tris = np.frombuffer(f.read(n * 50), dtype=np.uint8).reshape(n, 50)

verts = np.zeros((n, 3, 3))
for i in range(n):
    verts[i] = np.array(struct.unpack("<9f", tris[i, 12:48])).reshape(3, 3)

def snap(p, nd=4):
    return (round(p[0], nd), round(p[1], nd))

def slice_loops(z):
    segs = []
    for t in verts:
        d = t[:, 2] - z
        if (d > 1e-9).all() or (d < -1e-9).all():
            continue
        pts = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            da, db = d[a], d[b]
            if (da > 0) == (db > 0) or da == db:
                continue
            s = da / (da - db)
            pts.append(t[a, :2] + s * (t[b, :2] - t[a, :2]))
        for k in range(0, len(pts) - 1, 2):
            p, q = snap(pts[k]), snap(pts[k + 1])
            if p != q:
                segs.append((p, q))
    if not segs:
        return []
    adj = defaultdict(set)
    for i, (p, q) in enumerate(segs):
        adj[p].add(i)
        adj[q].add(i)
    used = [False] * len(segs)
    loops = []
    for i0 in range(len(segs)):
        if used[i0]:
            continue
        p0, cur = segs[i0]
        used[i0] = True
        loop = [p0]
        guard = 0
        while cur != p0 and guard < 100000:
            loop.append(cur)
            guard += 1
            nxt = None
            for j in adj[cur]:
                if used[j]:
                    continue
                a, b = segs[j]
                nxt = b if a == cur else a
                used[j] = True
                break
            if nxt is None:
                break
            cur = nxt
        if cur == p0 and len(loop) >= 8:
            loops.append(loop)
    return loops

def loop_area(L):
    a = 0.0
    for i in range(len(L)):
        x1, y1 = L[i]; x2, y2 = L[(i + 1) % len(L)]
        a += x1 * y2 - x2 * y1
    return a / 2

def fit_circle(pts):
    P = np.array(pts)
    A = np.c_[2 * P[:, 0], 2 * P[:, 1], np.ones(len(P))]
    b = (P ** 2).sum(1)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy, c = sol
    r = math.sqrt(max(c + cx * cx + cy * cy, 0))
    resid = np.abs(np.hypot(P[:, 0] - cx, P[:, 1] - cy) - r).max()
    return cx, cy, r, resid

out = {}
probe_z = [0.05, 2, 4, 6, 8, 12, 16, 20, 24, 28, 32, 34, 36.3]
for z in probe_z:
    loops = slice_loops(z)
    if not loops:
        out[z] = []
        continue
    entries = []
    for L in loops:
        area = abs(loop_area(L))
        cx, cy, r, resid = fit_circle(L)
        is_circ = (r < 25) and (resid < 0.02 * max(r, 1e-6)) and len(set(L)) > 8
        if is_circ and area < 120:
            entries.append({"type": "hole", "d": round(2 * r, 3),
                            "cx": round(cx, 3), "cy": round(cy, 3), "area": round(area, 2)})
        else:
            xs = [p[0] for p in L]; ys = [p[1] for p in L]
            entries.append({"type": "outline", "area": round(area, 1),
                            "bbox": [round(min(xs), 3), round(min(ys), 3),
                                     round(max(xs), 3), round(max(ys), 3)],
                            "verts": sorted(set(L))})
    out[z] = entries
    print(f"\n=== z={z} ===")
    for e in sorted(entries, key=lambda e: -e["area"]):
        if e["type"] == "hole":
            print(f"  HOLE d={e['d']:.3f} @ ({e['cx']:.3f}, {e['cy']:.3f})")
        else:
            print(f"  OUTLINE area={e['area']} bbox={e['bbox']}")
            print(f"     verts={e['verts']}")

with open(r"D:\GitHub\Bipedal_Robot\tmp\sw_cam\stl_geometry2.json", "w") as f:
    json.dump(out, f, indent=1)
print("\nsaved")
