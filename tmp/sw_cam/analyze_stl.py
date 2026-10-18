"""Extract exact geometry from 05_01_TI_R_008.STL for CAM programming.

Strategy: slice the mesh at many Z levels; at each level collect closed loops of
the cross-section; classify loops (outer vs holes); fit circles to small loops;
report rectilinear step profile per Z slab.
Units: mm.
"""
import struct, math, json
from collections import defaultdict

import numpy as np

SRC = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\CAM\05_01_TI_R_008\05_01_TI_R_008.stl"

with open(SRC, "rb") as f:
    header = f.read(80)
    n = struct.unpack("<I", f.read(4))[0]
    tris = np.frombuffer(f.read(n * 50), dtype=np.uint8).reshape(n, 50)

verts = np.zeros((n, 3, 3), dtype=np.float64)
for i in range(n):
    v = struct.unpack("<9f", tris[i, 12:48])   # skip normal
    verts[i] = np.array(v).reshape(3, 3)

print(f"tris: {n}, bbox: {verts.reshape(-1,3).min(0)} .. {verts.reshape(-1,3).max(0)}")

def slice_loops(z):
    """Return list of closed polylines (list of (x,y)) at plane z."""
    segs = []
    for t in verts:
        d = t[:, 2] - z
        if (d > 0).all() or (d < 0).all():
            continue
        pts = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            da, db = d[a], d[b]
            if (da > 0) == (db > 0):
                continue
            if da == db:
                continue
            s = da / (da - db)
            pts.append(t[a, :2] + s * (t[b, :2] - t[a, :2]))
        if len(pts) == 2:
            segs.append((tuple(pts[0]), tuple(pts[1])))
        elif len(pts) > 2 and len(pts) % 2 == 0:
            for k in range(0, len(pts), 2):
                segs.append((tuple(pts[k]), tuple(pts[k + 1])))
    if not segs:
        return []
    # join segments into loops
    adj = defaultdict(list)
    for i, (p, q) in enumerate(segs):
        adj[p].append((i, q))
        adj[q].append((i, p))
    used = [False] * len(segs)
    loops = []
    for i in range(len(segs)):
        if used[i]:
            continue
        start = segs[i][0]
        loop = [start]
        cur = segs[i][1]
        used[i] = True
        while True:
            loop.append(cur)
            nxt = None
            for (j, other) in adj[cur]:
                if not used[j]:
                    nxt = (j, other)
                    break
            if nxt is None:
                break
            used[nxt[0]] = True
            cur = nxt[1]
            if cur == start:
                break
        loops.append(loop)
    return loops

def loop_area(loop):
    a = 0.0
    for i in range(len(loop) - 1):
        x1, y1 = loop[i]
        x2, y2 = loop[i + 1]
        a += x1 * y2 - x2 * y1
    return a / 2.0

def fit_circle(pts):
    P = np.array(pts)
    A = np.c_[2 * P[:, 0], 2 * P[:, 1], np.ones(len(P))]
    b = (P ** 2).sum(1)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy, c = sol
    r = math.sqrt(max(c + cx * cx + cy * cy, 0))
    resid = np.abs(np.hypot(P[:, 0] - cx, P[:, 1] - cy) - r).max()
    return cx, cy, r, resid

def roundness(loop):
    """1.0 = perfect circle."""
    cx, cy, r, resid = fit_circle(loop)
    if r <= 0:
        return 0, (cx, cy, r, resid)
    return max(0.0, 1 - resid / r), (cx, cy, r, resid)

# --- find flat Z levels ---
zvals = verts[:, :, 2].reshape(-1)
levels = []
for z in np.unique(np.round(zvals, 3)):
    frac = np.mean(np.abs(zvals - z) < 1e-6)
    # count triangles fully in a horizontal plane at z
    cnt = np.sum(np.all(np.abs(verts[:, :, 2] - z) < 1e-6, axis=1))
    if cnt > 20:
        levels.append((float(z), int(cnt)))
print("flat levels (z, tri count):")
for z, c in levels:
    print(f"  {z:9.3f}  {c}")

# --- slice at informative levels ---
zs = sorted(set([z for z, _ in levels] + [-16.0, 0.0, 5.0, 10.0, 15.0, 18.5]))
out = {}
for z in zs:
    loops = slice_loops(z)
    if not loops:
        continue
    info = []
    for L in loops:
        area = abs(loop_area(L))
        if area < 1e-3:
            continue
        ro, circ = roundness(L)
        cx, cy, r, resid = circ
        entry = {"area": round(area, 2), "n": len(L)}
        if ro > 0.97 and area < 100:
            entry["circle"] = {"cx": round(cx, 3), "cy": round(cy, 3),
                               "r": round(r, 3), "d": round(2 * r, 3)}
        else:
            xs = [p[0] for p in L]; ys = [p[1] for p in L]
            entry["bbox"] = [round(min(xs), 3), round(min(ys), 3),
                             round(max(xs), 3), round(max(ys), 3)]
            entry["verts"] = sorted(set((round(x, 3), round(y, 3)) for x, y in L))
        info.append(entry)
    out[z] = info
    print(f"\n=== z = {z} : {len(loops)} loops ===")
    for e in sorted(info, key=lambda e: -e["area"]):
        if "circle" in e:
            print(f"  HOLE d={e['circle']['d']:.3f} @ ({e['circle']['cx']:.3f}, {e['circle']['cy']:.3f})")
        else:
            print(f"  LOOP area={e['area']} bbox={e['bbox']}")
            print(f"        verts={e['verts']}")

with open(r"D:\GitHub\Bipedal_Robot\tmp\sw_cam\stl_geometry.json", "w") as f:
    json.dump({str(k): v for k, v in out.items()}, f, indent=1)
print("\nsaved stl_geometry.json")
