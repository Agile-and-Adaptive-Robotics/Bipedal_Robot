"""Bisect Z transitions using the clean slicer from analyze_stl2."""
import struct, math, json
from collections import defaultdict
import numpy as np

SRC = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\CAM\05_01_TI_R_008\05_01_TI_R_008.stl"
with open(SRC, "rb") as f:
    f.read(80); n = struct.unpack("<I", f.read(4))[0]
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
        adj[p].add(i); adj[q].add(i)
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
            loop.append(cur); guard += 1
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

def state(z):
    loops = slice_loops(z)
    outer = sorted(abs(loop_area(L)) for L in loops if abs(loop_area(L)) > 120)[-1] if loops else 0
    return round(sum(abs(loop_area(L)) for L in loops if abs(loop_area(L)) > 120), 1)

def bisect(zlo, zhi, tol=0.002):
    slo, shi = state(zlo), state(zhi)
    print(f"  state({zlo})={slo}  state({zhi})={shi}")
    while zhi - zlo > tol:
        zm = (zlo + zhi) / 2
        sm = state(zm)
        if abs(sm - slo) < abs(sm - shi):
            zhi = zm
        else:
            zlo = zm
    return (zlo + zhi) / 2

print("bisecting rail+pads -> ? transitions")
for lo, hi in [(12.0, 16.0), (20.0, 24.0)]:
    zc = bisect(lo, hi)
    print(f"  transition near z = {zc:.4f}")
    print(f"    just below: {state(zc-0.01)}, just above: {state(zc+0.01)}")

# also check lower transitions 0..12 for hidden steps in the rail
for z in [0.5, 3, 6, 9, 11.5, 12.5, 14.5, 15.5]:
    print(f"  z={z:6.2f} area={state(z)}")
