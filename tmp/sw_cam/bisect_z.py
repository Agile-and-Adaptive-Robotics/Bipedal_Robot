"""Bisect exact Z levels where the cross-section changes (pad bottoms, top-plate bottom)."""
import struct, math
from collections import defaultdict
import numpy as np

SRC = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\CAM\05_01_TI_R_008\05_01_TI_R_008.stl"
with open(SRC, "rb") as f:
    f.read(80); n = struct.unpack("<I", f.read(4))[0]
    tris = np.frombuffer(f.read(n * 50), dtype=np.uint8).reshape(n, 50)
verts = np.zeros((n, 3, 3))
for i in range(n):
    verts[i] = np.array(struct.unpack("<9f", tris[i, 12:48])).reshape(3, 3)

def footprint(z):
    """Return set of triangle indices above z whose... simpler: total cross-section area + count of loops by bbox."""
    segs = 0
    ymin, ymax = 1e9, -1e9
    for t in verts:
        d = t[:, 2] - z
        if (d > 1e-9).all() or (d < -1e-9).all():
            continue
        # triangle crosses: sample its xy extent contribution
        segs += 1
    return segs

def area_at(z):
    """Total planar area of cross-section via rasterization-free approach: use the
    an|0 estimator = sum over triangles of their clipped polygon areas."""
    total = 0.0
    for t in verts:
        d = t[:, 2] - z
        pos = d > 0
        if pos.all() or (~pos).all():
            continue
        # clip triangle to plane z: collect intersection segment
        pts = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            da, db = d[a], d[b]
            if (da > 0) == (db > 0):
                continue
            s = da / (da - db)
            pts.append(t[a, :2] + s * (t[b, :2] - t[a, :2]))
        if len(pts) >= 2:
            # degenerate segment contributes nothing to area; skip
            pass
    return total

# Simpler robust approach: sample solid occupancy on a 2D grid at height z using
# ray parity along Y for each grid column X — expensive but exact enough.
def solid_area(z, nx=380, ny=679):
    xs = np.linspace(0.001, 37.829, nx)
    ys = np.linspace(0.001, 67.899, ny)
    # for each column x, find inside-intervals by crossing count of triangle edges in XY projected? no —
    # use z-crossing segments: a point (x,y) is inside the solid at height z iff the vertical ray
    # from (x,y) upward crosses the surface an odd number of times ABOVE z... compute per column x:
    inside = np.zeros((nx, ny), dtype=bool)
    # gather all triangles that could cover grid: compute barycentric containment in XY of triangles
    # whose plane z-range brackets z — approximate with point-in-2D-triangle of the FULL projected mesh
    # parity: count intersections of vertical line at (x,y) with ALL triangles at heights > z
    # Precompute per-triangle 2D bbox
    tri2d = verts[:, :, :2]
    tb0 = tri2d.min(1); tb1 = tri2d.max(1)
    for ix, x in enumerate(xs):
        cand = np.where((tb0[:, 0] <= x) & (tb1[:, 0] >= x))[0]
        for iy, y in enumerate(ys):
            cnt = 0
            for ti in cand:
                t = verts[ti]
                if t[:, 2].min() <= z <= t[:, 2].max() or t[:, 2].min() > z:
                    pass
                # exact: intersect vertical line with triangle plane
                p0, p1, p2 = t
                # barycentric in projection
                det = (p1[1]-p2[1])*(p0[0]-p2[0]) + (p2[0]-p1[0])*(p0[1]-p2[1])
                if abs(det) < 1e-12:
                    continue
                l0 = ((p1[1]-p2[1])*(x-p2[0]) + (p2[0]-p1[0])*(y-p2[1]))/det
                l1 = ((p2[1]-p0[1])*(x-p2[0]) + (p0[0]-p2[0])*(y-p2[1]))/det
                l2 = 1-l0-l1
                if l0 < -1e-9 or l1 < -1e-9 or l2 < -1e-9:
                    continue
                zh = l0*p0[2] + l1*p1[2] + l2*p2[2]
                if zh > z + 1e-9:
                    cnt += 1
            inside[ix, iy] = (cnt % 2) == 1
    return inside.sum() * (37.829/nx) * (67.899/ny)

import sys
for z in [float(a) for a in sys.argv[1:]] or [13.0, 14.0, 15.0, 21.0, 22.0, 23.0]:
    a = solid_area(z, nx=76, ny=136)  # coarse grid for speed
    print(f"z={z:7.3f}  solid area ~= {a:9.1f} mm^2")
