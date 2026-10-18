"""Wall-guard check: no cutting point may violate keep-out geometry.

OP10 (top side): below facing depth (-0.070), tool center must stay >= 0.245 in
outside the part outline rectangle (tool radius 0.250, 0.005 tolerance).
OP20 (bottom side, flipped): at depth reaching part height h = -Z*25.4, tool center
must stay >= 0.245 from the solid bodies present at that level:
  rail  y in [0, 5.00]    for h <= 13.475
  padA  y in [0, 30.00]   for h <= 24.000
  padB  y in [48.65,67.90] for h <= 24.000
  plate all y             for h <= 36.425
(X span of all bodies: 0..37.83.)
Distances measured in inches; XY in G54 part-center coords (mm converted).
"""
import math, os, re, sys

CAM = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\CAM\05_01_TI_R_008"
MM = 25.4
PW, PYLEN = 37.83 / MM, 67.90 / MM           # part half-extents inches
HX, HY = PW / 2, PYLEN / 2
R = 0.250

def rect_clear(px, py, x0, x1, y0, y1):
    """Distance from point to rectangle interior; 0 if inside."""
    dx = max(x0 - px, 0, px - x1)
    dy = max(y0 - py, 0, py - y1)
    return math.hypot(dx, dy)

def parse_cutpoints(path):
    pts = []
    pos = {"X": 0.0, "Y": 0.0, "Z": 2.0}
    motion = 0
    for ln, raw in enumerate(open(path).read().splitlines(), 1):
        line = raw.split("(")[0].strip()
        if not line:
            continue
        if line.startswith(("G00", "G01")):
            motion = 0 if line.startswith("G00") else 1
        d = dict(re.findall(r"([XYZ])(-?\d+\.?\d*)", line))
        if d:
            for k in "XYZ":
                if k in d:
                    pos[k] = float(d[k])
            if motion == 1 and pos["Z"] < -0.001:
                pts.append((ln, pos["X"], pos["Y"], pos["Z"]))
    return pts

def op10_guard():
    pts = parse_cutpoints(os.path.join(CAM, "05_01_TI_R_008_OP10_TOP_3AXIS.nc"))
    bad = []
    for ln, x, y, z in pts:
        if z > -0.070:            # facing depth band — allowed over the whole top
            continue
        d = rect_clear(x, y, -HX, HX, -HY, HY)
        if d < R - 0.006:
            bad.append((ln, x, y, z, d))
    print(f"OP10 wall guard: {len(pts)} cutting points checked below facing depth")
    for b in bad[:10]:
        print(f"   VIOLATION line {b[0]}: X{b[1]:.4f} Y{b[2]:.4f} Z{b[3]:.4f} wall clearance {b[4]:.4f}")
    return not bad

def op20_guard():
    pts = parse_cutpoints(os.path.join(CAM, "05_01_TI_R_008_OP20_BOTTOM_3AXIS.nc"))
    bad = []
    # flipped coords: cutting depth z reaches part height h = -z*25.4 mm measured from part z0
    z_pada = 13.475 / MM; z_padb = 14.615 / MM; z_plate = 24.0 / MM
    eps = 0.0015  # absorb 4-dp depth rounding in the .nc
    y_rail1, y_pada1 = 5.0 / MM - PYLEN / 2, 30.0 / MM - PYLEN / 2
    y_padb0 = 48.65 / MM - PYLEN / 2
    for ln, x, y, z in pts:
        depth = -z
        if depth <= 0:
            continue
        bodies = []
        if depth > 0:
            bodies.append((-HY, y_rail1))     # rail occupies full 0..13.475 band
        if depth > z_pada + eps:
            bodies.append((-HY, y_pada1))     # pad A occupies 13.475..24
        if depth > z_padb + eps:
            bodies.append((y_padb0, HY))      # pad B occupies 14.615..24
        if depth > z_plate + eps:
            bodies.append((-HY, HY))          # plate occupies 24..36.425
        ok = True
        for (y0, y1) in bodies:
            if rect_clear(x, y, -HX, HX, y0, y1) < R - 0.006:
                ok = False
        if not ok:
            bad.append((ln, x, y, z))
    print(f"OP20 keep-out guard: {len(pts)} cutting points checked")
    for b in bad[:10]:
        print(f"   VIOLATION line {b[0]}: X{b[1]:.4f} Y{b[2]:.4f} Z{b[3]:.4f}")
    return not bad

ok = op10_guard() & op20_guard()
print("WALL GUARD:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
