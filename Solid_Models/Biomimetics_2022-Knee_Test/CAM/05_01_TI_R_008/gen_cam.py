"""Generate CAM programs for 05_01_TI_R_008 (tibia plate) from extracted geometry.

Program frame (all three files):
  G54 X0 Y0 = center of finished part footprint in XY (part X center = 18.915 mm,
             part Y center = 33.95 mm in STL frame).
  Z0 = finished TOP face of the part as fixtured for that op.
Imperial output (G20), 4 decimals. Fanuc/HAAS-compatible RS-274 dialect.

Part geometry (mm, STL frame; x 0..37.83, y 0..67.90, z 0..36.425):
  top plate : z 24.000 .. 36.425, full 37.83 x 67.90
  pad A     : z 13.475 .. 24.000, y 0 .. 30.00
  pad B     : z 14.615 .. 24.000, y 48.65 .. 67.90
  base rail : z  0.000 .. 13.475, y 0 .. 5.00
  holes thru:
    D7.00 @ (18.915, 17.500)  depth 22.950 (plate+pad A)
    D6.60 @ (18.915, 57.400)  depth 21.815 (plate+pad B)
    D6.60 @ (28.214, 39.325)  depth 12.425 (plate only)
    D6.60 @ (10.414, 39.325)  depth 12.425 (plate only)
"""
import math, json, os

OUT = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\CAM\05_01_TI_R_008"

MM = 25.4
def IN(mm):
    return round(mm / MM, 4)

# ---------- geometry (mm, part frame centered) ----------
PW, PYLEN, PZ = 37.83, 67.90, 36.425          # width, length, height
XC, YC = PW / 2.0, PYLEN / 2.0               # 18.915, 33.95
Z_PLATE_BOT = 24.000
Z_PADA_BOT, PADB_Y0 = 13.475, 48.65
Z_PADB_BOT, RAIL_Y1, PADA_Y1 = 14.615, 5.00, 30.00
HOLES = [  # (dia, x_mm, y_mm, depth_mm, tag)
    (7.00, 18.915, 17.500, 22.950, "D7.00"),
    (6.60, 18.915, 57.400, 21.815, "D6.60-B"),
    (6.60, 28.214, 39.325, 12.425, "D6.60-M2"),
    (6.60, 10.414, 39.325, 12.425, "D6.60-M1"),
]
def hx(x_mm): return IN(x_mm - XC)
def hy(y_mm): return IN(y_mm - YC)

# ---------- stock & process ----------
STOCK = (1.500, 2.750, 1.500)                 # in: X, Y, Z (6061-T6 bar)
FACE_ALLOW = 0.060                            # in, total facing off top
SAFE_Z = 1.0                                  # in above Z0
CLEAR_Z = 0.5

# tools
T = {
    "EM12":  {"n": "1/2 3FL carbide flat EM",   "d_in": 0.500, "rpm": 5200, "feed": 45.0, "plunge": 12.0},
    "SPOT":  {"n": "1/2 90deg spot drill",     "d_in": 0.500, "rpm": 3400, "feed": 5.0,  "plunge": 5.0},
    "DRG":   {"n": "Letter G drill (0.261in / 6.60mm)", "d_in": 0.261, "rpm": 4200, "feed": 8.0, "plunge": 8.0},
    "DRJ":   {"n": "Letter J drill (0.277in / 7.04mm)", "d_in": 0.277, "rpm": 4000, "feed": 8.0, "plunge": 8.0},
    "CHAM":  {"n": "1/2 90deg chamfer mill",   "d_in": 0.500, "rpm": 5200, "feed": 20.0, "plunge": 10.0},
}

def header(f, title, tools_used, notes):
    f.write("(%)\n")
    f.write(f"(PROGRAM: {title})\n")
    f.write("(PART: 05_01_TI_R_008  tibia plate, right)\n")
    f.write("(MATERIAL: assumed 6061-T6 aluminum — verify before cutting)\n")
    f.write("(UNITS: INCH. G54 X0Y0 = part footprint center, Z0 = part top face as fixtured)\n")
    f.write("(STOCK: 1.500 x 2.750 x 1.500 in 6061-T6 bar)\n")
    for n in notes:
        f.write(f"({n})\n")
    f.write("(TOOLS:)\n")
    for k in tools_used:
        t = T[k]
        f.write(f"(  T{t['slot']} {t['n']}  S{t['rpm']} F{t['feed']})\n")
    f.write("%\n")

def abs_g90(f):
    f.write("G90 G94 G17 G20 G40 G49 G80\n")

def profile_rect(f, dx, dy, extra, z, feed):
    """Contour rectangle (dx,dy = part dims) with tool center offset (tool radius + extra) outside walls."""
    x0 = -(dx / 2 + 0.250 + extra)
    x1 = (dx / 2 + 0.250 + extra)
    y0 = -(dy / 2 + 0.250 + extra)
    y1 = (dy / 2 + 0.250 + extra)
    f.write(f"G00 Z{SAFE_Z:.4f}\n")
    f.write(f"G00 X{x0:.4f} Y{y0:.4f}\n")
    f.write(f"G01 Z{z:.4f} F{T['EM12']['plunge']:.1f}\n")
    f.write(f"G01 X{x0:.4f} Y{y1:.4f} F{feed:.1f}\n")
    f.write(f"G01 X{x1:.4f} Y{y1:.4f}\n")
    f.write(f"G01 X{x1:.4f} Y{y0:.4f}\n")
    f.write(f"G01 X{x0:.4f} Y{y0:.4f}\n")
    f.write(f"G00 Z{SAFE_Z:.4f}\n")

# ================= OP10: TOP SIDE (3-axis) =================
T["EM12"]["slot"], T["SPOT"]["slot"] = 1, 2
T["DRG"]["slot"], T["DRJ"]["slot"], T["CHAM"]["slot"] = 3, 4, 5

p10 = os.path.join(OUT, "05_01_TI_R_008_OP10_TOP_3AXIS.nc")
with open(p10, "w") as f:
    header(f, "OP10 TOP — face, holes, profile (3-AXIS)", ["EM12", "SPOT", "DRG", "DRJ", "CHAM"],
           ["FIXTURE: stock glued/bolted to SACRIFICIAL PLATE (profile runs 0.020 into it).",
            "No parallels under the footprint. TOP face up, long edge (2.750) along Y.",
            "Stock top sits 0.066 above finished Z0 (two 0.033 facing passes).",
            "T1 faces top, spot+drill all 4 holes, chamfer, then full-depth profile.",
            "HOLES (in, G54 center): D7.00 @ (" + f"{hx(18.915):+.4f},{hy(17.5):+.4f})" +
            " | D6.60 @ (" + f"{hx(18.915):+.4f},{hy(57.4):+.4f}), (" +
            f"{hx(28.214):+.4f},{hy(39.325):+.4f}), ({hx(10.414):+.4f},{hy(39.325):+.4f})",
            "PECK: deep holes >3xD use G83."])
    abs_g90(f)
    # --- T1 face + profile ---
    f.write("(--- T1 1/2 EM: FACE TOP ---)\n")
    f.write("T1 M06\n")
    f.write(f"S{T['EM12']['rpm']} M03\n")
    f.write("G54\n")
    f.write("M08\n")
    # facing: X-direction lanes across the full stock width; part Y half = 1.3364 must be covered
    # total stock over-height: 1.500 - 1.4340 = 0.066 -> passes at -0.033 and -0.066
    ypass = [-1.09, -0.72, -0.35, 0.02, 0.39, 0.76, 1.13]
    zpass = [-0.033, -0.066]
    for z in zpass:
        for yy in ypass:
            f.write(f"G00 Z{SAFE_Z:.4f}\n")
            f.write(f"G00 X{-STOCK[0]/2 - 0.15:.4f} Y{yy:.4f}\n")
            f.write(f"G01 Z{z:.4f} F{T['EM12']['plunge']:.1f}\n")
            f.write(f"G01 X{STOCK[0]/2 + 0.15:.4f} F{T['EM12']['feed']:.1f}\n")
        f.write(f"G00 Z{SAFE_Z:.4f}\n")
    # profile outer rectangle (part 37.83 x 67.90 + tool comp left in coordinates: cut ON outline)
    # rough steps 0.150 then finish pass full depth; tabs at Y sides mid-height
    part_dx, part_dy = IN(PW), IN(PYLEN)
    prof_off_r = IN(0.5 * 1.0 + 0.010)   # rough: tool center 0.260 from wall (r=0.250 + 0.010 stock)
    prof_off_f = IN(0.5)                 # finish: tool center exactly r from wall
    depth_total = IN(PZ) + 0.020         # 0.020 into sacrificial/parrel zone
    tab_z = -IN(PZ) + 0.060              # tab top
    f.write("(--- T1: PROFILE ROUGH (0.010 wall stock, 0.150 steps) ---)\n")
    z = -0.150
    while z >= -depth_total + 1e-9:
        zuse = max(z, -depth_total)
        profile_rect(f, part_dx, part_dy, 0.010, zuse, T["EM12"]["feed"])
        z -= 0.150
    f.write("(--- T1: PROFILE FINISH full depth (tool center r outside final wall) ---)\n")
    tx0, tx1 = -part_dx / 2 - 0.250, part_dx / 2 + 0.250
    ty0, ty1 = -part_dy / 2 - 0.250, part_dy / 2 + 0.250
    f.write(f"G00 Z{SAFE_Z:.4f}\n")
    f.write(f"G00 X{tx0:.4f} Y{ty0 - 0.15:.4f}\n")
    f.write(f"G01 Z{-depth_total:.4f} F{T['EM12']['plunge']:.1f}\n")
    f.write(f"G01 Y{ty0:.4f} F{T['EM12']['feed']:.1f}\n")
    f.write(f"G01 X{tx1:.4f}\n")
    f.write(f"G01 Y{ty1:.4f}\n")
    f.write(f"G01 X{tx0:.4f}\n")
    f.write(f"G01 Y{ty0:.4f}\n")
    f.write(f"G00 Z{SAFE_Z:.4f}\n")
    f.write("(Part held by adhesive on sacrificial plate — no tabs needed; see README)\n")
    f.write("M09\nM05\n")

    # --- T2 spot ---
    f.write("(--- T2 SPOT all 4 holes 0.060 deep ---)\n")
    f.write("T2 M06\n")
    f.write(f"S{T['SPOT']['rpm']} M03\nG54\nM08\n")
    for d, x, y, dep, tag in HOLES:
        f.write(f"G00 X{hx(x):.4f} Y{hy(y):.4f}\n")
        f.write(f"G81 Z-0.060 R0.100 F{T['SPOT']['feed']:.1f} (SPOT {tag})\n")
        f.write("G80\n")
    f.write("M09\nM05\n")

    # --- T3 letter G: three 6.60 holes ---
    f.write("(--- T3 DRILL 3x D6.60 (letter G 0.261 / metric 6.60) ---)\n")
    f.write("T3 M06\n")
    f.write(f"S{T['DRG']['rpm']} M03\nG54\nM08\n")
    for d, x, y, dep, tag in HOLES:
        if d != 6.60:
            continue
        zdrill = -(IN(dep) + 0.040)
        f.write(f"G00 X{hx(x):.4f} Y{hy(y):.4f}\n")
        f.write(f"G83 Z{zdrill:.4f} R0.100 Q0.150 F{T['DRG']['feed']:.1f} ({tag} depth {dep}mm)\n")
        f.write("G80\n")
    f.write("M09\nM05\n")

    # --- T4 letter J: one 7.00 hole ---
    f.write("(--- T4 DRILL D7.00 (letter J 0.277 / metric 7.00) ---)\n")
    f.write("T4 M06\n")
    f.write(f"S{T['DRJ']['rpm']} M03\nG54\nM08\n")
    d, x, y, dep, tag = HOLES[0]
    f.write(f"G00 X{hx(x):.4f} Y{hy(y):.4f}\n")
    f.write(f"G83 Z{-(IN(dep)+0.040):.4f} R0.100 Q0.150 F{T['DRJ']['feed']:.1f} ({tag})\n")
    f.write("G80\n")
    f.write("M09\nM05\n")

    # --- T5 chamfer (optional) ---
    f.write("(--- T5 CHAMFER 0.010-0.015 at hole mouths only ---)\n")
    f.write("(Outer top edge: chamfer manually or with a tool-specific post — omitted here)\n")
    f.write("T5 M06\n")
    f.write(f"S{T['CHAM']['rpm']} M03\nG54\nM08\n")
    for d, x, y, dep, tag in HOLES:
        f.write(f"G00 X{hx(x):.4f} Y{hy(y):.4f}\n")
        f.write(f"G81 Z-0.035 R0.100 F{T['CHAM']['feed']:.1f} (CHAM {tag})\n")
        f.write("G80\n")
    f.write(f"G00 Z{SAFE_Z:.4f}\n")
    f.write(f"G00 X{tx0:.4f} Y{ty0:.4f}\n")
    f.write("G01 Z-0.020 F10.0\n")
    f.write(f"G01 X{tx1:.4f} Y{ty0:.4f} F{T['CHAM']['feed']:.1f}\n")
    f.write(f"G01 X{tx1:.4f} Y{ty1:.4f}\n")
    f.write(f"G01 X{tx0:.4f} Y{ty1:.4f}\n")
    f.write(f"G01 X{tx0:.4f} Y{ty0:.4f}\n")
    f.write(f"G00 Z{SAFE_Z:.4f}\n")
    f.write("M09\nM05\nG00 Z2.0\nM30\n%\n")
print("wrote", p10)

# ================= OP20: BOTTOM SIDE (3-axis) =================
p20 = os.path.join(OUT, "05_01_TI_R_008_OP20_BOTTOM_3AXIS.nc")
with open(p20, "w") as f:
    header(f, "OP20 BOTTOM — 3 stepped pockets (3-AXIS)", ["EM12"],
           ["FLIP part: bottom (stepped) face UP. Z0 = surface that touches parallels",
            "= former part z=0 plane. Indicate X-Y off the two machined long side walls.",
            "POCKETS (part mm frame): P1 y5..30 deep 13.475 | P2 y30..48.65 deep 24.000",
            "                    | P3 y48.65..67.9 deep 14.615. All span full width.",
            "Leave 0.010 floor/wall stock on rough; finish floors + step walls.",
            "OP10 profiled the outer rectangle already — do NOT cut outside x +/-0.7415."])
    abs_g90(f)
    f.write("T1 M06\n")
    f.write(f"S{T['EM12']['rpm']} M03\nG54\nM08\n")
    r = 0.25                     # tool radius in
    pockets = [  # (y0_mm, y1_mm, depth_mm, name)
        (RAIL_Y1, PADA_Y1, Z_PADA_BOT, "P1 beside rail"),
        (PADA_Y1, PADB_Y0, Z_PLATE_BOT, "P2 middle gap"),
        (PADB_Y0, PYLEN, Z_PADB_BOT, "P3 beyond pad B"),
    ]
    for y0m, y1m, depm, name in pockets:
        yA, yB = hy(y0m) + r + 0.010, hy(y1m) - r - 0.010      # rough wall line + 0.010 stock
        xA, xB = -IN(PW) / 2 + r + 0.010, IN(PW) / 2 - r - 0.010
        depth = IN(depm)
        f.write(f"(--- POCKET {name}: y {y0m}..{y1m} mm, depth {depm} mm ---)\n")
        z = -0.200
        # zig-zag lanes (stepover 0.200) per level
        nlane = max(1, int(math.ceil((yB - yA) / 0.200)))
        while z >= -depth - 1e-9:
            zuse = max(z, -depth)
            f.write(f"G00 Z{SAFE_Z:.4f}\n")
            f.write(f"G00 X{xA:.4f} Y{yA:.4f}\n")
            f.write(f"G01 Z{zuse:.4f} F{T['EM12']['plunge']:.1f}\n")
            for i in range(nlane + 1):
                yy = yA + (yB - yA) * i / nlane
                if i % 2 == 0:
                    f.write(f"G01 X{xB:.4f} Y{yy:.4f} F{T['EM12']['feed']:.1f}\n")
                else:
                    f.write(f"G01 X{xA:.4f} Y{yy:.4f} F{T['EM12']['feed']:.1f}\n")
            f.write(f"G00 Z{SAFE_Z:.4f}\n")
            z -= 0.200
        # finish pass: floor cut + walls one spring pass at final depth, full wall height
        yAf, yBf = hy(y0m) + r, hy(y1m) - r
        xAf, xBf = -IN(PW) / 2 + r, IN(PW) / 2 - r
        f.write(f"(finish {name})\n")
        f.write(f"G00 Z{SAFE_Z:.4f}\n")
        f.write(f"G00 X{xAf:.4f} Y{yAf:.4f}\n")
        f.write(f"G01 Z{-depth:.4f} F{T['EM12']['plunge']:.1f}\n")
        f.write(f"G01 X{xBf:.4f} Y{yAf:.4f} F{T['EM12']['feed']*0.7:.1f}\n")
        f.write(f"G01 X{xBf:.4f} Y{yBf:.4f}\n")
        f.write(f"G01 X{xAf:.4f} Y{yBf:.4f}\n")
        f.write(f"G01 X{xAf:.4f} Y{yAf:.4f}\n")
        f.write(f"G00 Z{SAFE_Z:.4f}\n")
    f.write("M09\nM05\nG00 Z2.0\nM30\n%\n")
print("wrote", p20)

# ================= 2-AXIS (XY CNC, manual Z) =================
p2a = os.path.join(OUT, "05_01_TI_R_008_2AXIS_XY.nc")
with open(p2a, "w") as f:
    header(f, "2-AXIS XY + MANUAL Z — holes + contours with depth schedule", [],
           ["MACHINE: 2-axis (XY) CNC, manual quill Z, imperial. Spindle per README table.",
            "Z IS MANUAL: run each contour pass, then lower Z by the listed step at the",
            "[SET Z] checkpoints. Program commands Z only as comments.",
            "Work from the same 1.500 x 2.750 x 1.500 stock; do OP10-equivalent first",
            "(drill from top), then flip and run pocket contours.",
            "HOLES: set quill depth stop per table; positions below are G54 XY (part center).",
            "TOOL for contours: 1/2 in 3FL carbide EM; hole drills: letter G (3x) + J (1x)."])
    f.write("G90 G94 G17 G20\n")
    f.write("(===== STAGE 1 — DRILLING, part TOP up =====)\n")
    f.write("(Set spot depth 0.060 | G-hole depth stop 0.940 (12.425mm plate holes) or 0.905 for the two deep ones +0.040)\n")
    f.write("(DRILL DEPTH TABLE (from top face): D7.00 hole: 0.940 in | D6.60@B: 0.900 in | D6.60 mid: 0.530 in)\n")
    for d, x, y, dep, tag in HOLES:
        f.write(f"(MOVE {tag}: X{hx(x):+.4f} Y{hy(y):+.4f}  quill depth {IN(dep)+0.040:.4f})\n")
        f.write(f"G00 X{hx(x):.4f} Y{hy(y):.4f}\n")
        f.write("(  [MANUAL] spot, then peck-drill to depth stop)\n")
        f.write("G04 P1.0\n")
    f.write("(===== STAGE 2 — OUTER PROFILE, top still up =====)\n")
    f.write("(Contour rectangle: tool center at part outline + 0.250 (1/2 EM radius)\n")
    f.write("(X from -0.9947 to +0.9947, Y from -1.3236 to +1.3236)\n")
    tx0, tx1 = -IN(PW) / 2 - 0.25, IN(PW) / 2 + 0.25
    ty0, ty1 = -IN(PYLEN) / 2 - 0.25, IN(PYLEN) / 2 + 0.25
    prof_depth = IN(PZ) + 0.020
    f.write(f"(DEPTH SCHEDULE: 0.100 per pass, final -{prof_depth:.3f} (0.020 into sacrificial)\n")
    npass_prof = int(math.ceil(prof_depth / 0.100))
    for p in range(npass_prof):
        zd = min((p + 1) * 0.100, prof_depth)
        f.write(f"(--- pass {p+1}: [SET Z] = -{zd:.3f} ---)\n")
        f.write(f"G00 X{tx0:.4f} Y{ty0:.4f}\n")
        f.write(f"G01 X{tx0:.4f} Y{ty1:.4f} F45.0\n")
        f.write(f"G01 X{tx1:.4f} Y{ty1:.4f}\n")
        f.write(f"G01 X{tx1:.4f} Y{ty0:.4f}\n")
        f.write(f"G01 X{tx0:.4f} Y{ty0:.4f}\n")
        f.write("([MANUAL] raise clear, lower Z 0.100, rerun)\n" if p < 14 else "(done — profile complete)\n")
    f.write("(===== STAGE 3 — POCKETS, flip bottom up =====)\n")
    f.write("(Same re-zero as 3-axis OP20. Three rectangles, tool center inset 0.250+0.010 rough:\n")
    pockets = [
        (RAIL_Y1, PADA_Y1, Z_PADA_BOT, "P1 beside rail"),
        (PADA_Y1, PADB_Y0, Z_PLATE_BOT, "P2 middle gap"),
        (PADB_Y0, PYLEN, Z_PADB_BOT, "P3 beyond pad B"),
    ]
    for y0m, y1m, depm, name in pockets:
        yA, yB = hy(y0m) + 0.26, hy(y1m) - 0.26
        xA, xB = -IN(PW) / 2 + 0.26, IN(PW) / 2 - 0.26
        nlan = max(1, int(math.ceil((yB - yA) / 0.200)))
        f.write(f"(--- {name}: Y {hy(y0m):+.4f}..{hy(y1m):+.4f}, total depth {IN(depm):.4f} ---)\n")
        ztot = IN(depm)
        npass = int(math.ceil(ztot / 0.100))
        for p in range(npass):
            f.write(f"(pass {p+1} [SET Z] = -{min((p+1)*0.100, ztot):.4f})\n")
            for i in range(nlan + 1):
                yy = yA + (yB - yA) * i / nlan
                if i % 2 == 0:
                    f.write(f"G00 X{xA:.4f} Y{yy:.4f}\n")
                    f.write(f"G01 X{xB:.4f} Y{yy:.4f} F45.0\n")
                else:
                    f.write(f"G00 X{xB:.4f} Y{yy:.4f}\n")
                    f.write(f"G01 X{xA:.4f} Y{yy:.4f} F45.0\n")
        f.write(f"(finish walls: contour X{xA:.4f}..{xB:.4f} Y{hy(y0m)+0.25:.4f}..{hy(y1m)-0.25:.4f} at final Z)\n")
    f.write("M30\n%\n")
print("wrote", p2a)

# ---------- machine-readable geometry summary for the folder ----------
geo = {
    "part_frame_mm": {"x": [0, PW], "y": [0, PYLEN], "z": [0, PZ]},
    "top_plate": {"z": [Z_PLATE_BOT, PZ], "xy": [0, 0, PW, PYLEN]},
    "pad_A": {"z": [Z_PADA_BOT, Z_PLATE_BOT], "y": [0, PADA_Y1]},
    "pad_B": {"z": [Z_PADB_BOT, Z_PLATE_BOT], "y": [PADB_Y0, PYLEN]},
    "base_rail": {"z": [0, Z_PADA_BOT], "y": [0, RAIL_Y1]},
    "holes_mm": [{"d": d, "x": x, "y": y, "depth": dep, "tag": t} for d, x, y, dep, t in HOLES],
    "drawing_crosscheck": {"22.95": "plate top to pad A bottom", "21.81": "plate top to pad B bottom",
                            "18.65": "middle gap width (y30..48.65)", "19.25": "pad B depth",
                            "10.50": "pad A step (24.0-13.475)", "17.80": "mid hole spacing"},
    "source": "STL slice extraction 2026-10-16; verified against drawing 05_01_TI_R_008.SLDDRW dims",
}
with open(os.path.join(OUT, "part_geometry.json"), "w") as fjson:
    json.dump(geo, fjson, indent=1)
print("wrote part_geometry.json")
