"""Verify generated G-code: bounds, depths, hole positions vs part geometry."""
import re, json, math, os, sys

CAM = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\CAM\05_01_TI_R_008"
geo = json.load(open(os.path.join(CAM, "part_geometry.json")))
MM = 25.4
PW, PYLEN, PZ = 37.83, 67.90, 36.425
XC, YC = PW / 2, PYLEN / 2
HOLES = [(h["x"], h["y"]) for h in geo["holes_mm"]]
depths = {
    "P1": 13.475 / MM, "P2": 24.000 / MM, "P3": 14.615 / MM,
}

def parse(path):
    lines = open(path).read().splitlines()
    events = []
    pos = {"X": 0.0, "Y": 0.0, "Z": 2.0}
    motion = 0
    for ln, raw in enumerate(lines, 1):
        line = raw.split("(")[0].strip()
        if not line or line.startswith("%"):
            continue
        m = re.match(r"[NT]\d+", line)
        if line.startswith(("T", "S", "M", "G04", "G80")) and not line.startswith(("G00", "G01", "G81", "G83")):
            continue
        if line.startswith("G81") or line.startswith("G83"):
            d = dict(re.findall(r"([XYZQRF])(-?\d+\.?\d*)", line))
            if "X" in d or "Y" in d:
                pos.update({k: float(v) for k, v in d.items() if k in "XY"})
            events.append({"ln": ln, "type": "cycle", "X": pos["X"], "Y": pos["Y"],
                           "Z": float(d.get("Z", 0)), "line": raw.strip()})
            continue
        if line.startswith(("G00", "G01")):
            motion = 0 if line.startswith("G00") else 1
            d = dict(re.findall(r"([XYZF])(-?\d+\.?\d*)", line))
            new = dict(pos)
            new.update({k: float(v) for k, v in d.items() if k in "XYZ"})
            events.append({"ln": ln, "type": "cut" if (motion == 1 and new["Z"] < -0.001) else "move",
                           **new})
            pos = new
            continue
        if line.startswith("G54") or line.startswith("G90") or line.startswith("M"):
            continue
        d = dict(re.findall(r"([XYZF])(-?\d+\.?\d*)", line))
        if d:
            new = dict(pos)
            new.update({k: float(v) for k, v in d.items() if k in "XYZ"})
            events.append({"ln": ln, "type": "cut" if (motion == 1 and new["Z"] < -0.001) else "move", **new})
            pos = new
    return events

def check(path, zmin_allowed, xy_margin, expect_holes):
    name = os.path.basename(path)
    ev = parse(path)
    errors, warns = [], []
    zmin = min(e["Z"] for e in ev if "Z" in e)
    if zmin < zmin_allowed - 0.005:
        errors.append(f"Z below allowed floor: {zmin:.4f} < {zmin_allowed:.4f}")
    # XY envelope applies to CUTTING points only (rapids at safe Z may stage outside stock)
    cuts = [e for e in ev if e["type"] == "cut"]
    if cuts:
        xmax = max(abs(e["X"]) for e in cuts)
        ymax = max(abs(e["Y"]) for e in cuts)
    else:
        xmax = ymax = 0.0
    if xmax > PW / 2 / MM + xy_margin or ymax > PYLEN / 2 / MM + xy_margin:
        errors.append(f"XY exceeds stock+tool margin: X {xmax:.4f} Y {ymax:.4f} "
                      f"(allowed {PW/2/MM+xy_margin:.4f} / {PYLEN/2/MM+xy_margin:.4f})")
    # hole cycle positions
    got = sorted(set((round(e["X"], 3), round(e["Y"], 3)) for e in ev if e["type"] == "cycle"))
    want = sorted(set((round((x - XC) / MM, 3), round((y - YC) / MM, 3)) for x, y in expect_holes))
    for g in got:
        if g not in want:
            warns.append(f"cycle at non-hole position {g}")
    missing = [w for w in want if w not in got]
    if missing:
        errors.append(f"missing hole positions in cycles: {missing}")
    ncuts = sum(1 for e in ev if e["type"] == "cut")
    print(f"\n== {name}")
    print(f"   events {len(ev)} | cutting moves {ncuts} | Zmin {zmin:.4f} (allowed {zmin_allowed:.4f})")
    print(f"   |X|max {xmax:.4f}  |Y|max {ymax:.4f} | cycle positions: {len(got)}")
    for e in errors:
        print("   ERROR:", e)
    for w in warns:
        print("   WARN :", w)
    return not errors

ok = True
ok &= check(os.path.join(CAM, "05_01_TI_R_008_OP10_TOP_3AXIS.nc"),
            -(PZ / MM) - 0.020 - 0.005, 0.45, HOLES)
ok &= check(os.path.join(CAM, "05_01_TI_R_008_OP20_BOTTOM_3AXIS.nc"),
            -(24.000 / MM) - 0.005, 0.30, [])
ok &= check(os.path.join(CAM, "05_01_TI_R_008_2AXIS_XY.nc"), 0.0, 0.45, [])

# 2-axis file: hole positions are G00 targets (manual drilling)
txt2 = open(os.path.join(CAM, "05_01_TI_R_008_2AXIS_XY.nc")).read()
g00xy = set()
for m in re.finditer(r"G00 X(-?\d+\.\d+) Y(-?\d+\.\d+)", txt2):
    g00xy.add((round(float(m.group(1)), 3), round(float(m.group(2)), 3)))
want = set((round((x - XC) / MM, 3), round((y - YC) / MM, 3)) for x, y in HOLES)
hole_ok = want <= g00xy
print("\n2-axis G00 hole positions present:", hole_ok, "" if hole_ok else f"missing {want - g00xy}")
ok &= hole_ok

# 2-axis: per-section [SET Z] schedules must get monotonically deeper within a section
txt = txt2
sections = re.split(r"===== STAGE|\(--- P", txt)
mono = True
detail = []
for si, sec in enumerate(sections):
    zs = [float(v) for v in re.findall(r"\[SET Z\] = (-?\d+\.\d+)", sec)]
    if not zs:
        continue
    ok_sec = all(zs[i+1] <= zs[i] + 1e-9 for i in range(len(zs) - 1)) and zs[0] > zs[-1]
    detail.append((si, zs[0], zs[-1], len(zs), ok_sec))
    mono &= ok_sec
for d in detail:
    print(f"   section {d[0]}: {d[1]} -> {d[2]} in {d[3]} steps, per-section monotonic: {d[4]}")
print(f"2-axis [SET Z] sections checked: {len(detail)}; all per-section monotonic: {mono}")
print("\nOVERALL:", "PASS" if ok and mono else "FAIL")
