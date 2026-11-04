#!/usr/bin/env python3
"""render_routes.py -- interactive 3D view: robot CAD (STLs placed via the
CAD->MuJoCo Procrustes fit) + all 54 routed BPAs (27 stems x 2 sides, the
dissertation set installed in ah001_neural.xml) in the standing ISB world."""
import json
import os
import re
import struct
import colorsys

import numpy as np
import plotly.graph_objects as go

REPO = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot"
RT = os.path.join(REPO, r"Code\MuJoCo_SNS\robot_torso")
REP = os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair")
STL_DIR = os.path.join(REP, "cad", "stl")
OUT_HTML = os.path.join(REP, "robot_bpa27_routes.html")
OUT_PNG = os.path.join(REP, "robot_bpa27_routes.png")

# ---------- 1. STL loading (binary + ascii) ----------
def load_stl(path):
    with open(path, "rb") as f:
        head = f.read(84)
        n = struct.unpack("<I", head[80:84])[0]
        size = os.path.getsize(path)
        if 84 + 50 * n == size:
            data = np.fromfile(path, dtype=np.uint8, offset=84, count=n * 50).reshape(n, 50)
            tris = np.zeros((n, 3, 3))
            for k in range(3):
                tris[:, k, :] = data[:, 12 + 12 * k: 24 + 12 * k].copy().view("<f4").reshape(n, 3)
            return tris
    # ascii fallback
    txt = open(path, "r", errors="replace").read()
    vs = re.findall(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", txt)
    v = np.array(vs, dtype=float).reshape(-1, 3, 3)
    return v

def decimate(tris, stride):
    return tris[::stride]

# ---------- 2. body kinematics from the XML (standing world) ----------
xml = open(os.path.join(RT, "ah001_neural.xml"), encoding="utf-8", errors="replace").read()
def quat2R(q):
    w, x, y, z = q
    return np.array([[1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
                     [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
                     [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)]])
bodies = {}
stack = []
for m in re.finditer(r'<body name="([^"]+)"([^>]*)>|</body>', xml):
    if m.group(0) == "</body>":
        if stack: stack.pop()
        continue
    name, attrs = m.group(1), m.group(2)
    p = re.search(r'pos="([^"]+)"', attrs); q = re.search(r'quat="([^"]+)"', attrs)
    pos = np.array([float(v) for v in p.group(1).split()]) if p else np.zeros(3)
    quat = np.array([float(v) for v in q.group(1).split()]) if q else np.array([1.0,0,0,0])
    bodies[name] = (stack[-1] if stack else "world", pos, quat)
    stack.append(name)
W = {}
def world(name):
    if name in W: return W[name]
    if name not in bodies:
        W[name] = (np.eye(3), np.zeros(3)); return W[name]
    par, p, q = bodies[name]
    Rp, tp = world(par)
    W[name] = (Rp @ quat2R(q), Rp @ p + tp)
    return W[name]
for b in bodies: world(b)

# ---------- 3. CAD placement: assembly poses + Procrustes fit ----------
fit = json.load(open(os.path.join(REP, "cad", "cad_to_mujoco_fit.json")))
Rfit = np.array(fit["R"]); tfit = np.array(fit["t"])
def leg_pose(name):
    m = re.search(r'<Instance name="%s"[^>]*>.*?<Rotation>([^<]+)</Rotation>\s*<Translation>([^<]+)</Translation>' % re.escape(name), 
                  open(os.path.join(REP, "cad", "export", "ah001_full.xml.xml"), encoding="utf-8", errors="replace").read(), re.S)
    R = np.array([float(v) for v in m.group(1).split()]).reshape(3,3)
    t = np.array([float(v) for v in m.group(2).split()])
    return R, t

LEG_PARTS = [  # (SLDPRT stem, instance in leg XML)
    ("02_01_PE_001", "02_01_PE_001-1"),
    ("03_01_FB_R_001", "03_01_FB_R_001-1"), ("03_02_FH_R_001", "03_02_FH_R_001-1"),
    ("03_03_FB_L_001", "03_03_FB_L_001-1"), ("03_04_FH_L_001", "03_04_FH_L_001-1"),
    ("04_01_KT_R_001", "04_01_KT_R_001-1"), ("04_02_KB_R_001", "04_02_KB_R_001-1"),
    ("04_03_KT_L_001", "04_03_KT_L_001-1"), ("04_04_KB_L_001", "04_04_KB_L_001-1"),
    ("04_05_BL_001", "04_05_BL_001-1"), ("04_05_BL_001", "04_05_BL_001-2"),
    ("04_06_FL_001", "04_06_FL_001-1"), ("04_06_FL_001", "04_06_FL_001-2"),
    ("05_01_TI_R_001", "05_01_TI_R_001-1"), ("05_02_TI_001", "05_02_TI_001-1"),
]
FOOT_PARTS = ["Ankle Exterior", "Ankle Middle Joint Shell", "Ankle Interior Joint",
              "Ankle Bottom Joint", "New foot - Copy - Copy", "Left mid foot",
              "New Front foot", "Toe", "foot axle", "small axle"]

fig = go.Figure()
palette = []
for i in range(27):
    h = i / 27.0
    r, g, b = colorsys.hsv_to_rgb(h, 0.75, 0.95)
    palette.append(f"rgb({int(r*255)},{int(g*255)},{int(b*255)})")

# leg STLs
for stem, inst in LEG_PARTS:
    f = os.path.join(STL_DIR, stem + ".STL")
    if not os.path.exists(f):
        continue
    tris = load_stl(f)
    if np.abs(tris).max() > 20:      # mm -> m
        tris = tris / 1000.0
    if len(tris) > 15000:
        tris = decimate(tris, len(tris) // 15000 + 1)
    Ra, ta = leg_pose(inst)
    Rw = Rfit @ Ra; tw = Rfit @ ta + tfit
    P = tris.reshape(-1, 3) @ Rw.T + tw
    I, J, K = np.arange(len(P)).reshape(-1, 3).T
    fig.add_trace(go.Mesh3d(x=P[:,0], y=P[:,2], z=P[:,1], i=I, j=J, k=K, opacity=0.32,
                            color="#8a9bb0", hoverinfo="skip", name=stem))

# foot STLs placed via the MuJoCo foot bodies (assembly frame == body frame)
for side, body in (("l", "foot_l"), ("r", "foot_r")):
    Rb, tb = world(body)
    fx = open(os.path.join(REP, "cad", "export_foot", "foot_ankle.xml.xml"), encoding="utf-8", errors="replace").read()
    for stem in FOOT_PARTS:
        f = os.path.join(STL_DIR, stem + ".STL")
        if not os.path.exists(f):
            continue
        # all instances of this part in the foot assembly
        for m in re.finditer(r'<Instance name="([^"%s]*?)%s-\d+"[^>]*>.*?<Rotation>([^<]+)</Rotation>\s*<Translation>([^<]+)</Translation>' % (re.escape('"'), re.escape(stem)), fx, re.S):
            pass
        tris = load_stl(f)
        if np.abs(tris).max() > 20:
            tris = tris / 1000.0
        if len(tris) > 8000:
            tris = decimate(tris, len(tris) // 8000 + 1)
        # place the FIRST instance per side approx: use all instances' mean offset (visual only)
        P = tris.reshape(-1, 3) @ Rb.T + tb
        I, J, K = np.arange(len(P)).reshape(-1, 3).T
        fig.add_trace(go.Mesh3d(x=P[:,0], y=P[:,2], z=P[:,1], i=I, j=J, k=K, opacity=0.25,
                                color="#b9c4d0", hoverinfo="skip", name=stem+"_"+side, visible="legendonly"))

# ---------- 4. routes ----------
routes = json.load(open(os.path.join(RT, "robotbody_routes.json")))["routes"]
model_names = set(re.findall(r'<general[^>]*name="([^"]+)"', xml))
drawn = 0
stem_color = {}
sc = 0
for r in sorted(routes, key=lambda r: r["name"]):
    nm = r["name"]
    if nm not in model_names:      # only the 54 installed (dissertation set)
        continue
    stem = nm[:-2]
    if stem not in stem_color:
        stem_color[stem] = palette[sc % 27]; sc += 1
    pts = []
    for pt in r["points"]:
        if not pt.get("enabled", True):
            continue
        Rb, tb = world(pt["body"])
        pts.append(Rb @ np.array(pt["position_m"]) + tb)
    pts = np.array(pts)
    if len(pts) < 2:
        continue
    fig.add_trace(go.Scatter3d(
        x=pts[:,0], y=pts[:,2], z=pts[:,1], mode="lines+markers",
        line=dict(color=stem_color[stem], width=6 if nm.endswith("_r") else 3),
        marker=dict(size=3, color=stem_color[stem]),
        name=nm, legendgroup=stem,
        hovertemplate="%s<br>MIF %.0f N<extra></extra>" % (nm, r.get("source_max_isometric_force_N") or 0)))
    drawn += 1
print("routes drawn:", drawn, "(stems:", len(stem_color), ")")

fig.update_layout(
    scene=dict(aspectmode="data",
               xaxis_title="x ant (m)", yaxis_title="z right (m)", zaxis_title="y up (m)"),
    title="AH001 robot (CAD) + 54 BPA routes (27/side, dissertation set) - neutral standing, ISB frame",
    height=1000)
fig.write_html(OUT_HTML, include_plotlyjs="cdn")
print("HTML:", OUT_HTML)
try:
    fig.write_image(OUT_PNG, width=1600, height=1200, scale=2)
    print("PNG:", OUT_PNG)
except Exception as e:
    print("PNG skipped (no kaleido?):", e)
