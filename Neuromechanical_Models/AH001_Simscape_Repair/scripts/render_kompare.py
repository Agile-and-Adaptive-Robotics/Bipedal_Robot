#!/usr/bin/env python3
"""render_kompare.py -- robot CAD + dissertation 27-route set (KompareXXVII,
side=selected) in the standing ISB world. Route points: route_points.csv
(frame-local), placed via the ah001_neural.xml body chain. Survey sub-frames
(calcnA/calcnSub/femurB/tibia_rigid/pelvisShift/torsoA/torsoB) approximated by
their parent frames (disclosed)."""
import csv
import json
import os
import re
import struct
from collections import defaultdict

import colorsys
import numpy as np
import plotly.graph_objects as go

REPO = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot"
K27 = os.path.join(REPO, r"Code\Matlab\Mesh_Optimization\Results\KompareXXVII")
RT = os.path.join(REPO, r"Code\MuJoCo_SNS\robot_torso")
REP = os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair")
STL_DIR = os.path.join(REP, "cad", "stl")

# ---------- body chain from the XML ----------
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

# survey frame -> MuJoCo body (left-referenced robotbody; sub-frames by parent)
FRAME2BODY = {
    "pelvis": "pelvis", "pelvisShift": "pelvis",
    "femur": "thigh_l", "femurB": "thigh_l",
    "tibia": "shank_l", "tibia_rigid": "shank_l",
    "calcn": "foot_l", "calcnA": "foot_l", "calcnSub": "foot_l", "toes": "foot_l",
    "torso": "lumbar_torso", "torsoA": "lumbar_torso", "torsoB": "lumbar_torso",
}

# ---------- selected configs ----------
cfg = {}
for row in csv.DictReader(open(os.path.join(K27, "route_configurations.csv"))):
    if row["side"] == "selected":
        cfg[row["site"]] = row

# ---------- route points (selected) ----------
by_site = defaultdict(list)
for row in csv.DictReader(open(os.path.join(K27, "route_points.csv"))):
    if row["side"] != "selected":
        continue
    by_site[row["site"]].append(row)

palette = []
for i in range(27):
    r, g, b = colorsys.hsv_to_rgb(i / 27.0, 0.75, 0.95)
    palette.append(f"rgb({int(r*255)},{int(g*255)},{int(b*255)})")

fig = go.Figure()
npts_total = 0
for si, (site, rows) in enumerate(sorted(by_site.items(), key=lambda kv: int(kv[1][0]["actuator_id"]))):
    rows.sort(key=lambda r: int(r["point_index"]))
    body_seen = set()
    P = []
    for r in rows:
        b = FRAME2BODY.get(r["frame"], "pelvis")
        body_seen.add(r["frame"])
        Rb, tb = world(b)
        P.append(Rb @ np.array([float(r["x_m"]), float(r["y_m"]), float(r["z_m"])]) + tb)
    P = np.array(P)
    if len(P) < 2:
        continue
    c = cfg.get(site, {})
    diam = c.get("diameter_saved_mm", "?"); nb = c.get("nBPA", "?"); G = c.get("G_survey", "?")
    rest = c.get("resting_length_m", ""); tend = c.get("tendon_length_m", "")
    hov = (f"{site} ({c.get('muscle','')})<br>{nb}x{diam}mm, G={G}<br>"
           f"rest {float(rest)*1000:.0f} mm, tendon {float(tend)*1000:.0f} mm"
           if rest and tend else site)
    fig.add_trace(go.Scatter3d(
        x=P[:,0], y=P[:,2], z=P[:,1], mode="lines+markers",
        line=dict(color=palette[si % 27], width=7),
        marker=dict(size=[6] + [3]*(len(P)-2) + [6]),
        name=f"{site} {nb}x{diam} G{G}", legendgroup=site,
        hovertemplate=hov + "<extra></extra>"))
    npts_total += len(P)
print("routes drawn:", len(by_site), "points:", npts_total)

# ---------- CAD STLs (if present) ----------
fit = json.load(open(os.path.join(REP, "cad", "cad_to_mujoco_fit.json")))
Rfit = np.array(fit["R"]); tfit = np.array(fit["t"])
leg_xml = open(os.path.join(REP, "cad", "export", "ah001_full.xml.xml"), encoding="utf-8", errors="replace").read()
LEG = [("02_01_PE_001", "02_01_PE_001-1"),
       ("03_01_FB_R_001", "03_01_FB_R_001-1"), ("03_02_FH_R_001", "03_02_FH_R_001-1"),
       ("03_03_FB_L_001", "03_03_FB_L_001-1"), ("03_04_FH_L_001", "03_04_FH_L_001-1"),
       ("04_01_KT_R_001", "04_01_KT_R_001-1"), ("04_02_KB_R_001", "04_02_KB_R_001-1"),
       ("04_03_KT_L_001", "04_03_KT_L_001-1"), ("04_04_KB_L_001", "04_04_KB_L_001-1"),
       ("04_05_BL_001", "04_05_BL_001-1"), ("04_05_BL_001", "04_05_BL_001-2"),
       ("04_06_FL_001", "04_06_FL_001-1"), ("04_06_FL_001", "04_06_FL_001-2"),
       ("05_01_TI_R_001", "05_01_TI_R_001-1"), ("05_02_TI_001", "05_02_TI_001-1")]

def load_stl(path):
    with open(path, "rb") as f:
        head = f.read(84)
        n = struct.unpack("<I", head[80:84])[0]
        if 84 + 50 * n == os.path.getsize(path):
            data = np.fromfile(path, dtype=np.uint8, offset=84, count=n*50).reshape(n, 50)
            t = np.zeros((n, 3, 3))
            for k in range(3):
                t[:, k, :] = data[:, 12+12*k:24+12*k].copy().view("<f4").reshape(n, 3)
            return t
    txt = open(path, "r", errors="replace").read()
    v = np.array(re.findall(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", txt), dtype=float).reshape(-1, 3, 3)
    return v

nstl = 0
for stem, inst in LEG:
    f = os.path.join(STL_DIR, stem + ".STL")
    if not os.path.exists(f) or os.path.getsize(f) < 200:
        continue
    tris = load_stl(f)
    if np.abs(tris).max() > 20: tris = tris / 1000.0
    if len(tris) > 12000: tris = tris[:: len(tris)//12000 + 1]
    m = re.search(r'<Instance name="%s"[^>]*>.*?<Rotation>([^<]+)</Rotation>\s*<Translation>([^<]+)</Translation>' % re.escape(inst), leg_xml, re.S)
    if not m: continue
    Ra = np.array([float(v) for v in m.group(1).split()]).reshape(3,3)
    ta = np.array([float(v) for v in m.group(2).split()])
    Rw = Rfit @ Ra; tw = Rfit @ ta + tfit
    P = tris.reshape(-1, 3) @ Rw.T + tw
    I, J, K = np.arange(len(P)).reshape(-1, 3).T
    fig.add_trace(go.Mesh3d(x=P[:,0], y=P[:,2], z=P[:,1], i=I, j=J, k=K,
                            opacity=0.28, color="#93a2b4", hoverinfo="skip",
                            name=stem, showlegend=False))
    nstl += 1
print("CAD parts rendered:", nstl)


# ---------- feet: foot-assembly STLs placed rel. Ankle Exterior, on the MuJoCo foot bodies ----------
fx = open(os.path.join(REP, "cad", "export_foot", "foot_ankle.xml.xml"), encoding="utf-8", errors="replace").read()
def foot_pose(name):
    m = re.search(r'<Instance name="%s"[^>]*>.*?<Rotation>([^<]+)</Rotation>\s*<Translation>([^<]+)</Translation>' % re.escape(name), fx, re.S)
    if not m: return None
    return (np.array([float(v) for v in m.group(1).split()]).reshape(3,3),
            np.array([float(v) for v in m.group(2).split()]))
FOOT_INST = [("Ankle Exterior", "Ankle Exterior-1"),
             ("Ankle Middle Joint Shell", "Ankle Middle Joint Shell-1"),
             ("Ankle Interior Joint", "Ankle Interior Joint-1"),
             ("Ankle Bottom Joint", "Ankle Bottom Joint-1"),
             ("New foot - Copy - Copy", "New foot - Copy - Copy-1"),
             ("Left mid foot", "Left mid foot-1"),
             ("New Front foot", "New Front foot-1"),
             ("Toe", "Toe-1"), ("Toe", "Toe-2"), ("Toe", "Toe-3"), ("Toe", "Toe-4"), ("Toe", "Toe-5")]
Re_, te_ = foot_pose("Ankle Exterior-1")
for body_name, mirror in (("foot_l", 1), ("foot_r", -1)):
    Rb, tb = world(body_name)
    for stem, inst in FOOT_INST:
        f = os.path.join(STL_DIR, stem + ".STL")
        pr = foot_pose(inst)
        if not os.path.exists(f) or pr is None:
            continue
        Rp, tp = pr
        Rrel = Re_.T @ Rp; trel = Re_.T @ (tp - te_)
        tris = load_stl(f)
        if np.abs(tris).max() > 20: tris = tris / 1000.0
        if len(tris) > 8000: tris = tris[:: len(tris)//8000 + 1]
        Rw = Rb @ Rrel; tw = Rb @ trel + tb
        if mirror < 0:
            P = tris.reshape(-1,3) @ Rw.T + tw
            P[:, 0] = -P[:, 0]   # mirror across the sagittal plane for the right foot (visual)
            I, J, K = np.arange(len(P)).reshape(-1,3).T
        else:
            P = tris.reshape(-1,3) @ Rw.T + tw
            I, J, K = np.arange(len(P)).reshape(-1,3).T
        fig.add_trace(go.Mesh3d(x=P[:,0], y=P[:,2], z=P[:,1], i=I, j=J, k=K,
                                opacity=0.25, color="#aebccb", hoverinfo="skip",
                                name=stem+"_"+body_name, showlegend=False))

fig.update_layout(
    scene=dict(aspectmode="data",
               xaxis_title="x ant (m)", yaxis_title="z right (m)", zaxis_title="y up (m)"),
    title="AH001 CAD + KompareXXVII selected 27-route BPA set (dissertation Table 4.12 designs)",
    height=1050)
out_html = os.path.join(REP, "robot_bpa27_kompare.html")
fig.write_html(out_html, include_plotlyjs="cdn")
print("HTML:", out_html)
try:
    fig.write_image(os.path.join(REP, "robot_bpa27_kompare.png"), width=1700, height=1300, scale=2)
    print("PNG written")
except Exception as e:
    print("PNG skipped:", str(e)[:80])
