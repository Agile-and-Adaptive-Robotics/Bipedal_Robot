#!/usr/bin/env python3
"""render_kompare_anim.py -- MODULAR articulated BPA view: knee-flexion sweep
with per-route diameter expansion D(L) (bpa_braid law ported to python).

Modularity contract (for re-optimized routes): point the ROUTE_CSV env/arg at
the new KompareXXVII-style route_points.csv (fields: site, side, point_index,
frame, x_m/y_m/z_m) + CONFIG_CSV (route_configurations.csv) -- nothing else
changes. Frames mapped: pelvis/pelvisShift->pelvis, femur/femurB->thigh_l,
tibia/tibia_rigid->shank_l, calcn*/toes->foot_l, torso*->lumbar_torso.
Articulation: knee hinge (axis from the XML joint) rotating the shank+foot
subtree; hip hinge rotating the whole leg. BPA: D(L)=sqrt(b^2-L^2)/(pi*N),
N solved from the Kompare diameter/contraction endpoints per size class.
Output: animated HTML with slider + per-frame PNGs."""
import csv
import json
import math
import os
import re
import struct
import sys
from collections import defaultdict

import numpy as np
import plotly.graph_objects as go

REPO = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot"
K27 = os.path.join(REPO, r"Code\Matlab\Mesh_Optimization\Results\KompareXXVII")
RT = os.path.join(REPO, r"Code\MuJoCo_SNS\robot_torso")
REP = os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair")
STL_DIR = os.path.join(REP, "cad", "stl")
ROUTE_CSV = sys.argv[1] if len(sys.argv) > 1 else os.path.join(K27, "route_points.csv")
CONFIG_CSV = sys.argv[2] if len(sys.argv) > 2 else os.path.join(K27, "route_configurations.csv")

# ---------------- model chain + joint hinges ----------------
xml = open(os.path.join(RT, "ah001_neural.xml"), encoding="utf-8", errors="replace").read()
def quat2R(q):
    w, x, y, z = q
    return np.array([[1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
                     [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
                     [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)]])
def axis_angle_R(axis, ang):
    a = np.asarray(axis, float); a = a / np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * (K @ K)

TREE = {}   # body -> (parent, pos, quat, {joint: (anchor_body_rel, axis)})
stack = []
for m in re.finditer(r'<body name="([^"]+)"([^>]*)>|</body>', xml):
    if m.group(0) == "</body>":
        if stack: stack.pop()
        continue
    name, attrs = m.group(1), m.group(2)
    p = re.search(r'pos="([^"]+)"', attrs); q = re.search(r'quat="([^"]+)"', attrs)
    pos = np.array([float(v) for v in p.group(1).split()]) if p else np.zeros(3)
    quat = np.array([float(v) for v in q.group(1).split()]) if q else np.array([1.0,0,0,0])
    TREE[name] = (stack[-1] if stack else "world", pos, quat, {})
    stack.append(name)
# joint axes: hinge joints carry axis= (z default); ball joints have none
# (skip ball joints in the sweep -- this animation only articulates hinges).
JOINTS = {}
for jm in re.finditer(r'<joint name="([^"]+)"[^>]*?axis="([^"]+)"[^>]*?(?:pos="([^"]+)")?[^>]*/?>', xml):
    jname, axis, jpos = jm.group(1), jm.group(2), jm.group(3)
    JOINTS[jname] = (np.array([float(v) for v in axis.split()]),
                     np.array([float(v) for v in jpos.split()]) if jpos else np.zeros(3))
JOINT2BODY = {}
stack2 = []
for m in re.finditer(r'<body name="([^"]+)"[^>]*>|</body>|<joint name="([^"]+)"[^>]*?type="hinge"', xml):
    if m.group(0) == "</body>":
        if stack2: stack2.pop()
    elif m.group(1) is not None:
        stack2.append(m.group(1))
    elif m.group(2) and stack2:
        JOINT2BODY[m.group(2)] = stack2[-1]
BJ = {}
for jn, b in JOINT2BODY.items():
    if jn not in JOINTS: continue
    ax, jp = JOINTS[jn]
    BJ.setdefault(b, []).append((ax, jp))

def fk(angles):
    """angles: {body_name: rad} applied to that body's FIRST hinge (axis/anchor
    in body-local frame), rotating the body and its subtree."""
    W = {}
    def rec(name):
        if name in W: return W[name]
        par, pos, quat, _ = TREE[name]
        Rp, tp = rec(par) if par in TREE else (np.eye(3), np.zeros(3))
        Rl = Rp @ quat2R(quat)
        tl = Rp @ pos + tp
        if name in angles and name in BJ:
            ax, anch = BJ[name][0]
            Rj = axis_angle_R(Rl @ ax, angles[name])
            # rotate about the anchor: shift, rotate, shift back
            c = tl + Rl @ anch
            tl = c + Rj @ (tl - c)
            Rl = Rj @ Rl
        W[name] = (Rl, tl)
        return W[name]
    for b in TREE: rec(b)
    return W

FRAME2BODY = {"pelvis": "pelvis", "pelvisShift": "pelvis", "femur": "thigh_l",
              "femurB": "thigh_l", "tibia": "shank_l", "tibia_rigid": "shank_l",
              "calcn": "foot_l", "calcnA": "foot_l", "calcnSub": "foot_l",
              "toes": "foot_l", "torso": "lumbar_torso", "torsoA": "lumbar_torso",
              "torsoB": "lumbar_torso"}

# ---------------- BPA law (bpa_braid port) ----------------
SIZE_KMAX_OD = {10: (0.17, 0.0225), 20: (0.25, 0.038), 40: (0.32, 0.075)}
def braid(L, L0, D0):
    Kmax, Dk = SIZE_KMAX_OD.get(int(D0 * 1000), (0.25, 1.9 * D0))
    Lmin = L0 * (1 - Kmax)
    Lc = min(max(L, Lmin), L0)
    piN = math.sqrt((L0**2 - Lmin**2) / max(Dk**2 - D0**2, 1e-12))
    b2 = L0**2 + (piN * D0) ** 2
    D = math.sqrt(max(b2 - Lc**2, 0)) / piN
    return D, Lc

# ---------------- routes + configs (modular inputs) ----------------
cfg = {}
for row in csv.DictReader(open(CONFIG_CSV)):
    if row["side"] == "selected":
        cfg[row["site"]] = row
by_site = defaultdict(list)
for row in csv.DictReader(open(ROUTE_CSV)):
    if row["side"] != "selected": continue
    by_site[row["site"]].append(row)
for s in by_site:
    by_site[s].sort(key=lambda r: int(r["point_index"]))

# neutral lengths L0 (neutral standing pose)
W0 = fk({})
def route_world(W, site):
    P = []
    for r in by_site[site]:
        Rb, tb = W[FRAME2BODY.get(r["frame"], "pelvis")]
        P.append(Rb @ np.array([float(r["x_m"]), float(r["y_m"]), float(r["z_m"])]) + tb)
    return np.array(P)
L0 = {}
for site in by_site:
    P = route_world(W0, site)
    L0[site] = float(np.sum(np.linalg.norm(np.diff(P, axis=0), axis=1)))

# ---------------- geometry: STLs ----------------
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
    return np.array(re.findall(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", txt), dtype=float).reshape(-1, 3, 3)

fitj = json.load(open(os.path.join(REP, "cad", "cad_to_mujoco_fit.json"))) if os.path.exists(os.path.join(REP, "cad", "cad_to_mujoco_fit.json")) else None
import json
fit = json.load(open(os.path.join(REP, "cad", "cad_to_mujoco_fit.json")))
Rfit = np.array(fit["R"]); tfit = np.array(fit["t"])
leg_xml = open(os.path.join(REP, "cad", "export", "ah001_full.xml.xml"), encoding="utf-8", errors="replace").read()
LEG = [("02_01_PE_001","02_01_PE_001-1","pelvis"),
       ("03_03_FB_L_001","03_03_FB_L_001-1","thigh_l"),("03_04_FH_L_001","03_04_FH_L_001-1","thigh_l"),
       ("04_03_KT_L_001","04_03_KT_L_001-1","thigh_l"),("04_04_KB_L_001","04_04_KB_L_001-1","thigh_l"),
       ("04_05_BL_001","04_05_BL_001-2","shank_l"),("04_06_FL_001","04_06_FL_001-2","thigh_l"),
       ("05_02_TI_001","05_02_TI_001-1","shank_l")]
# NOTE: CAD part -> MuJoCo-body association is approximate (visual); at neutral
# the Procrustes fit places parts correctly; under articulation the moving
# parts follow their mapped MuJoCo body's rotation about the same joints.

# ---------------- frames: knee sweep 0..-90 deg ----------------
angles_list = [0, -15, -30, -45, -60, -75, -90]
frames = []
frame_data = []
for adeg in angles_list:
    W = fk({"shank_l": math.radians(adeg)})
    fd = []
    for site in sorted(by_site, key=lambda s: int(by_site[s][0]["actuator_id"])):
        P = route_world(W, site)
        segs = np.linalg.norm(np.diff(P, axis=0), axis=1)
        L = float(np.sum(segs))
        c = cfg.get(site, {})
        D0 = float(c.get("diameter_saved_mm", 20)) / 1000.0
        D, Lc = braid(L, L0[site], D0)
        fd.append((site, P, D, L, L0[site], c))
    frame_data.append((adeg, W, fd))

# ---------------- figure with animation frames ----------------
import colorsys
pal = []
for i in range(27):
    r, g, b = colorsys.hsv_to_rgb(i / 27, 0.75, 0.95)
    pal.append(f"rgb({int(r*255)},{int(g*255)},{int(b*255)})")
sites_sorted = sorted(by_site, key=lambda s: int(by_site[s][0]["actuator_id"]))
color_of = {s: pal[i % 27] for i, s in enumerate(sites_sorted)}

fig = go.Figure()
first = True
for fi, (adeg, W, fd) in enumerate(frame_data):
    visible = (fi == 0)
    # CAD meshes (left leg only for animation lightness)
    for stem, inst, muj_body in LEG:
        f = os.path.join(STL_DIR, stem + ".STL")
        if not os.path.exists(f): continue
        m = re.search(r'<Instance name="%s"[^>]*>.*?<Rotation>([^<]+)</Rotation>\s*<Translation>([^<]+)</Translation>' % re.escape(inst), leg_xml, re.S)
        if not m: continue
        Ra = np.array([float(v) for v in m.group(1).split()]).reshape(3,3)
        ta = np.array([float(v) for v in m.group(2).split()])
        tris = load_stl(f)
        if np.abs(tris).max() > 20: tris = tris / 1000.0
        if len(tris) > 6000: tris = tris[:: len(tris)//6000 + 1]
        # neutral world placement then apply body rotation delta
        Rw0 = Rfit @ Ra; tw0 = Rfit @ ta + tfit
        Rb0, tb0 = W0[muj_body]; Rb1, tb1 = W[muj_body]
        # delta rotation of the body between neutral and this angle
        Rd = Rb1 @ Rb0.T
        pivot = tb1
        P = tris.reshape(-1, 3) @ Rw0.T + tw0
        P = (Rd @ (P - pivot).T).T + pivot
        I, J, K = np.arange(len(P)).reshape(-1, 3).T
        fig.add_trace(go.Mesh3d(x=P[:,0], y=P[:,2], z=P[:,1], i=I, j=J, k=K,
                                opacity=0.22, color="#93a2b4", hoverinfo="skip",
                                visible=visible, showlegend=False))
    # routes as tubes (radius = D/2)
    for site, P, D, L, L0s, c in fd:
        pts = P
        fig.add_trace(go.Scatter3d(
            x=pts[:,0], y=pts[:,2], z=pts[:,1], mode="lines+markers",
            line=dict(color=color_of[site], width=max(4, min(16, D*400))),
            marker=dict(size=[6]+[3]*(len(pts)-2)+[6], color=color_of[site]),
            name=site, legendgroup=site, visible=visible, showlegend=(fi == 0),
            hovertemplate=(f"{site}<br>L {L*1000:.0f} mm (L0 {L0s*1000:.0f})<br>"
                           f"D {D*1000:.1f} mm<extra></extra>")))
    frames.append(go.Frame(
        name=f"{adeg}",
        traces=list(range(len(fig.data) - len(fd) - len(LEG), len(fig.data))),
        layout=dict(title=f"Knee {adeg} deg — BPA diameters expand as routes shorten")))

def frame_args(idx):
    return [go.Frame(name=f"{a}") for a in angles_list], dict(
        buttons=[dict(label="&#9654;", method="animate",
                      args=[None, dict(frame=dict(duration=700, redraw=True), fromcurrent=True)]),
                 dict(label="&#10074;&#10074;", method="animate",
                      args=[[None], dict(mode="immediate", frame=dict(duration=0))])],
        slider=dict(steps=[dict(method="animate", label=f"{a}°",
                                args=[[f"{a}"], dict(mode="immediate,",
                                     frame=dict(duration=0, redraw=True))]) for a in angles_list]))
fig.frames = [f for f in frames]
fig.update_layout(
    scene=dict(aspectmode="data", xaxis_title="x ant (m)", yaxis_title="z right (m)", zaxis_title="y up (m)"),
    title="AH001 left leg + KompareXXVII 27-route BPAs — knee sweep, D(L) expansion",
    height=1050,
    updatemenus=[dict(type="buttons", direction="left", x=0.02, y=1.06,
                      buttons=[dict(label="Play", method="animate",
                                    args=[None, dict(frame=dict(duration=700, redraw=True), fromcurrent=True)]),
                               dict(label="Pause", method="animate",
                                    args=[[None], dict(mode="immediate", frame=dict(duration=0))])])],
    sliders=[dict(active=0, x=0.12, len=0.85, currentvalue=dict(prefix="knee: "),
                  steps=[dict(method="animate", label=f"{a}°", args=[[f"{a}"],
                          dict(mode="immediate", frame=dict(duration=0, redraw=True))]) for a in angles_list])])
out = os.path.join(REP, "robot_bpa27_knee_sweep.html")
fig.write_html(out, include_plotlyjs="cdn")
print("animated HTML:", out)
print("neutral->90deg route length changes (top 5):")
fd0, fdn = frame_data[0][2], frame_data[-1][2]
deltas = []
for (s1, _, _, L1, _, _), (s2, _, D2, L2, L02, _) in zip(fd0, fdn):
    deltas.append((s2, (L2 - L1) * 1000, D2 * 1000))
for s, dl, dd in sorted(deltas, key=lambda x: x[1])[:5]:
    print(f"  {s}: dL {dl:+.0f} mm, D@90deg {dd:.1f} mm")
