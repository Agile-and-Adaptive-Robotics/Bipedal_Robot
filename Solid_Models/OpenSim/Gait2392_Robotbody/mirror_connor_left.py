"""
Build ConnorBipedal_symmetric.osim from ConnorBipedal.osim by mirroring the
RIGHT muscle paths onto the LEFT side (Ben, 2026-09-22: "copy them to the
left"). The ConnorBipedal models are the UNTOUCHED master's-era robot (repo
import May 2023); the RH side carries the custom BPA/muscle routings matching
the 27-actuator set of RobotPAMCalculationOptimization.m (Connor Morrow's
PamData machinery). The LH side still carries pristine stock gait2392 geometry.

Same mirror rules as repair_robotbody_muscles.py (verified against
gait2392_simbody.osim, whose left side is exactly the mirror of its right):
  - PathPoint / ConditionalPathPoint location: (x, y, z) -> (x, y, -z)
  - ConditionalPathPoint: socket_coordinate _r -> _l; 'range' copied unchanged
  - MovingPathPoint: x/y SimmSplines copied unchanged; z_location <y> negated;
    socket_*_coordinate _r -> _l
  - socket_parent_frame and point names: _r -> _l

Differences from repair_robotbody_muscles.py:
  - source is ConnorBipedal.osim; output is a NEW file
    (ConnorBipedal_symmetric.osim) — the originals are untouched
  - appliesForce flags are NOT touched (Ben's 16 right-side test disables
    stay; geometry/volume work does not depend on them)

Run:  python mirror_connor_left.py
"""

import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).parent
SRC = HERE / "ConnorBipedal.osim"
DST = HERE / "ConnorBipedal_symmetric.osim"

PTAGS = ("PathPoint", "ConditionalPathPoint", "MovingPathPoint")


def load(path):
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    return ET.parse(path, parser=parser).getroot()


def swap(text):
    """_r -> _l in a socket path or point name."""
    return text.replace("_r", "_l") if text else text


def fmt(v):
    v = float(v)
    if v == 0.0:
        v = 0.0  # normalize -0.0
    return repr(v)


def mirror_geom_path(gp):
    """Deep-copy a right-side GeometryPath into a mirrored left-side one."""
    gp2 = ET.fromstring(ET.tostring(gp))
    for pp in gp2.iter():
        if pp.tag in PTAGS:
            name = pp.get("name") or ""
            head, sep, tail = name.partition("-P")
            if sep and head.endswith("_r"):
                pp.set("name", head[:-2] + "_l-P" + tail)
            for s in pp.findall("socket_parent_frame"):
                s.text = swap(s.text)
            loc = pp.find("location")
            if loc is not None and loc.text:
                x, y, z = (float(v) for v in loc.text.split())
                loc.text = f"{fmt(x)} {fmt(y)} {fmt(-z)}"
            for s in pp.findall("socket_coordinate"):
                s.text = swap(s.text)
            for s in ("socket_x_coordinate", "socket_y_coordinate",
                      "socket_z_coordinate"):
                for el in pp.findall(s):
                    el.text = swap(el.text)
            if pp.tag == "MovingPathPoint":
                zl = pp.find("z_location/SimmSpline")
                if zl is not None:
                    yt = zl.find("y")
                    if yt is not None and yt.text:
                        vals = [float(v) for v in yt.text.split()]
                        yt.text = " ".join(fmt(-v) for v in vals)
    return gp2


def geom_path(muscle):
    for child in muscle:
        if child.tag == "GeometryPath":
            return child
    raise KeyError("GeometryPath not found")


def muscles(root):
    return {m.get("name"): m for m in root.iter("Thelen2003Muscle")}


rb = load(SRC)
M = muscles(rb)

n_mirrored = 0
for name, m in M.items():
    if not name.endswith("_r"):
        continue
    l = name[:-1] + "l"
    ml = M[l]
    old = geom_path(ml)
    idx = list(ml).index(old)
    ml.remove(old)
    ml.insert(idx, mirror_geom_path(geom_path(m)))
    n_mirrored += 1

# flags deliberately untouched (see docstring)
n_off = sum(1 for m in M.values()
            if (m.find("appliesForce") is not None
                and (m.find("appliesForce").text or "").strip() == "false"))

tree = ET.ElementTree(rb)
ET.indent(tree, space="\t")
tree.write(DST, encoding="UTF-8", xml_declaration=True)
print(f"wrote {DST} (mirrored {n_mirrored} right->left GeometryPaths; "
      f"{n_off} appliesForce=false flags kept as-is)")

# ---------------- validation ----------------
def paths_of(root):
    out = {}
    for m in root.iter("Thelen2003Muscle"):
        pts = []
        for pp in m.iter():
            if pp.tag in PTAGS:
                loc = pp.find("location")
                locv = tuple(float(v) for v in loc.text.split()) \
                    if loc is not None and loc.text else None
                body = None
                for s in pp.findall("socket_parent_frame"):
                    body = (s.text or "").strip()
                pts.append((pp.tag, body, locv))
        out[m.get("name")] = pts
    return out


def zmirror(pts):
    out = []
    for (t, b, l) in pts:
        b2 = swap(b) if b else b
        l2 = None if l is None else (l[0], l[1], -l[2])
        out.append((t, b2, l2))
    return out


fixed = load(DST)
P = paths_of(fixed)
bad = [n for n in P if n.endswith("_r") and P[n[:-1] + "l"] != zmirror(P[n])]
short = [n for n, p in P.items() if len(p) < 2]
print(f"non-mirrored pairs: {len(bad)} {bad}")
print(f"muscles with <2 path points: {len(short)} {short}")

src_pts = paths_of(load(SRC))
unchanged = sum(1 for n in P if P[n] == src_pts[n])
print(f"paths identical to source (flags/params untouched model-wide): "
      f"{unchanged}/{len(P)} checked by point comparison")
assert not bad and not short, "validation failed"
print("VALIDATION OK")
