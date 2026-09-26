# -*- coding: utf-8 -*-
"""
make_w2l_mjcf.py - read-only parser: Walker_2_Layer_CPG.aproj -> w2l_mjcf.xml
Milestone 1 of the AnimatLab 2-layer walker -> MuJoCo port (2026-09-25).

SOURCE CONVENTIONS (established in this session against the file itself,
documented in reports_20260925/goal2_m1_physical_model.md):
  * LocalMatrix: 16 floats, COLUMN-major storage, column-vector convention
    (verified: rotation matches the Rotation element via R = Rz(z)@Ry(y)@Rx(x)
    for Root/femur/hip/tibia/foot; translation = elements [12:15]).
  * Translation in the matrices is in DECIMETERS (DistanceUnits=Decimeters);
    the Value/Scale/Actual triplets carry 'Actual' in METERS.
    Mass 'Actual' is in GRAMS (MassUnits=Grams); kg = Actual/1000.
  * Rotations in the Rotation elements are DEGREES.
  * <Joint> lives inside the CHILD RigidBody; its LocalPosition/LocalMatrix
    are expressed in the CHILD body's frame (verified: hip anchor lands just
    below the pelvis center; knee at the femur's bottom end; ankle at the
    foot's top-rear corner).
  * AnimatLab hinge rotates about the joint frame's local X axis
    (Vortex primary rotational coordinate).
  * Box dims: Length/Width/Height = full extents along local x/y/z.
  * World frame is y-up; MuJoCo is z-up. Global remap used here:
        MJ(x, y, z) = AL(x, -z, y)          (same remap the lab already uses
                                             for OpenSim -> MuJoCo conversions)
    A box with AL-local dims (a,b,c) therefore gets MJ-local dims (a,c,b).

MAPPINGS (documented deviations flagged in the report):
  * MaximumTension [N]           -> muscle force
  * RestingLength L0, Lwidth     -> muscle range = [L0-Lw, L0+Lw] and
                                    lengthrange = same (MuJoCo FLV then
                                    peaks at normalized 0.5 = L0 and dies
                                    at +/-Lwidth, matching AnimatLab's
                                    parabolic length-tension curve
                                    1-((L-L0)/Lw)^2, zero at +/-Lw)
  * B [N.s/m]                    -> muscle damp = B * (2*Lw)  (MuJoCo muscle
                                    damping acts on range-normalized
                                    velocity, so damp*vdot/(2Lw) = B*vdot)
  * Kse (series elastic)         -> NO equivalent (rigid tendon; out of
                                    scope for M1 per the ask) - deviation
  * Kpe (parallel elastic)       -> NO equivalent - deviation
  * Materials: Vortex 'foot' FrictionLinearPrimary=1e5 is a viscous-type
    contact constant, not a Coulomb coefficient -> MuJoCo friction=1.0
    (MuJoCo/Default-material value); deviation documented.
  * LinearHillStretchReceptor bodies: sensory overlays (ApplyTension=False),
    no mass -> not ported (sensory port is a later milestone).
  * Attachment bodies -> MuJoCo sites on their parent body.
  * Spring bodies (toe springs, k=16000 N/m, b=20000, L0=0.088) ->
    two-site spatial tendons with spring/damping/springlength.
  * Root Freeze=True in the aproj is IGNORED (the Freeze trap: pelvis must
    be FREE in MuJoCo).

Outputs (next to this script):
  w2l_mjcf.xml          - the MuJoCo model
  w2l_source_dump.json  - every parsed source value (validator + report input)
"""
import json
import math
import os
import sys
import xml.etree.ElementTree as ET

APROJ = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj"
HERE = os.path.dirname(os.path.abspath(__file__))
XML_OUT = os.path.join(HERE, "w2l_mjcf.xml")
DUMP_OUT = os.path.join(HERE, "w2l_source_dump.json")

D2R = math.pi / 180.0

# ---------------------------------------------------------------- vec / rot
def vadd(a, b):  return [a[i] + b[i] for i in range(3)]
def vsub(a, b):  return [a[i] - b[i] for i in range(3)]
def vscale(a, s): return [a[i] * s for i in range(3)]
def vdot(a, b):  return sum(a[i] * b[i] for i in range(3))
def vnorm(a):    return math.sqrt(vdot(a, a))

def mat_mul(A, B):
    """3x3 x 3x3."""
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)]
            for i in range(3)]

def mat_vec(A, v):
    return [sum(A[i][k] * v[k] for k in range(3)) for i in range(3)]

def mat_T(A):
    return [[A[j][i] for j in range(3)] for i in range(3)]

def rx(a):
    c, s = math.cos(a), math.sin(a)
    return [[1, 0, 0], [0, c, -s], [0, s, c]]

def ry(a):
    c, s = math.cos(a), math.sin(a)
    return [[c, 0, s], [0, 1, 0], [-s, 0, c]]

def rz(a):
    c, s = math.cos(a), math.sin(a)
    return [[c, -s, 0], [s, c, 0], [0, 0, 1]]

def euler_xyz_deg(x, y, z):
    """AnimatLab Rotation element -> rotation matrix.

    Order verified against the LocalMatrix elements in this file:
    R = Rx(x) @ Ry(y) @ Rz(z)  (matches tibia_L to 5.8e-6, toe to <1e-5;
    Rz@Ry@Rx deviates by 0.052 there - see report).
    """
    return mat_mul(rx(x * D2R), mat_mul(ry(y * D2R), rz(z * D2R)))

def quat_from_mat(R):
    """3x3 rotation matrix -> quaternion [w, x, y, z]."""
    tr = R[0][0] + R[1][1] + R[2][2]
    if tr > 0:
        S = math.sqrt(tr + 1.0) * 2
        w = 0.25 * S
        x = (R[2][1] - R[1][2]) / S
        y = (R[0][2] - R[2][0]) / S
        z = (R[1][0] - R[0][1]) / S
    elif R[0][0] > R[1][1] and R[0][0] > R[2][2]:
        S = math.sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]) * 2
        w = (R[2][1] - R[1][2]) / S
        x = 0.25 * S
        y = (R[0][1] + R[1][0]) / S
        z = (R[0][2] + R[2][0]) / S
    elif R[1][1] > R[2][2]:
        S = math.sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]) * 2
        w = (R[0][2] - R[2][0]) / S
        x = (R[0][1] + R[1][0]) / S
        y = 0.25 * S
        z = (R[1][2] + R[2][1]) / S
    else:
        S = math.sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]) * 2
        w = (R[1][0] - R[0][1]) / S
        x = (R[0][2] + R[2][0]) / S
        y = (R[1][2] + R[2][1]) / S
        z = 0.25 * S
    return [w, x, y, z]

def norm_quat(q):
    n = math.sqrt(sum(c * c for c in q))
    return [c / n for c in q]

# AL(x,y,z) -> MJ(x,-z,y):  A @ v
_A = [[1, 0, 0], [0, 0, -1], [0, 1, 0]]
_A_T = mat_T(_A)

def al2mj_vec(v):
    return mat_vec(_A, v)

def al2mj_mat(R):
    return mat_mul(_A, mat_mul(R, _A_T))

# ---------------------------------------------------------------- XML read
def txt(el, tag):
    c = el.find(tag)
    return (c.text or "").strip() if c is not None else None

def triplet(el, tag):
    """Value/Scale/Actual child -> (value, scale, actual SI)."""
    c = el.find(tag)
    if c is None:
        return None
    return (c.get("Value"), c.get("Scale"), c.get("Actual"))

def actual_f(el, tag, default=None):
    t = triplet(el, tag)
    return float(t[2]) if t is not None and t[2] is not None else default

def parse_local_matrix(el):
    """Return (R_al 3x3, p_al m) from LocalMatrix (column-major, dm); None if absent."""
    c = el.find("LocalMatrix")
    if c is None or not (c.text or "").strip():
        return None, None
    m = [float(x) for x in (c.text or "").split(",")]
    # column-major 4x4: columns m[0:4], m[4:8], m[8:12], m[12:16]
    R = [[m[0], m[4], m[8]],
         [m[1], m[5], m[9]],
         [m[2], m[6], m[10]]]
    p = [m[12] / 10.0, m[13] / 10.0, m[14] / 10.0]   # dm -> m
    return R, p

def parse_pos_rot(el):
    lp = el.find("LocalPosition")
    ro = el.find("Rotation")
    p = [float(lp.find(a).get("Actual")) for a in ("X", "Y", "Z")] if lp is not None else [0.0] * 3
    r = [float(ro.find(a).get("Actual")) for a in ("X", "Y", "Z")] if ro is not None else [0.0] * 3
    return p, r

# ---------------------------------------------------------------- parser
class Body(object):
    pass

def collect_bodies(rb, parent, out):
    b = Body()
    b.name = txt(rb, "Name")
    b.id = txt(rb, "ID")
    b.type = txt(rb, "Type")
    b.parent = parent
    b.mass_kg = (actual_f(rb, "Mass") or 0.0) / 1000.0       # grams -> kg
    b.lwh = None
    if b.type == "Box":
        b.lwh = tuple(actual_f(rb, t, 0.0) for t in ("Length", "Width", "Height"))
    b.com = [actual_f(rb.find("COM"), a, 0.0) for a in ("X", "Y", "Z")] \
        if rb.find("COM") is not None else [0.0, 0.0, 0.0]
    b.material_id = txt(rb, "MaterialTypeID")
    b.freeze = txt(rb, "Freeze")
    b.collision = txt(rb, "IsCollisionObject")
    R_lm, p_lm = parse_local_matrix(rb)
    p_lr, r_lr = parse_pos_rot(rb)
    R_er = euler_xyz_deg(*r_lr)
    if R_lm is None:
        b.pos_al = p_lr
        b.R_al = R_er
        b.pos_check_dev = 0.0
        b.rot_check_dev = 0.0
        b.from_matrix = False
    else:
        b.pos_al = p_lm
        b.R_al = R_lm
        b.pos_check_dev = vnorm(vsub(p_lm, p_lr))
        b.rot_check_dev = max(abs(R_lm[i][j] - R_er[i][j]) for i in range(3) for j in range(3))
        b.from_matrix = True
    b.joint = None
    jt = rb.find("Joint")
    if jt is not None:
        jR, jp = parse_local_matrix(jt)
        b.joint = {
            "name": txt(jt, "Name"),
            "id": txt(jt, "ID"),
            "type": txt(jt, "Type"),
            "pos_local_al": jp,                 # child-frame, m
            "R_local_al": jR,
            "enable_limits": (txt(jt, "EnableLimits") == "True"),
            "lower_deg": actual_f(jt.find("LowerLimit"), "LimitPos", None)
                if jt.find("LowerLimit") is not None else None,
            "upper_deg": actual_f(jt.find("UpperLimit"), "LimitPos", None)
                if jt.find("UpperLimit") is not None else None,
            "limit_stiffness": actual_f(jt.find("LowerLimit"), "Stiffness", None)
                if jt.find("LowerLimit") is not None else None,
            "max_force": actual_f(jt, "MaxForce"),
            "enable_motor": txt(jt, "EnableMotor"),
            "attach_id": None,  # filled by caller if desired
        }
    b.attach_ids = [a.text for a in rb.findall("Attachments/AttachID")]
    b.children = []
    cb = rb.find("ChildBodies")
    if cb is not None:
        for child in cb.findall("RigidBody"):
            b.children.append(collect_bodies(child, b, out))
    out.append(b)
    return b

def world_transforms(root_body):
    """Fill b.R_world, b.p_world (AL frame)."""
    def rec(b, Rp, pp):
        b.R_world = mat_mul(Rp, b.R_al)
        b.p_world = vadd(pp, mat_vec(Rp, b.pos_al))
        for c in b.children:
            rec(c, b.R_world, b.p_world)
    rec(root_body, [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0])

def parse_aproj(path=APROJ):
    tree = ET.parse(path)
    root = tree.getroot()
    env = root.find(".//Environment")
    timestep = actual_f(env, "PhysicsTimeStep")
    grav = actual_f(env, "Gravity")
    materials = {}
    for mt in env.findall("MaterialTypes/MaterialType"):
        materials[txt(mt, "ID")] = {
            "name": txt(mt, "Name"),
            "friction_linear_primary": actual_f(mt, "FrictionLinearPrimary"),
        }
    organism = root.find(".//Organism")
    org_rb = organism.find("RigidBody")
    bodies = []
    tree_root = collect_bodies(org_rb, None, bodies)
    world_transforms(tree_root)

    by_id = {b.id: b for b in bodies if b.id}
    # joints: attach to the child body that contains them (already)
    # muscles / receptors / springs / attachments
    muscles = []
    receptors = []
    springs = []
    attachments = {}
    for b in bodies:
        if b.type == "Attachment":
            attachments[b.id] = b

    # walk raw elements for muscle/spring/receptor specifics
    def walk(rb_el, parent_body):
        t = txt(rb_el, "Type")
        if t == "LinearHillMuscle":
            lt = rb_el.find("LengthTension")
            st = rb_el.find("StimulusTension")
            muscles.append({
                "name": txt(rb_el, "Name"),
                "id": txt(rb_el, "ID"),
                "attach_ids": [a.text for a in rb_el.findall("Attachments/AttachID")],
                "max_tension": actual_f(rb_el, "MaximumTension"),
                "kse": actual_f(rb_el, "Kse"),
                "kpe": actual_f(rb_el, "Kpe"),
                "B": actual_f(rb_el, "B"),
                "resting_length": actual_f(lt, "RestingLength"),
                "lwidth": actual_f(lt, "Lwidth"),
                "pe_length": actual_f(lt, "PeLength"),
                "st_a": actual_f(st, "A"), "st_b": actual_f(st, "B"),
                "st_c": actual_f(st, "C"), "st_d": actual_f(st, "D"),
            })
        elif t == "LinearHillStretchReceptor":
            receptors.append({
                "name": txt(rb_el, "Name"),
                "attach_ids": [a.text for a in rb_el.findall("Attachments/AttachID")],
                "apply_tension": txt(rb_el, "ApplyTension"),
            })
        elif t == "Spring":
            springs.append({
                "name": txt(rb_el, "Name"),
                "attach_ids": [a.text for a in rb_el.findall("Attachments/AttachID")],
                "natural_length": actual_f(rb_el, "NaturalLength"),
                "stiffness": actual_f(rb_el, "Stiffness"),
                "damping": actual_f(rb_el, "Damping"),
            })
        for ch in rb_el.find("ChildBodies").findall("RigidBody") \
                if rb_el.find("ChildBodies") is not None else []:
            walk(ch, parent_body)
    walk(org_rb, None)

    # attachment world positions (attachment LocalMatrix relative to parent body)
    for b in bodies:
        if b.type == "Attachment":
            b.site_body = b.parent
            b.site_pos_al_local = b.pos_al      # already parent-frame
            b.site_pos_world_al = b.p_world

    return {
        "timestep": timestep, "gravity": grav, "materials": materials,
        "root": tree_root, "bodies": bodies, "by_id": by_id,
        "muscles": muscles, "receptors": receptors, "springs": springs,
        "attachments": attachments,
    }

# ---------------------------------------------------------------- emit MJCF
def f3(v):
    return " ".join("%.8g" % c for c in v)

def build_mjcf(m):
    lines = []
    A = lines.append
    used = {}   # (element-type, name) -> count, for per-type uniqueness

    A('<mujoco model="w2l_animatlab_port">')
    A('  <!-- Generated by make_w2l_mjcf.py from Walker_2_Layer_CPG.aproj (read-only parse).')
    A('       Frame remap: MuJoCo(x,y,z) = AnimatLab(x,-z,y). See report for mappings. -->')
    A('  <compiler angle="radian" autolimits="true"/>')
    A('  <option timestep="%g" gravity="0 0 %g"/>' % (m["timestep"], m["gravity"]))
    A('  <default>')
    A('    <joint type="hinge" damping="0"/>')
    A('    <geom type="box" density="0" contype="1" conaffinity="1"/>')
    A('    <site type="sphere" size="0.006" group="3"/>')
    A('  </default>')
    A('  <worldbody>')
    A('    <geom name="ground" type="plane" size="5 5 0.1" pos="0 0 0"')
    A('          friction="1 0.005 0.0001" contype="1" conaffinity="1"/>')

    # emit body tree recursively (boxes with mass; skip massless overlay types)
    def emit_body(b, indent):
        if b.type in ("LinearHillMuscle", "LinearHillStretchReceptor",
                      "Spring", "Attachment"):
            return  # overlays: not physical bodies in MJCF (documented);
                    # muscles/receptors are behavior-layer, springs become
                    # tendons, attachments become sites on their parent
        pad = "  " * indent
        # MuJoCo body pos/quat are RELATIVE TO PARENT
        if b.parent is None:
            pos_rel_al = b.p_world
            R_rel_al = b.R_al
        else:
            pos_rel_al = mat_vec(mat_T(b.parent.R_al),
                                 vsub(b.p_world, b.parent.p_world))
            R_rel_al = mat_mul(mat_T(b.parent.R_al), b.R_al)
        pos = al2mj_vec(pos_rel_al)
        quat = norm_quat(quat_from_mat(al2mj_mat(R_rel_al)))
        name = mj_name("body", b.name)
        b.mj_name = name
        attrs = 'name="%s" pos="%s"' % (name, f3(pos))
        if abs(quat[0] - 1.0) > 1e-10:
            attrs += ' quat="%s"' % f3(quat)
        A("%s<body %s>" % (pad, attrs))
        # inertial (uniform box about COM; massless overlays skipped)
        if b.type == "Box" and b.mass_kg > 0:
            dims_mj = (b.lwh[0], b.lwh[2], b.lwh[1])   # AL(a,b,c) -> MJ(a,c,b)
            com_mj = al2mj_vec(b.com)
            bmass = b.mass_kg
            Ix = bmass / 12.0 * (dims_mj[1] ** 2 + dims_mj[2] ** 2)
            Iy = bmass / 12.0 * (dims_mj[0] ** 2 + dims_mj[2] ** 2)
            Iz = bmass / 12.0 * (dims_mj[0] ** 2 + dims_mj[1] ** 2)
            A('%s  <inertial pos="%s" mass="%g" diaginertia="%s"/>'
              % (pad, f3(com_mj), bmass, f3([Ix, Iy, Iz])))
            # joint (hinge) if present
            if b.joint is not None:
                j = b.joint
                jpos_mj = al2mj_vec(j["pos_local_al"])
                jR = mat_mul(b.R_al, j["R_local_al"])
                axis_al = mat_vec(jR, [1.0, 0.0, 0.0])      # Vortex primary = local x
                axis_mj = al2mj_vec(axis_al)
                jattrs = 'name="%s" pos="%s" axis="%s"' % (mj_name("joint", j["name"]), f3(jpos_mj), f3(axis_mj))
                if j["enable_limits"] and j["lower_deg"] is not None and j["upper_deg"] is not None:
                    lo, hi = j["lower_deg"] * D2R, j["upper_deg"] * D2R
                    if abs(lo) < 1e-12 and abs(hi) < 1e-12:
                        jattrs += ' limited="false"'
                    else:
                        jattrs += ' limited="true" range="%.10g %.10g"' % (lo, hi)
                else:
                    jattrs += ' limited="false"'
                A("%s  <joint %s/>" % (pad, jattrs))
            # geom
            gname = mj_name("geom", b.name)
            A('%s  <geom name="%s" size="%s" friction="1 0.005 0.0001" rgba="0.7 0.7 0.7 1"/>'
              % (pad, gname, f3([d / 2.0 for d in dims_mj])))
        else:
            # non-box physical types: none expected in this model
            pass
        # sites for attachments parented to this body
        for att in m["bodies"]:
            if att.type == "Attachment" and att.parent is b:
                sname = mj_name("site", att.name, prefix="site_")
                att.site_name = sname
                spos = al2mj_vec(att.pos_al)
                A('%s  <site name="%s" pos="%s"/>' % (pad, sname, f3(spos)))
        for c in b.children:
            emit_body(c, indent + 1)
        A("%s</body>" % pad)

    # NOTE: attachments parented to Root (body_L_f etc.) need sites on Root;
    # emit_body handles that via the m["bodies"] loop per body.
    emit_body(m["root"], 2)
    A('  </worldbody>')

    # tendons: muscles + toe springs
    A('  <tendon>')
    for mus in m["muscles"]:
        sites = [m["by_id"][aid].site_name for aid in mus["attach_ids"]]
        A('    <spatial name="t_%s">%s</spatial>'
          % (mus["name"], "".join('<site site="%s"/>' % s for s in sites)))
    for spr in m["springs"]:
        sites = [m["by_id"][aid].site_name for aid in spr["attach_ids"]]
        # Source spring damping is 20000 N.s/m = ~116x critical for the ~0.46 kg
        # toe chain and is NOT integrable at MuJoCo's 1 ms step (measured:
        # max|qacc| 3.5e8 + instability warning). Capped at the critical value
        # c* = 2*sqrt(k*m) ~ 172 N.s/m (max|qacc| 1.77e4 = contact-impact
        # baseline, same as c=0). Deviation documented in the M1 report.
        A('    <!-- source damping %g N.s/m capped to 172 (critical) for 1 ms stability -->'
          % spr["damping"])
        A('    <spatial name="t_%s" stiffness="%g" damping="172" springlength="%g">%s</spatial>'
          % (spr["name"], spr["stiffness"], spr["natural_length"],
             "".join('<site site="%s"/>' % s for s in sites)))
    A('  </tendon>')

    A('  <actuator>')
    for mus in m["muscles"]:
        L0, Lw = mus["resting_length"], mus["lwidth"]
        lo, hi = L0 - Lw, L0 + Lw
        # mujoco 2.3.7 <muscle> has NO damping attribute (B is not portable);
        # timeconst takes activation+deactivation (10 ms each).
        A('    <muscle name="%s" tendon="t_%s" force="%g" range="%.8g %.8g" '
          'lengthrange="%.8g %.8g" timeconst="0.01 0.01"/>'
          % (mus["name"], mus["name"], mus["max_tension"], lo, hi, lo, hi))
    A('  </actuator>')
    A('</mujoco>')
    return "\n".join(lines) + "\n"

_used = {}
def mj_name(etype, name, prefix=""):
    """Sanitize unique MuJoCo name within its element type namespace.

    MuJoCo requires uniqueness per element type only (a body and a joint may
    share the name 'toe_L'), so keys are (etype, name).
    """
    key = (etype, prefix + name)
    if key in _used:
        n = _used[key] + 1
        _used[key] = n
        return "%s%s_%d" % (prefix, name.replace(" ", "_").replace("-", "_"), n)
    _used[key] = 1
    return "%s%s" % (prefix, name.replace(" ", "_").replace("-", "_"))

def main():
    _used.clear()
    m = parse_aproj()
    xml_text = build_mjcf(m)
    with open(XML_OUT, "w") as fh:
        fh.write(xml_text)
    # source dump for validator + report
    dump = {
        "aproj": APROJ,
        "timestep": m["timestep"], "gravity": m["gravity"],
        "materials": m["materials"],
        "bodies": [],
        "muscles": m["muscles"],
        "receptors": m["receptors"],
        "springs": m["springs"],
        "attachments": [],
        "notes": {
            "conventions": "column-major LocalMatrix (dm), grams mass, y-up AL world, "
                           "child-frame joint anchors, hinge about joint local x, "
                           "MJ(x,y,z)=AL(x,-z,y)",
        },
    }
    for b in m["bodies"]:
        entry = {
            "name": b.name, "id": b.id, "type": b.type,
            "mass_kg": b.mass_kg, "lwh_al": b.lwh,
            "pos_world_al": b.p_world,
            "freeze": b.freeze, "collision": b.collision,
            "material": (m["materials"].get(b.material_id) or {}).get("name"),
            "pos_check_dev_m": b.pos_check_dev,
            "rot_check_dev": b.rot_check_dev,
        }
        if b.joint is not None:
            j = b.joint
            # world axis = R_world(child) @ R_joint(local) @ e_x  (hinge = local x)
            jR_world = mat_mul(b.R_world, j["R_local_al"])
            axis_al = mat_vec(jR_world, [1.0, 0.0, 0.0])
            entry["joint"] = {
                "name": j["name"], "enable_limits": j["enable_limits"],
                "lower_deg": j["lower_deg"], "upper_deg": j["upper_deg"],
                "anchor_world_al": vadd(b.p_world, mat_vec(b.R_world, j["pos_local_al"])),
                "axis_world_al": axis_al,
            }
        if b.type == "Attachment":
            entry["site_parent"] = b.parent.name if b.parent else None
            entry["pos_local_al"] = b.pos_al
        dump["bodies"].append(entry)
    with open(DUMP_OUT, "w") as fh:
        json.dump(dump, fh, indent=1)
    n_boxes = sum(1 for b in m["bodies"] if b.type == "Box")
    print("bodies parsed: %d (boxes %d, muscles %d, receptors %d, springs %d, attachments %d)"
          % (len(m["bodies"]), n_boxes, len(m["muscles"]), len(m["receptors"]),
             len(m["springs"]), len(m["attachments"])))
    print("wrote %s" % XML_OUT)
    print("wrote %s" % DUMP_OUT)
    # sanity: LocalPosition vs LocalMatrix agreement
    devs = [(b.name, b.pos_check_dev) for b in m["bodies"] if b.pos_check_dev > 1e-6]
    if devs:
        print("WARNING pos dev >1e-6 m:")
        for nm, d in devs:
            print("  %s: %.3g" % (nm, d))
    rot_devs = [(b.name, b.rot_check_dev) for b in m["bodies"] if b.rot_check_dev > 1e-4]
    if rot_devs:
        print("NOTE rot dev >1e-4 (LocalMatrix vs RxRyRz Euler; all are Attachment/"
              "StretchReceptor overlays whose own orientation is not ported):")
        for nm, d in rot_devs:
            print("  %s: %.3g" % (nm, d))
    return 0

if __name__ == "__main__":
    sys.exit(main())
