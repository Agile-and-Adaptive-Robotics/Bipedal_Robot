#!/usr/bin/env python3
"""mjcf2urdf.py - convert a MyoConverter-output MuJoCo MJCF (cvt3 XML) to URDF
for Simscape Multibody smimport.

Handles the exact subset used by MyoConverter Gait2392 output:
  * bodies with pos/quat (w x y z) transforms
  * joints: hinge (default type) / slide / multiple joints per body
  * MuJoCo joint anchor (pos) offset != body origin -> URDF can't put the
    revolute axis off the child-frame origin, so we insert the exact MuJoCo
    kinematic chain: parent -fixed- <name>_pre -revolute@jpos- <name>_anchor
    -fixed@-jpos- <name>
  * inertial: pos/quat + diaginertia; missing inertial -> tiny point mass
    (MyoConverter pathpoint bodies are massless)
  * geoms: mesh assets -> visual+collision; plane/other types skipped
  * sites, tendons, muscles, equality: ignored (SNS/BPA layer added later)

Usage: python mjcf2urdf.py <input.xml> <output.urdf> <robot_name>
Stdlib only (ElementTree + math)."""
import sys
import xml.etree.ElementTree as ET
import math


def f3(s, default=(0.0, 0.0, 0.0)):
    if s is None:
        return tuple(default)
    v = [float(x) for x in s.split()]
    if len(v) == 1:
        return (v[0], v[0], v[0])
    return tuple(v)


def f4(s):
    if s is None:
        return (1.0, 0.0, 0.0, 0.0)
    v = [float(x) for x in s.split()]
    return tuple(v)


def quat_to_rpy(q):
    """MuJoCo quat (w,x,y,z) -> URDF rpy (fixed-axis XYZ = ZYX intrinsic)."""
    w, x, y, z = q
    n = math.sqrt(w*w + x*x + y*y + z*z)
    if n < 1e-12:
        return (0.0, 0.0, 0.0)
    w, x, y, z = w/n, x/n, y/n, z/n
    roll = math.atan2(2*(w*x + y*z), 1 - 2*(x*x + y*y))
    siny = 2*(w*y - z*x)
    siny = max(-1.0, min(1.0, siny))
    pitch = math.asin(siny)
    yaw = math.atan2(2*(w*z + x*y), 1 - 2*(y*y + z*z))
    return (roll, pitch, yaw)


def origin_xml(xyz, rpy, indent):
    o = ET.Element('origin')
    o.set('xyz', '%.12g %.12g %.12g' % xyz)
    o.set('rpy', '%.12g %.12g %.12g' % rpy)
    return o


def add_origin(parent, xyz, rpy):
    parent.append(origin_xml(xyz, rpy, 0))


def sub3(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])


def neg3(a):
    return (-a[0], -a[1], -a[2])


def main():
    src, dst, name = sys.argv[1], sys.argv[2], sys.argv[3]
    tree = ET.parse(src)
    root = tree.getroot()

    # mesh assets: name -> file (forward slashes for URDF)
    meshes = {}
    for m in root.iter('mesh'):
        fn = m.get('file').replace('\\', '/')
        meshes[m.get('name')] = fn

    world = root.find('worldbody')
    if world is None:
        raise SystemExit('no worldbody')

    # find the root skeletal body: first <body> in worldbody that is not the
    # ground plane body. MyoConverter puts a 'ground' body with a plane geom.
    root_body = None
    for b in world.findall('body'):
        g = b.find('geom')
        if g is not None and g.get('type') == 'plane':
            continue
        root_body = b
        break
    if root_body is None:
        raise SystemExit('no root body found')

    urdf = ET.Element('robot', {'name': name})
    world_link = ET.SubElement(urdf, 'link', {'name': 'world'})

    stats = {'links': 0, 'joints': 0, 'meshes': 0, 'massless': 0}

    def link(nm):
        L = ET.SubElement(urdf, 'link', {'name': nm})
        stats['links'] += 1
        return L

    def add_inertial(L, body):
        inert = body.find('inertial')
        it = ET.SubElement(L, 'inertial')
        if inert is None:
            stats['massless'] += 1
            MV = ET.SubElement(it, 'mass')
            MV.set('value', '1e-6')
            add_origin(it, (0, 0, 0), (0, 0, 0))
            i = ET.SubElement(it, 'inertia')
            for k in ('ixx', 'ixy', 'ixz', 'iyy', 'iyz', 'izz'):
                i.set(k, '1e-9')
            return
        MV = ET.SubElement(it, 'mass')
        MV.set('value', inert.get('mass'))
        iq = f4(inert.get('quat'))
        add_origin(it, f3(inert.get('pos')), quat_to_rpy(iq))
        di = f3(inert.get('diaginertia'), (1e-6, 1e-6, 1e-6))
        i = ET.SubElement(it, 'inertia')
        i.set('ixx', '%.12g' % di[0])
        i.set('iyy', '%.12g' % di[1])
        i.set('izz', '%.12g' % di[2])
        i.set('ixy', '0'); i.set('ixz', '0'); i.set('iyz', '0')

    def add_geoms(L, body):
        for g in body.findall('geom'):
            gt = g.get('type', 'sphere')
            if gt == 'mesh':
                mf = meshes.get(g.get('mesh'))
                if mf is None:
                    continue
                for tag in ('visual', 'collision'):
                    V = ET.SubElement(L, tag)
                    add_origin(V, f3(g.get('pos')), quat_to_rpy(f4(g.get('quat'))))
                    G = ET.SubElement(V, 'geometry')
                    M = ET.SubElement(G, 'mesh', {'filename': mf})
                    stats['meshes'] += 1
            # plane/sphere/box/capsule collision geoms skipped: decoration

    def add_joint(nm, jtype, parent_l, child_l, xyz, rpy, axis, lim=None):
        J = ET.SubElement(urdf, 'joint', {'name': nm, 'type': jtype})
        add_origin(J, xyz, rpy)
        J.append(ET.Element('parent', {'link': parent_l}))
        J.append(ET.Element('child', {'link': child_l}))
        if jtype in ('revolute', 'prismatic', 'continuous'):
            A = ET.SubElement(J, 'axis', {'xyz': '%.12g %.12g %.12g' % axis})
        if lim is not None:
            Le = ET.SubElement(J, 'limit')
            Le.set('lower', '%.12g' % lim[0])
            Le.set('upper', '%.12g' % lim[1])
            Le.set('effort', '1e6')
            Le.set('velocity', '1e3')
        stats['joints'] += 1
        return J

    def convert_body(body, parent_link, parent_name):
        bname = body.get('name')
        bpos = f3(body.get('pos'))
        brpy = quat_to_rpy(f4(body.get('quat')))

        joints = body.findall('joint')
        # chain: parent_link -> (per joint: link_i) -> body frame link
        cur_link = parent_link
        cur_name = parent_name
        anchor_off = (0.0, 0.0, 0.0)   # -jpos correction applied at body link
        for ji, j in enumerate(joints):
            jtype_in = j.get('type', 'hinge')   # MuJoCo default = hinge
            jname = j.get('name', '%s_j%d' % (bname, ji))
            jpos = f3(j.get('pos'))
            axis = f3(j.get('axis'), (0, 0, 1))
            limited = j.get('limited', 'true')
            rng = f3(j.get('range'), (0, 0)) if j.get('range') else (0.0, 0.0)
            if jtype_in == 'hinge':
                jtype = 'continuous' if limited in ('false', '0') else 'revolute'
            elif jtype_in == 'slide':
                jtype = 'prismatic' if limited in ('false', '0') else 'prismatic'
            else:
                raise SystemExit('unsupported joint type %s on %s' % (jtype_in, jname))
            lim = None
            if jtype in ('revolute', 'prismatic') and rng != (0.0, 0.0):
                lim = rng
            if jtype == 'prismatic':
                lim = lim or (-1e3, 1e3)

            # joint anchor sits at jpos in the CHILD frame; URDF revolute axis
            # passes through the child-frame origin, so rebuild MuJoCo's exact
            # chain when jpos != 0.
            need_anchor = any(abs(v) > 1e-12 for v in jpos)
            first = (ji == 0)
            if need_anchor:
                # URDF can't hang a revolute axis off the child-frame origin,
                # so rebuild MuJoCo's exact chain:
                #   parent -fixed@bodypose- <pre> -revolute@jpos- <anchor>
                # and the body link later hangs at -jpos from <anchor>.
                pre_nm = '%s_%s_pre' % (bname, jname)
                L_pre = link(pre_nm)
                add_inertial(L_pre, ET.Element('dummy'))  # tiny mass
                JF = ET.SubElement(urdf, 'joint', {'name': pre_nm + '_fix', 'type': 'fixed'})
                add_origin(JF, bpos if first else (0, 0, 0),
                           brpy if first else (0, 0, 0))
                JF.append(ET.Element('parent', {'link': cur_name}))
                JF.append(ET.Element('child', {'link': pre_nm}))
                stats['joints'] += 1
                anc_nm = '%s_%s_anchor' % (bname, jname)
                L_anc = link(anc_nm)
                add_inertial(L_anc, ET.Element('dummy'))
                add_joint(jname + '_anchor', jtype, pre_nm, anc_nm, jpos, (0, 0, 0), axis, lim)
                cur_link, cur_name = L_anc, anc_nm
                anchor_off = neg3(jpos)
            else:
                # joint carries the body transform if it is the first one
                nxt_nm = '%s_%s_f' % (bname, jname)
                L_nxt = link(nxt_nm)
                add_inertial(L_nxt, ET.Element('dummy'))
                add_joint(jname, jtype, cur_name, nxt_nm,
                          bpos if first else (0, 0, 0),
                          brpy if first else (0, 0, 0), axis, lim)
                cur_link, cur_name = L_nxt, nxt_nm
            bpos, brpy = (0, 0, 0), (0, 0, 0)   # transform consumed by chain

        # final link = actual body frame (with -jpos correction if anchored)
        L_body = link(bname)
        add_inertial(L_body, body)
        add_geoms(L_body, body)
        if cur_link is not parent_link:
            JF = ET.SubElement(urdf, 'joint', {'name': bname + '_frame_fix', 'type': 'fixed'})
            add_origin(JF, anchor_off, (0, 0, 0))
            JF.append(ET.Element('parent', {'link': cur_name}))
            JF.append(ET.Element('child', {'link': bname}))
            stats['joints'] += 1
        else:
            # body with NO joints: fixed to parent at the body pose
            JF = ET.SubElement(urdf, 'joint', {'name': bname + '_fix', 'type': 'fixed'})
            add_origin(JF, bpos, brpy)
            JF.append(ET.Element('parent', {'link': cur_name}))
            JF.append(ET.Element('child', {'link': bname}))
            stats['joints'] += 1

        for child in body.findall('body'):
            convert_body(child, L_body, bname)

    # root body attaches straight to world; its own joint chain (for gait2392:
    # 3 pelvis slides + 3 pelvis hinges) provides the 6 DOF. A floating joint
    # would double-count DOF against the pelvis chain.
    rb_joints = root_body.findall('joint')
    if not rb_joints:
        # ungrounded root: hang it on a floating joint
        L_float = link(root_body.get('name') + '_floatbase')
        add_inertial(L_float, ET.Element('dummy'))
        JF = ET.SubElement(urdf, 'joint', {'name': root_body.get('name') + '_floating', 'type': 'floating'})
        add_origin(JF, (0, 0, 0), (0, 0, 0))
        JF.append(ET.Element('parent', {'link': 'world'}))
        JF.append(ET.Element('child', {'link': L_float.get('name')}))
        stats['joints'] += 1
        convert_body(root_body, L_float, L_float.get('name'))
    else:
        convert_body(root_body, world_link, 'world')

    ET.ElementTree(urdf).write(dst, encoding='utf-8', xml_declaration=True)
    print('mjcf2urdf: %s -> %s | links=%d joints=%d meshrefs=%d massless=%d'
          % (src, dst, stats['links'], stats['joints'], stats['meshes'], stats['massless']))


if __name__ == '__main__':
    main()
