"""Build the current AH001 CAD geometry with manually audited joint frames.

Preserve both four-bar knees. The URDF mesh package is an intermediate export;
MJCF is authoritative because URDF cannot represent the knee loop closures.
"""
import json
from pathlib import Path
import struct
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from add_torso import add_torso, values

HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent/'spinal/reports_20261027_robot_control'
# CAD and simulation share the requested OpenSim/ISB frame.
Q = np.eye(3)
_cad_pose=json.loads((AUDIT/'robot_mass_geometry.json').read_text())
_pelvis_component=next(c for c in _cad_pose['components'] if '/02_01_PE_' in c['name'])
_pelvis_transform=np.array(_pelvis_component['transform'])
CAD_R=_pelvis_transform[:9].reshape(3,3)
CAD_T=_pelvis_transform[9:12]

def cad_to_isb(point):
    return CAD_R@(np.array(point)-CAD_T)


def inertia_for_group(names, components, densities,frame_rotation=None,frame_translation=None):
    frame_rotation=CAD_R if frame_rotation is None else frame_rotation
    frame_translation=-CAD_R@CAD_T if frame_translation is None else frame_translation
    bodies = []
    for name in names:
        comp = components[name]
        t = np.array(comp['transform'])
        rotation = frame_rotation @ t[:9].reshape(3, 3).T
        for raw in comp['mass_properties_density_one']:
            raw = np.array(raw)
            density = densities[name]
            mass = raw[3]*density
            center = rotation @ raw[:3] + frame_rotation @ t[9:12]+frame_translation
            # SOLIDWORKS reports products of inertia; tensors use their negatives.
            xx, yy, zz, xy, xz, yz = raw[6:12]*density
            tensor = rotation @ np.array([[xx, -xy, -xz], [-xy, yy, -yz], [-xz, -yz, zz]]) @ rotation.T
            bodies.append((mass, center, tensor))
    total = sum(x[0] for x in bodies)
    center = sum(m*c for m,c,_ in bodies)/total
    tensor = sum(i + m*(np.dot(c-center,c-center)*np.eye(3)-np.outer(c-center,c-center)) for m,c,i in bodies)
    assert np.linalg.eigvalsh(tensor).min() > 0
    return total, center, tensor


def localize_mesh(source, output, origin, rotation=None,translation=None):
    if rotation is None:
        rotation=CAD_R
        translation=-CAD_R@CAD_T if translation is None else translation
    else:
        rotation=np.asarray(rotation)
        translation=np.zeros(3) if translation is None else translation
    raw = bytearray(source.read_bytes())
    count = struct.unpack_from('<I', raw, 80)[0]
    if len(raw) != 84 + 50*count:
        raise ValueError('Expected binary STL: '+str(source))
    triangles = np.ndarray(count, dtype=np.dtype([('normal','<f4',(3,)), ('vertices','<f4',(3,3)), ('attr','<u2')]), buffer=raw, offset=84)
    triangles['normal'] = triangles['normal'] @ rotation.T
    world = triangles['vertices'].astype(float) @ rotation.T+translation
    # SW2URDF exports metres. Reject accidental millimetre output.
    if np.max(np.abs(world)) > 3:
        raise ValueError('Mesh scale is inconsistent with metre CAD coordinates.')
    triangles['vertices'] = world - origin
    if np.linalg.det(rotation)<0:
        triangles['vertices'][:,[1,2]]=triangles['vertices'][:,[2,1]]
    files = []
    # MuJoCo 2.3.7 limits each binary STL asset to 200,000 faces. Partition
    # triangles without decimation so the CAD geometry remains unchanged.
    for start in range(0, count, 190000):
        chunk = triangles[start:start+190000]
        path = output.with_name(output.stem+'_%02d.stl' % len(files))
        path.write_bytes(bytes(80)+struct.pack('<I',len(chunk))+chunk.tobytes())
        files.append(path.name)
    return [world.min(axis=(0,1)).tolist(), world.max(axis=(0,1)).tolist()], files


def main():
    audit = json.loads((AUDIT/'robot_geometry_audit.json').read_text())
    cad = json.loads((AUDIT/'robot_mass_geometry.json').read_text())
    components = {x['name']:x for x in cad['components']}
    densities = {x['name']:x['density_kg_m3'] for x in audit['components']}
    group_of = {name:group for group,names in audit['rigid_groups'].items() for name in names}
    connections = {}
    for joint in audit['joints']:
        a, b = group_of[joint['a']], group_of[joint['b']]
        connections[frozenset((a,b))] = joint
    origins = {'pelvis':np.zeros(3)}
    links = []
    closures = []
    for side in ('r','l'):
        thigh, back, front, shank = [x+'_'+side for x in ('thigh','back','front','shank')]
        for parent,child,kind in [('pelvis',thigh,'ball'), (thigh,back,'hinge'), (back,shank,'hinge'), (thigh,front,'hinge')]:
            joint = connections[frozenset((parent,child))]
            if joint['separation_m'] > 1e-6:
                raise ValueError('CAD/import anchor mismatch.')
            origins[child] = cad_to_isb(joint['anchor_a'])
            links.append(dict(parent=parent, child=child, kind=kind,
                              axis=CAD_R@np.array(joint['axis_a'])))
        joint = connections[frozenset((front,shank))]
        closures.append(dict(front=front, shank=shank, anchor=cad_to_isb(joint['anchor_a'])))
    mesh_dir = HERE/'meshes'
    mesh_dir.mkdir(exist_ok=True)
    root = ET.Element('mujoco', model='AH001 four-bar lower body')
    ET.SubElement(root, 'compiler', angle='radian', inertiafromgeom='auto', meshdir='meshes')
    ET.SubElement(root, 'option', timestep='0.002', integrator='implicitfast', gravity='0 -9.81 0',
                  solver='Newton', iterations='100', jacobian='dense')
    ET.SubElement(ET.SubElement(root,'visual'),'global',offwidth='1200',offheight='1000')
    ET.SubElement(root.find('visual'),'headlight',ambient='.6 .6 .6',diffuse='.6 .6 .6')
    asset = ET.SubElement(root, 'asset')
    world = ET.SubElement(root, 'worldbody')
    ET.SubElement(world,'light',pos='0 3 2',dir='0 -1 0')
    ET.SubElement(world,'geom',name='ground',type='plane',quat='.7071067811865476 -.7071067811865476 0 0',size='3 3 .1',contype='2',conaffinity='1',rgba='.85 .85 .85 1')
    pelvis = ET.SubElement(world,'body',name='pelvis',pos='0 .99 0')
    ET.SubElement(pelvis,'freejoint',name='root')
    bodies = {'pelvis':pelvis}
    for link in links:
        body = ET.SubElement(bodies[link['parent']], 'body', name=link['child'],
                             pos=values(origins[link['child']]-origins[link['parent']]))
        attrs = dict(name=link['child']+'_joint',type=link['kind'],damping='1' if link['kind']=='ball' else '.15',armature='0')
        if link['kind'] == 'hinge':
            attrs['axis'] = values(link['axis'])
        ET.SubElement(body,'joint',**attrs)
        bodies[link['child']] = body
    properties = {}
    for name,names in audit['rigid_groups'].items():
        mass, center, tensor = inertia_for_group(names,components,densities)
        body = bodies[name]
        ET.SubElement(body,'inertial',mass=str(mass),pos=values(center-origins[name]),
                      fullinertia=values([tensor[0,0],tensor[1,1],tensor[2,2],tensor[0,1],tensor[0,2],tensor[1,2]]))
        source = HERE/'ah001_cad_meshes/meshes'/(name+'.STL')
        bounds, files = localize_mesh(source,mesh_dir/(name+'.stl'),origins[name])
        for index, filename in enumerate(files):
            meshname = name+'_%02d' % index
            ET.SubElement(asset,'mesh',name=meshname,file=filename)
            ET.SubElement(body,'geom',name=meshname+'_cad',type='mesh',mesh=meshname,
                          contype='0',conaffinity='0',
                          friction='.8 .005 .0001',rgba='.25 .3 .35 1')
        properties[name] = dict(mass_kg=mass, center_of_mass=center.tolist(), inertia_com=tensor.tolist(), mesh_bounds=bounds)
    from foot_geometry import add_feet
    foot_audit=add_feet(root,bodies,origins,properties,AUDIT,HERE)
    for side in ('r','l'):links.append(dict(parent='shank_'+side,child='foot_'+side,kind='fixed',axis=np.zeros(3)))
    pelvis.set('pos',values([0,-.0002-min(x['mesh_bounds'][0][1] for x in foot_audit.values()),0]))
    from robot_kinematics import camera_axes
    camera_position=np.array([1.25,1.25,1.65]);target_position=np.array([-.06,.64,0.])
    camera_rotation=camera_axes(camera_position,target_position)
    ET.SubElement(world,'camera',name='isb_preview',pos=values(camera_position),xyaxes=values(camera_rotation[:,:2].T.reshape(-1)))
    equality = ET.SubElement(root,'equality')
    for closure in closures:
        # Verified against compiled eq_data: anchor is local to body1.
        ET.SubElement(equality,'connect',name='knee_loop_'+closure['front'][-1],
                      body1=closure['front'],body2=closure['shank'],
                      anchor=values(closure['anchor']-origins[closure['front']]),solref='.01 1')
    ET.indent(root,space='  ')
    lower = HERE/'ah001_lower.xml'
    ET.ElementTree(root).write(lower,encoding='utf-8',xml_declaration=True)
    model = mujoco.MjModel.from_xml_path(str(lower))
    data = mujoco.MjData(model)
    mujoco.mj_forward(model,data)
    for index,closure in enumerate(closures):
        expected = closure['anchor']-origins[closure['front']]
        if not np.allclose(model.eq_data[index,:3], expected, atol=1e-7):
            raise AssertionError('Loop anchor frame mismatch: '+str(model.eq_data[index,:3])+' vs '+str(expected))
    torso_spec=json.loads((HERE/'torso_spec.json').read_text())
    hinge=np.array(torso_spec['lumbar_hinge_position_in_pelvis_m'])
    torso_audit = add_torso(lower, HERE/'ah001_with_torso.xml', 'pelvis', hinge, torso_spec)
    output = dict(material_and_hardware_assumptions=audit['warning'], groups=properties,
                  lumbar_hinge_pelvis_frame_m=hinge.tolist(),
                  hinge_location_status='upper rear center approximation, user approved adjustable',
                  torso=torso_audit, compiled_version=mujoco.__version__,
                  frame='x anterior, y superior, z right',feet=foot_audit,
                  constraints_at_reference=data.efc_pos.tolist(),
                  status='Geometry and mass build; current routed control model and run assessments are separate artifacts')
    (HERE/'ah001_build_audit.json').write_text(json.dumps(output,indent=2))
    # A manually articulated URDF companion for inspection/export tools. Its
    # tree omits two cut knee edges, recorded in the sidecar. Use MJCF for physics.
    urdf = ET.Element('robot',name='ah001_manual_tree')
    urdf.append(ET.Comment('Knee loop closures are in ah001_lower.xml; this URDF is an open tree. Hip spherical joints use three serial passive continuous axes.'))
    for name, prop in properties.items():
        link = ET.SubElement(urdf,'link',name=name)
        inertial = ET.SubElement(link,'inertial')
        ET.SubElement(inertial,'origin',xyz=values(np.array(prop['center_of_mass'])-origins[name]),rpy='0 0 0')
        ET.SubElement(inertial,'mass',value=str(prop['mass_kg']))
        tensor=np.array(prop['inertia_com'])
        ET.SubElement(inertial,'inertia',**{key:str(tensor[a,b]) for key,a,b in [('ixx',0,0),('iyy',1,1),('izz',2,2),('ixy',0,1),('ixz',0,2),('iyz',1,2)]})
        for mesh in root.findall("./asset/mesh"):
            if mesh.get('name').startswith(name+'_'):
                visual=ET.SubElement(link,'visual')
                ET.SubElement(visual,'origin',xyz='0 0 0',rpy='0 0 0')
                ET.SubElement(ET.SubElement(visual,'geometry'),'mesh',filename='meshes/'+mesh.get('file'))
    for item in links:
        parent,child=item['parent'],item['child']
        position=origins[child]-origins[parent]
        axes=[np.eye(3)[i] for i in range(3)] if item['kind']=='ball' else [item['axis']]
        for index,axis in enumerate(axes):
            target=child if index==len(axes)-1 else child+'_axis%d' % index
            if target != child:
                ET.SubElement(urdf,'link',name=target)
            joint=ET.SubElement(urdf,'joint',name=child+'_passive_%d' % index,type='fixed' if item['kind']=='fixed' else 'continuous')
            ET.SubElement(joint,'parent',link=parent)
            ET.SubElement(joint,'child',link=target)
            ET.SubElement(joint,'origin',xyz=values(position if index==0 else [0,0,0]),rpy='0 0 0')
            ET.SubElement(joint,'axis',xyz=values(axis))
            ET.SubElement(joint,'limit',effort='0',velocity='0')
            parent=target
    ET.indent(urdf,space='  ')
    ET.ElementTree(urdf).write(HERE/'ah001_manual_tree.urdf',encoding='utf-8',xml_declaration=True)
    (HERE/'ah001_loop_closures.json').write_text(json.dumps([dict(body1=x['front'],body2=x['shank'],anchor_body1_m=(x['anchor']-origins[x['front']]).tolist()) for x in closures],indent=2))
    print(json.dumps(torso_audit,indent=2))


if __name__ == '__main__':
    main()
