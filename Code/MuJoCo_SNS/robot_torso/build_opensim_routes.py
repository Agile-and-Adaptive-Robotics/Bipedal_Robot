"""Compile robotbody routes and 27 reference tension actuators per side."""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from add_torso import values

HERE=Path(__file__).resolve().parent


def build(specification,output,include_reference=False):
    configuration=json.loads(Path(specification).read_text())
    tree=ET.parse(HERE/'ah001_with_torso.xml')
    root=tree.getroot()
    old=root.find('tendon')
    old_sites={s.get('site') for s in old.findall('.//site')} if old is not None else set()
    if old is not None:
        root.remove(old)
    for body in root.findall('.//body'):
        for site in list(body.findall('site')):
            if site.get('name') in old_sites:
                body.remove(site)
    bodies={x.get('name'):x for x in root.findall('.//body')}
    tendons=ET.SubElement(root,'tendon')
    actuators=ET.SubElement(root,'actuator')
    routes_used=[r for r in configuration['routes'] if include_reference or r.get('actuator_assignment') is not None]
    for route in routes_used:
        enabled=[]
        for point in route['points']:
            if point['body'] not in bodies:
                raise ValueError('Unknown target body '+point['body'])
            position=np.asarray(point['position_m'],dtype=float)
            if position.shape!=(3,) or not np.isfinite(position).all():
                raise ValueError('Invalid site position '+point['id'])
            ET.SubElement(bodies[point['body']],'site',name=point['id'],pos=values(position),size='.0015',group='4')
            # Allocate every site in the path's fixed storage. Runtime membership
            # compacts active sites within this allocation without touching others.
            enabled.append(point)
        if len(enabled)<2:
            raise ValueError('Route requires at least two enabled points: '+route['name'])
        spatial=ET.SubElement(tendons,'spatial',name=route['name'],width='.001',stiffness='0',damping='0',rgba='.7 .2 .15 .6',group='2')
        for point in enabled:
            ET.SubElement(spatial,'site',site=point['id'])
        assignment=route.get('actuator_assignment')
        if assignment is not None:
            capacity=float(assignment['reference_capacity_N'])
            ET.SubElement(actuators,'general',name=route['name'],tendon=route['name'],
                          dyntype='filter',dynprm='.03',gaintype='fixed',gainprm=str(-capacity),biastype='none',
                          ctrllimited='true',ctrlrange='0 1',forcelimited='true',forcerange=values([-capacity,0]))
    root.append(ET.Comment('Robotbody right paths mirrored left; moving/conditional paths require DynamicRoutes updates. Actuators are first-order reference tensions capped by the existing 27-actuator design MIFs, not calibrated BPA models.'))
    ET.indent(tree,space='  ')
    tree.write(output,encoding='utf-8',xml_declaration=True)
    model=mujoco.MjModel.from_xml_path(str(output))
    data=mujoco.MjData(model)
    mujoco.mj_forward(model,data)
    from dynamic_routes import DynamicRoutes
    dynamic=DynamicRoutes(model,configuration)
    dynamic.update(data)
    assert model.nu==54 and np.isfinite(data.ten_length).all()
    assert np.isclose(model.body_mass.sum(),75.,rtol=0,atol=1e-9)
    rows=np.flatnonzero(data.efc_type==int(mujoco.mjtConstraint.mjCNSTR_EQUALITY))
    jacobian=data.efc_J.reshape(data.nefc,model.nv)[rows]
    _,singular,vh=np.linalg.svd(jacobian,full_matrices=True)
    rank=int(np.sum(singular>1e-8))
    nullspace=vh[rank:].T
    route_audits=[]
    for index,route in enumerate(routes_used):
        analytic=-data.ten_J[index].copy()
        projected=analytic@nullspace
        route_audits.append(dict(name=route['name'],length_m=float(data.ten_length[index]),
                                 source_neutral_length_m=route['source_neutral_length_m'],
                                 source_max_isometric_force_N=route['source_max_isometric_force_N'],
                                 effective_moment_norm_m=float(np.linalg.norm(projected)),
                                 mechanical_status='no joint leverage' if np.linalg.norm(projected)<1e-10 else 'has joint leverage; sign/range audit required'))
    output=Path(output)
    output.with_suffix('.audit.json').write_text(json.dumps(dict(model=str(output),routes=route_audits,
        actuator_count=model.nu,tendon_count=model.ntendon,no_joint_leverage=sum(x['effective_moment_norm_m']<1e-10 for x in route_audits),
        total_mass_kg=float(model.body_mass.sum()),loop_rank=rank,
        dynamic_path_policy=configuration['dynamic_path_policy'],force_policy=configuration['force_policy']),indent=2))
    print(json.dumps(dict(model=str(output),routes=model.nu,total_mass_kg=float(model.body_mass.sum()),
                         no_joint_leverage=sum(x['effective_moment_norm_m']<1e-10 for x in route_audits)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--routes',type=Path,default=HERE/'robotbody_routes.json')
    parser.add_argument('--output',type=Path,default=HERE/'ah001_neural.xml')
    parser.add_argument('--include-reference-routes',action='store_true',help='Include the other 38 passive source routes for complete audits')
    args=parser.parse_args()
    build(args.routes,args.output,args.include_reference_routes)
