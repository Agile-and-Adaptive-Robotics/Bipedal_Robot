"""Create an editable transfer specification; never overwrite user adjustments.

Transfers right-side robotbody routes and their explicit left mirrors. Runtime
moving and conditional definitions are carried into the editable specification.
"""
import argparse
import json
from pathlib import Path

import mujoco
import numpy as np

from build_ah001 import Q,cad_to_isb

HERE=Path(__file__).resolve().parent
AUDIT=HERE.parent/'spinal/reports_20261027_robot_control'


def basis(proximal,distal):
    up=np.array(proximal)-np.array(distal)
    up/=np.linalg.norm(up)
    forward=np.array([1.,0,0])-up*up[0]
    forward/=np.linalg.norm(forward)
    return np.column_stack((forward,up,np.cross(forward,up)))


def alignment(source_proximal,source_distal,target_proximal,target_distal,transverse=1.):
    sp,sd,tp,td=map(np.array,(source_proximal,source_distal,target_proximal,target_distal))
    length_scale=float(np.linalg.norm(tp-td)/np.linalg.norm(sp-sd))
    matrix=basis(tp,td)@np.diag([transverse,length_scale,transverse])@basis(sp,sd).T
    return dict(matrix=matrix.tolist(),offset=(tp-matrix@sp).tolist(),
                source_proximal_m=sp.tolist(),source_distal_m=sd.tolist(),
                target_proximal_m=tp.tolist(),target_distal_m=td.tolist(),
                longitudinal_scale=length_scale,transverse_scale=transverse)


def apply(mapping,point):
    return np.array(mapping['matrix'])@point+mapping['offset']


def distal_shank_center(side,cad):
    # Use the tibia's own coordinate frame to find its distal end from the
    # exported CAD mesh, then take the center of that end's transverse bounds.
    component=next(x for x in cad['components'] if ('05_01_TI_R' if side=='r' else '05_02_TI') in x['name'])
    transform=np.array(component['transform'])
    rotation=transform[:9].reshape(3,3).T
    path=HERE/'ah001_cad_meshes/meshes'/('shank_'+side+'.STL')
    raw=path.read_bytes()
    count=int(np.frombuffer(raw,dtype='<u4',count=1,offset=80)[0])
    faces=np.frombuffer(raw,dtype=np.dtype([('normal','<f4',(3,)),('vertices','<f4',(3,3)),('attr','<u2')]),count=count,offset=84)
    local=(faces['vertices'].reshape(-1,3).astype(float)-transform[9:12])@rotation
    end=local[local[:,1]<local[:,1].min()+.0005]
    center=(end.min(axis=0)+end.max(axis=0))/2
    center[1]=local[:,1].min()
    return cad_to_isb(rotation@center+transform[9:12])


def main(destination):
    destination=Path(destination)
    if destination.exists():
        raise FileExistsError('Preserve existing adjustments; choose a new --output file.')
    source=json.loads((HERE/'robotbody_source_inventory.json').read_text())
    cad=json.loads((AUDIT/'robot_mass_geometry.json').read_text())
    audit=json.loads((AUDIT/'robot_geometry_audit.json').read_text())
    spec=json.loads((HERE/'torso_spec.json').read_text())
    model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_with_torso.xml'))
    data=mujoco.MjData(model)
    mujoco.mj_forward(model,data)
    pelvis_id=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'pelvis')
    pelvis_world=data.xpos[pelvis_id].copy()
    body_origins={mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_BODY,i):data.xpos[i]-pelvis_world for i in range(1,model.nbody)}
    source_joint=lambda name: Q@np.array(source['joints'][name]['child_frame_origin_ground_m'])
    target_hip={side:body_origins['thigh_'+side] for side in ('r','l')}
    source_hip={side:source_joint('hip_'+side) for side in ('r','l')}
    hinge=np.array(spec['lumbar_hinge_position_in_pelvis_m'])
    hip_width=np.linalg.norm(target_hip['l']-target_hip['r'])/np.linalg.norm(source_hip['l']-source_hip['r'])
    mappings={'pelvis':alignment(source_joint('back'),(source_hip['l']+source_hip['r'])/2,
                                 hinge,(target_hip['l']+target_hip['r'])/2,hip_width)}
    target_names={'pelvis':'pelvis','torso':'lumbar_torso'}
    foot_mappings={}
    for side in ('r','l'):
        joint_anchors=[cad_to_isb(x['anchor_a']) for x in audit['joints']
                       if x['kind']!='SphericalJoint' and ('_R_' if side=='r' else '_L_') in (x['a']+x['b'])]
        assert len(joint_anchors)==4
        target_knee=np.mean(joint_anchors,axis=0)
        target_ankle=distal_shank_center(side,cad)
        mappings['femur_'+side]=alignment(source_hip[side],source_joint('knee_'+side),target_hip[side],target_knee)
        mappings['tibia_'+side]=alignment(source_joint('knee_'+side),source_joint('ankle_'+side),target_knee,target_ankle)
        target_names['femur_'+side]='thigh_'+side
        target_names['tibia_'+side]='shank_'+side
        foot_mappings[side]=dict(source_ankle_m=source_joint('ankle_'+side).tolist(),
                                 target_distal_center_m=target_ankle.tolist(),
                                 target_basis=basis(target_knee,target_ankle).tolist(),
                                 source_basis=basis(source_joint('knee_'+side),source_joint('ankle_'+side)).tolist(),
                                 transverse_scale=.05,radial_limit_m=.015)
        for foot in ('talus','calcn','toes'):
            target_names[foot+'_'+side]='foot_'+side
    torso_points=[np.array(p['location_ground_m']) for m in source['muscles'] for p in m['points'] if p['body']=='torso']
    source_hinge=source_joint('back')
    offsets=np.array([Q@p-source_hinge for p in torso_points])
    transverse_scale=(spec['spine_square_width_m']/2)/np.max(np.abs(offsets[:,[0,2]]))
    vertical_scale=spec['spine_length_m']/np.max(offsets[:,1])
    torso_matrix=np.diag([transverse_scale,vertical_scale,transverse_scale])
    mappings['torso']=dict(matrix=torso_matrix.tolist(),offset=(hinge-torso_matrix@source_hinge).tolist(),
                           transverse_scale=transverse_scale,longitudinal_scale=vertical_scale,
                           policy='Compress torso guides onto the square spine bar; preserve relative vertical order')
    selected={}
    names={m['name'] for m in source['muscles']}
    for actuator in source['actuators_27']:
        seed=next(n for n in actuator['robot_seed_preference'] if n in names)
        for side in ('r','l'):
            selected[seed[:-1]+side]=dict(actuator_id_per_side=actuator['id'],master_name=actuator['master_name'],
                                        reference_capacity_N=actuator['master_mif_n'],primary_dof=actuator['primary_dof'].replace('_r','_'+side))
    routes=[]
    source_length_errors=[]
    for muscle in source['muscles']:
        route=dict(name=muscle['name'],source_max_isometric_force_N=muscle['max_isometric_force_N'],
                   source_optimal_fiber_length_m=muscle['optimal_fiber_length_m'],
                   source_tendon_slack_length_m=muscle['tendon_slack_length_m'],
                   source_neutral_length_m=muscle['source_neutral_path_length_m'],
                   force_model='unassigned',points=[],source_neutral_moment_arms_m=muscle['source_neutral_moment_arms_m'])
        route['actuator_assignment']=selected.get(muscle['name'])
        enabled_source_points=[]
        for index,point in enumerate(muscle['points']):
            body=point['body']
            original=Q@np.array(point['location_ground_m'])
            if body.startswith(('talus_','calcn_','toes_')):
                mapping=foot_mappings[body[-1]]
                target=np.array(mapping['target_distal_center_m'])+(original-mapping['source_ankle_m'])
                policy='source foot point transferred to supplied CAD foot; ankle and foot joints welded'
            else:
                target=apply(mappings[body],original)
                policy='body landmark affine transfer'
            target_body=target_names[body]
            route['points'].append(dict(id='%s_%02d' % (muscle['name'],index),source_point_index=index,
                source_point_name=point['name'],source_point_type=point['type'],source_body=body,
                body=target_body,position_m=(target-body_origins[target_body]).tolist(),
                enabled=point['neutral_enabled'],transfer_policy=policy,
                source_motion_status='static' if point['type']=='PathPoint' else 'runtime source function or conditional membership',
                motion_functions=point.get('motion_functions',[]),
                conditional_coordinate=(point['condition_coordinate'] or '').split('/')[-1],
                conditional_range_rad=[float(x) for x in (point['condition_range'] or '').split()],
                source_location_parent_m=point['location_parent_m'],
                motion_to_target_matrix=(np.array(mappings.get(body,{'matrix':np.eye(3)})['matrix'])@np.array(source['bodies'][body]['rotation_to_ground'])).tolist()))
            if point['neutral_enabled']:
                enabled_source_points.append(original)
        reconstructed=sum(np.linalg.norm(b-a) for a,b in zip(enabled_source_points,enabled_source_points[1:]))
        source_length_errors.append(abs(reconstructed-muscle['source_neutral_path_length_m']))
        active_bodies=sorted({p['body'] for p in route['points'] if p['enabled']})
        route['crossed_bodies']=active_bodies
        route['mechanical_status']='no joint leverage: all points on one rigid body' if len(active_bodies)==1 else 'requires moment-arm audit'
        routes.append(route)
    if max(source_length_errors)>1e-8:
        raise AssertionError('Source neutral path reconstruction mismatch: '+str(max(source_length_errors)))
    result=dict(source=source['source'],source_sha256=source['sha256'],source_neutral_coordinate_values=source['coordinates'],
                coordinate_frame='target body-local metres; editable positions are authoritative in build_opensim_routes.py',
                source_neutral_length_reconstruction_max_error_m=max(source_length_errors),
                body_mappings=mappings,foot_mappings=foot_mappings,
                dynamic_path_policy='Exact source cubic moving-point functions and conditional membership updated at runtime; retargeted geometry remains an engineering adaptation',
                force_policy='zero-gain interfaces only; source human Fmax values are provenance, not robot actuator capacity',
                routes=routes)
    destination.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(routes=len(routes),points=sum(len(r['points']) for r in routes),
        actuators=sum(r['actuator_assignment'] is not None for r in routes),
        no_joint_leverage=sum(len(r['crossed_bodies'])==1 for r in routes),
        source_length_reconstruction_max_error_m=max(source_length_errors),output=str(destination)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'robotbody_routes.json')
    main(parser.parse_args().output)
