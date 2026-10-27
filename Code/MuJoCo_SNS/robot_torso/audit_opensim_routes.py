"""Check transferred route lengths and moment arms on the constrained robot."""
import json
from pathlib import Path

import mujoco
import numpy as np
from PIL import Image

HERE=Path(__file__).resolve().parent


def knee_pose(model,data,side,angle,seed=None):
    joints=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,x+'_'+side+'_joint') for x in ('back','shank','front')]
    addresses=model.jnt_qposadr[joints]
    signs=model.jnt_axis[joints,1]
    assert np.allclose(np.abs(signs),1.,atol=1e-9)
    eq=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_EQUALITY,'knee_loop_'+side)
    first,second=model.eq_obj1id[eq],model.eq_obj2id[eq]
    def residual(q):
        data.qpos[addresses]=q
        mujoco.mj_kinematics(model,data)
        a=data.xpos[first]+data.xmat[first].reshape(3,3)@model.eq_data[eq,:3]
        b=data.xpos[second]+data.xmat[second].reshape(3,3)@model.eq_data[eq,3:6]
        return np.r_[(a-b)[[0,2]],q[:2]@signs[:2]-angle]
    q=np.zeros(3) if seed is None else seed.copy()
    for _ in range(40):
        r=residual(q)
        if np.linalg.norm(r)<1e-11:
            residual(q)
            mujoco.mj_forward(model,data)
            return q
        h=1e-7
        jac=np.column_stack([(residual(q+np.eye(3)[i]*h)-residual(q-np.eye(3)[i]*h))/(2*h) for i in range(3)])
        update=np.linalg.solve(jac,-r)
        # Backtracking avoids jumping between four-bar assembly branches.
        for factor in (1.,.5,.25,.125,.0625,.03125):
            candidate=q+factor*update
            if np.linalg.norm(residual(candidate))<np.linalg.norm(r):
                q=candidate
                break
        else:
            raise ValueError('Knee branch could not reach angle %.6f' % angle)
    raise ValueError('Knee solve did not converge')


def pose(model,side,angle,kind):
    data=mujoco.MjData(model)
    if kind=='knee':
        knee_pose(model,data,side,angle)
    else:
        joint=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,
                              'thigh_'+side+'_joint' if kind=='hip' else 'lumbar_flexion')
        tangent=np.zeros(model.nv)
        tangent[model.jnt_dofadr[joint]+(1 if kind=='hip' else 0)]=1.
        mujoco.mj_integratePos(model,data.qpos,tangent,angle)
        mujoco.mj_forward(model,data)
    return data


def main():
    config=json.loads((HERE/'opensim_routes.json').read_text())
    model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_with_opensim_routes.xml'))
    neutral=pose(model,'r',0.,'lumbar')
    assert model.nu==len(config['routes'])==92
    source_length_error=config['source_neutral_length_reconstruction_max_error_m']
    assert source_length_error<1e-10
    # A control pulse must produce exactly zero active forces at these interfaces.
    neutral.ctrl[:]=1.
    mujoco.mj_forward(model,neutral)
    assert np.array_equal(neutral.actuator_force,np.zeros(model.nu))
    results={}
    largest_difference=0.
    for kind,side,limit in [('hip','r',30),('hip','l',30),('knee','r',15),('knee','l',15),('lumbar','r',20)]:
        label=kind if kind=='lumbar' else kind+'_'+side
        samples=[]
        for degrees in np.linspace(-limit,limit,13):
            theta=np.deg2rad(degrees)
            try:
                middle=pose(model,side,theta,kind)
                plus=pose(model,side,theta+1e-6,kind)
                minus=pose(model,side,theta-1e-6,kind)
            except (ValueError,np.linalg.LinAlgError) as error:
                samples.append(dict(degrees=float(degrees),reachable=False,error=str(error)))
                continue
            assert np.isfinite(middle.ten_length).all()
            velocity=np.zeros(model.nv)
            mujoco.mj_differentiatePos(model,velocity,2e-6,minus.qpos,plus.qpos)
            direct=-middle.ten_J@velocity
            finite=-(plus.ten_length-minus.ten_length)/(2e-6)
            error=float(np.max(np.abs(direct-finite)))
            largest_difference=max(largest_difference,error)
            assert error<2e-6
            rows=np.flatnonzero(middle.efc_type==int(mujoco.mjtConstraint.mjCNSTR_EQUALITY))
            assert np.max(np.abs(middle.efc_pos[rows]))<1e-9
            samples.append(dict(degrees=float(degrees),reachable=True,
                lengths_m=middle.ten_length.tolist(),moment_arms_m=finite.tolist(),
                maximum_derivative_difference_m=error))
        results[label]=samples
    route_reports=[]
    for index,route in enumerate(config['routes']):
        excursions={}
        sign_changes=[]
        neutral_sign_disagreements=[]
        for label,samples in results.items():
            valid=[s for s in samples if s['reachable']]
            arms=np.array([s['moment_arms_m'][index] for s in valid])
            lengths=np.array([s['lengths_m'][index] for s in valid])
            excursions[label]=dict(length_min_m=float(lengths.min()),length_max_m=float(lengths.max()),
                                    moment_arm_min_m=float(arms.min()),moment_arm_max_m=float(arms.max()))
            if arms.min() < -.001 and arms.max() > .001:
                sign_changes.append(label)
            source_coord='lumbar_extension' if label=='lumbar' else ('hip_flexion_' if label.startswith('hip') else 'knee_angle_')+label[-1]
            expected=-route['source_neutral_moment_arms_m'][source_coord]
            center=next(s for s in valid if abs(s['degrees'])<1e-12)['moment_arms_m'][index]
            if abs(expected)>.001 and abs(center)>.001 and expected*center<0:
                neutral_sign_disagreements.append(dict(coordinate=label,source_sign_mapped_m=expected,target_moment_arm_m=center))
        route_reports.append(dict(name=route['name'],excursions=excursions,
                                   sign_changes_over_sampled_ranges=sign_changes,
                                   neutral_sign_disagreements=neutral_sign_disagreements,
                                   joint_leverage=any(max(abs(x['moment_arm_min_m']),abs(x['moment_arm_max_m']))>1e-8 for x in excursions.values())))
    lumbar=next(s for s in results['lumbar'] if abs(s['degrees'])<1e-12)
    es_indices=[i for i,r in enumerate(config['routes']) if r['name'].startswith('ercspn_')]
    es_arm_sum=sum(lumbar['moment_arms_m'][i] for i in es_indices)
    joint=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'lumbar_flexion')
    holding=float(neutral.qfrc_bias[model.jnt_dofadr[joint]])
    equal_force=holding/es_arm_sum if abs(es_arm_sum)>1e-12 else None
    summary=dict(routes=92,source_path_length_reconstruction_max_error_m=source_length_error,
                 moment_arm_derivative_max_error_m=largest_difference,
                 sampled_configurations=sum(len(x) for x in results.values()),
                 reachable_configurations=sum(s['reachable'] for x in results.values() for s in x),
                 source_sign_disagreement_routes=[r['name'] for r in route_reports if r['neutral_sign_disagreements']],
                 routes_with_sign_changes=[r['name'] for r in route_reports if r['sign_changes_over_sampled_ranges']],
                 no_sagittal_leverage_routes=[r['name'] for r in route_reports if not r['joint_leverage']],
                 erector_pair_signed_moment_arm_sum_m=es_arm_sum,
                 equal_erector_force_each_for_neutral_hold_N=equal_force,
                 total_mass_kg=float(model.body_mass.sum()),all_interface_forces_zero=True,
                 status='Neutral-path engineering transfer; moving points and conditional membership frozen; actuator capacity and neural control unvalidated')
    (HERE/'opensim_route_motion_audit.json').write_text(json.dumps(dict(summary=summary,routes=route_reports,samples=results),indent=2))
    model.tendon_rgba[:,3]=.1
    for index,route in enumerate(config['routes']):
        if route['name'].startswith(('vas_lat_','semimem_','iliopsoas_','glut_max1_','ercspn_','intobl_','extobl_')):
            model.tendon_rgba[index]=[.9,.25,.15,.85]
            model.tendon_width[index]=.002
    camera=mujoco.MjvCamera()
    mujoco.mjv_defaultCamera(camera)
    camera.lookat[:]=[-.06,0,.64]
    camera.distance=2.05
    camera.azimuth=140
    camera.elevation=-12
    renderer=mujoco.Renderer(model,height=1000,width=1200)
    renderer.update_scene(neutral,camera)
    Image.fromarray(renderer.render()).save(HERE/'ah001_opensim_routes_preview.png')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
