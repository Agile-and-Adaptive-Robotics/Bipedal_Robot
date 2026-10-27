"""Fit small mirrored attachment corrections to restore source neutral actions.

This changes engineering attachment geometry, never muscle strength or neural
signals. Every changed point and before/after moment arm is saved for review.
"""
import json
from pathlib import Path
import mujoco
import numpy as np
from scipy.optimize import minimize
from dynamic_routes import DynamicRoutes
from robot_kinematics import pose

HERE=Path(__file__).resolve().parent
def main():
    path=HERE/'robotbody_routes.json';configuration=json.loads(path.read_text())
    if configuration.get('neutral_action_corrections'):
        raise RuntimeError('Corrections already applied; inspect current verification rather than applying twice.')
    backup=HERE.parent/'spinal/reports_20261027_robot_control/robotbody_routes_before_sign_fixes.json'
    backup.write_text(json.dumps(configuration,indent=2))
    routes={r['name']:r for r in configuration['routes']}
    model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_neural.xml'));dynamic=DynamicRoutes(model,configuration)
    neutral=pose(model,dynamic,'hip','r',0.)
    tangents={}
    for kind in ('hip','knee'):
        for side in ('r','l'):
            plus=pose(model,dynamic,kind,side,1e-6);minus=pose(model,dynamic,kind,side,-1e-6)
            velocity=np.zeros(model.nv);mujoco.mj_differentiatePos(model,velocity,2e-6,minus.qpos,plus.qpos)
            tangents[kind+'_'+side]=velocity
    corrections=[]
    for base,indexes in [('add_mag2',[0,1]),('grac',[0,1,2,3]),('semimem',[0,1])]:
        constraints=[]
        for side in ('r','l'):
            route=routes[base+'_'+side]
            tid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_TENDON,route['name'])
            for kind in ('hip','knee'):
                name=('hip_flexion_' if kind=='hip' else 'knee_angle_')+side
                expected=route['source_neutral_moment_arms_m'].get(name,0.)
                if abs(expected)>.001:
                    constraints.append((tid,kind+'_'+side,np.sign(expected),max(.0015,.2*abs(expected)),expected))
        originals={}
        for side in ('r','l'):
            for index in indexes:
                point=routes[base+'_'+side]['points'][index]
                assert point['source_point_type']=='PathPoint'
                sid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,point['id'])
                originals[sid]=model.site_pos[sid].copy()
        def apply(offset):
            for side in ('r','l'):
                for k,index in enumerate(indexes):
                    point=routes[base+'_'+side]['points'][index]
                    sid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,point['id'])
                    model.site_pos[sid]=originals[sid]+np.r_[offset[2*k:2*k+2],0.]
            dynamic.update(neutral)
        def margins(offset,full=False):
            apply(offset)
            actual=np.array([dynamic.moment_arms[tid]@tangents[label] for tid,label,*_ in constraints])
            return actual if full else np.array([value*sign-minimum for value,(_,_,sign,minimum,_) in zip(actual,constraints)])
        zero=np.zeros(2*len(indexes));before=margins(zero,True)
        result=minimize(lambda x:float(x@x),zero,method='SLSQP',bounds=[(-.06,.06)]*len(zero),
                        constraints=[dict(type='ineq',fun=margins)],options=dict(maxiter=150,ftol=1e-13))
        after=margins(result.x,True)
        if not result.success or margins(result.x).min() < -1e-7:
            raise AssertionError((base,result.message,margins(result.x)))
        points=[]
        for side in ('r','l'):
            for k,index in enumerate(indexes):
                point=routes[base+'_'+side]['points'][index]
                delta=np.r_[result.x[2*k:2*k+2],0.]
                point['position_m']=(np.array(point['position_m'])+delta).tolist()
                point['neutral_action_adjustment_m']=delta.tolist()
                points.append(dict(id=point['id'],body=point['body'],adjustment_m=delta.tolist()))
        corrections.append(dict(route_pair=base,method='minimum squared x/y attachment offsets; same correction on mirrored sides',
            maximum_point_displacement_m=float(max(np.linalg.norm(result.x[2*k:2*k+2]) for k in range(len(indexes)))),
            before_moment_arms_m=before.tolist(),after_moment_arms_m=after.tolist(),
            constraints=[dict(coordinate=label,expected_source_m=expected,required_signed_minimum_m=minimum) for _,label,_,minimum,expected in constraints],points=points))
    configuration['neutral_action_corrections']=corrections
    path.write_text(json.dumps(configuration,indent=2))
    (HERE/'route_sign_corrections.json').write_text(json.dumps(corrections,indent=2))
    print(json.dumps([{k:v for k,v in c.items() if k!='points'} for c in corrections],indent=2))
if __name__=='__main__':main()
