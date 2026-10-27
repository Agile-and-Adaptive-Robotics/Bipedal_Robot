"""Verify path motion, conditional switching, full Jacobians, and source parity."""
import json,os,sys
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import mujoco
import numpy as np
import opensim as osim
from dynamic_routes import DynamicRoutes,spline_value_derivative
from robot_kinematics import pose

HERE=Path(__file__).resolve().parent
def main():
    config=json.loads((HERE/'robotbody_routes.json').read_text())
    source=json.loads((HERE/'robotbody_source_inventory.json').read_text())
    native=osim.Model(str(HERE.parents[2]/source['source']));native.setUseVisualizer(False);state=native.initSystem()
    coordinates=native.getCoordinateSet()
    for i in range(coordinates.getSize()):coordinates.get(i).setValue(state,0.,False)
    native.realizePosition(state)
    maximum_function_error=0.;condition_checks=0;function_checks=0
    for muscle in source['muscles']:
        if not muscle['name'].endswith('_r'):continue
        native_points=native.getMuscles().get(muscle['name']).getGeometryPath().getPathPointSet()
        for index,point in enumerate(muscle['points']):
            if point['type']=='MovingPathPoint':
                native_point=osim.MovingPathPoint.safeDownCast(native_points.get(index))
                for axis,f in zip('xyz',point['motion_functions']):
                    function=getattr(native_point,'get_'+axis+'_location')()
                    x=osim.Vector(1,0.)
                    for value in np.linspace(-1.8,.12,41):
                        x.set(0,float(value));reference=float(function.calcValue(x))
                        actual,_=spline_value_derivative(value,np.array(f['knots']),np.array(f['coefficients']),np.array(f['endpoints']))
                        maximum_function_error=max(maximum_function_error,abs(reference-actual));function_checks+=1
            if point['type']=='ConditionalPathPoint':
                native_point=osim.ConditionalPathPoint.safeDownCast(native_points.get(index))
                name=point['condition_coordinate'].split('/')[-1];lo,hi=map(float,point['condition_range'].split())
                for value in (lo-.001,lo+.001,(lo+hi)/2,hi-.001,hi+.001):
                    coordinates.get(name).setValue(state,float(value),False);native.realizePosition(state)
                    assert bool(native_point.isActive(state))==(lo<=value<=hi);condition_checks+=1
                coordinates.get(name).setValue(state,0.,False)
    assert maximum_function_error<1e-10
    model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_neural.xml'));dynamic=DynamicRoutes(model,config)
    results={};maximum_derivative_error=0.;switch_counts={};boundary_events=[]
    for kind,side,limit in [('hip','r',30),('hip','l',30),('knee','r',35),('knee','l',35),('lumbar','r',20)]:
        label=kind if kind=='lumbar' else kind+'_'+side;samples=[]
        for degrees in np.linspace(-limit,limit,15):
            angle=np.deg2rad(degrees)
            try:
                middle=pose(model,dynamic,kind,side,angle);jac=dynamic.length_jacobian.copy();active=dynamic.active_counts.copy()
                plus=pose(model,dynamic,kind,side,angle+1e-6);length_plus=plus.ten_length.copy();active_plus=dynamic.active_counts.copy()
                minus=pose(model,dynamic,kind,side,angle-1e-6);length_minus=minus.ten_length.copy();active_minus=dynamic.active_counts.copy()
            except (ValueError,np.linalg.LinAlgError) as error:
                samples.append(dict(degrees=float(degrees),reachable=False,error=str(error)));continue
            velocity=np.zeros(model.nv);mujoco.mj_differentiatePos(model,velocity,2e-6,minus.qpos,plus.qpos)
            direct=-jac@velocity;finite=-(length_plus-length_minus)/(2e-6)
            stable=active_plus==active_minus
            for index in np.flatnonzero(~stable):
                boundary_events.append(dict(coordinate=label,degrees=float(degrees),route=config['routes'][index]['name'],
                    length_change_across_two_microradians_m=float(length_plus[index]-length_minus[index]),
                    selected_actuator=config['routes'][index]['actuator_assignment'] is not None))
            error=float(np.max(np.abs(direct[stable]-finite[stable])));maximum_derivative_error=max(maximum_derivative_error,error)
            assert error<2e-6,(label,degrees,error)
            rows=np.flatnonzero(middle.efc_type==int(mujoco.mjtConstraint.mjCNSTR_EQUALITY))
            assert np.max(np.abs(middle.efc_pos[rows]))<1e-9
            samples.append(dict(degrees=float(degrees),reachable=True,lengths_m=middle.ten_length.tolist(),
                                moment_arms_m=direct.tolist(),active_point_counts=active.tolist()))
        results[label]=samples
    reports=[]
    for index,route in enumerate(config['routes']):
        disagreements=[]
        for label,samples in results.items():
            zero=next(s for s in samples if s['reachable'] and abs(s['degrees'])<1e-10)
            source_coordinate='lumbar_extension' if label=='lumbar' else ('hip_flexion_' if label.startswith('hip') else 'knee_angle_')+label[-1]
            expected=route['source_neutral_moment_arms_m'].get(source_coordinate,0.)
            actual=zero['moment_arms_m'][index]
            if abs(expected)>.001 and abs(actual)>.001 and expected*actual<0:
                disagreements.append(dict(coordinate=label,expected=expected,actual=actual))
        counts={s['active_point_counts'][index] for samples in results.values() for s in samples if s['reachable']}
        if len(counts)>1:switch_counts[route['name']]=sorted(counts)
        reports.append(dict(name=route['name'],neutral_sign_disagreements=disagreements,
                            selected_actuator=route['actuator_assignment'] is not None))
    summary=dict(frame='x anterior, y superior, z right',actuators=model.nu,
                 source_function_checks=function_checks,source_conditional_checks=condition_checks,
                 source_function_maximum_error_m=maximum_function_error,
                 full_moment_arm_maximum_error_m=maximum_derivative_error,
                 reachable_poses=sum(s['reachable'] for ss in results.values() for s in ss),
                 conditional_switches=switch_counts,
                 conditional_boundary_events=boundary_events,
                 reversed_routes=[r['name'] for r in reports if r['neutral_sign_disagreements']])
    (HERE/'dynamic_route_verification.json').write_text(json.dumps(dict(summary=summary,routes=reports,samples=results),indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
