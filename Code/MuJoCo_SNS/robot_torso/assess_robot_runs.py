import json
from pathlib import Path
import numpy as np
from run_neural_robot import rotation_error
HERE=Path(__file__).resolve().parent
def main():
    reference=np.load(HERE/'runs/stand_20s_reference.npz');fast=np.load(HERE/'runs/stand_20s_exact.npz')
    differences={key:float(np.max(np.abs(reference[key]-fast[key]))) for key in reference.files}
    window=fast['time']>=15
    frame=fast['qpos'][window]
    # Root qpos: translation, free quaternion, then hips/knee mechanism and ball.
    import mujoco
    model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_neural.xml'))
    joint=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'lumbar_ball');address=model.jnt_qposadr[joint]
    tilt=np.rad2deg([np.linalg.norm(rotation_error(q[3:7])) for q in frame])
    lumbar=np.rad2deg([np.linalg.norm(rotation_error(q[address:address+4])) for q in frame])
    metrics=dict(final_five_seconds_knee_excursion_deg=np.rad2deg(np.ptp(fast['knee_coordinate'][window],axis=0)).tolist(),
        final_five_seconds_pelvis_position_range_m=np.ptp(frame[:,:3],axis=0).tolist(),
        maximum_pelvis_tilt_deg=float(max(tilt)),maximum_lumbar_angle_deg=float(max(lumbar)),
        mean_ground_loads_N=np.mean(fast['contact_loads'][window].sum(axis=2),axis=0).tolist(),
        mean_vertical_support_fraction=float(np.mean(fast['rig_force'][window,1])/(75*9.81)),
        mean_rig_torque_Nm=np.mean(np.abs(fast['rig_torque'][window]),axis=0).tolist())
    result=dict(full_twenty_second_backend_parity=all(x==0 for x in differences.values()),maximum_differences=differences,
                steady_standing_metrics=metrics)
    if (HERE/'runs/step_20s.npz').exists():
        step=np.load(HERE/'runs/step_20s.npz');walking=step['time']>=5
        rg=step['rg'][walking]
        differences_rg=np.column_stack((rg[:,1]-rg[:,0],rg[:,3]-rg[:,2]))
        cycles=[int(np.sum((x[1:]>0)&(x[:-1]<=0))) for x in differences_rg.T]
        ground=step['contact_loads'][walking].sum(axis=2)
        contacts=[int(np.sum((x[1:]>10)&(x[:-1]<=10))) for x in ground.T]
        result['stepping_attempt']=dict(rg_difference_excursion_mV=np.ptp(differences_rg,axis=0).tolist(),
            rg_flexor_crossings=cycles,raw_contact_threshold_crossings=contacts,
            knee_excursions_after_drive_deg=np.rad2deg(np.ptp(step['knee_coordinate'][walking],axis=0)).tolist(),
            displacement_after_drive_m=float(step['qpos'][walking][-1,0]-step['qpos'][walking][0,0]),
            status='Contact threshold crossings alone do not establish gait; phase, lifting, and progression require review')
    (HERE/'controller_assessment.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    config=json.loads((HERE/'robotbody_routes.json').read_text())
    audit=json.loads((HERE/'dynamic_route_verification.json').read_text())
    neutral=next(s for s in audit['samples']['knee_r'] if abs(s['degrees'])<1e-10)
    for i,route in enumerate(config['routes']):
        if route['actuator_assignment'] and route['name'].endswith('_r'):
            print(route['name'], 'knee arm mm',round(1000*neutral['moment_arms_m'][i],3))
if __name__=='__main__':main()
