"""Measure uncrossed-leg geometry and actuator changes in recorded physics."""
import argparse,json
from run_neural_robot import *

def assess(path):
    trace=np.load(path);summary=json.loads(path.with_suffix('.json').read_text())
    model=mujoco.MjModel.from_xml_path(summary.get('model_path',str(HERE/'ah001_neural.xml')));data=mujoco.MjData(model)
    pelvis=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'pelvis')
    feet=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'foot_'+s) for s in ('r','l')]
    start=min(15.,max(2.,trace['time'][-1]-5.))
    heading=[];position=[];samples=[]
    for index in range(0,len(trace['time']),5):
        data.qpos[:]=trace['qpos'][index];mujoco.mj_kinematics(model,data)
        R=data.xmat[pelvis].reshape(3,3).T
        heading.append([np.rad2deg(np.arctan2((R@data.xmat[f].reshape(3,3)[:,0])[2],(R@data.xmat[f].reshape(3,3)[:,0])[0])) for f in feet])
        position.append([(R@(data.xpos[f]-data.xpos[pelvis]))[2] for f in feet]);samples.append(trace['time'][index])
    position=np.array(position);heading=np.array(heading);samples=np.array(samples)
    window=samples>=start;spacing=position[:,0]-position[:,1]
    activation={}
    for name in ('glut_max3_r','glut_max3_l','tfl_r','tfl_l','rect_fem_r','sar_r'):
        aid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_ACTUATOR,name)
        activation[name]=float(trace['activation'][trace['time']>=start,aid].mean())
    result=dict(run=path.stem,failed=summary['failed'],duration_s=summary['duration_simulated_s'],window_start_s=start,
        reversed_ankle_order_fraction=float(np.mean(spacing[window]<0)),
        full_run_reversed_ankle_order_fraction=float(np.mean(spacing<0)),
        minimum_ankle_spacing_m=float(spacing[window].min()),median_ankle_spacing_m=float(np.median(spacing[window])),
        median_heading_deg=np.median(heading[window],axis=0).tolist(),maximum_abs_heading_deg=np.max(abs(heading[window]),axis=0).tolist(),
        mean_abductor_activation=activation,
        mean_support_fraction=float(trace['rig_force'][trace['time']>=start,1].mean()/(75*9.81)))
    path.with_name(path.stem+'_abduction_assessment.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('runs',type=Path,nargs='+')
    for path in p.parse_args().runs:assess(path)
