"""Evaluate assisted robot stepping without counting contact chatter as steps."""
import argparse,json
from pathlib import Path
from run_neural_robot import rotation_error
import numpy as np
import mujoco
from scipy.signal import find_peaks

HERE=Path(__file__).resolve().parent

def bouts(mask,time,minimum=.06):
    edges=np.diff(np.r_[False,mask,False].astype(int))
    return [(float(time[a]),float(time[b-1]+.01)) for a,b in zip(np.flatnonzero(edges==1),np.flatnonzero(edges==-1)) if time[b-1]+.01-time[a]>=minimum]

def assess(path):
    trace=np.load(path);summary=json.loads(path.with_suffix('.json').read_text())
    options=summary['run_options']
    start=options['stand_time']+2
    if options.get('support_cap_final') is not None:
        start=max(start,options['support_wean_start']+options['support_wean_duration'])
    keep=trace['time']>=start
    if keep.sum()<20:return {'run':path.stem,'failed':summary['failed'],'assessment':'Failed before stepping assessment window'}
    t=trace['time'][keep];q=trace['qpos'][keep]
    rg=trace['rg'][keep];difference=rg[:,[1,3]]-rg[:,[0,2]]
    peaks=[find_peaks(x,prominence=3,distance=30)[0] for x in difference.T]
    correlation=float(np.corrcoef(difference.T)[0,1])
    loads=trace['contact_loads'][keep].sum(axis=2)
    swings=[bouts(x<5,t) for x in loads.T]
    contacts=[bouts(x>10,t) for x in loads.T]
    # Track actual foot geometry at 20 Hz; no guessed ankle-height proxy.
    model=mujoco.MjModel.from_xml_path(summary.get('model_path',str(HERE/'ah001_neural.xml')));data=mujoco.MjData(model)
    foot_geoms=[]
    for side in ('r','l'):
        foot_geoms.append([i for i in range(model.ngeom) if (mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_GEOM,i) or '').startswith('foot_'+side+'_')])
    clearance=[];ankle_lateral_order=[]
    pelvis=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'pelvis')
    feet=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'foot_'+s) for s in ('r','l')]
    for row in q[::5]:
        data.qpos[:]=row;mujoco.mj_forward(model,data)
        pelvis_rotation=data.xmat[pelvis].reshape(3,3).T
        ankle_lateral_order.append(float((pelvis_rotation@(data.xpos[feet[0]]-data.xpos[feet[1]]))[2]))
        heights=[]
        for ids in foot_geoms:
            minimum=np.inf
            for gid in ids:
                mid=model.geom_dataid[gid]
                if mid<0:continue
                vertices=model.mesh_vert[model.mesh_vertadr[mid]:model.mesh_vertadr[mid]+model.mesh_vertnum[mid]]
                height=vertices@data.geom_xmat[gid].reshape(3,3)[1]+data.geom_xpos[gid,1]
                minimum=min(minimum,float(height.min()))
            heights.append(minimum)
        clearance.append(heights)
    clearance=np.array(clearance)
    clearance_time=t[::5]
    apex=[find_peaks(clearance[:,side],height=.005,prominence=.004,distance=12)[0] for side in range(2)]
    apex_times=[clearance_time[index] for index in apex]
    lr_phases=[]
    for a,b in zip(apex_times[0][:-1],apex_times[0][1:]):
        inside=apex_times[1][(apex_times[1]>a)&(apex_times[1]<b)]
        if len(inside)==1:lr_phases.append(float((inside[0]-a)/(b-a)))
    lifted_bouts=[]
    for side in range(2):
        lifted_bouts.append([(a,b) for a,b in swings[side]
            if np.any((clearance_time>=a)&(clearance_time<b)) and
            float(clearance[(clearance_time>=a)&(clearance_time<b),side].max())>.005])
    lumbar=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'lumbar_ball');a=model.jnt_qposadr[lumbar]
    displacement=float(q[-1,0]-q[0,0]);duration=float(t[-1]-t[0])
    result=dict(run=path.stem,failed=summary['failed'],window_start_s=start,window_duration_s=duration,
        rg_bursts=[len(x) for x in peaks],rg_difference_excursion_mV=np.ptp(difference,axis=0).tolist(),rg_lr_correlation=correlation,
        sustained_swing_bouts=swings,sustained_stance_bouts=contacts,lifted_swing_bouts=lifted_bouts,
        foot_apex_times_s=[x.tolist() for x in apex_times],left_apex_phases_within_right_cycles=lr_phases,
        peak_foot_clearance_m=clearance.max(axis=0).tolist(),
        knee_excursion_deg=np.rad2deg(np.ptp(trace['knee_coordinate'][keep],axis=0)).tolist(),
        forward_displacement_m=displacement,mean_forward_speed_m_s=displacement/duration,
        actual_support_fraction=float(trace['rig_force'][keep,1].mean()/(75*9.81)),
        maximum_actual_support_fraction=float(trace['rig_force'][keep,1].max()/(75*9.81)),
        mean_absolute_horizontal_rig_force_N=np.mean(np.abs(trace['rig_force'][keep][:,[0,2]]),axis=0).tolist(),
        mean_absolute_rig_torque_Nm=np.mean(np.abs(trace['rig_torque'][keep]),axis=0).tolist(),
        minimum_recorded_activation=float(trace['activation'][keep].min()),
        minimum_right_minus_left_ankle_spacing_m=float(min(ankle_lateral_order)),
        reversed_ankle_order_fraction=float(np.mean(np.array(ankle_lateral_order)<0)),
        maximum_pelvis_tilt_deg=float(max(np.rad2deg(np.linalg.norm(rotation_error(row[3:7]))) for row in q)),
        maximum_lumbar_rotation_deg=float(max(np.rad2deg(np.linalg.norm(rotation_error(row[a:a+4]))) for row in q)))
    # These are deliberately separate necessary gates, not an all-purpose gait score.
    result['gates']=dict(completed=summary['failed'] is None,bilateral_rhythm=min(result['rg_bursts'])>=3,
        antiphase=correlation<-.5,bilateral_sustained_swings=min(map(len,swings))>=3,
        bilateral_foot_lifting=min(map(len,lifted_bouts))>=3,forward_progression=displacement>.1)
    coverage=len(lr_phases)/max(1,len(apex_times[0])-1)
    result['matched_foot_cycle_fraction']=coverage
    result['gates']['repeated_alternating_foot_apices']=len(lr_phases)>=3 and coverage>=.8 and all(.2<phase<.8 for phase in lr_phases)
    result['gates']['uncrossed_ankle_order']=min(ankle_lateral_order)>=0
    stepping=all(value for key,value in result['gates'].items() if key!='forward_progression')
    result['assessment']='Supported stepping candidate; forward walking not established' if stepping else 'Supported motion; stepping gates incomplete'
    path.with_name(path.stem+'_assessment.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('runs',type=Path,nargs='+')
    for run in parser.parse_args().runs:assess(run)
