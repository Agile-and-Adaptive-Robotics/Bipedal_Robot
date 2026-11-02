"""Bounded backend-specific Hunt recruitment fit; no production mutation.

The target is the accepted native-surrogate air kinematics, not imposed motion.
Every candidate is an ordinary live forward simulation, passing both original
air gates before ground promotion. Physical lengths and afferent definitions
are fixed: the controller's existing muscle-specific midpoint/half-range is
derived from routed lengths along the source IK trajectory (prepare_targets).
Changing that normalization would invalidate its separately fitted Hunt B fit.
"""
import os,sys,json,subprocess,hashlib,argparse
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
HERE=Path(__file__).resolve().parent;S=HERE.parent
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\physics')
BASE=ROOT/'calibration_air_baselines';MYO=r'D:\Anaconda\envs\myo\python.exe'
SHARED=S/'hunt_training/runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/zero_hyperextension_air_v1/shared.json'
_JOINT_ADDRESSES=None

def joint_addresses():
    global _JOINT_ADDRESSES
    if _JOINT_ADDRESSES is None:
        import mujoco
        sys.path.insert(0,str(S));import runner as R
        m=mujoco.MjModel.from_xml_path(str(R.MODEL))
        names=[j+'_'+side for side in ('r','l') for j in ('hip_flexion','knee_angle','ankle_angle')]
        _JOINT_ADDRESSES=(names,[int(m.joint(name).qposadr[0]) for name in names])
    return _JOINT_ADDRESSES

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def load_result(folder):
    b=json.loads((folder/'BACKEND_RESULT.json').read_text());r=json.loads((folder/'controller.result.json').read_text())
    return b,r

def run_stage(folder,stage,parameters,physics='compliant'):
    if (folder/'BACKEND_RESULT.json').exists():return load_result(folder)
    folder.mkdir(parents=True,exist_ok=True);cfg=json.loads(SHARED.read_text());cfg['diagnose_swing']=False
    cfg['motor_scales']={group+'_'+side:float(parameters[group]) for group in ('hip_flex','knee_flex','ankle_df') for side in ('r','l')}
    target=folder/'INPUT.json';target.write_text(json.dumps(cfg,indent=2))
    with (folder.parent/'simulation_output.log').open('a',encoding='utf-8') as log:
        p=subprocess.run([MYO,str(HERE/'live_tendon_trial.py'),'--protocol','hunt-air'+str(stage),'--physics',physics,'--config',str(target),'--out',str(folder)],stdout=log,stderr=subprocess.STDOUT,timeout=180,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    if not (folder/'BACKEND_RESULT.json').exists():raise RuntimeError('Missing simulation result, exit='+str(p.returncode))
    return load_result(folder)

def errors(folder,stage):
    # Backend qpos is actual physical state in radians, not runner's degree q.
    ref=np.load(BASE/('native_s'+str(stage))/'backend_trace.npz')
    z=np.load(folder/'backend_trace.npz');names,ids=joint_addresses()
    rt=np.asarray(ref['t']);t=np.asarray(z['t']);keep=(t>=5.)&(t<=13.);t=t[keep];q=z['qpos'][keep][:,ids];rq=ref['qpos'][:,ids]
    scale=np.maximum(np.ptp(rq[(rt>=5)&(rt<=13)],axis=0),np.radians(10.))
    # One common small time shift; never shift the legs independently.
    best=None
    for shift in np.arange(-.2,.2001,.01):
        target=np.column_stack([np.interp(t+shift,rt,rq[:,i]) for i in range(6)])
        row=(float(np.mean(((q-target)/scale)**2)),float(shift),np.sqrt(np.mean((q-target)**2,axis=0)))
        if best is None or row[0]<best[0]:best=row
    return dict(normalized_mse=best[0],common_time_shift_s=best[1],joint_rms_deg=np.degrees(best[2]).tolist(),joint_order=names)

def handle_request(request):
    manifest=request['manifest'];folder=Path(manifest['output_dir'])
    if manifest.get('calibration_kind')=='source_loaded_ground':return handle_loaded_request(request)
    if request['op']=='preflight':
        for stage in (1,2):
            for mode in ('native','compliant'):
                b,r=load_result(BASE/(mode+'_s'+str(stage)))
                if b['guard'] is not None or not r['assessment']['accepted']:raise RuntimeError('Unchanged baseline air gate failed: '+mode+str(stage))
        return dict(ok=True,ready_for_optimization=True,scope='Shared bilateral PF-to-MN group recruitment only; fixed source lengths, physical laws, timing, reflexes, initializer, pose, and gates.',source_config_sha256=sha(SHARED),driver_sha256=sha(HERE/'live_tendon_trial.py'),force_backend_sha256=sha(HERE/'thelen_dynamic.py'))
    if request['op']!='evaluate':raise ValueError('Unknown request')
    p=request['params'];trial=folder/request['evaluation_id'];rows=[];failed=False
    for stage in (1,2):
        # Reuse already measured identity baseline without rerunning.
        case=BASE/('compliant_s'+str(stage)) if all(float(p[k])==1. for k in p) else trial/('stage'+str(stage))
        b,r=load_result(case) if (case/'BACKEND_RESULT.json').exists() else run_stage(case,stage,p)
        accepted=b['guard'] is None and r['assessment']['accepted'];row=dict(stage=stage,accepted=accepted,path=str(case),assessment=r['assessment'],guard=b['guard']);rows.append(row)
        if not accepted:failed=True;break
        row['tracking']=errors(case,stage)
    value=1000.+100.*(2-len(rows)) if failed else float(np.mean([r['tracking']['normalized_mse'] for r in rows]))+.0001*sum((float(v)-1.)**2 for v in p.values())
    return dict(ok=True,value=value,diagnostics=dict(air_pair_accepted=not failed,stages=rows,parameters_tied=True,ground_run=False))

def handle_loaded_request(request):
    """Bounded stance-recruitment correction after observed loaded collapse."""
    root=ROOT/'hunt_balance_v1';manifest=request['manifest'];folder=Path(manifest['output_dir'])
    bases={1:root/'pt_calf_air1',2:root/'pt_calf_air2',3:root/'pt_calf_stand_vertical_free',4:root/'force_support_0p5'}
    if request['op']=='preflight':
        for stage in (1,2,3):
            b,r=load_result(bases[stage])
            if b['guard'] is not None or not r['assessment']['accepted']:raise RuntimeError('Prior stages must pass before load calibration')
        baseline=json.loads((bases[4]/'BACKEND_RESULT.json').read_text())
        return dict(ok=True,ready_for_optimization=True,baseline_loaded_guard=baseline['guard'],scope='Only shared knee_ext and ankle_pf PF-to-MN scales; fixed source dimensions, POST0.5, PT-only, all original gates, actual force-controlled0.5BW.')
    if request['op']!='evaluate':raise ValueError('Unknown loaded-calibration operation')
    p=request['params']
    if set(p)!={'knee_ext','ankle_pf'} or any(not 1.<=float(v)<=2.5 for v in p.values()):raise ValueError('Declared bilateral stance gains only')
    rows=[];accepted=False;value=2000.
    for stage in (1,2,3,4):
        case=bases[stage] if all(float(v)==1. for v in p.values()) else folder/request['evaluation_id']/('stage'+str(stage))
        if not (case/'BACKEND_RESULT.json').exists():
            case.mkdir(parents=True,exist_ok=True);cfg=json.loads((bases[stage]/'controller.json').read_text());cfg['stage']=stage
            for group,v in p.items():
                for side in ('r','l'):cfg['motor_scales'][group+'_'+side]=float(v)
            target=case/'INPUT.json';target.write_text(json.dumps(cfg,indent=2));protocol='hunt-air'+str(stage) if stage<3 else 'hunt-stand' if stage==3 else 'hunt-ground'
            with (case.parent/'simulation_output.log').open('a',encoding='utf-8') as log:
                subprocess.run([MYO,str(HERE/'live_tendon_trial.py'),'--protocol',protocol,'--physics','compliant','--source-task-basis','--config',str(target),'--out',str(case)],stdout=log,stderr=subprocess.STDOUT,timeout=180,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        b=json.loads((case/'BACKEND_RESULT.json').read_text());r=json.loads((case/'controller.result.json').read_text()) if (case/'controller.result.json').exists() else None
        ok=b['guard'] is None and r is not None and r['assessment']['accepted'];row=dict(stage=stage,path=str(case),accepted=ok,guard=b['guard'],assessment=None if r is None else r['assessment']);rows.append(row)
        if not ok:
            if r is None:
                z=np.load(case/'backend_trace.npz');row['last_physical_time_s']=float(z['t'][-1]);value=2000.-10.*row['last_physical_time_s']
            elif stage==3:value=1000.+1000.*r['assessment'].get('sway_m',1.)
            elif stage<3:value=1000.-r['assessment'].get('score',-1000.)
            else:
                swings=[f['sustained_swings'] for f in r.get('foot_clearance',{}).get('feet',[])];row['sustained_swings']=swings
                value=500.+100.*sum(max(0,3-n) for n in swings)+100.*r['assessment'].get('cycle_phase_error',.5)
            break
        if stage==4:
            swings=[f['sustained_swings'] for f in r['foot_clearance']['feet']];row['sustained_swings']=swings;row['actual_assistance_bw']=r['ground_contacts']['mean_harness_vertical_bw'];row['contact_phase_error']=r['ground_contacts']['phase_error']
            accepted=True;value=-100.-5.*sum(swings)+20.*row['contact_phase_error']+10.*row['actual_assistance_bw']
    return dict(ok=True,value=float(value),diagnostics=dict(all_four_stages_accepted=accepted,stages=rows,parameters_tied=True,source_lengths_fixed=True))

def diagnose():
    import mujoco
    sys.path.insert(0,str(S));import runner as R
    m=mujoco.MjModel.from_xml_path(str(R.MODEL));joints={j+'_'+s:int(m.joint(j+'_'+s).qposadr[0]) for s in ('r','l') for j in ('hip_flexion','knee_angle','ankle_angle')}
    rows={}
    for mode in ('native','compliant'):
        for stage in (1,2):
            folder=BASE/(mode+'_s'+str(stage));b,r=load_result(folder);z=np.load(folder/'backend_trace.npz');t=z['t'];keep=(t>=5)&(t<=13);names=z['muscle_names'].tolist();key=mode+'_s'+str(stage)
            row=dict(assessment=r['assessment'],joint_ranges_deg={name:float(np.degrees(np.ptp(z['qpos'][keep,adr]))) for name,adr in joints.items()},muscles={})
            for family in ('psoas','iliacus','bifemsh','semimem','semiten','bifemlh','rect_fem','soleus','med_gas','lat_gas','tib_ant','ext_dig','ext_hal'):
                ids=[names.index(family+'_'+s) for s in ('r','l')];a=z['activation'][keep][:,ids];f=z['force'][keep][:,ids]
                item=dict(mean_activation=a.mean(axis=0).tolist(),mean_delivered_force_N=f.mean(axis=0).tolist(),p95_delivered_force_N=np.percentile(f,95,axis=0).tolist())
                if mode=='compliant':
                    comp=z['fiber_force_components'][keep][:,:,ids];item['mean_active_passive_damping_N']=comp.mean(axis=0).tolist()
                row['muscles'][family]=item
            if mode=='compliant':row['component_sum_max_error_N']=float(np.max(abs(z['fiber_force_components'].sum(axis=1)-z['force'])))
            rows[key]=row
    (BASE/'BASELINE_DIAGNOSIS.json').write_text(json.dumps(rows,indent=2));print(json.dumps({k:{'knee_cycles':v['assessment']['knee_cycles'],'ranges':v['joint_ranges_deg'],'component_error':v.get('component_sum_max_error_N')} for k,v in rows.items()},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--diagnose',action='store_true');p.add_argument('--make-manifest',action='store_true');a=p.parse_args()
    if a.diagnose:diagnose()
    if a.make_manifest:
        folder=ROOT/'hunt_motor_calibration_v1';folder.mkdir(parents=True,exist_ok=True)
        cfg=dict(adapter_path=str(Path(__file__).resolve()),output_dir=str(folder),max_trials=12,wall_budget_s=1800,trial_timeout_s=360,seed=3101031,n_startup_trials=6,fixed_parameters={},
            search_space={k:dict(low=lo,high=hi) for k,lo,hi in [('hip_flex',1.,2.5),('knee_flex',1.,2.5),('ankle_df',.75,2.5)]},
            initial_candidates=[dict(hip_flex=h,knee_flex=k,ankle_df=d) for h,k,d in [(1.,1.,1.),(1.5,1.5,1.5),(2.,2.,1.),(1.5,2.,2.)]])
        target=folder/'INPUT_MANIFEST.json'
        if target.exists():raise RuntimeError('Do not overwrite a study manifest')
        target.write_text(json.dumps(cfg,indent=2));print(target)
