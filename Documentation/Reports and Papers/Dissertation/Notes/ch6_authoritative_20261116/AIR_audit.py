"""Recompute bounded saved-recording checks; never integrate or fit a model."""
from pathlib import Path
import hashlib
import json
import time
import importlib.util
import numpy as np
import mujoco

OUT=Path(__file__).resolve().parent
T=Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104')
W=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\LiDeng_W2L_20261104')
sources={}
def bind(p):
    p=Path(p);sources[str(p)]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size};return p
def js(p):return json.loads(bind(p).read_text(encoding='utf-8-sig'))
def nz(p):return np.load(bind(p),allow_pickle=False)
start=time.perf_counter()
summary={}
air=T/'hunt_pose_air_validation_v3'
complete=js(air/'COMPLETE.json')
cfg=js(T/'ordered_output_fits_v5/hunt/C/MANIFEST.json')['configuration']
summary['command_conditions']={'sustained_normalized_FE_stage1':0.,'sustained_normalized_FE_stage2':[cfg['command_F'],cfg['command_E']], 'common_initial_seed':'+0.4 E_r and F_l for first 50 neural calls at 2 ms (=100 ms); initialization, not a stride clock','stage1_sensory_enabled':False,'stage2_sensory_enabled':True,'both_contact_enabled':False,'POST_BAL_FUNC':'clamped off; rg_tonic_current=0; old DRIVE disabled by Runtime.prepare_inputs'}
summary['air']={}
cycle_path=bind(T/'source_snapshot_final/controller_profiles_20261104/dominant_cycles.py')
spec=importlib.util.spec_from_file_location('saved_dominant_cycles',cycle_path)
cycles=importlib.util.module_from_spec(spec)
try:
    spec.loader.exec_module(cycles)
    summary['cycle_detector_self_check_passed']=cycles.self_check()['passed']
except PermissionError as error:
    summary['cycle_recalculation_blocked']=str(error)
    cycles=None
for stage in (1,2):
    d=air/f'stage{stage}';a=nz(d/'body.assistance.npz');z=nz(d/'capture/aligned_trace.npz');m=mujoco.MjModel.from_binary_path(str(bind(d/'capture/model.mjb')));meta=js(d/'capture/ALIGNED_TRACE_METADATA.json')
    pelvis=m.body('pelvis').id
    legacy_bw=float(m.body_mass.sum()*9.81);subtree_bw=float(m.body_subtreemass[pelvis]*9.81)
    cols=a['columns'].tolist();s=a['samples'];late=s[:,0]>=5
    force=float(s[late,1].mean())
    t=z['t_state'];mask=t>=5
    anatomical_mean=float(z['harness_wrench_anatomical'][mask,1].mean())
    jn=z['joint_names'].tolist();qadr=z['joint_qpos_addresses'];knee_columns=[qadr[jn.index('knee_angle_'+side)]for side in ('r','l')]
    summary['air'][str(stage)]={'metadata':meta,'sample_columns':cols,'legacy_total_mass_kg':float(m.body_mass.sum()),'pelvis_subtree_mass_kg':float(m.body_subtreemass[pelvis]),'difference_ground_mass_kg':float(m.body('ground').mass[0]),'legacy_BW_N':legacy_bw,'pelvis_subtree_BW_N':subtree_bw,'legacy_support_window':'body.assistance.samples[:,0] >= 5 s (5.000 through15.998 s, 5500 samples)','support_N':force,'legacy_support_percent':100*force/legacy_bw,'same_force_subtree_percent':100*force/subtree_bw,'aligned_force_mean_N_tstate_ge5':anatomical_mean,'aligned_support_percent_tstate_ge5':100*anatomical_mean/subtree_bw,'foot_force_max_abs_N':float(np.max(np.abs(z['foot_force_anatomical_N']))),'aligned_knee_max_deg_RL':np.degrees(z['qpos'][:,knee_columns].max(axis=0)).tolist(),'adopted':complete[stage-1]}
    keep=(t>=5)&(t<14)
    hipcols=[qadr[jn.index('hip_flexion_'+side)]for side in ('r','l')]
    if cycles is not None:
        recomputed=cycles.measure(t[keep],-np.degrees(z['qpos'][keep][:,knee_columns]),np.degrees(z['qpos'][keep][:,hipcols]))
        summary['air'][str(stage)]['principal_cycles_recomputed_now']=recomputed
    summary['air'][str(stage)]['cycle_window']='5 <= t_state < 14 s; knee sign inverted only for flexion-peak detection, raw coordinates retained'

base_dir=T/'hunt_layer_stimulus_2nA_v2/Baseline';base=nz(base_dir/'capture/aligned_trace.npz');bn=base['neural_names'].tolist()
summary['stimuli']=[]
for label,root in [('silent',T/'hunt_layer_stimulus_2nA_v2'),('active',T/'hunt_active_phase_stimulus_v3')]:
    for layer in ['RG','PF','MN','PF_inhibition']:
        d=root/layer;config=js(d/'body.json');z=nz(d/'capture/aligned_trace.npz');names=z['neural_names'].tolist();pulse=z['stimulus_nA'];on=np.where(np.any(np.nan_to_num(pulse)!=0,axis=1))[0];onset=float(z['t_neural_input'][on[0]]);after=(z['t_neural']>=onset)&(z['t_neural']<onset+3.4);before=z['t_neural']<onset
        r={'phase':label,'layer':layer,'target':config['neural_stimulus']['cell'],'onset_s':onset,'last_input_s':float(z['t_neural_input'][on[-1]]),'pulse_samples':len(on),'current_nA':config['neural_stimulus']['current_nA'],'response_window_s':[onset,onset+3.4],'response_end_exclusive':True,'pre_pulse_max_abs_mV':0.,'maximum_absolute_differences_mV':{},'peak_times_s':{}}
        for cell in ['RG_F_r','PF_HIP-F_r','MN_iliacus_r']:
            y=z['neural_voltage'][:,names.index(cell)];ref=np.interp(z['t_neural'],base['t_neural'],base['neural_voltage'][:,bn.index(cell)]);delta=np.abs(y-ref)
            r['pre_pulse_max_abs_mV']=max(r['pre_pulse_max_abs_mV'],float(delta[before].max()))
            r['maximum_absolute_differences_mV'][cell]=float(delta[after].max());r['peak_times_s'][cell]=float(z['t_neural'][after][np.argmax(delta[after])])
        pf=z['neural_voltage'][:,names.index('PF_HIP-F_r')];during=(z['t_neural']>=onset)&(z['t_neural']<onset+.1)
        r['PF_voltage_range_during_pulse_mV']=[float(pf[during].min()),float(pf[during].max())];r['PF_transmission_range_during_pulse']=np.clip(np.array(r['PF_voltage_range_during_pulse_mV'])/5.,0,1).tolist()
        summary['stimuli'].append(r)

fresh=js(W/'FRESH_VERIFICATION.json');wd=W/'fresh_repeat/W2L_VirtualWalker_Newer/repeat';z=nz(wd/'trajectory.npz');result=fresh['result'];bw=result['robot_mass_kg']*9.81;mask=z['t']>=1
summary['W2L108']={'robot_mass_kg':result['robot_mass_kg'],'BW_N':bw,'joint_names':z['joint_names'].tolist(),'minimum_saved_angles_deg_by_column':z['angles'][mask].min(axis=0).tolist(),'maximum_saved_angles_deg_by_column':z['angles'][mask].max(axis=0).tolist(),'result_keys':list(result),'result_knee_min_deg':result.get('anatomical_knee_min_deg'),'assessment':fresh['assessment'],'source_model_xml_sha256':hashlib.sha256(bind(wd/'model.xml').read_bytes()).hexdigest(),'window_after1':[float(z['t'][mask][0]),float(z['t'][mask][-1])],'upward_percent_after1':float(100*sum(np.maximum(z[k][mask,1],0).mean()for k in ['rail_force_world','guide_force_world','harness_force'])/bw),'signed_net_percent_whole':float(100*sum(z[k][:,1].mean()for k in ['rail_force_world','guide_force_world','harness_force'])/bw),'repeat_verified_now':None,'stored_parity':fresh['all_numeric_arrays_byte_exact']}
trial=W/f"body_trials/W2L_VirtualWalker_Newer/rustuna108_{fresh['selected_trial']}"
if not trial.exists():
    trial=next(p.parent for p in (W/'body_trials').rglob('trajectory.npz') if 'trial_0010_' in str(p.parent))
selected=nz(trial/'trajectory.npz')
summary['W2L108']['selected_path']=str(trial);summary['W2L108']['repeat_verified_now']={k:bool(np.array_equal(z[k],selected[k],equal_nan=True))for k in z.files if np.issubdtype(z[k].dtype,np.number)}
noise=js(T/'hunt_rg_noise_v1/COMPLETE.json');summary['noise']={'cases':[dict(label=r['label'],sigma_nA=r['sigma_nA'],seed=r['seed'],accepted=r['accepted'])for r in noise['cases']],'stored_repeat':noise['cases'][-1].get('fresh_repeat_max_abs_difference')}
original_noise=nz(T/'hunt_rg_noise_v1/0.025_41001/capture/aligned_trace.npz');repeat_noise=nz(T/'hunt_rg_noise_v1/repeat/capture/aligned_trace.npz')
summary['noise']['repeat_verified_now']={k:float(np.max(np.abs(original_noise[k]-repeat_noise[k])))for k in ['qpos','neural_voltage','stimulus_nA']}
for rel in ['source_snapshot_final/controller_profiles_20261104/body_probe.py','source_snapshot_final/controller_profiles_20261104/profile.py','source_snapshot_final/controller_profiles_20261104/stimulus_figures.py','source_snapshot_final/existing_dependencies/hunt_training__physical_eval.py']:
    bind(T/rel)
bind(W/'ADAPTER_SOURCE.py')
summary['elapsed_s']=time.perf_counter()-start;summary['sources']=sources
(OUT/'AIR_RESULTS.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps({'air':{s:{k:v for k,v in r.items()if k not in ['metadata','adopted']}for s,r in summary['air'].items()},'stimuli':summary['stimuli'],'W2L108':{k:v for k,v in summary['W2L108'].items()if k not in ['assessment']},'elapsed_s':summary['elapsed_s']},indent=2))
