"""Reusable predeclared reserve-to-valve reference experiment runner."""
import argparse,json,hashlib,subprocess,sys
from pathlib import Path
import numpy as np
from normalized_prepared_transfer import ROOT,OUT,F,configuration,SpikeInputForward
from inverse_feasibility_policy import install
from pressure_valve_neural import ValveReserveForward,install_valve_adapter
from pressure_reserve_neural import CHANNELS
from reserve_neural_runner import CONDITIONS
from gait_standing_objective import score
DEST=ROOT/'pressure_reserve_valve'

def run(cap,preview,condition,encoding,tag,duration=20.,reserve_target=0.,output_dir=None,actuator_configuration=None,preparation_configuration=None,common_release_validator=None,sensor_gain=1.,mechanical_push_bw=None):
 if '.' in tag:raise ValueError('Use filesystem-safe tags with p instead of decimal points')
 if not .5<=sensor_gain<=2.:raise ValueError('Declared reserve sensor gain range is0.5..2')
 cfg=configuration();p=cfg['global_parameters'];out=Path(output_dir) if output_dir is not None else DEST;out.mkdir(parents=True,exist_ok=True);prefix=out/('w2lvar_'+tag)
 if prefix.with_suffix('.comparison.json').exists():raise FileExistsError(prefix)
 context={'active':False};module,original,instances=install('nearest');install_valve_adapter(module.NativeBPAForce,module,context);wrappers=[];encoders=[]
 def transform(inner,names):
  if encoding=='spike':inner=SpikeInputForward(inner,names,population_size=16);encoders.append(inner)
  wrapper=ValveReserveForward(inner,names,context,cap,preview,reserve_target);wrapper.force_scale=1./sensor_gain;wrappers.append(wrapper);return wrapper
 actual=Path(actuator_configuration) if actuator_configuration is not None else ROOT/'virtual_bpa_design'/CONDITIONS[condition];preparation=Path(preparation_configuration) if preparation_configuration is not None else ROOT/'pinned_preparation/PREPARATION_REFERENCE.json';push=(.05 if condition=='pulse' else 0.) if mechanical_push_bw is None else mechanical_push_bw
 try:r=F.run('w2lvar',out,duration=duration,tag=tag,gain=p['pf_mn_gain'],ib_group_gain=p['ib_group_gain'],standing_bias_scale=p['standing_bias_scale'],pathway_group_scales=cfg['pathway_group_scales'],capture_neural_inputs=True,bpa_configuration=actual,compiled_transform=transform,pinned_actuator_settling=True,settling_reference_configuration=preparation,mechanical_push_bw=push)
 finally:module.NativeBPAForce=original
 wrapper=wrappers[0];bank=instances[0].bank;data=context['data']
 np.savez_compressed(prefix.with_suffix('.relay.npz'),samples=np.asarray(wrapper.logs),channels=np.array(CHANNELS),routes=np.array(wrapper.routes),pressure_cap_kpa=cap,preview_s=preview,neural_dt_s=.002)
 np.savez_compressed(prefix.with_suffix('.valve.npz'),samples=np.asarray(bank.valve_logs),channels=np.array(['original_nominal_reference_kpa','neural_requested_increment_kpa','applied_increment_kpa','final_reference_kpa']),routes=np.array(wrapper.routes),physics_dt_s=.002)
 np.savez_compressed(prefix.with_suffix('.final_neural_state.npz'),**wrapper.export_state())
 final=dict(pressure=bank.pressure,command=bank.command,history=bank.history,force=bank.force,flags=bank.flags,step_count=np.array(bank.step_count),direction_state=instances[0].direction_state,nominal_force=instances[0].nominal_force,policy_steps=np.array(bank.policy_steps),valve_steps=np.array(bank.valve_steps),qpos=data.qpos,qvel=data.qvel,act=data.act,ctrl=data.ctrl,qacc_warmstart=data.qacc_warmstart,physics_time=np.array(data.time))
 np.savez_compressed(prefix.with_suffix('.final_actuator_mechanics_state.npz'),**final)
 baseline=OUT/('w2lvar_bpa_'+encoding+'_v1');common={}
 if common_release_validator is not None:common=common_release_validator(prefix.with_suffix('.release_state.npz'))
 else:
  with np.load(prefix.with_suffix('.release_state.npz')) as a,np.load(baseline.with_suffix('.release_state.npz')) as b:
   for k in b.files:
    key='transformed_inner_'+k[len('transformed_'):] if encoding=='spike' and k.startswith('transformed_') else k;common[k]=bool(np.array_equal(a[key],b[k]))
 if not all(common.values()):raise ValueError('Common preparation changed')
 valve=np.asarray(bank.valve_logs);report=dict(pressure_cap_kpa=cap,preview_s=preview,reserve_target_N=reserve_target,sensor_gain=sensor_gain,condition=condition,encoding=encoding,result=str(prefix.with_suffix('.result.json')),assessment=score(r,prefix)[1],common_release_exact=common,circuit=wrapper.circuit(),maximum_applied_increment_kpa=float(np.max(valve[:,:,2])) if len(valve) else 0.,integrated_applied_increment_kpa_s=float(np.sum(valve[:,:,2])*.002) if len(valve) else 0.,complete_neural_restore_verified=wrapper.restore_verified,spiking=encoders[0].metadata() if encoders else None,source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ('pressure_valve_neural.py','pressure_reserve_neural.py','inverse_feasibility_policy.py')},scope='Distinct explicitly neural valve-reference output map. Original NN supplies nominal reference; two reserve neurons add bounded reference only. Actual pressure/force/geometry and forward guard unchanged; no direct MN current or state projection.')
 prefix.with_suffix('.comparison.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(tag=tag,assessment=report['assessment'],moving_s=r['completed_samples']*.002,maximum_increment=report['maximum_applied_increment_kpa'],restore=wrapper.restore_verified,failure=r['failure']),indent=2))
 return report

def fit():
 DEST.mkdir(exist_ok=True)
 if not (DEST/'ZERO_GAIN_PARITY.json').exists():raise ValueError('Zero-gain physical parity gate required')
 plan=dict(pressure_caps_kpa=[.5,1.,2.],previews_s=[.1],reserve_targets_N=[0.,1.],budget=6,duration_s=20.,training='nominal graded only, existing late pressure-domain failure',selection='accepted first; smallest pressure cap then smallest reserve target; no held-out selection',relay_current_cap_nA=10.,relay_tau_s=.006,pressure_filter_tau_s=.01,reserve_setpoint_scope='Explicit neural error setpoint:1N equivalent is2.5micrometers at400kN/m, not an imposed actual force/pressure floor')
 path=DEST/'PREDECLARED_FIT.json'
 if path.exists() and json.loads(path.read_text())!=plan:raise ValueError('Plan changed')
 path.write_text(json.dumps(plan,indent=2),encoding='utf-8')
 for cap in plan['pressure_caps_kpa']:
  for reserve_target in plan['reserve_targets_N']:
   preview=.1;tag=f'fit_p{cap:g}_r{reserve_target:g}'.replace('.','p');target=DEST/('w2lvar_'+tag+'.comparison.json')
   if not target.exists():
    with (DEST/(tag+'.log')).open('w',encoding='utf-8') as log:process=subprocess.run([sys.executable,str(Path(__file__)),'--cap',str(cap),'--preview',str(preview),'--reserve-target',str(reserve_target),'--tag',tag],stdout=log,stderr=subprocess.STDOUT)
    if process.returncode:raise RuntimeError(tag)
   c=json.loads(target.read_text());print(json.dumps(dict(tag=tag,assessment=c['assessment'])),flush=True)
 winners=[json.loads(p.read_text()) for p in DEST.glob('w2lvar_fit_*.comparison.json') if json.loads(p.read_text())['assessment']['accepted']]
 if winners:(DEST/'SELECTED.json').write_text(json.dumps(sorted(winners,key=lambda c:(c['pressure_cap_kpa'],c['reserve_target_N']))[0],indent=2),encoding='utf-8')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--cap',type=float,default=0.);p.add_argument('--preview',type=float,default=.1);p.add_argument('--reserve-target',type=float,default=0.);p.add_argument('--condition',choices=CONDITIONS,default='nominal');p.add_argument('--encoding',choices=['graded','spike'],default='graded');p.add_argument('--tag',default='zero_parity');p.add_argument('--duration',type=float,default=20.);p.add_argument('--fit',action='store_true');a=p.parse_args();fit() if a.fit else run(a.cap,a.preview,a.condition,a.encoding,a.tag,a.duration,a.reserve_target)
