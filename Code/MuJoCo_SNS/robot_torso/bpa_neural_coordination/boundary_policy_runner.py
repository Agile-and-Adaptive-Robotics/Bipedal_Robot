"""Reusable frozen-controller interface-policy comparison; one physics run per call."""
import argparse,json,hashlib,subprocess,sys
from pathlib import Path
import numpy as np
from normalized_prepared_transfer import F,ROOT,OUT,configuration,SpikeInputForward
from inverse_feasibility_policy import install,CHANNELS
from gait_standing_objective import score

def run(policy,encoding,duration,repeat):
 gate=ROOT/'inverse_boundary_policy/BOUNDARY_UNIT_CHECKS.json'
 unit=json.loads(gate.read_text())
 if not unit['forward_domain_unchanged']:raise ValueError('Unit gate missing')
 for path,expected in unit.get('source_sha256',{}).items():
  if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=expected:raise ValueError('Boundary source changed after unit gate: '+path)
 cfg=configuration();p=cfg['global_parameters'];out=gate.parent;tag=f'{policy}_{encoding}_{duration:g}s_{repeat}';prefix=out/('w2lvar_'+tag)
 if prefix.with_suffix('.comparison.json').exists():raise FileExistsError('Completed case preserved: '+str(prefix))
 module,original,instances=install(policy);wrappers=[]
 def transform(inner,names):
  obj=SpikeInputForward(inner,names,population_size=16);wrappers.append(obj);return obj
 try:
  r=F.run('w2lvar',out,duration=duration,tag=tag,gain=p['pf_mn_gain'],ib_group_gain=p['ib_group_gain'],standing_bias_scale=p['standing_bias_scale'],pathway_group_scales=cfg['pathway_group_scales'],capture_neural_inputs=True,bpa_configuration=ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json',compiled_transform=transform if encoding=='spike' else None,pinned_actuator_settling=True,settling_reference_configuration=ROOT/'pinned_preparation/PREPARATION_REFERENCE.json')
 finally:module.NativeBPAForce=original
 if len(instances)!=1:raise ValueError('Expected one pressure adapter')
 bank=instances[0].bank;events=np.asarray(bank.policy_events).reshape(-1,len(CHANNELS)+2)
 np.savez_compressed(prefix.with_suffix('.policy.npz'),events=events,channels=np.array(['policy_step','route_index']+CHANNELS),route_names=np.array(bank.names),last_policy=bank.last_policy,policy_steps=bank.policy_steps)
 if wrappers:
  np.savez_compressed(prefix.with_suffix('.spikes.npz'),counts=np.asarray(wrappers[0].logs),names=np.array([wrappers[0].names[i] for i in wrappers[0].indices]),neural_step_s=.002,cells_per_population=16)
  np.savez_compressed(prefix.with_suffix('.final_spiking_state.npz'),**wrappers[0].export_state())
 base=OUT/('w2lvar_bpa_'+encoding+'_v1');checks={}
 with np.load(prefix.with_suffix('.release_state.npz')) as a,np.load(base.with_suffix('.release_state.npz')) as b:
  if set(a.files)!=set(b.files):raise ValueError('Release state field mismatch')
  checks={k:bool(np.array_equal(a[k],b[k])) for k in a.files}
 if not all(checks.values()):raise ValueError('Common preparation state changed')
 report=dict(policy=policy,encoding=encoding,result=str(prefix.with_suffix('.result.json')),assessment=score(r,prefix)[1],policy_events=len(events),policy_steps=bank.policy_steps,legacy_release_state_exact=checks,checkpoint_sha256=hashlib.sha256((ROOT/'group_diagnosis/NORMALIZED_STANDING_CHECKPOINT.json').read_bytes()).hexdigest(),spiking=wrappers[0].metadata() if wrappers else None,scope='Frozen neural parameters/geometry, unchanged forward domain/actual pressure lag. Policy activates only after matched preparation. Static feasible-reference saturation is not a force override and cannot guarantee dynamic-domain safety. No saturation-error neural input added.')
 prefix.with_suffix('.comparison.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(case=tag,assessment=report['assessment'],events=len(events),free_s=r['completed_samples']*.002,failure=r['failure']),indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--policy',choices=['abort','nearest']);p.add_argument('--encoding',default='graded',choices=['graded','spike']);p.add_argument('--duration',type=float,default=20.);p.add_argument('--repeat',default='v1');p.add_argument('--batch',action='store_true');a=p.parse_args()
 if a.batch:
  out=ROOT/'inverse_boundary_policy'
  for encoding in ('graded','spike'):
   for policy,repeat in (('abort','v1'),('nearest','v1'),('nearest','fresh')):
    tag=f'{policy}_{encoding}_{a.duration:g}s_{repeat}';path=out/('w2lvar_'+tag+'.comparison.json')
    if not path.exists():
     with (out/(tag+'.log')).open('w',encoding='utf-8') as log:process=subprocess.run([sys.executable,str(Path(__file__)),'--policy',policy,'--encoding',encoding,'--duration',str(a.duration),'--repeat',repeat],stdout=log,stderr=subprocess.STDOUT)
     if process.returncode:raise RuntimeError('Case process failed: '+tag)
    c=json.loads(path.read_text());r=json.loads(Path(c['result']).read_text());print(json.dumps(dict(case=tag,assessment=c['assessment'],moving_s=r['completed_samples']*.002)),flush=True)
 else:
  if a.policy is None:p.error('--policy or --batch is required')
  run(a.policy,a.encoding,a.duration,a.repeat)
