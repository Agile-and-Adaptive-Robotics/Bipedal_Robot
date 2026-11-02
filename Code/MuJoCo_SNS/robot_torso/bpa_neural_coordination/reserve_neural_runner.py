"""Bounded reusable reserve-relay experiments; fixed model and recorded variants."""
import argparse,json,hashlib,subprocess,sys
from pathlib import Path
import numpy as np
from normalized_prepared_transfer import ROOT,OUT,F,configuration,SpikeInputForward
from inverse_feasibility_policy import install
from pressure_reserve_neural import ReserveForward,reserve_adapter,CHANNELS
from gait_standing_objective import score
DEST=ROOT/'pressure_reserve_neural'
CONDITIONS={'nominal':'VIRTUAL_GAIT_BPA_CONFIGURATION.json','delay':'HELDOUT_DELAY.json','map':'HELDOUT_MAPERROR.json','harsh':'HELDOUT_SLOW_WEAK_MEMORY.json','pulse':'VIRTUAL_GAIT_BPA_CONFIGURATION.json'}

def run(cap,preview,condition,encoding,tag,duration=20.):
 if '.' in tag:raise ValueError('Use filesystem-safe tags without decimal points (0p05), because suffix-based result paths would collide')
 cfg=configuration();p=cfg['global_parameters'];DEST.mkdir(exist_ok=True);prefix=DEST/('w2lvar_'+tag)
 if prefix.with_suffix('.comparison.json').exists():raise FileExistsError(prefix)
 context={'active':False};module,original,instances=install('nearest');reserve_adapter(module.NativeBPAForce,module,context);wrappers=[];spikes=[]
 def transform(inner,names):
  if encoding=='spike':inner=SpikeInputForward(inner,names,population_size=16);spikes.append(inner)
  wrapper=ReserveForward(inner,names,context,cap,preview);wrappers.append(wrapper);return wrapper
 try:r=F.run('w2lvar',DEST,duration=duration,tag=tag,gain=p['pf_mn_gain'],ib_group_gain=p['ib_group_gain'],standing_bias_scale=p['standing_bias_scale'],pathway_group_scales=cfg['pathway_group_scales'],capture_neural_inputs=True,bpa_configuration=ROOT/'virtual_bpa_design'/CONDITIONS[condition],compiled_transform=transform,pinned_actuator_settling=True,settling_reference_configuration=ROOT/'pinned_preparation/PREPARATION_REFERENCE.json',mechanical_push_bw=.05 if condition=='pulse' else 0.)
 finally:module.NativeBPAForce=original
 wrapper=wrappers[0];np.savez_compressed(prefix.with_suffix('.relay.npz'),samples=np.asarray(wrapper.logs),channels=np.array(CHANNELS),routes=np.array(wrapper.routes),current_cap_nA=cap,preview_s=preview,neural_dt_s=.002)
 np.savez_compressed(prefix.with_suffix('.final_neural_state.npz'),**wrapper.export_state())
 bank=instances[0].bank;np.savez_compressed(prefix.with_suffix('.policy.npz'),events=np.asarray(bank.policy_events).reshape(-1,11),last_policy=bank.last_policy,policy_steps=bank.policy_steps)
 baseline=OUT/('w2lvar_bpa_'+encoding+'_v1');common={}
 with np.load(prefix.with_suffix('.release_state.npz')) as a,np.load(baseline.with_suffix('.release_state.npz')) as b:
  for k in b.files:
   # A reserve wrapper around sensory-spiking adds one explicit prefix layer.
   key='transformed_inner_'+k[len('transformed_'):] if encoding=='spike' and k.startswith('transformed_') else k
   common[k]=bool(np.array_equal(a[key],b[k]))
 if not all(common.values()):raise ValueError('Preparation changed before relay activation: '+str([k for k,v in common.items() if not v]))
 samples=np.asarray(wrapper.logs);active=samples[:,:,5];report=dict(cap_nA=cap,preview_s=preview,condition=condition,encoding=encoding,result=str(prefix.with_suffix('.result.json')),assessment=score(r,prefix)[1],legacy_complete_release_common_exact=common,circuit=wrapper.circuit(),maximum_added_current_nA=float(active.max()),integrated_added_current_nA_s=float(active.sum()*.002),nonzero_current_samples=int(np.count_nonzero(active)),complete_restore_selftest=wrapper.restore_verified,checkpoint_sha256=hashlib.sha256((ROOT/'group_diagnosis/NORMALIZED_STANDING_CHECKPOINT.json').read_bytes()).hexdigest(),spiking=spikes[0].metadata() if spikes else None,scope='New explicit two-relay engineering architecture. Frozen underlying neural parameters/geometry; actual pressure evolves normally. Same nearest-static reference policy. No pressure projection/force override/hidden MN voltage edit.')
 prefix.with_suffix('.comparison.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(tag=tag,assessment=report['assessment'],moving_s=r['completed_samples']*.002,max_current=report['maximum_added_current_nA'],restore_verified=wrapper.restore_verified,failure=r['failure']),indent=2))

def fit():
 DEST.mkdir(exist_ok=True)
 if not (DEST/'ZERO_GAIN_PARITY.json').exists():raise ValueError('Exact zero-gain gate required')
 plan={'caps_nA':[.05,.1,.2],'preview_s':[.05,.1],'duration_s':20.,'training':'nominal only, existing late pectineus failure','condition':'nominal','encoding':'graded','selection':'accepted only, then smallest current cap, then longest preview; held-outs not used for selection','budget':6,'current_relay':'2mirrored cells,6ms,1N sensor scale;10ms pressure derivative filter;40mV excitatory reversal'}
 path=DEST/'PREDECLARED_FIT.json'
 if path.exists() and json.loads(path.read_text())!=plan:raise ValueError('Fit plan changed')
 path.write_text(json.dumps(plan,indent=2),encoding='utf-8')
 for cap in plan['caps_nA']:
  for preview in plan['preview_s']:
   tag=f'fit_c{cap:g}_h{preview:g}'.replace('.','p');target=DEST/('w2lvar_'+tag+'.comparison.json')
   legacy=DEST/'w2lvar_fit_c0.05_h0.comparison.json'
   if cap==.05 and preview==.05 and legacy.exists():
    saved=json.loads(legacy.read_text())
    if saved['cap_nA']!=cap or saved['preview_s']!=preview:raise ValueError('Legacy first-trial metadata mismatch')
    target=legacy
   if not target.exists():
    with (DEST/(tag+'.log')).open('w',encoding='utf-8') as log:process=subprocess.run([sys.executable,str(Path(__file__)),'--cap',str(cap),'--preview',str(preview),'--tag',tag],stdout=log,stderr=subprocess.STDOUT)
    if process.returncode:raise RuntimeError('Trial process failed: '+tag)
   c=json.loads(target.read_text());print(json.dumps(dict(tag=tag,assessment=c['assessment'])),flush=True)
 winners=[json.loads(p.read_text()) for p in DEST.glob('w2lvar_fit_*.comparison.json') if json.loads(p.read_text())['assessment']['accepted']]
 if winners:
  best=sorted(winners,key=lambda c:(c['cap_nA'],-c['preview_s']))[0];(DEST/'SELECTED.json').write_text(json.dumps(best,indent=2),encoding='utf-8')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--cap',type=float,default=0.);p.add_argument('--preview',type=float,default=.1);p.add_argument('--condition',choices=CONDITIONS,default='nominal');p.add_argument('--encoding',choices=['graded','spike'],default='graded');p.add_argument('--tag',default='zero_parity');p.add_argument('--duration',type=float,default=20.);p.add_argument('--fit',action='store_true');a=p.parse_args();fit() if a.fit else run(a.cap,a.preview,a.condition,a.encoding,a.tag,a.duration)
