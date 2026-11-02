"""Frozen normalized-w2lvar biological/pressure and graded/sensory-spiking comparison."""
import sys,json,argparse,hashlib,subprocess
from pathlib import Path
import numpy as np
import gait_standing_family as F
from gait_standing_objective import score
sys.path.insert(0,str(F.SPINAL/'virtual_walker_20261031'))
from spiking_sns_inputs import SpikeInputForward
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');OUT=ROOT/'normalized_prepared_transfer';CHECKPOINT=ROOT/'group_diagnosis/NORMALIZED_STANDING_CHECKPOINT.json'

def configuration():
 cfg=json.loads(CHECKPOINT.read_text())
 for name,expected in cfg['expected_source_sha256'].items():
  path=Path(F.__file__) if name==Path(F.__file__).name else F.SPINAL/name
  if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise ValueError('Controller source changed; revalidate: '+str(path))
 for path,expected in cfg['expected_dependency_sha256'].items():
  if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=expected:raise ValueError('Dependency changed; revalidate: '+path)
 return cfg

def run(case):
 cfg=configuration();OUT.mkdir(exist_ok=True);allowed={'bio_graded','bio_spike','bpa_graded','bpa_spike','bpa_graded_fresh','bpa_spike_fresh','bpa_graded_replay','bpa_spike_replay','bio_graded_fresh','bio_spike_fresh','harsh_graded','harsh_spike','harsh_graded_replay','harsh_spike_replay','push_graded','push_spike','push_graded_replay','push_spike_replay','rot0p8_graded','rot0p8_spike','rot0p8_graded_fresh','rot0p8_spike_fresh','rot0p8_hold20_graded','rot0p8_hold20_spike','hold20_graded','hold20_spike'}
 if case not in allowed:raise ValueError(case)
 biological=case.startswith('bio');spike='spike' in case;wrappers=[]
 def transform(inner,names):
  obj=SpikeInputForward(inner,names,population_size=16);wrappers.append(obj);return obj
 actual=None if biological else ROOT/'virtual_bpa_design'/('HELDOUT_SLOW_WEAK_MEMORY.json' if case.startswith('harsh') else 'VIRTUAL_GAIT_BPA_CONFIGURATION.json')
 nominal=OUT/('w2lvar_bpa_'+('spike' if spike else 'graded')+'_v1')
 replay=nominal.with_suffix('.inputs.npz') if case.endswith('replay') else None
 if replay and not replay.exists():raise FileNotFoundError(replay)
 p=cfg['global_parameters'];r=F.run('w2lvar',OUT,duration=20. if 'hold20' in case else 8.,tag=case+'_v1',gain=p['pf_mn_gain'],ib_group_gain=p['ib_group_gain'],standing_bias_scale=p['standing_bias_scale'],pathway_group_scales=cfg['pathway_group_scales'],capture_neural_inputs=True,bpa_configuration=actual,compiled_transform=transform if spike else None,pinned_actuator_settling=True,settling_reference_configuration=ROOT/'pinned_preparation/PREPARATION_REFERENCE.json' if actual else None,sensory_replay=replay,mechanical_push_bw=.05 if case.startswith('push') else 0.,rig_rotational_scale=.8 if case.startswith('rot0p8') else 1.)
 prefix=OUT/('w2lvar_'+case+'_v1');report=dict(case=case,result=str(prefix.with_suffix('.result.json')),assessment=score(r,prefix)[1],spiking=wrappers[0].metadata() if wrappers else None,checkpoint=str(CHECKPOINT),checkpoint_sha256=hashlib.sha256(CHECKPOINT.read_bytes()).hexdigest(),scope='Frozen normalized-w2lvar standing configuration. Same0.6s pinned preparation policy; actuator/encoder types naturally acquire different internal states. All complete-state equality claims are within a given nominal/live-replay/fresh treatment, not across biological versus pressure states. Virtual actuator design unchanged, no postfailure resizing.')
 if wrappers:
  np.savez_compressed(prefix.with_suffix('.spikes.npz'),counts=np.asarray(wrappers[0].logs),names=np.array([wrappers[0].names[i] for i in wrappers[0].indices]),neural_step_s=.002,cells_per_population=16)
  np.savez_compressed(prefix.with_suffix('.final_spiking_state.npz'),**wrappers[0].export_state())
 prefix.with_suffix('.comparison.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

def sequence():
 OUT.mkdir(exist_ok=True)
 for case in ('bio_graded','bio_spike','bpa_graded','bpa_spike'):
  path=OUT/('w2lvar_'+case+'_v1.comparison.json')
  if not path.exists():
   with (OUT/(case+'.log')).open('w',encoding='utf-8') as log:process=subprocess.run([sys.executable,str(Path(__file__)),'--case',case],stdout=log,stderr=subprocess.STDOUT)
   if process.returncode:raise RuntimeError(case+' failed; inspect preserved log')
  report=json.loads(path.read_text());print(json.dumps(dict(case=case,assessment=report['assessment'])),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case');a=p.parse_args();run(a.case) if a.case else sequence()
