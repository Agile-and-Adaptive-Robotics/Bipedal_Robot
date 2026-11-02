"""Frozen matched preparation policy; original central circuit remains graded."""
import sys,json,argparse
import numpy as np
from pathlib import Path
import gait_standing_family as F
from gait_standing_objective import score
sys.path.insert(0,str(F.SPINAL/'virtual_walker_20261031'))
from spiking_sns_inputs import SpikeInputForward
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')
VERSION='v2' # v1 preserves a failed release-snapshot scope check; no physics result.

def main(case):
 allowed={'bpa_spike','bpa_spike_fresh','bpa_spike_replay','bio_spike','bio_spike_fresh','harsh_spike','harsh_spike_replay','graded_bpa_current'}
 if case not in allowed:raise ValueError(case)
 out=ROOT/'prepared_spiking';out.mkdir(exist_ok=True);wrappers=[]
 def transform(inner,names):
  wrapper=SpikeInputForward(inner,names,population_size=16);wrappers.append(wrapper);return wrapper
 params=json.loads((ROOT/'rustuna_syn6_v1/BEST.json').read_text())['trial']['params']
 config=ROOT/'virtual_bpa_design'/('HELDOUT_SLOW_WEAK_MEMORY.json' if case.startswith('harsh') else 'VIRTUAL_GAIT_BPA_CONFIGURATION.json')
 if case.startswith('bio'):config=None
 replay=out/('syn6_bpa_spike_'+VERSION+'.inputs.npz') if case.endswith('replay') else None
 if replay and not replay.exists():raise FileNotFoundError(replay)
 r=F.run('syn6',out,duration=8.,tag=case+'_'+VERSION,gain=params['pf_mn_gain'],ib_group_gain=params['ib_group_gain'],standing_bias_scale=params['standing_bias_scale'],bpa_configuration=config,sensory_replay=replay,compiled_transform=None if case=='graded_bpa_current' else transform,pinned_actuator_settling=True,settling_reference_configuration=ROOT/'pinned_preparation/PREPARATION_REFERENCE.json' if config else None)
 prefix=out/('syn6_'+case+'_'+VERSION);result=dict(case=case,result=str(prefix.with_suffix('.result.json')),assessment=score(r,prefix)[1],spiking=wrappers[0].metadata() if wrappers else None,scope='Frozen central syn6 gains; sensory-only spiking substitution; 0.6s pinned preparation; unchanged baseline rig. Model/actuator types naturally acquire different preparation states. Common-state equality applies within each nominal/live-replay pair, not across biological versus BPA physics.')
 if wrappers:
  np.savez_compressed(prefix.with_suffix('.spikes.npz'),counts=np.asarray(wrappers[0].logs),names=np.array([wrappers[0].names[i] for i in wrappers[0].indices]),neural_step_s=.002,cells_per_population=16)
  np.savez_compressed(prefix.with_suffix('.final_spiking_state.npz'),**wrappers[0].export_state())
 prefix.with_suffix('.comparison.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(case=case,assessment=result['assessment'],spiking_cells=result['spiking']['actual_spiking_cells'] if wrappers else 0,failure=r['failure']),indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case',required=True);main(p.parse_args().case)
