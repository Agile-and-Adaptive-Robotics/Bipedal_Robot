"""Matched opt-in sensory-spiking transfer of the calibrated syn6 standing case."""
import os,sys,argparse,json,contextlib,hashlib
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['AARL_NET']='syn6';os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
from spiking_sns_inputs import SpikeInputForward
HERE=Path(__file__).resolve().parent;COORD=HERE.parent.parent/'robot_torso/bpa_neural_coordination';sys.path.insert(0,str(COORD))
import gait_standing_family as F
from gait_standing_objective import score
p=argparse.ArgumentParser();p.add_argument('--mode',choices=['graded','hybrid'],required=True);p.add_argument('--tag',default='sensory_v1');p.add_argument('--bpa-configuration');p.add_argument('--verify-state',action='store_true');p.add_argument('--out',default=r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\syn6_standing');a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
source=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\rustuna_syn6_v1\BEST.json');parameters=json.loads(source.read_text())['trial']['params'];saved=[];state_check={}
def transform(inner,names):
 c=SpikeInputForward(inner,names,population_size=16)
 if a.verify_state:
  initial=c.export_state();inputs=[np.full(inner.c.input_connectivity.shape[1],.8+.7*np.sin(j*.41)) for j in range(32)]
  first=np.array([c.forward(u).copy() for u in inputs]);end=c.export_state();c.import_state(initial)
  second=np.array([c.forward(u).copy() for u in inputs]);again=c.export_state()
  state_check.update(output_byte_parity=first.tobytes()==second.tobytes(),state_byte_parity={k:bool(v.shape==again[k].shape and v.dtype==again[k].dtype and v.tobytes()==again[k].tobytes()) for k,v in end.items()},actual_spikes_during_probe=int(c.total_counts.sum()),scope='32 actual sparse-backend steps; restore initial state after verification, including encoder, Ca gates, caches, counters, and history. External mechanics and runner counters are not owned here.')
  assert state_check['output_byte_parity'] and all(state_check['state_byte_parity'].values());c.import_state(initial)
 saved.append(c);return c
tag=a.tag+'_'+a.mode;prefix=out/('syn6_'+tag)
if prefix.with_suffix('.SPIKING.json').exists():raise FileExistsError('Preserving completed trial; use a new tag')
with (out/(tag+'.log')).open('w',encoding='utf-8') as log,contextlib.redirect_stdout(log):
 r=F.run('syn6',out,duration=8.,tag=tag,gain=parameters['pf_mn_gain'],standing_bias_scale=parameters['standing_bias_scale'],ib_group_gain=parameters['ib_group_gain'],bpa_configuration=a.bpa_configuration,compiled_transform=transform if a.mode=='hybrid' else None)
 value,assessment=score(r,prefix)
record=dict(state_snapshot_verification=state_check,mode=a.mode,source_parameters=str(source),source_parameters_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),parameters=parameters,result=r,assessment=assessment,scope='Same calibrated syn6 standing controller, with declared virtual BPA pressure dynamics' if a.bpa_configuration else 'Same calibrated syn6 standing controller and biological force law; input-only sensory spike encoding only. Not the Hunt controller, stepping, or harness-free balance.')
if a.bpa_configuration:record.update(bpa_configuration=str(Path(a.bpa_configuration).resolve()),bpa_configuration_sha256=hashlib.sha256(Path(a.bpa_configuration).read_bytes()).hexdigest(),scope='Sensory-spiking syn6 with declared virtual BPA pressure/force dynamics on Gait2392. Assumed dimensions and actuator dynamics, not identified CAD hardware. Rig remains; standing only, not the Hunt controller.')
if saved:
 c=saved[-1];record['spiking']=c.metadata();np.savez_compressed(out/(tag+'.spikes.npz'),spike_counts=np.array(c.logs),names=record['spiking']['spiking_source_names'])
prefix.with_suffix('.SPIKING.json').write_text(json.dumps(record,indent=2));print(json.dumps(dict(mode=a.mode,assessment=assessment,upward_BW=r['mean_upward_rig_support_bw'],downward_BW=r['mean_downward_rig_force_bw'],net_BW=r['mean_net_rig_force_bw'],spiking_cells=record.get('spiking',{}).get('actual_spiking_cells',0)),indent=2))
