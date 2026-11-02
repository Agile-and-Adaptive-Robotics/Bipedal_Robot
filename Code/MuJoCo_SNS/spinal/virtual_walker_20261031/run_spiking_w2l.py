"""Bounded external-output graded vs actual sensory-spiking W2L comparison."""
import argparse,contextlib,hashlib,json
from pathlib import Path
import numpy as np
from import_virtual_walker import HERE
import native_kernel as K
from spiking_sensory import BaseController,SpikingSensoryController
from run_virtual_walker import run
from assess_trial import assess
p=argparse.ArgumentParser();p.add_argument('--output-root',default=r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking');p.add_argument('--population-size',type=int,default=32);p.add_argument('--duration',type=float,default=10.);p.add_argument('--tag',default='sensory_v1');p.add_argument('--case',choices=['both','graded','hybrid'],default='both');p.add_argument('--basin-seeds',type=int,default=0);p.add_argument('--ablate-sensory-outputs',action='store_true');a=p.parse_args();out=Path(a.output_root);out.mkdir(parents=True,exist_ok=True)
cfg=json.loads((HERE/'GENTLE_HIP_0P6_FLEX_0P5.json').read_text());expected=cfg.pop('expected_source_sha256');assert hashlib.sha256(Path(cfg['source']).read_bytes()).hexdigest()==expected
rows=[];jobs=[('hybrid',s) for s in range(a.basin_seeds)] if a.basin_seeds else [(m,None) for m in (['graded','hybrid'] if a.case=='both' else [a.case])]
for mode,seed in jobs:
 candidate=dict(cfg);encoder=dict(rate_max=200.,tau_m=.005,refractory_s=.001,tau_syn_scale=1.);factors={}
 if seed is not None:
  rng=np.random.default_rng(100031+seed)
  for key in ('pf_mn_recruitment_scale','toe_mn_gmax','hip_pf_mn_scale','hip_flexor_pf_mn_scale'):
   factors[key]=float(rng.uniform(.99,1.01));candidate[key]*=factors[key]
  for key in encoder:factors[key]=float(rng.uniform(.99,1.01));encoder[key]*=factors[key]
 tag=a.tag+('' if seed is None else '_s'+str(seed))+'_'+mode;folder=out/'W2L_VirtualWalker_Newer'/tag
 fingerprint=hashlib.sha256(json.dumps(dict(mode=mode,source=expected,config=candidate,encoder=encoder,population=a.population_size,duration=a.duration,ablate_outputs=a.ablate_sensory_outputs),sort_keys=True).encode()).hexdigest()
 if (folder/'SPIKING_COMPARISON.json').exists():
  prior=json.loads((folder/'SPIKING_COMPARISON.json').read_text());assert prior.get('experiment_fingerprint')==fingerprint,'Existing result differs or predates fingerprint; choose a new tag';rows.append(prior);continue
 K.NativeController=BaseController if mode=='graded' else SpikingSensoryController;SpikingSensoryController.population_size=a.population_size;SpikingSensoryController.ablate_outputs=a.ablate_sensory_outputs
 for key,value in encoder.items():setattr(SpikingSensoryController,key,value)
 with (out/(a.tag+'.log')).open('a',encoding='utf-8') as log,contextlib.redirect_stdout(log):r=run(newer=True,duration=a.duration,tag=tag,output_root=out,**candidate);gate=assess(folder)
 row=dict(mode=mode,experiment_fingerprint=fingerprint,seed=seed,perturbation_factors=factors,source_sha256=expected,result=r,assessment=gate,config=candidate,encoder_parameters=encoder,source_logical_channels=108,source_chemical_edges=216,scope='Graded source comparator vs explicitly sensory-spiking hybrid; RG/PF/IN/MN remain source graded. Sixteen input-only sensory relay cells replaced by LIF populations; no shadow graded encoder drives their outputs.')
 if mode=='hybrid':
  c=SpikingSensoryController.last_instance;c.save(folder);counts=np.array([x[1] for x in c.logs]);row.update(spiking_populations=16,cells_per_population=a.population_size,actual_spiking_cells=16*a.population_size,remaining_graded_cells=92,actual_total_cells=92+16*a.population_size,spike_count_per_population=counts.sum(axis=0).tolist(),spiking_names=[c.names[i] for i in c.encoded],rate_max_hz=c.rate_max,source_matched_synaptic_filter_tau_s=c.tau_syn.tolist(),lif_rest_reset_mV=-70.,lif_threshold_mV=-50.,lif_tau_s=c.tau_m,refractory_s=c.refractory_s,ablation_edges=getattr(c,'ablated_edges',[]))
 else:
  reference=HERE/'runs/W2L_VirtualWalker_Newer/canonical108_ground_hip0p6_flex0p5_fresh_v1/trajectory.npz'
  if a.duration==10.:
   x=np.load(reference);y=np.load(folder/'trajectory.npz');row['existing_graded_trajectory_byte_parity']={k:x[k].tobytes()==y[k].tobytes() for k in x.files if x[k].dtype.kind in 'fiu'};assert all(row['existing_graded_trajectory_byte_parity'].values())
 (folder/'SPIKING_COMPARISON.json').write_text(json.dumps(row,indent=2));rows.append(row);(out/(a.tag+'_SUMMARY.json')).write_text(json.dumps(rows,indent=2));print(json.dumps(dict(case=mode,gate=gate['alternating_supported_stepping_gate'],swings={s:gate['legs'][s]['unloaded_clear_swings'] for s in 'LR'},alternation=gate['alternating_transition_fraction'],net_assistance_BW=r['mean_net_vertical_assistance_BW'],sim_per_wall=r['simulated_s_per_wall_s'],failure=r['failure']),indent=2),flush=True)
K.NativeController=BaseController
