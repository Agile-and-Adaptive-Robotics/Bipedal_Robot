"""Gated neural parity, air, then supported-body tests of adaptive spiking RG/PF."""
import os,sys,argparse,json,hashlib,contextlib
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
from scipy.signal import find_peaks
from import_virtual_walker import HERE
import native_kernel as K
from spiking_central_runtime import SpikingCentralController
from recruitment import apply_recruitment
from run_virtual_walker import run
from assess_trial import assess
p=argparse.ArgumentParser();p.add_argument('--mode',choices=['neural','air','supported','drop'],required=True);p.add_argument('--feedback-off',action='store_true');p.add_argument('--duration',type=float,default=10.);p.add_argument('--tag',required=True);p.add_argument('--out',default=r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking');p.add_argument('--configuration',default=r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\CENTRAL_NEURAL_WINNER.json');p.add_argument('--perturbation-seed',type=int);p.add_argument('--rail-height-offset',type=float);p.add_argument('--belt-motor-end-s',type=float);p.add_argument('--reference-calculations',action='store_true');p.add_argument('--sensory-spiking',action='store_true');p.add_argument('--audit-state',action='store_true');p.add_argument('--release-at',type=float,default=2.);p.add_argument('--assessment-start',type=float,default=1.);p.add_argument('--initial-qpos');p.add_argument('--pool-inherited-hip-physics',action='store_true');p.add_argument('--balance-configuration');p.add_argument('--belt-speed',type=float);p.add_argument('--clock-scale',type=float,default=1.);a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);config=json.loads(Path(a.configuration).read_text());SpikingCentralController.configuration=config;SpikingCentralController.reference_calculations=a.reference_calculations;SpikingCentralController.sensory_spiking=a.sensory_spiking;SpikingCentralController.audit_state=a.audit_state;assert not(a.sensory_spiking and (a.mode=="neural" or a.reference_calculations)),"Combined encoders currently use explicit air/ground checks, not the central-only neural reference"
SpikingCentralController.clock_scale=a.clock_scale
assert a.mode!='neural' or a.clock_scale==1.,'Neural parity mode is pinned to source clock'
SpikingCentralController.balance_configuration=None if a.balance_configuration is None else json.loads(Path(a.balance_configuration).read_text())
cfg=json.loads((HERE/'GENTLE_HIP_0P6_FLEX_0P5.json').read_text());expected=cfg.pop('expected_source_sha256');assert hashlib.sha256(Path(cfg['source']).read_bytes()).hexdigest()==expected
if a.pool_inherited_hip_physics:cfg['pool_inherited_hip_physics']=True
if a.initial_qpos is not None:
 initial=np.load(a.initial_qpos)['qpos'];assert initial.shape==(17,) and np.isfinite(initial).all();cfg['initial_qpos']=initial.tolist()
if a.rail_height_offset is not None:cfg['rail_height_offset']=a.rail_height_offset
if a.belt_motor_end_s is not None:cfg['belt_motor_end_s']=a.belt_motor_end_s
if a.belt_speed is not None:cfg['belt_speed_m_s']=a.belt_speed
if a.mode=='drop':cfg['release_at']=a.release_at
factors={}
if a.perturbation_seed is not None:
 from hip_knee_interlock import OPTIONAL_KEYS
 assert a.mode!='neural';rng=np.random.default_rng(a.perturbation_seed)
 for key in config:
  if key=='rg_reset_mV' or key in OPTIONAL_KEYS:continue
  factors[key]=float(rng.uniform(.99,1.01));config[key]*=factors[key]
 for key in ('pf_mn_recruitment_scale','toe_mn_gmax','hip_pf_mn_scale','hip_flexor_pf_mn_scale'):factors[key]=float(rng.uniform(.99,1.01));cfg[key]*=factors[key]
 if 'rg_reset_mV' in config:factors['rg_reset_mV']=float(rng.uniform(.99,1.01));config['rg_reset_mV']*=factors['rg_reset_mV']
 for key in OPTIONAL_KEYS:
  if key in config:factors[key]=float(rng.uniform(.99,1.01));config[key]*=factors[key]
 if SpikingCentralController.balance_configuration is not None and SpikingCentralController.balance_configuration.get('pathway_kind')=='contralateral_load_termination':
  for key in ('termination_gmax_uS','pulse_s','on_load_BW'):
   factors['pathway.'+key]=float(rng.uniform(.99,1.01));SpikingCentralController.balance_configuration[key]*=factors['pathway.'+key]
 if SpikingCentralController.balance_configuration is not None and SpikingCentralController.balance_configuration.get('pathway_kind')=='phase_stance_recruitment':
  from stance_recruitment import StanceRecruitmentRelay
  for key in StanceRecruitmentRelay.learned_keys(SpikingCentralController.balance_configuration):
   factors['pathway.'+key]=float(rng.uniform(.99,1.01));SpikingCentralController.balance_configuration[key]*=factors['pathway.'+key]
if a.mode=='neural':
 folder=out/a.tag;folder.mkdir(parents=True,exist_ok=True);assert not (folder/'NEURAL_PARITY.json').exists();bank=json.loads((out/'W2L_VirtualWalker_Newer/IMPORT_AUDIT.json').read_text())['muscles'];c=SpikingCentralController(cfg['source'],bank,ankle_pf_crossed=cfg['ankle_pf_crossed']);apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale']);c.apply_output_gain();trace=[];events=[]
 for tick in range(round(a.duration/c.dt)):
  t=tick*c.dt;c.applied[:]=0.
  for i,start,stop,amplitude in c.stimuli:
   if start<=t<stop:c.applied[i]+=amplitude
  c.tick(c.v,c.m,c.h,c.newv,c.gs,c.gi,c.lv,c.lm,c.lh,c.ref,c.adapt,c.decoded,c.counts)
  if tick%10==9:trace.append(c.v.copy());events.append(c.counts.copy());c.counts[:]=0
 trace=np.array(trace);events=np.array(events);reference=np.load(out/'central_pf_phase_v1/BEST_NEURAL_TRACE.npz');diff=float(np.max(abs(trace-reference['v']))) if trace.shape==reference['v'].shape else None;parity=dict(source_sha256=expected,configuration=config,max_voltage_difference_mV=diff,spike_count_byte_parity=bool(events.shape==reference['spike_counts'].shape and events.tobytes()==reference['spike_counts'].tobytes()),scope='One-tick runtime kernel vs fitted neural-only integrator; not mechanical validation.');np.savez_compressed(folder/'runtime_neural.npz',v=trace,spike_counts=events,names=c.names);(folder/'NEURAL_PARITY.json').write_text(json.dumps(parity,indent=2));print(json.dumps(parity,indent=2));assert diff is not None and diff<1e-10 and parity['spike_count_byte_parity'];sys.exit(0)
folder=out/'W2L_VirtualWalker_Newer'/a.tag
if (folder/'SPIKING_COMPARISON.json').exists():raise FileExistsError('Preserving completed run; use a new tag')
if a.mode=='air':
 assert len(cfg['initial_qpos'])==17
 if a.initial_qpos is not None:cfg['fixed_root_pose']=cfg['initial_qpos'][:7]
 cfg['initial_qpos']=cfg['initial_qpos'][7:]
K.NativeController=SpikingCentralController
try:
 with (out/(a.tag+'.log')).open('w',encoding='utf-8') as log,contextlib.redirect_stdout(log):r=run(newer=True,mode=a.mode,duration=a.duration,feedback=not a.feedback_off,reference=a.reference_calculations,tag=a.tag,output_root=out,**cfg);assessment=assess(folder,start_s=a.assessment_start)
finally:K.NativeController=SpikingCentralController.__bases__[0]
c=SpikingCentralController.last_instance;c.save(folder);r.update(controller_kind='AdaptiveCaLIF_RG_PF_hybrid',logical_source_channels=108,actual_spiking_cells=c.metadata()["actual_spiking_cells"],actual_graded_cells=c.metadata()["remaining_graded_cells"],spike_driven_source_edges=int(np.isin(c.edges[:,0].astype(int),np.r_[c.encoded,c.sensors.indices] if c.sensors is not None else c.encoded).sum()),legacy_census_field_scope='neurons and nonspiking_edges describe the imported logical source graph, not runtime cell types.',scope=c.metadata()['scope']);
if c.balance is not None:r.update(controller_kind=c.metadata()['architecture'],runtime_logical_channels=len(c.v),spike_driven_source_edges=r['spike_driven_source_edges']+int(np.isin(c.edges[:,0].astype(int),c.balance.indices).sum()),legacy_census_field_scope=f'Source108/216 preserved. Runtime proxy graph{len(c.v)}/{len(c.edges)} adds two actual LIF populations and{len(c.balance.edge_specs)} chemical edges; actual cell counts are separately reported.')
r['neural_clock_multiplier']=a.clock_scale;r['neural_physical_substep_s']=c.dt*a.clock_scale
(folder/'RESULT.json').write_text(json.dumps(r,indent=2));record=dict(initial_qpos_source=None if a.initial_qpos is None else dict(path=str(Path(a.initial_qpos).resolve()),sha256=hashlib.sha256(Path(a.initial_qpos).read_bytes()).hexdigest()),rail_release_at_s=a.release_at if a.mode=="drop" else None,reference_calculations=a.reference_calculations,reference_scope='Scalar Python neural and Hill equations; identical MuJoCo C engine, preallocation, and physics' if a.reference_calculations else 'Numba-compiled neural and Hill equations',source_sha256=expected,configuration=config,source_configuration=cfg,feedback=not a.feedback_off,mode=a.mode,result=r,assessment=assessment,spiking=c.metadata(),perturbation_seed=a.perturbation_seed,perturbation_factors=factors,implementation_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ('spiking_central_runtime.py','run_spiking_central.py','native_kernel.py','run_virtual_walker.py','recruitment.py','balance_pathway.py','spiking_sensory.py','hip_knee_interlock.py')})
record['balance_configuration_source']=None if a.balance_configuration is None else dict(path=str(Path(a.balance_configuration).resolve()),sha256=hashlib.sha256(Path(a.balance_configuration).read_bytes()).hexdigest())
if a.mode=='air':
 z=np.load(folder/'trajectory.npz');mask=z['t']>=1.;times={};legs={}
 for side,col in [('L',1),('R',4)]:
  peaks,_=find_peaks(z['angles'][mask,col],prominence=5,distance=round(.3/r['physics_dt']));times[side]=z['t'][mask][peaks];legs[side]=dict(knee_cycles=len(peaks),period_s=float(np.median(np.diff(times[side]))) if len(peaks)>1 else None,knee_range_deg=float(np.ptp(z['angles'][mask,col])))
 errors=[]
 for x,y in zip(times['L'][:-1],times['L'][1:]):
  v=times['R'][(times['R']>=x)&(times['R']<y)]
  if len(v):errors.append(abs(float((v[0]-x)/(y-x))-.5))
 phase=float(np.mean(errors)) if errors else .5;passed=r['failure'] is None and all(v['knee_cycles']>=3 and v['period_s'] is not None and .7<=v['period_s']<=1.4 for v in legs.values()) and phase<=.15
 record['air_gate']=dict(passed=passed,legs=legs,knee_peak_antiphase_error=phase,scope='Fixed elevated pelvis, all collisions disabled, same joint pose after removing free-root coordinates. Ground load/clearance gate intentionally not applied.')
(folder/'SPIKING_COMPARISON.json').write_text(json.dumps(record,indent=2));print(json.dumps(dict(mode=a.mode,air_gate=record.get('air_gate'),ground_gate=assessment['alternating_supported_stepping_gate'],failure=r['failure'],spiking_cells=c.metadata()['actual_spiking_cells']),indent=2))
