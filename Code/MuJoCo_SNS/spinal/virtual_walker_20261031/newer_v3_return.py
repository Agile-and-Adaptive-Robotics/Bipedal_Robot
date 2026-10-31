from source_inputs import historical_source
"""Bounded deafferented→afferented→ground return of explicit V3 ablation."""
import json
import numpy as np
from scipy.signal import find_peaks
from import_virtual_walker import HERE,REPO
from run_virtual_walker import run
from assess_trial import assess
source=historical_source('W2L_VirtualWalker_Newer_GUI1nA');base=HERE/'runs/W2L_VirtualWalker_Newer_GUI1nA';rows=[]
def full_cycles(t,y):
 lo=float(y.min());span=float(np.ptp(y));high=lo+.65*span;low=lo+.25*span;armed=False;events=[]
 if span<10:return np.array(events)
 for ti,yi in zip(t,y):
  if yi<=low:armed=True
  elif armed and yi>=high:events.append(ti);armed=False
 return np.array(events)
cases=[('air_deafferented','air',False,0.,0.,None,8.),('air_afferented','air',True,0.,0.,None,8.),('ground_guides','supported',True,0.,0.,.03,10.),('source_default_regression','supported',True,None,1.,None,10.)]
for label,mode,feedback,tonic,v3,clearance,duration in cases:
 tag='newer_v3_diag_v1_'+label;folder=base/tag;params=dict(symmetric_rg_ext_tonic=tonic,v3_output_scale=v3,neural_time_scale=2.,guide_clearance=clearance,sagittal_harness=0.,rail_height_offset=-.02,yaw_perturb_nm=0.)
 result=run(source=source,mode=mode,feedback=feedback,duration=duration,tag=tag,**params);z=np.load(folder/'trajectory.npz');mask=z['t']>=1.;names=z['neuron_names'].tolist();metrics={};pk={}
 for name in ('L RG ext','L RG flx','R RG ext','R RG flx'):
  v=z['neuro'][mask,names.index(name)];p,_=find_peaks(v,prominence=1,distance=round(.15/result['physics_dt']));times=z['t'][mask][p];pk[name]=times;metrics[name]=dict(range_mV=float(np.ptp(v)),prominent_peaks=len(p),median_period_s=float(np.median(np.diff(times))) if len(times)>1 else None)
 phase=None;period=metrics['L RG ext']['median_period_s']
 if period and len(pk['L RG ext']) and len(pk['R RG ext']):phase=float(np.median([((t-pk['L RG ext'][np.argmin(abs(pk['L RG ext']-t))])/period)%1 for t in pk['R RG ext']]))
 report=dict(result=result,rg=metrics,RG_ext_LR_phase_fraction=phase)
 if mode=='air':
  cycles=[full_cycles(z['t'][mask],z['angles'][mask,k]) for k in (1,4)];period=float(np.median(np.diff(cycles[0]))) if len(cycles[0])>1 else None;cyclephase=float(np.median([((t-cycles[0][np.argmin(abs(cycles[0]-t))])/period)%1 for t in cycles[1]])) if period and len(cycles[1]) else None;report.update(full_excursion_cycles=[len(c) for c in cycles],full_excursion_phase_fraction=cyclephase,air_gate=bool(result['failure'] is None and min(len(c) for c in cycles)>=4 and cyclephase is not None and .35<=cyclephase<=.65 and np.max(abs(z['foot_vertical']))==0))
 else:report['ground']=assess(folder)
 if label=='source_default_regression':
  old=np.load(base/'newer_guide_pair_v1_t2_cNone/trajectory.npz');diff={k:float(np.max(abs(old[k]-z[k]))) for k in ('q','neuro','forces','root','foot_force_world')};assert max(diff.values())==0;report['default1_bit_identical_max_abs_differences']=diff
 (folder/'V3_DIAGNOSTIC.json').write_text(json.dumps(report,indent=2));(folder/'PARAMETERS.json').write_text(json.dumps(params,indent=2));rows.append(dict(label=label,folder=str(folder),report=report));(HERE/'runs/NEWER_V3_RETURN_V1.json').write_text(json.dumps(rows,indent=2))
