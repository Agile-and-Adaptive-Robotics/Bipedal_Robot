"""One positive bilateral V3 scale candidate, no gain sweep or new model family."""
import json
import numpy as np
from scipy.signal import find_peaks
from source_inputs import historical_source,HERE
from run_virtual_walker import run
from assess_trial import assess
source=historical_source('W2L_VirtualWalker_Newer_GUI1nA');base=HERE/'runs/W2L_VirtualWalker_Newer_GUI1nA';rows=[]
params=dict(symmetric_rg_ext_tonic=0.,v3_output_scale=.01,neural_time_scale=2.,sagittal_harness=0.,rail_height_offset=-.02,yaw_perturb_nm=0.)
for label,mode,feedback,guide,duration in [('air_deafferented','air',False,None,8.),('air_afferented','air',True,None,8.),('ground','supported',True,.03,10.)]:
 tag='positive_v3_0p01_v1_'+label;folder=base/tag;r=run(source=source,mode=mode,feedback=feedback,guide_clearance=guide,duration=duration,tag=tag,**params);z=np.load(folder/'trajectory.npz');mask=z['t']>=1.;names=z['neuron_names'].tolist();rg={}
 for name in ['L RG ext','L RG flx','R RG ext','R RG flx']:
  v=z['neuro'][mask,names.index(name)];pk,_=find_peaks(v,prominence=1.,distance=round(.15/r['physics_dt']));rg[name]=dict(peaks=len(pk),range_mV=float(np.ptp(v)),period_s=float(np.median(np.diff(z['t'][mask][pk]))) if len(pk)>1 else None)
 report=dict(result=r,RG=rg,source_scope='Immutable GUI-authoritative96/204 historical snapshot; proposed calibration for canonical NEWER, not original80 and not yet silently written to canonical.')
 if mode=='air':
  events=[]
  for k in (1,4):
   y=z['angles'][mask,k];low=float(y.min())+.25*np.ptp(y);high=float(y.min())+.65*np.ptp(y);armed=False;ev=[]
   if np.ptp(y)>=10:
    for t,value in zip(z['t'][mask],y):
     if value<=low:armed=True
     elif armed and value>=high:ev.append(t);armed=False
   events.append(np.array(ev))
  period=float(np.median(np.diff(events[0]))) if len(events[0])>1 else None;phase=float(np.median([((t-events[0][np.argmin(abs(events[0]-t))])/period)%1 for t in events[1]])) if period and len(events[1]) else None;report.update(knee_full_cycles=[len(e) for e in events],knee_phase=phase,air_gate=bool(r['failure'] is None and min(len(e) for e in events)>=4 and phase is not None and .35<=phase<=.65 and np.max(abs(z['foot_vertical']))==0))
 else:report['ground']=assess(folder)
 (folder/'PARAMETERS.json').write_text(json.dumps(dict(contact_si=False,pf_knee_flex_scale=1.,ankle_pf_crossed=False,guide_clearance=guide,belt_motor_end_s=10.,**params),indent=2));(folder/'POSITIVE_V3_REPORT.json').write_text(json.dumps(report,indent=2));rows.append(dict(label=label,folder=str(folder),report=report));(HERE/'runs/POSITIVE_V3_0P01_V1.json').write_text(json.dumps(rows,indent=2))
