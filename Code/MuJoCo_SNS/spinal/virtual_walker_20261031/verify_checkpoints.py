"""Read completed fresh-process trials; no hidden retuning or pooling of settings."""
import json
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
from import_virtual_walker import HERE
from assess_trial import assess
base=HERE/'runs/W2L_VirtualWalker';out={}
def full_cycles(t,y):
 lo=float(y.min());span=float(np.ptp(y));high=lo+.65*span;low=lo+.25*span;armed=False;events=[]
 if span<10:return np.array(events)
 for ti,yi in zip(t,y):
  if yi<=low:armed=True
  elif armed and yi>=high:events.append(ti);armed=False
 return np.array(events)
for tag in ['timing2_fresh_ground','timing2_air_deafferented','timing2_air_afferented','timing2_dt1ms']:
 folder=base/tag
 if not (folder/'RESULT.json').exists():continue
 r=json.loads((folder/'RESULT.json').read_text());z=np.load(folder/'trajectory.npz');t=z['t'];mask=t>=1.;dt=r['physics_dt'];metrics=dict(result=r)
 if r['mode']=='air':
  pk=[]
  for k in [1,4]:p,_=find_peaks(z['angles'][mask,k],prominence=5,distance=round(.25/dt));pk.append(t[mask][p])
  period=float(np.median(np.diff(pk[0]))) if len(pk[0])>1 else None;phases=[]
  if period:
   for q in pk[1]:phases.append(float(((q-pk[0][np.argmin(abs(pk[0]-q))])/period)%1))
  phase=float(np.median(phases)) if phases else None
  full=[full_cycles(t[mask],z['angles'][mask,k]) for k in [1,4]];full_period=float(np.median(np.diff(full[0]))) if len(full[0])>1 else None;full_phases=[]
  if full_period:
   for q in full[1]:full_phases.append(float(((q-full[0][np.argmin(abs(full[0]-q))])/full_period)%1))
  full_phase=float(np.median(full_phases)) if full_phases else None
  metrics.update(local_maxima_count=[len(x) for x in pk],local_peak_period_s=period,local_peak_phase_fraction=phase,full_excursion_cycles=[len(x) for x in full],full_excursion_period_s=full_period,full_excursion_phase_fraction=full_phase,measurement_note='Full excursion reaches65% of its recorded angle range then returnsbelow25%; range must exceed10deg. Local maxima retained separately because small recoil is not a complete flexioncycle.',air_gate=bool(r['failure'] is None and min(len(x) for x in full)>=4 and full_phase is not None and .35<=full_phase<=.65 and np.max(abs(z['foot_vertical']))==0))
 else:metrics['assessment']=assess(folder)
 out[tag]=metrics
if 'timing2_fresh_ground' in out:
 first=np.load(base/'timing_v1_scale_2/trajectory.npz');repeat=np.load(base/'timing2_fresh_ground/trajectory.npz');out['fresh_repeat_max_abs_diff']={k:float(np.max(abs(first[k]-repeat[k]))) for k in ['q','angles','neuro','forces','root','rail_vertical','foot_vertical','foot_clearance']}
(HERE/'runs/CHECKPOINT_VERIFICATION.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:({q:v[q] for q in ('air_gate','local_maxima_count','full_excursion_cycles','full_excursion_period_s','full_excursion_phase_fraction') if q in v} if isinstance(v,dict) else v) for k,v in out.items()},indent=2))
