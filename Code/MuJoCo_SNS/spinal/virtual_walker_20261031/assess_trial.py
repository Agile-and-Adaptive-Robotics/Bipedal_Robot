"""Separate oscillation, supported clearance, and unrestricted walking claims."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
def assess(folder):
 folder=Path(folder);z=np.load(folder/'trajectory.npz');r=json.loads((folder/'RESULT.json').read_text());t=z['t'];dt=r['physics_dt'];mask=t>=1.;angles=z['angles'];weight=41.919002*9.81;report={};events=[]
 for side,knee in [('L',1),('R',4)]:
  si=0 if side=='L' else 1;indices=np.flatnonzero(mask);peaks,_=find_peaks(angles[mask,knee],prominence=5,distance=round(.2/dt));peak_times=t[indices[peaks]]
  clear=np.min(z['foot_clearance'][:,2*si:2*si+2],axis=1);load=z['foot_vertical'][:,si]
  swing=mask&(clear>.001)&(load<.025*weight);edges=np.diff(np.r_[False,swing,False].astype(int));starts=np.flatnonzero(edges==1);stops=np.flatnonzero(edges==-1);qualified=[]
  for a,b in zip(starts,stops):
   other_stance=float(np.mean(z['foot_vertical'][a:b,1-si]>.05*weight))
   if (b-a)*dt>=.08 and clear[a:b].max()>=.01 and other_stance>=.5:qualified.append(dict(start_s=float(t[a]),end_s=float(t[min(b,len(t)-1)]),peak_clearance_m=float(clear[a:b].max()),other_stance_fraction=other_stance));events.append((float(t[a]),side))
  report[side]=dict(knee_peak_count=len(peaks),period_median_s=float(np.median(np.diff(peak_times))) if len(peaks)>1 else None,unloaded_clear_swings=len(qualified),swing_intervals=qualified,loaded_fraction=float(np.mean(load[mask]>.025*weight)),max_whole_foot_clearance_m=float(clear[mask].max()),min_anatomical_knee_flexion_deg=float(angles[:,knee].min()))
 events.sort();alternating=sum(a[1]!=b[1] for a,b in zip(events[:-1],events[1:]));minimum=min(report[x]['unloaded_clear_swings'] for x in report);fraction=alternating/max(1,len(events)-1)
 out=dict(source=str(folder.resolve()),evaluation_start_s=1.,criteria='Whole-foot clearance>1mm for80ms, peak≥10mm, load<.025BW; opposite foot bears>.05BW for≥half interval. ≥3 swings per leg and≥80% alternating transitions. Knee peaks separately reported.',legs=report,alternating_swing_transitions=alternating,alternating_transition_fraction=fraction,alternating_supported_stepping_gate=minimum>=3 and alternating>=5 and fraction>=.8 and r['failure'] is None,walking_claim=False,reason='Passive support rails bear substantial weight; this gate does not establish independent walking.')
 (folder/'ASSESSMENT.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('folder');a=p.parse_args();assess(a.folder)
