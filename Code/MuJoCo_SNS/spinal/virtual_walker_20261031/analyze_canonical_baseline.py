"""Declare stable baseline cycles and phase-matched stimulation before pulse runs."""
import json
import numpy as np
from scipy.signal import find_peaks
from import_virtual_walker import HERE
base=HERE/'runs/W2L_VirtualWalker_Newer/canonical108_air_baseline_v2';z=np.load(base/'trajectory.npz');r=json.loads((base/'RESULT.json').read_text());names=z['neuron_names'].tolist();assert len(names)==108 and r['failure'] is None;t=z['t'];dt=r['physics_dt'];mask=t>=1.;rg=z['neuro'][:,names.index('L RG ext')];peaks,_=find_peaks(rg,prominence=1.,distance=round(.2/dt));peaks=peaks[t[peaks]>=1.];stable=None
for k in range(len(peaks)-4):
 p=peaks[k:k+5];periods=np.diff(t[p[:4]]);amps=np.array([[np.ptp(z['angles'][p[j]:p[j+1],col]) for col in (1,4)] for j in range(3)]);cvp=float(np.std(periods)/np.mean(periods));cva=np.std(amps,axis=0)/np.mean(amps,axis=0)
 if cvp<=.05 and np.max(cva)<=.10 and np.min(amps)>=10:
  region=np.arange(p[3],p[4]);targets=['L RG ext','L Hip PF ext','L Hip MN ext'];activity=[]
  for n in targets:
   v=z['neuro'][mask,names.index(n)];lo=float(np.percentile(v,5));hi=float(np.percentile(v,95));activity.append((z['neuro'][region,names.index(n)]-lo)/(hi-lo))
  score=np.min(np.array(activity),axis=0);idx=region[int(np.argmax(score))];stable=dict(first_cycle_start_s=float(t[p[0]]),three_cycles_completed_s=float(t[p[3]]),periods_s=periods.tolist(),period_CV=cvp,knee_excursions_deg=amps.tolist(),knee_excursion_CV=cva.tolist(),pulse_t0_s=float(t[idx]),baseline_period_s=float(np.median(periods)),targets=targets,minimum_target_normalized_activity=float(np.max(score)),target_voltage_mV={n:float(z['neuro'][idx,names.index(n)]) for n in targets});break
metrics={}
for name in ('L RG ext','L RG flx','R RG ext','R RG flx'):
 v=z['neuro'][mask,names.index(name)];pk,_=find_peaks(v,prominence=1.,distance=round(.2/dt));metrics[name]=dict(peaks=len(pk),range_mV=float(np.ptp(v)),median_period_s=float(np.median(np.diff(t[mask][pk]))) if len(pk)>1 else None)
out=dict(source_result=str(base/'RESULT.json'),actual_neurons=len(names),source_sha256=r['source_sha256'],stability_rule='At least3 completed RG cycles, interval CV≤5%, bilateral knee excursion CV≤10%, each knee excursion≥10deg. Matched pulse onset maximizes common normalized activity of the3 chosen extensor targets during next baseline cycle.',stable=stable,RG=metrics,min_knee_flexion_deg=np.min(z['angles'][:,[1,4]],axis=0).tolist(),ground_contact_max_N=float(np.max(abs(z['foot_vertical']))),duration_s=r['duration_completed_s']);(base/'BASELINE_ANALYSIS.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
