"""Matched ankle-crossing sensory-source diagnostics with unchanged mechanics."""
import argparse,json,subprocess,sys
from pathlib import Path
import numpy as np
from scipy.signal import periodogram
from dominant_cycles import measure
HERE=Path(__file__).resolve().parent
def main():
 p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--preserve-knee-crossing',action='store_true');p.add_argument('--expression',type=float,default=.2913677891057415);p.add_argument('--fine-preserved',action='store_true');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);rows=[]
 cases=[('baseline',1.,1.),('Ia01',.1,1.),('Ia03',.3,1.),('II01',1.,.1),('II03',1.,.3),('both01',.1,.1)]
 if a.preserve_knee_crossing:cases=[('baseline',1.,1.),('both01',.1,.1),('both003',.03,.03),('Ia01',.1,1.)]
 if a.fine_preserved:
  if not a.preserve_knee_crossing:raise ValueError('Fine ladder requires preserved knee-crossing afferents')
  cases=[('both01',.1,.1),('both015',.15,.15),('both02',.2,.2)]
 (a.output/'MANIFEST.json').write_text(json.dumps(dict(fit=str(a.fit),cases=cases,expression=a.expression,preserve_knee_crossing=a.preserve_knee_crossing,global_local_II_scale=.2114065546411311,scope='Shared L/R local outputs from ankle-crossing muscle afferents; knee-crossing sources excluded when preserve_knee_crossing is true; RG/PF feedback and mechanical laws unchanged',quality_screen='Report fraction of ankle-angle spectral power above3Hz; flag>0.15 as engineering waveform-quality failure, in addition to movement/cadence gates'),indent=2))
 for label,ia,ii in cases:
  folder=a.output/label
  with(a.output/(label+'.log')).open('w',encoding='utf-8')as log:
   job=subprocess.run([sys.executable,str(HERE/'body_probe.py'),'--reflex-fit',str(a.fit),'--stage','2','--expression',str(a.expression),'--local-ii-scale','.2114065546411311','--ankle-ia-scale',str(ia),'--ankle-ii-scale',str(ii),'--audit-electrical','--output',str(folder)]+(['--preserve-knee-crossing']if a.preserve_knee_crossing else[]),stdout=log,stderr=subprocess.STDOUT,timeout=240,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
  row=dict(case=label,ankle_Ia_scale=ia,ankle_II_scale=ii,exit_code=job.returncode,folder=str(folder))
  if(folder/'body.result.json').exists():
   r=json.loads((folder/'body.result.json').read_text());z=np.load(folder/'capture/aligned_trace.npz');t=z['t_state'];keep=(t>=5.)&(t<14.);names=z['joint_names'].tolist();q=z['qpos'][:,z['joint_qpos_addresses']];principal=measure(t[keep],-np.degrees(q[keep][:,[names.index('knee_angle_'+s)for s in('r','l')]]),np.degrees(q[keep][:,[names.index('hip_flexion_'+s)for s in('r','l')]]));spectra=[]
   for side in('r','l'):
    f,power=periodogram(q[keep,names.index('ankle_angle_'+side)],fs=1/np.median(np.diff(t)),window='hann');spectra.append(dict(side=side,dominant_Hz=float(f[1:][np.argmax(power[1:])]),power_above3Hz_fraction=float(power[f>3].sum()/max(power[1:].sum(),1e-30))))
   audit=json.loads((folder/'ELECTRICAL_AUDIT.json').read_text());row.update(principal=principal,legacy=r['assessment'],ankle_spectra=spectra,measured=r['measured'],maximum_live_Euler_load_ratio=max(audit.get('maximum_Euler_load_ratio')or[0]),waveform_quality_pass=all(v['power_above3Hz_fraction']<=.15 for v in spectra))
  elif(folder/'STOP.json').exists():row['stop']=json.loads((folder/'STOP.json').read_text())
  rows.append(row);(a.output/'STATUS.json').write_text(json.dumps(rows,indent=2))
 (a.output/'COMPLETE.json').write_text(json.dumps(rows,indent=2))
if __name__=='__main__':main()
