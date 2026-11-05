"""Bounded short-stride sensory feedback trials on the slow candidate."""
import argparse,json,contextlib
from pathlib import Path
from run_virtual_walker import run
from shuffle_optimization_20261104 import measure
def main():
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--duty-ladder',action='store_true');p.add_argument('--heel-ladder',action='store_true');p.add_argument('--heel-ablation',action='store_true');a=p.parse_args();a.root.mkdir(parents=True,exist_ok=False);parent=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\W2L_short_slow_shuffle_20261104');best=json.loads((parent/'BEST.json').read_text())['trial'];cfg=best['diagnostics']['configuration'];rows=[]
 cases=[('baseline',None)]+[(f'limit_{lower}_{upper}_{gain}',dict(lower_deg=lower,upper_deg=upper,gain_nA_deg=gain))for lower,upper in((-8.,8.),(-5.,10.),(-5.,5.))for gain in(.1,.3,.6)]
 if a.duty_ladder:cases=[('baseline',None)]+[(f'duty_{duration}_{clock}_{gain}',None if gain==0 else dict(lower_deg=-8.,upper_deg=8.,gain_nA_deg=gain))for duration,clock in((1.5,1.8),(2.,1.6),(2.5,1.5))for gain in(0.,.3,.6)]
 if a.heel_ladder:cases=[('baseline',None)]+[(f'heel_{dwell}_{gain}',None)for dwell in(0.,.06,.12)for gain in(1.,2.,4.)]
 if a.heel_ablation:cases=[('baseline',None),('ablate_contra_heel',None),('static_0.25',None),('static_0.5',None)]
 (a.root/'MANIFEST.json').write_text(json.dumps(dict(source_trial=best,cases=cases,scope='Soft neural feedback through actual inhibitory hip PF relays; no mechanical limits or force clipping'),indent=2))
 for label,limit in cases:
  c=dict(cfg)
  if label.startswith('duty_'):
   fields=label.split('_');c['rg_ext_duration_scale']=float(fields[1]);c['neural_time_scale']=float(fields[2])
  if label.startswith('heel_'):
   fields=label.split('_');c['heel_handover']=dict(dwell_s=float(fields[1]),gain=float(fields[2]))
  if label=='ablate_contra_heel':c['heel_handover']=dict(mode='ablate')
  if label.startswith('static_'):c['heel_handover']=dict(mode='tonic_scale',gain=float(label.split('_')[1]))
  with(a.root/(label+'.log')).open('w',encoding='utf-8')as log,contextlib.redirect_stdout(log):
   run(newer=True,duration=16.,tag=label,output_root=a.root/'trials',hip_pf_limits=limit,**c);folder=a.root/'trials/W2L_VirtualWalker_Newer'/label;result=measure(folder)
  rows.append(dict(case=label,configuration=dict(c,hip_pf_limits=limit),folder=str(folder),**result));(a.root/'STATUS.json').write_text(json.dumps(rows,indent=2))
 (a.root/'COMPLETE.json').write_text(json.dumps(rows,indent=2));print(json.dumps(dict(trials=len(rows),passed=sum(r['passed']for r in rows))))
if __name__=='__main__':main()
