"""Bounded positive-gain sensitivity test after fixed PF readout fitting.

No zero-conductance edges, new cells, or biomechanical changes. Distinguishes
inherited reflex magnitude from lack of feedforward representational capacity.
"""
import argparse, hashlib, json, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'virtual_walker_20261031'))
from rustuna_transfer import Worker

def main():
 p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--relay-outputs',action='store_true');a=p.parse_args()
 a.output.mkdir(parents=True,exist_ok=False)
 cases=[dict(local_Ia=g,local_Ib=g,local_II=g)for g in(.0001,.001,.01,.03)]
 cases += [dict(local_Ia=.001,local_Ib=g,local_II=.001)for g in(.1,1.)]
 if a.relay_outputs:
  cases=[dict(local_Ia=.1,local_Ib=.1,local_II=.01,relay_Ib_output=ib,relay_II_output=ii)for ib,ii in((1.,1.),(.1,1.),(1.,.1),(.1,.1),(.03,.03),(.01,.01),(.001,.001),(.1,.01))]
 preflight=dict(op='preflight',fit=str(a.fit),fit_sha256=hashlib.sha256((a.fit/'FINAL.json').read_bytes()).hexdigest())
 manifest=dict(**preflight,cases=cases,scope='Offline sensitivity diagnostic, mirrored positive gains; accepted requires original D waveform gate and C peak preservation',source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest()for f in(Path(__file__),HERE/'motor_adapter.py',HERE/'readout_keys.py',HERE/'replay_adapter.py',HERE/'pose_replay.py')})
 (a.output/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
 def guard(r):
  d=r['diagnostics'];d['waveform_gate_passed']=d['accepted'];d['peak_preservation_gate']=d['mn_peak_rmse']<=.25;d['accepted']=bool(d['accepted']and d['peak_preservation_gate']);return r
 w=Worker(HERE/'motor_adapter.py',a.output);rows=[]
 try:
  ready=w.request(preflight,180.);(a.output/'PREFLIGHT.json').write_text(json.dumps(ready,indent=2))
  if not ready['ready_for_optimization']:raise ValueError('C gate failed')
  for i,parameters in enumerate(cases):
   row=dict(trial=i,parameters=parameters,**guard(w.request(dict(op='evaluate',params=parameters),180.)));rows.append(row)
   (a.output/'STATUS.json').write_text(json.dumps(rows,indent=2))
 finally:w.close()
 best=min(rows,key=lambda r:(not r['diagnostics']['accepted'],r['value']))
 fresh=Worker(HERE/'motor_adapter.py',a.output)
 try:
  fresh.request(preflight,180.)
  repeat=guard(fresh.request(dict(op='capture',params=best['parameters'],path=str(a.output/'BEST_TRACE.npz')),180.))
 finally:fresh.close()
 matched=abs(repeat['value']-best['value'])<1e-10
 (a.output/'FRESH_REPEAT.json').write_text(json.dumps(dict(matched=matched,**repeat),indent=2))
 (a.output/'FINAL.json').write_text(json.dumps(dict(accepted=bool(best['diagnostics']['accepted']and repeat['diagnostics']['accepted']and matched),best=best,parameters=best['parameters'],fresh_process_matched=matched,scope=manifest['scope']),indent=2))
 print(json.dumps(dict(accepted=best['diagnostics']['accepted'],parameters=best['parameters'],diagnostics=best['diagnostics'])),flush=True)
if __name__=='__main__':main()
