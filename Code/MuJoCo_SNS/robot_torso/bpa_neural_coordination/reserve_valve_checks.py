"""Restartable frozen winner repeat/ablation and withheld valve-feedback checks."""
import json,subprocess,sys,argparse
from pathlib import Path
import numpy as np
from reserve_valve_runner import DEST
HERE=Path(__file__).resolve().parent

def equal(a,b):
 with np.load(a) as x,np.load(b) as y:
  assert set(x.files)==set(y.files)
  return {k:bool(np.array_equal(x[k],y[k])) for k in x.files}

def audit():
 selected=json.loads((DEST/'SELECTED.json').read_text());base=Path(selected['result'].removesuffix('.result.json'));report={};fresh=DEST/'w2lvar_winner_fresh'
 if fresh.with_suffix('.comparison.json').exists():
  suffixes=['.release_state.npz','.npz','.motor_state.npz','.pressure.npz','.inputs.npz','.apparatus.npz','.apparatus_extended.npz','.relay.npz','.valve.npz','.final_neural_state.npz','.final_actuator_mechanics_state.npz'];checks={s:equal(base.with_suffix(s),fresh.with_suffix(s)) for s in suffixes}
  if not all(v for g in checks.values() for v in g.values()):raise ValueError('Fresh complete-state repeat differs')
  report['fresh_all_numeric_exact']=checks
 for condition in ('delay','map','pulse','harsh'):
  path=DEST/('w2lvar_heldout_'+condition)
  if path.with_suffix('.comparison.json').exists():
   checks=equal(base.with_suffix('.release_state.npz'),path.with_suffix('.release_state.npz'))
   if not all(checks.values()):raise ValueError('Withheld test common release differs: '+condition)
   report[condition+'_complete_release_exact']=checks
 zero=DEST/'w2lvar_zero_current'
 if zero.with_suffix('.comparison.json').exists():
  legacy=DEST.parent/'inverse_boundary_policy/w2lvar_nearest_graded_20s_v1';checks={s:equal(legacy.with_suffix(s),zero.with_suffix(s)) for s in ('.motor_state.npz','.pressure.npz','.inputs.npz','.apparatus.npz','.apparatus_extended.npz')}
  if not all(v for g in checks.values() for v in g.values()):raise ValueError('Latest zero-output physical parity failed')
  report['latest_zero_output_physical_exact']=checks
 (DEST/'COMPLETE_STATE_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report

def run_case(tag,cap,preview,reserve,condition='nominal',encoding='graded'):
 path=DEST/('w2lvar_'+tag+'.comparison.json')
 if not path.exists():
  with (DEST/(tag+'.log')).open('w',encoding='utf-8') as log:proc=subprocess.run([sys.executable,str(HERE/'reserve_valve_runner.py'),'--cap',str(cap),'--preview',str(preview),'--reserve-target',str(reserve),'--condition',condition,'--encoding',encoding,'--tag',tag],stdout=log,stderr=subprocess.STDOUT)
  if proc.returncode:raise RuntimeError(tag)
 c=json.loads(path.read_text());print(json.dumps(dict(tag=tag,assessment=c['assessment'])),flush=True);return c

def main(phase):
 selected=json.loads((DEST/'SELECTED.json').read_text());cap=selected['pressure_cap_kpa'];preview=selected['preview_s'];reserve=selected['reserve_target_N']
 if phase=='fresh':
  run_case('zero_current',0.,preview,reserve);fresh=run_case('winner_fresh',cap,preview,reserve)
  if not fresh['assessment']['accepted']:raise ValueError('Fresh winner failed')
  audit()
 elif phase=='heldout':
  gates=audit()
  if not all(k in gates for k in ('fresh_all_numeric_exact','latest_zero_output_physical_exact')):raise ValueError('Fresh/zero gates required')
  plan=dict(conditions=['delay','map','pulse','harsh'],duration_s=20.,frozen_selection=selected['result'],selection_does_not_change_after_failures=True,paired_controls='same model/condition with valve output gain zero')
  p=DEST/'PREDECLARED_HELDOUT.json'
  if p.exists() and json.loads(p.read_text())!=plan:raise ValueError('Heldout plan changed')
  p.write_text(json.dumps(plan,indent=2),encoding='utf-8')
  for condition in plan['conditions']:
   run_case('heldout_'+condition,cap,preview,reserve,condition);run_case('off_'+condition,0.,preview,reserve,condition)
  audit()
 elif phase=='spike':
  if 'fresh_all_numeric_exact' not in audit():raise ValueError('Fresh nominal gate required')
  run_case('sensory_spike_nominal',cap,preview,reserve,encoding='spike')
 else:audit()

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--phase',choices=['fresh','heldout','spike','audit'],required=True);main(p.parse_args().phase)
