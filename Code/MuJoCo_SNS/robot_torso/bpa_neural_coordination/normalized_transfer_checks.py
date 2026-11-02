"""Restartable matched-transfer gates; one owned physics worker at a time."""
import argparse,json,subprocess,sys
from pathlib import Path
import numpy as np
from normalized_prepared_transfer import OUT
HERE=Path(__file__).resolve().parent

def equal(a,b):
 with np.load(a) as x,np.load(b) as y:
  if set(x.files)!=set(y.files):raise ValueError('State field mismatch: '+str(set(x.files)^set(y.files)))
  return {k:bool(np.array_equal(x[k],y[k])) for k in x.files}

def audit():
 report={}
 for neural in ('graded','spike'):
  base=OUT/('w2lvar_bpa_'+neural+'_v1')
  for suffix_case in ('fresh','replay'):
   other=OUT/('w2lvar_bpa_'+neural+'_'+suffix_case+'_v1')
   if not other.with_suffix('.comparison.json').exists():continue
   result={}
   for suffix in ('.release_state.npz','.npz','.apparatus.npz','.apparatus_extended.npz','.motor_state.npz','.pressure.npz','.inputs.npz'):
    result[suffix]=equal(base.with_suffix(suffix),other.with_suffix(suffix))
   if neural=='spike':
    for suffix in ('.spikes.npz','.final_spiking_state.npz'):result[suffix]=equal(base.with_suffix(suffix),other.with_suffix(suffix))
   if not all(v for g in result.values() for v in g.values()):raise ValueError('Matched nominal repeat/replay failed')
   report[neural+'_'+suffix_case]=dict(all_equal=True,comparisons=result)
  for kind in ('harsh','push'):
   for replay in ('','_replay'):
    other=OUT/('w2lvar_'+kind+'_'+neural+replay+'_v1')
    if other.with_suffix('.comparison.json').exists():
     result=equal(base.with_suffix('.release_state.npz'),other.with_suffix('.release_state.npz'))
     if not all(result.values()):raise ValueError('Held-out release states do not match')
     report[kind+'_'+neural+replay]=dict(release_all_equal=True,comparisons=result)
 for neural in ('graded','spike'):
  base=OUT/('w2lvar_bio_'+neural+'_v1');other=OUT/('w2lvar_bio_'+neural+'_fresh_v1')
  if other.with_suffix('.comparison.json').exists():
   result={s:equal(base.with_suffix(s),other.with_suffix(s)) for s in ('.release_state.npz','.npz','.apparatus.npz','.motor_state.npz','.inputs.npz')}
   if neural=='spike':result['.final_spiking_state.npz']=equal(base.with_suffix('.final_spiking_state.npz'),other.with_suffix('.final_spiking_state.npz'))
   if not all(v for group in result.values() for v in group.values()):raise ValueError('Biological fresh state repeat failed')
   report['bio_'+neural+'_fresh']=dict(all_equal=True,comparisons=result)
 for neural in ('graded','spike'):
  base=OUT/('w2lvar_rot0p8_'+neural+'_v1');other=OUT/('w2lvar_rot0p8_'+neural+'_fresh_v1')
  if other.with_suffix('.comparison.json').exists():
   suffixes=['.release_state.npz','.npz','.apparatus.npz','.apparatus_extended.npz','.motor_state.npz','.pressure.npz','.inputs.npz']
   if neural=='spike':suffixes+=['.spikes.npz','.final_spiking_state.npz']
   result={s:equal(base.with_suffix(s),other.with_suffix(s)) for s in suffixes}
   if not all(v for group in result.values() for v in group.values()):raise ValueError('Reduced-guidance fresh repeat failed')
   report['rot0p8_'+neural+'_fresh']=dict(all_equal=True,comparisons=result)
 for neural in ('graded','spike'):
  base=OUT/('w2lvar_bpa_'+neural+'_v1');other=OUT/('w2lvar_hold20_'+neural+'_v1')
  if other.with_suffix('.comparison.json').exists():
   release=equal(base.with_suffix('.release_state.npz'),other.with_suffix('.release_state.npz'));prefix={}
   for suffix in ('.motor_state.npz','.pressure.npz','.inputs.npz','.apparatus.npz'):
    with np.load(base.with_suffix(suffix)) as a,np.load(other.with_suffix(suffix)) as b:
     for key in a.files:
      if key in ('time','activation','excitation','engine_force_cache_N','samples','inputs','neuro'):prefix[suffix+':'+key]=bool(np.array_equal(a[key],b[key][:len(a[key])]))
   if neural=='spike':
    with np.load(base.with_suffix('.spikes.npz')) as a,np.load(other.with_suffix('.spikes.npz')) as b:prefix['spike_counts']=bool(np.array_equal(a['counts'],b['counts'][:len(a['counts'])]))
   if not all(release.values()) or not all(prefix.values()):raise ValueError('Long-hold common state or short-run prefix differs')
   report['hold20_'+neural]=dict(release_all_equal=True,short_run_prefix_exact=True,release=release,prefix=prefix)
 (OUT/'COMMON_STATE_AND_REPEAT_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report

def main(phase):
 if phase=='audit':print(json.dumps(audit(),indent=2));return
 if phase=='fresh':cases=('bio_graded_fresh','bio_spike_fresh','bpa_graded_fresh','bpa_spike_fresh','bpa_graded_replay','bpa_spike_replay')
 elif phase=='heldout':
  gates=audit()
  if not all(k in gates for k in ('graded_fresh','spike_fresh','graded_replay','spike_replay')):raise ValueError('Complete nominal repeat/replay gates first')
  cases=('harsh_graded','harsh_graded_replay','harsh_spike','harsh_spike_replay')
 elif phase=='push':
  gates=audit()
  if not all(k in gates for k in ('graded_fresh','spike_fresh','graded_replay','spike_replay')):raise ValueError('Complete nominal gates required')
  cases=('push_graded','push_graded_replay','push_spike','push_spike_replay')
 elif phase=='guidance':
  gates=audit()
  if not all(k in gates for k in ('graded_fresh','spike_fresh','graded_replay','spike_replay')):raise ValueError('Complete nominal gates required')
  cases=('rot0p8_graded','rot0p8_spike')
 elif phase=='guidance_validate':
  for neural in ('graded','spike'):
   base=OUT/('w2lvar_rot0p8_'+neural+'_v1');target=json.loads(base.with_suffix('.comparison.json').read_text());r=json.loads(Path(target['result']).read_text());nom=json.loads((OUT/('w2lvar_bpa_'+neural+'_v1.result.json')).read_text())
   if not target['assessment']['accepted']:raise ValueError('Reduced-guidance baseline failed')
   if not np.array_equal(r['standing_activation_targets'],nom['standing_activation_targets']):raise ValueError('Guide changed static neural targets')
  cases=('rot0p8_graded_fresh','rot0p8_spike_fresh','rot0p8_hold20_graded','rot0p8_hold20_spike')
 elif phase=='hold20':
  gates=audit()
  if not all(k in gates for k in ('graded_fresh','spike_fresh','graded_replay','spike_replay')):raise ValueError('Complete nominal gates required')
  cases=('hold20_graded','hold20_spike')
 else:raise ValueError(phase)
 for case in cases:
  target=OUT/('w2lvar_'+case+'_v1.comparison.json')
  if not target.exists():
   with (OUT/(case+'.log')).open('w',encoding='utf-8') as log:r=subprocess.run([sys.executable,str(HERE/'normalized_prepared_transfer.py'),'--case',case],stdout=log,stderr=subprocess.STDOUT)
   if r.returncode:raise RuntimeError(case)
  row=json.loads(target.read_text());print(json.dumps(dict(case=case,assessment=row['assessment'])),flush=True)
 audit()

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--phase',required=True);main(p.parse_args().phase)
