"""Audit all serialized release states, not just visible graded neuron voltages."""
import json
from pathlib import Path
import numpy as np
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')
OUT=ROOT/'prepared_spiking'

def equal(a,b,common=False):
 with np.load(a) as za,np.load(b) as zb:
  ka=set(za.files);kb=set(zb.files)
  if not common and ka!=kb:raise ValueError('Snapshot fields differ: '+str(ka^kb))
  rows={k:bool(np.array_equal(za[k],zb[k])) for k in sorted(ka&kb)}
  return dict(all_equal=all(rows.values()),fields=len(rows),checks=rows,scope='Common pre-existing fields only' if common else 'All fields, including exact field-set identity')

def main():
 pairs=[('fresh','bpa_spike','bpa_spike_fresh'),('nominal_replay','bpa_spike','bpa_spike_replay'),('harsh_live','bpa_spike','harsh_spike'),('harsh_replay','bpa_spike','harsh_spike_replay'),('bio_fresh','bio_spike','bio_spike_fresh')];checks={}
 for label,a,b in pairs:
  pa=OUT/('syn6_'+a+'_v2');pb=OUT/('syn6_'+b+'_v2')
  if not pb.with_suffix('.comparison.json').exists():continue
  row={'release':equal(pa.with_suffix('.release_state.npz'),pb.with_suffix('.release_state.npz'))}
  if label in ('fresh','nominal_replay','bio_fresh'):
   for suffix in ('.npz','.apparatus.npz','.apparatus_extended.npz'):
    row[suffix]=equal(pa.with_suffix(suffix),pb.with_suffix(suffix))
   if label!='bio_fresh':row['.pressure.npz']=equal(pa.with_suffix('.pressure.npz'),pb.with_suffix('.pressure.npz'))
   ma=json.loads(pa.with_suffix('.comparison.json').read_text());mb=json.loads(pb.with_suffix('.comparison.json').read_text());row['final_spike_counts_equal']=ma['spiking']['spike_counts']==mb['spiking']['spike_counts']
  checks[label]=row
 graded=OUT/'syn6_graded_bpa_current_v2'
 if graded.with_suffix('.result.json').exists():
  original=ROOT/'pinned_preparation/syn6_prepared_nominal_v1';checks['default_graded_unchanged']={s:equal(original.with_suffix(s),graded.with_suffix(s),common=True) for s in ('.npz','.apparatus.npz','.pressure.npz','.release_state.npz')}
 nominal=OUT/'syn6_bpa_spike_v2.release_state.npz'
 with np.load(nominal) as z:
  fields=z.files;required=['transformed_encoder_'+k for k in ('v','ref','decoded','spikes','counts','total_counts','target','steps','logs','rng_state','names','indices','ports','parameters')]
  if not all(k in fields for k in required):raise ValueError('Missing required encoder state')
  checks['encoder_coverage']=dict(required_fields=required,all_present=True,explicit_empty_rng_state=z['transformed_encoder_rng_state'].size==0,spikes_during_preparation=int(z['transformed_encoder_total_counts'].sum()),fields=len(fields))
 result=dict(machine='easteregg2',checks=checks,scope='Equality applies within an encoder/actuator treatment under common preparation. Across biological/BPA or graded/spiking treatments, pose and preparation policy match but dynamics naturally differ. Source constants are fixed by saved configuration and canonical source ordering.')
 (OUT/'COMMON_STATE_AUDIT.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
