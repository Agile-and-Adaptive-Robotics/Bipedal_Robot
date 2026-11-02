"""Confirm two held-out live-feedback successes in separate Python processes."""
import json,subprocess,sys
from pathlib import Path
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')
HERE=Path(__file__).resolve().parent

def main():
 cases=[('maperror','pressure_forward_recheck','maperror_live_equilibrium_v3','maperror_fresh_equilibrium_v3'),('delay','pressure_series_recheck','delay_live_series_v2','delay_fresh_series_v2')]
 for label,folder,source,tag in cases:
  out=ROOT/folder;result=out/('syn6_'+tag+'.result.json')
  if not result.exists():
   command=[sys.executable,str(HERE/'gait_standing_family.py'),'--family','syn6','--out',str(out),'--tag',tag,'--parameters-json',str(ROOT/'rustuna_syn6_v1/BEST.json'),'--bpa-configuration',str(ROOT/('virtual_bpa_design/HELDOUT_'+label.upper()+'.json'))]
   with (ROOT/(tag+'.log')).open('w',encoding='utf-8') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
 import numpy as np
 from gait_standing_objective import score
 results=[]
 for label,folder,source,tag in cases:
  base=ROOT/folder/('syn6_'+source);fresh=ROOT/folder/('syn6_'+tag)
  r=json.loads(fresh.with_suffix('.result.json').read_text());assessment=score(r,fresh)[1];comparisons={}
  for suffix,keys in (('.npz',('qfull','ctrl','neuro')),('.apparatus.npz',('samples','neuro')),('.pressure.npz',('samples',))):
   a=np.load(base.with_suffix(suffix));b=np.load(fresh.with_suffix(suffix))
   for key in keys:
    x,y=a[key],b[key];comparisons[suffix+':'+key]=dict(bit_exact=x.shape==y.shape and np.array_equal(x,y),max_abs_difference=float(np.max(abs(x-y))) if x.shape==y.shape else None)
  results.append(dict(case=label,source_result=str(base.with_suffix('.result.json')),fresh_result=str(fresh.with_suffix('.result.json')),assessment=assessment,fresh_bit_exact=all(x['bit_exact'] for x in comparisons.values()),comparisons=comparisons))
 report=dict(machine='easteregg2',results=results,scope='Fresh repeats of fixed held-out virtual BPA standing cases; no new tuning, hardware identification, free-standing, or gait claim.')
 (ROOT/'BPA_HELDOUT_FRESH_CHECKS.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
