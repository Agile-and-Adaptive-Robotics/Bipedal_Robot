"""Restartable single-worker frozen comparison; no parameter tuning."""
import sys,subprocess,json
from pathlib import Path
from run_prepared_spiking import ROOT,VERSION
HERE=Path(__file__).resolve().parent
cases=('graded_bpa_current','bpa_spike_fresh','bpa_spike_replay','bio_spike','bio_spike_fresh','harsh_spike','harsh_spike_replay')
out=ROOT/'prepared_spiking'
for case in cases:
 p=out/('syn6_'+case+'_'+VERSION+'.comparison.json')
 if p.exists():print('Reuse '+case,flush=True);continue
 with (out/(case+'_'+VERSION+'.log')).open('w',encoding='utf-8') as log:
  completed=subprocess.run([sys.executable,str(HERE/'run_prepared_spiking.py'),'--case',case],stdout=log,stderr=subprocess.STDOUT)
 if completed.returncode:raise RuntimeError(case+' execution failed; inspect preserved log')
 result=json.loads(p.read_text());print(json.dumps(dict(case=case,assessment=result['assessment'])),flush=True)
