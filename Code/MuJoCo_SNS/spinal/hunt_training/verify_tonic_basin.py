"""Fresh-process basin around the shared two-stage tonic-PF candidate.

Checks each new tuning scalar at +/-1%, and two seeded simultaneous perturbations
of PF/BI/MN/reflex edge weights. Intrinsic and expected-afferent fits stay frozen.
"""
import os,sys,json,shutil,subprocess
from pathlib import Path
import numpy as np
H=Path(__file__).resolve().parent;candidate=H/'runs/tonic_abduction_20261030_s1_v1/0p02'
out=candidate/'basin'
if out.exists():raise SystemExit('Existing basin retained')
out.mkdir();summary={}
# The three saved controller scalars are the tonic increment, commissural gain,
# and biarticular pair gain. Probe both signs one scalar at a time.
cases=[(knob,sign) for knob in ('tonic','commissural','biarticular') for sign in (-1,1)]+[('edge_jitter',seed) for seed in (301030,301031)]
for knob,sign in cases:
    name=knob+'_'+str(sign).replace('-','minus');root=out/name;weights=root/'weights';weights.mkdir(parents=True)
    for filename in ('pf_weights.json','bi_weights.json','mn_weights.json','reflex_weights.json','timing_weights.json','circuit_options.json'):shutil.copyfile(candidate/'weights'/filename,weights/filename)
    config=dict(biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True)
    if knob=='tonic':
        path=weights/'mn_weights.json';z=json.loads(path.read_text())
        for change in z['experimental_tonic_abductor_support']['changes']:
            z['edges'][change['edge']]=change['original']+change['added_conductance']*(1+sign*.01)
        path.write_text(json.dumps(z,indent=2))
    elif knob=='edge_jitter':
        rng=np.random.default_rng(sign)
        for filename in ('pf_weights.json','bi_weights.json','mn_weights.json','reflex_weights.json'):
            path=weights/filename;z=json.loads(path.read_text());z['edges']={key:value*float(rng.uniform(.99,1.01)) for key,value in z['edges'].items()};path.write_text(json.dumps(z,indent=2))
    else:config[knob+'_gain' if knob=='commissural' else 'biarticular_pair_gain']*=1+sign*.01
    path=root/'config.json';path.write_text(json.dumps(config,indent=2))
    with (root/'pair.log').open('w') as log:
        code=subprocess.run([sys.executable,str(H/'check_shared_air.py'),'--weights',str(weights),'--config',str(path),'--output-directory',str(root/'air_pair')],stdout=log,stderr=subprocess.STDOUT).returncode
    summary[name]=json.loads((root/'air_pair/SUMMARY.json').read_text()) if code==0 else dict(accepted=False,exit=code)
    (out/'SUMMARY.json').write_text(json.dumps(summary,indent=2));print(name,summary[name]['accepted'],flush=True)
shutil.copyfile(Path(__file__),out/Path(__file__).name)
print('BASIN_COMPLETE',sum(r['accepted'] for r in summary.values()),len(summary),flush=True)
