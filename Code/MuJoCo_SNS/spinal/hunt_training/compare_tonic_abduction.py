"""Separate tonic hip-abductor support from global phasic actuator scaling.

Adds a small conductance from the existing gait-gated tonic PF cell to
hip-abductor MNs. No positional sensor, new oscillator, target clamp, or
body restraint is added. This differs from ZCode A's scale-8 feedback sweep.
"""
import os,sys,json,shutil,subprocess
from pathlib import Path
H=Path(__file__).resolve().parent;sys.path.insert(0,str(H.parent))
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import build_network_w2lvar as W
source=H/'runs/native_joint_pf_mono_20261030_v1/weights';out=H/'runs/tonic_abduction_20261030_s1_v1'
if out.exists():raise SystemExit('Existing experiment retained')
out.mkdir();summary={}
for increment in (.02,.04,.08):
    root=out/str(increment).replace('.','p');weights=root/'weights';weights.mkdir(parents=True)
    for name in ('pf_weights.json','bi_weights.json','mn_weights.json','reflex_weights.json','timing_weights.json','circuit_options.json'):shutil.copyfile(source/name,weights/name)
    value=json.loads((weights/'mn_weights.json').read_text());changes=[]
    for key,gain in value['edges'].items():
        pre,post=key.split(' -> ')
        if pre=='PF_HUNT_TONIC' and W.classify(post[3:]).groups[0]=='hip_abd':
            value['edges'][key]=gain+increment;changes.append(dict(edge=key,original=gain,added_conductance=increment))
    assert changes
    value['experimental_tonic_abductor_support']=dict(changes=changes,scope='Added conductance; not a guaranteed activation floor or a physiological coactivation index')
    (weights/'mn_weights.json').write_text(json.dumps(value,indent=2))
    with (root/'pair.log').open('w') as log:
        code=subprocess.run([sys.executable,str(H/'check_shared_air.py'),'--weights',str(weights),'--output-directory',str(root/'air_pair')],stdout=log,stderr=subprocess.STDOUT).returncode
    summary[str(increment)]=json.loads((root/'air_pair/SUMMARY.json').read_text()) if code==0 else dict(accepted=False,exit=code)
    (out/'SUMMARY.json').write_text(json.dumps(summary,indent=2));print(increment,summary[str(increment)]['accepted'],flush=True)
shutil.copyfile(Path(__file__),out/Path(__file__).name)
