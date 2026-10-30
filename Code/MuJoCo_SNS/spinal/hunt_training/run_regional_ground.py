"""Bounded regional-contact stage-4 diagnostics, accelerated and fresh-process."""
import sys,json,subprocess,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/tonic_abduction_20261030_s1_v1/0p02'
out=root/'regional_ground_v1'
if out.exists():raise SystemExit('Existing evidence retained')
out.mkdir();shutil.copyfile(Path(__file__),out/Path(__file__).name)
base=json.loads((root/'ground_initial/force_support70.json').read_text())
ib=root/'stance_support/ib_0p3/weights'
plans={
 'regional_no_cross':{},
 'regional_tonic_cross':dict(cross_contact_gain=6.),
 'regional_event_cross':dict(cross_contact_gain=6.,cross_contact_event_amplitude=2.),
 'regional_event_ib':dict(cross_contact_gain=6.,cross_contact_event_amplitude=2.,weights_directory=str(ib)),
 'regional_event_ib_timing_off':dict(cross_contact_gain=6.,cross_contact_event_amplitude=2.,weights_directory=str(ib),timing_afferents=False)}
def run(item):
    name,changes=item;cfg={**base,'regional_contacts':True,**changes};p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
    with p.with_suffix('.log').open('w') as log:
        code=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
    r=json.loads(p.with_suffix('.result.json').read_text()) if code==0 else dict(exit_code=code)
    return name,r
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
    for name,r in pool.map(run,plans.items()):
        rows[name]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(name,r.get('assessment',{}).get('accepted'),flush=True)
