"""Test the available contralateral heel excitation and a truly absent task branch."""
import sys,json,subprocess,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/tonic_abduction_20261030_s1_v1/0p02';out=root/'contact_pathway'
if out.exists():raise SystemExit('Existing contact-pathway tests retained')
out.mkdir();ground=json.loads((root/'ground_initial/force_support70.json').read_text());spring=json.loads((root/'ground_initial/spring_support.json').read_text());stand=json.loads((root/'standing/task_low.json').read_text())
cases=[('task_physically_absent',{**stand,'balance':None}),
       *[(f'contra_heel_{gain:g}',{**ground,'cross_contact_gain':gain,'export_graph':True}) for gain in (1.,4.,12.)],
       ('spring_contra_heel12',{**spring,'cross_contact_gain':12.,'export_graph':True})]
def run(item):
    name,cfg=item;path=out/(name+'.json');path.write_text(json.dumps(cfg,indent=2))
    with path.with_suffix('.log').open('w') as log:
        code=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(path)],stdout=log,stderr=subprocess.STDOUT).returncode
    return name,json.loads(path.with_suffix('.result.json').read_text()) if code==0 else dict(accepted=False,exit=code)
with ThreadPoolExecutor(max_workers=2) as pool:
    rows={}
    for name,result in pool.map(run,cases):
        rows[name]=result;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(name,result.get('assessment',{}).get('accepted',False),flush=True)
shutil.copyfile(Path(__file__),out/Path(__file__).name)
