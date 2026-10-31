"""Common PF stance/swing conductances as an explicit functional diagnostic.

Air training never established stance-load force. These additions are declared
engineering candidates, not fitted or literature-prescribed conductances.
"""
import sys,json,shutil,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from bilateral import assert_tied
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'stance_swing_pf_v1';out.mkdir()
base=json.loads((root/'shared_force_coupling_v1/shared_c1p2.json').read_text());source=Path(base['weights_directory']);plans={}
for strength in (.15,.3,.6):
 name='stance_'+str(strength).replace('.','p');W=out/name/'weights';W.mkdir(parents=True)
 for p in source.glob('*.json'):shutil.copyfile(p,W/p.name)
 p=W/'mn_weights.json';z=json.loads(p.read_text());changes=[]
 for side in ('r','l'):
  for pre,act,amount in [('PF_KNEE-E_',a,strength) for a in ('vas_med','vas_int','vas_lat')]+[('PF_HIP-E_',a,strength*.5) for a in ('glut_max1','glut_max2','glut_max3')]+[('PF_FUNC_ANKLE_Z_P_','soleus',.15)]:
   edge=pre+side+' -> MN_'+act+'_'+side;assert edge in z['edges'],edge;z['edges'][edge]+=amount;changes.append(dict(edge=edge,increment=amount))
 assert_tied(z['edges']);z['functional_stance_recruitment']=dict(changes=changes,scope='Engineering PF stance-load recruitment; not learned from ground targets.')
 p.write_text(json.dumps(z,indent=2));plans[name]={**base,'weights_directory':str(W)}
def run(item):
 name,cfg=item;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,plans.items()):
  rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment',{}).get('accepted'),flush=True)
