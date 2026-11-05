"""Diagnose a failed peak fit before expanding its parameter search."""
import argparse,json,os,sys,time
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
from replay_adapter import Replay
from muscle_map import classify

def main():
    p=argparse.ArgumentParser();p.add_argument('--block',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((a.block/'MANIFEST.json').read_text());r=Replay(cfg);best=json.loads((a.block/'BEST.json').read_text())['parameters'];rows=[]
    groups=('hip_ext','hip_flex','knee_ext','knee_flex','ankle_pf','ankle_df')
    scenarios=[('best',best)]+[('shared_'+str(g),{**best,**{'ff_'+k:g for k in groups}})for g in (1.,3.,10.,30.,100.)]
    for tag,params in scenarios:
        result=r.evaluate(params,time.time()+180.);keep=r.t>=2*r.T
        per_act=[]
        for i,act in enumerate(r.acts):
            if not r.selected[i]:continue
            ids=[j for j,(pre,post)in enumerate(zip(r.c.pre,r.c.post))if post==r.mn[i] and r.names[pre].startswith('PF_')and not r.names[pre].startswith('PF_IN')]
            per_act.append(dict(actuator=act,group=classify(act).groups[0],target_peak=float(r.desired[keep,i].max()),MN_peak=float(r.output[keep,i].max()),active_PF_MN_edges=int(sum(r.c.gmax[ids]>0)),maximum_PF_MN_gmax_uS=float(r.c.gmax[ids].max())if ids else 0.))
        rows.append(dict(case=tag,parameters=params,**result,per_actuator=per_act));(a.output/'RESULTS.json').write_text(json.dumps(rows,indent=2))
        print(json.dumps(dict(case=tag,**result)),flush=True)
    (a.output/'MANIFEST.json').write_text(json.dumps(cfg,indent=2))
if __name__=='__main__':main()
