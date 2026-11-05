"""Factorial intrinsic/commissural check; do not hide a structural failure in TPE."""
import argparse,copy,json,os,sys,time
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
from replay_adapter import Replay
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--seed',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    cfg=json.loads(a.manifest.read_text());seed=json.loads(a.seed.read_text());rows=[]
    import params
    old=dict(nap=dict(params.NAP),tau_h=params.TAU['rg_nap_h'],mutual_inhibition=params.G['rg_mutual_inh'])
    for core,spec in (('original',old),('accepted_Hunt_RG',seed)):
        for c1 in (1.,.03):
            config=copy.deepcopy(cfg);config['controller_profile']['rg_design']={k:copy.deepcopy(spec[k])for k in ('nap','tau_h','mutual_inhibition')};config['controller_profile']['rg_design']['c1_gain']=c1
            r=Replay(config);result=r.evaluate({},time.time()+180.);tag=core+'_c1_'+str(c1)
            rows.append(dict(case=tag,**result));(a.output/'RESULTS.json').write_text(json.dumps(rows,indent=2))
            np.savez_compressed(a.output/(tag+'.npz'),t=r.t,rg=r.rg_trace)
            print(json.dumps(dict(case=tag,**result)),flush=True)
    (a.output/'MANIFEST.json').write_text(json.dumps(dict(source_manifest=str(a.manifest),seed=str(a.seed),old=old),indent=2))
if __name__=='__main__':main()
