"""Select and independently reconstruct an already evaluated feasible C fit."""
import argparse,json,time,sys,os
from pathlib import Path
for key in('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
HERE=Path(__file__).resolve().parent;sys.path[:0]=[str(HERE),str(HERE.parent)]
import numpy as np
from pose_replay import ContractReplay
from readout_keys import multiplier
def main():
 p=argparse.ArgumentParser();p.add_argument('--search',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);manifest=json.loads((a.search/'MANIFEST.json').read_text());rows=json.loads((a.search/'ITERATIONS.json').read_text());feasible=[r for r in rows if r['diagnostics']['accepted']];best=min(feasible,key=lambda r:r['training_RMSE']);r=ContractReplay(manifest['configuration'])
 for i,(pre,post)in enumerate(zip(r.c.pre,r.c.post)):
  source,target=r.names[pre],r.names[post]
  if target.startswith('MN_')and source.startswith('PF_')and not source.startswith('PF_IN'):r.base_g[i]*=multiplier(source,target,best)
 fixed=json.loads((Path(manifest['source_timing_block'])/'FINAL.json').read_text())['parameters'];repeat=r.evaluate(fixed,time.time()+180.);matched=abs(repeat['value']-best['value'])<1e-10;manifest.update(source_search=str(a.search),selection='Minimum training waveform error among previously evaluated candidates passing the original C peak/rhythm gates; no gate relaxed');(a.output/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));(a.output/'BEST.json').write_text(json.dumps(best,indent=2));(a.output/'FINAL.json').write_text(json.dumps(dict(accepted=bool(matched and repeat['diagnostics']['accepted']),fresh_reconstruction_matched=matched,best=best,repeat=repeat,scope=manifest['scope']),indent=2));np.savez_compressed(a.output/'BEST_TRACE.npz',t=r.t,rg=r.rg_trace,pf=r.pf_trace,mn=r.output,desired=r.desired,acts=np.array(r.acts),rg_names=np.array([r.names[i]for i in r.rg]),pf_names=np.array([r.names[i]for i in r.pf]));print(json.dumps(dict(iteration=best['iteration'],matched=matched,diagnostics=repeat['diagnostics'])))
if __name__=='__main__':main()
