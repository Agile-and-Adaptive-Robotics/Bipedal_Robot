"""Audit actual conductance-loaded Euler step ratios and raw neuron voltages.

For an ungated isolated compartment with frozen conductances, dt*(Gm+Gsyn)/C
above1 permits overshoot and above2 violates explicit-Euler linear stability.
This is not a complete recurrent-network stability proof. No circuit is changed.
"""
import argparse,json,hashlib,os,sys,time
from pathlib import Path
for key in('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));sys.path.insert(0,str(HERE.parent))
import numpy as np
import motor_adapter as M
def main():
 p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);final=json.loads((a.fit/'FINAL.json').read_text());manifest=json.loads((a.fit/'MANIFEST.json').read_text());peak=Path(manifest['fit'])
 ready=M.handle_request(dict(op='preflight',fit=str(peak),fit_sha256=hashlib.sha256((peak/'FINAL.json').read_bytes()).hexdigest()))
 if not ready['ready_for_optimization']:raise RuntimeError('Frozen peak fit does not reproduce')
 r=M.INSTANCE;c=r.c;core=r.core;max_ratio=np.zeros(len(core.V));over1=np.zeros(len(core.V),dtype=int);over2=over1.copy();lo=np.full(len(core.V),np.inf);hi=-lo.copy();steps=0;samples=0;forward=c.forward
 def observed(u):
  nonlocal steps,samples
  if steps*r.dt>=2*r.T:
   g=c.gmax*np.clip((core.V[c.pre]-c.elo)/c.erange,0.,1.);total=np.bincount(c.post,weights=g,minlength=len(core.V));ratio=core.time_factor_membrane*(core.g_m+total);np.maximum(max_ratio,ratio,out=max_ratio);over1[:]+=ratio>1.;over2[:]+=ratio>2.;samples+=1
  result=forward(u)
  if steps*r.dt>=2*r.T:np.minimum(lo,core.V,out=lo);np.maximum(hi,core.V,out=hi)
  steps+=1;return result
 c.forward=observed;result=M.handle_request(dict(op='evaluate',params=final['parameters']));rows=[]
 for i,name in enumerate(r.names):
  if name.startswith(('MN_','PF_','RG_','II','Ia','Ib','RC_')):rows.append(dict(neuron=name,maximum_Euler_load_ratio=float(max_ratio[i]),fraction_above1=float(over1[i]/samples),fraction_above2=float(over2[i]/samples),minimum_raw_voltage_mV=float(lo[i]),maximum_raw_voltage_mV=float(hi[i]),ionic_channels_present=bool(np.any(core.g_ion[:,i]))))
 rows.sort(key=lambda x:x['maximum_Euler_load_ratio'],reverse=True);out=dict(family=r.manifest['family'],source_D_fit=str(a.fit),model_dt_s=r.dt,samples=samples,replay=result,rows=rows,scope=__doc__);(a.output/'RESULT.json').write_text(json.dumps(out,indent=2));print(json.dumps(dict(family=out['family'],rmse=result['diagnostics']['mn_command_rmse'],top=rows[:8]),indent=2))
if __name__=='__main__':main()
