from pathlib import Path
import subprocess,time,json,numpy as np,argparse
parser=argparse.ArgumentParser();parser.add_argument('--variant-id',default='SIContact');args=parser.parse_args();assert args.variant_id.replace('_','').isalnum()
R=Path(__file__).resolve().parents[3];name='W2L_VirtualWalker_Newer_'+args.variant_id;W=R/('Neuromechanical_Models/'+name+'_20261031');exe=Path(r'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe')
t=time.perf_counter();p=subprocess.run([str(exe),str(W/(name+'_Standalone.asim'))],cwd=exe.parent,capture_output=True,text=True,timeout=80);(W/'native_smoke.log').write_text(p.stdout+'\n'+p.stderr);print(p.returncode,p.stdout[-1800:],p.stderr[-500:]);result={'exit_code':p.returncode,'wall_seconds':time.perf_counter()-t}
for filename in ['Virtual walker verification.txt','SI contact neural check.txt']:
 f=W/filename
 if f.exists() and f.stat().st_size:
  a=np.loadtxt(f,skiprows=1);a=a[(a[:,1]<=10)&np.any(a[:,2:]!=0,axis=1)];result[filename]={'rows':len(a),'finite':bool(np.isfinite(a).all()),'signal_min':float(a[:,2:].min()),'signal_max':float(a[:,2:].max())};np.savez_compressed(W/(f.stem+'.npz'),data=a)
  if filename.startswith('Virtual'):
   labels=json.loads((W/'diagnostic_columns.json').read_text());result['root_y_m']={k:float(v) for k,v in [('first',a[0,labels.index('root_y')]),('min',a[:,labels.index('root_y')].min()),('max',a[:,labels.index('root_y')].max()),('last',a[-1,labels.index('root_y')])]}
result['scope']='Passive wing/rail support; finite root height does not prove free walking. SI companion follows editable project; not identical to previous native208 graph.'
(W/'NATIVE_SMOKE.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
