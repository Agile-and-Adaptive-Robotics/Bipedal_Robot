from pathlib import Path
import subprocess,json,time,numpy as np,argparse
parser=argparse.ArgumentParser();parser.add_argument('--variant-id',default='Guided');args=parser.parse_args();assert args.variant_id.replace('_','').isalnum();name='W2L_VirtualWalker_'+args.variant_id
R=Path(__file__).resolve().parents[3];W=R/('Neuromechanical_Models/'+name+'_20261031');exe=Path(r'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe');t=time.perf_counter();p=subprocess.run([str(exe),str(W/(name+'_Standalone.asim'))],cwd=exe.parent,capture_output=True,text=True,timeout=60);(W/'run.log').write_text(p.stdout+'\n'+p.stderr);print(p.returncode,p.stdout[-1400:]);report={'exit':p.returncode,'wall_seconds':time.perf_counter()-t}
for filename,labelfile in [('Guide diagnostics.txt','guide_columns.json'),('Virtual walker verification.txt','diagnostic_columns.json')]:
 f=W/filename
 if f.exists() and f.stat().st_size:
  z=np.loadtxt(f,skiprows=1);z=z[(z[:,1]<=10)&np.any(z[:,2:]!=0,axis=1)];labels=json.loads((W/labelfile).read_text());np.savez_compressed(W/(f.stem+'.npz'),data=z,labels=labels);report[filename]={'samples':len(z),'finite':bool(np.isfinite(z).all()),'signals':{name:{'first':float(z[0,i]),'min':float(z[:,i].min()),'max':float(z[:,i].max()),'last':float(z[-1,i])} for i,name in enumerate(labels) if name.startswith(('root_','carriage_')) or 'ContactCount' in name}}
(W/'NATIVE_GUIDE_SMOKE.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2)[:6500])
