from pathlib import Path
import subprocess, numpy as np,json,time,sys
folder='W2L_VirtualWalker_Newer_20261031' if '--newer' in sys.argv else 'W2L_VirtualWalker_20261031'
name='W2L_VirtualWalker_Newer' if '--newer' in sys.argv else 'W2L_VirtualWalker'
w=Path('Neuromechanical_Models')/folder;w=w.resolve();exe=Path(r'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe')
t=time.time();p=subprocess.run([str(exe),str(w/(name+'_Standalone.asim'))],cwd=exe.parent,capture_output=True,text=True,timeout=80)
(w/'runner.log').write_text(p.stdout+'\n'+p.stderr,encoding='utf-8')
print(p.returncode,time.time()-t,p.stdout[-1500:],p.stderr[-300:])
f=w/'Virtual walker verification.txt'
if f.exists():
 a=np.loadtxt(f,skiprows=1);labels=json.loads((w/'diagnostic_columns.json').read_text());a=a[(a[:,1]<=10)&np.any(a[:,2:]!=0,axis=1)]
 np.savez_compressed(w/'verification.npz',data=a,labels=labels)
 summary={x:{'first':float(a[0,i]),'min':float(a[:,i].min()),'max':float(a[:,i].max()),'last':float(a[-1,i])} for i,x in enumerate(labels) if x in ['root_x','root_y','root_z','Left hip wing_y','Right hip wing_y'] or 'contacts' in x or 'angle' in x or 'MN' in x}
 (w/'TEST_RESULTS.json').write_text(json.dumps(dict(exit_code=p.returncode,wall_seconds=time.time()-t,samples=len(a),signals=summary),indent=2),encoding='utf-8');print(json.dumps(summary,indent=2)[:6000])
