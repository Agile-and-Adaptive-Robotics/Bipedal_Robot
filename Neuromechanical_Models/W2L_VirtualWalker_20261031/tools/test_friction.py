from pathlib import Path
from lxml import etree as E
import subprocess,copy,uuid,json,numpy as np
w=Path('Neuromechanical_Models/W2L_VirtualWalker_20261031').resolve();src=w/'W2L_VirtualWalker_Standalone.asim';exe=Path(r'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe');labels=json.loads((w/'diagnostic_columns.json').read_text());summary=[]
for mu in [.05,.2,.5]:
 d=w/'build'/f'friction_{mu}';d.mkdir(parents=True,exist_ok=True);r=E.parse(str(src));mts=r.find('.//MaterialTypes');m=copy.deepcopy(mts[0]);mid=str(uuid.uuid4());m.find('ID').text=mid;m.find('Name').text='Rail trial material'
 for x in m:
  if x.tag.startswith('FrictionLinear'):x.text=str(mu)
 mts.append(m)
 for b in r.iter('RigidBody'):
  if 'hip wing' in b.findtext('Name','') or 'support rail' in b.findtext('Name',''):b.find('MaterialTypeID').text=mid
 r.find('ProjectPath').text=str(d)+'\\';fn=d/'trial.asim';r.write(str(fn),encoding='utf-8',xml_declaration=True)
 p=subprocess.run([str(exe),str(fn)],cwd=exe.parent,capture_output=True,text=True,timeout=80);(d/'runner.log').write_text(p.stdout+p.stderr)
 a=np.loadtxt(d/'Virtual walker verification.txt',skiprows=1);a=a[(a[:,1]<=10)&np.any(a[:,2:]!=0,axis=1)];np.savez_compressed(d/'signals.npz',data=a,labels=labels)
 y=a[:,labels.index('root_y')];fall=a[y<.75,1];s=dict(mu=mu,exit=p.returncode,root_y_min=float(y.min()),root_y_end=float(y[-1]),first_below_75cm=float(fall[0]) if len(fall) else None,x_end=float(a[-1,labels.index('root_x')]),z_end=float(a[-1,labels.index('root_z')]));summary.append(s);print(json.dumps(s),flush=True)
(w/'FRICTION_CHECK.json').write_text(json.dumps(summary,indent=2))
