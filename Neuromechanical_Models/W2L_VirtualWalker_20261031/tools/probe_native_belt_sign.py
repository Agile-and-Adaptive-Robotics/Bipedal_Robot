from pathlib import Path
from lxml import etree as E
import copy,uuid,json,subprocess,numpy as np,hashlib
R=Path(__file__).resolve().parents[3];S=R/'Neuromechanical_Models/W2L_VirtualWalker_20261031/W2L_VirtualWalker_Standalone.asim';W=S.parent/'build/native_belt_sign';W.mkdir(exist_ok=True)
tree=E.parse(str(S));root=tree.getroot();root.find('ProjectPath').text=str(W)+'\\';root.find('SimEndTime').text='2.01';belt=next(b for b in root.iter('RigidBody') if b.findtext('Name')=='Treadmill Belt');joint=belt.find('Joint');charts=root.find('DataCharts');chart=copy.deepcopy(charts[0]);charts.clear();charts.append(chart)
def put(n,k,v):
 z=n.find(k)
 if z is None:z=E.SubElement(n,k)
 z.text=str(v)
for k,v in [('Name','Belt sign probe'),('OutputFilename','Belt sign probe.txt'),('EndTime',2),('CollectInterval',.001)]:put(chart,k,v)
cols=chart.find('DataColumns');cols.clear();labels=['TimeSlice','Time']
for label,target,typ in [('belt_world_X',belt.findtext('ID'),'WorldPositionX'),('belt_velocity_X',belt.findtext('ID'),'BodyLinearVelocityX'),('slider_position',joint.findtext('ID'),'JointPosition'),('slider_velocity',joint.findtext('ID'),'JointActualVelocity'),('slider_command',joint.findtext('ID'),'JointSetVelocity')]:
 c=E.SubElement(cols,'DataColumn')
 for k,v in [('ID',str(uuid.uuid5(uuid.NAMESPACE_URL,'belt_sign_native/'+label))),('Type','DataColumn'),('ColumnName',label),('TargetID',target),('DataType',typ)]:put(c,k,v)
 labels.append(label)
file=W/'probe.asim';tree.write(str(file),encoding='utf-8',xml_declaration=True);exe=Path(r'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe');q=subprocess.run([str(exe),str(file)],cwd=exe.parent,capture_output=True,text=True,timeout=20);(W/'run.log').write_text(q.stdout+'\n'+q.stderr);report={'source':str(S),'source_sha256':hashlib.sha256(S.read_bytes()).hexdigest(),'exit':q.returncode,'note':'Onlycopy outputpath/endtime/chartchanged; sourcephysics/neural/motorcommand untouched.'};print(q.returncode,q.stdout[-900:])
if q.returncode==0:
 a=np.loadtxt(W/'Belt sign probe.txt',skiprows=1);a=a[(a[:,1]<=2)&np.any(a[:,2:]!=0,axis=1)];np.savez_compressed(W/'probe.npz',data=a,labels=labels);sel=a[:,1]>=.6;report['samples']=len(a);report['belt_world_X_initial_final']=[float(a[0,2]),float(a[-1,2])];report['belt_world_displacement_m']=float(a[-1,2]-a[0,2]);report['post_0p6s_mean']={name:float(a[sel,j].mean()) for j,name in enumerate(labels) if j>=3};report['world_vs_joint_position_max_error_m']=float(np.max(np.abs((a[:,2]-a[0,2])+a[:,4])));report['world_velocity_plus_joint_velocity_max_error']=float(np.max(np.abs(a[:,3]+a[:,5])))
(W/'BELT_SIGN.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
