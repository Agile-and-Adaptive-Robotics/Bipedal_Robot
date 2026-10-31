from pathlib import Path
from lxml import etree as E
import copy,uuid,json,subprocess,time,hashlib,numpy as np
BASE=Path('Neuromechanical_Models');EXE=Path(r'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe')
def put(n,k,v):
 q=n.find(k)
 if q is None:q=E.SubElement(n,k)
 q.text=str(v);return q
for name in ['W2L_Geometry_20261031','W2L_Geometry_Newer_20261031']:
 w=BASE/name;out=w/'verification';out.mkdir(exist_ok=True);t=E.parse(str(w/'W2L_Geometry_Standalone.asim'));r=t.getroot();windows=r.find('WindowMgr/Windows');windows.clear();put(r,'ProjectPath',str(out.resolve())+'\\');put(r,'SimEndTime','10.1');put(r,'PlaybackControlMode','0')
 for tx in r.iter('Texture'):
  if tx.text and tx.text.endswith('Treadmill_checker.bmp'):tx.text=str((w/'Treadmill_checker.bmp').resolve())
 bodies={b.findtext('Name'):b for b in r.findall('.//RigidBody') if b.findtext('Type')=='Box'};proto=next(b for b in r.findall('.//RigidBody') if b.findtext('Type')=='Attachment');probes=[]
 for leg in ['L','R']:
  for namebody,pos,label in [('femur_'+leg,(0,-.21,0),'femur_distal_'+leg),('tibia_'+leg,(0,0,-.21),'tibia_proximal_'+leg),('tibia_'+leg,(0,0,.21),'tibia_distal_'+leg),('foot_'+leg,(0,0,0),'foot_center_'+leg),('foot_'+leg,(.1,0,0),'foot_forward_'+leg)]:
   probe=copy.deepcopy(proto)
   for n in probe.iter('ID'):n.text=str(uuid.uuid4())
   put(probe,'Name',label);probe.find('Position').attrib.update(dict(zip('xyz',map(str,pos))))
   bodies[namebody].find('ChildBodies').append(probe);probes.append((probe.findtext('ID'),label))
 for bodyname in ['foot_L','foot_R','toe_L','toe_R']:
  body=bodies[bodyname];size=np.array([float(body.findtext('Length')),float(body.findtext('Height')),float(body.findtext('Width'))])
  for idx,signs in enumerate([(x,y,z) for x in [-1,1] for y in [-1,1] for z in [-1,1]]):
   probe=copy.deepcopy(proto)
   for n in probe.iter('ID'):n.text=str(uuid.uuid4())
   label='corner_'+bodyname+'_'+str(idx);put(probe,'Name',label);probe.find('Position').attrib.update(dict(zip('xyz',map(str,size*np.array(signs)/2))));body.find('ChildBodies').append(probe);probes.append((probe.findtext('ID'),label))
 charts=r.find('DataCharts');charts.clear();c=E.SubElement(charts,'DataChart')
 for k,v in dict(ID=str(uuid.uuid4()),Type='FileChart',Name='Geometry verification',OutputFilename='geometry.txt',SetStartEndTime='True',Enabled='True',CollectTimeWindow='-1',AlwaysActive='False',StartTime='0',EndTime='10',CollectInterval='.01').items():put(c,k,v)
 cs=E.SubElement(c,'DataColumns');signals=[]
 for target,label in probes+[(bodies['Root'].findtext('ID'),'root')]:
  for a in 'XYZ':signals.append((target,'WorldPosition'+a,label+'_'+a))
 for j in r.iter('Joint'):
  if j.findtext('Type')=='Hinge':signals.append((j.findtext('ID'),'JointRotation','angle_'+j.findtext('Name')))
 for target,typ,label in signals:
  col=E.SubElement(cs,'DataColumn')
  for k,v in dict(ID=str(uuid.uuid4()),ColumnName=label,Type='DataColumn',TargetID=target,DataType=typ).items():put(col,k,v)
 f=out/'geometry_test.asim';t.write(str(f),encoding='utf-8',xml_declaration=True);start=time.time()
 p=subprocess.run([str(EXE),str(f.resolve())],cwd=str(EXE.parent),capture_output=True,text=True,timeout=45)
 (out/'runner.log').write_text(p.stdout+'\n'+p.stderr,encoding='utf-8')
 if p.returncode:raise RuntimeError(p.stdout[-500:])
 a=np.loadtxt(out/'geometry.txt',skiprows=1);a=a[(a[:,1]<=10)&np.any(a[:,2:]!=0,axis=1)];labels=['TimeSlice','Time']+[s[2] for s in signals];report={'exit_code':p.returncode,'wall_seconds':time.time()-start,'end_time':float(a[-1,1]),'rows':len(a),'root_y_range_m':[float(a[:,labels.index('root_Y')].min()),float(a[:,labels.index('root_Y')].max())],'legs':{}}
 for leg in ['L','R']:
  f=a[:,[labels.index('femur_distal_'+leg+'_'+q) for q in 'XYZ']];t=a[:,[labels.index('tibia_proximal_'+leg+'_'+q) for q in 'XYZ']];gap=np.linalg.norm(f-t,axis=1)
  report['legs'][leg]={'native_initial_femur_distal_world_m':f[0].tolist(),'native_initial_tibia_proximal_world_m':t[0].tolist(),'initial_gap_m':float(gap[0]),'maximum_dynamic_gap_m':float(gap.max())}
 meta=json.loads((w/'GEOMETRY_REPORT.json').read_text())['models']['asim'];cornercols=[i for i,label in enumerate(labels) if label.startswith('corner_') and label.endswith('_Y')];minsole=float(a[0,cornercols].min());report['native_initial_minimum_solid_sole_y_m']=minsole;report['floor_top_y_m']=meta['floor_top_y_m'];report['native_initial_minimum_solid_sole_gap_m']=minsole-meta['floor_top_y_m'];report['joint_anatomical_ranges_deg']={j['joint']:[float(np.rad2deg(a[:,labels.index('angle_'+j['joint'])].min())+j['new_anatomical_zero_deg']),float(np.rad2deg(a[:,labels.index('angle_'+j['joint'])].max())+j['new_anatomical_zero_deg'])] for j in meta['hinge_limit_changes']};report['joint_angle_range_rad']={label:[float(a[:,i].min()),float(a[:,i].max())] for i,label in enumerate(labels) if label.startswith('angle_')}
 np.savez_compressed(out/'geometry.npz',data=a,labels=labels);(w/'NATIVE_TEST.json').write_text(json.dumps(report,indent=2));print(name,json.dumps(report),flush=True)
