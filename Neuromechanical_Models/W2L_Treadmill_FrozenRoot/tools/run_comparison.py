"""Matched W2L mechanics; saved legacy and bilateral controllers, without retuning."""
from pathlib import Path
import xml.etree.ElementTree as E,copy,uuid,subprocess,json,concurrent.futures,hashlib
import numpy as np

REPO=Path(__file__).resolve().parents[3]
OUT=REPO/'Neuromechanical_Models/W2L_Comparison_20261030';OUT.mkdir(exist_ok=True)
EXE=Path(r'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe')
LEG=REPO/'Neuromechanical_Models/W2L_Treadmill_FrozenRoot/W2L_Treadmill_FrozenRoot_Standalone.asim'
NEW=REPO/'Neuromechanical_Models/Walker_2_Layer_CPG_BilateralRG/Walker_2_Layer_CPG_BilateralRG_Standalone_modern.asim'
def put(n,k,v):
 c=n.find(k)
 if c is None:c=E.SubElement(n,k)
 c.text=str(v);return c
def uid():return str(uuid.uuid4())
def add_stim(r,fields):
 s=E.SubElement(r.find('ExternalStimuli'),'Stimulus')
 for k,v in dict(ID=uid(),Name='trial stimulus',AlwaysActive='False',Enabled='True',**fields).items():put(s,k,v)
 return s
def make(model,condition,name='baseline',pulse=None):
 folder=OUT/model/condition/name;folder.mkdir(parents=True,exist_ok=True)
 r=E.parse(LEG if model=='legacy' else NEW).getroot()
 legacy=E.parse(LEG).getroot()
 env=r.find('Environment');env.remove(env.find('Structures'));env.append(copy.deepcopy(legacy.find('Environment/Structures')))
 ground=next(s for s in env.find('Structures') if s.findtext('Name')=='Ground')
 ground.find('Position').set('y','-2')
 ground.find('RigidBody/Position').set('y','-2')
 windows=r.find('WindowMgr/Windows')
 if windows is not None:windows.clear()
 org=env.find('.//Organism');original_root=org.find('RigidBody');org.remove(original_root);org.append(copy.deepcopy(legacy.find('.//Organism/RigidBody')))
 root=org.find('RigidBody');put(root,'Freeze','True' if condition in ['air','treadmill'] else 'False')
 if model=='bilateral':
  # The new controller retains all its saved II/contact adapters and RG connections.
  for s in list(r.find('ExternalStimuli')):
   if s.findtext('Type')=='MotorVelocity':r.find('ExternalStimuli').remove(s)
  motor=next(s for s in legacy.find('ExternalStimuli') if s.findtext('Type')=='MotorVelocity')
  r.find('ExternalStimuli').append(copy.deepcopy(motor))
 for s in r.find('ExternalStimuli'):
  if s.findtext('Type')=='MotorVelocity':put(s,'Equation','0.1' if condition=='treadmill' else '0')
 if condition=='air':
  path=next(s for s in env.find('Structures') if s.findtext('Name')=='Path')
  path.find('Position').set('y','-2');path.find('RigidBody/Position').set('y','-2')
 mass_g=sum(float(b.findtext('Mass','0')) for b in root.iter('RigidBody') if b.findtext('IsCollisionObject')=='True' and b.findtext('IsContactSensor')!='True')
 force=0.5*mass_g/1000*abs(float(env.findtext('Gravity','-9.8')))
 if condition=='partial50':
  s=add_stim(r,dict(Type='ForceInput',StructureID=org.findtext('ID'),BodyID=root.findtext('ID'),StartTime='0',EndTime='10',ForceX='0',ForceY=str(force),ForceZ='0',TorqueX='0',TorqueY='0',TorqueZ='0'))
  E.SubElement(s,'RelativePosition',x='0',y='0',z='0')
 if pulse:
  target,start,duration,amps=pulse
  nid=next(n.findtext('ID') for n in r.iter('Neuron') if n.findtext('Name')==target)
  add_stim(r,dict(Type='Current',CurrentType='Tonic',OrganismID=org.findtext('ID'),TargetNodeID=nid,StartTime=start,EndTime=start+duration,CycleOnDuration=duration,CycleOffDuration='0',BurstOnDuration=duration,BurstOffDuration='0',CurrentOn=amps,CurrentOff='0',CurrentBurstOff='0'))
 put(r,'ProjectPath',str(folder)+'\\');put(r,'SimEndTime','10.1');put(r,'PlaybackControlMode','0')
 charts=r.find('DataCharts');charts.clear();chart=E.SubElement(charts,'DataChart')
 for k,v in dict(ID=uid(),Type='FileChart',Name='Signals',OutputFilename='signals.txt',SetStartEndTime='True',CollectTimeWindow='-1',AlwaysActive='False',StartTime='0',EndTime='10',CollectInterval='.001').items():put(chart,k,v)
 cols=E.SubElement(chart,'DataColumns');signals=[]
 for n in r.iter('Neuron'):signals.append((n.findtext('ID'),'MembraneVoltage',n.findtext('Name')))
 for n in root.iter('Joint'):
  if n.findtext('Type')=='Hinge':signals.append((n.findtext('ID'),'JointRotation','angle_'+n.findtext('Name')))
 for axis in 'XYZ':signals.append((root.findtext('ID'),'WorldPosition'+axis,'root_'+axis.lower()))
 for n in root.iter('RigidBody'):
  if n.findtext('IsContactSensor')=='True':signals.append((n.findtext('ID'),'ContactCount','contact_'+n.findtext('Name')))
 belt=next(b for b in env.iter('RigidBody') if b.findtext('Name')=='Treadmill Belt')
 put(belt,'Texture',str(LEG.parent/'Treadmill_checker.bmp'))
 signals.append((belt.findtext('ID'),'WorldPositionX','belt_x'))
 for target,typ,label in signals:
  c=E.SubElement(cols,'DataColumn')
  for k,v in dict(ID=uid(),ModuleName='',Type='DataColumn',ColumnName=label,TargetID=target,DataType=typ).items():put(c,k,v)
 sim=folder/'trial.asim';E.ElementTree(r).write(sim,encoding='utf-8',xml_declaration=True)
 meta=dict(model=model,condition=condition,name=name,pulse=pulse,declared_dynamic_mass_kg=mass_g/1000,support_force_N=force if condition=='partial50' else 0,neurons=sum(1 for n in r.iter('Neuron')),source=str(LEG if model=='legacy' else NEW),labels=['TimeSlice','Time']+[s[2] for s in signals])
 (folder/'protocol.json').write_text(json.dumps(meta,indent=2));return sim
def run(sim):
 folder=sim.parent;p=subprocess.run([str(EXE),str(sim)],cwd=EXE.parent,capture_output=True,text=True,timeout=80)
 (folder/'runner.log').write_text(p.stdout+'\n'+p.stderr)
 if p.returncode or 'Error occurred' in p.stdout:raise RuntimeError(str(sim)+' '+p.stdout[-600:])
 a=np.loadtxt(folder/'signals.txt',skiprows=1);a=a[(a[:,1]<=10)&np.any(a[:,2:]!=0,axis=1)]
 meta=json.loads((folder/'protocol.json').read_text());np.savez_compressed(folder/'signals.npz',data=a,labels=meta['labels'])
 (folder/'signals.txt').unlink();print('completed',str(folder.relative_to(OUT)),flush=True)
 return sim
def onsets(a,col):
 v=a[:,col];t=a[:,1];cut=(np.quantile(v[t>1],.1)+np.quantile(v[t>1],.9))/2
 ids=np.where((v[1:]>=cut)&(v[:-1]<cut))[0]+1
 ids=ids[t[ids]>1];return t[ids]
if __name__=='__main__':
 jobs=[make(m,c) for m in ['legacy','bilateral'] for c in ['air','treadmill','drop','partial50']]
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:list(ex.map(run,jobs))
 jobs=[];plan=[]
 for m in ['legacy','bilateral']:
  for c in ['air','treadmill','drop','partial50']:
   for label,target in [('deng2019_PF','L Hip PF ext'),('deng2019_RG','L RG ext')]:jobs.append(make(m,c,label,(target,2,.1,-10e-9)))
   if c not in ['air','treadmill']:continue
   z=np.load(OUT/m/c/'baseline/signals.npz');labels=list(z['labels']);a=z['data'];ts=onsets(a,labels.index('L RG ext'));ts=ts[(ts>3)&(ts<7)]
   if len(ts)<4:plan.append(dict(model=m,condition=c,prc='undefined: insufficient baseline cycles'));continue
   period=float(np.median(np.diff(ts)));cv=float(np.std(np.diff(ts))/np.mean(np.diff(ts)));anchor=float(ts[-1]);plan.append(dict(model=m,condition=c,period=period,period_cv=cv,phase_reference='L RG ext adaptive-midrange upward crossing',stationary=cv<.1))
   for sign in [-1,1]:
    for phi in np.arange(.05,1,.1):
     name=f'deng2022_{sign:+d}_phase{phi:.2f}';jobs.append(make(m,c,name,('L RG flx IN',anchor+float(phi)*period,.05*period,sign*.5e-9)))
 (OUT/'phase_plan.json').write_text(json.dumps(plan,indent=2))
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:list(ex.map(run,jobs))
 (OUT/'COMPLETE.json').write_text(json.dumps(dict(baselines=8,perturbed=len(jobs),sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [LEG,NEW]}),indent=2))
