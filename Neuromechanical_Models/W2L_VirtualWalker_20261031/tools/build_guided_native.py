"""Separate passive-guide AnimatLab companion; source robot and neural objects preserved."""
from pathlib import Path
from historical_build_guard import require_historical_regeneration
require_historical_regeneration()
from lxml import etree as E
import copy,uuid,json,hashlib,shutil,argparse
R=Path(__file__).resolve().parents[3];S=R/'Neuromechanical_Models/W2L_VirtualWalker_20261031'
parser=argparse.ArgumentParser();parser.add_argument('--clearance',type=float,default=.01);parser.add_argument('--rail-drop',type=float,default=.02);parser.add_argument('--variant-id',default='Guided');parser.add_argument('--neural-clock',type=float,default=1.);parser.add_argument('--source-folder',type=Path,default=S);parser.add_argument('--source-name',default='W2L_VirtualWalker');args=parser.parse_args();S=args.source_folder;c=args.clearance;assert .001<=c<=.1;assert args.neural_clock>0;assert args.variant_id.replace('_','').isalnum();NAME='W2L_VirtualWalker_'+args.variant_id;W=R/('Neuromechanical_Models/'+NAME+'_20261031');W.mkdir(exist_ok=True)
P=E.XMLParser(strip_cdata=False)
def uid(s):return str(uuid.uuid5(uuid.NAMESPACE_URL,'ben-w2l-native-guided-20261031/'+s))
def put(n,k,v):
 z=n.find(k)
 if z is None:z=E.SubElement(n,k)
 z.text=str(v);return z
def num(n,k,v,g):
 z=put(n,k,'' if g else v)
 if g:z.attrib.update({'Value':str(v),'Scale':'None','Actual':str(v)})
def vec(n,k,v,g):
 z=n.find(k)
 if z is None:z=E.SubElement(n,k)
 if g:
  for a,x in zip('XYZ',v):num(z,a,x,True)
 else:z.attrib.update(dict(zip('xyz',map(str,v))))
def xyz(n,k,g):
 z=n.find(k);return [float(z.find(a).get('Actual')) for a in 'XYZ'] if g else [float(z.get(a,'0')) for a in 'xyz']
def position(n,v,g):
 vec(n,'LocalPosition' if g else 'Position',v,g)
 if g:put(n,'LocalMatrix',','.join(map(str,[1,0,0,0,0,1,0,0,0,0,1,0,*[x*10 for x in v],1])))
def reroll(n,key):
 for j,q in enumerate(n.iter('ID')):q.text=uid(key+'/'+str(j))
def box(template,key,name,pos,dims,mass,gui,collide=True):
 n=copy.deepcopy(template)
 for q in list(n):
  if q.tag in ('ChildBodies','Joint'):n.remove(q)
 reroll(n,key);put(n,'ID',uid(key));put(n,'Name',name);position(n,pos,gui);vec(n,'Rotation',[0,0,0],gui);vec(n,'COM',[0,0,0],gui)
 for a,x in zip(('Length','Height','Width'),dims):num(n,a,x,gui)
 for a,x in [('Mass',mass),('Density',mass/(dims[0]*dims[1]*dims[2]*1e6))]:num(n,a,x,gui)
 for a,x in [('Freeze','False'),('IsVisible','True'),('IsContactSensor','False'),('IsCollisionObject',str(collide)),('Texture','')]:put(n,a,x)
 return n
source_hashes={};reports=[]
for gui,ext in [(True,'.aproj'),(False,'_Standalone.asim')]:
 src=S/(args.source_name+ext);source_hashes[str(src)]=hashlib.sha256(src.read_bytes()).hexdigest();tree=E.parse(str(src),P);root=tree.getroot();env=root.find('.//Environment');physical=next(b for b in root.iter('RigidBody') if b.findtext('Name')=='Root');px,py,pz=xyz(physical,'LocalPosition' if gui else 'Position',gui)
 neural=[E.tostring(n,method='c14n') for n in root.findall('.//NervousSystem')]
 belt=next(b for b in root.iter('RigidBody') if b.findtext('Name')=='Treadmill Belt');joint_template=belt.find('Joint');assert joint_template.findtext('Type')=='Prismatic'
 rail=next(b for b in root.iter('RigidBody') if b.findtext('Name')=='Left support rail');rp=xyz(rail,'LocalPosition' if gui else 'Position',gui);rail_top=rp[1]+.02-args.rail_drop
 # Shift both rails and their enclosing structure positions, without touching robot placement.
 for b in [b for b in root.iter('RigidBody') if b.findtext('Name') in ('Left support rail','Right support rail')]:
  bp=xyz(b,'LocalPosition' if gui else 'Position',gui);bp[1]-=args.rail_drop;position(b,bp,gui);st=b.getparent();sp=xyz(st,'LocalPosition' if gui else 'Position',gui);sp[1]-=args.rail_drop;position(st,sp,gui)
 height=.3;bottom=rail_top+.005;center=bottom+height/2;cap_inner=py-.04+.0125+c;half_x=.06+c+.02;origin=[px,cap_inner+.01,.3];com=[0,rail_top+.15-origin[1],-.3]
 geom=[]
 for side in (-1,1):
  geom.append(('side_'+str(side),[0,center,side*(.4+c+.01)],[2*half_x,height,.02]))
  geom.append(('cap_'+str(side),[0,cap_inner+.01,side*.3],[2*half_x,.02,.2]))
  for front in (-1,1):geom.append(('post_'+str(side)+'_'+str(front),[front*(.06+c+.01),center,side*.3],[.02,height,.2]))
 # The left cap itself is the dynamic compound body's root shape. This avoids
 # introducing a hidden central collision box inside the robot pelvis.
 root_geom=next(q for q in geom if q[0]=='cap_1')
 # The existing frozen left horizontal rail is a valid Vortex joint anchor.
 # Reusing it avoids adding an extra fixed collision box to the apparatus.
 anchor_pos=xyz(rail,'LocalPosition' if gui else 'Position',gui);carrier_pos=[a-b for a,b in zip(origin,anchor_pos)]
 carrier=box(physical,'carrier','Guide carriage cap',carrier_pos,root_geom[2],2000,gui);vec(carrier,'COM',com,gui)
 mat=rail.findtext('MaterialTypeID');put(carrier,'MaterialTypeID',mat)
 joint=copy.deepcopy(joint_template);reroll(joint,'joint');put(joint,'ID',uid('joint'));put(joint,'Name','Passive guide X slider');put(joint,'EnableMotor','False');position(joint,[0,0,0],gui);vec(joint,'Rotation',[0,0,0],gui)
 # Vortex requires strictly positive capacity parameters even with motor off.
 for k,v in [('MaxForce',10000),('MaxVelocity',.2)]:num(joint,k,v,gui)
 for item,value in [('LowerLimit',-1.9),('UpperLimit',1.9)]:num(joint.find(item),'LimitPos',value,gui)
 if joint.find('Friction') is not None:put(joint.find('Friction'),'Enabled','False')
 for q in list(joint):
  if q.tag.startswith('Relaxation'):put(q,'Enabled','False')
 carrier.insert(0,joint);parts=E.SubElement(carrier,'ChildBodies')
 for key,world,dims in geom:
  if key==root_geom[0]:continue
  local=[world[0],world[1]-origin[1],world[2]-origin[2]];b=box(physical,key,'Guide '+key,local,dims,.001,gui);put(b,'MaterialTypeID',mat);parts.append(b)
 anchor_children=rail.find('ChildBodies')
 if anchor_children is None:anchor_children=E.SubElement(rail,'ChildBodies')
 anchor_children.append(carrier)
 assert physical.find('Joint') is None and physical.findtext('Freeze')=='False'
 assert neural==[E.tostring(n,method='c14n') for n in root.findall('.//NervousSystem')]
 if args.neural_clock!=1:
  neuronal=[n for n in root.findall('.//NervousSystem//Node') if '.Neurons.' in (n.findtext('ClassName') or '')] if gui else list(root.iter('Neuron'))
  for n in neuronal:
   for path in (['TimeConstant','CaActivation/TimeConstant','CaDeactivation/TimeConstant'] if gui else ['TimeConst','CaActivation/TimeConstant','CaDeactivation/TimeConstant']):
    z=n.find(path)
    if z is not None:
     if gui:z.set('Value',str(float(z.get('Value'))*args.neural_clock));z.set('Actual',str(float(z.get('Actual'))*args.neural_clock))
     else:z.text=str(float(z.text)*args.neural_clock)
 if gui:
  put(root,'ProjectName',NAME);put(root,'SimulationFile',NAME+'.asim');put(root.find('Simulation'),'ProjectPath',str(W)+'\\')
  for n in root.findall('.//NervousSystem//Node'):
   if n.find('DiagramXml') is not None:E.fromstring(n.findtext('DiagramXml').encode())
 else:
  put(root,'ProjectPath',str(W)+'\\');charts=root.find('DataCharts');chart=copy.deepcopy(charts[0]);put(chart,'ID',uid('chart'));put(chart,'Name','Guide diagnostics');put(chart,'OutputFilename','Guide diagnostics.txt');cols=chart.find('DataColumns');cols.clear();labels=['TimeSlice','Time']
  targets=[('root',physical.findtext('ID')),('carriage',uid('carrier')),('belt',belt.findtext('ID'))]+[(key,uid('carrier') if key=='cap_1' else uid(key)) for key,_,_ in geom]+[('left_wing',b.findtext('ID')) for b in root.iter('RigidBody') if b.findtext('Name')=='Left hip wing']
  for label,target in targets:
   for axis in 'XYZ':
    name=label+'_pos_'+axis;cq=E.SubElement(cols,'DataColumn')
    for k,v in [('ID',uid(name)),('Type','DataColumn'),('ColumnName',name),('TargetID',target),('DataType','WorldPosition'+axis)]:put(cq,k,v)
    labels.append(name)
   for typ in (['RotationX','RotationY','RotationZ','Mass']+[base+axis for base in ('BodyAngularVelocity','BodyLinearVelocity','BodyForce','BodyTorque') for axis in 'XYZ'] if label in ('root','carriage','belt') else ['ContactCount']):
    name=label+'_'+typ;cq=E.SubElement(cols,'DataColumn')
    for k,v in [('ID',uid(name)),('Type','DataColumn'),('ColumnName',name),('TargetID',target),('DataType',typ)]:put(cq,k,v)
    labels.append(name)
  charts.append(chart);(W/'guide_columns.json').write_text(json.dumps(labels,indent=2))
 tree.write(str(W/(NAME+ext)),encoding='utf-8',xml_declaration=True)
 reports.append({'gui':gui,'carriage_frame_world_origin':origin,'carriage_COM_local_m':com,'source_prismatic_rotation':[0,0,0],'native_axis_note':'copied unpowered original floor slider, source signed axis -X; symmetric travel limits','geometry':[{'name':k,'center_relative_initial_root_x_world_yz':v,'dimensions_xyz':d} for k,v,d in geom]})
shutil.copy2(W/(NAME+'_Standalone.asim'),W/(NAME+'.asim'))
for f in S.glob('*.aform'):shutil.copy2(f,W/f.name)
shutil.copy2(S/'diagnostic_columns.json',W/'diagnostic_columns.json')
report={'source_hashes':source_hashes,'clearance_m':c,'horizontal_rail_lowering_m':args.rail_drop,'slider_travel_m':[-1.9,1.9],'motor_enabled':False,'root_frozen':False,'root_joint_added':False,'neural_clock':args.neural_clock,'neural_clock_scope':'Scales neuronal membrane TimeConst and Ca activation/deactivation TimeConstant only, matching restricted graded MuJoCo kernel. No physics/muscle timestep or timeconstants changed. Accommodation/AHP spike-related times remain source values because the graded port has no spikes. No gain fitting.','nominal_parent_mass_kg':2,'seven_compound_children_declared_mass_kg_each':1e-6,'inertia_note':'Native inferred compound inertia differs from MuJoCo explicit principal inertia; carriage rotation is constrained by apparatus slider, robot remains free. COM frame adapter is explicit.','contact_note':'Reuses native rail material linearfriction.05; MuJoCo solref/solimp and friction torsion/rolling are not claimed equivalent. ContactCount on ordinary colliders is logged but may not report true reaction contacts; no aggregate guide force is fabricated.','geometry':reports}
(W/'BUILD_GUIDE.json').write_text(json.dumps(report,indent=2));(W/'.gitignore').write_text('*.aproj\n*.asim\n*.aform\n*.npz\n*.txt\n*.log\n!README.txt\n');print(json.dumps({k:v for k,v in report.items() if k!='geometry'},indent=2))
