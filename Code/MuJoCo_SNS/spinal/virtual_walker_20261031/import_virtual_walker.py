"""Direct, y-up mechanical import of the specified AnimatLab virtual walkers.

Keeps source projects untouched. Import validation precedes controller fitting.
The old W2L converter deliberately is not reused: it drops support structures.
"""
from __future__ import annotations
import argparse, hashlib, json, math, os, re, shutil
from pathlib import Path
import xml.etree.ElementTree as E
os.environ.setdefault('CONDA_PREFIX',r'D:\Anaconda\envs\myo')
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
import mujoco
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3]
FRAME=np.diag([1.,1.,-1.])  # Source L is +Z; anatomical target L is -Z.
def number(e,k,d=0.):return float(e.findtext(k,str(d)))
def vec(e,k):
 x=e.find(k);return np.zeros(3) if x is None else np.array([float(x.get(c,'0')) for c in 'xyz'])
def rot(v):
 x,y,z=v;cx,sx=math.cos(x),math.sin(x);cy,sy=math.cos(y),math.sin(y);cz,sz=math.cos(z),math.sin(z)
 return np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]])@np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]])@np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]])
def quat(R):
 q=np.empty(4);mujoco.mju_mat2Quat(q,np.array(R,order='C').ravel());return q
def fmt(x):return ' '.join(f'{v:.14g}' for v in np.asarray(x).ravel())
def safe(s):return re.sub(r'[^a-zA-Z0-9_]+','_',s)
def build(newer=False,source_override=None):
 name='W2L_VirtualWalker'+('_Newer' if newer else '');source=REPO/'Neuromechanical_Models'/(name+'_20261031')/(name+'_Standalone.asim')
 if source_override:source=Path(source_override).resolve();name=source.stem.replace('_Standalone','')
 dest=HERE/'runs'/name;dest.mkdir(parents=True,exist_ok=True)
 raw=source.read_bytes();sha=hashlib.sha256(raw).hexdigest();snap=dest/f'source_{sha[:12]}.asim'
 if not snap.exists():snap.write_bytes(raw)
 root=E.fromstring(raw);org=root.find('.//Organism');rb=org.find('RigidBody')
 mj=E.Element('mujoco',model=name);E.SubElement(mj,'compiler',angle='radian',autolimits='true',inertiafromgeom='false')
 E.SubElement(mj,'statistic',extent='3',center=fmt(FRAME@vec(rb,'Position')+np.array([0.,-.45,0.])))
 E.SubElement(mj,'option',timestep='.002',gravity='0 -9.81 0',integrator='implicitfast',iterations='100',tolerance='1e-10')
 visual=E.SubElement(mj,'visual');E.SubElement(visual,'global',azimuth='90',elevation='-15')
 default=E.SubElement(mj,'default');E.SubElement(default,'geom',friction='1 .005 .0001',solref='.006 1',solimp='.95 .999 .001',condim='3')
 world=E.SubElement(mj,'worldbody');E.SubElement(world,'light',pos='0 4 2',dir='0 -1 -.2',directional='true')
 allb={};joints=[];muscles=[];receptors=[];springs=[];mass=0.;geometry=[];sensors=[];sites={};warnings=[]
 def visit(b,container,worldp,worldR,isroot=False,static=False):
  nonlocal mass
  ident=b.findtext('ID');typ=b.findtext('Type');nameb=b.findtext('Name');p=FRAME@vec(b,'Position');R=FRAME@rot(vec(b,'Rotation'))@FRAME;wp=worldp+worldR@p;wR=worldR@R
  allb[ident]=dict(name=nameb,type=typ,wp=wp,R=wR,element=b)
  if typ=='Attachment':
   sn='site_'+safe(nameb)+'_'+ident[:8];E.SubElement(container,'site',name=sn,pos=fmt(p),size='.002',rgba='.1 .7 .2 1');sites[ident]=sn;allb[ident]['site']=sn;return
  if typ in ('LinearHillMuscle','LinearHillStretchReceptor','Spring'):
   if typ=='LinearHillMuscle':muscles.append(b)
   elif typ=='LinearHillStretchReceptor':receptors.append(b)
   elif typ=='Spring':springs.append(b)
   return
  if typ not in ('Box','Plane'):warnings.append(f'Unsupported body type: {nameb}: {typ}');return
  sensor=b.findtext('IsContactSensor')=='True'
  node=E.SubElement(container,'body',name=safe(nameb),pos=fmt(p),quat=fmt(quat(R)))
  allb[ident]['mjbody']=safe(nameb)
  if isroot and not static:E.SubElement(node,'freejoint',name='root')
  if isroot and not static:E.SubElement(node,'camera',name='anatomical_oblique',pos='2 .1 2',xyaxes='.70710678 0 -.70710678 -.134974 .981606 -.134974',mode='trackcom',fovy='35')
  jt=b.find('Joint')
  if jt is not None:
   jtype=jt.findtext('Type');jn=jt.findtext('Name');jr=FRAME@rot(vec(jt,'Rotation'))@FRAME;axis=jr@np.array([-1. if jtype=='Prismatic' else 1.,0,0]);lo=number(jt.find('LowerLimit'),'LimitPos') if jt.find('LowerLimit') is not None else 0.;hi=number(jt.find('UpperLimit'),'LimitPos') if jt.find('UpperLimit') is not None else 0.
   attrs=dict(name=safe(jn),type='slide' if jtype=='Prismatic' else 'hinge',pos=fmt(FRAME@vec(jt,'Position')),axis=fmt(axis),damping='0',armature='.00001',limited='false')
   if jt.findtext('EnableLimits')=='True' and hi>lo:attrs.update(limited='true',range=fmt([lo,hi]),solreflimit='.004 1',solimplimit='.999 .9999 .001')
   E.SubElement(node,'joint',**attrs);joints.append(dict(name=jn,type=jtype,axis_local=axis.tolist(),axis_world=(wR@axis).tolist(),range=[lo,hi],enable_limits=attrs['limited'],source_id=jt.findtext('ID')))
  if typ=='Box':
   # AnimatLab dimensions are x=Length, y=Height, z=Width.
   dims=np.array([number(b,'Length'),number(b,'Height'),number(b,'Width')]);kg=number(b,'Mass')/1000
   if (not static or jt is not None) and kg>0:
    inertia=kg/12*np.array([dims[1]**2+dims[2]**2,dims[0]**2+dims[2]**2,dims[0]**2+dims[1]**2]);E.SubElement(node,'inertial',pos=fmt(FRAME@vec(b,'COM')),mass=f'{kg:.14g}',diaginertia=fmt(inertia));mass+=kg
   rail='support rail' in nameb;wing='hip wing' in nameb
   guide=nameb.startswith('Guide ');color='.2 .8 .85 1' if wing else ('.6 .6 .6 .3' if rail or guide else ('.1 .5 .8 1' if nameb.endswith('_L') else '.9 .65 .15 1'))
   coll=b.findtext('IsCollisionObject')!='False' and not sensor
   E.SubElement(node,'geom',name='geom_'+safe(nameb),type='box',size=fmt(dims/2),rgba=color,contype='1' if coll else '0',conaffinity='1' if coll else '0',friction='.05 .005 .0001' if wing or rail or guide else '1 .005 .0001')
   geometry.append(dict(name=nameb,dims_xyz_m=dims.tolist(),world_center=wp.tolist(),world_R=wR.tolist(),collision=coll,mass_kg=0 if static else kg,source_sensor=sensor))
   if sensor:sensors.append(nameb)
  elif typ=='Plane':E.SubElement(node,'geom',name='geom_'+safe(nameb),type='plane',size='10 10 .1',quat='0 1 0 0',rgba='.8 .8 .8 1')
  for child in b.findall('ChildBodies/RigidBody'):visit(child,node,wp,wR,static=static)
 visit(rb,world,np.zeros(3),np.eye(3),isroot=True)
 # Source structure Position is duplicated in its root body's Position; use once.
 for structure in root.findall('.//Environment/Structures/Structure'):
  b=structure.find('RigidBody')
  if b is not None:visit(b,world,np.zeros(3),np.eye(3),isroot=True,static=True)
 tendon=E.SubElement(mj,'tendon');actuator=E.SubElement(mj,'actuator');muscle_meta=[]
 for b in muscles+receptors:
  nameb=b.findtext('Name');ids=[a.text for a in b.findall('Attachments/AttachID')];tn='route_'+safe(nameb)
  t=E.SubElement(tendon,'spatial',name=tn,width='.0015',rgba='.1 .65 .2 .7')
  for aid in ids:E.SubElement(t,'site',site=sites[aid])
  if b.findtext('Type')=='LinearHillMuscle':E.SubElement(actuator,'general',name=safe(nameb),tendon=tn,gear='-1',gainprm='1',gaintype='fixed',biastype='none',dyntype='none',ctrllimited='true',ctrlrange=fmt([0,number(b,'MaximumTension')]))
  lt=b.find('LengthTension');st=b.find('StimulusTension');muscle_meta.append(dict(name=nameb,id=b.findtext('ID'),type=b.findtext('Type'),attachment_ids=ids,initial_route_length=sum(float(np.linalg.norm(allb[x]['wp']-allb[y]['wp'])) for x,y in zip(ids[:-1],ids[1:])),maximum_tension=number(b,'MaximumTension'),kse=number(b,'Kse'),kpe=number(b,'Kpe'),B=number(b,'B'),resting_length=number(lt,'RestingLength'),lwidth=number(lt,'Lwidth'),pe_length=number(lt,'PeLength'),min_pe_length=number(lt,'MinPeLength'),ib_discharge=number(b,'IbDischarge',100),ia_discharge=number(b,'IaDischarge',100),ii_discharge=number(b,'IIDischarge',100),stimulus=[number(st,k) for k in 'ABCD']))
 for b in springs:
  ids=[a.text for a in b.findall('Attachments/AttachID')];t=E.SubElement(tendon,'spatial',name='spring_'+safe(b.findtext('Name')),stiffness=str(number(b,'Stiffness')),damping=str(number(b,'Damping')),springlength=str(number(b,'NaturalLength')))
  for aid in ids:E.SubElement(t,'site',site=sites[aid])
 # Explicit source exclusions, including sensor-parent overlaps.
 for stim in root.findall('.//ExternalStimuli/Stimulus'):
  if stim.findtext('Type')=='MotorVelocity':
   source_joint=next((x for x in joints if x['source_id']==stim.findtext('JointID')),None)
   if source_joint:E.SubElement(actuator,'velocity',name='motor_'+safe(source_joint['name']),joint=safe(source_joint['name']),kv='10000',forcelimited='true',forcerange='-10000 10000')
 contact=E.SubElement(mj,'contact');seen=set()
 for pair in root.findall('.//CollisionExclusionPair'):
  refs=[x.text for x in pair if x.text in allb]
  if len(refs)==2 and all('mjbody' in allb[x] for x in refs):
   names=tuple(sorted(allb[x]['mjbody'] for x in refs))
   if names not in seen:E.SubElement(contact,'exclude',body1=names[0],body2=names[1]);seen.add(names)
 path=dest/(name+'.xml');E.indent(mj);E.ElementTree(mj).write(path,encoding='utf-8',xml_declaration=True)
 model=mujoco.MjModel.from_xml_path(str(path));data=mujoco.MjData(model);mujoco.mj_forward(model,data)
 lengths=[float(data.ten_length[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_TENDON,'route_'+safe(x['name']))]) for x in muscle_meta]
 errors=[abs(a-x['initial_route_length']) for a,x in zip(lengths,muscle_meta)]
 body_errors=[]
 for b in allb.values():
  if 'mjbody' in b:
   idx=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,b['mjbody']);body_errors.append(float(np.max(abs(data.xpos[idx]-b['wp']))))
 report=dict(machine='easteregg2',source=str(source),source_sha256=sha,source_snapshot=str(snap),xml=str(path),axes='x anterior, y superior, z right; gravity [0,-9.81,0]',joints=joints,muscles=muscle_meta,geometries=geometry,total_dynamic_mass_kg=mass,compiled_dynamic_mass_kg=float(np.sum(model.body_mass)),nq=model.nq,nv=model.nv,nu=model.nu,route_initial_length_max_error_m=max(errors),body_initial_position_max_error_m=max(body_errors),warnings=warnings,deviations=['Vortex contact law replaced with explicit MuJoCo friction/soft contact; coefficients are engineering assumptions except source rail sliding friction .05.','Contact-sensor overlay boxes are non-force-producing visual geometry; contact feedback must be computed from physical foot contact.','Hill force update is external, not native MuJoCo muscle; no neural/controller equivalence is claimed by mechanical import.','Anatomical knee neutral and sign must be audited before training; source joint coordinates are displacements from initialized pose.'])
 report['frame_adapter']=dict(source='+X anterior, +Y superior, +Z source-left',target='+X anterior, +Y superior, +Z right',matrix=FRAME.tolist(),rotation_rule='Rtarget=S Rsource S',hinge_axis_rule='Axial-vector det(S) S a; native source hinge axis is negative joint-frame X')
 report['source_gravity_m_s2']=number(root.find('.//Environment'),'Gravity');report['selected_gravity_m_s2']=-9.81
 report['deviations'].extend(['Source gravity −9.8m/s² uses the requested target−9.81m/s²; direction remains superior-axis negative.','Finite treadmill velocity servo kv10000Ns/m and±10000N replaces Vortex velocity motor; source drive .1m/s in native positive slide coordinate maps to backward world X.','Physics target2ms versus native1ms; neural substeps retain native.2ms.'])
 (dest/'IMPORT_AUDIT.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:report[k] for k in ('xml','nq','nv','nu','route_initial_length_max_error_m','body_initial_position_max_error_m','total_dynamic_mass_kg','compiled_dynamic_mass_kg','warnings')},indent=2))
 assert max(errors)<1e-10 and max(body_errors)<1e-10
 return path
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--newer',action='store_true');p.add_argument('--source');a=p.parse_args();build(a.newer,a.source)
