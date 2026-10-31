"""Opt-in physical deadband guides and read-only whole-robot measurements."""
import math
import xml.etree.ElementTree as E
import numpy as np
from import_virtual_walker import mujoco,fmt

def add_guides(root,clearance):
 assert .001<=clearance<=.1
 world=root.find('worldbody');pelvis=world.find("body[@name='Root']");rootp=np.array(list(map(float,pelvis.get('pos').split())));left=world.find("body[@name='Left_support_rail']");railp=np.array(list(map(float,left.get('pos').split())));rail_top=railp[1]+.02
 carriage=E.SubElement(world,'body',name='GuideCarriage',pos=fmt([rootp[0],0,0]));E.SubElement(carriage,'joint',name='GuideCarriageSlide',type='slide',axis='1 0 0',limited='true',range='-1.9 1.9',damping='0',solreflimit='.004 1',solimplimit='.999 .9999 .001')
 E.SubElement(carriage,'inertial',pos=fmt([0,rail_top+.15,0]),mass='2',diaginertia='.14 .2 .08')
 bottom=rail_top+.005;height=.3;center=bottom+height/2;roof_inner=rootp[1]-.04+.0125+clearance;half_x=.06+clearance+.02
 def box(name,p,size):E.SubElement(carriage,'geom',name='guide_'+name,type='box',pos=fmt(p),size=fmt(size),friction='.05 .005 .0001',rgba='.55 .6 .68 .4',solref='.006 1',solimp='.95 .999 .001',condim='3')
 for side in (-1,1):
  box('side_'+str(side),[0,center,side*(.4+clearance+.01)],[half_x,height/2,.01])
  box('cap_'+str(side),[0,roof_inner+.01,side*.3],[half_x,.01,.1])
  for front in (-1,1):box('post_'+str(side)+'_'+str(front),[front*(.06+clearance+.01),center,side*.3],[.01,height/2,.1])
 return dict(side_foreaft_clearance_at_neutral_m=clearance,initial_vertical_cap_clearance_m=clearance,cap_inner_height_m=float(roof_inner),post_bottom_height_m=float(bottom),carriage_mass_kg=2.,carriage_constraint='Passive unpowered frictionless x slider, finite±1.9m; apparatus held in y/z/rotation, robot remains6DoF and interacts only by contact.',source_rail_span_unchanged=True,vertical_gap_note='Initial wings are above horizontal rails. Cap gap grows as robot settles; measured supported-state gap must be reported, not equated to side clearance.')

class Measurements:
 shapes={'pelvis_quaternion':(4,),'pelvis_euler_XYZ_deg':(3,),'pelvis_omega_world':(3,),'pelvis_rotational_energy_J':(), 'com_position':(3,),'com_velocity':(3,),'rail_force_world':(3,),'rail_torque_about_pelvis':(3,),'guide_force_world':(3,),'guide_torque_about_pelvis':(3,),'guide_contact_count':(),'foot_force_world':(2,3),'foot_region_force_world':(4,3),'foot_cop_world':(2,3),'foot_body_position':(2,3),'guide_carriage_x':(),'guide_carriage_vx':(),'guide_carriage_kinetic_J':()}
 def __init__(self,m,rootid):
  self.m=m;self.rootid=rootid;self.robot=set()
  for bid in range(1,m.nbody):
   p=bid
   while p and p!=rootid:p=int(m.body_parentid[p])
   if p==rootid:self.robot.add(bid)
  self.rails={g for g in range(m.ngeom) if 'support_rail' in (mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g) or '')};self.guides={g for g in range(m.ngeom) if (mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g) or '').startswith('guide_')};self.ground={mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,b) for b in ['Treadmill_Belt','WalkingPath','GroundPlane']};self.feet={mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,b):i for i,b in enumerate(['foot_L','toe_L','foot_R','toe_R'])};self.carriage=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'GuideCarriageSlide');self.force=np.empty(6);self.velocity=np.empty(6)
  self.guides|={g for g in range(m.ngeom) if (mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g) or '').startswith('geom_Guide_')}
  if self.carriage<0:self.carriage=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'Passive_guide_X_slider')
 def sample(self,d):
  m=self.m;r={k:np.zeros(s) for k,s in self.shapes.items()};mujoco.mj_objectVelocity(m,d,mujoco.mjtObj.mjOBJ_BODY,self.rootid,self.velocity,0);omega=self.velocity[:3];R=d.xmat[self.rootid].reshape(3,3);Ri=d.ximat[self.rootid].reshape(3,3);principal_omega=Ri.T@omega
  r['pelvis_quaternion']=d.xquat[self.rootid].copy();r['pelvis_euler_XYZ_deg']=np.degrees([math.atan2(R[2,1],R[2,2]),math.asin(np.clip(-R[2,0],-1,1)),math.atan2(R[1,0],R[0,0])]);r['pelvis_omega_world']=omega.copy();r['pelvis_rotational_energy_J']=.5*float(np.sum(m.body_inertia[self.rootid]*principal_omega**2));r['com_position']=d.subtree_com[self.rootid].copy();mujoco.mj_subtreeVel(m,d);r['com_velocity']=d.subtree_linvel[self.rootid].copy();cop_sum=np.zeros((2,3));cop_weight=np.zeros(2)
  for c in range(d.ncon):
   contact=d.contact[c];g1,g2=int(contact.geom1),int(contact.geom2);b1,b2=int(m.geom_bodyid[g1]),int(m.geom_bodyid[g2]);mujoco.mj_contactForce(m,d,c,self.force);force=contact.frame.reshape(3,3).T@self.force[:3];couple=contact.frame.reshape(3,3).T@self.force[3:]
   # Contact frame normal points geom1→geom2; positive force acts on geom2.
   if (g1 in self.rails or g2 in self.rails) and (b1 in self.robot or b2 in self.robot):
    f=force if b2 in self.robot else -force;torque=couple if b2 in self.robot else -couple;r['rail_force_world']+=f;r['rail_torque_about_pelvis']+=np.cross(contact.pos-d.xipos[self.rootid],f)+torque
   if (g1 in self.guides or g2 in self.guides) and (b1 in self.robot or b2 in self.robot):
    f=force if b2 in self.robot else -force;torque=couple if b2 in self.robot else -couple;r['guide_force_world']+=f;r['guide_torque_about_pelvis']+=np.cross(contact.pos-d.xipos[self.rootid],f)+torque;r['guide_contact_count']+=1
   foot=b1 if b1 in self.feet else b2 if b2 in self.feet else None
   if foot is not None and ((b1 in self.ground) or (b2 in self.ground)):
    index=self.feet[foot];side=index//2;f=force if foot==b2 else -force;r['foot_region_force_world'][index]+=f;r['foot_force_world'][side]+=f;positive=max(0.,f[1]);cop_sum[side]+=contact.pos*positive;cop_weight[side]+=positive
  for side in range(2):r['foot_cop_world'][side]=cop_sum[side]/cop_weight[side] if cop_weight[side]>1e-9 else np.full(3,np.nan)
  for body,index in self.feet.items():
   if index in (0,2):r['foot_body_position'][index//2]=d.xpos[body]
  if self.carriage>=0:
   q=m.jnt_qposadr[self.carriage];v=m.jnt_dofadr[self.carriage];direction=float(m.jnt_axis[self.carriage,0]);r['guide_carriage_x']=float(d.qpos[q]*direction);r['guide_carriage_vx']=float(d.qvel[v]*direction);r['guide_carriage_kinetic_J']=float(d.qvel[v]**2)
  return r
