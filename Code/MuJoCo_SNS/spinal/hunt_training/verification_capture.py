"""Read-only, preallocated same-run neural/mechanical verification telemetry.

This module never advances physics or neurons. Call sample at an explicitly
declared runner phase, and retain distinct state/neural timestamps when needed.
"""
from pathlib import Path
import json,hashlib
import numpy as np

FRAMES={
 'legacy_z_up':np.array([[1.,0.,0.],[0.,0.,1.],[0.,-1.,0.]]),
 'anatomical_y_up':np.eye(3),
}

def requested_cells():
 names=[]
 for s in ('r','l'):
  names += ['RG_E_'+s,'RG_F_'+s]
  names += [x+s for x in ('PF_BI_JOINT_RECT_','PF_BI_JOINT_HAM_',
                          'PF_IN_BI_PAIR_RECT_','PF_IN_BI_PAIR_HAM_')]
  names += ['MN_'+a+'_'+s for a in ('iliacus','bifemlh','vas_int','bifemsh','soleus','tib_ant')]
 return names

class AlignedCapture:
 def __init__(self,model,network,max_samples,frame='legacy_z_up',cells=None,
              sample_phase='caller must declare',voltage_label='model voltage (mV)',
              provenance=None,region_split=.5):
  import mujoco
  self.mj=mujoco;self.model=model;self.network=network;self.rotation=FRAMES[frame]
  self.capacity=int(max_samples);assert self.capacity>1;self.count=0
  requested=requested_cells()if cells is None else list(cells)
  self.cells=[n for n in requested if n in network.idx]
  self.missing=[n for n in requested if n not in network.idx]
  assert len(self.cells)==len(set(self.cells));self.indices=np.array([network.idx[n]for n in self.cells],int)
  self.joints=[j+'_'+s for s in ('r','l')for j in ('hip_flexion','knee_angle','ankle_angle')]
  self.qadr=np.array([int(model.joint(n).qposadr[0])for n in self.joints])
  self.calcn=[model.body('calcn_'+s).id for s in ('r','l')]
  self.toes=[model.body('toes_'+s).id for s in ('r','l')]
  self.owners={model.body(b+'_'+s).id:i for i,s in enumerate(('r','l'))for b in ('talus','calcn','toes')}
  self.pelvis=model.body('pelvis').id;self.region_split=float(region_split);assert 0<self.region_split<1
  N=self.capacity
  self.time=np.empty(N);self.neural_time=np.empty(N);self.voltage=np.empty((N,len(self.cells)))
  self.qpos=np.empty((N,model.nq));self.qvel=np.empty((N,model.nv))
  self.activation=np.empty((N,model.na));self.control=np.empty((N,model.nu))
  self.foot_force=np.empty((N,2,3));self.regional_vertical=np.empty((N,2,2))
  self.harness=np.empty((N,6));self.support=np.empty((N,3));self.root=np.empty((N,3))
  self.stimulus=np.full((N,len(self.cells)),np.nan);self.force=np.empty(6)
  self.metadata={'schema':'aarl_aligned_verification_v1','native_frame':frame,
   'display_frame':'x anterior, y superior, z right','native_to_anatomical':self.rotation.tolist(),
   'sample_phase':sample_phase,'voltage_label':voltage_label,'side_order':['r','l'],
   'region_order':['hindfoot','forefoot'],'region_split':self.region_split,
   'regional_sensor_definition':'Physical ground contact locations projected from calcaneus toward toes; virtual regions, not calibrated heel/toe pads.',
   'force_definition':'Once-counted force on each foot, transformed from contact frame into anatomical world axes. Regional trace is true superior component, not normal-force proxy.',
   'harness_definition':'Pelvis translational/rotational joint spring-damper wrench unless caller supplies actual apparatus wrench; explicit support kept separate.',
   'missing_requested_cells':self.missing,'provenance':provenance or {},
   'stimulus_policy':'NaN means not logged; never synthesize an inferred pulse.',
   'no_model_or_neural_steps':True}

 def sample(self,data,state_time_s=None,neural_time_s=None,harness_wrench_native=None,
            explicit_support_native=None,stimulus_nA=None):
  k=self.count
  if k>=self.capacity:raise RuntimeError('Verification capture capacity exceeded; refusing silent truncation')
  t=float(data.time if state_time_s is None else state_time_s)
  if k and t<=self.time[k-1]:raise ValueError('State timestamps must increase monotonically')
  self.time[k]=t;self.neural_time[k]=t if neural_time_s is None else neural_time_s
  self.voltage[k]=self.network.compiled.V[self.indices]
  self.qpos[k]=data.qpos;self.qvel[k]=data.qvel;self.activation[k]=data.act;self.control[k]=data.ctrl
  self.root[k]=self.rotation@data.xpos[self.pelvis]
  loads=np.zeros((2,3));regions=np.zeros((2,2))
  for ci in range(data.ncon):
   c=data.contact[ci];a,b=int(c.geom1),int(c.geom2)
   ga=self.model.geom_type[a]==self.mj.mjtGeom.mjGEOM_PLANE
   gb=self.model.geom_type[b]==self.mj.mjtGeom.mjGEOM_PLANE
   if ga==gb:continue
   body=int(self.model.geom_bodyid[b if ga else a])
   while body and body not in self.owners:body=int(self.model.body_parentid[body])
   if not body:continue
   side=self.owners[body];self.mj.mj_contactForce(self.model,data,ci,self.force)
   world=c.frame.reshape(3,3).T@self.force[:3]
   if gb:world=-world # mj_contactForce result acts on geom2; select force on foot.
   anatomical=self.rotation@world;loads[side]+=anatomical
   origin=data.xpos[self.calcn[side]];direction=data.xpos[self.toes[side]]-origin
   fraction=float((c.pos-origin)@direction/max(direction@direction,1e-12))
   region=int(body==self.toes[side] or fraction>=self.region_split)
   regions[side,region]+=anatomical[1]
  self.foot_force[k]=loads;self.regional_vertical[k]=regions
  if harness_wrench_native is None:
   wrench=np.zeros(6)
   for names,offset in ((('pelvis_tx','pelvis_ty','pelvis_tz'),0),
                        (('pelvis_tilt','pelvis_list','pelvis_rotation'),3)):
    for name in names:
     j=self.model.joint(name);ad=int(j.qposadr[0]);vd=int(j.dofadr[0])
     value=-self.model.jnt_stiffness[j.id]*(data.qpos[ad]-self.model.qpos_spring[ad])-self.model.dof_damping[vd]*data.qvel[vd]
     wrench[offset:offset+3]+=value*data.xaxis[j.id]
  else:wrench=np.asarray(harness_wrench_native,dtype=float)
  assert wrench.shape==(6,)
  self.harness[k,:3]=self.rotation@wrench[:3];self.harness[k,3:]=self.rotation@wrench[3:]
  self.support[k]=np.nan if explicit_support_native is None else self.rotation@explicit_support_native
  if stimulus_nA is not None:
   stimulus=np.asarray(stimulus_nA)
   self.stimulus[k]=stimulus[self.indices]if len(stimulus)==len(self.network.compiled.V)else stimulus
  self.count+=1

 def save(self,folder,model_filename='model.mjb'):
  out=Path(folder);out.mkdir(parents=True,exist_ok=True);n=self.count;assert n>1
  modelpath=out/model_filename
  self.mj.mj_saveModel(self.model,str(modelpath),None)
  meta={**self.metadata,'samples':n,'capacity':self.capacity,'body_weight_N':float(self.model.body_subtreemass[self.pelvis]*9.81),
        'body_weight_definition':'Pelvis subtree mass times9.81; excludes separate fixed ground/apparatus bodies.',
        'model_file':modelpath.name,'model_sha256':hashlib.sha256(modelpath.read_bytes()).hexdigest(),
        'actuator_names':[self.model.actuator(i).name for i in range(self.model.nu)],
        'mujoco_version':self.mj.__version__}
  np.savez_compressed(out/'aligned_trace.npz',t_state=self.time[:n],t_neural=self.neural_time[:n],
   neural_voltage=self.voltage[:n],neural_names=np.array(self.cells),qpos=self.qpos[:n],qvel=self.qvel[:n],
   activation=self.activation[:n],control=self.control[:n],joint_names=np.array(self.joints),joint_qpos_addresses=self.qadr,
   foot_force_anatomical_N=self.foot_force[:n],regional_vertical_N=self.regional_vertical[:n],
   harness_wrench_anatomical=self.harness[:n],explicit_support_anatomical_N=self.support[:n],
   root_anatomical_m=self.root[:n],stimulus_nA=self.stimulus[:n])
  (out/'ALIGNED_TRACE_METADATA.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
  return out/'aligned_trace.npz'
