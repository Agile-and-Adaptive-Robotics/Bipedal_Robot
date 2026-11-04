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

SPARSE_CURRENT_STATUS={0:'not_requested',1:'valid_actual_sparse_step',
 2:'configured_backend_not_sparse',3:'not_the_actual_sparse_forward_binding',
 4:'actual_applied_buffer_unavailable_for_non_onehot_input',
 5:'no_actual_sparse_step_or_stale_V_last_buffer',6:'nonfinite_current_fields',
 7:'actual_step_update_does_not_close'}
SPARSE_CURRENT_FIELDS=('neural_V_last_mV','neural_applied_current_nA',
 'neural_chemical_current_nA','neural_leak_current_nA','neural_bias_current_nA',
 'neural_ionic_current_nA','neural_synaptic_conductance_uS','neural_synaptic_gE_nA',
 'neural_membrane_conductance_uS','neural_membrane_capacitance_numeric',
 'neural_rest_mV','neural_time_factor_membrane','neural_update_residual_mV')


def sparse_current_snapshot(network,indices,backend=None):
 """Observe the last actual sparse step, never reconstruct a caller/stimulus.

 Native overrides retain stale sparse workspaces. Exact binding, V_last alias,
 and the measured update closure are therefore required, not just class labels.
 """
 if backend is None:return dict(status=0,available=False,fields=None)
 if backend not in ('sparse','compact'):return dict(status=2,available=False,fields=None)
 from sparse_pairwise_eval import SparsePairwiseForward
 expected_forward=SparsePairwiseForward.forward
 if backend=='compact':
  from compact_pairwise_eval import CompactPairwiseForward
  expected_forward=CompactPairwiseForward.forward
 w=network.compiled;forward=w.forward;forward_scope='direct_sparse_forward'
 if backend=='compact':forward_scope='direct_compact_forward'
 def actual_sparse_binding(method):
  return getattr(method,'__self__',None)is w and getattr(method,'__func__',None)is expected_forward
 recognized=actual_sparse_binding(forward)
 if isinstance(w,SparsePairwiseForward)and not recognized and getattr(forward,'__closure__',None):
  from swing_diagnostics import SwingDiagnostics
  # Recognize the actual audited recorder, not an arbitrary function with a
  # convenient name. Its captured native must still be this sparse instance.
  if any(forward.__code__ is code for code in SwingDiagnostics.__init__.__code__.co_consts):
   closure={name:cell.cell_contents for name,cell in zip(forward.__code__.co_freevars,forward.__closure__)}
   recognized=(closure.get('c')is w and closure.get('net')is network and
               getattr(closure.get('self'),'net',None)is network and actual_sparse_binding(closure.get('native')))
   if recognized:forward_scope='audited_SwingDiagnostics_delegating_to_actual_sparse'
 if not isinstance(w,SparsePairwiseForward)or not recognized:
  return dict(status=3,available=False,fields=None)
 if not w.input_onehot:return dict(status=4,available=False,fields=None)
 c=w.c
 if not np.shares_memory(c.V_last,w.last_voltage):return dict(status=5,available=False,fields=None)
 ids=np.asarray(indices,dtype=int);last=c.V_last[ids].copy()
 applied=w.i_app[ids].copy();sg=w.sums[ids,0].copy();sge=w.sums[ids,1].copy()
 chemical=sge-last*sg;gm=c.g_m[ids].copy();rest=c.V_rest[ids].copy()
 leak=-gm*(last-rest);bias=c.i_b[ids].copy();factor=c.time_factor_membrane[ids].copy()
 ionic=np.zeros(len(ids))
 if c.gated:
  # Algebraic evaluation from the actual post-update gates, no gate integration.
  a=1/(1+c.k_a[:,ids]*np.exp(c.slope_a[:,ids]*(c.e_a[:,ids]-last)))
  channels=c.g_ion[:,ids]*(a**c.pow_a[:,ids])*(c.b_gate[:,ids]**c.pow_b[:,ids])*(c.c_gate[:,ids]**c.pow_c[:,ids])*(c.e_ion[:,ids]-last)
  ionic=np.sum(channels,axis=0)
 predicted=last+factor*(leak+bias+chemical+applied+ionic)
 residual=c.V[ids]-predicted
 fields=dict(zip(SPARSE_CURRENT_FIELDS,(last,applied,chemical,leak,bias,ionic,sg,sge,gm,c.c_m[ids].copy(),rest,factor,residual)))
 if not all(np.isfinite(value).all()for value in fields.values()):return dict(status=6,available=False,fields=None)
 if np.any(abs(residual)>1e-11):return dict(status=7,available=False,fields=None,maximum_closure_error_mV=float(abs(residual).max()))
 return dict(status=1,available=True,fields=fields,dt_s=float(c.dt),forward_scope=forward_scope,maximum_closure_error_mV=float(abs(residual).max())if len(ids)else 0.)

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
  self.time=np.empty(N);self.neural_time=np.empty(N);self.neural_input_time=np.empty(N)
  self.physics_clock=np.empty(N);self.force_interval_start=np.empty(N);self.voltage=np.empty((N,len(self.cells)))
  self.qpos=np.empty((N,model.nq));self.qvel=np.empty((N,model.nv))
  self.activation=np.empty((N,model.na));self.control=np.empty((N,model.nu))
  self.foot_force=np.empty((N,2,3));self.regional_vertical=np.empty((N,2,2))
  self.harness=np.empty((N,6));self.support=np.empty((N,3));self.root=np.empty((N,3))
  self.stimulus=np.full((N,len(self.cells)),np.nan);self.force=np.empty(6)
  self.load_work=np.zeros((2,3));self.region_work=np.zeros((2,2));self.wrench_work=np.zeros(6)
  self.current_buffers=None;self.current_status=None
  self.metadata={'schema':'aarl_aligned_verification_v1','native_frame':frame,
   'display_frame':'x anterior, y superior, z right','native_to_anatomical':self.rotation.tolist(),
   'sample_phase':sample_phase,'voltage_label':voltage_label,'side_order':['r','l'],
   'region_order':['hindfoot','forefoot'],'region_split':self.region_split,
   'regional_sensor_definition':'Physical ground contact locations projected from calcaneus toward toes; virtual regions, not calibrated heel/toe pads.',
   'force_definition':'Once-counted force on each foot, transformed from contact frame into anatomical world axes. Regional trace is true superior component, not normal-force proxy.',
   'harness_definition':'Pelvis generalized passive-force wrench projected on solver-stage root-joint axes unless caller supplies apparatus wrench; confirm only apparatus spring/damper acts on those coordinates. Explicit support kept separate.',
   'kinematic_evaluation_phase':'Contact frames/body xpos are from the latest solver evaluation; qpos/qvel are final state. No extra forward call is made to retime those intermediate quantities.',
   'missing_requested_cells':self.missing,'provenance':provenance or {},
   'stimulus_policy':'NaN means not logged; never synthesize an inferred pulse.',
   'no_model_or_neural_steps':True,'capture_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

 def sample(self,data,state_time_s=None,neural_time_s=None,harness_wrench_native=None,
            explicit_support_native=None,stimulus_nA=None,neural_input_time_s=None,
            force_interval_start_s=None,sparse_current_backend=None):
  k=self.count
  if k>=self.capacity:raise RuntimeError('Verification capture capacity exceeded; refusing silent truncation')
  t=float(data.time if state_time_s is None else state_time_s)
  if k and t<=self.time[k-1]:raise ValueError('State timestamps must increase monotonically')
  self.time[k]=t;self.neural_time[k]=t if neural_time_s is None else neural_time_s
  self.neural_input_time[k]=t if neural_input_time_s is None else neural_input_time_s
  self.physics_clock[k]=data.time
  self.force_interval_start[k]=t if force_interval_start_s is None else force_interval_start_s
  self.voltage[k]=self.network.compiled.V[self.indices]
  self.qpos[k]=data.qpos;self.qvel[k]=data.qvel;self.activation[k]=data.act;self.control[k]=data.ctrl
  self.root[k]=self.rotation@data.xpos[self.pelvis]
  loads=self.load_work;regions=self.region_work;loads[:]=0.;regions[:]=0.
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
   wrench=self.wrench_work;wrench[:]=0.
   for names,offset in ((('pelvis_tx','pelvis_ty','pelvis_tz'),0),
                        (('pelvis_tilt','pelvis_list','pelvis_rotation'),3)):
    for name in names:
     j=self.model.joint(name);vd=int(j.dofadr[0]);value=float(data.qfrc_passive[vd])
     wrench[offset:offset+3]+=value*data.xaxis[j.id]
  else:wrench=np.asarray(harness_wrench_native,dtype=float)
  assert wrench.shape==(6,)
  self.harness[k,:3]=self.rotation@wrench[:3];self.harness[k,3:]=self.rotation@wrench[3:]
  self.support[k]=np.nan if explicit_support_native is None else self.rotation@explicit_support_native
  if stimulus_nA is not None:
   stimulus=np.asarray(stimulus_nA)
   self.stimulus[k]=stimulus[self.indices]if len(stimulus)==len(self.network.compiled.V)else stimulus
  if sparse_current_backend is not None:
   if self.current_buffers is None:
    self.current_buffers={name:np.full((self.capacity,len(self.cells)),np.nan)for name in SPARSE_CURRENT_FIELDS}
    self.current_status=np.zeros(self.capacity,dtype=np.uint8)
    self.metadata['current_telemetry']={
     'schema':'aarl_optional_actual_sparse_currents_v1','status_codes':SPARSE_CURRENT_STATUS,
     'source':'Last actual configured sparse/compact i_app and sums with actual V_last; exact forward binding and update closure checked.',
     'applied_current':'All actually mapped external inputs, including direct POST when present; not the separately logged stimulus.',
     'chemical_current':'sums_gE - V_last*sums_g from the last real sparse forward call.',
     'ionic_current':'Algebraic evaluation using actual last voltage and already-updated gates; no new gate/neural update.',
     'update':'V_new=V_last+time_factor*(leak+bias+chemical+applied+ionic)',
     'unavailable':'NaN fields plus non-valid status; native/stale workspaces are never represented as zero.',
     'closure_atol_mV':1e-11,'clock':'Same recorded neural output time; V_last is the actual preceding neural state.',
     'activation_scope':'activation is actual MuJoCo data.act; control is data.ctrl. Legacy main NPZ act may contain ctrl and is not used to infer either currents or activation.'}
   observed=sparse_current_snapshot(self.network,self.indices,sparse_current_backend)
   self.current_status[k]=observed['status']
   if observed['available']:
    for name,value in observed['fields'].items():self.current_buffers[name][k]=value
    self.metadata['current_telemetry']['neural_step_dt_s']=observed['dt_s']
    self.metadata['current_telemetry']['forward_scope']=observed['forward_scope']
  self.count+=1

 def _saved_current_arrays(self,n):
  if self.current_buffers is None:return {}
  return {**{name:value[:n]for name,value in self.current_buffers.items()},
          'neural_current_status':self.current_status[:n],
          'neural_current_available':self.current_status[:n]==1}

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
   t_neural_input=self.neural_input_time[:n],physics_clock_s=self.physics_clock[:n],
   force_interval_start=self.force_interval_start[:n],force_interval_end=self.time[:n],
   neural_voltage=self.voltage[:n],neural_names=np.array(self.cells),qpos=self.qpos[:n],qvel=self.qvel[:n],
   activation=self.activation[:n],control=self.control[:n],joint_names=np.array(self.joints),joint_qpos_addresses=self.qadr,
   foot_force_anatomical_N=self.foot_force[:n],regional_vertical_N=self.regional_vertical[:n],
   harness_wrench_anatomical=self.harness[:n],explicit_support_anatomical_N=self.support[:n],
   root_anatomical_m=self.root[:n],stimulus_nA=self.stimulus[:n],**self._saved_current_arrays(n))
  (out/'ALIGNED_TRACE_METADATA.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
  return out/'aligned_trace.npz'
