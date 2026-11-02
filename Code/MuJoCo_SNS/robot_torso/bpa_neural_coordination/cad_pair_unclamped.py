"""Explicit paired recruitment normalization; pair POST clamp is discarded."""
import json,hashlib,argparse,sys,shutil
from pathlib import Path
import numpy as np
import cad_pair_pf_audit as A
import cad_mixed_pressure_followup as K
C=A.C;J=A.J;L=A.L;OUT=A.OUT/'normalized_unclamped';HERE=Path(__file__).resolve().parent

def weights():
 report=json.loads((A.OUT/'STANCE_PATH_AUDIT.json').read_text());gain=report['conditional_PF_and_IBEXC_common_gain'];return dict(PF_E1=.2*gain,PF_E2=.3*gain,IBEXC_knee_ext=.6*gain)

def numerical_state(inner):
 d={}
 for prefix,obj in [('core_',inner.c),('sparse_',inner)]:
  for key,value in obj.__dict__.items():
   if isinstance(value,np.ndarray) and value.dtype.kind!='O':d[prefix+key]=value.copy()
   elif isinstance(value,(bool,int,float)):d[prefix+key]=np.asarray(value)
 return d

class FreePair:
 def __init__(self,inner,names,ports,enabled=False):
  self.inner=inner;self.enabled=bool(enabled);self.names=list(names);self.ids=np.array([names.index('MN_vas_med_'+s) for s in ('r','l')]);self.ports=np.array([ports.index('POST_vas_med_'+s) for s in ('r','l')]);self.buffer=np.zeros(len(ports));self.steps=0;self.logs=[];self.reset_done=False;self.absolute=weights();self.edits=[]
  assert not inner.c.spiking and not np.any(inner.c.g_ion[:,self.ids]),'PairMNs must be ordinary graded neurons'
  for post,port in zip(self.ids,self.ports):assert inner.c.input_connectivity[post,port]==1 and np.count_nonzero(inner.c.input_connectivity[:,port])==1
  if self.enabled:
   for side,post in zip(('r','l'),self.ids):
    for prefix,value in self.absolute.items():
     pre=names.index(prefix+'_'+side);expected={'PF_E1':.2,'PF_E2':.3,'IBEXC_knee_ext':.6}[prefix];assert inner.c.g_max_non[post,pre]==expected
     hit=(inner.post==post)&(inner.pre==pre);assert np.sum(hit)==1
     self.edits.append(dict(source=names[pre],target=names[post],original_uS=expected,absolute_uS=value));inner.c.g_max_non[post,pre]=value;inner.gmax[hit]=value
 def __getattr__(self,key):return getattr(self.inner,key)
 def forward(self,x):
  if not self.enabled:result=self.inner.forward(x);self.steps+=1;return result
  c=self.inner.c
  if not self.reset_done:c.V[self.ids]=c.V_rest[self.ids];self.reset_done=True
  np.copyto(self.buffer,x);requested=self.buffer[self.ports].copy();self.buffer[self.ports]=0.;before=c.V[self.ids].copy();gmax=c.g_max_non[self.ids];lo=c.e_lo[self.ids];hi=c.e_hi[self.ids];g=np.clip(gmax*(c.V[None,:]-lo)/np.where(hi==lo,1.,hi-lo),0,gmax);syn=np.sum(g*c.del_e[self.ids],axis=1)-before*np.sum(g,axis=1);leak=-c.g_m[self.ids]*(before-c.V_rest[self.ids]);tonic=c.i_b[self.ids].copy();app=c.input_connectivity[self.ids]@self.buffer;expected=before+c.time_factor_membrane[self.ids]*(syn+leak+tonic+app);result=self.inner.forward(self.buffer);after=c.V[self.ids].copy();error=after-expected
  if np.max(np.abs(error))>1e-12:raise ValueError('Pair differential equation mismatch')
  self.logs.append(np.column_stack((before,after,requested,np.zeros(2),syn,leak,tonic,app,error)));self.steps+=1;return result
 def export_state(self):
  d=numerical_state(self.inner);d.update(pair_parameters=np.array([self.enabled,*self.absolute.values()]),pair_ids=self.ids.copy(),pair_ports=self.ports.copy(),pair_buffer=self.buffer.copy(),pair_steps=np.array(self.steps),pair_reset_done=np.array(self.reset_done),pair_logs=np.asarray(self.logs).reshape(-1,2,9),pair_rng_state=np.empty(0,dtype=np.uint8));return d
 def import_state(self,d):
  assert np.array_equal(d['pair_parameters'],[self.enabled,*self.absolute.values()]) and len(d['pair_rng_state'])==0
  for prefix,obj in [('core_',self.inner.c),('sparse_',self.inner)]:
   for key,value in d.items():
    if key.startswith(prefix):
     name=key[len(prefix):];old=getattr(obj,name);setattr(obj,name,value.copy() if isinstance(old,np.ndarray) else type(old)(value))
  self.buffer[:]=d['pair_buffer'];self.steps=int(d['pair_steps']);self.reset_done=bool(d['pair_reset_done']);self.logs=[r.copy() for r in d['pair_logs']]

def network():
 R=C.R;m=C.mujoco.MjModel.from_xml_path(str(C.ROBOT/'ah001_neural.xml'));acts=[m.actuator(i).name for i in range(m.nu)];R.params.G.update(c1_gain=.15,v3_gain=0.,heel_rge=.15,toe_rge=.05);R.params.TAU['rg_nap_h']=.35;n=R.bn.build(acts,dt=m.opt.timestep,interleg=True);return n,sorted(n.idx,key=n.idx.get)

def checks():
 OUT.mkdir(exist_ok=True)
 with np.load(A.OUT/'TEACHER_NEURAL.npz') as z:inputs=z['inputs'][:96]
 n1,names=network();n2,_=network();old=C.R.ExactFastForward(n1.compiled);off=FreePair(C.R.ExactFastForward(n2.compiled),names,n2.inputs,False);default=True
 for x in inputs[:64]:default=bool(default and np.array_equal(old.forward(x),off.forward(x)))
 default=bool(default and all(np.array_equal(v,numerical_state(off.inner)[k]) for k,v in numerical_state(old).items()))
 n3,_=network();on=FreePair(C.R.ExactFastForward(n3.compiled),names,n3.inputs,True)
 for x in inputs[:64]:on.forward(x)
 saved=on.export_state();first=[on.forward(x).copy() for x in inputs[64:96]];end=on.export_state();on.import_state(saved);second=[on.forward(x).copy() for x in inputs[64:96]];again=on.export_state();state=bool(np.array_equal(first,second) and all(np.array_equal(end[k],again[k]) for k in end));logs=np.asarray(on.logs);fit=json.loads((A.OUT/'STANCE_PATH_AUDIT.json').read_text());r=dict(default_off_all_outputs_and_core_arrays_exact=default,enabled_continuation_all_states_exact=state,actual_pair_POST_current_max_nA=float(np.max(np.abs(logs[:,:,3]))),MN_equation_max_error_mV=float(np.max(np.abs(logs[:,:,8]))),pair_initialized_at_rest_not_target=on.reset_done,absolute_shared_weights_uS=on.absolute,edited_connections=on.edits,offline_pressure_domain_valid=fit['conditional_pressure_failure'] is None,offline_pressure_peak_kpa=fit['conditional_pressure_peak_kpa'],offline_force_RMS_error_N=fit['conditional_force_RMS_error_N'],passed=default and state and np.max(np.abs(logs[:,:,8]))<=1e-12 and fit['conditional_pressure_failure'] is None and fit['conditional_pressure_peak_kpa']<620 and max(fit['conditional_force_RMS_error_N'])<45.3,scope='Defaultoff electricalparity on actualteacherinputs; enabledfullstatecontinuation in SNSkernel, not wholeCADrunnerrestart. PairdirectPOSTzero, no synaptic cancellation or instantaneous voltage target afterone restinitialization.');(OUT/'OFFLINE_AND_STATE_GATE.json').write_text(json.dumps(r,indent=2));return r

def run(mode):
 gate=checks();assert gate['passed'];OUT.mkdir(exist_ok=True)
 for name in ('ASSUMED_PAIR_CONFIGURATION.json','STATIC_AUDIT.json'):shutil.copyfile(A.OLD/name,OUT/name)
 (OUT/'source_snapshot').mkdir(exist_ok=True);C.OUT=OUT;L.OUT=OUT;oldbuild=C.R.bn.build;holder=[];neurons=[]
 def build(*args,**kwargs):
  n=oldbuild(*args,**kwargs);neurons.append(n);names=sorted(n.idx,key=n.idx.get);step=n.step
  def forward(x):
   if not holder:wrapper=FreePair(n.compiled,names,n.inputs,True);n.compiled=wrapper;holder.append(wrapper)
   return step(x)
  n.step=forward;return n
 C.R.bn.build=build;oldcompile=L.compile_opt_in
 def compiled():
  fn,parser,code,edits=oldcompile();parser.set_defaults(activation_floor=.05);return fn,parser,code,edits
 L.compile_opt_in=compiled;protocol=dict(machine='easteregg2',phase='Explicit selectedpair recruitment normalization',absolute_shared_weights_uS=weights(),PF_and_IBEXC_only=True,other52_original_hybrid=True,pair_POST_clamp_applied=False,pair_target_initialization=False,actual_pair_activation_floor=False,teaching_target_normalized=.05,duration_s=2.,nominal_support_bw=.9,other_parameters_frozen=True,known_nominal_only_not_validation=True,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(C.__file__),Path(L.__file__),A.OUT/'TEACHER_NEURAL.npz',A.OUT/'STANCE_PATH_AUDIT.json')});path=OUT/'RECRUITMENT_PREREGISTRATION.json'
 if path.exists():assert json.loads(path.read_text())==protocol,'Frozenprotocolchanged'
 else:path.write_text(json.dumps(protocol,indent=2))
 shutil.copyfile(Path(__file__),OUT/'source_snapshot/cad_pair_unclamped.py')
 summary=L.run(mode);wrapper=holder[0];folder=OUT/mode;np.savez_compressed(folder/'PAIR_NEURAL_STATE.npz',**wrapper.export_state());nlog=np.asarray(wrapper.logs);np.savez_compressed(folder/'PAIR_NEURAL_LOG.npz',samples=nlog,channels=np.array(['Vbefore','Vafter','discarded_external_POST','applied_external_POST','synaptic_current','leak','intrinsic_bias','applied_input','differential_error']),physical_time=np.maximum(0,(np.arange(len(nlog))-100)*.002),sample_is_warm=np.arange(len(nlog))<100);summary.update(scope=('Two freeSNS motorneurons drive pressureBPAs; other52 retainoriginalhybridreferencecontext' if mode=='mixed' else 'Two freeSNS motorneurons drive referenceforce actuators; other52 retainoriginalhybridreferencecontext'),free_pair=True,actual_POST_current_max_nA=float(np.max(np.abs(nlog[:,:,3]))),equation_error_max_mV=float(np.max(np.abs(nlog[:,:,8]))),pair_voltage_range_mV=np.column_stack((nlog[:,:,1].min(axis=0),nlog[:,:,1].max(axis=0))).tolist(),physical_pair_voltage_range_mV=np.column_stack((nlog[100:,:,1].min(axis=0),nlog[100:,:,1].max(axis=0))).tolist(),pair_current_clamp_removed=True,actual_pair_floor_imposed=False);(folder/'UNCLAMPED_SUMMARY.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));return summary

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['reference','mixed']);a=p.parse_args()
 if a.mode:run(a.mode)
 else:print(json.dumps(checks(),indent=2))
