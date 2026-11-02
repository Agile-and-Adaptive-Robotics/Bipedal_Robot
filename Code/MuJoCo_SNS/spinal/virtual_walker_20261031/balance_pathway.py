"""Explicit whole-body balance experiments; source baseline remains unchanged."""
import json,itertools
from pathlib import Path
import numpy as np
from numba import njit

@njit(cache=True)
def chemical_current_integral(v,edges,indices,integral,dt):
 for k in range(len(indices)):
  e=edges[indices[k]];pre=int(e[0]);post=int(e[1]);activity=min(1.,max(0.,(v[pre]+60.)/20.));integral[k]+=e[2]*activity*(e[5]-v[post])*dt

class BodyBalanceRelay:
 """Two actual LIF populations and phase-gated excitatory chemical pathways.

 Signed orientation error recruits opposite hip actions in stance and swing.
 Existing Ia cells provide the source circuit's reciprocal inhibition.
 """
 fields=('kp_nA_per_rad','kd_s','target_pitch_deg','delay_s','max_current_nA','mn_gmax_uS','ia_gmax_uS','swing_gain')
 def __init__(self,c,configuration):
  from spiking_sensory import encode_step
  self.c=c;self.cfg=configuration;assert set(configuration)-{'capture_weight','functional_projection','ankle_gain_ratio','capture_margin_m'}==set(self.fields);assert all(np.isfinite(v) for v in configuration.values());assert configuration['delay_s']>=0 and 0<configuration['max_current_nA']<=8
  assert configuration.get('functional_projection',0) in (0,1) and 0<=configuration.get('capture_weight',0)<=1
  assert all(configuration[k]>=0 for k in ('kp_nA_per_rad','kd_s','mn_gmax_uS','ia_gmax_uS','swing_gain'));self.encode=encode_step;self.indices=np.arange(len(c.v),len(c.v)+2);self.edge_specs=[]
  for key in ('v','newv','m','h','gs','gi','applied','tonic','last_test_impulse_nA_s'):
   old=getattr(c,key);setattr(c,key,np.r_[old,[-60.,-60.] if key in ('v','newv') else [0.,0.]])
  c.names.extend(['Body posterior balance relay','Body anterior balance relay']);c.params=np.vstack([c.params,np.tile([-60.,25.,1.,0.,0.,1.,1.,0.,1.,1.,0.],(2,1))])
  rows=[]
  for pop,sign in enumerate(('posterior','anterior')):
   for side in 'LR':
    for phase in ('stance','swing'):
     if configuration.get('functional_projection',0)==0:joint='Hip';action='flx' if (sign=='posterior')==(phase=='stance') else 'ext';scale=1.
     else:joint='Ank' if phase=='stance' else 'Knee';action=('ext' if sign=='posterior' else 'flx') if phase=='stance' else ('flx' if sign=='posterior' else 'ext');scale=configuration.get('ankle_gain_ratio',4.) if phase=='stance' else 1.
     for target,kind in [(side+' '+joint+' MN '+action,'MN'),(side+' '+joint+' '+action+' Ia','Ia')]:
      i=c.names.index(target);assert c.params[i,2]==1.;rows.append([self.indices[pop],i,0.,-60.,-40.,-10.]);self.edge_specs.append(dict(index=len(c.edges)+len(rows)-1,source=c.names[self.indices[pop]],target=target,side=side,phase=phase,kind=kind,source_label_role='Actual inhibitory Ia cell' if kind=='Ia' else 'Source-action-verified motor neuron',configured_gmax_uS=configuration['mn_gmax_uS' if kind=='MN' else 'ia_gmax_uS']*scale))
  c.edges=np.vstack([c.edges,rows]);self.v=np.tile(np.linspace(.01,.99,16),(2,1));self.ref=np.zeros((2,16));self.decoded=np.zeros(2);self.spikes=np.zeros((2,16),dtype=np.int8);self.counts=np.zeros(2,dtype=np.int64);self.total_counts=self.counts.copy();self.target=np.zeros(2);self.current=np.zeros(2);self.gates=np.zeros(2);self.tau=np.full(2,.025);self.delay=np.zeros(max(1,round(configuration['delay_s']/.002)));self.pointer=np.array(0,dtype=np.int64);self.measurement=np.zeros(4);self.integral=np.zeros(16);self.logs=[]
  self.encode(self.target,self.v.copy(),self.ref.copy(),self.decoded.copy(),self.spikes.copy(),self.counts.copy(),c.dt,200.,.005,self.tau,.001)
  self.edge_indices=np.array([s['index'] for s in self.edge_specs],dtype=np.int64);chemical_current_integral(c.v,c.edges,self.edge_indices,self.integral.copy(),c.dt)
 def observe(self,model,data,previous,time,interval):
  from import_virtual_walker import mujoco
  root=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Root');R=data.xmat[root].reshape(3,3);pitch=float(np.arctan2(R[0,1],R[1,1]));vel=np.zeros(6);mujoco.mj_objectVelocity(model,data,mujoco.mjtObj.mjOBJ_BODY,root,vel,0);rate=-float(vel[2]);forces=np.zeros(2) if previous is None else previous['foot_force_world'][:,1];mass=sum(model.body_mass[b] for b in range(1,model.nbody) if self._robot_descendant(model,b,root));bw=mass*9.81;desired=np.clip(forces/(.2*bw),0,1);self.gates+=(desired-self.gates)*(-np.expm1(-interval/.02))
  error=np.radians(self.cfg['target_pitch_deg'])-pitch-self.cfg['kd_s']*rate
  if self.cfg.get('capture_weight',0)>0:
   capture_error=0.
   if previous is not None and np.any(forces>.05*bw):
    belt=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Treadmill_Belt');mujoco.mj_objectVelocity(model,data,mujoco.mjtObj.mjOBJ_BODY,belt,vel,0);g=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,'geom_Treadmill_Belt');top=data.geom_xpos[g,1]+np.abs(data.geom_xmat[g].reshape(3,3)[1])@model.geom_size[g];height=max(.1,previous['com_position'][1]-top);capture=previous['com_position'][0]+(previous['com_velocity'][0]-vel[3])/np.sqrt(9.81/height);bounds=[]
    for k,name in enumerate(('foot_L','toe_L','foot_R','toe_R')):
     if forces[k//2]<=.05*bw or previous['foot_region_force_world'][k,1]<=.01*bw:continue
     g=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+name);extent=np.abs(data.geom_xmat[g].reshape(3,3)[0])@model.geom_size[g];bounds.extend([data.geom_xpos[g,0]-extent,data.geom_xpos[g,0]+extent])
    if bounds:
     margin=min(self.cfg.get('capture_margin_m',.005),(max(bounds)-min(bounds))*.2);target=np.clip(capture,min(bounds)+margin,max(bounds)-margin);capture_error=(target-capture)/height
   error=(1-self.cfg['capture_weight'])*error+self.cfg['capture_weight']*capture_error
  drive=self.cfg['max_current_nA']*np.tanh(self.cfg['kp_nA_per_rad']*error/self.cfg['max_current_nA']) if np.any(forces>.05*bw) else 0.;i=int(self.pointer);delayed=drive if self.cfg['delay_s']==0 else float(self.delay[i]);self.delay[i]=drive;self.pointer[...] = (i+1)%len(self.delay);self.current[:]=[max(0.,delayed),max(0.,-delayed)];self.target[:]=self.current/20.;self.measurement[:]=[pitch,rate,error,float(np.any(forces>.05*bw))]
  for spec in self.edge_specs:
   side=0 if spec['side']=='L' else 1;F=self.c.names.index(spec['side']+' RG flx');E=self.c.names.index(spec['side']+' RG ext');f=max(0.,self.c.v[F]-self.c.params[F,0]);e=max(0.,self.c.v[E]-self.c.params[E,0]);gate=self.gates[side] if spec['phase']=='stance' else (1-self.gates[side])*f/(f+e+1e-12)*self.cfg['swing_gain'];self.c.edges[spec['index'],2]=spec['configured_gmax_uS']*gate
 @staticmethod
 def _robot_descendant(model,body,root):
  while body and body!=root:body=int(model.body_parentid[body])
  return body==root
 def prepare(self):self.counts[:]=0;self.integral[:]=0
 def advance(self):
  self.encode(self.target,self.v,self.ref,self.decoded,self.spikes,self.counts,self.c.dt,200.,.005,self.tau,.001);self.proxy=-60.+20.*self.decoded;self.c.v[self.indices]=self.proxy
  chemical_current_integral(self.c.v,self.c.edges,self.edge_indices,self.integral,self.c.dt)
 def finish(self,time):
  self.total_counts+=self.counts;self.logs.append((time,self.counts.copy(),self.decoded.copy(),self.current.copy(),self.measurement.copy(),self.gates.copy(),self.integral.copy()))
 def export_state(self):
  state={key:getattr(self,key).copy() for key in ('v','ref','decoded','spikes','counts','total_counts','target','current','gates','delay','pointer','measurement','integral')}
  for i,key in enumerate(('time','counts','decoded','current','measurement','gates','integral')):state['log_'+key]=np.array([row[i] for row in self.logs])
  state['configuration']=np.array(json.dumps(self.cfg,sort_keys=True));state['proxy_present']=np.array(hasattr(self,'proxy'))
  if hasattr(self,'proxy'):state['proxy']=self.proxy.copy()
  return state
 def import_state(self,state):
  assert str(state['configuration'])==json.dumps(self.cfg,sort_keys=True)
  for key in ('v','ref','decoded','spikes','counts','total_counts','target','current','gates','delay','pointer','measurement','integral'):np.copyto(getattr(self,key),state[key])
  self.logs=[tuple(float(values[0]) if i==0 else values[i].copy() for i in range(7)) for values in zip(*(state['log_'+key] for key in ('time','counts','decoded','current','measurement','gates','integral')))]
  if bool(state['proxy_present']):self.proxy=state['proxy'].copy()
  elif hasattr(self,'proxy'):del self.proxy
 def save(self,folder):
  keys=('t','counts','decoded','input_nA','body_measurement','stance_gate','synaptic_integral_nA_s');np.savez_compressed(folder/'body_balance_spikes.npz',**{key:np.array([row[i] for row in self.logs]) for i,key in enumerate(keys)});(folder/'BODY_BALANCE_ARCHITECTURE.json').write_text(json.dumps(dict(configuration=self.cfg,edges=self.edge_specs,spike_counts=self.total_counts.tolist(),actual_added_spiking_cells=32,added_logical_populations=2,scope='New engineering balance architecture. Actual LIF populations drive chemical conductances; reciprocal inhibition uses inherited Ia cells. No joint torques or muscle force targets applied. Source graph itself unchanged; runtime graph adds2 proxy channels and16 excitatory edges. Configured cue uses pitch/rate or belt-relative extrapolated COM versus loaded support, with explicit delay; ideal simulator-state estimate, not physiological sensory fusion; foot force gate uses previous physics sample plus20ms filtering.'),indent=2))

class LoadTransferTerminationRelay:
 """Spike-encoded, force-qualified strikes briefly terminate opposite PF-E bursts."""
 learned=('termination_gmax_uS','pulse_s','on_load_BW')
 def __init__(self,c,configuration):
  from spiking_sensory import encode_step
  self.c=c;self.cfg=configuration;self.encode=encode_step;assert configuration['pathway_kind']=='contralateral_load_termination'
  assert set(configuration)-{'phase_conditioned'}==set(self.learned)|{'pathway_kind','off_load_BW','load_dwell_s','unload_dwell_s','refractory_s','input_current_nA'};assert configuration.get('phase_conditioned',0) in (0,1)
  # Nominal search limit2, with explicit headroom for learned-gain Â±1% validation.
  assert 0<=configuration['termination_gmax_uS']<=2.05 and .02<=configuration['pulse_s']<=.25 and 0<configuration['off_load_BW']<configuration['on_load_BW']<.8
  assert configuration['load_dwell_s']>=.02 and configuration['unload_dwell_s']>=.04 and configuration['refractory_s']>=.2 and 0<configuration['input_current_nA']<=4
  self.indices=np.arange(len(c.v),len(c.v)+2);self.edge_specs=[]
  for key in ('v','newv','m','h','gs','gi','applied','tonic','last_test_impulse_nA_s'):
   old=getattr(c,key);setattr(c,key,np.r_[old,[-60.,-60.] if key in ('v','newv') else [0.,0.]])
  c.names.extend(['L force-qualified load-transfer relay','R force-qualified load-transfer relay']);c.params=np.vstack([c.params,np.tile([-60.,25.,1.,0.,0.,1.,1.,0.,1.,1.,0.],(2,1))]);rows=[]
  for pop,source_side,target_side in [(0,'L','R'),(1,'R','L')]:
   for joint in ('Hip','Knee'):
    target=target_side+' '+joint+' PF flx IN';post=c.names.index(target);pf=c.names.index(target_side+' '+joint+' PF ext');down=[e for e in c.edges if int(e[0])==post and int(e[1])==pf];assert len(down)==1 and down[0][5]<-60
    rows.append([self.indices[pop],post,configuration['termination_gmax_uS'],-60.,-40.,-10.]);self.edge_specs.append(dict(index=len(c.edges)+len(rows)-1,source=c.names[self.indices[pop]],target=target,source_foot=source_side,opposite_leg=target_side,configured_gmax_uS=configuration['termination_gmax_uS'],reversal_mV=-10.,inherited_inhibitory_target=c.names[pf],inherited_inhibitory_reversal_mV=float(down[0][5])))
  c.edges=np.vstack([c.edges,rows]);self.v=np.tile(np.linspace(.01,.99,16),(2,1));self.ref=np.zeros((2,16));self.decoded=np.zeros(2);self.spikes=np.zeros((2,16),dtype=np.int8);self.counts=np.zeros(2,dtype=np.int64);self.total_counts=self.counts.copy();self.target=np.zeros(2);self.current=np.zeros(2);self.tau=np.full(2,.025);self.armed=np.zeros(2,dtype=np.bool_);self.high_dwell=np.zeros(2);self.low_dwell=np.zeros(2);self.pulse_until=np.full(2,-1.);self.last_event=np.full(2,-10.);self.event_counts=np.zeros(2,dtype=np.int64);self.measurement=np.zeros(8);self.phase_gate=np.ones(2);self.integral=np.zeros(4);self.logs=[];self.edge_indices=np.array([e['index'] for e in self.edge_specs],dtype=np.int64)
  self.phase_indices=np.array([[c.names.index(s+' RG flx'),c.names.index(s+' RG ext')] for s in ('R','L')],dtype=np.int64);self.encode(self.target,self.v.copy(),self.ref.copy(),self.decoded.copy(),self.spikes.copy(),self.counts.copy(),c.dt,200.,.005,self.tau,.001);chemical_current_integral(c.v,c.edges,self.edge_indices,self.integral.copy(),c.dt)
 def observe_forces(self,force_BW,time,interval):
  events=np.zeros(2)
  for side in range(2):
   load=max(0.,float(force_BW[side]));self.low_dwell[side]=self.low_dwell[side]+interval if load<self.cfg['off_load_BW'] else 0.
   if self.low_dwell[side]>=self.cfg['unload_dwell_s']:self.armed[side]=True
   self.high_dwell[side]=self.high_dwell[side]+interval if load>=self.cfg['on_load_BW'] and self.armed[side] else 0.
   if self.high_dwell[side]+1e-12>=self.cfg['load_dwell_s'] and time-self.last_event[side]>=self.cfg['refractory_s']:
    events[side]=1.;self.event_counts[side]+=1;self.last_event[side]=time;self.pulse_until[side]=time+self.cfg['pulse_s'];self.armed[side]=False;self.high_dwell[side]=0.
  self.current[:]=self.cfg['input_current_nA']*(time<self.pulse_until);self.target[:]=self.current/20.;self.measurement[:6]=np.r_[force_BW,events,np.maximum(0.,self.pulse_until-time)]
 def observe(self,model,data,previous,time,interval):
  from import_virtual_walker import mujoco
  root=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Root');bw=sum(model.body_mass[b] for b in range(1,model.nbody) if BodyBalanceRelay._robot_descendant(model,b,root))*9.81;forces=np.zeros(2) if previous is None else previous['foot_force_world'][:,1]/bw;self.observe_forces(forces,time,interval)
 def prepare(self):self.counts[:]=0;self.integral[:]=0
 def advance(self):
  for source in range(2):
   f,e=self.phase_indices[source];F=max(0.,self.c.v[f]-self.c.params[f,0]);E=max(0.,self.c.v[e]-self.c.params[e,0]);self.phase_gate[source]=F/(F+E+1e-12) if self.cfg.get('phase_conditioned',0) else 1.
  for spec in self.edge_specs:self.c.edges[spec['index'],2]=spec['configured_gmax_uS']*self.phase_gate[0 if spec['source_foot']=='L' else 1]
  self.measurement[6:]=self.phase_gate;self.encode(self.target,self.v,self.ref,self.decoded,self.spikes,self.counts,self.c.dt,200.,.005,self.tau,.001);self.proxy=-60.+20.*self.decoded;self.c.v[self.indices]=self.proxy;chemical_current_integral(self.c.v,self.c.edges,self.edge_indices,self.integral,self.c.dt)
 def finish(self,time):self.total_counts+=self.counts;self.logs.append((time,self.counts.copy(),self.decoded.copy(),self.current.copy(),self.measurement.copy(),self.armed.copy(),self.integral.copy()))
 def export_state(self):
  fields=('v','ref','decoded','spikes','counts','total_counts','target','current','armed','high_dwell','low_dwell','pulse_until','last_event','event_counts','measurement','phase_gate','integral');state={k:getattr(self,k).copy() for k in fields};state['configuration']=np.array(json.dumps(self.cfg,sort_keys=True));state['proxy_present']=np.array(hasattr(self,'proxy'))
  if hasattr(self,'proxy'):state['proxy']=self.proxy.copy()
  for i,k in enumerate(('time','counts','decoded','current','measurement','armed','integral')):state['log_'+k]=np.array([row[i] for row in self.logs])
  return state
 def import_state(self,state):
  assert str(state['configuration'])==json.dumps(self.cfg,sort_keys=True)
  for k in ('v','ref','decoded','spikes','counts','total_counts','target','current','armed','high_dwell','low_dwell','pulse_until','last_event','event_counts','measurement','phase_gate','integral'):np.copyto(getattr(self,k),state[k])
  self.logs=[tuple(float(v[0]) if i==0 else v[i].copy() for i in range(7)) for v in zip(*(state['log_'+k] for k in ('time','counts','decoded','current','measurement','armed','integral')))]
  if bool(state['proxy_present']):self.proxy=state['proxy'].copy()
  elif hasattr(self,'proxy'):del self.proxy
 def save(self,folder):
  keys=('t','counts','decoded','input_nA','load_event_measurement','armed','synaptic_integral_nA_s');np.savez_compressed(folder/'load_termination_spikes.npz',**{k:np.array([r[i] for r in self.logs]) for i,k in enumerate(keys)});(folder/'LOAD_TERMINATION_ARCHITECTURE.json').write_text(json.dumps(dict(configuration=self.cfg,edges=self.edge_specs,event_counts=self.event_counts.tolist(),actual_added_spiking_cells=32,spike_counts=self.total_counts.tolist(),scope='New event-limited contralateral normal-load sensory pathway. Four actual spike-driven excitatory axons reach inherited PF-flexor INs, which inhibit opposite PF-extensors. No direct RG/MN/muscle command. Measured normal load is delayed one physics sample; sustained-load and unload-rearming hysteresis reject impacts/chatter. Pulse expires despite sustained contact. Initial loads do not fire until the source foot has unloaded. Source graph108/216 remains unchanged; runtime110/220 and736spiking/80graded cells.'),indent=2))

def probe_load_termination(root,phase_conditioned=False):
 """Isolated event/current/state checks before any coupled optimization."""
 import inspect
 import native_kernel as K
 from spiking_central_runtime import SpikingCentralController as C
 from recruitment import apply_recruitment
 root=Path(root);out=root/'load_termination_v1';out=out/'isolated_phase_v1' if phase_conditioned else out;out.mkdir(exist_ok=True,parents=True);reference=root/'W2L_VirtualWalker_Newer/pace_0p2_supported10_v1';base=json.loads((reference/'SPIKING_COMPARISON.json').read_text());bank=json.loads((reference.parent/'IMPORT_AUDIT.json').read_text())['muscles'];z=np.load(reference/'trajectory.npz');allowed=set(inspect.signature(C.__bases__[0]).parameters)-{'source','bank'};settings=dict(pathway_kind='contralateral_load_termination',termination_gmax_uS=.5,pulse_s=.1,on_load_BW=.25,off_load_BW=.025,load_dwell_s=.04,unload_dwell_s=.05,refractory_s=.25,input_current_nA=4.,phase_conditioned=int(phase_conditioned))
 def make(gain):
  C.configuration=base['configuration'].copy();C.balance_configuration=dict(settings,termination_gmax_uS=gain);C.clock_scale=1.;C.sensory_spiking=True;C.reference_calculations=False;C.audit_state=False;c=C(base['source_configuration']['source'],bank,**{k:v for k,v in base['source_configuration'].items() if k in allowed});cfg=base['source_configuration'];apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale']);c.apply_output_gain();return c
 dt=.002;c=make(.5);pulse=[]
 for t in np.arange(0,2.,dt):c.balance.observe_forces(np.array([.6 if t>=.5 else 0.,0.]),float(t),dt);pulse.append(c.balance.current.copy())
 pulse=np.array(pulse);latch=dict(event_counts=c.balance.event_counts.tolist(),pulse_duration_s=float(np.sum(pulse[:,0]>0)*dt),late_current_nA=float(pulse[-1,0]),scope='One sustained load step after initial unloading: one self-expiring pulse, never repeated during continuous contact.');assert latch['event_counts']==[1,0] and abs(latch['pulse_duration_s']-.1)<=dt and latch['late_current_nA']==0
 rows=[];traces={};state_audit=None;bw=base['result']['robot_mass_kg']*9.81
 for gain in (0.,.5,1.):
  c=make(gain);zeros=np.zeros(len(bank));trace=[]
  for n in range(round(4./dt)):
   t=n*dt;c.balance.observe_forces(z['foot_vertical'][n]/bw,t,dt);c.step(t,zeros,zeros,zeros,{},dt,False)
   if gain==.5 and n==round(2./dt):
    initial=c.export_state();snapshots=[]
    for repeat in range(2):
     c.import_state(initial)
     for j in range(32):
      k=n+j+1;tt=k*dt;c.balance.observe_forces(z['foot_vertical'][k]/bw,tt,dt);c.step(tt,zeros,zeros,zeros,{},dt,False)
     snapshots.append(c.export_state())
    state_audit={key:a.shape==snapshots[1][key].shape and a.tobytes()==snapshots[1][key].tobytes() for key,a in snapshots[0].items()};assert all(state_audit.values());c.import_state(initial)
   trace.append(c.v.copy())
  folder=out/('isolated_g'+str(gain).replace('.','p'));folder.mkdir(exist_ok=True);c.save(folder);trace=np.array(trace);traces[str(gain)]=trace;caps=[]
  for spec in c.balance.edge_specs:
   target=c.names.index(spec['target']);pre=c.names.index(spec['target'].removesuffix(' IN'));edge=next(e for e in c.edges[:216] if int(e[0])==pre and int(e[1])==target);caps.append(float(edge[2]*abs(edge[5]-c.params[target,0])))
  current=np.array([row[6] for row in c.balance.logs])/dt;rows.append(dict(gain=gain,path=str(folder),source_full_PF_to_IN_current_caps_nA=caps,actual_added_current_peak_nA=np.max(abs(current),axis=0).tolist(),current_cap_ratio=float(np.max(np.abs(current)/np.array(caps))),event_counts=c.balance.event_counts.tolist(),new_spikes=c.balance.total_counts.tolist(),PF_E_spikes={name:int(c.total_counts[p]) for p,i in enumerate(c.encoded) if ' PF ext' in (name:=c.names[i])}));np.savez_compressed(folder/'isolated_neural.npz',t=np.arange(len(trace))*dt,v=trace,names=c.names)
 report=dict(source_sha256=base['source_sha256'],configuration=settings,source_graph=[108,216],experimental_runtime_graph=[110,220],experimental_cell_counts=dict(spiking=736,graded=80),sustained_contact_test=latch,rows=rows,complete_state_replay=state_audit,scope='New spike-encoded normal-force event is replayed from actual supported-reference foot loads into an isolated deafferented neural circuit; no physical integration. Actual chemical current into graded PF IN measured every neural substep. Source PFâ†’IN full-drive at rest supplies the declared current envelope. No direct RG or motor input and no claim of grounded stepping. Full event timers/hysteresis, LIF/filter state and histories included in continuation audit.');(out/'ISOLATED_LOAD_TRANSFER.json').write_text(json.dumps(report,indent=2),encoding='utf-8');(out/'CONFIGURATION.json').write_text(json.dumps(settings,indent=2));return report

def plot_termination_timing(root):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 root=Path(root);folder=root/'W2L_VirtualWalker_Newer/loadterm_zero_ground_v1';z=np.load(folder/'trajectory.npz');s=np.load(folder/'load_termination_spikes.npz');r=json.loads((folder/'RESULT.json').read_text());names=z['neuron_names'].tolist();muscles=z['muscle_names'].tolist();roles={x['actual_motor_neuron']:x['muscle'] for x in json.loads((root/'W2L_VirtualWalker_Newer/pace_0p2_supported10_v1/STANCE_LOAD_DIAGNOSIS.json').read_text())['muscles']};bw=r['robot_mass_kg']*9.81;plt.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42});fig,axes=plt.subplots(5,1,figsize=(8,8.6),sharex=True);t=z['t'];colours=('#4477AA','#EE7733')
 for side,col,c,style in [('Left',0,'.15','-'),('Right',1,'#AA3377','--')]:axes[0].plot(t,z['foot_vertical'][:,col]/bw,c=c,ls=style,label=side+' foot',lw=1)
 for joint,c,style in [('Hip',colours[0],'-'),('Knee',colours[1],'--')]:
  name='L '+joint+' PF ext';axes[1].plot(t,np.clip((z['neuro'][:,names.index(name)]+60)/20,0,1),c=c,ls=style,label=joint+' extensor',lw=1);axes[2].plot(t,z['neuro'][:,names.index('L '+joint+' MN ext')],c=c,ls=style,label=joint,lw=1);axes[3].plot(t,z['forces'][:,muscles.index(roles['L '+joint+' MN ext'])],c=c,ls=style,label=joint,lw=1)
 axes[2].axhline(-70,c=colours[0],ls=':',lw=.7,label='Hip activation midpoint');axes[2].axhline(-50,c=colours[1],ls=':',lw=.7,label='Knee activation midpoint')
 for k,label,c,style in [(0,'Hip',colours[0],'-'),(1,'Knee',colours[1],'--'),(2,'Ankle','.35',':')]:axes[4].plot(t,z['angles'][:,k],label=label,c=c,ls=style,lw=1)
 events=t[s['load_event_measurement'][:,3]>.5];events=events[(events>=1.3)&(events<=2.8)]
 for ax in axes:
  for event in events:ax.axvline(event,c='#CC6677',ls='-.',lw=.85)
  ax.grid(alpha=.2);ax.set_xlim(1.3,2.8);ax.legend(loc='upper right',fontsize=10,ncol=2,frameon=False)
 for ax,label in zip(axes,['Foot normal load\n(body weight)','Left PF-E output\n(normalized)','Left extensor MN\n(mV)','Left extensor\ntension (N)','Left joint angle\n(deg)']):ax.set_ylabel(label)
 axes[1].set_ylim(-.03,1.1);axes[4].set_xlabel('Time (s)');fig.suptitle('Qualified opposite-foot load follows extensor-burst termination',fontsize=11);fig.tight_layout(rect=(0,.045,1,.97));fig.text(.06,.015,'Passing reference; new axons disconnected. Red lines: right-foot load-qualified cue. No independent-walking claim.',fontsize=10);out=root/'load_termination_v1';fig.savefig(out/'TIMING_AUDIT.png',dpi=200);fig.savefig(out/'TIMING_AUDIT.pdf');plt.close(fig);(out/'TIMING_AUDIT_ALT.txt').write_text('Five aligned measured panels show foot normal forces, left hip/knee extensor spike-filtered PF output, motor voltage, actual muscle tension, and joint motion. Opposite-foot qualified-load events occur after the left leg has already unloaded and its extensor PF output has decayed. This is the unchanged passing reference with new conductances zero, not evidence of improved balance.',encoding='utf-8')

def write_termination_checkpoint(root):
 root=Path(root);out=root/'load_termination_v1';w=root/'W2L_VirtualWalker_Newer';cases=[]
 for tag in ('loadterm_air_deaff_v1','loadterm_air_aff_v1','loadterm_zero_ground_v1','loadterm_pilot_ground_v1'):
  f=w/tag;c=json.loads((f/'SPIKING_COMPARISON.json').read_text());r=c['result'];z=np.load(f/'trajectory.npz');sp=np.load(f/'load_termination_spikes.npz');F=np.stack([z[k][:,1] for k in ('rail_force_world','guide_force_world','harness_force')],axis=1);bw=r['robot_mass_kg']*9.81;cases.append(dict(path=str(f),air_gate=c.get('air_gate'),ground_gate=c['assessment']['alternating_supported_stepping_gate'],swings=[c['assessment']['legs'][s]['unloaded_clear_swings'] for s in 'LR'],net_assistance_BW=r['mean_net_vertical_assistance_BW'],upward_channel_BW=float(np.maximum(F,0).sum(1).mean()/bw),downward_channel_BW=float(np.maximum(-F,0).sum(1).mean()/bw),progress_relative_belt_m=r['root_forward_relative_to_belt_m'],actual_cell_counts=[r['actual_spiking_cells'],r['actual_graded_cells']],actual_event_counts=json.loads((f/'LOAD_TERMINATION_ARCHITECTURE.json').read_text())['event_counts'],added_spikes=sp['counts'].sum(0).tolist(),added_current_peak_2ms_average_nA=float(np.max(abs(sp['synaptic_integral_nA_s']/r['physics_dt']))),state_audit=c['spiking']['state_snapshot_verification']))
 ledger=out/'rustuna_v1/TRIALS.jsonl';trials=[json.loads(s) for s in ledger.read_text().splitlines()] if ledger.exists() else [];report=dict(machine='easteregg2',source_baseline=str(w/'pace_0p2_supported10_v1'),source_hash=json.loads((w/'pace_0p2_supported10_v1/RESULT.json').read_text())['source_sha256'],source_census=[108,216],experimental_census=[110,220],configuration=str(out/'PILOT_CONFIGURATION.json'),cases=cases,zero_parity=json.loads((out/'ZERO_PARITY.json').read_text()),isolated_probe=str(out/'isolated_phase_v1/ISOLATED_LOAD_TRANSFER.json'),timing=json.loads((out/'BASELINE_EVENT_TIMING.json').read_text()),early_acceptance_timing=json.loads((out/'EARLY_ACCEPTANCE_TIMING.json').read_text()),rustuna=dict(completed_trials=len(trials),planned_max_trials=8,rows=trials,status='Stopped on software constructor-bound error during the third evaluation; two completed trial records valid. Subsequent timing audit found the cue largely postdates actual burst termination, so no resumed gain sweep.',source_snapshot=str(out/'rustuna_v1/ADAPTER_SOURCE.py')),bound_recheck=json.loads((out/'BOUND_RECHECK.json').read_text()),validation_seeds_not_used=[100048,*range(100061,100069)],candidate_promoted=False,scope='Prototype input/state/current behavior verified, not a demonstrated solution to support transfer. No robustness or independent-walking claim. Retain704spiking/80graded baseline. No source model or manuscript changed. The asymmetric existing static LP is a feasibility witness, not a ready symmetric neural teaching pattern.');(out/'CHECKPOINT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report

def fit_current_budget(folder,output):
 """Constrain new synaptic gains using actual source phase/current budgets.

 Direct-MN small-signal response is a fit diagnostic, not a force allocator.
 Reciprocal-Ia network return effects remain for coupled validation.
 """
 from import_virtual_walker import mujoco
 import native_kernel as K
 from recruitment import apply_recruitment
 import inspect
 folder=Path(folder);output=Path(output);output.mkdir(parents=True,exist_ok=True);record=json.loads((folder/'SPIKING_COMPARISON.json').read_text());cfg=record['source_configuration'];z=np.load(folder/'trajectory.npz');bank=json.loads((folder.parent/'IMPORT_AUDIT.json').read_text())['muscles'];allowed=set(inspect.signature(K.NativeController).parameters)-{'source','bank'};c=K.NativeController(cfg['source'],bank,**{k:v for k,v in cfg.items() if k in allowed});apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale'])
 for row in record['spiking']['output_gain_edits']:
  matches=[e for e in c.edges if c.names[int(e[0])]==row['source'] and c.names[int(e[1])]==row['target']];assert len(matches)==1;matches[0][2]=row['effective_gmax_uS']
 model=mujoco.MjModel.from_xml_path(str(folder/'model.xml'));data=mujoco.MjData(model);params=np.load(folder/'POOLED_FORCE_PARAMETERS.npz')['pooled'];outputs={c.names[pre]:(pre,index,coeff) for pre,index,coeff in c.outputs};bw=record['result']['robot_mass_kg']*9.81;rows=[]
 for gmax in (.15,.3,.6,1.):
  responses=[]
  for side,col in [('L',0),('R',1)]:
   mask=(z['t']>=1)&(z['t']<=8)&(z['foot_vertical'][:,col]>.1*bw)&(z['foot_vertical'][:,1-col]<.025*bw);candidates=np.flatnonzero(mask);samples=candidates[np.linspace(0,len(candidates)-1,min(12,len(candidates))).astype(int)]
   for n in samples:
    v=z['neuro'][n,:108];pre,idx,coeff=outputs[side+' Hip MN flx'];ia=c.names.index(side+' Hip flx Ia');active=np.clip((v[c.edges[:,0].astype(int)]-c.edges[:,3])/(c.edges[:,4]-c.edges[:,3]),0,1);incoming=c.edges[:,1].astype(int)==pre;conductance=float(np.sum(c.edges[incoming,2]*active[incoming]));extra=gmax*.1;current=extra*(-10.-v[pre]);delta=current/(1.+conductance+extra);vm=K.polynomial(coeff,v[pre]*.001);newvm=K.polynomial(coeff,(v[pre]+delta)*.001);data.qpos[:]=z['q'][n];mujoco.mj_kinematics(model,data);mujoco.mj_comPos(model,data);mujoco.mj_tendon(model,data);tid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_TENDON,'route_'+bank[idx]['name']);j=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'hip_'+side);moment=-data.ten_J.reshape(model.ntendon,model.nv)[tid,model.jnt_dofadr[j]];assert moment>0;l0,width,kse,kpe,damping,maximum,A,B,C,D,*_=params[idx];fl=max(0,1-((data.ten_length[tid]-l0)/width)**2);force_delta=(max(0,B/(1+np.exp(C*(A-newvm)))+D)-max(0,B/(1+np.exp(C*(A-vm)))+D))*fl/(1+kpe/kse);responses.append(dict(side=side,time_s=float(z['t'][n]),source_MN_mV=float(v[pre]),source_incoming_g_uS=conductance,added_current_nA=float(current),predicted_direct_MN_delta_mV=float(delta),predicted_direct_hip_torque_delta_Nm=float(moment*force_delta),reciprocal_Ia_current_at_2nA_input_nA=float(gmax*.5*.1*(-10-v[ia]))))
  torques=[row['predicted_direct_hip_torque_delta_Nm'] for row in responses];rows.append(dict(mn_gmax_uS=gmax,ia_gmax_uS=gmax*.5,median_direct_torque_delta_Nm=float(np.median(torques)),maximum_MN_current_at_4nA_input_nA=float(max(row['added_current_nA'] for row in responses)*2),responses=responses))
 eligible=[row for row in rows if 2<=row['median_direct_torque_delta_Nm']<=8 and row['maximum_MN_current_at_4nA_input_nA']<=20];chosen=min(eligible,key=lambda row:abs(row['median_direct_torque_delta_Nm']-4)) if eligible else rows[-1];assert chosen['maximum_MN_current_at_4nA_input_nA']<=20
 configuration=dict(kp_nA_per_rad=40.,kd_s=.15,target_pitch_deg=float(-np.degrees(2*np.arctan2(record['result']['initial_qpos'][6],record['result']['initial_qpos'][3]))),delay_s=.04,max_current_nA=4.,mn_gmax_uS=chosen['mn_gmax_uS'],ia_gmax_uS=chosen['ia_gmax_uS'],swing_gain=.25);report=dict(source=str(folder),source_sha256=record['source_sha256'],configuration=configuration,candidates=rows,scope='Actual saved phase-specific source conductance and motor voltages. Moderate2nA relay input corresponds to0.1 source activation;4nA cap gives0.2. Direct MN current/voltage change estimated with existing source conductance, then source Hill force-length law and actual moment arm. This ignores dynamic/reciprocal network return, so coupled validation is mandatory. No force targets or joint torques applied. Gain fit aims for2â€“8Nm small-signal stance correction and <=20nA per MN, not a physiological identification.')
 (output/'BODY_BALANCE_CURRENT_FIT.json').write_text(json.dumps(report,indent=2));(output/'BODY_BALANCE_CANDIDATE.json').write_text(json.dumps(configuration,indent=2));return report

def fit_functional_current_budget(folder,output):
 """Calibrate new ankle/clearance projection to existing source PF capacity."""
 from import_virtual_walker import mujoco
 import native_kernel as K
 import inspect
 folder=Path(folder);output=Path(output);output.mkdir(parents=True,exist_ok=True);record=json.loads((folder/'SPIKING_COMPARISON.json').read_text());z=np.load(folder/'trajectory.npz');cfg=record['source_configuration'];bank=json.loads((folder.parent/'IMPORT_AUDIT.json').read_text())['muscles'];allowed=set(inspect.signature(K.NativeController).parameters)-{'source','bank'};c=K.NativeController(cfg['source'],bank,**{k:v for k,v in cfg.items() if k in allowed});PF={row['target']:row['effective_gmax_uS'] for row in record['spiking']['output_gain_edits']}
 from recruitment import apply_recruitment
 apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale'])
 for row in record['spiking']['output_gain_edits']:
  for e in c.edges:
   if c.names[int(e[0])]==row['source'] and c.names[int(e[1])]==row['target']:e[2]=row['effective_gmax_uS']
 model=mujoco.MjModel.from_xml_path(str(folder/'model.xml'));data=mujoco.MjData(model);parameters=np.load(folder/'POOLED_FORCE_PARAMETERS.npz')['pooled'];outputs={c.names[pre]:(pre,index,coeff) for pre,index,coeff in c.outputs};bw=record['result']['robot_mass_kg']*9.81;rows=[];limits={}
 for joint,action,phase,gmax,expected_sign in [('Ank','ext','stance',4.,1),('Ank','flx','stance',4.,-1),('Knee','flx','swing',1.,1),('Knee','ext','swing',1.,-1)]:
  for side,col in [('L',0),('R',1)]:
   name=side+' '+joint+' MN '+action;pre,idx,coeff=outputs[name];mask=(z['t']>=1)&(z['t']<=8)
   if phase=='stance':mask&=(z['foot_vertical'][:,col]>.1*bw)&(z['foot_vertical'][:,1-col]<.025*bw)
   else:
    F=c.names.index(side+' RG flx');E=c.names.index(side+' RG ext');mask&=(z['foot_vertical'][:,col]<.025*bw)&(z['foot_vertical'][:,1-col]>.1*bw)&(z['neuro'][:,F]>z['neuro'][:,E])
   indices=np.flatnonzero(mask);indices=indices[np.linspace(0,len(indices)-1,min(16,len(indices))).astype(int)];responses=[]
   for n in indices:
    v=z['neuro'][n,:108];active=np.clip((v[c.edges[:,0].astype(int)]-c.edges[:,3])/(c.edges[:,4]-c.edges[:,3]),0,1);incoming=c.edges[:,1].astype(int)==pre;gs=np.sum(c.edges[incoming,2]*active[incoming]);extra=gmax*.1;current=extra*(-10-v[pre]);delta=current/(1+gs+extra);vm=K.polynomial(coeff,v[pre]*.001);newvm=K.polynomial(coeff,(v[pre]+delta)*.001);data.qpos[:]=z['q'][n];mujoco.mj_kinematics(model,data);mujoco.mj_comPos(model,data);mujoco.mj_tendon(model,data);tid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_TENDON,'route_'+bank[idx]['name']);j=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,('ankle_' if joint=='Ank' else 'knee_')+side);moment=-data.ten_J.reshape(model.ntendon,model.nv)[tid,model.jnt_dofadr[j]];assert moment*expected_sign>0;l0,width,kse,kpe,damping,maximum,A,B,C,D,*_=parameters[idx];fl=max(0,1-((data.ten_length[tid]-l0)/width)**2);change=(max(0,B/(1+np.exp(C*(A-newvm)))+D)-max(0,B/(1+np.exp(C*(A-vm)))+D))*fl/(1+kpe/kse);responses.append(dict(time_s=float(z['t'][n]),MN_mV=float(v[pre]),added_current_at_2nA_input_nA=float(current),direct_predicted_torque_change_Nm=float(moment*change)))
   full_source_current=PF[name]*(-10-c.params[pre,0]);limits.setdefault(joint,[]).append(float(1.25*full_source_current));rows.append(dict(target=name,phase=phase,anatomical_moment_sign=expected_sign,new_maximum_gmax_uS=gmax,effective_g_at_4nA_mean_input_uS=gmax*.2,source_PF_full_gmax_uS=PF[name],source_PF_full_current_at_rest_nA=float(full_source_current),median_direct_torque_change_at_2nA_input_Nm=float(np.median([x['direct_predicted_torque_change_Nm'] for x in responses])),samples=responses))
 candidate=json.loads((output.parent/'body_balance_v1/BODY_BALANCE_CANDIDATE.json').read_text());candidate.update(mn_gmax_uS=1.,ia_gmax_uS=.5,swing_gain=1.,capture_weight=1.,functional_projection=1,ankle_gain_ratio=4.,capture_margin_m=.005);report=dict(source=str(folder),source_sha256=record['source_sha256'],rows=rows,current_limits_nA={joint:min(values) for joint,values in limits.items()},candidate=candidate,scope='Same actual LIF encoder and source20mV activation window. At4nA mean input, ankle configured4uS supplies0.8uS effective conductance, within source0.812uS full PF input; knee configured1uS supplies0.2uS, below both source knee PF gains. Current envelope uses125% of the smaller antagonist full-PF source current at rest, allowing measured spike-filter ripple; this is a declared new joint-specific calibration, not the prior hip20nA envelope. Local direct-MN estimates exclude reciprocal network return and require coupled validation. Mechanical signs checked from actual tendon Jacobians. No torque/force targets applied online.')
 (output/'FUNCTIONAL_CURRENT_FIT.json').write_text(json.dumps(report,indent=2));(output/'CAPTURE_BALANCE_CANDIDATE.json').write_text(json.dumps(candidate,indent=2));return report

def fit_task_actions(folder,output):
 """Offline muscle-force sensitivities, never an online force allocator."""
 from import_virtual_walker import mujoco
 import native_kernel as K
 folder=Path(folder);output=Path(output);output.mkdir(parents=True,exist_ok=True);z=np.load(folder/'trajectory.npz');r=json.loads((folder/'RESULT.json').read_text());model=mujoco.MjModel.from_xml_path(str(folder/'model.xml'));bank=json.loads((folder.parent/'IMPORT_AUDIT.json').read_text())['muscles'];record=json.loads((folder/'SPIKING_COMPARISON.json').read_text());c=K.NativeController(record['source_configuration']['source'],bank);root=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Root');act=[(i,mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_ACTUATOR,b['name'])) for i,b in enumerate(bank) if b['type']=='LinearHillMuscle'];motor=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_ACTUATOR,'motor_Treadmill_Slide');dt=r['physics_dt'];rows=[]
 for g in range(model.ngeom):
  if 'support_rail' in (mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_GEOM,g) or ''):model.geom_contype[g]=0;model.geom_conaffinity[g]=0
 for time in (1.55,1.85,2.05,2.2):
  n=int(np.argmin(abs(z['t']-time)));q=z['q'][n-1];v=np.empty(model.nv);mujoco.mj_differentiatePos(model,v,dt,z['q'][n-2],q);data=mujoco.MjData(model);data.qpos[:]=q;data.qvel[:]=v
  for i,a in act:data.ctrl[a]=z['forces'][n,i]
  data.ctrl[motor]=r.get('belt_command_m_s',.1);mujoco.mj_forward(model,data);baseline=data.qacc.copy();warm=data.qacc_warmstart.copy();J=np.empty((3,model.nv));rot=np.empty_like(J);mujoco.mj_jacSubtreeCom(model,data,J,root);mujoco.mj_jacBody(model,data,np.empty_like(J),rot,root);feet=[]
  for side in 'LR':
   jac=np.empty_like(J);mujoco.mj_jacBody(model,data,jac,np.empty_like(J),mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'foot_'+side));feet.append(jac[0].copy()-J[0])
  for i,a in act:
   old=float(data.ctrl[a]);data.ctrl[a]=min(old+5.,model.actuator_ctrlrange[a,1]);increment=float(data.ctrl[a]-old)
   if increment<=0:continue
   data.qacc_warmstart[:]=warm;mujoco.mj_forward(model,data);delta=(data.qacc-baseline)/increment;name=next(c.names[pre] for pre,target,coeff in c.outputs if target==i);rows.append(dict(time_s=float(z['t'][n]),motor_neuron=name,physical_route=bank[i]['name'],source_tension_N=old,force_increment_N=increment,COM_ax_per_N=float(J[0]@delta),COM_ay_per_N=float(J[1]@delta),pelvis_anterior_pitch_accel_per_N=float(-rot[2]@delta),left_foot_relative_COM_ax_per_N=float(feet[0]@delta),right_foot_relative_COM_ax_per_N=float(feet[1]@delta),saved_foot_load_N=z['foot_vertical'][n].tolist()));data.ctrl[a]=old
 report=dict(path=str(folder),source_sha256=r['source_sha256'],rows=rows,scope='Offline local forward-dynamics sensitivity to +5N in existing muscle tension, at actual saved phase poses/velocities with horizontal rails removed and guides retained. No integration, training force target, inverse allocator, or direct joint torque. Constraints re-solved after each perturbation. Instantaneous action may vary with contact mode; this does not model neural/actuator delay or prove recovery.');(output/'TASK_ACTION_SENSITIVITY.json').write_text(json.dumps(report,indent=2));return report

def run_balance_case(stage,configuration,manifest,tag):
 """Persistent-worker coupled evaluation, with frozen baseline neural/physics."""
 import hashlib
 import native_kernel as K
 from spiking_central_runtime import SpikingCentralController as C
 from run_virtual_walker import run
 from assess_trial import assess
 baseline=json.loads(Path(manifest['baseline_record']).read_text());C.configuration=baseline['configuration'].copy();C.sensory_spiking=True;C.reference_calculations=False;C.audit_state=False;C.clock_scale=1.;C.balance_configuration=configuration;cfg=baseline['source_configuration'].copy();cfg.update(newer=True,mode='supported' if stage=='ground' else 'drop',duration=10. if stage=='ground' else 6.,feedback=True,release_at=2.,tag=tag,output_root=Path(manifest['output_dir']));previous=K.NativeController;K.NativeController=C
 try:r=run(**cfg)
 finally:K.NativeController=previous
 folder=Path(manifest['output_dir'])/'W2L_VirtualWalker_Newer'/tag;c=C.last_instance;c.save(folder);assessment=assess(folder,start_s=1. if stage=='ground' else 2.5);meta=c.metadata();r.update(controller_kind=meta['architecture'],logical_source_channels=108,runtime_logical_channels=len(c.v),actual_spiking_cells=meta['actual_spiking_cells'],actual_graded_cells=meta['remaining_graded_cells'],scope=meta['scope']);(folder/'RESULT.json').write_text(json.dumps(r,indent=2));(folder/'BODY_BALANCE_CONFIGURATION.json').write_text(json.dumps(configuration,indent=2));(folder/'SPIKING_COMPARISON.json').write_text(json.dumps(dict(source_sha256=baseline['source_sha256'],configuration=C.configuration,source_configuration=cfg,spiking=meta,result=r,assessment=assessment,mode=cfg['mode'],rail_release_at_s=2. if stage!='ground' else None,balance_configuration=configuration,implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,default=str));return folder,r,assessment

def added_current_peak(folder,dt):
 z=np.load(folder/'body_balance_spikes.npz');edges=json.loads((folder/'BODY_BALANCE_ARCHITECTURE.json').read_text())['edges'];targets=sorted(set(e['target'] for e in edges));currents=np.column_stack([z['synaptic_integral_nA_s'][:,[i for i,e in enumerate(edges) if e['target']==target]].sum(1)/dt for target in targets]);return float(np.max(np.abs(currents)))

def current_budget_ratio(folder,dt,limits):
 if limits is None:return added_current_peak(folder,dt)/20.
 z=np.load(folder/'body_balance_spikes.npz');edges=json.loads((folder/'BODY_BALANCE_ARCHITECTURE.json').read_text())['edges'];ratios=[]
 for target in set(e['target'] for e in edges):
  group='Ank' if ' Ank ' in target else 'Knee' if ' Knee ' in target else 'Hip';current=z['synaptic_integral_nA_s'][:,[i for i,e in enumerate(edges) if e['target']==target]].sum(1)/dt;ratios.append(float(np.max(abs(current)))/limits[group])
 return max(ratios)

def handle_request(request):
 """Rustuna adapter: search only new balance-pathway gains, never source timing."""
 import hashlib
 manifest=request['manifest'];op=request['op'];baseline=json.loads(Path(manifest['baseline_record']).read_text());expected=baseline['source_sha256'];assert hashlib.sha256(Path(baseline['source_configuration']['source']).read_bytes()).hexdigest()==expected
 if manifest.get('objective_kind')=='motor_output':return handle_motor_request(request,baseline)
 if manifest.get('objective_kind')=='load_termination':return handle_termination_request(request,baseline)
 if op=='preflight':
  fit=json.loads(Path(manifest['current_fit']).read_text());parity=json.loads(Path(manifest['parity_record']).read_text());air=[json.loads(Path(p).read_text())['air_gate']['passed'] for p in manifest['air_records']];ready=baseline['assessment']['alternating_supported_stepping_gate'] and fit['source_sha256']==expected and all(air) and all(all(v.values()) for v in parity['checks'].values()) and parity['state_audit']['all_arrays_byte_identical'];return dict(ok=True,ready_for_optimization=bool(ready),reason='Frozen source baseline, both air gates, exact default/zero-path controls, source-current fit and state continuation verified. New balance release still fails and is the tuning target.',source_sha256=expected,scope='Only three new symmetric task-balance parameters; source neural timing/recruitment, valid home, pooled physical law and guides fixed.')
 assert op=='evaluate';sample=request['params'];assert set(sample)=={'mn_gmax_uS','kp_nA_per_rad','swing_gain'};configuration=json.loads(Path(manifest['balance_configuration']).read_text());configuration.update(sample);configuration['ia_gmax_uS']=.5*configuration['mn_gmax_uS'];identifier=request['evaluation_id'];ground,r,a=run_balance_case('ground',configuration,manifest,'balance_'+identifier+'_ground');peak=added_current_peak(ground,r['physics_dt']);diagnostics=dict(ground_path=str(ground),configuration=configuration,ground_gate=a['alternating_supported_stepping_gate'],ground_swings=[a['legs'][s]['unloaded_clear_swings'] for s in 'LR'],ground_support_BW=r['mean_net_vertical_assistance_BW'],ground_belt_progress_m=r['root_forward_relative_to_belt_m'],actual_added_current_peak_nA=peak,source_sha256=expected)
 ratio=current_budget_ratio(ground,r['physics_dt'],manifest.get('current_limits_nA'));diagnostics['new_current_budget_ratio']=ratio
 if not a['alternating_supported_stepping_gate'] or r['root_forward_relative_to_belt_m']<=0 or ratio>1:
  value=20000+200*(10-r['duration_completed_s'])+1000*max(0,ratio-1)-sum(diagnostics['ground_swings']);diagnostics['reason']='Supported gate, forward progression, or declared joint-specific new-current envelope failed; no release promotion.';return dict(ok=True,value=float(value),diagnostics=diagnostics)
 release,r,a=run_balance_case('release',configuration,manifest,'balance_'+identifier+'_release');z=np.load(release/'trajectory.npz');mask=z['t']>=2.;bw=r['robot_mass_kg']*9.81;channels=np.stack([z[key][:,1] for key in ('rail_force_world','guide_force_world','harness_force')],axis=1);up=float(np.maximum(channels[mask],0).sum(1).mean()/bw);down=float(np.maximum(-channels[mask],0).sum(1).mean()/bw);peak=max(peak,added_current_peak(release,r['physics_dt']));duration=r['duration_completed_s'];passed=a['alternating_supported_stepping_gate'] and r['failure'] is None and duration>=5.999 and up<=.05 and peak<=20
 ratio=max(ratio,current_budget_ratio(release,r['physics_dt'],manifest.get('current_limits_nA')));passed=a['alternating_supported_stepping_gate'] and r['failure'] is None and duration>=5.999 and up<=.05 and ratio<=1;diagnostics.update(release_path=str(release),release_duration_s=duration,release_failure=r['failure'],release_gate=a['alternating_supported_stepping_gate'],release_swings=[a['legs'][s]['unloaded_clear_swings'] for s in 'LR'],release_upward_channel_BW=up,release_downward_channel_BW=down,actual_added_current_peak_nA=peak,new_current_budget_ratio=ratio,physical_pass=bool(passed));value=(100*up+5*down) if passed else 1000+200*(6-duration)-sum(diagnostics['release_swings'])+2000*max(0,ratio-1)+1000*max(0,up-.05);return dict(ok=True,value=float(value),diagnostics=diagnostics)

def run_motor_case(mode,overrides,manifest,tag,feedback=True,pathway_configuration=None,perturbation_seed=None):
 """Source motor-output cases with an optional explicitly declared neural relay."""
 from scipy.signal import find_peaks
 import hashlib,native_kernel as K
 from spiking_central_runtime import SpikingCentralController as C
 from run_virtual_walker import run
 from assess_trial import assess
 base=json.loads(Path(manifest['baseline_record']).read_text());C.configuration=dict(base['configuration'],**overrides);C.sensory_spiking=True;C.balance_configuration=None if pathway_configuration is None else pathway_configuration.copy();C.clock_scale=1.;C.reference_calculations=False;C.audit_state=False;cfg=base['source_configuration'].copy();cfg.update(newer=True,mode=mode,duration=6. if mode=='air' else 10.,feedback=feedback,tag=tag,output_root=Path(manifest['output_dir']));factors={}
 if perturbation_seed is not None:
  from hip_knee_interlock import OPTIONAL_KEYS
  rng=np.random.default_rng(perturbation_seed)
  for key in C.configuration:
   if key=='rg_reset_mV' or key in OPTIONAL_KEYS:continue
   factors[key]=float(rng.uniform(.99,1.01));C.configuration[key]*=factors[key]
  for key in ('pf_mn_recruitment_scale','toe_mn_gmax','hip_pf_mn_scale','hip_flexor_pf_mn_scale'):factors[key]=float(rng.uniform(.99,1.01));cfg[key]*=factors[key]
  if 'rg_reset_mV' in C.configuration:factors['rg_reset_mV']=float(rng.uniform(.99,1.01));C.configuration['rg_reset_mV']*=factors['rg_reset_mV']
  for key in OPTIONAL_KEYS:
   if key in C.configuration:factors[key]=float(rng.uniform(.99,1.01));C.configuration[key]*=factors[key]
  if C.balance_configuration is not None:
   from stance_recruitment import StanceRecruitmentRelay
   learned=StanceRecruitmentRelay.learned_keys(C.balance_configuration) if C.balance_configuration.get('pathway_kind')=='phase_stance_recruitment' else LoadTransferTerminationRelay.learned
   for key in learned:factors['pathway.'+key]=float(rng.uniform(.99,1.01));C.balance_configuration[key]*=factors['pathway.'+key]
 if mode=='air':cfg['fixed_root_pose']=cfg['initial_qpos'][:7];cfg['initial_qpos']=cfg['initial_qpos'][7:]
 previous=K.NativeController;K.NativeController=C
 try:r=run(**cfg)
 finally:K.NativeController=previous
 folder=Path(manifest['output_dir'])/'W2L_VirtualWalker_Newer'/tag;c=C.last_instance;c.save(folder);assessment=assess(folder,start_s=1.);z=np.load(folder/'trajectory.npz');air=None
 if mode=='air':
  mask=z['t']>=1.;times={};legs={}
  for side,col in [('L',1),('R',4)]:
   peaks,_=find_peaks(z['angles'][mask,col],prominence=5,distance=round(.3/r['physics_dt']));times[side]=z['t'][mask][peaks];legs[side]=dict(knee_cycles=len(peaks),period_s=float(np.median(np.diff(times[side]))) if len(peaks)>1 else None,knee_range_deg=float(np.ptp(z['angles'][mask,col])))
  errors=[]
  for left,right in zip(times['L'][:-1],times['L'][1:]):
   other=times['R'][(times['R']>=left)&(times['R']<right)]
   if len(other):errors.append(abs(float((other[0]-left)/(right-left))-.5))
  phase=float(np.mean(errors)) if errors else .5;air=dict(passed=bool(r['failure'] is None and all(x['knee_cycles']>=3 and x['period_s'] is not None and .7<=x['period_s']<=1.4 for x in legs.values()) and phase<=.15),legs=legs,knee_peak_antiphase_error=phase)
 meta=c.metadata();r.update(controller_kind=meta['architecture'],actual_spiking_cells=meta['actual_spiking_cells'],actual_graded_cells=meta['remaining_graded_cells']);record=dict(source_sha256=base['source_sha256'],configuration=C.configuration,source_configuration=cfg,spiking=meta,result=r,assessment=assessment,air_gate=air,mode=mode,feedback=feedback,perturbation_seed=perturbation_seed,perturbation_factors=factors,balance_configuration=C.balance_configuration,implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),additional_implementation_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ('hip_knee_interlock.py','spiking_central_runtime.py','run_virtual_walker.py')});(folder/'RESULT.json').write_text(json.dumps(r,indent=2));(folder/'SPIKING_COMPARISON.json').write_text(json.dumps(record,indent=2,default=str));return folder,r,assessment,air,c

def handle_motor_request(request,baseline):
 """One shared knee-flexor output scalar; freeze all source timing and mechanics."""
 manifest=request['manifest']
 if request['op']=='preflight':
  parity=json.loads(Path(manifest['output_checkpoint']).read_text())['default_path_parity'];fit=json.loads(Path(manifest['current_fit']).read_text());ready=baseline['assessment']['alternating_supported_stepping_gate'] and bool(parity) and all(all(x.values()) for x in parity.values()) and fit['source_sha256']==baseline['source_sha256'];return dict(ok=True,ready_for_optimization=bool(ready),reason='Passing physical reference, source-current fit, positive knee guard, and exact default-path parity recorded. One shared output scalar; air gates are repeated before each ground evaluation.')
 assert request['op']=='evaluate';params=request['params'];assert set(params)=={'knee_flexor_pf_output_multiplier'};assert 1<=params['knee_flexor_pf_output_multiplier']<=2.6;tag='motor_'+request['evaluation_id'];diagnostics=dict(parameters=params,air={},source_sha256=baseline['source_sha256'])
 for kind,feedback in [('deaff',False),('aff',True)]:
  folder,r,a,air,c=run_motor_case('air',params,manifest,tag+'_air_'+kind,feedback);diagnostics['air'][kind]=dict(path=str(folder),**air)
  if not air['passed']:diagnostics['physical_pass']=False;return dict(ok=True,value=float(20000-sum(x['knee_cycles'] for x in air['legs'].values())),diagnostics=diagnostics)
 folder,r,a,_,c=run_motor_case('supported',params,manifest,tag+'_ground');z=np.load(folder/'trajectory.npz');bw=r['robot_mass_kg']*9.81;F=np.stack([z[key][:,1] for key in ('rail_force_world','guide_force_world','harness_force')],axis=1);up=float(np.maximum(F,0).sum(1).mean()/bw);down=float(np.maximum(-F,0).sum(1).mean()/bw);foot=z['foot_vertical'].sum(1)/bw;p95=float(np.percentile(foot,95));peak=float(foot.max());ratios=[]
 for edit in c.output_edits:
  if edit['target'] not in ('L Knee MN flx','R Knee MN flx'):continue
  pre=c.names.index(edit['source']);post=c.names.index(edit['target']);e=next(e for e in c.edges if int(e[0])==pre and int(e[1])==post);current=e[2]*np.clip((z['neuro'][:,pre]-e[3])/(e[4]-e[3]),0,1)*(e[5]-z['neuro'][:,post]);cap=edit['before_gmax_uS']*abs(e[5]-c.params[post,0]);ratios.append(float(np.max(abs(current))/cap))
 ratio=max(ratios);passed=bool(a['alternating_supported_stepping_gate'] and r['failure'] is None and r['root_forward_relative_to_belt_m']>0 and ratio<=1);diagnostics.update(ground_path=str(folder),ground_gate=a['alternating_supported_stepping_gate'],physical_pass=passed,swings=[a['legs'][s]['unloaded_clear_swings'] for s in 'LR'],upward_channel_BW=up,downward_channel_BW=down,net_assistance_BW=r['mean_net_vertical_assistance_BW'],belt_relative_progress_m=r['root_forward_relative_to_belt_m'],knee_min_deg=r['anatomical_knee_min_deg'],foot_p95_BW=p95,foot_peak_BW=peak,PF_current_to_unattenuated_source_cap_ratio=ratio)
 value=100*up+.05*p95+.005*peak if passed else 10000-sum(diagnostics['swings'])+1000*max(0,ratio-1);return dict(ok=True,value=float(value),diagnostics=diagnostics)

def handle_termination_request(request,baseline):
 """Three shared load-event parameters, scored across prespecified training seeds."""
 manifest=request['manifest'];iso=json.loads(Path(manifest['isolated_probe']).read_text());parity=json.loads(Path(manifest['zero_parity']).read_text())
 if request['op']=='preflight':
  air=[json.loads(Path(f).read_text())['air_gate']['passed'] for f in manifest['air_records']];ready=baseline['assessment']['alternating_supported_stepping_gate'] and all(air) and parity['all_original_numerical_arrays_exact'] and parity['original_names_equal'] and all(all(d.values()) for d in parity['spike_checks'].values()) and all(iso['complete_state_replay'].values());return dict(ok=True,ready_for_optimization=bool(ready),reason='Passing baseline, both real air gates with zero new spikes, source-current envelope, complete event-state replay, and zero-path numerical/spike parity verified. Trials reuse force-free air invariance; new path cannot activate without load. Nominal plus two declared training seeds are evaluated;100048 and100061â€“68 not used for fitting.')
 assert request['op']=='evaluate';params=request['params'];assert set(params)==set(LoadTransferTerminationRelay.learned);configuration=json.loads(Path(manifest['pathway_configuration']).read_text());configuration.update(params);caps=np.array(iso['rows'][0]['source_full_PF_to_IN_current_caps_nA']);rows=[];identifier=request['evaluation_id']
 for seed in (None,*manifest['training_seeds']):
  tag='loadterm_'+identifier+('_nominal' if seed is None else '_s'+str(seed));folder,r,a,_,c=run_motor_case('supported',{},manifest,tag,True,configuration,seed);z=np.load(folder/'trajectory.npz');s=np.load(folder/'load_termination_spikes.npz');bw=r['robot_mass_kg']*9.81;F=np.stack([z[key][:,1] for key in ('rail_force_world','guide_force_world','harness_force')],axis=1);up=float(np.maximum(F,0).sum(1).mean()/bw);down=float(np.maximum(-F,0).sum(1).mean()/bw);foot=z['foot_vertical'].sum(1)/bw;ratio=float(np.max(np.abs(s['synaptic_integral_nA_s']/r['physics_dt'])/caps));passed=bool(a['alternating_supported_stepping_gate'] and r['failure'] is None and r['root_forward_relative_to_belt_m']>0 and ratio<=1);row=dict(path=str(folder),seed=seed,gate=a['alternating_supported_stepping_gate'],physical_pass=passed,swings=[a['legs'][side]['unloaded_clear_swings'] for side in 'LR'],upward_channel_BW=up,downward_channel_BW=down,net_assistance_BW=r['mean_net_vertical_assistance_BW'],belt_relative_progress_m=r['root_forward_relative_to_belt_m'],knee_min_deg=r['anatomical_knee_min_deg'],event_counts=c.balance.event_counts.tolist(),added_spikes=c.balance.total_counts.tolist(),added_current_peak_nA=float(np.max(np.abs(s['synaptic_integral_nA_s']/r['physics_dt']))),current_cap_ratio=ratio,foot_p95_BW=float(np.percentile(foot,95)),foot_peak_BW=float(foot.max()),failure=r['failure']);rows.append(row)
  if not passed:return dict(ok=True,value=float(10000-sum(row['swings'])+1000*max(0,ratio-1)),diagnostics=dict(configuration=configuration,cases=rows,physical_pass=False,reason='Air-compatible architecture failed an actual grounded/current/forward/kinematic gate; no support promotion.'))
 ups=np.array([r['upward_channel_BW'] for r in rows]);value=100*(ups.mean()+.5*ups.max())+.05*np.mean([r['foot_p95_BW'] for r in rows])+.005*max(r['foot_peak_BW'] for r in rows);return dict(ok=True,value=float(value),diagnostics=dict(configuration=configuration,cases=rows,physical_pass=True,mean_upward_BW=float(ups.mean()),worst_upward_BW=float(ups.max()),scope='All three nominal/training cases must pass first. Mean plus half-worst upward support, with small physical-impact penalties; not a metabolic objective. No validation seed used.'))

def landing_geometry(folder,output):
 """Actual force-qualified landing geometry and pace, without simulation."""
 from import_virtual_walker import mujoco
 folder=Path(folder);z=np.load(folder/'trajectory.npz');r=json.loads((folder/'RESULT.json').read_text());model=mujoco.MjModel.from_xml_path(str(folder/'model.xml'));data=mujoco.MjData(model);bw=r['robot_mass_kg']*9.81;dt=r['physics_dt'];belt=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Treadmill_Belt');joint=int(model.body_jntadr[belt]);adr=int(model.jnt_qposadr[joint]);data.qpos[:]=z['q'][0];mujoco.mj_kinematics(model,data);axis=(data.xmat[belt].reshape(3,3)@model.jnt_axis[joint])[0];beltx=axis*(z['q'][:,adr]-z['q'][0,adr]);beltv=np.gradient(beltx,z['t']);events=[];names=z['neuron_names'].tolist();periods={}
 for side in 'LR':
  F=names.index(side+' RG flx');E=names.index(side+' RG ext');on=z['neuro'][:,F]>z['neuro'][:,E];starts=np.flatnonzero(on[1:]&~on[:-1])+1;times=z['t'][starts];times=times[times>=1];periods[side]=float(np.median(np.diff(times))) if len(times)>1 else None
 for side,col in [('L',0),('R',1)]:
  on=False;high=low=0;landing=None;segments=[]
  for n,force in enumerate(z['foot_vertical'][:,col]/bw):
   if not on:
    high=high+1 if force>=.05 else 0
    if high>=round(.04/dt):landing=n-high+1;on=True;low=0
   else:
    low=low+1 if force<.025 else 0
    if low>=round(.05/dt):segments.append((landing,n-low+1));on=False;high=0
  if on:segments.append((landing,len(z['t'])-1))
  for n,end in segments:
   if z['t'][n]<1 or end-n<round(.04/dt):continue
   data.qpos[:]=z['q'][n];mujoco.mj_kinematics(model,data);bounds=[]
   for region,name in enumerate(('foot_'+side,'toe_'+side)):
    k=2*col+region
    if z['foot_region_force_world'][n:min(end,n+round(.04/dt)),k,1].mean()<.01*bw:continue
    g=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+name);extent=np.abs(data.geom_xmat[g].reshape(3,3)[0])@model.geom_size[g];bounds.extend([data.geom_xpos[g,0]-extent,data.geom_xpos[g,0]+extent])
   if not bounds:continue
   floor=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,'geom_Treadmill_Belt');top=data.geom_xpos[floor,1]+model.geom_size[floor,1];height=max(.1,z['com_position'][n,1]-top);omega=np.sqrt(9.81/height);relative_velocity=z['com_velocity'][n,0]-beltv[n];com=z['com_position'][n,0];capture=com+relative_velocity/omega;gap=max(0,min(bounds)-com);events.append(dict(side=side,time_s=float(z['t'][n]),stance_until_s=float(z['t'][end]),stance_duration_s=float(z['t'][end]-z['t'][n]),foot_frame_ahead_COM_m=float(z['foot_body_position'][n,col,0]-com),loaded_rear_ahead_COM_m=float(min(bounds)-com),capture_margin_m=float(min(capture-min(bounds),max(bounds)-capture)),COM_relative_belt_velocity_m_s=float(relative_velocity),relative_velocity_to_put_XCOM_at_loaded_rear_m_s=float(omega*gap),constant_velocity_catch_time_s=float(gap/relative_velocity) if relative_velocity>0 else None,forward_foot_frame_separation_m=float(z['foot_body_position'][n,col,0]-z['foot_body_position'][n,1-col,0]),belt_x_m=float(beltx[n])))
 events.sort(key=lambda x:x['time_s']);same={};steps=[]
 for row in events:
  side=row['side']
  if side in same:
   old=same[side];row['same_leg_stride_relative_belt_m']=float((z['foot_body_position'][np.argmin(abs(z['t']-row['time_s'])),0 if side=='L' else 1,0]-z['foot_body_position'][np.argmin(abs(z['t']-old['time_s'])),0 if side=='L' else 1,0])-(row['belt_x_m']-old['belt_x_m']))
  same[side]=row
 for first,second in zip(events[:-1],events[1:]):
  if first['side']!=second['side']:steps.append(second['time_s']-first['time_s'])
 fields=('stance_duration_s','foot_frame_ahead_COM_m','loaded_rear_ahead_COM_m','capture_margin_m','COM_relative_belt_velocity_m_s','relative_velocity_to_put_XCOM_at_loaded_rear_m_s','forward_foot_frame_separation_m');medians={key:float(np.median([row[key] for row in events])) if events else None for key in fields};report=dict(path=str(folder),source_sha256=r['source_sha256'],RG_period_s=periods,median_alternating_contact_interval_s=float(np.median(steps)) if steps else None,median_actual_belt_speed_m_s=float(np.median(-beltv[z['t']>=1])),medians=medians,events=events,scope='Force-qualified landing: normal load >=5%BW for40ms; unload <2.5%BW for50ms. Events after1s. Geometry uses actually loaded foot/toe box projections. Foot-frame separation is not a labelled heel-marker step length. XCoM uses belt-relative velocity and LIPM full-gravity approximation; catch time assumes constant velocity, not a stability proof under rails. No neural/physical parameters changed.');Path(output).write_text(json.dumps(report,indent=2));return report

def summarize_campaign(root):
 root=Path(root);w=root/'W2L_VirtualWalker_Newer';baseline=w/'central_combined_pooledhome_ground_v2';reference=np.load(baseline/'trajectory.npz');searches={};pace=[];controls={}
 for kind in ('body_balance_v1','capture_balance_v1'):
  folder=root/kind/'rustuna_v1'
  if not (folder/'SUMMARY.json').exists():continue
  s=json.loads((folder/'SUMMARY.json').read_text());rows=[json.loads(line) for line in (folder/'TRIALS.jsonl').read_text().splitlines()];searches[kind]=dict(path=str(folder),source_snapshot=str(folder/'ADAPTER_SOURCE.py'),completed_trials=s['completed_trials'],wall_s=s['this_process_wall_s'],physical_passes=sum(bool(row['diagnostics'].get('physical_pass')) for row in rows),best_trial=s['best']['number'],best_release_duration_s=s['best']['diagnostics'].get('release_duration_s'),best_configuration=s['best']['diagnostics']['configuration'],all_ground_gates_pass=sum(bool(row['diagnostics']['ground_gate']) for row in rows),scope='No release candidate promoted. Tiny changes in time-to-fall are not recovery.')
 for tag in ('body_balance_default_path_fresh_v1','body_balance_zero_path_v1','body_balance_matched_zero_v2','capture_balance_zero_control_v1'):
  f=w/tag/'trajectory.npz'
  if not f.exists():continue
  z=np.load(f);checks={}
  for key in reference.files:
   a,b=reference[key],z[key]
   if a.shape!=b.shape:
    if a.ndim==2 and b.shape[1]==a.shape[1]+2:b=b[:,:a.shape[1]]
    elif a.ndim==1 and len(b)==len(a)+2:b=b[:len(a)]
   checks[key]=a.shape==b.shape and a.tobytes()==b.tobytes()
  controls[tag]=dict(path=str(f.parent),all_original_arrays_exact=all(checks.values()),array_checks=checks)
 for f in sorted(w.glob('pace_*/RESULT.json')):
  r=json.loads(f.read_text());a=json.loads((f.parent/'ASSESSMENT.json').read_text());z=np.load(f.parent/'trajectory.npz');bw=r['robot_mass_kg']*9.81;F=np.stack([z[key][:,1] for key in ('rail_force_world','guide_force_world','harness_force')],axis=1);pace.append(dict(path=str(f.parent),belt_command_m_s=r['belt_command_m_s'],mode=r['mode'],duration_s=r['duration_completed_s'],failure=r['failure'],stepping_gate=a['alternating_supported_stepping_gate'],cleared_swings=[a['legs'][s]['unloaded_clear_swings'] for s in 'LR'],upward_channel_BW=float(np.maximum(F,0).sum(1).mean()/bw),downward_channel_BW=float(np.maximum(-F,0).sum(1).mean()/bw),net_assistance_BW=r['mean_net_vertical_assistance_BW'],belt_relative_progress_m=r['root_forward_relative_to_belt_m'],minimum_knees_deg=r['anatomical_knee_min_deg'],slide_limits=r['slide_travel_checks']))
 fresh=w/'pace_0p2_supported10_fresh_v2/trajectory.npz';first=w/'pace_0p2_supported10_v1/trajectory.npz';repeat={}
 if fresh.exists() and first.exists():
  a,b=np.load(first),np.load(fresh);repeat={key:a[key].shape==b[key].shape and a[key].tobytes()==b[key].tobytes() for key in a.files}
 report=dict(machine='easteregg2',preserved_baseline=str(root/'POOLED_PHYSICS_CHECKPOINT.json'),landing_geometry=str(root/'LANDING_GEOMETRY_BASELINE.json'),body_balance_searches=searches,disconnected_controls=controls,pace_conditions=pace,pace_0p2_10s_fresh_array_parity=repeat,independent_walking_achieved=False,scope='Separate experimental architectures and explicit environment conditions. No neural gain improvement, support-system weaning, or independent balance inferred from the belt-speed comparison. Default0.1 baseline preserved; physical force/clearance gates unchanged. Full release still fails. Signed support is positive/negative per-time net channels, not per-contact decomposition.')
 basin=[row for row in pace if Path(row['path']).name.startswith('pace_0p2_basin_')]
 report['pace_0p2_basin']=dict(count=len(basin),passes=sum(row['stepping_gate'] for row in basin),seeds=list(range(100041,100049)),scope='Same sixteen learned scalar perturbations as the preserved0.1m/s basin; paired environment comparison, no fitting.',rows=basin)
 waveform=root/'NEURAL_CLOCK_WAVEFORM_PARITY.json';plan=root/'NEURAL_CLOCK_RELEASE_PLAN.json'
 if waveform.exists() and plan.exists():
  cadence=dict(waveform_parity=json.loads(waveform.read_text()),release_protocol=json.loads(plan.read_text()),supported=[],phase_matched_release=[])
  for tag in ('cadence_fast_belt_0p1_ground_v1','cadence_fast_belt_0p2_ground_v1'):
   folder=w/tag;c=json.loads((folder/'SPIKING_COMPARISON.json').read_text());r=c['result'];a=c['assessment'];cadence['supported'].append(dict(path=str(folder),neural_clock_multiplier=r['neural_clock_multiplier'],belt_speed_m_s=r['belt_command_m_s'],stepping_gate=a['alternating_supported_stepping_gate'],cleared_swings=[a['legs'][s]['unloaded_clear_swings'] for s in 'LR'],net_assistance_BW=r['mean_net_vertical_assistance_BW'],belt_relative_progress_m=r['root_forward_relative_to_belt_m']))
  for row in cadence['release_protocol']['conditions']:
   folder=w/row['tag'];c=json.loads((folder/'SPIKING_COMPARISON.json').read_text());r=c['result'];a=c['assessment'];z=np.load(folder/'trajectory.npz');ref=np.load(Path(row['reference'])/'trajectory.npz');n=int(np.sum(z['t']<row['release_at_s']));checks={key:z[key][:n].shape==ref[key][:n].shape and z[key][:n].tobytes()==ref[key][:n].tobytes() for key in z.files if z[key].ndim and z[key].shape[0]==len(z['t']) and key in ref.files and ref[key].ndim and ref[key].shape[0]==len(ref['t'])};dt=float(np.median(np.diff(z['t'])));actual_release=np.ceil((row['release_at_s']-1e-10)/dt)*dt
   cadence['phase_matched_release'].append(dict(path=str(folder),clock_multiplier=row['clock_scale'],belt_speed_m_s=row['belt_speed'],requested_release_s=row['release_at_s'],actual_release_grid_s=float(actual_release),duration_completed_s=r['duration_completed_s'],survival_after_release_s=float(r['duration_completed_s']-actual_release),stepping_gate=a['alternating_supported_stepping_gate'],cleared_swings=[a['legs'][s]['unloaded_clear_swings'] for s in 'LR'],failure=r['failure'],supported_reference_prefix_exact=all(checks.values()),prefix_array_checks=checks))
  cadence['scope']='Whole neural biological clock, not RG-only kinetics or a performance speed-up. Deafferented matched-biological-time108-neuron waveforms are exactly unchanged. Faster passes air/ground but not release; slower fails air and is not advanced. Release timing is predeclared from supported landings, not chosen after outcomes. Different histories mean this is not identical mechanical release state.';report['cadence']=cadence
 if (root/'motor_output_fit_v1').exists():
  summarize_motor_outputs(root);report['motor_output_continuation']=str(root/'motor_output_fit_v1/MOTOR_OUTPUT_CHECKPOINT.json')
 (root/'BALANCE_CHECKPOINT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');record=root/'BALANCE_PATHWAY_RECORD.md';text=record.read_text(encoding='utf-8').split('<!-- GENERATED_BALANCE_CHECKPOINT -->')[0];text+='\n<!-- GENERATED_BALANCE_CHECKPOINT -->\n## Current durable checkpoint\n\nFull machine-readable summary: `D:\\Users\\Ben Bolen\\NeuroCampaign_20261031\\spiking\\BALANCE_CHECKPOINT.json`. Pitch-only search:24 trials, best2.728 s fall; capture/functional projection:16 trials, best2.790 s fall. Neither completes a single qualifying postrelease alternating sequence. The frozen baseline falls at2.776 s; no new architecture is promoted.\n\nForce-qualified landing geometry at the source0.1m/s belt: RG period0.916s, alternating contact interval0.458s, median loaded rear edge95.2mm ahead of CoM, XCoM marginâˆ’79.3mm. Actual belt-relative CoM speed0.0536m/s versus an LIPM estimate0.328m/s to put XCoM at that rear edge. Stance lasts0.585s; catch-up from the median gap/speed takes about1.78s. Interfoot separation is only30mm; this is mainly a body-behind-feet offset. Stance ankles are aroundâˆ’6.3 degrees plantarflexion, not saturated at the3-degree DF limit.\n\nMatched6s speed conditions with the exact frozen704/80 controller:0.1m/s gives6/6 swings and40.89%BW net assistance;0.2 gives6/5 and33.99%;0.3 gives0/0 and74.82%, failing the actual stepping gate. Downward channel reaction is about0.01%BW at all three, so the0.2 reduction is not downward-cap cancellation. No slide limits are approached in these6s trials. All rail-release probes still fail. This is an environment-condition finding, not neural training or harness weaning.\n\nThe0.2m/s10s case passes10/9 unloaded, cleared swings at36.4754%BW net assistance (up36.4835%, down0.0081%), with+0.517m belt-relative progress. Its fresh process repeats every saved trajectory array exactly. The treadmill slider reaches1.893m, inside its2m limit; this duration is not extended without explicitly handling the finite arena. Closed-loop held-out Â±1% checks of the same16 learned scalars pass'+str(sum(row['stepping_gate'] for row in basin))+'/'+str(len(basin))+'; net assistance spans'+(f"{100*min(row['net_assistance_BW'] for row in basin):.2f}â€“{100*max(row['net_assistance_BW'] for row in basin):.2f}%BW" if basin else 'not yet measured')+'. The original0.1m/s checkpoint remains preserved.\n\nWhole-neural-clock controls keep all currents, conductances, reset/adaptation amplitudes, decoder normalization, and physical Hill/MuJoCo clock fixed. At matched biological times2â€“4s, all108 deafferented neuronal waveforms, including every MN, are byte-exact (max0mV difference) for both5/6 faster and1.25 slower clocks. The faster case passes both air stages and7/6 supported swings at each belt, but requires42.92%BW at0.1 and37.36%BW at0.2, rather than reducing assistance. Slower timing fails the physical air-knee gate despite the unchanged motor waveform, so it is not promoted to ground. This is a mechanical timing response, not neural amplitude suppression.\n\nPhase-matched release is selected in advance from the second right force-qualified landing after one clock-scaled second,2.5% of an RG period before that landing. At belts0.1/0.2, normal-clock survival after actual release is0.776/0.740s; faster-clock survival0.762/0.798s. All four fail the released stepping gate. Every recorded time-series array before release matches its own supported reference exactly. Side guides remain, so even these failed probes are not tests of independent3D walking. Exact timing, rounding to the2ms physics grid, and per-array checks are saved in the checkpoint. No cadence winner or release success is claimed.\n\nSource-mechanics output fitting continues separately in `motor_output_fit_v1/MOTOR_OUTPUT_RECORD.md`; completed exact cases and source gains are in `motor_output_fit_v1/MOTOR_OUTPUT_CHECKPOINT.json`. No output candidate replaces the preserved passing configuration without a fresh measured improvement.\n';record.write_text(text,encoding='utf-8');return report

def summarize_motor_outputs(root):
 root=Path(root);w=root/'W2L_VirtualWalker_Newer';rows=[];reference=w/'pace_0p2_supported10_v1';base=json.loads((reference/'SPIKING_COMPARISON.json').read_text());optional=('ankle_pf_output_multiplier','ankle_dorsiflexor_pf_output_multiplier','hip_extensor_pf_output_multiplier','hip_extensor_pf_activation_floor','knee_flexor_pf_output_multiplier');defaults={}
 for f in sorted(w.glob('motor_*/SPIKING_COMPARISON.json')):
  c=json.loads(f.read_text());r=c['result'];a=c['assessment'];z=np.load(f.parent/'trajectory.npz');bw=r['robot_mass_kg']*9.81;F=np.stack([z[k][:,1] for k in ('rail_force_world','guide_force_world','harness_force')],axis=1);row=dict(path=str(f.parent),mode=c['mode'],feedback=c['feedback'],source_sha256=c['source_sha256'],actual_spiking_cells=r['actual_spiking_cells'],actual_graded_cells=r['actual_graded_cells'],configuration=c['configuration'],changed_output_settings={k:c['configuration'][k] for k in optional if k in c['configuration']},effective_outputs=c['spiking']['output_gain_edits'],air_gate=None if c.get('air_gate') is None else c['air_gate']['passed'],ground_gate=a['alternating_supported_stepping_gate'],cleared_swings=[a['legs'][s]['unloaded_clear_swings'] for s in 'LR'],net_assistance_BW=r['mean_net_vertical_assistance_BW'],upward_channel_BW=float(np.maximum(F,0).sum(1).mean()/bw),downward_channel_BW=float(np.maximum(-F,0).sum(1).mean()/bw),belt_relative_progress_m=r['root_forward_relative_to_belt_m'],minimum_knees_deg=r['anatomical_knee_min_deg'],failure=r['failure'],duration_s=r['duration_completed_s'],implementation_sha256=c['implementation_sha256']);rows.append(row)
  if f.parent.name.startswith('motor_output_default_fresh'):
   ref=np.load(reference/'trajectory.npz');defaults[f.parent.name]={key:ref[key].shape==z[key].shape and ref[key].tobytes()==z[key].tobytes() for key in ref.files}
 report=dict(machine='easteregg2',preserved_baseline=str(reference),baseline_net_support_BW=base['result']['mean_net_vertical_assistance_BW'],source_current_fit=str(reference/'MOTOR_OUTPUT_CURRENT_FIT.json'),actual_chemical_currents=str(root/'motor_output_fit_v1/ACTUAL_CHEMICAL_CURRENTS.json'),mechanical_sensitivities=str(root/'motor_output_fit_v1/TASK_ACTION_SENSITIVITY.json'),default_path_parity=defaults,cases=rows,scope='Source-based output primitive calibration only. Timing, force laws, geometry, initial pose, physical guides, belt0.2, and afferent maps remain frozen. Shared anatomical output gains/thresholds are explicitly recorded. Air gates precede ground. No candidate is promoted solely from joint cycles, force capacity, or lower net support; a fresh repeat and actual clearance/unloading/alternation/forward progression remain required. Not full Hunt training or rail-free walking.')
 study=root/'motor_output_fit_v1/rustuna_v1/SUMMARY.json';fresh=root/'motor_output_fit_v1/KNEE_FRESH_PARITY.json';basin=[]
 for f in sorted(w.glob('motor_knee_basin_*/SPIKING_COMPARISON.json')):
  c=json.loads(f.read_text());basin.append(dict(path=str(f.parent),seed=c['perturbation_seed'],scalar_count=len(c['perturbation_factors']),factors=c['perturbation_factors'],gate=c['assessment']['alternating_supported_stepping_gate'],alternation=c['assessment']['alternating_transition_fraction'],cleared_swings=[c['assessment']['legs'][s]['unloaded_clear_swings'] for s in 'LR'],net_assistance_BW=c['result']['mean_net_vertical_assistance_BW'],belt_relative_progress_m=c['result']['root_forward_relative_to_belt_m'],failure=c['result']['failure']))
 if study.exists():report['rustuna']=dict(path=str(study),source_snapshot=str(study.parent/'ADAPTER_SOURCE.py'),summary=json.loads(study.read_text()))
 if fresh.exists():report['fresh_parity']=json.loads(fresh.read_text())
 report['candidate_basin']=dict(passes=sum(row['gate'] for row in basin),count=len(basin),rows=basin)
 control=w/'motor_knee_baseline17_seed48_v1/SPIKING_COMPARISON.json';failure=w/'motor_knee_basin_100048_v1/SPIKING_COMPARISON.json'
 if control.exists() and failure.exists():
  a,b=json.loads(control.read_text()),json.loads(failure.read_text());air={}
  for mode in ('deaff','aff'):
   f=w/('motor_knee_seed48_air_'+mode+'_v1/SPIKING_COMPARISON.json')
   if f.exists():air[mode]=json.loads(f.read_text())['air_gate']
  report['failed_seed_diagnosis']=dict(baseline_control=str(control),candidate=str(failure),all17_factor_maps_identical=a['perturbation_factors']==b['perturbation_factors'],baseline_ground_gate=a['assessment']['alternating_supported_stepping_gate'],candidate_ground_gate=b['assessment']['alternating_supported_stepping_gate'],candidate_air=air,interpretation='Both isolated air gates pass and the matched baseline passes ground. The new knee output introduces load-transfer/alternation sensitivity, not an intrinsic RG latch. No gate relaxation.')
 report['output_candidate_promoted']=False;report['preferred_configuration']=str(root/'CENTRAL_RESET_REST_CANDIDATE.json');report['promotion_status']='Retain the verified original0.2m/s controller. The nominal knee2.3 support benefit repeats exactly but the17-scalar basin passes only7/8, so the fitted output remains experimental; no release advancement.'
 report['all_motor_current_budgets']=str(reference/'ACTUAL_MOTOR_CHEMICAL_CURRENTS.json');(root/'motor_output_fit_v1/MOTOR_OUTPUT_CHECKPOINT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report

def diagnose(folder):
 from import_virtual_walker import mujoco
 import native_kernel as K
 folder=Path(folder);z=np.load(folder/'trajectory.npz');r=json.loads((folder/'RESULT.json').read_text());model=mujoco.MjModel.from_xml_path(str(folder/'model.xml'));data=mujoco.MjData(model)
 bank=json.loads((folder.parent/'IMPORT_AUDIT.json').read_text())['muscles'];controller=K.NativeController(json.loads((folder/'SPIKING_COMPARISON.json').read_text())['source_configuration']['source'],bank)
 N=len(z['t']);bw=r['robot_mass_kg']*9.81;root=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Root');belt=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Treadmill_Belt');signs=np.array(list(itertools.product((-1,1),repeat=3)))
 geometries=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+name) for name in ('foot_L','toe_L','foot_R','toe_R')];floor=next(g for g in range(model.ngeom) if model.geom_bodyid[g]==belt)
 loaded_bounds=np.full((N,2),np.nan);active_bounds=loaded_bounds.copy();beltx=np.empty(N);height=np.empty(N);pitch=np.empty(N);torques=np.empty((N,6));qprev=np.array(r['initial_qpos']);act=[i for i,b in enumerate(bank) if b['type']=='LinearHillMuscle'];tids=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_TENDON,'route_'+bank[i]['name']) for i in act];jointnames=['hip_L','knee_L','ankle_L','hip_R','knee_R','ankle_R'];dofs=[int(model.jnt_dofadr[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,j)]) for j in jointnames]
 for n,q in enumerate(z['q']):
  data.qpos[:]=qprev;mujoco.mj_kinematics(model,data);mujoco.mj_comPos(model,data);mujoco.mj_tendon(model,data);torques[n]=(-z['forces'][n,act,None]*data.ten_J.reshape(model.ntendon,model.nv)[np.ix_(tids,dofs)]).sum(0)
  data.qpos[:]=q;mujoco.mj_kinematics(model,data);mujoco.mj_comPos(model,data);beltx[n]=data.xpos[belt,0];top=(data.geom_xpos[floor]+(signs*model.geom_size[floor])@data.geom_xmat[floor].reshape(3,3).T)[:,1].max();height[n]=z['com_position'][n,1]-top;R=data.xmat[root].reshape(3,3);pitch[n]=np.arctan2(R[0,1],R[1,1]);full=[];low=[]
  for k,g in enumerate(geometries):
   if z['foot_region_force_world'][n,k,1]<=.01*bw or z['foot_vertical'][n,k//2]<=.05*bw:continue
   vertices=data.geom_xpos[g]+(signs*model.geom_size[g])@data.geom_xmat[g].reshape(3,3).T;full.extend(vertices[:,0]);minimum=vertices[:,1].min();low.extend(vertices[vertices[:,1]<=minimum+.002,0])
  if full:loaded_bounds[n]=[min(full),max(full)];active_bounds[n]=[min(low),max(low)]
  qprev=q
 beltv=np.gradient(beltx,z['t']);omega=np.sqrt(9.81/np.maximum(.1,height));capture=z['com_position'][:,0]+(z['com_velocity'][:,0]-beltv)/omega;com=z['com_position'][:,0];margin=np.minimum(com-active_bounds[:,0],active_bounds[:,1]-com);cp_margin=np.minimum(capture-active_bounds[:,0],active_bounds[:,1]-capture)
 release=json.loads((folder/'SPIKING_COMPARISON.json').read_text()).get('rail_release_at_s');release=2. if release is None else release
 phase=np.where(z['neuro'][:,controller.names.index('L RG ext')]>z['neuro'][:,controller.names.index('L RG flx')],'L_RG_E','L_RG_F');support=np.where(np.all(z['foot_vertical']>.05*bw,axis=1),'double',np.where(z['foot_vertical'][:,0]>.05*bw,'left',np.where(z['foot_vertical'][:,1]>.05*bw,'right','none')))
 def stats(mask):
  if not mask.any():return None
  return dict(samples=int(mask.sum()),pitch_mean_deg=float(np.degrees(pitch[mask]).mean()),pitch_rate_mean_deg_s=float(np.degrees(-z['pelvis_omega_world'][mask,2]).mean()),COM_velocity_relative_belt_mean_m_s=float((z['com_velocity'][mask,0]-beltv[mask]).mean()),COM_loaded_margin_mean_m=float(np.nanmean(margin[mask])) if np.isfinite(margin[mask]).any() else None,capture_loaded_margin_mean_m=float(np.nanmean(cp_margin[mask])) if np.isfinite(cp_margin[mask]).any() else None,COM_behind_support_fraction=float(np.mean(com[mask]<active_bounds[mask,0])),capture_behind_support_fraction=float(np.mean(capture[mask]<active_bounds[mask,0])),foot_load_BW=(z['foot_vertical'][mask].mean(0)/bw).tolist(),joint_muscle_torque_mean_Nm=torques[mask].mean(0).tolist(),rail_BW=float(z['rail_force_world'][mask,1].mean()/bw),guide_BW=float(z['guide_force_world'][mask,1].mean()/bw))
 windows={}
 for lo,hi in [(-.2,0),(0,.15),(.15,.35),(.35,.55),(.55,1.)]:windows[f'{lo:+.2f}_{hi:+.2f}']=stats((z['t']>=release+lo)&(z['t']<release+hi))
 phases={p+'_'+s:stats((z['t']>=release)&(phase==p)&(support==s)) for p in ('L_RG_E','L_RG_F') for s in ('left','right','double','none')}
 result=dict(path=str(folder),foot_load_units='BW',release_s=release,source_sha256=r['source_sha256'],joint_order=jointnames,windows_after_release_s=windows,neural_support_phase=phases,scope='Saved actual COM/velocity, foot region forces, and source muscle tensions. Tendon Jacobians use pre-step pose; geometry footprint uses saved post-step pose. Loaded region threshold1%BW and leg threshold5%BW. Active footprint is2mm lowest vertex band of loaded foot/toe boxes, a geometric approximation rather than exact force-patch boundary. Capture point uses moving-belt-relative COM velocity and linear inverted-pendulum approximation, not an independent dynamic stability proof. All torques are actual muscle moments, never applied targets.')
 (folder/'BALANCE_DIAGNOSIS.json').write_text(json.dumps(result,indent=2));np.savez_compressed(folder/'BALANCE_MECHANICS.npz',t=z['t'],loaded_bounds=loaded_bounds,active_bounds=active_bounds,COM_x=com,capture_x=capture,COM_margin=margin,capture_margin=cp_margin,pitch_rad=pitch,belt_velocity_x=beltv,joint_torque_Nm=torques,phase=phase,support=support);return result

if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument('folders',nargs='+');parser.add_argument('--fit-budget',action='store_true');parser.add_argument('--fit-task-actions',action='store_true');parser.add_argument('--fit-functional-budget',action='store_true');parser.add_argument('--landing',action='store_true');parser.add_argument('--summary',action='store_true');parser.add_argument('--output');args=parser.parse_args()
 for folder in args.folders:print(json.dumps(summarize_campaign(folder) if args.summary else landing_geometry(folder,args.output) if args.landing else fit_functional_current_budget(folder,args.output) if args.fit_functional_budget else fit_task_actions(folder,args.output) if args.fit_task_actions else fit_current_budget(folder,args.output) if args.fit_budget else diagnose(folder),indent=2))
