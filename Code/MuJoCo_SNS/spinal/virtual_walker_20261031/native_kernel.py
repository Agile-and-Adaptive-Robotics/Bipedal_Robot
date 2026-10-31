"""Source-parameter, nonspiking IntegrateFireSim/Hill port with compiled sparse edges.

Equations independently transcribed from AnimatLabPublicSource's Neuron.cpp,
IntegrateFireModule.cpp, LinearHillMuscle.cpp, and gain/receptor classes.
This deliberately does NOT substitute the generic calibrated SNS W2L circuit.
"""
import math, json, xml.etree.ElementTree as E
import numpy as np
from numba import njit
from import_virtual_walker import number

@njit(cache=True)
def neural_steps(v,m,h,newv,gs,gi,params,edges,tonic,applied,dt_ms,steps):
 for tick in range(steps):
  gs[:]=0.;gi[:]=0.
  for k in range(len(edges)):
   pre=int(edges[k,0]);post=int(edges[k,1]);g=edges[k,2];threshold=edges[k,3];sat=edges[k,4]
   if v[pre]<=threshold:g=0.
   elif v[pre]<sat:g*= (v[pre]-threshold)/(sat-threshold)
   gs[post]+=g;gi[post]+=g*(edges[k,5]-params[post,0])
  for i in range(len(v)):
   rest,tau,size,gmax,vm,sm,tm,vh,sh,th,eca=params[i]
   gca=gmax*m[i]*h[i];ica=gca*(eca-rest)
   if gmax>0:
    z=math.exp(-sm*(v[i]-vm));inf=1/(1+z);t=inf*math.sqrt(z)*tm;m[i]+=(inf-m[i])*(-math.expm1(-dt_ms/t))
    z=math.exp(-sh*(v[i]-vh))*.5;inf=1/(1+z);t=inf*math.sqrt(z)*th;h[i]+=(inf-h[i])*(-math.expm1(-dt_ms/t))
   total=1.+gs[i]+gca;decay=math.exp(-total*dt_ms/tau)
   current=(tonic[i]+applied[i])/size+gi[i]+ica
   newv[i]=rest+(v[i]-rest)*decay+current*(1-decay)/total
  v[:]=newv

@njit(cache=True)
def hill_step(length,velocity,vm,params,internal,force,ia,ii,ib,dt):
 for i in range(len(length)):
  l0,width,kse,kpe,damping,maximum,A,B,C,D,pe,minpe,ia_gain,ii_gain,ib_gain=params[i]
  tl=max(0.,1.-((length[i]-l0)/width)**2)
  active=max(0.,B/(1.+math.exp(min(700.,C*(A-vm[i]))))+D)*tl
  internal[i]+=dt*kse/damping*(kpe*(length[i]-l0)+damping*velocity[i]-(1+kpe/kse)*internal[i]+active)
  force[i]=min(maximum,max(0.,internal[i]))
  series=(1.-pe*.01)*l0+force[i]/kse
  parallel=length[i]-series;minimum=minpe*.01*l0
  if parallel<minimum:parallel=minimum;series=length[i]-minimum
  ia[i]=ia_gain*series;ii[i]=ii_gain*parallel;ib[i]=ib_gain*force[i]

def polynomial(coeff,x):return ((coeff[0]*x+coeff[1])*x+coeff[2])*x+coeff[3]

class NativeController:
 def __init__(self,source,bank,contact_si=False,pf_knee_flex_scale=1.,ankle_pf_crossed=False,symmetric_rg_ext_tonic=None,neural_time_scale=1.,v3_output_scale=1.):
  root=E.parse(source).getroot();module=next(x for x in root.iter('NeuralModule') if x.findtext('ModuleName')=='IntegrateFireSim')
  assert module.findtext('HH')=='False' and module.findtext('Cd')=='False'
  nodes=list(module.iter('Neuron'));self.names=[n.findtext('Name') for n in nodes];self.ids={n.findtext('ID'):i for i,n in enumerate(nodes)}
  self.params=np.array([[number(n,'RestingPot'),number(n,'TimeConst'),number(n,'Size'),number(n,'GMaxCa'),number(n.find('CaActivation'),'MidPoint'),number(n.find('CaActivation'),'Slope'),number(n.find('CaActivation'),'TimeConstant'),number(n.find('CaDeactivation'),'MidPoint'),number(n.find('CaDeactivation'),'Slope'),number(n.find('CaDeactivation'),'TimeConstant'),number(module,'CaEquilPot')] for n in nodes])
  assert neural_time_scale>0;self.params[:,[1,6,9]]*=neural_time_scale
  assert all(number(n,'Noise')==0 and number(n,'InitialThresh')>=50 for n in nodes),'Port restricted to these graded/no-spike circuits'
  assert not list(module.iter('IonChannel')),'Additional ion channels need a separately validated port'
  self.v=self.params[:,0].copy();self.m=np.zeros(len(nodes));self.h=self.m.copy()
  for i,n in enumerate(nodes):
   if number(n,'GMaxCa')>0 and n.findtext('BurstInitAtBottom')=='True':self.v[i]+=7.408;self.m[i]=.0945;self.h[i]=.0208
  self.newv=self.v.copy();self.gs=np.zeros(len(nodes));self.gi=self.gs.copy();self.applied=self.gs.copy();self.tonic=np.array([number(n,'TonicStimulus') for n in nodes])
  if symmetric_rg_ext_tonic is not None:
   for side in ('L','R'):self.tonic[self.names.index(side+' RG ext')]=symmetric_rg_ext_tonic
  types={x.findtext('ID'):x for x in module.findall('Synapses/NonSpikingSynapses/SynapseType')};edges=[]
  for x in module.iter('Connexion'):
   assert x.findtext('Type')=='1' and number(x,'Delay')==0,'Spiking/electrical/delayed synapse not supported by graded port'
   t=types[x.findtext('SynapseTypeID')];edges.append([self.ids[x.findtext('SourceID')],self.ids[x.findtext('TargetID')],number(t,'SynAmp'),number(t,'ThreshV'),number(t,'SaturateV'),number(t,'Equil')])
  self.edges=np.array(edges);self.tuned_edges=[]
  assert v3_output_scale>=0.;self.v3_output_edits=[]
  for k,e in enumerate(self.edges):
   if (self.names[int(e[0])],self.names[int(e[1])]) in [('L V3','R RG ext IN'),('R V3','L RG ext IN')]:self.v3_output_edits.append(dict(source=self.names[int(e[0])],target=self.names[int(e[1])],source_gmax=float(e[2]),scale=v3_output_scale));self.edges[k,2]*=v3_output_scale
  if v3_output_scale!=1.:assert len(self.v3_output_edits)==2,'This diagnostic requires exactly two mirrored V3 output edges'
  for k,edge in enumerate(self.edges):
   pre,post=int(edge[0]),int(edge[1])
   if any(self.names[pre]==s+' Knee PF flx' and self.names[post]==s+' Knee MN flx' for s in ('L','R')):
    self.tuned_edges.append(dict(source=self.names[pre],target=self.names[post],source_gmax=float(edge[2]),scale=pf_knee_flex_scale));self.edges[k,2]*=pf_knee_flex_scale
  assert len(self.tuned_edges)==2,'Expected one shared bilateral PF→knee-flexor MN parameter'
  self.ankle_phase_edits=[]
  if ankle_pf_crossed:
   for k,edge in enumerate(self.edges):
    pre,post=int(edge[0]),int(edge[1]);pre_name=self.names[pre];post_name=self.names[post]
    if ' Knee PF ' in pre_name and ' Ank MN ' in post_name and not pre_name.endswith(' IN') and not post_name.endswith(' RE'):
     target=post_name.replace(' ext',' flx') if post_name.endswith(' ext') else post_name.replace(' flx',' ext');self.ankle_phase_edits.append(dict(source=pre_name,old_target=post_name,new_target=target,gmax_uS=float(edge[2])));self.edges[k,1]=self.names.index(target)
   assert len(self.ankle_phase_edits)==4,'Expected four bilateral knee-PF→ankle-MN drive edges'
  self.dt=number(module,'TimeStep');self.stimuli=[]
  self.test_pulses=[];self.last_test_impulse_nA_s=np.zeros(len(nodes))
  for s in root.findall('.//ExternalStimuli/Stimulus'):
   if s.findtext('Type')=='Current' and s.findtext('Enabled')=='True':
    assert s.findtext('CurrentType')=='Tonic';self.stimuli.append((self.ids[s.findtext('TargetNodeID')],number(s,'StartTime'),number(s,'EndTime'),number(s,'CurrentOn')*1e9))
  bankids={b['id']:i for i,b in enumerate(bank)};self.ibgain=np.array([b['ib_discharge'] for b in bank]);self.outputs=[];self.inputs=[]
  bodies={b.findtext('ID'):b.findtext('Name') for b in root.iter('RigidBody')}
  for a in root.iter('Adapter'):
   if a.findtext('Enabled')=='False':continue
   g=a.find('Gain');assert g.findtext('Type')=='Polynomial' and g.findtext('UseLimits')=='False'
   coeff=np.array([number(g,k) for k in 'ABCD'])
   if a.findtext('SourceDataType')=='ContactCount' and contact_si:coeff[:3]*=1e-9
   if a.findtext('Type')=='NodeToPhysical':
    assert a.findtext('SourceDataType')=='MembraneVoltage';self.outputs.append((self.ids[a.findtext('SourceID')],bankids[a.findtext('TargetID')],coeff))
   else:
    assert a.findtext('TargetDataType')=='ExternalCurrent';sid=a.findtext('SourceID');kind=a.findtext('SourceDataType');self.inputs.append((self.ids[a.findtext('TargetID')],bankids.get(sid,-1),kind,coeff,bodies.get(sid,'')))
 def step(self,time,ia,ii,ib,contacts,interval,feedback=True):
  self.applied[:]=0
  self.last_test_impulse_nA_s[:]=0.
  if feedback:
   for target,idx,kind,coeff,body in self.inputs:
    if kind=='Ia':value=ia[idx]
    elif kind=='II':value=ii[idx]
    elif kind=='Ib':value=ib[idx]
    elif kind=='Tension':value=ib[idx]/self.ibgain[idx]
    elif kind in ('ContactCount','Contacts'):value=contacts.get(body,0.)
    else:raise ValueError('Unported sensory signal '+kind)
    self.applied[target]+=polynomial(coeff,value)*1e9
  count=round(interval/self.dt);assert abs(count*self.dt-interval)<1e-12
  for sub in range(count):
   stimtime=time+sub*self.dt
   for idx,start,end,amplitude in self.stimuli:
    if start<=stimtime<end:self.applied[idx]+=amplitude
   for idx,start,end,amplitude in self.test_pulses:
    if start<=stimtime<end:self.applied[idx]+=amplitude;self.last_test_impulse_nA_s[idx]+=amplitude*self.dt
   neural_steps(self.v,self.m,self.h,self.newv,self.gs,self.gi,self.params,self.edges,self.tonic,self.applied,self.dt*1000,1)
   for idx,start,end,amplitude in self.test_pulses:
    if start<=stimtime<end:self.applied[idx]-=amplitude
   for idx,start,end,amplitude in self.stimuli:
    if start<=stimtime<end:self.applied[idx]-=amplitude
  assert np.all(self.v<40),'A spike-capable trajectory needs the spiking kernel, not this graded port'
  return self.v
 def motor_voltages(self,vm):
  for pre,target,coeff in self.outputs:vm[target]=polynomial(coeff,self.v[pre]*.001)
 def rg_current_terms(self):
  out=[]
  for post,name in enumerate(self.names):
   if ' RG ' not in name and name not in ('L V3','R V3','L c1','R c1'):continue
   assert self.params[post,2]==1.,'Native-equivalent nA diagnostic restricted to these Size1 source cells'
   incoming=[]
   for e in self.edges:
    if int(e[1])!=post:continue
    pre=int(e[0]);g=float(e[2]*np.clip((self.v[pre]-e[3])/(e[4]-e[3]),0,1));incoming.append(dict(source=self.names[pre],g_active=g,reversal_mV=float(e[5]),current_nA=float(g*(e[5]-self.v[post]))))
   leak=float(self.params[post,0]-self.v[post]);ca=float(self.params[post,3]*self.m[post]*self.h[post]*(self.params[post,10]-self.v[post]));tonic=float(self.tonic[post]);afferent=float(self.applied[post]);syn=sum(x['current_nA'] for x in incoming);out.append(dict(target=name,voltage_mV=float(self.v[post]),leak_nA=leak,calcium_nA=ca,tonic_nA=tonic,physical_applied_nA=afferent,synaptic_nA=syn,net_nA=leak+ca+tonic+afferent+syn,incoming=sorted(incoming,key=lambda x:abs(x['current_nA']),reverse=True)))
  return out

def bank_parameters(bank):
 return np.array([[b[k] for k in ('resting_length','lwidth','kse','kpe','B','maximum_tension')]+b['stimulus']+[b[k] for k in ('pe_length','min_pe_length','ia_discharge','ii_discharge','ib_discharge')] for b in bank])
