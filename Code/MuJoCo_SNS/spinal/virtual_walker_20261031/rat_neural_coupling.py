"""Rat-specific source receptor/adaptor coupling; no generic W2L substitution.

The strict Hunt-B gains are loaded separately and never fitted here. Source
bounded polynomial gain limits, disabled adapters, intrinsic states, and SI to
nA conversion are retained. No unsupported source signal is silently dropped.
"""
from pathlib import Path
import json,xml.etree.ElementTree as E
import numpy as np
from numba import njit
from native_kernel import neural_steps


@njit(cache=True)
def replay_steps(tape, params, edges, tonic, v0, m0, h0, dt_ms):
    """Preallocated source-kernel replay; timing does not depend on physics."""
    v = v0.copy(); m = m0.copy(); h = h0.copy()
    out = np.zeros((len(tape), len(v))); new = v.copy()
    gs = np.zeros(len(v)); gi = gs.copy()
    for k in range(len(tape)):
        neural_steps(v, m, h, new, gs, gi, params, edges, tonic, tape[k], dt_ms, 1)
        out[k] = v
    return out


@njit(cache=True)
def bounded_polynomial(coeff,limits,value):
    if limits[0]>0:
        if value<limits[1]:return limits[3]
        if value>limits[2]:return limits[4]
    return ((coeff[0]*value+coeff[1])*value+coeff[2])*value+coeff[3]


@njit(cache=True)
def receptor_values(length,force,params,ia,ii,ib,series,parallel):
    for i in range(len(length)):
        p=params[i];se=(1.-p[10]*.01)*p[0]+force[i]/p[2]
        pe=length[i]-se;minimum=p[11]*.01*p[0]
        if pe<minimum:pe=minimum;se=length[i]-minimum
        series[i]=se;parallel[i]=pe;ia[i]=p[12]*se;ii[i]=p[13]*pe;ib[i]=p[14]*force[i]


@njit(cache=True)
def adapter_currents(targets,bodies,kinds,coeff,limits,ia,ii,ib,current):
    current[:]=0.
    for i in range(len(targets)):
        value=ia[bodies[i]]if kinds[i]==0 else ii[bodies[i]]if kinds[i]==1 else ib[bodies[i]]
        current[targets[i]]+=bounded_polynomial(coeff[i],limits[i],value)*1e9


@njit(cache=True)
def motor_voltages(v,pre,body,coeff,limits,out):
    out[:]=-.15
    for i in range(len(pre)):
        out[body[i]]=bounded_polynomial(coeff[i],limits[i],v[pre[i]]*.001)


def all_source_hill_bank(source):
    root=E.parse(source).getroot();bank=[]
    for b in root.iter('RigidBody'):
        typ=b.findtext('Type')
        if typ not in('LinearHillMuscle','LinearHillStretchReceptor'):continue
        def val(key,default=None):return float(b.findtext(key,default))
        p=[val('LengthTension/RestingLength'),val('LengthTension/Lwidth'),val('Kse'),val('Kpe'),val('B'),val('MaximumTension')]+[val('StimulusTension/'+q)for q in'ABCD']+[val('LengthTension/PeLength'),val('LengthTension/MinPeLength'),val('IaDischarge','100'),val('IIDischarge','100'),val('IbDischarge','100')]
        assert b.findtext('Enabled')!='False','Disabled source muscles require their exact sensor reset semantics'
        if typ=='LinearHillStretchReceptor':assert b.findtext('ApplyTension','False')=='False','Force-producing receptor not implemented'
        bank.append({'name':b.findtext('Name'),'id':b.findtext('ID'),'type':typ,'params':p})
    return bank,np.array([x['params']for x in bank])


class RatController:
    def __init__(self,source,strict_fit):
        root=E.parse(source).getroot();module=next(m for m in root.iter('NeuralModule')if m.findtext('ModuleName')=='IntegrateFireSim');nodes=list(module.iter('Neuron'));self.names=[n.findtext('Name')for n in nodes];ids={n.findtext('ID'):i for i,n in enumerate(nodes)};self.dt=float(module.findtext('TimeStep'))
        assert module.findtext('HH')=='False'and module.findtext('Cd')=='False';assert len(nodes)==132
        self.params=np.array([[float(n.findtext(q))for q in('RestingPot','TimeConst','Size','GMaxCa','CaActivation/MidPoint','CaActivation/Slope','CaActivation/TimeConstant','CaDeactivation/MidPoint','CaDeactivation/Slope','CaDeactivation/TimeConstant')]+[float(module.findtext('CaEquilPot'))]for n in nodes]);self.tonic=np.array([float(n.findtext('TonicStimulus'))for n in nodes]);assert all(float(n.findtext('Noise'))==0 and float(n.findtext('InitialThresh'))>=50 for n in nodes)
        types={t.findtext('ID'):t for t in module.findall('Synapses/NonSpikingSynapses/SynapseType')};edges=[];self.edge_ids=[]
        for c in module.iter('Connexion'):
            assert c.findtext('Type')=='1'and float(c.findtext('Delay','0'))==0;t=types[c.findtext('SynapseTypeID')];edges.append([ids[c.findtext('SourceID')],ids[c.findtext('TargetID')]]+[float(t.findtext(q))for q in('SynAmp','ThreshV','SaturateV','Equil')]);self.edge_ids.append(c.findtext('ID'))
        self.source_edges=np.array(edges);assert len(edges)==270;fit=np.load(strict_fit);assert self.names==fit['neuron_names'].tolist();assert np.array_equal(self.source_edges,fit['source_edges']);self.edges=fit['fitted_edges'].copy();changed=np.flatnonzero(np.any(self.edges!=self.source_edges,axis=1));assert len(changed)==14
        self.v=self.params[:,0].copy();self.m=np.zeros(len(nodes));self.h=self.m.copy()
        for i,n in enumerate(nodes):
            if self.params[i,3]>0 and n.findtext('BurstInitAtBottom')=='True':self.v[i]+=7.408;self.m[i]=.0945;self.h[i]=.0208
        self.newv=self.v.copy();self.gs=np.zeros(len(nodes));self.gi=self.gs.copy();self.applied=self.gs.copy();self.bank,self.bank_params=all_source_hill_bank(source);bids={b['id']:i for i,b in enumerate(self.bank)};ins=[];outs=[]
        for a in root.iter('Adapter'):
            if a.findtext('Enabled')=='False':continue
            g=a.find('Gain');assert g.findtext('Type')=='Polynomial';coeff=[float(g.findtext(q))for q in'ABCD'];limited=g.findtext('UseLimits')=='True';limits=[1. if limited else 0.]+([float(g.findtext(q))for q in('LowerLimit','UpperLimit','LowerOutput','UpperOutput')]if limited else[0.,0.,0.,0.])
            if a.findtext('Type')=='NodeToPhysical':
                assert a.findtext('SourceDataType')=='MembraneVoltage';assert a.findtext('TargetDataType')=='MembraneVoltage';outs.append((ids[a.findtext('SourceID')],bids[a.findtext('TargetID')],coeff,limits,a.findtext('ID')))
            else:
                kind=a.findtext('SourceDataType');assert kind in('Ia','II','Ib'),kind;assert a.findtext('TargetDataType')=='ExternalCurrent';bi=bids[a.findtext('SourceID')]
                if kind in('Ia','II'):assert self.bank[bi]['type']=='LinearHillStretchReceptor'
                ins.append((ids[a.findtext('TargetID')],bi,('Ia','II','Ib').index(kind),coeff,limits,a.findtext('ID')))
        assert len(ins)==44 and len(outs)==16
        self.input_targets=np.array([x[0]for x in ins],dtype=np.int64);self.input_bodies=np.array([x[1]for x in ins],dtype=np.int64);self.input_kinds=np.array([x[2]for x in ins],dtype=np.int64);self.input_coeff=np.array([x[3]for x in ins]);self.input_limits=np.array([x[4]for x in ins]);self.output_pre=np.array([x[0]for x in outs],dtype=np.int64);self.output_bodies=np.array([x[1]for x in outs],dtype=np.int64);self.output_coeff=np.array([x[2]for x in outs]);self.output_limits=np.array([x[3]for x in outs]);self.inputs=ins;self.outputs=outs
        # Exact strict-B teaching run cleared external pulses; retain that
        # configuration, rather than silently reinstate source perturbations.
        self.excluded_source_stimuli=[{'id':s.findtext('ID'),'name':s.findtext('Name'),'type':s.findtext('Type'),'enabled':s.findtext('Enabled')}for s in root.findall('.//ExternalStimuli/Stimulus')if s.findtext('Enabled')=='True'];self.changed_gain_ids=[self.edge_ids[k]for k in changed]

    def step(self,applied):
        self.applied[:]=applied
        neural_steps(self.v,self.m,self.h,self.newv,self.gs,self.gi,self.params,self.edges,self.tonic,self.applied,self.dt*1000,1)
        return self.v

    def feedback(self,ia,ii,ib):
        adapter_currents(self.input_targets,self.input_bodies,self.input_kinds,self.input_coeff,self.input_limits,ia,ii,ib,self.applied)
        return self.applied

    def voltages(self,out):motor_voltages(self.v,self.output_pre,self.output_bodies,self.output_coeff,self.output_limits,out)
