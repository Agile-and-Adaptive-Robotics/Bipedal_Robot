"""Experimental axis-specific and biarticular PF microcircuits for Gait2392.

The w2lvar RG and reflex pathways are retained. Extra axis pairs have two
excitatory PF cells and two inhibitory relays. Dedicated biarticular pairs
are stance/swing pattern channels, not a claim that every muscle has a
single anatomical antagonist. No MN current cancellation is introduced.
"""
from dataclasses import dataclass,field
from pathlib import Path
import json
import numpy as np
import build_network as BN
import build_network_w2lvar as W
from mechanics import JOINTS,CoupledArms,BI_FAMILIES
from params import TAU,G,E_HI
from sns_toolbox import backends

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'reports_20261028_functional_pf'
DEFAULTS=dict(axes=False,biarticular=False,pose_dependent=True,axis_gain=.1,
              axis_feedback=1.,trunk_gain=.1,bi_gain=1.,bi_afferent=.03,update_steps=5)

@dataclass
class FunctionalNetwork(W.W2LVarnet):
    options:dict=field(default_factory=dict)
    routing:dict=field(default_factory=dict)
    extra_pairs:dict=field(default_factory=dict)
    legacy_pf_budget:dict=field(default_factory=dict)
    extra_edges:list=field(default_factory=list)
    live:object=None

    def _pair(self,n,label,side,rg=False):
        key=label+'_'+side
        cells=[f'PF_FUNC_{label}_{phase}_{side}' for phase in ('P','N')]
        relays=[f'PF_IN_FUNC_{label}_{phase}_{side}' for phase in ('P','N')]
        for name in cells+relays:self._add(name,TAU['pf'],n)
        for i in (0,1):
            n.add_connection(BN._syn(W.GAINS['w2l_pf_recip'],True),cells[i],relays[i])
            n.add_connection(BN._syn(W.GAINS['w2l_pf_recip'],False),relays[i],cells[1-i])
            if rg:
                n.add_connection(BN._syn(W.GAINS['w2l_rg_to_pf'],True),f'RG_{"E" if i==0 else "F"}_{side}',cells[i])
            else:
                n.add_input(cells[i]);self.inputs.append('FUNC_'+key+'_'+str(i))
        self.extra_pairs[key]=cells

    def _build_pf_w2l(self,n,side):
        super()._build_pf_w2l(n,side)
        if self.options['axes']:
            for label in ('HIP_X','HIP_Y','FOOT_X'):self._pair(n,label,side)
            if side==self.sides[0]:
                for label in ('LUM_X','LUM_Y','LUM_Z'):self._pair(n,label,'mid')
        if self.options['biarticular']:
            self._pair(n,'ANKLE_Z',side,rg=True)  # P plantarflexion, N dorsiflexion
            for family in BI_FAMILIES:
                if any(v['family']==family and a.endswith('_'+side) for a,v in self.routing['actuator_mapping'].items()):
                    self._pair(n,'BI_'+family,side,rg=True)
            # The old shared toe pathway suppressed knee flexion too. With
            # separate ankle PF, toe inhibition belongs on dorsiflexion alone.
            names=[p['name'] for p in n.populations]
            for c in n.connections:
                if names[c['source']]=='TOEDF_'+side and names[c['destination']]=='PF_KNEE-F_'+side:
                    c['destination']=names.index(self.extra_pairs['ANKLE_Z_'+side][1])

    def _edge(self,n,source,target,g,record=None):
        # Tiny nonzero compile-time reservations let pose-dependent sign changes
        # enable either direction without changing compiled topology mid-trial.
        n.add_connection(BN._syn(max(g,1e-12),True),source,target)
        if record is not None:self.extra_edges.append((source,target,*record))

    def _wire_muscle_w2l(self,n,act,mi):
        super()._wire_muscle_w2l(n,act,mi)
        names=[p['name'] for p in n.populations];mn='MN_'+act
        old=[c for c in n.connections if names[c['destination']]==mn and names[c['source']].startswith('PF_') and not names[c['source']].startswith('PF_IN')]
        budget=sum(float(c['params']['max_conductance']) for c in old)
        self.legacy_pf_budget[act]=budget
        info=self.routing['actuator_mapping'][act];family=info['family']
        if self.options['biarticular'] and family:
            for c in old:n.connections.remove(c)
            cells=self.extra_pairs['BI_'+family+'_'+mi.side]
            weights=self.routing['phase_weights'][act]
            for cell,phase in zip(cells,('stance','swing')):
                self._edge(n,cell,mn,self.options['bi_gain']*budget*weights[phase])
            # Retain autogenic feedback, now reaching this muscle's own PF pair.
            for aff in ('Ia','II'):
                n.add_connection(BN._syn(self.options['bi_afferent'],True),aff+'_'+act,cells[1])
            # Phase gate of Ia reciprocal inhibition follows the new swing PF.
            for c in n.connections:
                if names[c['destination']]=='IaIN_'+act and names[c['source']].startswith('PF_'):
                    c['source']=names.index(cells[1])
        elif self.options['biarticular'] and mi.groups[0] in ('ankle_pf','ankle_df'):
            for c in old:n.connections.remove(c)
            cells=self.extra_pairs['ANKLE_Z_'+mi.side]
            self._edge(n,cells[0 if mi.groups[0]=='ankle_pf' else 1],mn,budget)
        if self.options['axes']:
            side=mi.side;leg=0 if side=='r' else 7
            if mi.groups[0].startswith('trunk'):
                entries=[('LUM_X_mid',15),('LUM_Y_mid',16),('LUM_Z_mid',14)]
            else:entries=[('HIP_X_'+side,leg+1),('HIP_Y_'+side,leg+2),('FOOT_X_'+side,leg+5)]
            for key,row in entries:
                for sign,cell in zip((1.,-1.),self.extra_pairs[key]):
                    self._edge(n,cell,mn,1e-12,record=(act,row,sign,key))

    def bind_mujoco(self,model,data):
        self.live=(model,data);self.arm_model=CoupledArms(model)
        self.act_order=[__import__('mujoco').mj_id2name(model,__import__('mujoco').mjtObj.mjOBJ_ACTUATOR,i) for i in range(model.nu)]
        self.act_index={a:i for i,a in enumerate(self.act_order)}
        self.goal=data.qpos[self.arm_model.addresses].copy();self.step_count=0
        self.ports={name:i for i,name in enumerate(self.inputs)}
        self.change_log=[]
        samples=np.load(OUT/'routing_samples.npz')
        rho=np.concatenate((samples['ik_moment_arms'],samples['sweep_moment_arms']))
        self.row_scale=np.maximum(np.max(abs(rho),axis=(0,2)),.003)
        self._update_weights(force=True)

    def _update_weights(self,force=False):
        if not self.options['axes']:return
        if not force and not self.options['pose_dependent']:return
        _,data=self.live;rho=self.arm_model.update(data)
        c=self.compiled
        if not hasattr(self,'edge_locations'):
            self.edge_locations=[]
            for source,target,act,row,sign,key in self.extra_edges:
                post,pre=self.idx[target],self.idx[source]
                matches=np.flatnonzero((c.post==post)&(c.pre==pre))
                assert len(matches)==1,(source,target,matches)
                self.edge_locations.append((post,pre,int(matches[0]),self.act_index[act],row,sign,key))
        for post,pre,edge,ai,row,sign,key in self.edge_locations:
            gain=self.options['trunk_gain'] if key.startswith('LUM') else self.options['axis_gain']
            value=gain*np.clip(sign*rho[row,ai]/self.row_scale[row],0.,1.)
            c.c.g_max_non[post,pre]=value;c.gmax[edge]=value
        if self.step_count%100==0:
            self.change_log.append(dict(time_s=float(data.time),hip_flexion_r=float(data.qpos[self.arm_model.addresses[0]]),
                active_functional_synapses=int(np.count_nonzero(c.gmax[[r[2] for r in self.edge_locations]]))))

    def step(self,input_currents):
        if self.live is not None and self.options['axes']:
            _,data=self.live
            if self.step_count%self.options['update_steps']==0:self._update_weights()
            q=data.qpos[self.arm_model.addresses];v=data.qvel[self.arm_model.rows]
            demand=np.clip((self.goal-q)/.2-v/.8,-2.,2.)*self.options['axis_feedback']
            for key,cells in self.extra_pairs.items():
                if key.startswith(('BI_','ANKLE_Z')):continue
                side=key[-1];leg=0 if side=='r' else 7
                row=15 if key.startswith('LUM_X') else 16 if key.startswith('LUM_Y') else 14 if key.startswith('LUM_Z') else leg+(1 if key.startswith('HIP_X') else 2 if key.startswith('HIP_Y') else 5)
                for i,sign in enumerate((1.,-1.)):
                    input_currents[self.ports['FUNC_'+key+'_'+str(i)]]=max(0.,sign*demand[row])
        self.compiled.forward(input_currents);self.step_count+=1
        return self.compiled.V

def build(acts,dt=.002,interleg=True,options=None):
    options={**DEFAULTS,**(options or {})}
    if not options['axes'] and not options['biarticular']:
        return W.build(acts,dt=dt,interleg=interleg)
    for key,value in W.VARIANT_G.items():G[key]=float(value)
    muscles={a:W.classify(a) for a in acts};assert all(muscles.values())
    routing=json.loads((OUT/'routing_audit.json').read_text())
    net=FunctionalNetwork(muscles=muscles,sides=tuple(sorted({m.side for m in muscles.values()})),
                          interleg=interleg,options=options,routing=routing)
    net.compile(dt)
    return net
