"""Opt-in pure-input GLIF -> actual graded-network chemical conductance.

No NaP/Ca replacement. No graded presynaptic voltage/readout carrier. Original
compiled model is copied; only the copy loses the selected graded output axons.
The electrical input is constant during each physical host interval. Actual
within-interval LIF events and reset-synapse integrals are evaluated exactly.
Host neurons retain their existing update method, including graded ion gates.
"""
import copy,json,math,hashlib
import numpy as np
from table2_glif import ActivationDomain,derive_neuron,derive_synapse
from graded_profile import assert_campaign_allowed

class EventConductanceInputPort:
    def __init__(self,compiled,names,input_names,*,input_max_nA=1.,fmax_hz=500.,range_mV=20.,theta0_mV=1.,delta=.01,cells=16,enabled=True,outgoing_off=False):
        assert_campaign_allowed('unit_and_individual_pathway_probes')
        from sparse_pairwise_eval import SparsePairwiseForward
        raw=getattr(compiled,'c',compiled);self.names=list(names);self.indices=np.array([names.index(n) for n in input_names]);self.enabled=enabled;self.outgoing_off=outgoing_off;self.input_max=float(input_max_nA);self.dt=float(raw.dt);self.N=int(cells);self.time=0.;self.events=[];self.peaks=[];self.edge_specs=[];self.ports=[]
        if self.input_max<=0 or self.N<1:raise ValueError('Positive declared current domain/population required')
        self.activation=ActivationDomain(fmax_hz,range_mV,theta0_mV);self.design=derive_neuron(self.activation,0.,1.,.01);self.U=np.tile(np.linspace(0,theta0_mV,self.N,endpoint=False),(len(self.indices),1));self.total_counts=np.zeros(len(self.indices),dtype=np.int64)
        modified=copy.deepcopy(raw)
        for p,i in enumerate(self.indices):
            if np.any(raw.g_max_non[i]) or raw.i_b[i]!=0 or raw.V_rest[i]!=0:raise ValueError('Only zero-bias/rest pure physical-input cells are supported')
            if getattr(raw,'gated',False) and np.any(raw.g_ion[:,i]):raise ValueError('Ionic input cells are outside pure GLIF design')
            ports=np.flatnonzero(raw.input_connectivity[i])
            if len(ports)!=1 or np.count_nonzero(raw.input_connectivity[:,ports[0]])!=1:raise ValueError('Input port must uniquely address this cell')
            self.ports.append((int(ports[0]),float(raw.input_connectivity[i,ports[0]])))
            full_voltage=self.input_max/raw.g_m[i]
            for post in np.flatnonzero(raw.g_max_non[:,i]):
                fraction=float(np.clip((full_voltage-raw.e_lo[post,i])/(raw.e_hi[post,i]-raw.e_lo[post,i]),0,1));mean=float(raw.g_max_non[post,i]*fraction)
                s=derive_synapse(self.activation,1.,0.,float(raw.del_e[post,i]),delta)
                peak=mean/(s.tau_s*fmax_hz*(1-s.actual_delta_at_fmax));self.peaks.append(peak);self.edge_specs.append((p,int(post),float(raw.del_e[post,i]),s.tau_s,mean))
            if enabled:modified.g_max_non[:,i]=0.
        self.inner=SparsePairwiseForward(modified);self.gstate=np.zeros((len(self.edge_specs),self.N));self.gs=np.zeros(len(raw.V));self.ge=np.zeros(len(raw.V));self.counts=np.zeros(len(self.indices),dtype=np.int64)
        self.policy=dict(input_names=input_names,input_max_nA=input_max_nA,fmax_physical_hz=fmax_hz,range_mV=range_mV,theta0_mV=theta0_mV,delta=delta,cells=cells,enabled=enabled,outgoing_off=outgoing_off,clock='Physical seconds only; no global neural clock/cadence multiplier',target='Match each copied graded edge mean conductance at the declared physical-input maximum; Eq28 event reset/decay, not Eq29 spiking receiver correction',scope='New graded FSA profile plus pure-input spiking port; graded RG/PF/MN retained. Not the old standing winner or a fully Table2-verified controller.')
        self.parameter_hash=self._parameter_hash()
    def __getattr__(self,name):return getattr(self.inner,name)
    def _parameter_hash(self):
        h=hashlib.sha256()
        for key in ('c_m','g_m','V_rest','i_b','time_factor_membrane','g_max_non','del_e','e_lo','e_hi','input_connectivity'):h.update(np.ascontiguousarray(getattr(self.inner.c,key)).tobytes())
        return h.hexdigest()
    def verify_consumption(self):
        if self._parameter_hash()!=self.parameter_hash:raise ValueError('Underlying weights/parameters changed after event-port design; rebuild from the final intended profile')
        return dict(parameters_consumed=True,all_replaced_graded_axons_removed=bool(not self.enabled or not np.any(self.inner.c.g_max_non[:,self.indices])),physical_dt_s=self.dt)
    def _events(self,current):
        all_events=[];tau=self.design.tau_mem_s;theta=self.activation.theta0_mV
        for p,I in enumerate(current):
            inf=I+self.design.bias_nA;cell_events=[]
            for j in range(self.N):
                remaining=self.dt;elapsed=0.;v=self.U[p,j];times=[]
                while remaining>1e-15:
                    crossing=tau*math.log((inf-v)/(inf-theta)) if inf>theta and v<theta else (0. if v>=theta else math.inf)
                    if crossing>remaining:v=inf+(v-inf)*math.exp(-remaining/tau);break
                    elapsed+=crossing;remaining-=crossing;times.append(elapsed);v=0.
                    if len(times)>100:raise RuntimeError('Unexpected event saturation; no silent spike clipping')
                self.U[p,j]=v;cell_events.append(times)
                for t in times:self.events.append((self.time+t,p,j))
            all_events.append(cell_events)
        return all_events
    def forward(self,x=None):
        if not self.enabled:return self.inner.forward(x)
        if np.any(self.inner.c.g_max_non[:,self.indices]):raise ValueError('A replaced graded axon was reintroduced; no hidden graded carrier allowed')
        physical=np.array([np.asarray(x)[port]*weight for port,weight in self.ports]);fraction=physical/self.input_max
        if np.any(fraction<0) or np.any(fraction>1+1e-12):raise ValueError('Actual physical-input current exceeds declared design range; no silent clamp')
        events=self._events(fraction*self.activation.range_mV);self.gs[:]=0.;self.ge[:]=0.;self.counts[:]=0
        for p,e in enumerate(events):self.counts[p]=sum(len(t) for t in e)
        for k,(p,post,rev,tau,mean) in enumerate(self.edge_specs):
            total=0.
            for j in range(self.N):
                prev=0.;g=self.gstate[k,j];integral=0.
                for t in events[p][j]:integral+=g*tau*(-math.expm1(-(t-prev)/tau));g=self.peaks[k];prev=t
                integral+=g*tau*(-math.expm1(-(self.dt-prev)/tau));self.gstate[k,j]=g*math.exp(-(self.dt-prev)/tau);total+=integral/(self.N*self.dt)
            self.gs[post]+=total;self.ge[post]+=total*rev
        self.inner.forward(x)
        if not self.outgoing_off:self.inner.c.V+=self.inner.c.time_factor_membrane*(self.ge-self.gs*self.inner.last_voltage)
        if self.inner.output_onehot:np.take(self.inner.c.V,self.inner.oq,out=self.inner.output_values);self.inner.c.outputs=self.inner.output_values
        else:self.inner.c.outputs=self.inner.c.output_voltage_connectivity@self.inner.c.V
        self.time+=self.dt;self.total_counts+=self.counts;return self.inner.c.outputs
    def export_state(self):
        state=dict(policy=np.array(json.dumps(self.policy,sort_keys=True)),time=np.array(self.time),U=self.U.copy(),gstate=self.gstate.copy(),gs=self.gs.copy(),ge=self.ge.copy(),counts=self.counts.copy(),total_counts=self.total_counts.copy(),events=np.array(self.events).reshape(-1,3))
        for key in ('V','V_last','b_gate','c_gate','b_gate_last','c_gate_last','outputs'):
            if hasattr(self.inner.c,key):state['raw_'+key]=np.asarray(getattr(self.inner.c,key)).copy()
        for key in ('last_voltage','input_values','i_app','output_values','sums','acc','values'):state['sparse_'+key]=getattr(self.inner,key).copy()
        return state
    def import_state(self,state):
        if str(state['policy'])!=json.dumps(self.policy,sort_keys=True):raise ValueError('Port policy mismatch')
        self.time=float(state['time']);self.events=[tuple(v) for v in state['events']]
        for key in ('U','gstate','gs','ge','counts','total_counts'):np.copyto(getattr(self,key),state[key])
        for key in ('V','V_last','b_gate','c_gate','b_gate_last','c_gate_last','outputs'):
            if 'raw_'+key in state:setattr(self.inner.c,key,state['raw_'+key].copy())
        for key in ('last_voltage','input_values','i_app','output_values','sums','acc','values'):np.copyto(getattr(self.inner,key),state['sparse_'+key])
        self.inner.c.V_last=self.inner.last_voltage
        if self.inner.output_onehot:self.inner.c.outputs=self.inner.output_values
