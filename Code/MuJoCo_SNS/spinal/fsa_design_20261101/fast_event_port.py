"""Preallocated exact-event sensory port, physical-domain design v2.

Design-only: no curriculum/body authorization. The v1 reference is unchanged.
Default telemetry is counts and a bounded last event block, never a growing log.
"""
import math,json,hashlib,copy
import numpy as np
from numba import njit
if __package__:
    from .event_conductance_port import EventConductanceInputPort,default_port_profile
    from .graded_profile import assert_campaign_allowed,population_names,field_digest
else:
    from event_conductance_port import EventConductanceInputPort,default_port_profile
    from graded_profile import assert_campaign_allowed,population_names,field_digest

@njit(cache=True)
def event_conductance_tick(x,ports,input_limits,R,theta,cell_bias,tau_mem,dt,time,U,counts,total_counts,cell_counts,event_times,event_block,specs,peaks,gstate,gs,ge):
    """Actual threshold crossings and reset-G integrals; no fastmath or allocation."""
    for p in range(len(ports)):
        physical=x[int(ports[p,0])]*ports[p,1]
        if not math.isfinite(physical) or physical<0 or physical>input_limits[p]*(1+1e-12):raise ValueError('Physical input outside declared engineering domain')
    block=0
    for p in range(len(ports)):
        physical=x[int(ports[p,0])]*ports[p,1]
        applied=physical/input_limits[p]*R;counts[p]=0
        for j in range(U.shape[1]):
            inf=applied+cell_bias[p,j]
            remaining=dt;elapsed=0.;v=U[p,j];n=0
            while remaining>1e-15:
                crossing=tau_mem*math.log((inf-v)/(inf-theta)) if inf>theta and v<theta else (0. if v>=theta else math.inf)
                if crossing>remaining:v=inf+(v-inf)*math.exp(-remaining/tau_mem);break
                if n>=event_times.shape[2]:raise ValueError('Bounded event capacity exceeded; no silent clipping')
                elapsed+=crossing;remaining-=crossing;event_times[p,j,n]=elapsed;n+=1;v=0.
                event_block[block,0]=time+elapsed;event_block[block,1]=p;event_block[block,2]=j;block+=1
            U[p,j]=v;cell_counts[p,j]=n;counts[p]+=n
        total_counts[p]+=counts[p]
    for i in range(len(gs)):gs[i]=0.;ge[i]=0.
    for k in range(len(specs)):
        p=int(specs[k,0]);post=int(specs[k,1]);rev=specs[k,2];tau=specs[k,3];total=0.
        for j in range(U.shape[1]):
            previous=0.;g=gstate[k,j];integral=0.
            for z in range(cell_counts[p,j]):
                t=event_times[p,j,z];integral+=g*tau*(-math.expm1(-(t-previous)/tau));g=peaks[k];previous=t
            integral+=g*tau*(-math.expm1(-(dt-previous)/tau));gstate[k,j]=g*math.exp(-(dt-previous)/tau);total+=integral/(U.shape[1]*dt)
        gs[post]+=total;ge[post]+=total*rev
    return block

@njit(cache=True)
def add_event_current(V,factor,gs,ge,last):
    for i in range(len(V)):V[i]+=factor[i]*(ge[i]-gs[i]*last[i])

class FastEventConductanceInputPort(EventConductanceInputPort):
    def __init__(self,*args,capture_events=False,backend='numba',population_bias='homogeneous',input_limits_nA=None,**kwargs):
        if backend not in ('numba','reference') or not isinstance(capture_events,bool):raise ValueError('Explicit backend and capture switch required')
        if population_bias not in ('homogeneous','midpoint_uniform'):raise ValueError('Unknown analytically declared population bias design')
        source=args[0] if args else kwargs.get('compiled')
        if getattr(source,'inner',None) is not None:raise ValueError('Do not silently unwrap an existing neural encoder/central wrapper; explicit combined conductance/state coupling must be designed and verified separately')
        super().__init__(*args,**kwargs)
        # Deterministic initial scratch state; none of these buffers is a
        # presynaptic signal. Their first neural update overwrites them.
        self.inner.input_values.fill(0.);np.copyto(self.inner.last_voltage,self.inner.c.V_last)
        if hasattr(self.inner.c,'outputs'):np.copyto(self.inner.output_values,self.inner.c.outputs)
        else:self.inner.output_values.fill(0.);self.inner.c.outputs=self.inner.output_values
        for key in ('syn_ge_sum','syn_g_sum'):
            if hasattr(self.inner,key):getattr(self.inner,key).fill(0.)
        self.backend=backend;self.capture_events=capture_events;self.ports=np.asarray(self.ports,dtype=float);self.edge_specs=np.asarray(self.edge_specs,dtype=float).reshape(-1,5);self.peaks=np.asarray(self.peaks,dtype=float)
        self.input_limits=np.full(len(self.indices),self.input_max)
        if input_limits_nA is not None:
            if set(input_limits_nA)!=set(self.policy['input_names']):raise ValueError('Every selected input requires an explicit per-type limit')
            self.input_limits[:]=[input_limits_nA[n] for n in self.policy['input_names']]
            if not np.isfinite(self.input_limits).all() or np.any(self.input_limits<=0) or np.any(self.input_limits>self.input_max):raise ValueError('Per-type limits must lie inside the source linear window domain')
            for k,e in enumerate(self.edge_specs):
                scale=self.input_limits[int(e[0])]/self.input_max;self.edge_specs[k,4]*=scale;self.peaks[k]*=scale
        self.cell_bias=np.full(self.U.shape,self.design.bias_nA)
        if population_bias=='midpoint_uniform':
            biases=self.activation.theta0_mV*(np.arange(self.N)+.5)/self.N
            self.cell_bias[:]=biases;self.U[:]=self.cell_bias # Exact quiet equilibrium; no precharge/spikes.
            full=self.activation.range_mV+biases;rates=1/(self.design.tau_mem_s*np.log(full/(full-self.activation.theta0_mV)))
            rate_bound=float(rates.max());tau=-1/(rate_bound*math.log(self.policy['delta']));reset_mean=float(np.mean(tau*rates*(-np.expm1(-1/(tau*rates)))))
            self.edge_specs[:,3]=tau;self.peaks[:]=self.edge_specs[:,4]/reset_mean
            self.policy.update(population_bias_scope='Engineered mean-bias population extension, not Table2 per-cell bias. Mean bias equalsG*theta/2; every bias<threshold. Exact full-population reset-conductance normalization.',individual_bias_nA=biases.tolist(),full_input_individual_rate_Hz=rates.tolist(),full_input_mean_rate_Hz=float(rates.mean()),synapse_delta_max_rate_Hz=rate_bound)
        capacity=math.ceil(self.activation.fmax_hz*self.dt)+2
        self.cell_counts=np.zeros_like(self.U,dtype=np.int64);self.event_times=np.zeros((*self.U.shape,capacity));self.event_block=np.zeros((self.U.size*capacity,3));self.event_block_count=0
        self.policy.update(implementation='preallocated_exact_event_v2',population_bias=population_bias,input_limits_nA=dict(zip(self.policy['input_names'],self.input_limits.tolist())),physical_strength_preservation='Each endpointGtarget scaled by typed_max/source_linear_window_max; all physicalI slopes remain the sourceG/window slope, no motor gain change.',capture_events=capture_events,event_capacity_per_cell=capacity,logging='Counts+last bounded block by default; full event history only when explicitly capture_events=True')
        self.parameter_hash=self._parameter_hash()
    def _parameter_hash(self):
        h=hashlib.sha256(super()._parameter_hash().encode())
        if 'cell_bias' in self.__dict__:h.update(self.cell_bias.tobytes())
        if 'input_limits' in self.__dict__:h.update(self.input_limits.tobytes())
        return h.hexdigest()
    def forward(self,x=None):
        if not self.enabled:return self.inner.forward(x)
        for i in self.indices:
            if np.any(self.inner.c.g_max_non[:,i]):raise ValueError('A replaced graded axon was reintroduced')
        kernel=event_conductance_tick if self.backend=='numba' else event_conductance_tick.py_func
        self.event_block_count=kernel(np.asarray(x),self.ports,self.input_limits,self.activation.range_mV,self.activation.theta0_mV,self.cell_bias,self.design.tau_mem_s,self.dt,self.time,self.U,self.counts,self.total_counts,self.cell_counts,self.event_times,self.event_block,self.edge_specs,self.peaks,self.gstate,self.gs,self.ge)
        if self.capture_events:self.events.extend(map(tuple,self.event_block[:self.event_block_count]))
        self.inner.forward(x)
        if not self.outgoing_off:
            apply=add_event_current if self.backend=='numba' else add_event_current.py_func
            apply(self.inner.c.V,self.inner.c.time_factor_membrane,self.gs,self.ge,self.inner.last_voltage)
        if self.inner.output_onehot:np.take(self.inner.c.V,self.inner.oq,out=self.inner.output_values);self.inner.c.outputs=self.inner.output_values
        else:self.inner.c.outputs=self.inner.c.output_voltage_connectivity@self.inner.c.V
        self.time+=self.dt;return self.inner.c.outputs
    def last_events(self):return self.event_block[:self.event_block_count]
    def export_state(self):
        s=super().export_state();s.update(fast_version=np.array(2),fast_event_block_count=np.array(self.event_block_count),fast_cell_counts=self.cell_counts.copy(),fast_event_times=self.event_times.copy(),fast_event_block=self.event_block.copy(),fast_cell_bias=self.cell_bias.copy(),fast_input_limits=self.input_limits.copy());return s
    def import_state(self,s):
        if int(s['fast_version'])!=2:raise ValueError('Wrong fast state version')
        if not np.array_equal(s['fast_cell_bias'],self.cell_bias):raise ValueError('Population bias design mismatch')
        if not np.array_equal(s['fast_input_limits'],self.input_limits):raise ValueError('Physical input domain mismatch')
        super().import_state(s);self.event_block_count=int(s['fast_event_block_count'])
        for a,k in ((self.cell_counts,'fast_cell_counts'),(self.event_times,'fast_event_times'),(self.event_block,'fast_event_block')):np.copyto(a,s[k])

PHYSICAL_SCHEMA='aarl_exact_event_physical_domain_v2'
POPULATION_SCHEMA='aarl_exact_event_physical_population_bias_v2'
TYPED_SCHEMA='aarl_exact_event_canonical_typed_population_v3'
TIGHT_SCHEMA='aarl_exact_event_canonical_typed_tight_tail_v1'

def physical_domain_profile(family,all_afferents=True):
    p=default_port_profile(family,all_afferents);p.update(schema=PHYSICAL_SCHEMA,input_max_nA=5.,theta0_mV=.1)
    p['binding']=dict(kind='engineering_admissible_canonical_window',units='nA',range=[0.,5.],not_universal_maximum=True,caller_aff_scale=1.,gamma_coact=0.,vest_prop=0.,canonical_maximum_Ia=2.88,canonical_maximum_II=2.92,canonical_maximum_Ib=1.125,requested_normal_Ib=.225,Ib_zero_force=0.)
    p['capture_events']=False
    return p

def physical_population_profile(family,all_afferents=True):
    p=physical_domain_profile(family,all_afferents);p.update(schema=POPULATION_SCHEMA,population_bias='midpoint_uniform');p['binding']['carrier_rate_interpretation']='500Hz population-mean design target; exact heterogeneous per-cell full-input rates and maximum are separately recorded, and synapticdelta uses the actual maximum';return p

def canonical_typed_profile(family,all_afferents=True):
    p=physical_population_profile(family,all_afferents);p.update(schema=TYPED_SCHEMA,cells=128,input_limits_nA={n:dict(Ia=2.88,II=2.92,Ib=1.125)[n.split('_')[0]] for n in p['input_names']});p['binding'].update(kind='canonical_runner_context_per_type',range='Per-input limits; source linear window remains0..5mV',not_universal_maximum=True,count_resolution_design='N128,F500,W20ms: one-event normalIb modulation quantization1/(NFW)/(0.1/1.125)=0.879%; this is not a receiver error bound. Receiver5%RMS criterion is independently tested.');return p

def canonical_tight_tail_profile(family,all_afferents=True):
    p=canonical_typed_profile(family,all_afferents);p.update(schema=TIGHT_SCHEMA,delta=.001);p['binding']['synaptic_tail_design']='Declare0.1% residual at the fastest actual per-cell firing period. Finite-tail gain budgetdelta/(1-delta)=0.1001%, down from1.0101%; tau and resetpeak are rederived exactly, not multiplied by a fitted gain.';return p

def install_physical_port(spinal,profile,*,family):
    assert_campaign_allowed('optional_profile_implementation');expected=canonical_tight_tail_profile(family) if profile.get('schema')==TIGHT_SCHEMA else (canonical_typed_profile(family) if profile.get('schema')==TYPED_SCHEMA else (physical_population_profile(family) if profile.get('schema')==POPULATION_SCHEMA else physical_domain_profile(family)))
    if set(profile)!=set(expected) or profile['schema']!=expected['schema'] or profile['family']!=family:raise ValueError('Missing/unused physical-port schema fields')
    if profile['binding']!=expected['binding']:raise ValueError('Unknown physical input binding; derive a new declared profile')
    if profile['schema'] in (TYPED_SCHEMA,TIGHT_SCHEMA) and (profile['cells']!=128 or profile['input_limits_nA']!=expected['input_limits_nA']):raise ValueError('Typed profile requires the declared canonical bounds andN128 design; register alternatives separately')
    names=population_names(spinal);old=spinal.compiled;args={k:v for k,v in profile.items() if k not in ('schema','family','binding')};port=FastEventConductanceInputPort(old,names,**args)
    spinal.compiled=port;spinal.fsa_spiking_port_profile=copy.deepcopy(profile);spinal.fsa_spiking_port_ledger=dict(schema=profile['schema'],family=family,compiled_before_sha256=field_digest(getattr(old,'c',old)),parameter_sha256=port.parameter_hash,counts_only_default=not port.capture_events,input_binding=profile['binding'],training_allowed=False,consumption=port.verify_consumption());return port
