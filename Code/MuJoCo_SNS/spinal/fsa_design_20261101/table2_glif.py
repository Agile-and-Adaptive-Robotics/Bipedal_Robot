"""Analytical FSA design, Szczecinski et al.2020 Table2; no automatic probes.

Public units: seconds, Hz, relative-rest mV, nA, nF, and microSiemens.
nF/uS is milliseconds, so Cm_nF=1000*tau_s*G_uS.
This pure GLIF derivation does NOT cover NaP/Ca/other ionic bursters.
Existing controller implementations and fitted weights are never modified.
"""
from dataclasses import dataclass,asdict
import math
import numpy as np
from numba import njit

SOURCE=dict(doi='10.3389/fnbot.2020.577804',table_pdf_page=11,
            dynamics_pdf_page=5,mapping_pdf_page=8,
            eq23='tau_mem=R/(Fmax*theta_star)',
            eq26='tau_theta=graded_tau_mem*(1-m/2)')

def finite(name,value,positive=False):
    if not math.isfinite(value) or (positive and value<=0):
        raise ValueError(name+' must be finite'+(' and positive' if positive else ''))
    return float(value)

@dataclass(frozen=True)
class ActivationDomain:
    """Step1 choices must be supplied, not inferred from a successful trial."""
    fmax_hz:float
    range_mV:float
    theta0_mV:float
    def __post_init__(self):
        for name in ('fmax_hz','range_mV','theta0_mV'):finite(name,getattr(self,name),True)

@dataclass(frozen=True)
class NeuronDesign:
    activation:ActivationDomain
    m:float
    leak_uS:float
    graded_reference_tau_s:float
    tau_theta_s:float|None
    theta_star_approx_mV:float
    bias_nA:float
    tau_mem_s:float
    capacitance_nF:float
    compatibility:str
    approximation:str

def derive_neuron(domain,m,leak_uS,graded_reference_tau_s,*,ionic_conductance_uS=0.):
    """Literal Table2 steps2–5, with explicit high-rate approximation scope.

    The reference graded tau is an operating-point quantity. In a coupled
    conductance network Cm/(Gleak+Gsyn) varies with input; one reference tau
    does not establish uniform transient equivalence across all drives.
    """
    finite('m',m);finite('leak_uS',leak_uS,True)
    finite('graded_reference_tau_s',graded_reference_tau_s,True)
    if ionic_conductance_uS!=0:
        raise NotImplementedError('NaP/Ca/ionic neurons are outside pure-GLIF Table2; preserve graded RG or design a separately verified ionic-spiking model')
    if m>=2:raise ValueError('m<2 is required for positive theta_star and Table2 time constants')
    b=1/(1-m/2);theta=domain.theta0_mV*b
    tau=domain.range_mV/(domain.fmax_hz*theta)
    return NeuronDesign(domain,float(m),float(leak_uS),float(graded_reference_tau_s),
        None if m==0 else graded_reference_tau_s/b,theta,
        leak_uS*domain.theta0_mV/(2-m),tau,1000*leak_uS*tau,
        'pure_GLIF_constant_threshold' if m==0 else 'pure_GLIF_adaptive_threshold',
        'Eq16 theta*=theta0/(1-m/2) is a high-rate approximation; Eq23 endpoint and Eq26 transients require independent rate/step validation, not equality by assertion.')

@dataclass(frozen=True)
class SynapseDesign:
    fmax_hz:float
    tau_s:float
    target_delta:float
    actual_delta_at_fmax:float
    reversal_mV:float
    functional_gain:float
    target_mean_uS:float
    peak_table29_uS:float
    peak_finite_delta28_uS:float
    selected_peak_uS:float
    calibration_mode:str
    incremental_inhibition:bool

def derive_synapse(domain,post_leak_uS,gain,reversal_mV,delta,*,tau_s=None,
                   calibration_mode='table29',incremental_inhibition=False):
    """Steps6–7. Event G is RESET to peak, never incremented by peak.

    Eq29 prints normalized Gmem=1. Restoring a physical nonunit leak factor
    in Gbar follows the static equilibrium equation and is labelled here.
    Eq28 finite-delta correction and Eq29 asymptote are both exposed; the
    selected rule is explicit. Inhibition is an incremental signed pathway
    about a declared operating point, not a negative firing-rate target.
    """
    finite('post_leak_uS',post_leak_uS,True);finite('gain',gain);finite('reversal_mV',reversal_mV)
    if not 0<delta<1:raise ValueError('delta must lie strictly between0 and1')
    upper=-1/(domain.fmax_hz*math.log(delta));tau=upper if tau_s is None else finite('tau_s',tau_s,True)
    if tau>upper*(1+1e-12):raise ValueError('tau_s violates Eq27 maximum nonlinearity at Fmax')
    if gain<0 and not incremental_inhibition:raise ValueError('Negative gain requires explicit signed/baseline encoding')
    denominator=reversal_mV-gain*domain.range_mV
    if denominator==0:raise ValueError('Unbounded conductance at Es=k*R')
    mean=post_leak_uS*gain*domain.range_mV/denominator
    if mean<0:raise ValueError('Requested gain/reversal requires negative conductance')
    actual=math.exp(-1/(domain.fmax_hz*tau));asym=mean/(tau*domain.fmax_hz);exact=asym/(1-actual)
    if calibration_mode not in ('table29','finite_delta28'):raise ValueError('Choose table29 or finite_delta28 explicitly')
    return SynapseDesign(domain.fmax_hz,tau,float(delta),actual,float(reversal_mV),float(gain),
        mean,asym,exact,asym if calibration_mode=='table29' else exact,calibration_mode,bool(incremental_inhibition))

def constant_threshold_rate_hz(design,applied_nA,*,conductance_uS=0.,conductance_times_reversal_nA=0.):
    """Eq13 for m=0 only; not a simulation and not valid for ionic bursters."""
    if design.m!=0:raise NotImplementedError('Adaptive theta requires Eq14 or separately verified event dynamics')
    g=design.leak_uS+conductance_uS
    if g<=0:raise ValueError('Nonpositive total conductance')
    target=(applied_nA+design.bias_nA+conductance_times_reversal_nA)/g
    if target<=design.activation.theta0_mV:return 0.
    effective_tau=design.capacitance_nF/(1000*g)
    return -1/(effective_tau*math.log1p(-design.activation.theta0_mV/target))

def average_reset_synapse_uS(peak_uS,tau_s,rate_hz):
    """Eq18; a reset-to-Gmax synapse, not additive spike filtering."""
    if rate_hz<=0:return 0.
    return peak_uS*tau_s*rate_hz*(-math.expm1(-1/(rate_hz*tau_s)))

def analytic_passive_voltage_bound(leak_uS,rest_mV,bias_nA,input_current_nA_bounds,edges):
    """Conservative exact envelope, assuming independently bounded chemical G.

    Relative-rest equilibrium is a convex mixture of Ibase/Gleak and Erev.
    Missing physical input bounds yield an explicitly unresolved domain.
    This is NOT an invariant dynamic/ionic bound or a task design range.
    """
    if input_current_nA_bounds is None:return dict(resolved=False,reason='Physical/descending current domain not specified; observed traces are not design limits')
    lo,hi=input_current_nA_bounds
    if lo>hi:raise ValueError('Reversed input-current bounds')
    finite('leak_uS',leak_uS,True)
    values=[(bias_nA+lo)/leak_uS,(bias_nA+hi)/leak_uS]
    for peak,reversal in edges:
        if peak<0:raise ValueError('Negative conductance')
        if peak>0:values.append(reversal-rest_mV)
    return dict(resolved=True,relative_rest_equilibrium_mV=[min(values),max(values)],
                interpretation='Conservative equilibrium envelope; no claim all extrema are reachable; not a selected encoding domain or ionic bound')

class GLIFPopulationRuntime:
    """Design-probe runtime only. No physical/controller adapter is installed.

    Same voltage/threshold equations as paper Eqs8–10; no added refractory,
    SFA current, NaP channel, graded voltage carrier, or arbitrary rate decoder.
    Crossing is detected at a declared finite-step boundary, so convergence
    tests are mandatory before reuse. Initial phases are deterministic.
    Each pre-cell has its own RESET-and-decay synaptic state, not summed events.
    Population outputs are the average of actual per-cell conductances.
    """
    def __init__(self,neurons,edges,population_size,step_s,*,allow_design_probe=False):
        if not allow_design_probe:raise RuntimeError('Unverified design runtime: root authorization for pathway probes is required; global training remains blocked')
        self.neurons=list(neurons);self.edges=list(edges);self.N=int(population_size);self.dt=finite('step_s',step_s,True)
        if self.N<1:raise ValueError('Positive population size required')
        self.U=np.array([[d.activation.theta0_mV*j/self.N for j in range(self.N)] for d in neurons]);self.theta=np.array([[d.activation.theta0_mV]*self.N for d in neurons])
        self.gs=np.zeros((len(edges),self.N));self.events=np.zeros_like(self.U,dtype=np.bool_);self.total_events=np.zeros(len(neurons),dtype=np.int64);self.time_s=0.
        for pre,post,s in edges:
            if not(0<=pre<len(neurons) and 0<=post<len(neurons)):raise ValueError('Invalid edge index')
            if self.dt>s.tau_s/20:raise ValueError('Declare finer step: at least20 samples per synaptic decay; this alone is not a convergence proof')
        if self.dt>min(d.tau_mem_s for d in neurons)*1e-4:raise ValueError('Step exceeds the paper simulation policy min(tau_mem)*1e-4; no unvalidated coarse clock substitution')
        self.out_g=np.zeros(len(neurons));self.out_ge=np.zeros(len(neurons))
        self.numeric_neurons=np.array([[d.leak_uS,d.capacitance_nF,d.bias_nA,d.activation.theta0_mV,d.m,d.tau_theta_s or 1.] for d in neurons])
        self.numeric_edges=np.array([[pre,post,s.selected_peak_uS,s.tau_s,s.reversal_mV] for pre,post,s in edges]).reshape(-1,5)
    def conductance_ports(self):
        self.out_g[:]=0;self.out_ge[:]=0
        for k,(pre,post,s) in enumerate(self.edges):
            g=float(self.gs[k].mean());self.out_g[post]+=g;self.out_ge[post]+=g*s.reversal_mV
        return self.out_g,self.out_ge
    def step_reference(self,applied_nA):
        u=np.asarray(applied_nA,dtype=float)
        if u.shape!=(len(self.neurons),) or not np.isfinite(u).all():raise ValueError('Finite current per population required')
        for k,(_,_,s) in enumerate(self.edges):self.gs[k]*=math.exp(-self.dt/s.tau_s)
        g,ge=self.conductance_ports();self.events[:]=False
        for i,d in enumerate(self.neurons):
            a=1000*(d.leak_uS+g[i])/d.capacitance_nF;target=(u[i]+d.bias_nA+ge[i])/(d.leak_uS+g[i]);old=self.U[i].copy();ea=math.exp(-a*self.dt)
            self.U[i]=target+(old-target)*ea
            if d.m==0:self.theta[i]=d.activation.theta0_mV
            else:
                b=1/d.tau_theta_s;eb=math.exp(-b*self.dt);limit=d.activation.theta0_mV+d.m*target
                cross=d.m*b*(old-target)*(self.dt*ea if abs(a-b)<1e-12*max(a,b) else (ea-eb)/(b-a))
                self.theta[i]=limit+(self.theta[i]-limit)*eb+cross
            self.events[i]=self.U[i]>=self.theta[i];self.U[i,self.events[i]]=0.
        for k,(pre,_,s) in enumerate(self.edges):self.gs[k,self.events[pre]]=s.selected_peak_uS
        self.total_events+=self.events.sum(axis=1);self.time_s+=self.dt
        return self.events
    def step(self,applied_nA):
        u=np.asarray(applied_nA,dtype=float)
        if u.shape!=(len(self.neurons),) or not np.isfinite(u).all():raise ValueError('Finite current per population required')
        glif_tick(self.U,self.theta,self.gs,self.events,self.total_events,self.numeric_neurons,self.numeric_edges,u,self.out_g,self.out_ge,self.dt)
        self.time_s+=self.dt;return self.events
    def graded_target_currents_nA(self,graded_voltage_mV):
        """Actual spike-conductance port: sum Gs*(E−V), no proxy voltage."""
        g,ge=self.conductance_ports();return ge-g*np.asarray(graded_voltage_mV)
    def export_state(self):
        return dict(U=self.U.copy(),theta=self.theta.copy(),gs=self.gs.copy(),events=self.events.copy(),total_events=self.total_events.copy(),time_s=np.array(self.time_s),dt_s=np.array(self.dt))
    def import_state(self,state):
        if float(state['dt_s'])!=self.dt:raise ValueError('State clock mismatch')
        for key in ('U','theta','gs','events','total_events'):
            if np.shape(state[key])!=np.shape(getattr(self,key)):raise ValueError('State shape mismatch')
            np.copyto(getattr(self,key),state[key])
        self.time_s=float(state['time_s'])

@njit(cache=True)
def glif_tick(U,theta,gs,events,totals,parameters,edges,applied,out_g,out_ge,dt):
    """Preallocated frozen-conductance split step; no fastmath or RNG.

    Its reference/clock convergence checks are deliberately NOT claimed yet.
    Event reset occurs at the declared step end; at most one event/cell/step.
    """
    out_g[:]=0.;out_ge[:]=0.;events[:]=False
    N=U.shape[1]
    for e in range(len(edges)):
        pre,post,peak,tau,rev=edges[e];post=int(post);decay=math.exp(-dt/tau);mean=0.
        for j in range(N):gs[e,j]*=decay;mean+=gs[e,j]/N
        out_g[post]+=mean;out_ge[post]+=mean*rev
    for i in range(len(U)):
        leak,cap,bias,theta0,m,tau_theta=parameters[i];g=leak+out_g[i];a=1000*g/cap;target=(applied[i]+bias+out_ge[i])/g;ea=math.exp(-a*dt)
        for j in range(N):
            old=U[i,j];U[i,j]=target+(old-target)*ea
            if m==0:theta[i,j]=theta0
            else:
                b=1/tau_theta;eb=math.exp(-b*dt);limit=theta0+m*target
                response=dt*ea if abs(a-b)<1e-12*max(a,b) else (ea-eb)/(b-a)
                theta[i,j]=limit+(theta[i,j]-limit)*eb+m*b*(old-target)*response
            if U[i,j]>=theta[i,j]:events[i,j]=True;U[i,j]=0.;totals[i]+=1
    for e in range(len(edges)):
        pre=int(edges[e,0]);peak=edges[e,2]
        for j in range(N):
            if events[pre,j]:gs[e,j]=peak
