"""Hunt step A: train the undriven RG period, then verify in SNS-Toolbox.

The preallocated Numba kernel implements the compiled four-cell fixed-tau
SNS equations without fastmath. Its numerical error is measured explicitly;
the winner is replayed with the original SNS backend before acceptance.
"""
import os,sys,json,time,hashlib
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parent;S=H.parent;sys.path.insert(0,str(S))
import numpy as np
from numba import njit
from scipy.signal import find_peaks
from scipy.optimize import brentq
from sns_toolbox.networks import Network
from sns_toolbox.neurons import NonSpikingNeuronWithPersistentSodiumChannel
from build_network import _neu,_syn,SNS_NumpyFixedTau
from params import NAP
DT=.002

def build(tau_h,ginh=4.,nap_values=None,afferent_gains=None):
    nv={**NAP,**(nap_values or {})}
    n=Network(name='Hunt intrinsic RG')
    nap=NonSpikingNeuronWithPersistentSodiumChannel(membrane_capacitance=.05,
        membrane_conductance=1.,resting_potential=0.,bias=0.,
        g_ion=np.array([nv['g_ion']]),e_ion=np.array([nv['e_ion']]),
        k_m=np.array([nv['k_m']]),slope_m=np.array([nv['slope_m']]),e_m=np.array([nv['e_m']]),
        k_h=np.array([nv['k_h']]),slope_h=np.array([nv['slope_h']]),e_h=np.array([nv['e_h']]),
        tau_max_h=np.array([tau_h]))
    for name in ('E','F'):n.add_neuron(nap,name=name);n.add_input(name)
    for name in ('InE','InF'):n.add_neuron(_neu(.05),name=name)
    for a,b,exc in [('E','InE',True),('InE','F',False),('F','InF',True),('InF','E',False)]:
        n.add_connection(_syn(ginh,exc),a,b)
    if afferent_gains is not None:
        for i,target in enumerate(('F','E','F','E','E','E')):
            name='AFF'+str(i);n.add_neuron(_neu(.02),name=name);n.add_input(name)
            n.add_connection(_syn(float(afferent_gains[i]),True),name,target)
    for population in n.populations:n.add_output(population['name'])
    c=n.compile(dt=DT,backend='numpy');c.__class__=SNS_NumpyFixedTau
    return c

def arrays(c):
    assert np.all(c.pow_c==0)
    return tuple(np.ascontiguousarray(a,dtype=np.float64) for a in (
        c.g_max_non,c.del_e,c.e_lo,c.e_hi-c.e_lo,c.g_m,c.V_rest,c.i_b,
        c.time_factor_membrane,c.g_ion[0],c.e_ion[0],c.k_a[0],c.slope_a[0],c.e_a[0],
        c.k_b[0],c.slope_b[0],c.e_b[0],c.tau_max_b[0],c.pow_a[0],c.pow_b[0]))

@njit(cache=True)
def simulate(p,v,h,u):
    gm,rev,lo,span,leak,rest,bias,factor,ion,ei,ka,sa,ea,kb,sb,eb,th,pa,pb=p
    count=len(v);trace=np.empty((len(u),count));hh=np.empty_like(trace)
    old=np.empty(count)
    for t in range(len(u)):
        old[:]=v
        for i in range(count):
            gg=0.;ge=0.
            for j in range(count):
                if gm[i,j]==0.:continue
                saturation=min(1.,max(0.,(old[j]-lo[i,j])/span[i,j]))
                g=gm[i,j]*saturation;gg+=g;ge+=g*rev[i,j]
            a=1./(1.+ka[i]*np.exp(sa[i]*(ea[i]-old[i])))
            b=1./(1.+kb[i]*np.exp(sb[i]*(eb[i]-old[i])))
            h[i]+=.002*(b-h[i])/th[i]
            current=ion[i]*(a**pa[i])*(h[i]**pb[i])*(ei[i]-old[i])
            v[i]=old[i]+factor[i]*(-leak[i]*(old[i]-rest[i])+bias[i]+ge-old[i]*gg+u[t,i]+current)
        trace[t]=v;hh[t]=h
    return trace,hh

def measure(trace):
    t=np.arange(len(trace))*DT;keep=t>=10.;signal=trace[:,0]-trace[:,1]
    peaks=find_peaks(signal,prominence=2.,distance=int(.2/DT))[0];peaks=peaks[t[peaks]>=10]
    periods=np.diff(t[peaks]);T=float(np.median(periods)) if len(periods) else 0.
    return dict(period_s=T,cycles=len(peaks),period_cv=float(np.std(periods)/np.mean(periods)) if len(periods) else 9.,
        extensor_duty=float(np.mean(signal[keep]>0)),finite=bool(np.all(np.isfinite(trace))))

def main():
    out=H/'runs/step_A';out.mkdir(parents=True,exist_ok=True)
    source=np.load(S/'bsolve_out.npz',allow_pickle=True);target_T=float(np.median(np.diff(source['hs_r'])))*2
    u=np.zeros((int(40/DT),4));u[:50,0]=2.
    trials=[]
    candidates=[]
    for em in (-2.,-1.5,-1.,-.5,0.,.5,1.,1.5,2.):
        for eh in (2.,3.5,5.):
            for inh in (2.,4.,6.):
                nv=dict(NAP,e_m=em,e_h=eh)
                c=build(.6,inh,nv);v,h=simulate(arrays(c),c.V.copy(),c.b_gate[0].copy(),u)
                metric=measure(v)
                if metric['cycles']>=4 and metric['period_cv']<.05 and .2<metric['extensor_duty']<.8:
                    candidates.append((abs(metric['period_s']-target_T),inh,nv,metric))
    (out/'regime_scan.json').write_text(json.dumps([dict(error=a,ginh=b,nap=c,metrics=d) for a,b,c,d in candidates],indent=2))
    if not candidates:raise RuntimeError('No undriven oscillatory regime found; do not proceed to sensory fitting.')
    _,ginh,nv,_=min(candidates,key=lambda row:row[0])
    def evaluate(tau):
        c=build(tau,ginh,nv);v,h=simulate(arrays(c),c.V.copy(),c.b_gate[0].copy(),u)
        result=dict(tau_h=tau,**measure(v));trials.append(result)
        print(json.dumps(result),flush=True)
        return result['period_s']-target_T
    tau=float(brentq(evaluate,.25,1.5,xtol=1e-4));nv['tau_max_h']=tau;c=build(tau,ginh,nv)
    fast_v,fast_h=simulate(arrays(c),c.V.copy(),c.b_gate[0].copy(),u)
    native=np.empty_like(fast_v);native_h=np.empty_like(fast_h)
    for i in range(len(u)):
        c.forward(u[i,:2]);native[i]=c.V;native_h[i]=c.b_gate[0]
    voltage_error=float(np.max(abs(native-fast_v)));h_error=float(np.max(abs(native_h-fast_h)))
    fast=measure(fast_v);reference=measure(native)
    accepted=reference['finite'] and abs(reference['period_s']/target_T-1)<.02 and reference['period_cv']<.02 and .2<reference['extensor_duty']<.8
    # The implementation may differ by floating-point roundoff; never call it bit-exact.
    assert voltage_error<1e-7 and h_error<1e-8,(voltage_error,h_error)
    result=dict(step='A_intrinsic',target_stepping_period_s=target_T/2,target_undriven_period_s=target_T,
        tau_h=tau,mutual_inhibition=ginh,nap=nv,dt=DT,tonic_drive=0.,
        initialization='2 nA to RG-E for first 100 ms; no input during the scored 10-40 s window',
        reference=reference,fast=fast,accepted=accepted,voltage_max_abs_error=voltage_error,h_gate_max_abs_error=h_error,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),trials=trials)
    (out/'winner.json').write_text(json.dumps(result,indent=2));np.savez_compressed(out/'validation.npz',t=np.arange(len(u))*DT,fast_v=fast_v,native_v=native)
    print('STEP_A',json.dumps({k:v for k,v in result.items() if k not in ('trials','nap')}),flush=True)
if __name__=='__main__':main()
