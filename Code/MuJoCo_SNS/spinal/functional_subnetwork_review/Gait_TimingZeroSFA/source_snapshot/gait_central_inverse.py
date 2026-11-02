"""Derived spike-rate RG decoder; unit/held-out validation precedes any fitting.

The runtime decoder accepts observed spike rate and actual conductance/channel
states. It does not accept a graded target voltage or a direct input current.
"""
import json,math,hashlib,argparse,copy,time,inspect
from pathlib import Path
import numpy as np
from numba import njit
from scipy.integrate import quad
from scipy.optimize import brentq
import gait_central_meanfield as M
C=M.C;ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_central_inverse_v3')
NODES,WEIGHTS=np.polynomial.legendre.leggauss(24)
STATE={}

@njit(cache=True)
def activation(v,ka,sa,ea):return 1./(1.+ka*math.exp(sa*(ea-v)))

@njit(cache=True)
def nap_vector(v,current,conductance,gna,h,ka,sa,ea,ena):
    return current-conductance*v+gna*h*activation(v,ka,sa,ea)*(ena-v)

@njit(cache=True)
def passage_time(current,cm,conductance,gna,h,ka,sa,ea,ena,reset=2.5,threshold=4.):
    result=0.;half=(threshold-reset)/2.;center=(threshold+reset)/2.
    for j in range(len(NODES)):
        v=center+half*NODES[j];speed=nap_vector(v,current,conductance,gna,h,ka,sa,ea,ena)
        if speed<=0:return 1e9
        result+=WEIGHTS[j]/speed
    return cm*half*result

@njit(cache=True)
def passive_inverse(rate,cm,conductance,reset=2.5,threshold=4.,refractory=.001,clock_correction=0.):
    if rate<=0:return -4.
    passage=max(1e-6,1/rate-refractory-clock_correction);z=math.exp(-conductance*passage/cm)
    return (threshold-reset*z)/max(1e-12,1-z)

@njit(cache=True)
def inverse_nap_rate(rate,cm,conductance,gna,h,ka,sa,ea,ena,clock_correction=.0001):
    if rate<.1:return -4. # A silent subthreshold voltage is not identifiable from no spikes.
    if gna*h<1e-12:return passive_inverse(rate,cm,conductance,clock_correction=clock_correction)
    target=max(1e-6,1/rate-.001-clock_correction)
    critical=-1e30
    for j in range(33):
        v=2.5+1.5*j/32.;critical=max(critical,conductance*v-gna*h*activation(v,ka,sa,ea)*(ena-v))
    lo=critical+1e-10;hi=lo+cm*1.5/target+gna*abs(ena)+conductance*4+1.
    for _ in range(8):
        if passage_time(hi,cm,conductance,gna,h,ka,sa,ea,ena)<target:break
        hi=lo+2*(hi-lo)
    for _ in range(32):
        mid=(lo+hi)/2.
        if passage_time(mid,cm,conductance,gna,h,ka,sa,ea,ena)>target:lo=mid
        else:hi=mid
    current=(lo+hi)/2.;lo=4.;hi=max(9.,(current+gna*abs(ena)+10)/conductance)
    for _ in range(36):
        mid=(lo+hi)/2.
        if nap_vector(mid,current,conductance,gna,h,ka,sa,ea,ena)>0:lo=mid
        else:hi=mid
    return (lo+hi)/2.

@njit(cache=True)
def fixed_channel_rate(current,cm,conductance,gna,h,ka,sa,ea,ena,duration=3.,dt=.0002):
    v=2.5;ref=0.;events=0;first=-1.;last=-1.
    for step in range(round(duration/dt)):
        if ref>dt/2:ref=max(0.,ref-dt);v=2.5;continue
        ion=gna*h*activation(v,ka,sa,ea);g=conductance+ion;steady=(current+ion*ena)/g
        v=steady+(v-steady)*math.exp(-g*dt/cm)
        if v>=4.:
            v=2.5;ref=.001
            if step*dt>=1.:
                if first<0:first=step*dt
                events+=1;last=step*dt
    return (events-1)/(last-first) if events>2 else 0.

@njit(cache=True)
def fixed_population_rate(current,cm,conductance,gna,h,ka,sa,ea,ena):
    total=0.
    for cell in range(16):
        threshold=4*(1+.02*(cell/15.-.5));v=2.5;ref=0.;events=0;first=-1.;last=-1.
        for step in range(15000):
            if ref>.0001:ref=max(0.,ref-.0002);v=2.5;continue
            ion=gna*h*activation(v,ka,sa,ea);g=conductance+ion;steady=(current+ion*ena)/g
            v=steady+(v-steady)*math.exp(-g*.0002/cm)
            if v>=threshold:
                v=2.5;ref=.001
                if step>=5000:
                    if first<0:first=step*.0002
                    events+=1;last=step*.0002
        total+=(events-1)/(last-first) if events>2 else 0.
    return total/16.

def validate():
    ROOT.mkdir(parents=True,exist_ok=True)
    protocol=dict(machine='easteregg2',derivation='C*dV/dt=J−(gLeak+gSyn)*V+gNa*h*mInf(V)*(ENa−V). For fixedh, inferJ from actual firing period tref+C∫reset^threshold dV/F(V;J), then solve stable supra-threshold F(V;J)=0. A pure passiveLIF reduces to Vinf=(threshold−reset*exp(−G*(1/f−tref)/C))/(1−exp(...)). Inputcurrent and reversals supply independent unit ground truth, never runtime outputcarrier.',runtime_signals='Observed filtered spikeHz plus actual h/channel/conductance state; no original graded voltage, no inputcurrent argument. Sourcefastm nonlinear over the interspike interval is integrated rather than averaging gNa at a reset-biased soma.',clock_s=.0002,clock_inverse_correction_s=.0001,correction_scope='Midpoint of one200us threshold-crossing quantization bin; declared approximation, not a simulation timestep change.',identifiability='No-spike subthresholdpositive voltages are not inferable; zero-rate maps silent−4. Firing domain only. No invented positive activation from direct inputcurrent.',heldout=dict(passive_G_uS=[1.,1.75,3.],NaP_h=[.05,.2,.8],target_steady_mV=[4.1,4.56,5.25,6.25],cm=.05,gNa=12.,ENa=8.),gates=dict(passive_roundtrip_mV_max=1e-8,independent_scipy_integral_NaP_roundtrip_mV_max=.002,actual_200us_chemical_fraction_error_max=.05,actual_200us_hInf_error_max=.05),physical_trials=0,fit='No search permitted unless independently derived transfer passes the held-out firing drives. Frozen nominal comparison follows unitvalidation; sourcev1/v2 and all gates preserved.',source_v2=str(M.ROOT),script_sha256=C.X.digest(__file__))
    C.X.write(ROOT/'READOUT_PROTOCOL.json',protocol)
    passive=[]
    for g in (1.,1.75,3.):
        for target in (4.1,4.56,5.25,6.25):
            period=.001+.05/g*math.log((target-2.5)/(target-4.));inferred=passive_inverse(1/period,.05,g)
            passive.append(dict(g_total_uS=g,target_mV=target,inferred_mV=float(inferred),absolute_error_mV=abs(float(inferred)-target)))
    cases=[]
    for g in (1.,1.75,3.):
        for h in (.05,.2,.8):
            for target in (4.1,4.56,5.25,6.25):
                j=g*target-12*h/(1+math.exp(.8*(2-target)))*(8-target)
                vector=lambda v:j-g*v+12*h/(1+math.exp(.8*(2-v)))*(8-v)
                assert min(vector(v) for v in np.linspace(2.5,4.,200))>0
                exact_passage=.05*quad(lambda v:1/vector(v),2.5,4.,epsabs=1e-12,epsrel=1e-12)[0];exact_rate=1/(.001+exact_passage)
                inferred=inverse_nap_rate(exact_rate,.05,g,12.,h,1.,.8,2.,8.,0.)
                measured=fixed_channel_rate(j,.05,g,12.,h,1.,.8,2.,8.);decoded=inverse_nap_rate(measured,.05,g,12.,h,1.,.8,2.,8.)
                chemical_error=abs(np.clip(decoded/5,0,1)-np.clip(target/5,0,1));hinf=lambda v:1/(1+math.exp(2*(v-3.5)))
                poprate=fixed_population_rate(j,.05,g,12.,h,1.,.8,2.,8.);popdecoded=inverse_nap_rate(poprate,.05,g,12.,h,1.,.8,2.,8.)
                cases.append(dict(g_leak_plus_syn_uS=g,h=h,target_mV=target,effective_current_nA=j,ge_example_nA=8*(g-1),input_current_example_nA=j-8*(g-1),continuous_expected_rate_Hz=exact_rate,continuous_roundtrip_error_mV=float(abs(inferred-target)),actual_200us_rate_Hz=float(measured),actual_200us_decoded_mV=float(decoded),chemical_fraction_error=float(chemical_error),hInf_error=float(abs(hinf(decoded)-hinf(target))),heterogeneous16cell_rate_Hz=float(poprate),heterogeneous16cell_chemical_error=float(abs(np.clip(popdecoded/5,0,1)-np.clip(target/5,0,1))),heterogeneous16cell_hInf_error=float(abs(hinf(popdecoded)-hinf(target)))))
    passed=max(x['absolute_error_mV'] for x in passive)<1e-8 and max(x['continuous_roundtrip_error_mV'] for x in cases)<.002 and max(x['chemical_fraction_error'] for x in cases)<.05 and max(x['hInf_error'] for x in cases)<.05
    passed=passed and max(x['heterogeneous16cell_chemical_error'] for x in cases)<.05 and max(x['heterogeneous16cell_hInf_error'] for x in cases)<.05
    report=dict(derived_transfer_valid=bool(passed),passive=passive,heldout=cases,nonfiring_output_mV=float(inverse_nap_rate(0.,.05,2.,12.,.5,1.,.8,2.,8.)),maximum_continuous_roundtrip_error_mV=max(x['continuous_roundtrip_error_mV'] for x in cases),maximum_actual_clock_chemical_error=max(x['chemical_fraction_error'] for x in cases),maximum_actual_clock_hInf_error=max(x['hInf_error'] for x in cases),maximum_population_chemical_error=max(x['heterogeneous16cell_chemical_error'] for x in cases),maximum_population_hInf_error=max(x['heterogeneous16cell_hInf_error'] for x in cases),derivation_sha256=hashlib.sha256((inspect.getsource(inverse_nap_rate.py_func)+inspect.getsource(passage_time.py_func)+inspect.getsource(nap_vector.py_func)).encode()).hexdigest(),scope='Held-out fixedh/steadycurrent firing neurons and16cell heterogeneous populations, not dynamic fullnetwork or physical validation.',physical_trials=0)
    C.X.write(ROOT/'READOUT_VALIDATION.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('passive','heldout')}),flush=True);return report

@njit(cache=True)
def inverse_population_step(lv,h,ref,adapt,decoded,counts,params,gs,ge,applied,config,dt):
    for p in range(len(params)):
        cm,gm,rest,bias,gna,ka,sa,ea,kb,sb,eb,th,ena,is_rg=params[p]
        reset=2.5 if is_rg else 0.;rate=config[1]/cm;spikes=0
        h_voltage=inverse_nap_rate(decoded[p]*rate,cm,gm+gs[p],gna,np.mean(h[p]),ka,sa,ea,ena) if is_rg else 0.
        for cell in range(lv.shape[1]):
            old=lv[p,cell];adapt[p,cell]*=math.exp(-dt/.35)
            if is_rg:
                inf=1/(1+kb*np.exp(sb*(eb-h_voltage)));h[p,cell]+=(inf-h[p,cell])*dt/th
                m=activation(old,ka,sa,ea);gion=gna*m*h[p,cell]
            else:gion=0.
            if ref[p,cell]>dt/2:ref[p,cell]=max(0.,ref[p,cell]-dt);lv[p,cell]=reset;continue
            g=gm+gs[p]+gion;steady=(gm*rest+bias+applied[p]+ge[p]+gion*ena-adapt[p,cell])/g
            lv[p,cell]=steady+(old-steady)*math.exp(-g*dt/cm)
            threshold=(4. if is_rg else 1.)*(1.+(.02 if is_rg else .15)*(cell/(lv.shape[1]-1)-.5))
            if lv[p,cell]>=threshold:
                lv[p,cell]=reset;ref[p,cell]=.001;spikes+=1
                if is_rg:adapt[p,cell]+=config[0]
        decoded[p]=decoded[p]*math.exp(-dt/.025)+spikes/(lv.shape[1]*.025)*(1./rate if is_rg else config[2]);counts[p]+=spikes

class InverseCentral(M.MeanFieldCentral):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.decoded[self.is_rg]=0. # Empty spike history; preserve initial chemical silence, not an invented positive voltage.
    def proxy_voltage(self):
        v=super().proxy_voltage()
        for p,i in enumerate(self.indices):
            if not self.is_rg[p]:continue
            cm,gm,rest,bias,gna,ka,sa,ea,kb,sb,eb,th,ena,_=self.params[p]
            v[p]=inverse_nap_rate(self.decoded[p]*self.config[1]/cm,cm,gm+self.inner.sums[i,0],gna,float(self.h[p].mean()),ka,sa,ea,ena)
        return v
    def forward(self,x=None):
        if not self.enabled:return self.inner.forward(x)
        c=self.inner.c;c.V[self.indices]=self.proxy_voltage();self.inner.forward(x)
        if self.poison_shadow:c.V[self.indices]=1e6
        self.actual_input[:]=self.inner.i_app[self.indices];self.counts[:]=0
        gs=self.inner.sums[self.indices,0];ge=self.inner.sums[self.indices,1]
        for _ in range(10):inverse_population_step(self.lv,self.h,self.ref,self.adapt,self.decoded,self.counts,self.params,gs,ge,self.actual_input,self.config,.0002)
        c.V[self.indices]=self.proxy_voltage();self.total_counts+=self.counts;self.steps+=1
        if self.inner.output_onehot:np.take(c.V,self.inner.oq,out=self.inner.output_values);c.outputs=self.inner.output_values
        else:c.outputs=c.output_voltage_connectivity@c.V
        return c.outputs
    def export_state(self):
        state=C.CentralForward.export_state(self);state['inactivation_voltage_policy']=np.array('inverse_nap_spike_rate_envelope');state['source_population_params']=self.params.copy();state['source_bias']=self.inner.c.i_b.copy();return state
    def import_state(self,state):
        assert str(state['inactivation_voltage_policy'])=='inverse_nap_spike_rate_envelope';C.CentralForward.import_state(self,state);np.copyto(self.params,state['source_population_params']);np.copyto(self.inner.c.i_b,state['source_bias'])

def prepare_comparison():
    validation=json.loads((ROOT/'READOUT_VALIDATION.json').read_text());assert validation['derived_transfer_valid'],'No network experiment before heldout transfer validation'
    expected=hashlib.sha256((inspect.getsource(inverse_nap_rate.py_func)+inspect.getsource(passage_time.py_func)+inspect.getsource(nap_vector.py_func)).encode()).hexdigest();assert validation['derivation_sha256']==expected
    protocol=dict(machine='easteregg2',frozen_configuration=C.DEFAULT,mechanism='Only RG readout coordinate changes from linear rate map to derived nonlinear conductance/NaP-aware inverse; h still uses the actual spike-derived population envelope. Actualcell sourceNaP, reset/refractory/adaptation,PFdecoder,chemicalweights,200usclock unchanged.',initial_state='RG filter initialized with zero spikehistory. Old linearcodec encoded sourceV0 using a nonzero0.444fraction; no positive firingrate can encode subthresholdV0 in the inversecodec. Both initial chemicalgates are silent; rawinitialproxy−4 vs0 is explicitly different. Somamembrane/h states unchanged. No hidden graded target or directcurrent suppliesoutput.',RG_rate_normalization='The inherited RG_C_rate5 merely changes filter coordinate and algebraically cancels in Hz reconstruction; it is fixed, not an optimization knob.',unchanged_acceptance=json.loads((M.ROOT/'PROTOCOL.json').read_text())['unchanged_gates'],h_gates=dict(standing_mean_error_max=.15,driven_peak_phase_error_max=.15),controls='Fresh unchangedv2 must reproduce all frozen traces/events/somas/h. New decoder disabled/default+complete-state+shadow-poison checks; incompatiblepolicy state rejected; NumBa/scalar localcheck.',physical_gate='Both offline fixtures and h gates required, no fit unlessheldoutderivationvalid, no newphysicalcondition; no bodybydefault.',source_v2_adapter_sha256=C.X.digest(M.__file__),derivation_sha256=expected,script_sha256=C.X.digest(__file__))
    C.X.write(ROOT/'COMPARISON_PROTOCOL.json',protocol)
    inner,names=C.X.construct();M.STATE.update(source=inner,names=names,inputs={},reference={},graded_h={})
    path=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait2392_reflex_demo_v2\cases\syn6_graded_baseline.inputs.npz')
    with np.load(path) as z:M.STATE['inputs']['standing']=z['inputs'].copy()
    u=np.zeros((5000,M.STATE['inputs']['standing'].shape[1]));u[:,np.flatnonzero(inner.c.input_connectivity[names.index('DRIVE')])[0]]=2.5;M.STATE['inputs']['deaff']=u
    with np.load(M.SOURCE/'GRADED_REFERENCE.npz') as z:M.STATE['reference']={k:dict(trace=z[k].copy()) for k in ('standing','deaff')}
    with np.load(M.ROOT/'MEANFIELD_FROZEN.npz') as z:M.STATE['graded_h']={k:z[k+'_graded_h'].copy() for k in ('standing','deaff')}
    a=C.clone(inner);b=InverseCentral(C.clone(inner),names,enabled=False)
    for u in M.STATE['inputs']['standing'][:64]:assert np.array_equal(a.forward(u),b.forward(u))
    a=InverseCentral(C.clone(inner),names);b=InverseCentral(C.clone(inner),names,poison_shadow=True)
    for u in M.STATE['inputs']['standing'][:64]:assert np.array_equal(a.forward(u),b.forward(u))
    start=a.export_state()
    for u in M.STATE['inputs']['standing'][64:96]:a.forward(u)
    first=a.export_state();a.import_state(start)
    for u in M.STATE['inputs']['standing'][64:96]:a.forward(u)
    second=a.export_state();checks={k:bool(np.array_equal(v,second[k],equal_nan=True)) if v.dtype.kind in 'fc' else bool(np.array_equal(v,second[k])) for k,v in first.items()};assert all(checks.values())
    rejected=False
    try:a.import_state(M.MeanFieldCentral(C.clone(inner),names).export_state())
    except (AssertionError,KeyError):rejected=True
    assert rejected
    arrays=[a.lv.copy(),a.h.copy(),a.ref.copy(),a.adapt.copy(),a.decoded.copy(),a.counts.copy()];other=[v.copy() for v in arrays];gs=a.inner.sums[a.indices,0];ge=a.inner.sums[a.indices,1]
    inverse_population_step(*arrays,a.params,gs,ge,a.actual_input,a.config,.0002);inverse_population_step.py_func(*other,a.params,gs,ge,a.actual_input,a.config,.0002)
    errors={k:float(np.max(abs(x-y))) for k,x,y in zip(('soma','h','refractory','adaptation','spikefilter','counts'),arrays,other)};assert max(errors.values())<1e-11
    C.X.write(ROOT/'STATE_AUDIT.json',dict(default_disabled64_exact=True,discarded_shadow_poison64_exact=True,complete_continuation32_exact=True,incompatible_v2_policy_rejected=True,compiled_vs_scalar_one200us_step_errors=errors,arrays=checks,physical_trials=0));print(json.dumps(dict(prepared=True,unit_validation_passed=True,state_audits_passed=True)),flush=True)

def compare():
    prepare_comparison();results={};traces={}
    for label,klass in (('V2_CONTROL',M.MeanFieldCentral),('INVERSE_FROZEN',InverseCentral)):
        rows={};started=time.perf_counter();all_data={}
        original=M.MeanFieldCentral;M.MeanFieldCentral=klass
        try:
            for key,u in M.STATE['inputs'].items():
                a,_=M.replay(C.clone(M.STATE['source']),M.STATE['names'],u,C.DEFAULT,fixture=key=='deaff');r=C.metrics(a,M.STATE['reference'][key],key)
                if key=='standing':r['RG_h_mean_absolute_error']=np.mean(abs(a['h_mean'][500:,[0,1,8,9]]-M.STATE['graded_h'][key][500:,[0,1,8,9]]),axis=0).tolist()
                else:
                    hp=M.h_profile(a['trace'],a['h_mean']);hs=M.h_profile(M.STATE['reference'][key]['trace'],M.STATE['graded_h'][key]);r['h_peak_cycle_phase_error']=(np.full(4,.5) if hp is None else abs((np.argmax(hp,axis=1)/64-np.argmax(hs,axis=1)/64+.5)%1-.5)).tolist()
                rows[key]=r;all_data[key]=a
        finally:M.MeanFieldCentral=original
        if label=='V2_CONTROL':
            with np.load(M.ROOT/'MEANFIELD_FROZEN.npz') as z:
                for key,a in all_data.items():
                    for field in ('trace','events','somas','h_mean'):assert np.array_equal(a[field],z[key+'_'+field]),(key,field)
        passed=all(r['passed'] for r in rows.values()) and max(rows['standing']['RG_h_mean_absolute_error'])<=.15 and max(rows['deaff']['h_peak_cycle_phase_error'])<=.15
        result=dict(configuration=C.DEFAULT,passed_offline_gate=passed,fixtures=rows,wall_s=time.perf_counter()-started,physical_trials=0)
        C.X.write(ROOT/(label+'.json'),result);np.savez_compressed(ROOT/(label+'.npz'),**{k+'_'+f:a[f] for k,a in all_data.items() for f in ('trace','events','somas','h_mean')},**{k+'_graded_h':v for k,v in M.STATE['graded_h'].items()},configuration_json=json.dumps(C.DEFAULT),dt_s=.002)
        results[label]=result;print(json.dumps(dict(case=label,**result)),flush=True)
    (ROOT/'EXPERIMENT_RECORD.md').write_text('# Conductance-aware RG spike readout\n\nMachine: easteregg2. Derived nonlinear inverse validated before any networkcomparison: runtime uses spikeHz and actualconductance/channelstates, nevergradedshadowor directinputcurrent. Frozenv2control exactly reproduces its traces/events/somas/h. Same sourceweights, physicalcondition,200µsclock,reset/adaptation,meanfieldhmechanism, and defaultparameters. Decodercoordinate initialization differs explicitly to preserve initialchemicalsilence with an emptyspikehistory.\n\nInversefrozen offline gate: '+str(results['INVERSE_FROZEN']['passed_offline_gate'])+'. No gainsearch or bodytrial was launched. See READOUT_VALIDATION.json, COMPARISON_PROTOCOL.json, STATE_AUDIT.json, and INVERSE_FROZEN.json.\n',encoding='utf-8')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--validate',action='store_true');p.add_argument('--compare',action='store_true');a=p.parse_args()
    if a.validate:validate()
    if a.compare:compare()
