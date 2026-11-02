"""Offline-gated syn6 RG/PF population prototype; native MN/reflex cells stay graded.

NaP gates are copied from actual syn6 RGs, not the W2L calcium law. Only the
reset/refractory/adaptation + spike-filter population pattern is reused.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import sys,json,copy,time,hashlib,argparse,math
from pathlib import Path
import numpy as np
from numba import njit
from scipy.signal import find_peaks
import gait_spiking_reflex as X
F=X.F;ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_central_v1')
DEFAULT=dict(rg_adapt_jump_nA=.3,rg_readout_C_rate=5.,pf_rate_scale=1.)

@njit(cache=True)
def population_step(lv,h,ref,adapt,decoded,counts,params,gs,ge,applied,config,dt):
    for p in range(len(params)):
        cm,gm,rest,bias,gna,ka,sa,ea,kb,sb,eb,th,ena,is_rg=params[p]
        reset=2.5 if is_rg else 0.;rate=config[1]/cm;spikes=0
        for cell in range(lv.shape[1]):
            old=lv[p,cell];adapt[p,cell]*=math.exp(-dt/.35)
            if is_rg:
                inf=1/(1+kb*math.exp(sb*(eb-old)));h[p,cell]+=(inf-h[p,cell])*dt/th
                m=1/(1+ka*math.exp(sa*(ea-old)));gion=gna*m*h[p,cell]
            else:gion=0.
            if ref[p,cell]>dt/2:
                ref[p,cell]=max(0.,ref[p,cell]-dt);lv[p,cell]=reset;continue
            g=gm+gs[p]+gion;steady=(gm*rest+bias+applied[p]+ge[p]+gion*ena-adapt[p,cell])/g
            lv[p,cell]=steady+(old-steady)*math.exp(-g*dt/cm)
            threshold=(4. if is_rg else 1.)*(1.+(.02 if is_rg else .15)*(cell/(lv.shape[1]-1)-.5))
            if lv[p,cell]>=threshold:
                lv[p,cell]=reset;ref[p,cell]=.001;spikes+=1
                if is_rg:adapt[p,cell]+=config[0]
        decoded[p]=decoded[p]*math.exp(-dt/.025)+spikes/(lv.shape[1]*.025)*(1./rate if is_rg else config[2]);counts[p]+=spikes

class CentralForward:
    def __init__(self,inner,names,configuration=None,enabled=True,poison_shadow=False):
        self.inner=inner;self.names=list(names);self.enabled=enabled;self.poison_shadow=poison_shadow
        self.configuration=dict(DEFAULT if configuration is None else configuration);assert set(self.configuration)==set(DEFAULT)
        self.config=np.array([self.configuration[k] for k in DEFAULT]);assert np.all(np.isfinite(self.config)) and np.all(self.config>0)
        c=inner.c;self.indices=np.array([i for i,n in enumerate(names) if n.startswith(('RG_','PF_S'))],dtype=np.int64);assert len(self.indices)==16
        self.is_rg=np.array([names[i].startswith('RG_') for i in self.indices]);assert self.is_rg.sum()==4
        self.params=np.array([[c.c_m[i],c.g_m[i],c.V_rest[i],c.i_b[i],c.g_ion[0,i],c.k_a[0,i],c.slope_a[0,i],c.e_a[0,i],c.k_b[0,i],c.slope_b[0,i],c.e_b[0,i],c.tau_max_b[0,i],c.e_ion[0,i],self.is_rg[p]] for p,i in enumerate(self.indices)])
        assert np.all(c.pow_c[:,self.indices]==0);assert c.dt==.002
        self.lv=np.repeat(c.V[self.indices,None],16,axis=1)+np.linspace(-.01,.01,16)[None,:]
        self.h=np.repeat(c.b_gate[0,self.indices,None],16,axis=1);self.ref=np.zeros_like(self.lv);self.adapt=self.ref.copy()
        self.decoded=np.where(self.is_rg,(c.V[self.indices]+4.)/9.,0.);self.counts=np.zeros(16,dtype=np.int64);self.total_counts=self.counts.copy();self.steps=0;self.actual_input=np.zeros(16)
        population_step(self.lv.copy(),self.h.copy(),self.ref.copy(),self.adapt.copy(),self.decoded.copy(),self.counts.copy(),self.params,np.zeros(16),np.zeros(16),np.zeros(16),self.config,.0002)
    def __getattr__(self,name):return getattr(self.inner,name)
    def proxy_voltage(self):
        voltage=-4.+9.*self.decoded
        for p,i in enumerate(self.indices):
            if self.is_rg[p]:continue
            rate=self.decoded[p]
            if rate<.1:voltage[p]=-4.;continue
            conductance=self.params[p,1]+self.inner.sums[i,0];interval=max(.00005,1./rate-.001)
            z=math.exp(-conductance*interval/self.params[p,0]);voltage[p]=1./max(1e-8,1-z)
        return voltage
    def forward(self,x=None):
        if not self.enabled:return self.inner.forward(x)
        c=self.inner.c;c.V[self.indices]=self.proxy_voltage()
        self.inner.forward(x) # All outgoing chemical gates used the spike-derived proxy, never this shadow update.
        if self.poison_shadow:c.V[self.indices]=1e6
        self.actual_input[:]=self.inner.i_app[self.indices];self.counts[:]=0
        gs=self.inner.sums[self.indices,0];ge=self.inner.sums[self.indices,1]
        for _ in range(10):population_step(self.lv,self.h,self.ref,self.adapt,self.decoded,self.counts,self.params,gs,ge,self.actual_input,self.config,.0002)
        c.V[self.indices]=self.proxy_voltage();self.total_counts+=self.counts;self.steps+=1
        if self.inner.output_onehot:np.take(c.V,self.inner.oq,out=self.inner.output_values);c.outputs=self.inner.output_values
        else:c.outputs=c.output_voltage_connectivity@c.V
        return c.outputs
    def export_state(self):
        state=X.raw_snapshot(self.inner);state.update(central_version=np.array(1),central_config=np.array(json.dumps(self.configuration,sort_keys=True)),central_names=np.array(self.names),central_enabled=np.array(self.enabled),central_steps=np.array(self.steps))
        for key in ('lv','h','ref','adapt','decoded','counts','total_counts','actual_input'):state['central_'+key]=getattr(self,key).copy()
        return state
    def import_state(self,state):
        assert int(state['central_version'])==1 and str(state['central_config'])==json.dumps(self.configuration,sort_keys=True) and bool(state['central_enabled'])==self.enabled
        assert np.array_equal(state['central_names'],self.names);X.raw_restore(self.inner,state);self.steps=int(state['central_steps'])
        for key in ('lv','h','ref','adapt','decoded','counts','total_counts','actual_input'):np.copyto(getattr(self,key),state['central_'+key])

STATE={}
def clone(inner):return X.SparsePairwiseForward(copy.deepcopy(inner.c))
def run_replay(inner,names,inputs,config=None,graded=False,fixture=False):
    controller=inner if graded else CentralForward(inner,names,config);ids=np.array([i for i,n in enumerate(names) if n.startswith(('RG_','PF_S'))]);mn=np.array([names.index('MN_'+m+'_r') for m in ('iliacus','bifemlh','rect_fem','vas_int','soleus','tib_ant')])
    selected=np.r_[ids,mn];trace=np.empty((len(inputs),len(selected)));events=np.zeros((len(inputs),16),dtype=np.int64);somas=np.empty((len(inputs),16,2)) if not graded else np.empty((0,16,2))
    peak_adapt=0.;minimum_h=1.;maximum_h=0.
    for k,u in enumerate(inputs):
        if fixture:inner.c.i_b[names.index('RG_F_l')]=1. if 50<=k<100 else 0.
        if not graded and fixture:controller.params[1,3]=inner.c.i_b[names.index('RG_F_l')]
        controller.forward(u);trace[k]=inner.c.V[selected]
        if not graded:
            events[k]=controller.counts;somas[k]=controller.lv[:,:2];peak_adapt=max(peak_adapt,float(controller.adapt.max()));minimum_h=min(minimum_h,float(controller.h.min()));maximum_h=max(maximum_h,float(controller.h.max()))
    return dict(trace=trace,events=events,somas=somas,selected=selected,peak_adapt_nA=peak_adapt,h_range=[minimum_h,maximum_h]),controller

def prepare():
    ROOT.mkdir(parents=True,exist_ok=True)
    protocol=dict(machine='easteregg2',configuration=str(X.BEST),configuration_sha256=X.digest(X.BEST),source_diagram=str(X.DIAGRAM),source_diagram_sha256=X.digest(X.DIAGRAM),source_equations='Four actual fixed-tau NaP RG cells and twelve graded synergy PF cells replaced by16-cell reset/refractory populations. RG copies source NaP gates plus explicit shared spike-triggered adaptation; PF uses conductance LIF. Not calcium neurons; W2L supplies population/reset/adaptation implementation pattern only.',architecture='256 actual central LIF cells;778 remaining graded cells;16 logical source positions are spike-filter readout proxies. MNs, sensory/reflex cells and all original synaptic conductances unchanged, including fittedPF gain0.013344897908126326, notzero.',readout='Only spike events drive25ms exponential synaptic filters: Vproxy=−4+9*decoded, max-rate normalization C_rate/sourceCm. No graded central presynaptic voltages are retained. This offset/span and rate policy are new explicit engineering choices.',membrane='Relative-rest mV; threshold4; RGreset2.5,PFreset.5;1ms refractory,200µs physical neural substeps. Source RGh tau350ms retained. Adaptation tau350ms; PF adaptation0.',fixtures=['Actual1500-frame3s recorded syn6 graded standing input, including preexisting0.6s preparation; reference is tonic, not rhythmic.','Deafferented10s source network, DRIVE2.5nA, allotherinput0; explicit leftRG-F1nA startup0.1–0.2s. Source graded reference is almost bilateral in-phase, not accepted alternating stepping.'],gates=dict(default_disabled_exact_64steps=True,complete_state_exact_continuation=True,shadow_voltage_poison_exact=True,standing_fraction_RMSE_max=.15,standing_MN_RMSE_mV_max=.10,central_fraction_RMSE_max=.20,RG_period_relative_error_max=.20,RG_relative_phase_error_max=.10,finite=True,actual_spikes_active_pools=True),fit=dict(max_trials=24,shared_LR=True,keys={'rg_adapt_jump_nA':[.05,1.5],'rg_readout_C_rate':[2.,12.],'pf_readout_C_rate':[1.5,4.]}),physical_gate='No body integration until BOTH offline fixture gates and default/state/shadow audits pass. Standing success is not rhythmic gait; no new support or gain-zero substitution.',source_script=str(Path(__file__)))
    protocol['readout']='RG: spike-filter Vproxy=−4+9*decoded with C-normalized rate. PF: measured spikeHz filtered25ms then inverse passive-LIF rate equation Vinf=theta/(1−exp(−G*(1/f−tref)/Cm)); theta1/reset0, actual chemical+leakG only, no graded membrane/current target. Common PF rate scale learned. Negligible f<.1Hz maps to silent−4mV. Deterministic threshold heterogeneity ±1% RG/±7.5%PF mirrored across sides. These are explicit new engineering choices.'
    protocol['membrane']='Relative-rest mV; RGthreshold4/reset2.5,PFthreshold1/reset0 with declared heterogeneity;1ms refractory rounded at200µs; RGsourceNaPh tau350ms and adaptationtau350ms. PF noNaP/adaptation.'
    protocol['fit']['keys'].pop('pf_readout_C_rate');protocol['fit']['keys']['pf_rate_scale']=[.8,1.2]
    protocol['gates'].update(peak_adaptation_nA_max=100,slow_gate_range=[0,1]);protocol['superseded_probe']='INITIAL_LINEAR_READOUT_RESULT.json failed; PFthreshold4 couldnot recruit mixedPF. Preserved separately, no bodytrial.'
    X.write(ROOT/'PROTOCOL.json',protocol)
    (ROOT/'EXPERIMENT_RECORD.md').write_text('# Central-spiking syn6 prototype\n\nMachine: easteregg2. Offline feasibility only until the preregistered gates pass. Accepted standing reference is tonic RG-E with suppressed RG-F; it supplies no cadence target. A separately labelled deafferented-drive fixture is rhythmic but almost bilateral in-phase, so it is not alternating walking.\n\nActual source: `'+str(X.DIAGRAM)+'`; fixed standing configuration: `'+str(X.BEST)+'`. Four RG NaP and twelve graded PF cells are targeted. W2L reset/refractory/adaptation code pattern is reused; its calcium parameters and controller are not silently imported. See PROTOCOL.json for all architecture/decoder changes and frozen physical gate.\n',encoding='utf-8')
    original,names=X.construct();STATE['source']=original;STATE['names']=names
    saved=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait2392_reflex_demo_v2\cases\syn6_graded_baseline.inputs.npz')
    with np.load(saved) as z:standing=z['inputs'].copy()
    deaff=np.zeros((5000,standing.shape[1]));port=np.flatnonzero(original.c.input_connectivity[names.index('DRIVE')])[0];deaff[:,port]=2.5
    STATE['inputs']={'standing':standing,'deaff':deaff};STATE['reference']={}
    for key,u in STATE['inputs'].items():
        a,_=run_replay(clone(original),names,u,graded=True,fixture=key=='deaff');STATE['reference'][key]=a
    source=STATE['reference']['standing']['trace'];ids=STATE['reference']['standing']['selected'];i0=np.flatnonzero(ids==names.index('RG_E_l'))[0]
    X.write(ROOT/'SOURCE_AUDIT.json',dict(logical_source_neurons=len(names),central_names=[names[i] for i in ids[:16]],source_nonzero_PF_to_MN_gain=X.settings()['pf_mn_gain'],standing_RG_E_range_after1s=[float(source[500:,i0].min()),float(source[500:,i0].max())],recorded_input_path=str(saved),recorded_input_sha256=X.digest(saved),reference_scope='Standing tonic; deafferented source reference near in-phase, no locomotion or mechanics claim.'))
    np.savez_compressed(ROOT/'GRADED_REFERENCE.npz',standing=STATE['reference']['standing']['trace'],deaff=STATE['reference']['deaff']['trace'],names=np.array(names)[ids],dt_s=.002)
    a=clone(original);b=CentralForward(clone(original),names,enabled=False)
    for u in standing[:64]:assert np.array_equal(a.forward(u),b.forward(u))
    c=CentralForward(clone(original),names);d=CentralForward(clone(original),names,poison_shadow=True)
    for u in standing[:64]:assert np.array_equal(c.forward(u),d.forward(u))
    initial=c.export_state()
    for u in standing[64:96]:c.forward(u)
    first=c.export_state();c.import_state(initial)
    for u in standing[64:96]:c.forward(u)
    second=c.export_state();checks={k:bool(np.array_equal(v,second[k],equal_nan=True)) if v.dtype.kind in 'fc' else bool(np.array_equal(v,second[k])) for k,v in first.items()};assert all(checks.values())
    X.write(ROOT/'ELECTRICAL_STATE_AUDIT.json',dict(default_disabled64_exact=True,shadow_poison64_exact=True,continuation32_exact=True,arrays=checks,physical_simulations=0))
    print(json.dumps(dict(prepared=True,state_audits_pass=True,central_names=[names[i] for i in ids[:16]])),flush=True)

def metrics(a,ref,key):
    tail=slice(1000,None) if key=='deaff' else slice(500,None);fraction=np.clip(a['trace'][:,:16]/5.,0,1);target=np.clip(ref['trace'][:,:16]/5.,0,1)
    rmse=float(np.sqrt(np.mean((fraction[tail]-target[tail])**2)));mn=float(np.sqrt(np.mean((a['trace'][tail,16:]-ref['trace'][tail,16:])**2)));spikes=a['events'][tail].sum(axis=0);active=np.mean(target[tail],axis=0)>.1
    out=dict(finite=bool(np.isfinite(a['trace']).all() and np.isfinite(a['somas']).all()),chemical_fraction_RMSE=rmse,MN_RMSE_mV=mn,spikes=spikes.tolist(),all_source_active_populations_spike=bool(np.all(spikes[active]>0)),peak_adaptation_nA=a['peak_adapt_nA'],h_range=a['h_range']);gate=out['finite'] and out['all_source_active_populations_spike'] and a['peak_adapt_nA']<100 and min(a['h_range'])>=0 and max(a['h_range'])<=1 and rmse<(.2 if key=='deaff' else .15) and (key=='deaff' or mn<.1)
    if key=='deaff':
        cycles={};period_errors=[];phases=[]
        for i in (0,1,8,9):
            refp,_=find_peaks(target[tail,i],prominence=.3,distance=100);p,_=find_peaks(fraction[tail,i],prominence=.3,distance=100)
            target_period=float(np.median(np.diff(refp))*.002) if len(refp)>2 else None;period=float(np.median(np.diff(p))*.002) if len(p)>2 else None;error=abs(period/target_period-1) if target_period and period else 1.
            cycles[str(i)]=dict(count=len(p),source_count=len(refp),period_s=period,source_period_s=target_period,relative_error=error);period_errors.append(error)
        # Relative left/right RG-F phase is compared with the measured source, not assumed antiphase.
        for v in (target[tail],fraction[tail]):
            left,_=find_peaks(v[:,1],prominence=.3,distance=100);right,_=find_peaks(v[:,9],prominence=.3,distance=100);errs=[]
            for l,r in zip(left[:-1],left[1:]):
                choices=right[(right>=l)&(right<r)]
                if len(choices):errs.append((choices[0]-l)/(r-l))
            phases.append(float(np.angle(np.mean(np.exp(2j*np.pi*np.array(errs))))/(2*np.pi)%1) if errs else None)
        phase_error=abs((phases[1]-phases[0]+.5)%1-.5) if None not in phases else .5
        out.update(cycles=cycles,source_RG_phase=phases[0],spiking_RG_phase=phases[1],phase_error=phase_error);gate=gate and max(period_errors)<=.2 and phase_error<=.1 and all(c['count']>=3 for c in cycles.values())
    out['passed']=bool(gate);return out

def evaluate(config,save=False):
    if not STATE:prepare()
    rows={};data={};started=time.perf_counter()
    for key,u in STATE['inputs'].items():
        a,c=run_replay(clone(STATE['source']),STATE['names'],u,config,fixture=key=='deaff');rows[key]=metrics(a,STATE['reference'][key],key);data[key]=a
    passed=all(v['passed'] for v in rows.values());value=sum(v['chemical_fraction_RMSE'] for v in rows.values())+.2*rows['standing']['MN_RMSE_mV']
    value+=sum(min(2.,c['relative_error']) for c in rows['deaff']['cycles'].values())+2*rows['deaff']['phase_error']
    if save:
        np.savez_compressed(ROOT/'PROTOTYPE_TRACE.npz',**{k+'_'+field:a[field] for k,a in data.items() for field in ('trace','events','somas')},names=np.array(STATE['names'])[data['standing']['selected']],configuration_json=json.dumps(config),dt_s=.002)
        X.write(ROOT/'PROTOTYPE_RESULT.json',dict(configuration=config,passed_offline_gate=passed,fixtures=rows,objective=value,wall_s=time.perf_counter()-started,physical_simulations=0))
    return dict(ok=True,value=float(value),diagnostics=dict(passed_offline_gate=passed,fixtures=rows,wall_s=time.perf_counter()-started))

def handle_request(request):
    if request['op']=='preflight':
        if not STATE:prepare()
        return dict(ok=True,ready_for_optimization=True,reason='Offline source/default/state/shadow gates passed; bounded calibration only. Body gate remains closed.',source_sha256=X.digest(X.DIAGRAM))
    if request['op']=='evaluate':return evaluate(request['params'])
    raise ValueError(request['op'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--probe',action='store_true');p.add_argument('--manifest',action='store_true');a=p.parse_args();prepare()
    if a.probe:print(json.dumps(evaluate(DEFAULT,save=True)),flush=True)
    if a.manifest:
        manifest=dict(adapter_path=str(Path(__file__)),output_dir=str(ROOT/'rustuna24'),search_space={'rg_adapt_jump_nA':dict(low=.05,high=1.5),'rg_readout_C_rate':dict(low=2.,high=12.),'pf_rate_scale':dict(low=.8,high=1.2)},fixed_parameters={},seed=100101,max_trials=24,wall_budget_s=300,trial_timeout_s=90,initial_candidates=[DEFAULT,dict(rg_adapt_jump_nA=1.,rg_readout_C_rate=5.,pf_rate_scale=1.)])
        X.write(ROOT/'RUSTUNA_MANIFEST.json',manifest);print(json.dumps(dict(manifest=str(ROOT/'RUSTUNA_MANIFEST.json'))),flush=True)
