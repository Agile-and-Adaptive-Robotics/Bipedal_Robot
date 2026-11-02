"""One explicitly isolated NaP-inactivation signal comparison; no body by default."""
import json,copy,time,hashlib,argparse
from pathlib import Path
import numpy as np
from numba import njit
from scipy.signal import find_peaks
import gait_central_spiking as C
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_central_meanfield_v2')
SOURCE=C.ROOT;STATE={}

@njit(cache=True)
def meanfield_step(lv,h,ref,adapt,decoded,counts,params,gs,ge,applied,config,dt):
    for p in range(len(params)):
        cm,gm,rest,bias,gna,ka,sa,ea,kb,sb,eb,th,ena,is_rg=params[p]
        reset=2.5 if is_rg else 0.;rate=config[1]/cm;spikes=0
        # Sole mechanism change: h sees the actual spike-derived RG envelope.
        h_voltage=-4.+9.*decoded[p]
        for cell in range(lv.shape[1]):
            old=lv[p,cell];adapt[p,cell]*=np.exp(-dt/.35)
            if is_rg:
                inf=1/(1+kb*np.exp(sb*(eb-h_voltage)));h[p,cell]+=(inf-h[p,cell])*dt/th
                m=1/(1+ka*np.exp(sa*(ea-old)));gion=gna*m*h[p,cell]
            else:gion=0.
            if ref[p,cell]>dt/2:
                ref[p,cell]=max(0.,ref[p,cell]-dt);lv[p,cell]=reset;continue
            g=gm+gs[p]+gion;steady=(gm*rest+bias+applied[p]+ge[p]+gion*ena-adapt[p,cell])/g
            lv[p,cell]=steady+(old-steady)*np.exp(-g*dt/cm)
            threshold=(4. if is_rg else 1.)*(1.+(.02 if is_rg else .15)*(cell/(lv.shape[1]-1)-.5))
            if lv[p,cell]>=threshold:
                lv[p,cell]=reset;ref[p,cell]=.001;spikes+=1
                if is_rg:adapt[p,cell]+=config[0]
        decoded[p]=decoded[p]*np.exp(-dt/.025)+spikes/(lv.shape[1]*.025)*(1./rate if is_rg else config[2]);counts[p]+=spikes

class MeanFieldCentral(C.CentralForward):
    def forward(self,x=None):
        if not self.enabled:return self.inner.forward(x)
        c=self.inner.c;c.V[self.indices]=self.proxy_voltage();self.inner.forward(x)
        if self.poison_shadow:c.V[self.indices]=1e6
        self.actual_input[:]=self.inner.i_app[self.indices];self.counts[:]=0
        gs=self.inner.sums[self.indices,0];ge=self.inner.sums[self.indices,1]
        for _ in range(10):meanfield_step(self.lv,self.h,self.ref,self.adapt,self.decoded,self.counts,self.params,gs,ge,self.actual_input,self.config,.0002)
        c.V[self.indices]=self.proxy_voltage();self.total_counts+=self.counts;self.steps+=1
        if self.inner.output_onehot:np.take(c.V,self.inner.oq,out=self.inner.output_values);c.outputs=self.inner.output_values
        else:c.outputs=c.output_voltage_connectivity@c.V
        return c.outputs
    def export_state(self):
        state=super().export_state();state['inactivation_voltage_policy']=np.array('actual_spike_filtered_population_envelope');return state
    def import_state(self,state):
        assert str(state['inactivation_voltage_policy'])=='actual_spike_filtered_population_envelope';super().import_state(state)

def replay(inner,names,inputs,configuration=None,fixture=False,policy='meanfield'):
    controller=(MeanFieldCentral if policy=='meanfield' else C.CentralForward)(inner,names,configuration)
    selected=np.r_[controller.indices,[names.index('MN_'+m+'_r') for m in ('iliacus','bifemlh','rect_fem','vas_int','soleus','tib_ant')]]
    a=dict(trace=np.empty((len(inputs),len(selected))),events=np.empty((len(inputs),16),dtype=np.int64),somas=np.empty((len(inputs),16,2)),h_mean=np.empty((len(inputs),16)),selected=selected,peak_adapt_nA=0.,h_range=[1.,0.])
    for k,u in enumerate(inputs):
        if fixture:
            bias=1. if 50<=k<100 else 0.;inner.c.i_b[names.index('RG_F_l')]=bias;controller.params[1,3]=bias
        controller.forward(u);a['trace'][k]=inner.c.V[selected];a['events'][k]=controller.counts;a['somas'][k]=controller.lv[:,:2];a['h_mean'][k]=controller.h.mean(axis=1)
        a['peak_adapt_nA']=max(a['peak_adapt_nA'],float(controller.adapt.max()));a['h_range'][0]=min(a['h_range'][0],float(controller.h.min()));a['h_range'][1]=max(a['h_range'][1],float(controller.h.max()))
    return a,controller

def h_profile(trace,h):
    peaks,_=find_peaks(np.clip(trace[1000:,1]/5,0,1),prominence=.3,distance=100);peaks=peaks+1000
    rows=[]
    for l,r in zip(peaks[:-1],peaks[1:]):rows.append(np.array([np.interp(np.linspace(l,r,64,endpoint=False),np.arange(len(h)),h[:,i]) for i in (0,1,8,9)]))
    if not rows:return None
    return np.mean(rows,axis=0)

def prepare():
    ROOT.mkdir(parents=True,exist_ok=True);prov=json.loads((SOURCE/'rustuna24/PROVENANCE.json').read_text());assert C.X.digest(C.__file__)==prov['adapter_sha256']
    protocol=dict(machine='easteregg2',source_v1=str(SOURCE),frozen_parameters=C.DEFAULT,mechanism='Only NaP h_inf voltage changes: actual threshold/reset soma versus actual spike-filter population envelope. m_inf still uses actual soma. h kinetics350ms,200us clock, all weights/decoder/reset/refractory/adaptation/initialization unchanged. No original graded central voltage supplies a chemical output or h signal.',physiology_scope='Mean-field spike-envelope dependent NaP availability is an explicit new phenomenological population mechanism, not an assertion of literal single-cell membrane voltage. Chemical presynaptic voltages always derive from spikes.',comparison='Fresh soma-h control must exactly reproduce saved v1 nominal traces/events/somas. Meanfield-h case uses same frozen nominal values on standing+10s deafferented reference, not new fitting.',unchanged_gates=json.loads((SOURCE/'PROTOCOL.json').read_text())['gates'],extra_h_gates=dict(standing_RG_h_mean_error_max=.15,deaff_h_peak_phase_error_max=.15),physical_gate='No body trial unless both original offline gates plus h and state/default/shadow audits pass. Failedv1 untouched. No optimizer launched by default.',source_adapter_sha256=C.X.digest(C.__file__),script=str(Path(__file__)),script_sha256=C.X.digest(__file__))
    C.X.write(ROOT/'PROTOCOL.json',protocol)
    if not (ROOT/'EXPERIMENT_RECORD.md').exists():
        (ROOT/'EXPERIMENT_RECORD.md').write_text('# Mean-field NaP mechanism comparison\n\nMachine: easteregg2. Separate explicit change from preserved failed v1: only h_inf now reads the actual spike-filter population envelope; fast NaP activation still reads actual spike soma. Source graded central outputs are discarded and never supply h or chemical outputs. Frozen shared parameters, source weights, 200 µs clock, physical condition, and decoder stay unchanged. No body trial is allowed before both source replay and additional h/state gates pass.\n\n',encoding='utf-8')
    inner,names=C.X.construct();STATE.update(source=inner,names=names,inputs={},reference={},graded_h={})
    inputfile=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait2392_reflex_demo_v2\cases\syn6_graded_baseline.inputs.npz')
    with np.load(inputfile) as z:STATE['inputs']['standing']=z['inputs'].copy()
    deaff=np.zeros((5000,STATE['inputs']['standing'].shape[1]));deaff[:,np.flatnonzero(inner.c.input_connectivity[names.index('DRIVE')])[0]]=2.5;STATE['inputs']['deaff']=deaff
    with np.load(SOURCE/'GRADED_REFERENCE.npz') as z:
        STATE['reference']={k:dict(trace=z[k].copy()) for k in ('standing','deaff')}
    ids=np.array([i for i,n in enumerate(names) if n.startswith(('RG_','PF_S'))])
    for key,u in STATE['inputs'].items():
        graded=C.clone(inner);h=np.empty((len(u),16))
        for k,row in enumerate(u):
            if key=='deaff':graded.c.i_b[names.index('RG_F_l')]=1. if 50<=k<100 else 0.
            graded.forward(row);h[k]=graded.c.b_gate[0,ids]
            assert np.array_equal(graded.c.V[ids],STATE['reference'][key]['trace'][k,:16])
        STATE['graded_h'][key]=h
    a=C.clone(inner);b=MeanFieldCentral(C.clone(inner),names,enabled=False)
    for u in STATE['inputs']['standing'][:64]:assert np.array_equal(a.forward(u),b.forward(u))
    a=MeanFieldCentral(C.clone(inner),names);b=MeanFieldCentral(C.clone(inner),names,poison_shadow=True)
    for u in STATE['inputs']['standing'][:64]:assert np.array_equal(a.forward(u),b.forward(u))
    start=a.export_state()
    for u in STATE['inputs']['standing'][64:96]:a.forward(u)
    first=a.export_state();a.import_state(start)
    for u in STATE['inputs']['standing'][64:96]:a.forward(u)
    second=a.export_state();checks={k:bool(np.array_equal(v,second[k],equal_nan=True)) if v.dtype.kind in 'fc' else bool(np.array_equal(v,second[k])) for k,v in first.items()};assert all(checks.values())
    C.X.write(ROOT/'STATE_AUDIT.json',dict(disabled_default64_exact=True,discarded_graded_voltage_poison64_exact=True,complete_continuation32_exact=True,graded_reference_recheck_all_central_frames_exact=True,arrays=checks,physical_trials=0));print(json.dumps(dict(prepared=True,state_audits_pass=True)),flush=True)

def evaluate(configuration,save=False,control=False):
    if not STATE:prepare()
    rows={};data={};started=time.perf_counter();h_pass=True
    for key,u in STATE['inputs'].items():
        a,c=replay(C.clone(STATE['source']),STATE['names'],u,configuration,fixture=key=='deaff',policy='soma' if control else 'meanfield');m=C.metrics(a,STATE['reference'][key],key)
        if key=='standing':
            error=np.mean(abs(a['h_mean'][500:,[0,1,8,9]]-STATE['graded_h'][key][500:,[0,1,8,9]]),axis=0);m['RG_h_mean_absolute_error']=error.tolist();h_pass=h_pass and max(error)<=.15
        else:
            source=h_profile(STATE['reference'][key]['trace'],STATE['graded_h'][key]);pred=h_profile(a['trace'],a['h_mean']);error=np.full(4,.5) if source is None or pred is None else abs((np.argmax(pred,axis=1)/64-np.argmax(source,axis=1)/64+.5)%1-.5);m['h_peak_cycle_phase_error']=error.tolist();h_pass=h_pass and max(error)<=.15
        rows[key]=m;data[key]=a
    if control:
        with np.load(SOURCE/'PROTOTYPE_TRACE.npz') as z:
            for key,a in data.items():
                for f in ('trace','events','somas'):assert np.array_equal(a[f],z[key+'_'+f]),(key,f)
    passed=all(r['passed'] for r in rows.values()) and h_pass;value=sum(r['chemical_fraction_RMSE'] for r in rows.values())+.2*rows['standing']['MN_RMSE_mV']+2*rows['deaff']['phase_error']+sum(min(2.,r['relative_error']) for r in rows['deaff']['cycles'].values())
    value+=2*max(rows['standing']['RG_h_mean_absolute_error'])+2*max(rows['deaff']['h_peak_cycle_phase_error'])
    result=dict(configuration=configuration,policy='actual_soma' if control else 'spike_filtered_envelope',passed_offline_gate=passed,h_gates_pass=bool(h_pass),fixtures=rows,wall_s=time.perf_counter()-started,physical_trials=0)
    if save:
        tag='SOMA_CONTROL' if control else 'MEANFIELD_FROZEN';C.X.write(ROOT/(tag+'.json'),result)
        np.savez_compressed(ROOT/(tag+'.npz'),**{k+'_'+field:a[field] for k,a in data.items() for field in ('trace','events','somas','h_mean')},**{k+'_graded_h':v for k,v in STATE['graded_h'].items()},names=np.array(STATE['names'])[data['standing']['selected']],configuration_json=json.dumps(configuration),dt_s=.002)
    return dict(ok=True,value=float(value),diagnostics=result)

def handle_request(request):
    if request['op']=='preflight':
        assert (ROOT/'MEANFIELD_FROZEN.json').exists() and (ROOT/'SOMA_CONTROL.json').exists(),'Frozen mechanism comparison must precede fitting'
        if not STATE:prepare()
        return dict(ok=True,ready_for_optimization=True,reason='Independent frozen mechanism comparison completed before bounded calibration; body gate remains closed.')
    if request['op']=='evaluate':return evaluate(request['params'])
    raise ValueError(request['op'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compare',action='store_true');p.add_argument('--manifest',action='store_true');a=p.parse_args();prepare()
    if a.compare:
        print(json.dumps(evaluate(C.DEFAULT,save=True,control=True)),flush=True);print(json.dumps(evaluate(C.DEFAULT,save=True)),flush=True)
    if a.manifest:
        assert (ROOT/'MEANFIELD_FROZEN.json').exists(),'Frozen mechanism comparison must precede any new calibration'
        cfg=dict(adapter_path=str(Path(__file__)),output_dir=str(ROOT/'rustuna24'),search_space={'rg_adapt_jump_nA':dict(low=.005,high=.4,log=True),'rg_readout_C_rate':dict(low=2.,high=8.),'pf_rate_scale':dict(low=.9,high=1.1)},fixed_parameters={},seed=100102,max_trials=24,wall_budget_s=300,trial_timeout_s=90,initial_candidates=[C.DEFAULT,dict(rg_adapt_jump_nA=.02,rg_readout_C_rate=3.5,pf_rate_scale=1.)])
        C.X.write(ROOT/'RUSTUNA_MANIFEST.json',cfg)
        protocol=json.loads((ROOT/'PROTOCOL.json').read_text());protocol['postmechanism_calibration']=dict(max_trials=24,manifest=cfg,justification='Envelope h restores phase/timing but residual added adaptation and readout scaling still distort source fidelity. Fit shared positive smaller adaptation, RG rate normalization, and near-unity PF conversion; no force/neuronal source-weight change.',objective='Sum standing+driven chemicalfractionRMSE +0.2standing sampledMN-RMSE +sumRG relativeperioderrors +2RGphaseerror +2maxstandingherror +2maxhphaseerror',physical_gate_unchanged=True);C.X.write(ROOT/'PROTOCOL.json',protocol)
        text=(ROOT/'EXPERIMENT_RECORD.md').read_text(encoding='utf-8');text+='Frozen comparison completed: soma-h control exactly reproduces all saved v1 traces/events/somas; it remains latched in the driven fixture. Envelope-h restores bilateral cycles (2.482–2.633 s versus source2.132 s) and h peak phase (0–0.03125 cycles error), but exceeds the period/h/waveform gates. Therefore zero physical trials.\n\nA bounded 24-trial Rustuna calibration is now preregistered in `'+str(ROOT/'RUSTUNA_MANIFEST.json')+'`, using the same weights, roles, clock, decoder, and guards. Adaptation0.005–0.4 nA/spike (log), RG C-rate2–8, PF rate scale0.9–1.1, all mirrored. No new physical conditions. The objective explicitly includes h errors. Both source fixtures and h/state gates must pass before any body trial.\n';(ROOT/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8')
