"""Causal event-marked inverse decoder, independently gated before network use."""
import json,math,inspect,hashlib,argparse,copy,time
from pathlib import Path
import numpy as np
from numba import njit
import gait_central_inverse as Q
C=Q.C;M=Q.M;ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_dynamic_readout_v4')
TAU=.025;N=16

@njit(cache=True)
def support_weight(rate):
    # Exponential rate-estimator variance gives an effective history of2*tau.
    # Poisson>=2 support is a declared uncertainty model, not fabricated spikes.
    count=max(0.,rate)*N*(2*TAU)
    if count>60:return 1.
    return max(0.,min(1.,-math.expm1(-count)-count*math.exp(-count)))

@njit(cache=True)
def decode(rate,g_mark,h_mark,g_current,h_current,cm,gna,ka,sa,ea,ena,mode):
    if rate<.1:return -4.
    g=g_mark/rate if mode in (1,3) else g_current;h=h_mark/rate if mode in (1,3) else h_current
    raw=Q.inverse_nap_rate(rate,cm,max(.01,g),gna,min(1.,max(0.,h)),ka,sa,ea,ena)
    confidence=support_weight(rate) if mode in (2,3) else 1.
    return -4.+confidence*(raw+4.)

@njit(cache=True)
def transient_sim(g,h,current,initial,mask_events=False):
    dt=.0002;cm=.05;gna=12.;rate=0.;gm=0.;hm=0.;v=np.linspace(2.5,3.9,N);ref=np.zeros(N);teacher=initial;rows=np.empty((len(g)//10,14));counts=0;observed=0
    for step in range(len(g)):
        # Independent fourth-order reference ODE. It is never read by decoder.
        k1=Q.nap_vector(teacher,current[step],g[step],gna,h[step],1.,.8,2.,8.)/cm
        k2=Q.nap_vector(teacher+.5*dt*k1,current[step],g[step],gna,h[step],1.,.8,2.,8.)/cm
        k3=Q.nap_vector(teacher+.5*dt*k2,current[step],g[step],gna,h[step],1.,.8,2.,8.)/cm
        k4=Q.nap_vector(teacher+dt*k3,current[step],g[step],gna,h[step],1.,.8,2.,8.)/cm
        teacher+=dt*(k1+2*k2+2*k3+k4)/6.;spikes=0
        for cell in range(N):
            if ref[cell]>.0001:ref[cell]=max(0.,ref[cell]-dt);v[cell]=2.5;continue
            ion=gna*h[step]*Q.activation(v[cell],1.,.8,2.);total=g[step]+ion;steady=(current[step]+ion*8)/total
            v[cell]=steady+(v[cell]-steady)*math.exp(-total*dt/cm)
            if v[cell]>=4*(1+.02*(cell/(N-1)-.5)):v[cell]=2.5;ref[cell]=.001;spikes+=1
        used=0 if mask_events else spikes;decay=math.exp(-dt/TAU)
        rate=rate*decay+used/(N*TAU);gm=gm*decay+used*g[step]/(N*TAU);hm=hm*decay+used*h[step]/(N*TAU);counts+=spikes;observed+=used
        if step%10==9:
            j=step//10;rows[j,0]=(step+1)*dt;rows[j,1]=teacher;rows[j,2]=rate;rows[j,3]=g[step];rows[j,4]=h[step];rows[j,5]=current[step]
            for mode in range(4):rows[j,6+mode]=decode(rate,gm,hm,g[step],h[step],cm,gna,1.,.8,2.,8.,mode)
            rows[j,10]=counts;rows[j,11]=observed;rows[j,12]=support_weight(rate);rows[j,13]=v.mean();counts=0;observed=0
    return rows

def validate_transients():
    ROOT.mkdir(parents=True,exist_ok=True)
    old=json.loads((Q.ROOT/'READOUT_VALIDATION.json').read_text());assert old['derived_transfer_valid']
    digest=hashlib.sha256((inspect.getsource(Q.inverse_nap_rate.py_func)+inspect.getsource(Q.passage_time.py_func)+inspect.getsource(Q.nap_vector.py_func)).encode()).hexdigest();assert digest==old['derivation_sha256']
    protocol=dict(machine='easteregg2',source_inverse=str(Q.ROOT),source_derivation_sha256=digest,carrier='Only actual AP events update exponentialrate and event-marked channel-history states. Current/teacherVm have no decoderargument. No invented spikes or gradedshadow.',derivation='For r=sum expkernel/N, matchedG=sum G_at_AP*expkernel/N dividedbyr; matchedh analogous. Both histories use exactly the same25mskernel as rate, avoiding mixing past spikes with instantaneous channels. Two-event support weight1−(1+lambda)exp(−lambda),lambda=2*N*tau*r;2tau is equivalent countwindow of exponential estimator. Weight shrinks low-support inverse toward silent−4mV.',assumption='Poisson two-event support is an explicit regularization assumption, not the deterministic population spiking law or a claim of Bayesian posterior correctness.',clock_s=.0002,shared_parameters=dict(N=N,tau=TAU),modes=['v3instantaneous','event_marked_channels_only','support_shrinkage_only','combined_preregistered_candidate'],reference='IndependentRK4 nonspiking NaP equation with prescribedphysicalh/G/current; referenceVm never entersdecoder. Exact same observed actualspikes feedallfourmodes.',heldout='InhibitoryG step at h.05/.2; excitatoryG step; current-onlypulse; Gramp; h ramp. BaselineG1.75,Vm4.56;changes0.8–1.2s. No parameters fitted.',gates=dict(score_start_s=.4,observable_chemical_RMSE_max=.15,observable_hInf_RMSE_max=.15,inhibitory_off_lag_s_max=.08,zero_event_output_exact_minus4=True,settled_silence_s=.2),nonidentifiable_control='Stable positiveVm3.2 belowAPthreshold produces no spikes; decodedsilence expected. This explicitly documented unidentifiable domain is not declared voltage-fidelity passing.',acceptance='Combinedcandidate mustpassALLobservabletransients pluszero/silence controls before fullnetwork. No gainsearch unless independenttransients valid; originalstanding/rhythm gates unchanged. No bodybydefault.',script_sha256=C.X.digest(__file__))
    C.X.write(ROOT/'TRANSIENT_PROTOCOL.json',protocol)
    t=np.arange(15000)*.0002;cases={};arrays={};pulse=(t>=.8)&(t<1.2)
    for name,h0 in [('inhibitory_h005',.05),('inhibitory_h020',.2),('inhibitory_h080_strong',.8),('excitatory_step',.05),('current_pulse',.05),('conductance_ramp',.05),('channel_h_ramp',.05),('unidentifiable_positive',.05),('zero_events',.05)]:
        target=3.2 if name=='unidentifiable_positive' else 4.56;g=np.full(len(t),1.75);h=np.full(len(t),h0);base=1.75*target-12*h0*Q.activation(target,1.,.8,2.)*(8-target);current=np.full(len(t),base)
        if name.startswith('inhibitory'):
            added=5. if name.endswith('strong') else 3.;g[pulse]+=added;current[pulse]-=5*added
        elif name=='excitatory_step':g[pulse]+=1.5;current[pulse]+=12.
        elif name=='current_pulse':current[pulse]+=2.
        elif name=='conductance_ramp':ramp=np.clip((t-.8)/.2,0,1)*np.clip((1.4-t)/.2,0,1);g+=3*ramp;current-=15*ramp
        elif name=='channel_h_ramp':h+=.15*np.clip((t-.8)/.2,0,1)*np.clip((1.4-t)/.2,0,1)
        a=transient_sim(g,h,current,target,mask_events=name=='zero_events');arrays[name]=a;mask=a[:,0]>=.4;rows={}
        for mode in range(4):
            voltage=a[:,6+mode];teacher=a[:,1];chemical=np.clip(voltage/5,0,1);truth=np.clip(teacher/5,0,1);hi=lambda v:1/(1+np.exp(np.clip(2*(v-3.5),-700,700)))
            r=dict(chemical_RMSE=float(np.sqrt(np.mean((chemical[mask]-truth[mask])**2))),hInf_RMSE=float(np.sqrt(np.mean((hi(voltage[mask])-hi(teacher[mask]))**2))),finite=bool(np.isfinite(voltage).all()))
            if name.startswith('inhibitory'):
                interval=np.flatnonzero((a[:,0]>=.8)&(a[:,0]<1.2));expected=interval[np.flatnonzero(truth[interval]<.1)[0]];off=interval[np.flatnonzero(chemical[interval]<.1)[0]] if np.any(chemical[interval]<.1) else interval[-1];r['off_lag_s']=max(0.,float(a[off,0]-a[expected,0]))
            r['passed']=r['finite'] and r['chemical_RMSE']<=.15 and r['hInf_RMSE']<=.15 and r.get('off_lag_s',0)<=.08;rows[str(mode)]=r
        cases[name]=dict(modes=rows,actual_events=int(a[:,10].sum()),observed_events=int(a[:,11].sum()),observable=name not in ('zero_events','unidentifiable_positive'))
    zero=bool(np.all(arrays['zero_events'][:,6:10]==-4.));quiet=bool(arrays['unidentifiable_positive'][:,10].sum()==0 and np.all(arrays['unidentifiable_positive'][:,6:10]==-4.));passed=all(c['modes']['3']['passed'] for c in cases.values() if c['observable']) and zero and quiet
    report=dict(independent_transients_valid=bool(passed),cases=cases,zero_observed_events_silent=zero,unidentifiable_positive_reports_silence=quiet,physical_trials=0,optimization_trials=0,scope='Prescribed-channel scalarcomponentvalidation, not fullnetwork/physical acceptance.')
    C.X.write(ROOT/'TRANSIENT_VALIDATION.json',report);np.savez_compressed(ROOT/'TRANSIENT_TRACES.npz',**arrays,columns=['time_s','teacher_mV','filtered_rate_Hz','gLeakSyn_uS','actual_h','effective_current_nA','instant_inverse_mV','matched_only_mV','support_only_mV','combined_mV','actual_events_2ms','observed_events_2ms','support_weight','mean_spiking_soma_mV'])
    (ROOT/'EXPERIMENT_RECORD.md').write_text('# Dynamic spike-derived RG readout\n\nMachine: easteregg2. Independent transient validation only, no gain search or body trial. Event-marked conductance/h histories share the actual spike-rate kernel; a declared two-event support regularizer handles sparse-event ambiguity. Teacher voltage and inputcurrent neverenter the decoder. Actualspikes are preserved, with an explicit observed-event ablation control.\n\nCombined observable-transient gate: '+str(passed)+'. Full results: `'+str(ROOT/'TRANSIENT_VALIDATION.json')+'`. Stable positive subthreshold voltage without spikes is explicitly unidentifiable and not relabelled voltage fidelity. V1/v2/v3 unchanged.\n',encoding='utf-8')
    print(json.dumps(dict(independent_transients_valid=passed,zero_event_silence=zero,unidentifiable_control=quiet,combined={k:v['modes']['3'] for k,v in cases.items() if v['observable']})),flush=True);return report

@njit(cache=True)
def dynamic_population_step(lv,h,ref,adapt,decoded,counts,event_g,event_h,params,gs,ge,applied,config,dt):
    decay=math.exp(-dt/TAU)
    for p in range(len(params)):
        cm,gm,rest,bias,gna,ka,sa,ea,kb,sb,eb,th,ena,is_rg=params[p];reset=2.5 if is_rg else 0.;norm=config[1]/cm;spikes=0;spike_h=0.
        envelope=decode(decoded[p]*norm,event_g[p],event_h[p],gm+gs[p],np.mean(h[p]),cm,gna,ka,sa,ea,ena,3) if is_rg else 0.
        for cell in range(lv.shape[1]):
            old=lv[p,cell];adapt[p,cell]*=math.exp(-dt/.35)
            if is_rg:
                hi=1/(1+kb*np.exp(sb*(eb-envelope)));h[p,cell]+=(hi-h[p,cell])*dt/th;ion=gna*Q.activation(old,ka,sa,ea)*h[p,cell]
            else:ion=0.
            if ref[p,cell]>dt/2:ref[p,cell]=max(0.,ref[p,cell]-dt);lv[p,cell]=reset;continue
            g=gm+gs[p]+ion;steady=(gm*rest+bias+applied[p]+ge[p]+ion*ena-adapt[p,cell])/g
            lv[p,cell]=steady+(old-steady)*math.exp(-g*dt/cm);threshold=(4. if is_rg else 1.)*(1+(.02 if is_rg else .15)*(cell/(lv.shape[1]-1)-.5))
            if lv[p,cell]>=threshold:
                lv[p,cell]=reset;ref[p,cell]=.001;spikes+=1;spike_h+=h[p,cell]
                if is_rg:adapt[p,cell]+=config[0]
        decoded[p]=decoded[p]*decay+spikes/(lv.shape[1]*TAU)*(1/norm if is_rg else config[2]);counts[p]+=spikes
        if is_rg:event_g[p]=event_g[p]*decay+spikes*(gm+gs[p])/(lv.shape[1]*TAU);event_h[p]=event_h[p]*decay+spike_h/(lv.shape[1]*TAU)

class DynamicCentral(Q.InverseCentral):
    def __init__(self,*args,**kwargs):super().__init__(*args,**kwargs);self.event_g=np.zeros(16);self.event_h=np.zeros(16)
    def proxy_voltage(self):
        voltage=C.CentralForward.proxy_voltage(self)
        for p,i in enumerate(self.indices):
            if not self.is_rg[p]:continue
            cm,gm,rest,bias,gna,ka,sa,ea,kb,sb,eb,th,ena,_=self.params[p];voltage[p]=decode(self.decoded[p]*self.config[1]/cm,self.event_g[p],self.event_h[p],gm+self.inner.sums[i,0],float(self.h[p].mean()),cm,gna,ka,sa,ea,ena,3)
        return voltage
    def forward(self,x=None):
        if not self.enabled:return self.inner.forward(x)
        c=self.inner.c;c.V[self.indices]=self.proxy_voltage();self.inner.forward(x)
        if self.poison_shadow:c.V[self.indices]=1e6
        self.actual_input[:]=self.inner.i_app[self.indices];self.counts[:]=0;gs=self.inner.sums[self.indices,0];ge=self.inner.sums[self.indices,1]
        for _ in range(10):dynamic_population_step(self.lv,self.h,self.ref,self.adapt,self.decoded,self.counts,self.event_g,self.event_h,self.params,gs,ge,self.actual_input,self.config,.0002)
        c.V[self.indices]=self.proxy_voltage();self.total_counts+=self.counts;self.steps+=1
        if self.inner.output_onehot:np.take(c.V,self.inner.oq,out=self.inner.output_values);c.outputs=self.inner.output_values
        else:c.outputs=c.output_voltage_connectivity@c.V
        return c.outputs
    def export_state(self):
        state=C.CentralForward.export_state(self);state.update(inactivation_voltage_policy=np.array('event_marked_spike_history_support_v4'),observer_event_g=self.event_g.copy(),observer_event_h=self.event_h.copy(),observer_constants=np.array([N,TAU,2]),source_population_params=self.params.copy(),source_bias=self.inner.c.i_b.copy());return state
    def import_state(self,state):
        assert str(state['inactivation_voltage_policy'])=='event_marked_spike_history_support_v4' and np.array_equal(state['observer_constants'],[N,TAU,2]);C.CentralForward.import_state(self,state)
        for attr,key in ((self.params,'source_population_params'),(self.inner.c.i_b,'source_bias'),(self.event_g,'observer_event_g'),(self.event_h,'observer_event_h')):np.copyto(attr,state[key])

def prepare_network():
    validation=json.loads((ROOT/'TRANSIENT_VALIDATION.json').read_text());assert validation['independent_transients_valid'],'No fullnetwork before independenttransient pass'
    protocol=dict(machine='easteregg2',configuration=C.DEFAULT,observer='Causal event-marked leak/synapticG andNaPh with matched25msAPhistory, plus preregistered two-event supportshrinkage; alloutputs and h-envelope depend on realspikes, neveroriginalgradedVm/currenttarget.',unchanged='v3cellmechanics,NaP,reset/refractory/adaptation,weights,PFdecoder,200usclock,emptyspikehistory,sharedL/Rparameters; no fitting.',source_v3=str(Q.ROOT),source_v3_script_sha256=C.X.digest(Q.__file__),transient_validation_sha256=C.X.digest(ROOT/'TRANSIENT_VALIDATION.json'),acceptance='Originalstanding+rhythm and h gates unchanged; complete electrical audits required; no bodybydefault. Independent scalartransients pass is not fullnetwork acceptance.',script_sha256=C.X.digest(__file__))
    C.X.write(ROOT/'NETWORK_PROTOCOL.json',protocol)
    inner,names=C.X.construct();M.STATE.update(source=inner,names=names,inputs={},reference={},graded_h={})
    source=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait2392_reflex_demo_v2\cases\syn6_graded_baseline.inputs.npz')
    with np.load(source) as z:M.STATE['inputs']['standing']=z['inputs'].copy()
    u=np.zeros((5000,M.STATE['inputs']['standing'].shape[1]));u[:,np.flatnonzero(inner.c.input_connectivity[names.index('DRIVE')])[0]]=2.5;M.STATE['inputs']['deaff']=u
    with np.load(M.SOURCE/'GRADED_REFERENCE.npz') as z:M.STATE['reference']={k:dict(trace=z[k].copy()) for k in ('standing','deaff')}
    with np.load(M.ROOT/'MEANFIELD_FROZEN.npz') as z:M.STATE['graded_h']={k:z[k+'_graded_h'].copy() for k in ('standing','deaff')}
    a=C.clone(inner);b=DynamicCentral(C.clone(inner),names,enabled=False)
    for u in M.STATE['inputs']['standing'][:64]:assert np.array_equal(a.forward(u),b.forward(u))
    a=DynamicCentral(C.clone(inner),names);b=DynamicCentral(C.clone(inner),names,poison_shadow=True)
    for u in M.STATE['inputs']['standing'][:64]:assert np.array_equal(a.forward(u),b.forward(u))
    start=a.export_state()
    for u in M.STATE['inputs']['standing'][64:96]:a.forward(u)
    first=a.export_state();a.import_state(start)
    for u in M.STATE['inputs']['standing'][64:96]:a.forward(u)
    second=a.export_state();checks={k:bool(np.array_equal(v,second[k],equal_nan=True)) if v.dtype.kind in 'fc' else bool(np.array_equal(v,second[k])) for k,v in first.items()};assert all(checks.values())
    rejected=False
    try:a.import_state(Q.InverseCentral(C.clone(inner),names).export_state())
    except (AssertionError,KeyError):rejected=True
    assert rejected
    arrays=[a.lv.copy(),a.h.copy(),a.ref.copy(),a.adapt.copy(),a.decoded.copy(),a.counts.copy(),a.event_g.copy(),a.event_h.copy()];other=[v.copy() for v in arrays];gs=a.inner.sums[a.indices,0];ge=a.inner.sums[a.indices,1]
    dynamic_population_step(*arrays,a.params,gs,ge,a.actual_input,a.config,.0002);dynamic_population_step.py_func(*other,a.params,gs,ge,a.actual_input,a.config,.0002)
    errors={k:float(np.max(abs(x-y))) for k,x,y in zip(('soma','h','ref','adapt','spikefilter','counts','eventG','eventh'),arrays,other)};assert max(errors.values())<1e-11
    C.X.write(ROOT/'STATE_AUDIT.json',dict(disabled_default64_exact=True,shadow_poison64_exact=True,complete_continuation32_exact=True,incompatible_v3_state_rejected=True,compiled_scalar_200us_errors=errors,arrays=checks,physical_trials=0));print(json.dumps(dict(prepared=True,electrical_audits_passed=True)),flush=True)

def compare_network():
    prepare_network();results={}
    for label,klass in (('V3_CONTROL',Q.InverseCentral),('DYNAMIC_FROZEN',DynamicCentral)):
        rows={};data={};start=time.perf_counter();old=M.MeanFieldCentral;M.MeanFieldCentral=klass
        try:
            for key,u in M.STATE['inputs'].items():
                a,_=M.replay(C.clone(M.STATE['source']),M.STATE['names'],u,C.DEFAULT,fixture=key=='deaff');metrics=C.metrics(a,M.STATE['reference'][key],key)
                if key=='standing':metrics['RG_h_mean_absolute_error']=np.mean(abs(a['h_mean'][500:,[0,1,8,9]]-M.STATE['graded_h'][key][500:,[0,1,8,9]]),axis=0).tolist()
                else:
                    hp=M.h_profile(a['trace'],a['h_mean']);hs=M.h_profile(M.STATE['reference'][key]['trace'],M.STATE['graded_h'][key]);metrics['h_peak_cycle_phase_error']=(np.full(4,.5) if hp is None else abs((np.argmax(hp,axis=1)/64-np.argmax(hs,axis=1)/64+.5)%1-.5)).tolist()
                data[key]=a;rows[key]=metrics
        finally:M.MeanFieldCentral=old
        if label=='V3_CONTROL':
            with np.load(Q.ROOT/'INVERSE_FROZEN.npz') as z:
                for key,a in data.items():
                    for f in ('trace','events','somas','h_mean'):assert np.array_equal(a[f],z[key+'_'+f]),(key,f)
        passed=all(r['passed'] for r in rows.values()) and max(rows['standing']['RG_h_mean_absolute_error'])<=.15 and max(rows['deaff']['h_peak_cycle_phase_error'])<=.15
        result=dict(configuration=C.DEFAULT,passed_offline_gate=passed,fixtures=rows,wall_s=time.perf_counter()-start,physical_trials=0,optimization_trials=0)
        C.X.write(ROOT/(label+'.json'),result);np.savez_compressed(ROOT/(label+'.npz'),**{key+'_'+f:a[f] for key,a in data.items() for f in ('trace','events','somas','h_mean')},**{k+'_graded_h':v for k,v in M.STATE['graded_h'].items()},dt_s=.002,configuration_json=json.dumps(C.DEFAULT));results[label]=result;print(json.dumps(dict(case=label,**result)),flush=True)
    text=(ROOT/'EXPERIMENT_RECORD.md').read_text(encoding='utf-8');text+='\nFrozen fullnetwork comparison completed. V3control exactly matches all saved traces/events/somas/h. Combinedobserver gate: '+str(results['DYNAMIC_FROZEN']['passed_offline_gate'])+'. No search or bodytrial. See DYNAMIC_FROZEN.json and STATE_AUDIT.json.\n';(ROOT/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--validate-transients',action='store_true');p.add_argument('--compare-network',action='store_true');a=p.parse_args()
    if a.validate_transients:validate_transients()
    if a.compare_network:compare_network()
