"""Fixed telemetry-only inverse decoder diagnosis; no fit or body integration."""
import sys,json,numpy as np
from pathlib import Path
from scipy.signal import find_peaks
import gait_central_inverse as Q
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
R=Q.ROOT;M=Q.M;C=Q.C
COLS=['time_s','spikeHz','gLeakSyn_uS','ge_nA','J_actual_nA','J_from_spikes_nA','sampled_totalG_uS','sampled_ion_nA','readout_mV','soma_mean_mV','h','adapt_nA']
def load(p):
    with np.load(p) as z:return {k:z[k].copy() for k in z.files}
class TelemetryInverse(Q.InverseCentral):
    def __init__(self,*args,**kwargs):super().__init__(*args,**kwargs);self.telemetry=[]
    def forward(self,x=None):
        ans=super().forward(x);rows=[]
        for p in (0,1,8,9):
            i=self.indices[p];cm,gm,rest,bias,gna,ka,sa,ea,kb,sb,eb,th,ena,_=self.params[p];h=float(self.h[p].mean());g=gm+self.inner.sums[i,0];ge=self.inner.sums[i,1];rate=self.decoded[p]*self.config[1]/cm;v=self.inner.c.V[i];adapt=float(self.adapt[p].mean());actual=gm*rest+bias+self.actual_input[p]+ge-adapt
            inferred=g*v-gna*h*Q.activation(v,ka,sa,ea)*(ena-v) if rate>=.1 else np.nan
            activation=1/(1+ka*np.exp(sa*(ea-self.lv[p])));gt=g+gna*h*activation.mean();ion=gna*h*np.mean(activation*(ena-self.lv[p]))
            rows.append([(self.steps-1)*.002,rate,g,ge,actual,inferred,gt,ion,v,float(self.lv[p].mean()),h,adapt])
        self.telemetry.append(rows);return ans

if '--plot-only' not in sys.argv:
    protocol=(R/'COMPARISON_PROTOCOL.json').read_bytes();Q.prepare_comparison();(R/'COMPARISON_PROTOCOL.json').write_bytes(protocol)
    probe={};original=M.MeanFieldCentral;M.MeanFieldCentral=TelemetryInverse
    try:
        for key,u in M.STATE['inputs'].items():
            a,c=M.replay(C.clone(M.STATE['source']),M.STATE['names'],u,C.DEFAULT,fixture=key=='deaff')
            with np.load(R/'INVERSE_FROZEN.npz') as z:
                for f in ('trace','events','somas','h_mean'):assert np.array_equal(a[f],z[key+'_'+f]),(key,f)
            probe[key]=np.array(c.telemetry)
    finally:M.MeanFieldCentral=original
    np.savez_compressed(R/'INVERSE_TELEMETRY.npz',**probe,columns=COLS,RG_names=['RG_E_l','RG_F_l','RG_E_r','RG_F_r'])
z=load(R/'INVERSE_TELEMETRY.npz');new=load(R/'INVERSE_FROZEN.npz');old=load(R/'V2_CONTROL.npz');ref=load(M.SOURCE/'GRADED_REFERENCE.npz');report={}
def profile(a):
    v=np.clip(a[:,:16]/5,0,1);p,_=find_peaks(v[1000:,1],prominence=.3,distance=100);p=p+1000;cycles=[]
    for l,r in zip(p[:-1],p[1:]):cycles.append(np.array([np.interp(np.linspace(l,r,64,endpoint=False),np.arange(len(v)),v[:,i]) for i in range(16)]))
    return np.mean(cycles,axis=0) if cycles else None
for key in ('standing','deaff'):
    a=z[key];mask=a[:,0,0]>= (2. if key=='deaff' else 1.);metrics={}
    for k,name in enumerate(('RG_E_l','RG_F_l','RG_E_r','RG_F_r')):
        row=a[:,k];valid=mask&np.isfinite(row[:,5]);dg=np.abs(np.r_[0.,np.diff(row[:,2])]);err=np.abs(row[:,5]-row[:,4]);transition=valid&(dg>.01);quiet=valid&(dg<.0001)
        events=new[key+'_events'][:,(0,1,8,9)[k]];no_recent=np.convolve(events,np.ones(25,dtype=int),mode='full')[:len(events)]==0
        metrics[name]=dict(mean_abs_drive_reconstruction_error_nA=float(np.mean(err[valid])) if valid.any() else None,drive_reconstruction_error_p95_nA=float(np.quantile(err[valid],.95)) if valid.any() else None,mean_error_slow_conductance_nA=float(np.mean(err[quiet])) if quiet.any() else None,mean_error_fast_conductance_nA=float(np.mean(err[transition])) if transition.any() else None,spike_free50ms_but_positive_readout_duty=float(np.mean((no_recent&(row[:,8]>0))[mask])),filtered_rate_range_Hz=[float(row[mask,1].min()),float(row[mask,1].max())])
    report[key]=metrics
sourcep=profile(ref['deaff']);newp=profile(new['deaff_trace']);oldp=profile(old['deaff_trace'])
report['cycle_registered_chemical_RMSE']=dict(inverse=None if newp is None else float(np.sqrt(np.mean((newp-sourcep)**2))),v2_control=float(np.sqrt(np.mean((oldp-sourcep)**2))),scope='Diagnostic only; raw acceptance unchanged.')
def availability_profile(a,h):
    v=np.clip(a[:,:16]/5,0,1);p,_=find_peaks(h[1000:,0],prominence=.15,distance=250);p=p+1000
    cycles=[np.array([np.interp(np.linspace(l,r,64,endpoint=False),np.arange(len(v)),v[:,i]) for i in range(16)]) for l,r in zip(p[:-1],p[1:])]
    return (np.mean(cycles,axis=0) if cycles else None),p
sourcehp,sourcepeaks=availability_profile(ref['deaff'],new['deaff_graded_h']);predhp,predpeaks=availability_profile(new['deaff_trace'],new['deaff_h_mean'])
eventrate=np.convolve(new['deaff_events'][:,0]/(.002*16),np.ones(50)/50,mode='same');eventpeaks,_=find_peaks(eventrate[1000:],prominence=10.,distance=300);eventpeaks=eventpeaks+1000
report['separate_slowburst_diagnostics']=dict(source_h_period_s=float(np.median(np.diff(sourcepeaks))*.002) if len(sourcepeaks)>2 else None,inverse_h_period_s=float(np.median(np.diff(predpeaks))*.002) if len(predpeaks)>2 else None,actual_events100ms_density_peak_period_s=float(np.median(np.diff(eventpeaks))*.002) if len(eventpeaks)>2 else None,availability_registered_waveform_RMSE=float(np.sqrt(np.mean((sourcehp-predhp)**2))) if predhp is not None else None,scope='Posthoc diagnostic detector: hprominence0.15,minspacing0.5s; actual eventdensity100ms,prominence10Hz,minspacing0.6s. Retains original readout-peak gates; does not promote intrinsic cadence or waveform acceptance.')
validation=json.loads((R/'READOUT_VALIDATION.json').read_text());fresh=json.loads((R/'INVERSE_FROZEN.json').read_text())
checkpoint=dict(machine='easteregg2',derived_transfer_validation=validation,nominal_inverse=fresh,telemetry_full_array_parity=True,telemetry_diagnostics=report,physical_trials=0,optimization_trials=0,scope='Static-derived readout improves standing electrical fidelity but fails dynamic waveform/cadence gates. No body or fit promotion.')
C.X.write(R/'CENTRAL_CHECKPOINT.json',checkpoint)
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'svg.fonttype':'none','pdf.fonttype':42})
fig,axes=plt.subplots(6,2,figsize=(6.5,9.3),sharex='col');handles=[]
for j,key in enumerate(('standing','deaff')):
    t=np.arange(1,len(new[key+'_trace'])+1)*.002
    for row,col in ((0,0),(1,1)):
        for label,a,color,style in (('V2 linear',old[key+'_trace'],'#d55e00','--'),('Inverse',new[key+'_trace'],'#0072b2','-'),('Graded source',ref[key],'#555555',':')):
            line,=axes[row,j].plot(t,a[:,col],color=color,ls=style,lw=1.,label=label)
            if row==0 and j==0:handles.append(line)
    axes[2,j].plot(t,new[key+'_h_mean'][:,0],color='#0072b2');axes[2,j].plot(t,new[key+'_graded_h'][:,0],color='#555555',ls=':')
    a=z[key][:,0];axes[3,j].plot(t,a[:,1],color='#0072b2');density=np.convolve(new[key+'_events'][:,0]/(.002*16),np.ones(50)/50,mode='same');densityline,=axes[3,j].plot(t,density,color='#222222',ls=':',lw=.9);axes[4,j].plot(t,a[:,2],color='#0072b2')
    axes[5,j].plot(t,a[:,4],color='#009e73',lw=1.,label='Actual net drive');axes[5,j].plot(t,a[:,5],color='#cc79a7',lw=1.,ls='--',label='From spike rate')
    axes[0,j].set_title('Recorded standing\nFrozen parameters' if key=='standing' else 'Deafferented drive\nFrozen parameters',fontsize=10);axes[0,j].set_xlim((1.,3.) if key=='standing' else (2.,9.8));axes[5,j].set_xlabel('Neural time (s)')
    for ax in axes[:,j]:ax.grid(alpha=.18)
for row,label in enumerate(['RG-E readout\n(mV)','RG-F readout\n(mV)','RG-E h\n(availability)','RG-E spike\nrate(Hz)','RG-E leak+syn\n(µS)','RG-E net drive\n(nA)']):axes[row,0].set_ylabel(label)
fig.subplots_adjust(left=.15,right=.98,top=.92,bottom=.24,wspace=.32,hspace=.21)
drivehandles=axes[5,0].get_lines();fig.legend(handles+list(drivehandles)+[densityline],['V2 linear','Inverse','Graded source','Actual net drive','From spike rate','100ms event density'],loc='lower center',bbox_to_anchor=(.5,.125),ncol=3,frameon=False)
fig.text(.035,.025,'Static transfer passes held-out drives; dynamic source gates still fail.\nCurrent is recorded for validation only; it never supplies a presynaptic voltage.\nNaP output uses actual spikes plus channel/conductance state; all MNs stay graded.\n256 central spiking cells; no gain search or body trial.',fontsize=10)
for suffix in ('pdf','svg','png'):fig.savefig(R/('INVERSE_READOUT_DIAGNOSTIC.'+suffix),dpi=220,bbox_inches='tight',facecolor='white')
plt.close(fig)
text='# Conductance-aware RG spike readout\n\nMachine: easteregg2. **Derived static transfer passes, but frozen dynamic source gates fail. Zero fitting or body trials.**\n\n'
text+='For each actual firing population, observed spike rate and known leak/synaptic/NaP channel state determine an effective current by inverting its interspike-time integral; the stable NaP equilibrium then provides the readout. Runtime accepts no original graded voltage or input-current target. The actual current is logged only for independent validation. No-spike subthreshold positive voltage remains unidentifiable and maps silent−4 mV.\n\n'
text+='Independent tests include the analytical passive limit, SciPy-integrated nonlinear NaP firing curves, and 36 held-out actual200µs drives repeated for the16-cell threshold-heterogeneous population. Maximum continuous roundtrip error8.41e−9 mV; population chemical/h∞ errors0.002276/0.002197. The midpoint threshold-crossing clock correction is explicitly0.1 ms; integration clock remains0.2 ms. This validates fixedh/steady-current firing transfer, not arbitrary dynamics.\n\n'
text+='Frozen inverse greatly improves standing: chemical-fraction RMSE0.01866, sampled MN RMSE0.000306 mV, extensor-h error0.04535 versus v2error0.33633. Dynamic raw chemical RMSE0.58484 and peak-counter periods0.26–0.504 s still fail the source2.132 s gate. These peak counts are not promoted to successful burst timing. Phase-normalized diagnostic: '+str(report['cycle_registered_chemical_RMSE'])+'.\n\n'
text+='A telemetry-only fixed replay reproduces all stored traces/events/somas/h exactly. It logs actual conductance, filtered spike rate, known net drive, and rate-reconstructed drive, to distinguish a static transfer from a valid dynamic observer. Full diagnostics: `'+str(R/'CENTRAL_CHECKPOINT.json')+'`. Positive readout after a50ms spike-free interval is explicitly measured as retained filter history, not new spikes or a hidden graded signal. Additional slow-availability/event-density diagnostics, without changing original gates: '+str(report['separate_slowburst_diagnostics'])+'.\n\n'
text+='The inverse filter starts with empty spike history, preserving initial chemical silence; its initial raw proxy−4 differs from the old linear-codec source0. Actual soma/h initial states, weights, default shared parameters,200µs clock,reset/refractory/adaptation, and PF decoder stay fixed. Both h and central chemical output use the new spike-derived envelope.\n\n'
text+='Checks: disabled/default64-step equality, complete-state32-step continuation, discarded-graded-voltage poison64-step equality, incompatible-policy snapshot rejection, and local compiled/scalar parity pass. V2control full trace/event/soma/h arrays reproduce exactly. V1/v2 remain untouched.\n\n'
text+='Files and commands on easteregg2:\n\n'
for name in ('READOUT_PROTOCOL.json','READOUT_VALIDATION.json','COMPARISON_PROTOCOL.json','STATE_AUDIT.json','INVERSE_FROZEN.json','CENTRAL_CHECKPOINT.json','INVERSE_READOUT_DIAGNOSTIC.pdf'):
    text+='- `'+str(R/name)+'`\n'
text+='\nNumerical program: `'+str(Path(Q.__file__))+'`, with `D:\\Anaconda\\envs\\myo\\python.exe` and `--validate --compare`. This repeats only the frozen analytical/neural comparisons, not fitting or body integration. Telemetry/figure program: `'+str(Path(__file__))+'`; `--plot-only` uses saved telemetry without neural replay.\n'
(R/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8');print(json.dumps(dict(physical_trials=0,optimization_trials=0,phase_registered=report['cycle_registered_chemical_RMSE'],slowburst=report['separate_slowburst_diagnostics'],standing_drive=report['standing']['RG_E_l'],dynamic_drive=report['deaff']['RG_E_l'])),flush=True)
