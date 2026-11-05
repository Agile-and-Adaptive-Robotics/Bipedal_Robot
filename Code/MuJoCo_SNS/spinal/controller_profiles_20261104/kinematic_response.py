"""Measured joint response versus IK cycles; harmonic ratios are not a Bode FRF.

One common right-knee peak anchor preserves hip/knee/ankle and interleg phase.
Cycle normalization is for shape only; actual stride period is scored separately.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from dominant_cycles import measure

JOINTS=('hip_flexion_r','knee_angle_r','ankle_angle_r',
        'hip_flexion_l','knee_angle_l','ankle_angle_l')
GRID=np.arange(200)/200.

def describe(values):
 a=np.asarray(values,dtype=float)
 if not a.size:return dict(n=0,mean=None,SD=None,minimum=None,maximum=None)
 return dict(n=int(a.size),mean=float(a.mean()),SD=float(a.std(ddof=1))if a.size>1 else None,minimum=float(a.min()),maximum=float(a.max()))

def phase_stats(values):
 a=(np.asarray(values,dtype=float)+180.)%360.-180.
 if not a.size:return dict(n=0,circular_mean_deg=None,circular_SD_deg=None,mean_absolute_error_deg=None,maximum_absolute_error_deg=None)
 vector=np.mean(np.exp(1j*np.radians(a)));R=abs(vector)
 return dict(n=int(a.size),circular_mean_deg=float(np.degrees(np.angle(vector))),circular_SD_deg=float(np.degrees(np.sqrt(-2*np.log(np.clip(R,1e-15,1.))))),mean_absolute_error_deg=float(np.mean(abs(a))),maximum_absolute_error_deg=float(np.max(abs(a))))

def reference(path=None):
 if path is None:
  import kine_ref
  ref=kine_ref.load_reference()
  sources=[kine_ref.IK_MOT,kine_ref.GRF_MOT]
  identity='subject01_walk1: existing calibration IK/GRF reference'
 else:
  path=Path(path);z=np.load(path,allow_pickle=False);ref={}
  for side in('r','l'):
   ref[side]={joint:np.asarray(z[side+'_'+joint])for joint in('hip','knee','ankle')}
   ref['T_'+side]=float(z['T_'+side])
  ref['lag_rl']=float(z['lag_rl']);sources=[path];identity=path.name
 if not all(np.isfinite(ref['T_'+s])and ref['T_'+s]>0 for s in('r','l')):raise ValueError('Invalid reference period')
 for side in('r','l'):
  for joint in('hip','knee','ankle'):
   if not np.isfinite(ref[side][joint]).all():raise ValueError('Nonfinite reference')
  if np.ptp(ref[side]['knee'])<5.:raise ValueError('Degenerate reference knee cycle')
 # Stored left/right cycles each start at their OWN loading onset. Restore
 # their relative phase before applying a shared kinematic anchor.
 out=[]
 for name in JOINTS:
  side=name[-1];joint=name.split('_')[0];y=np.asarray(ref[side][joint])
  p=np.linspace(0.,1.,len(y));query=(GRID-(ref['lag_rl']if side=='l'else 0.))%1.
  out.append(np.interp(query,p,y))
 raw=np.array(out).T;anchor=int(np.argmin(raw[:,1]))/len(GRID)
 aligned=np.column_stack([np.interp((GRID+anchor)%1.,GRID,raw[:,i],period=1.)for i in range(6)])
 return dict(identity=identity,period_s=[ref['T_r'],ref['T_l']],cycle_deg=aligned,
             source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in sources},
             right_knee_anchor_fraction=anchor,reference_left_loading_phase=float(ref['lag_rl']))

def harmonics(y):
 # Complex cosine convention: y=mean+Re(C_k exp(i*2*pi*k*phase)).
 return np.array([2*np.mean((y-y.mean(axis=0))*np.exp(-2j*np.pi*k*GRID[:,None]),axis=0)for k in(1,2,3)])

def stage_policy(stage,safety,motion,cadence_ok,ground_ok,kinematic_loss,air_terms):
 air=stage in(1,2)
 passed=bool(safety and motion and(air or(cadence_ok and ground_ok)))
 loss=5.*air_terms['bilateral_phase_error']+air_terms['period_cv']+air_terms['normalized_cycle_variability']if air else kinematic_loss
 return passed,float(loss+(0. if passed else 1000.))

def compare(cycles,periods,ref):
 mean=np.mean(cycles,axis=0);target=ref['cycle_deg'];scale=np.maximum(np.ptp(target,axis=0),5.)
 rms=np.sqrt(np.mean((mean-target)**2,axis=0));bias=np.mean(mean-target,axis=0)
 amplitude=np.ptp(mean,axis=0);target_amplitude=np.ptp(target,axis=0)
 gains=amplitude/np.maximum(target_amplitude,1e-9);sim_h=harmonics(mean);ref_h=harmonics(target)
 rows=[];phase_cost=[];cycle_h=np.stack([harmonics(cycle)for cycle in cycles])
 for i,name in enumerate(JOINTS):
  hs=[]
  for k in range(3):
   usable=abs(ref_h[k,i])>=1. # Ignore unstable ratios of near-zero reference harmonics.
   ratio=sim_h[k,i]/ref_h[k,i]if usable else None
   phase_usable=usable and abs(sim_h[k,i])>=.1
   lag=-np.angle(ratio,deg=True)if phase_usable else None
   if k==0 and usable:phase_cost.append(abs(lag)/180. if phase_usable else 1.)
   hs.append(dict(order=k+1,reference_amplitude_deg=float(abs(ref_h[k,i])),simulation_amplitude_deg=float(abs(sim_h[k,i])),usable=bool(usable),gain=float(abs(ratio))if usable else None,gain_dB=float(20*np.log10(max(abs(ratio),1e-12)))if usable else None,lag_deg=float(lag)if usable else None))
  errors=cycles[:,:,i]-target[None,:,i];per_cycle_amplitude=np.ptp(cycles[:,:,i],axis=1)
  per_cycle_max=np.max(abs(errors),axis=1);worst=np.unravel_index(np.argmax(abs(errors)),errors.shape)
  usable_phase=(abs(ref_h[0,i])>=1.)&(abs(cycle_h[:,0,i])>=.1)
  phase_errors=-np.angle(cycle_h[usable_phase,0,i]/ref_h[0,i],deg=True)if abs(ref_h[0,i])>=1. else np.array([])
  rows.append(dict(joint=name,reference_ROM_deg=float(target_amplitude[i]),reference_minimum_deg=float(target[:,i].min()),reference_maximum_deg=float(target[:,i].max()),simulation_ROM_deg=float(amplitude[i]),ROM_gain=float(gains[i]),cycle_amplitude_deg=describe(per_cycle_amplitude),cycle_amplitude_error_deg=describe(per_cycle_amplitude-target_amplitude[i]),cycle_amplitude_gain=describe(per_cycle_amplitude/max(target_amplitude[i],1e-9)),cycle_amplitude_percent_error=describe(100*(per_cycle_amplitude/max(target_amplitude[i],1e-9)-1.)),fundamental_phase_error=phase_stats(phase_errors),mean_bias_deg=float(bias[i]),mean_absolute_error_deg=float(np.mean(abs(errors))),maximum_absolute_error_deg=float(abs(errors).max()),mean_cycle_maximum_absolute_error_deg=float(abs(mean[:,i]-target[:,i]).max()),cycle_maximum_absolute_error_deg=describe(per_cycle_max),worst_error_cycle_index=int(worst[0]),worst_error_cycle_phase_percent=float(100*GRID[worst[1]]),waveform_RMSE_deg=float(rms[i]),normalized_RMSE=float(rms[i]/scale[i]),cycle_variability_RMS_deg=float(np.sqrt(np.mean((cycles[:,:,i]-mean[:,i])**2))),harmonics=hs))
 cadence=float(np.mean(abs(np.log(np.asarray(periods)/ref['period_s']))))
 cycle_phase=[row['fundamental_phase_error']['mean_absolute_error_deg']/180. if row['fundamental_phase_error']['n']else 1. for row in rows if row['harmonics'][0]['usable']]
 terms=dict(cadence_log_error=cadence,waveform_normalized_RMSE=float(np.mean(rms/scale)),ROM_log_gain_error=float(np.mean(abs(np.log(np.maximum(gains,1e-6))))),mean_bias_normalized=float(np.mean(abs(bias)/scale)),fundamental_phase_error=float(np.mean(cycle_phase))if cycle_phase else 1.,worst_cycle_max_error_normalized=float(np.mean([row['maximum_absolute_error_deg']/scale[i]for i,row in enumerate(rows)])))
 weights=dict(cadence_log_error=2.,waveform_normalized_RMSE=1.,ROM_log_gain_error=.5,mean_bias_normalized=.5,fundamental_phase_error=.5,worst_cycle_max_error_normalized=.5)
 return dict(joints=rows,terms=terms,weights=weights,kinematic_loss=float(sum(terms[k]*weights[k]for k in terms))),mean

def score_recording(folder,ref=None,start=5.,end=14.):
 folder=Path(folder);z=np.load(folder/'capture/aligned_trace.npz');body=json.loads((folder/'body.result.json').read_text());t=z['t_state'];keep=(t>=start)&(t<end)
 names=z['joint_names'].tolist();q=np.column_stack([np.degrees(z['qpos'][:,int(z['joint_qpos_addresses'][names.index(j)])])for j in JOINTS]);tt=t[keep];qq=q[keep]
 principal=measure(tt,-qq[:,[1,4]],qq[:,[0,3]])if len(tt)>2 and np.isfinite(qq).all()else dict(accepted=False,reason='Missing or nonfinite samples')
 report=dict(schema='stage_specific_kinematics_v2',recording=str(folder),reference=ref['identity']if ref else None,scope='IK waveform response comparison, not input-output transfer function or Bode FRF',window_s=[start,end],phase_anchor='One shared right-knee maximum-flexion peak; all six joints retain relative phase. Left reference cycles shifted by reference loading phase.',principal=principal)
 report['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in(Path(__file__),HERE/'dominant_cycles.py',folder/'capture/aligned_trace.npz')}
 peaks=principal.get('principal_peak_times_s',[[]])[0];cycles=[];lengths=[]
 for a,b in zip(peaks[:-1],peaks[1:]):
  if .8<=b-a<=2.5:
   cycles.append(np.column_stack([np.interp(a+GRID*(b-a),tt,qq[:,i])for i in range(6)]));lengths.append(b-a)
 gate=body['assessment'];safety=all(gate.get(key,False)for key in('finite','complete','frontal_gate','knee_extension_gate'))
 report.update(safety_gates_passed=bool(safety),complete_cycles=len(cycles))
 if len(cycles)<2:
  report.update(valid=False,accepted=False,objective=10000.,reason='Fewer than two complete right-knee-anchored cycles; no favorable frozen/fast score')
  return report,None,None
 ref=reference()if ref is None else ref;report['reference']=ref['identity'];report['reference_source_sha256']=ref['source_sha256']
 periods=[]
 for p in principal['principal_peak_times_s']:periods.append(float(np.median(np.diff(p)))if len(p)>1 else 0.)
 result,mean=compare(np.asarray(cycles),periods,ref);report.update(result)
 from scipy.signal import periodogram
 for i,row in enumerate(report['joints']):
  row.update(window_minimum_deg=float(qq[:,i].min()),window_maximum_deg=float(qq[:,i].max()),window_total_amplitude_deg=float(np.ptp(qq[:,i])))
  f,power=periodogram(qq[:,i],fs=1./np.median(np.diff(tt)),window='hann');df=f[1]-f[0];total=power[1:].sum();high=f>3.;fraction=float(power[high].sum()/max(total,1e-30))
  row['frequency_content']=dict(dominant_spectral_Hz=float(f[1:][np.argmax(power[1:])])if total>1e-15 else None,frequency_bin_width_Hz=float(df),power_above3Hz_fraction=fraction,high_frequency_RMS_deg=float(np.sqrt(power[high].sum()*df)),strongest_peak_above3Hz=float(f[high][np.argmax(power[high])])if np.any(high)and fraction>1e-6 else None,scope='Hann periodogram of unwarped window; 3Hz split is diagnostic, not a universal physiological bound')
 frequency=[]
 for side,peaks_side,Tref in zip(('r','l'),principal['principal_peak_times_s'],ref['period_s']):
  intervals=np.diff(peaks_side);hz=1./intervals
  frequency.append(dict(side=side,period_s=describe(intervals),frequency_Hz=describe(hz),reference_period_s=float(Tref),reference_frequency_Hz=float(1./Tref),frequency_error_Hz=describe(hz-1./Tref),frequency_percent_error=describe(100*(hz*Tref-1.))))
 report['frequency_statistics']=frequency
 report['statistics_scope']='Descriptive within-run cycle statistics; SD is sample SD, not independent-subject uncertainty. Maximum waveform error includes every retained cycle, using one shared right-knee anchor; cadence error remains separate.'
 report['statistics_revision']='amplitude_frequency_cycle_phase_maxerror_v1'
 cadence_ok=all(abs(np.asarray(periods)/ref['period_s']-1.)<=.15)
 air=body['stage']in(1,2)
 # Ground IK is a diagnostic covariate in air, not an asserted air target.
 # Rank air coordination by measured phase stability and repeatability only.
 variability=float(np.mean([row['cycle_variability_RMS_deg']/max(row['simulation_ROM_deg'],5.)for row in result['joints']]))
 air_terms=dict(bilateral_phase_error=float(principal['bilateral_phase_error']),period_cv=float(principal['period_cv']),normalized_cycle_variability=variability)
 passed,objective=stage_policy(body['stage'],safety,principal['accepted'],cadence_ok,gate['accepted'],result['kinematic_loss'],air_terms)
 report.update(valid=True,accepted=passed,simulation_periods_s=periods,reference_periods_s=ref['period_s'],reference_cadence_within15percent=bool(cadence_ok),reference_cadence_gate_applied=not air,comparison_role='diagnostic_only_air'if air else'ground_reference_target',air_coordination_terms=air_terms if air else None,objective=objective,acceptance_scope='Air: existing motion/safety gates, no grounded IK gate or loss. Ground: existing ground assessment plus reference cadence; waveform loss ranks, no validated waveform cutoff yet.')
 return report,mean,ref['cycle_deg']

def self_check():
 phase=GRID[:,None];base=np.array([20.,30.,12.,20.,30.,12.]);y=base*np.cos(2*np.pi*phase+np.arange(6)[None,:]*.3)
 ref=dict(cycle_deg=y,period_s=[1.,1.]);perfect,_=compare(np.stack([y,y]),[1.,1.],ref)
 doubled,_=compare(np.stack([2*y,2*y]),[1.,1.],ref)
 lagged=base*np.cos(2*np.pi*(phase-.125)+np.arange(6)[None,:]*.3);delayed,_=compare(np.stack([lagged,lagged]),[1.,1.],ref)
 slow,_=compare(np.stack([y,y]),[2.,2.],ref)
 assert perfect['kinematic_loss']<1e-12
 assert abs(doubled['joints'][0]['harmonics'][0]['gain']-2.)<1e-12
 assert abs(delayed['joints'][0]['harmonics'][0]['lag_deg']-45.)<1e-10
 assert slow['kinematic_loss']>perfect['kinematic_loss']
 biased,_=compare(np.stack([y+10,y+10]),[1.,1.],ref)
 swapped=y.copy();swapped[:,3:]*=-1.;wrong_side,_=compare(np.stack([swapped,swapped]),[1.,1.],ref)
 assert biased['kinematic_loss']>perfect['kinematic_loss']
 assert wrong_side['kinematic_loss']>perfect['kinematic_loss']
 air_terms=dict(bilateral_phase_error=.05,period_cv=.02,normalized_cycle_variability=.01)
 air_a=stage_policy(2,True,True,False,False,99.,air_terms)
 air_b=stage_policy(2,True,True,True,True,0.,air_terms)
 assert air_a==air_b and air_a[0]
 assert not stage_policy(4,True,True,False,True,0.,air_terms)[0]
 assert not stage_policy(2,False,True,True,True,0.,air_terms)[0]
 assert abs(doubled['joints'][0]['cycle_amplitude_gain']['mean']-2.)<1e-12
 assert abs(delayed['joints'][0]['fundamental_phase_error']['circular_mean_deg']-45.)<1e-10
 assert abs(biased['joints'][0]['maximum_absolute_error_deg']-10.)<1e-12
 spike=y.copy();spike[20,0]+=100.;spiked,_=compare(np.stack([y,spike]),[1.,1.],ref)
 assert abs(spiked['joints'][0]['maximum_absolute_error_deg']-100.)<1e-12
 wrapped=phase_stats([179.,-179.])
 assert abs(abs(wrapped['circular_mean_deg'])-180.)<1e-10 and wrapped['circular_SD_deg']<1.01
 return dict(passed=True,perfect_loss=perfect['kinematic_loss'],double_amplitude_gain=doubled['joints'][0]['harmonics'][0]['gain'],delayed_phase_lag_deg=delayed['joints'][0]['harmonics'][0]['lag_deg'],half_cadence_loss=slow['kinematic_loss'])

def write_results(report,path):
 lines=['# Measured kinematic statistics','',f"Machine: easteregg2. Source: `{report['recording']}`.",'',f"Comparison role: **{report.get('comparison_role','invalid')}**. Grounded IK differences do not reject air stepping.",'',f"Window: {report['window_s']} s; {report['complete_cycles']} complete right-knee-anchored cycles.",'']
 if not report['valid']:
  lines += [report['reason']];path.write_text('\n'.join(lines)+'\n');return
 fmt=lambda v:'undefined'if v is None else f'{v:.3f}'
 stat=lambda d:fmt(d['mean'])+' ± '+fmt(d['SD'])
 lines += ['## Frequency and period','','| Side | Intervals n | Frequency mean ± SD (Hz) | Reference (Hz) | Frequency error (%) | Period mean ± SD (s) |','| --- | ---: | ---: | ---: | ---: | ---: |']
 for row in report['frequency_statistics']:lines.append(f"| {row['side'].upper()} | {row['period_s']['n']} | {stat(row['frequency_Hz'])} | {fmt(row['reference_frequency_Hz'])} | {fmt(row['frequency_percent_error']['mean'])} | {stat(row['period_s'])} |")
 lines += ['','## Total amplitude and extrema','','Amplitude means peak-to-peak excursion, not half-amplitude. Cycle SD is sample SD.','','| Joint | Reference amplitude (deg) | Whole-window amplitude (deg) | Per-cycle amplitude mean ± SD (deg) | Mean amplitude error (%) | Window min / max (deg) |','| --- | ---: | ---: | ---: | ---: | ---: |']
 for row in report['joints']:lines.append(f"| {row['joint']} | {fmt(row['reference_ROM_deg'])} | {fmt(row['window_total_amplitude_deg'])} | {stat(row['cycle_amplitude_deg'])} | {fmt(row['cycle_amplitude_percent_error']['mean'])} | {fmt(row['window_minimum_deg'])} / {fmt(row['window_maximum_deg'])} |")
 lines += ['','## Joint frequency content','','Leg cadence alone can hide fast ankle oscillations. These spectra use the unwarped 5–14s window; bin width is approximately 0.111Hz. The 3Hz split is diagnostic.','','| Joint | Dominant spectral frequency (Hz) | Strongest peak above 3Hz (Hz) | Power above 3Hz (%) | High-frequency amplitude RMS (deg) |','| --- | ---: | ---: | ---: | ---: |']
 for row in report['joints']:
  s=row['frequency_content'];lines.append(f"| {row['joint']} | {fmt(s['dominant_spectral_Hz'])} | {fmt(s['strongest_peak_above3Hz'])} | {fmt(100*s['power_above3Hz_fraction'])} | {fmt(s['high_frequency_RMS_deg'])} |")
 lines += ['','## Phase error and maximum error','','Phase is the fundamental lag versus reference, with one shared knee anchor. Positive means lag; circular statistics handle wraparound. Maximum error uses **every cycle**, not just the averaged waveform.','','| Joint | Phase cycles n | Circular phase mean ± SD (deg) | Maximum absolute phase error (deg) | Worst absolute waveform error (deg) | Per-cycle maximum error mean ± SD (deg) | Mean bias (deg) |','| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
 for row in report['joints']:
  phase=row['fundamental_phase_error'];lines.append(f"| {row['joint']} | {phase['n']} | {fmt(phase['circular_mean_deg'])} ± {fmt(phase['circular_SD_deg'])} | {fmt(phase['maximum_absolute_error_deg'])} | {fmt(row['maximum_absolute_error_deg'])} | {stat(row['cycle_maximum_absolute_error_deg'])} | {fmt(row['mean_bias_deg'])} |")
 lines += ['','## Additional waveform diagnostics','','| Joint | Mean absolute error (deg) | Mean-cycle RMSE (deg) |','| --- | ---: | ---: |']
 for row in report['joints']:lines.append(f"| {row['joint']} | {fmt(row['mean_absolute_error_deg'])} | {fmt(row['waveform_RMSE_deg'])} |")
 lines += ['',report['statistics_scope'],'','These few cycles are correlated observations from one run. They do not justify independent-sample significance tests or population confidence intervals.','',f"Reference: {report['reference']}. See KINEMATIC_RESPONSE.json for all unrounded values, gains, per-side frequency errors, and provenance."]
 path.write_text('\n'.join(lines)+'\n',encoding='utf-8')

def main():
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path);p.add_argument('--output',type=Path);p.add_argument('--reference',type=Path);p.add_argument('--database',type=Path);p.add_argument('--self-check',action='store_true');a=p.parse_args()
 checks=self_check()
 if a.self_check:print(json.dumps(checks));return
 if not a.run or not a.output:p.error('--run and --output required')
 a.output.mkdir(parents=True,exist_ok=True);report,mean,target=score_recording(a.run,reference(a.reference));report['metric_self_check']=checks
 (a.output/'KINEMATIC_RESPONSE.json').write_text(json.dumps(report,indent=2,allow_nan=False))
 write_results(report,a.output/'RESULTS.md')
 if a.database:
  rows=[]
  for path in sorted(a.database.glob('*.npz')):
   try:
    candidate,_,_=score_recording(a.run,reference(path))
    rows.append({key:candidate.get(key)for key in('reference','valid','accepted','reference_source_sha256','simulation_periods_s','reference_periods_s','reference_cadence_within15percent','comparison_role','kinematic_loss','terms','reason')})
   except(ValueError,KeyError)as exc:rows.append(dict(reference=path.name,excluded=True,reason=str(exc)))
  (a.output/'DATABASE_COMPARISON.json').write_text(json.dumps(dict(scope='Same fixed-command saved trajectory compared to every supplied reference; not a speed-adaptation test and no best-reference selection',references=rows),indent=2,allow_nan=False))
 if mean is not None:
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  plt.rcParams.update({'font.family':'Arial','font.size':10})
  fig,axes=plt.subplots(3,2,figsize=(10,7),sharex=True,constrained_layout=True)
  for i,name in enumerate(JOINTS):
   ax=axes[i%3,i//3];ax.plot(100*GRID,target[:,i],'--',color='#222222',label='IK reference');ax.plot(100*GRID,mean[:,i],color='#4477AA',label='Simulation');ax.set_title(name);ax.set_ylabel('OpenSim coordinate (deg)');ax.grid(alpha=.2)
  for ax in axes[-1]:ax.set_xlabel('Right-knee-anchored cycle (%)')
  axes[0,0].legend();fig.suptitle('Cycle shape comparison; actual cadence scored separately\nSimulation periods '+str([round(v,3)for v in report['simulation_periods_s']])+' s; reference '+str([round(v,3)for v in report['reference_periods_s']])+' s')
  for ext in('png','pdf'):fig.savefig(a.output/('IK_OVERLAY.'+ext),dpi=170)
  plt.close(fig)
  fig,axes=plt.subplots(3,2,figsize=(9,7),sharex=True,constrained_layout=True)
  for i,row in enumerate(report['joints']):
   hs=[h for h in row['harmonics']if h['usable']];style='-'if i<3 else'--';label=row['joint'][-1].upper();axes[i%3,0].plot([h['order']for h in hs],[h['gain_dB']for h in hs],style+'o',label=label);axes[i%3,1].plot([h['order']for h in hs],[h['lag_deg']for h in hs],style+'o',label=label)
  for i in range(3):
   axes[i,0].set_ylabel(('Hip','Knee','Ankle')[i]+' gain (dB)');axes[i,1].set_ylabel('Lag (deg)')
   for ax in axes[i]:ax.grid(alpha=.2);ax.legend();ax.set_xticks([1,2,3])
  for ax in axes[-1]:ax.set_xlabel('Stride harmonic order')
  fig.suptitle('Simulation/reference harmonic ratios\nDescriptive gait response, not a Bode transfer function')
  for ext in('png','pdf'):fig.savefig(a.output/('HARMONIC_GAIN_PHASE.'+ext),dpi=170)
  plt.close(fig)
 print(json.dumps(dict(run=str(a.run),valid=report['valid'],accepted=report['accepted'],objective=report['objective'],terms=report.get('terms'))))
if __name__=='__main__':main()
