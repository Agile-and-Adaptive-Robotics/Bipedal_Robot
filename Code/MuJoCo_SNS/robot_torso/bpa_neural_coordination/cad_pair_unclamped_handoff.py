"""Measured literal-SNS pair result; qualified comparison and editable graph."""
import json,hashlib,shutil,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cad_pair_unclamped as U
OUT=U.OUT;HERE=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def main():
 proof=read(OUT/'FRESH_REPEAT_AUDIT.json');assert proof['live_neural_continuation_both_exact'];pressure=OUT/'fresh/mixed';reference=OUT/'fresh/reference';r=read(pressure/'UNCLAMPED_SUMMARY.json');ref=read(reference/'UNCLAMPED_SUMMARY.json');assert r['failure'] is None and ref['failure'] is None
 with np.load(pressure/'ACTUATOR_BRIDGE.npz') as z:a={k:z[k] for k in z.files}
 with np.load(reference/'ACTUATOR_BRIDGE.npz') as z:b={k:z[k] for k in z.files}
 with np.load(pressure/'PAIR_NEURAL_LOG.npz') as z:n=z['samples'][~z['sample_is_warm']];nt=z['physical_time'][~z['sample_is_warm']]
 with np.load(reference/'PAIR_NEURAL_LOG.npz') as z:rn=z['samples'][~z['sample_is_warm']];rnt=z['physical_time'][~z['sample_is_warm']]
 with np.load(pressure/'PRESYNAPTIC_TRACE.npz') as z:pre=z['voltage'][100:];pt=z['time'][100:]
 with np.load(pressure/'cad.npz') as z:ground=z['contact_loads'].sum(axis=(1,2))
 m=U.C.mujoco.MjModel.from_xml_path(str(U.C.ROBOT/'ah001_neural.xml'));ids=[m.actuator(s).id for s in ('vas_med_r','vas_med_l')];weight=m.body_mass.sum()*9.81;t=a['t'];mask=t>=.5;F=a['force'][:,ids];target=np.clip(n[:,:,1]/5,0,1)*4530.;work=np.trapz(-F*a['velocity'][:,ids],t,axis=0)
 metrics=dict(machine='easteregg2',duration_s=r['duration_s'],pair_neural_clamp_current_max_nA=r['actual_POST_current_max_nA'],SNS_equation_max_error_mV=r['equation_error_max_mV'],shared_absolute_conductances_uS=U.weights(),physical_MN_voltage_range_mV=r['physical_pair_voltage_range_mV'],normalized_command_range=(np.array(r['physical_pair_voltage_range_mV'])/5).tolist(),pair_activation_floor_imposed=False,other52_floor=.05,pair_force_range_N=r['pair_force_range_N'],pair_mean_force_after_half_second_N=F[mask].mean(axis=0).tolist(),actual_pressure_peak_kpa=r['pressure_peak_kpa'],upward_rig_bw=r['mean_upward_rig_bw'],lateral_rig_rms_bw=r['lateral_rig_rms_bw'],rig_mean_absolute_torques_Nm=r['mean_abs_rig_torque_Nm'],ground_mean_bw=float(ground.mean()/weight),ground_contact_fraction_100Hz_above1N=float(np.mean(ground>1)),knee_range_deg=np.column_stack((np.rad2deg(a['knee']).min(axis=0),np.rad2deg(a['knee']).max(axis=0))).tolist(),actual_mesh_y_range_m=proof['actual_pressure_model_foot_mesh_y_range_m'],pair_mechanical_work_J=work.tolist(),matched_reference=dict(physical_MN_voltage_range_mV=ref['physical_pair_voltage_range_mV'],pair_force_range_N=ref['pair_force_range_N'],upward_rig_bw=ref['mean_upward_rig_bw']),fresh_reference_and_pressure_exact=True,live_SNS_restore_exact=True,whole_runner_restart=False,uncertainty_or_long_standing_tested=False,scope='TwoSNScontrolled assumedBPAs,52hybridreference routes, strongapparatus; knownnominal2s transient, not wholebodyautonomy/gait/hardwarecalibration.')
 (OUT/'RESULTS.json').write_text(json.dumps(metrics,indent=2),encoding='utf-8')
 # Correct inherited description fields only, retaining the original text.
 corrections=[]
 for folder in (OUT/'reference',OUT/'mixed',reference,pressure):
  path=folder/'cad.json';d=read(path);old=d.get('controller_scope');new='Selected2MN clamp currents disabled; six shared recruitment weights normalized. Other52retain original online torqueallocation/MNcurrentclamp.'
  if old!=new:
   d['legacy_adapter_controller_scope']=old;d['controller_scope']=new;d['legacy_base_controller_description']=d.get('controller');d['controller']='LayeredSNS with two free motorneurons;52other routes retain explicit hybrid allocator/currentclamp';d['pair_recruitment_protocol']=str(OUT/'RECRUITMENT_PREREGISTRATION.json');path.write_text(json.dumps(d,indent=2));corrections.append(dict(path=str(path),old_controller_scope=old,new_controller_scope=new,numerical_fields_changed=False))
 if corrections:(OUT/'METADATA_CORRECTIONS.json').write_text(json.dumps(corrections,indent=2))
 plt.rcParams.update({'font.family':'Arial','font.size':10,'legend.fontsize':10,'pdf.fonttype':42,'svg.fonttype':'none'});fig,axes=plt.subplots(4,2,figsize=(9.3,10),sharex=True,constrained_layout=True)
 for side,color in enumerate(('#4477AA','#EE7733')):
  ax=axes[0,side];ax.plot(pt,pre[:,side*3],color=color,label='PF E1');ax.plot(pt,pre[:,side*3+1],color='#777777',ls='--',label='PF E2');ax.plot(pt,pre[:,side*3+2],color='#AA3377',ls=':',label='RG/Ib integrator');ax.set_title(f'({chr(65+side)}) '+('Right' if side==0 else 'Left')+' neural presynaptic activity');ax.legend(loc='lower right')
  ax=axes[1,side];ax.plot(nt,n[:,side,1],color=color,label='Pressure case');ax.plot(rnt,rn[:,side,1],color='#777777',ls='--',label='Matched reference case');ax.axhline(.25,color='black',ls=':',lw=1,label='Offline teaching value');ax.set_title(f'({chr(67+side)}) Free MN, applied clamp current = 0');ax.legend(loc='upper right')
  ax=axes[2,side];ax.plot(t,F[:,side],color=color,label='Actual BPA tension');ax.plot(nt,target[:,side],color=color,ls=':',label='Neural force reference');ax.plot(b['t'],b['force'][:,ids[side]],color='#777777',ls='--',label='Matched reference actuator');ax.set_title(f'({chr(69+side)}) Neural command and actual tension');ax.legend(loc='lower right')
  ax=axes[3,side];ax.plot(t,np.rad2deg(a['knee'][:,side]),color=color,label='Pressure case');ax.plot(b['t'],np.rad2deg(b['knee'][:,side]),color='#777777',ls='--',label='Reference case');ax.axhline(0,color='black',ls=':',lw=1);ax.set_title(f'({chr(71+side)}) Knee flexion, negative');ax.set_xlabel('Time (s)');ax.legend(loc='lower right')
 for row,label in enumerate(('SNS ΔV (mV)','MN ΔV (mV)','Tension (N)','Knee angle (deg)')):
  axes[row,0].set_ylabel(label)
  for ax in axes[row]:ax.set_xlim(0,2);ax.grid(alpha=.2)
 fig.suptitle('Two free SNS motor neurons drive assumed pressure BPAs\n52 other routes remain hybrid; approximately 92% BW apparatus support',fontsize=12)
 for ext in ('png','pdf','svg'):fig.savefig(OUT/('UNCLAMPED_PAIR_NEURAL_MECHANICAL.'+ext),dpi=200)
 plt.close(fig)
 caption='Known nominal2s experiment. Only the selected vas_med motorneurons lose their engineeringPOSTcurrentclamp; they start atrest and thenfollow literalSNS equations. SixexistingPF/RG-Ibintegrator outputweights are scaled bilaterally; intrinsicbias/leak,POSTURE,and nativeafferents remain. Offline0.25mV is a teaching reference, not an applied target or activationfloor. The reference-force comparator has the SAMEfree-neuron weights/controller; its30ms activation-filtered actuators replace only the pressurepair interface. Neuralforcereference isMNvoltage/5mV times4530N masterreference, not measuredBPAcapacity or actualtension. Pressurestate/series-equilibrium dynamics are explicit, assumed, and unchanged. Bothcases repeat everyoriginalsavedarray, and a liveSNS snapshot-step-restore-repeat is exact. Other52routes keep onlineallocation/MNclamps; strongapparatus and−2.5mm visualfoot penetration remain. This is not a spiking, walking, independentbalance, calibratedhardware, pressure-error-relay, uncertainty-validation, or endurance result.'
 caption+=' Voltages use the builder’s relative SNS convention, not absolute physiological membrane potentials. PF activity is tonic in this standing condition; no gait cycles are claimed.'
 (OUT/'UNCLAMPED_PAIR_NEURAL_MECHANICAL_caption.txt').write_text(caption,encoding='utf-8');(OUT/'UNCLAMPED_PAIR_NEURAL_MECHANICAL_alt.txt').write_text('Eight aligned bilateral panels show actual PF and RG/Ib presynaptic activity, unclamped motor voltage, force references versus measured modeledBPA tension, and knee flexion. Reference-actuator curves differ from pressure-actuator curves; applied pair clamp current is zero.',encoding='utf-8')
 # Editable actual graph, checked against the saved numerical connectivity.
 graph=read(U.A.OLD/'ACTUAL_CORE_CONNECTOME.json');names=[r['name'] for r in graph['nodes']]
 with np.load(pressure/'PAIR_NEURAL_STATE.npz') as z:g=z['core_g_max_non']
 changes=[]
 for edge in graph['edges']:
  post=edge['target'];preid=edge['source'];value=float(g[post,preid])
  if value!=edge['gmax_uS']:changes.append(dict(source=names[preid],target=names[post],old=edge['gmax_uS'],new=value));edge['gmax_uS']=value
 assert len(changes)==6 and np.count_nonzero(g)==len(graph['edges'])
 for port in graph['inputs']:
  if port['name'] in ('POST_vas_med_r','POST_vas_med_l'):port.update(runtime_current_forced_zero=True,role='Original engineeringclamp transport disabled at runtime; existinginputcolumn retainedfor defaultparity')
 graph.update(scope='Actualunclamped selectedpair:262neurons/650synapses/230ports. OnlytwoMNrecruitments areliteralSNScontrolled;other52retainhybridcontext.',six_normalized_edges=changes,applied_pair_POST_current_nA=0.,parent_protocol=str(OUT/'RECRUITMENT_PREREGISTRATION.json'));graph['external_blocks']=[dict(name='Other52hybridcontext',role='Originaltorqueallocator and compensatingMNcurrents, explicitly not removed'),dict(name='Selectedpair',role='OriginalPOSTcurrentdiscarded; restedinitialization; subsequentV fromSNS equation; no activationfloor'),dict(name='Pressureinterface',role='MN/5mV→4530N reference→staticinverse→4msdelay→40/60mspressure→seriesforce'),dict(name='AdditiveRG/Ibintegrator',role='SourceIBEXC isnot amultiplicativeloadgate; itsoutputto pair calibratedwithPF'),dict(name='Apparatus',role='Measuredindependentforces/torques;~92%BW verticalsupport')];(OUT/'ACTUAL_UNCLAMPED_CONNECTOME.json').write_text(json.dumps(graph,indent=2))
 record=f'''# CAD pair controlled by literal SNS currents

Machine: **easteregg2**. Both selected motor neurons now drive the pressure-actuated pair without applied engineering clamp current. The known nominal2s pressure trial and its matched reference-force trial complete with no domain failure or knee hyperextension. Fresh processes reproduce every original saved array; live SNS state restoration is exact. This is a two-muscle neural recruitment result inside a strongly supported hybrid model, not autonomous whole-body control.

## What changed

Only six existing output conductances change, identically on both sides: PF_E1→MN = {U.weights()['PF_E1']:.12f} µS; PF_E2→MN = {U.weights()['PF_E2']:.12f} µS; IBEXC_knee_ext→MN = {U.weights()['IBEXC_knee_ext']:.12f} µS. One shared recruitment normalization preserves their original ratios. POSTURE, leak, intrinsicbias, other sensory weights, and the262-neuron/650-synapse core remain unchanged. The IBEXC source combines RG and Ib additively; it is not described as a multiplicative load gate.

Direct POST currents for these two MNs are discarded before entering the SNS; measured appliedcurrent is exactly0nA. Original target-based pair initialization is replaced by one rest initialization, followed by ordinary neural warmup. Thereafter membrane voltages follow the unchanged SNS differential equations, with maximum numericalcheck error {r['equation_error_max_mV']:.3g}mV. The old helper still calculates unused clamp values for compatibility; PAIR_NEURAL_LOG.npz records them separately from the zero appliedcurrent. No actual voltage, pressure, or force is projected to a target.

The original5% normalized command is an offline teaching value. It is **not** imposed as a floor on the freed pair: actual normalized commands range {metrics['normalized_command_range']}. The other52routes retain the existing5% floor and explicit online torqueallocation/MNcurrentclamp. Thus they and the external apparatus can compensate for pair-induced body motion; no whole-SNS necessity or independentbalance claim is made.

## Matched physical comparison

The all-reference comparator uses the same free-pair neural weights and inputs/policy, with30ms activation-filtered reference tension instead of pressure actuation. Its actual pair voltage ranges are {ref['physical_pair_voltage_range_mV']}mV. Pressure-case ranges are {r['physical_pair_voltage_range_mV']}mV. Their differences reflect actual coupled neural/mechanical feedback, not separate gain fits. This observation does not isolate feedback necessity or demonstrate pressure-model-error compensation.

Actual pressure-case force peaks are250.01/261.85N; mean forces after0.5s are {np.round(metrics['pair_mean_force_after_half_second_N'],3).tolist()}N. Peakpressure is {r['pressure_peak_kpa']:.3f}kPa under the unchanged620kPa ceiling. Pair mechanicalwork is {np.round(work,5).tolist()}J, not pneumatic/metabolic energy. Pressuremodel and reference scaling remain declared assumptions: single20mm sleeves, source-derived rest/fixed remote path,200kN/m seriesstiffness,40/60ms pressure response,4ms delay, no calibratedhardware assignment. Other52actuators remain reference forces.

Mean upward apparatusforce is {100*r['mean_upward_rig_bw']:.3f}%BW; lateral RMS is {100*r['lateral_rig_rms_bw']:.3f}%BW. Mean absolute apparatus torques are {np.round(r['mean_abs_rig_torque_Nm'],3).tolist()}Nm. Meanground load is {100*metrics['ground_mean_bw']:.3f}%BW, with contact above1N in {100*metrics['ground_contact_fraction_100Hz_above1N']:.1f}% of100Hz samples. Actualfoot-mesh minimum is {1000*metrics['actual_mesh_y_range_m'][0]:.3f}mm at2ms sampling. Knee ranges are {np.round(metrics['knee_range_deg'],3).tolist()}degrees, all below zero. This2s result includes startup and posture recovery/drift; longer stable standing, reduced assistance, locomotion, and uncertainty robustness remain untested.

## Verification and reproduction

OFFLINE_AND_STATE_GATE.json proves default-off SNS outputs/core arrays are bit-exact, enabled continuation restores all serialized neural states exactly, and frozen pressure replay remains feasible. The live fresh audit repeats snapshot→step→restore→repeat inside both physical contexts; actual original trajectories remain exact. It does not implement a complete restart API for the original CAD runner's local scheduling variables.

```powershell
& 'D:\\Anaconda\\envs\\myo\\python.exe' '{HERE / 'cad_pair_unclamped_followup.py'}'
```

The command reuses completed trials and checks hashes/arrays. Direct --mode commands are not resumable and should not overwrite retained evidence. RECRUITMENT_PREREGISTRATION.json fixes this known-nominal experiment; there was no post-test tuning or independentvalidation. The earlier PF-only and conditionaldiagnosis negatives remain in the parent auditfolder, and the hybrid5% baseline is untouched. METADATA_CORRECTIONS.json transparently records corrected controller-description fields in this branch's outputJSONs; no numerical fields changed. Source XML, original neural/robot runners, accepted Gait checkpoints, and manuscripts are untouched.

ACTUAL_UNCLAMPED_CONNECTOME.json is editable actual connectivity with sixweights checked against saved numericalarrays; two runtime-disabled POST ports and the other52hybrid blocks are explicit. UNCLAMPED_PAIR_NEURAL_MECHANICAL.png/.pdf/.svg plots the real recorded neural/force/knee activity. Full commands/configurations and pressure/force/rig/state evidence are retained locally with compact output.
'''
 record+='''\nTo run a new named study without overwriting retained evidence, use the guarded entry point below. It checks the frozen source/model hashes and requires a new empty outputfolder. The fixed duration is2s; --mode both repeats the matched reference and pressure cases sequentially.\n\n```powershell\n& 'D:\\Anaconda\\envs\\myo\\python.exe' 'D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\robot_torso\\bpa_neural_coordination\\run_saved_cad_pair.py' --output 'D:\\Users\\Ben Bolen\\NeuroCampaign_20261031\\cad_bpa\\my_new_free_pair_run' --mode both\n```\n\nVoltage traces use relative SNS units, not absolute physiological membrane potentials.\n'''
 for old,new in {'intrinsicbias':'intrinsic bias','appliedcurrent':'applied current','torqueallocation':'torque allocation','MNcurrentclamp':'MN current clamp','Meanground':'Mean ground','Actualfoot':'Actual foot','Peakpressure':'Peak pressure','seriesstiffness':'series stiffness','physicalwork':'physical work','neuralforcereference':'neural force reference','whole-SNS':'whole SNS','independentbalance':'independent balance','outputJSONs':'output JSONs','source XML':'source XML','%BW':'% BW'}.items():record=record.replace(old,new)
 lines=[]
 for line in record.splitlines():
  if 'D:\\' not in line and not line.startswith('& '):
   line=re.sub(r'(?<=[A-Za-z])(?=\d)',' ',line);line=re.sub(r'(?<=\d)(?!(?:e|E)[+-]?\d)(?=[A-Za-z])',' ',line)
  lines.append(line)
 record=('\n'.join(lines)+'\n').replace('easteregg 2','easteregg2').replace('Gait 2392','Gait2392')
 (OUT/'EXPERIMENT_RECORD.md').write_text(record,encoding='utf-8');snap=OUT/'source_snapshot';manifest={}
 for p in (Path(__file__),Path(U.__file__),HERE/'cad_pair_unclamped_followup.py',HERE/'run_saved_cad_pair.py',Path(U.C.__file__),Path(U.L.__file__),U.C.ROBOT/'run_neural_robot.py'):
  digest=hashlib.sha256(p.read_bytes()).hexdigest();dest=snap/(p.stem+'_'+digest[:8]+p.suffix)
  if not dest.exists():shutil.copyfile(p,dest)
  manifest[str(p)]=dict(sha256=digest,snapshot=str(dest))
 (snap/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(metrics,indent=2))
if __name__=='__main__':main()
