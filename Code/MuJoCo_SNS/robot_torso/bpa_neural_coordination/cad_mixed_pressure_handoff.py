"""Qualified CAD mixed-actuation results, plots, source bundle, and commands."""
import json,sys,hashlib,shutil,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cad_mixed_coactivation05 as J
BASE=J.OUT;HERE=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def main():
 proof=read(BASE/'FRESH_AND_MESH_AUDIT.json');case=BASE/'fresh_a/mixed';summary=read(case/'BRIDGE_SUMMARY.json');core=read(case/'cad.json');cfg=read(BASE/'ASSUMED_PAIR_CONFIGURATION.json');initial=read(BASE/'mixed/BRIDGE_SUMMARY.json');zero=read(BASE.parent/'mixed_pressure_pair/mixed/BRIDGE_SUMMARY.json');parity=read(BASE/'MOVING_AND_DEFAULT_PARITY.json')
 with np.load(case/'ACTUATOR_BRIDGE.npz') as z:a={k:z[k] for k in z.files}
 with np.load(case/'cad.npz') as z:ground=z['contact_loads'].sum(axis=(1,2));contact_fraction=float(np.mean(ground>1.));min_ground=float(ground.min());max_ground=float(ground.max())
 with np.load(BASE/'reference/ACTUATOR_BRIDGE.npz') as z:reference={k:z[k] for k in z.files}
 names=list(cfg['actuators']);m=J.C.mujoco.MjModel.from_xml_path(str(J.C.ROBOT/'ah001_neural.xml'));ids=[m.actuator(n).id for n in names];time=a['t'];start=time>=.5;weight=m.body_mass.sum()*9.81;force=a['force'][:,ids];mean_force=force[start].mean(axis=0);work=np.trapz(-force*a['velocity'][:,ids],time,axis=0);stats=dict(machine='easteregg2',trial_duration_s=2.,pressure_actuator_names=names,pressure_actuator_count=2,reference_actuator_count=52,base_SNS_cells=262,controller='Hybrid torqueallocator with explicit MN synaptic-current cancellation; not autonomous neural coordination',activation_floor=.05,post_startup_mean_pair_tension_N=mean_force.tolist(),pressure_peak_kpa=summary['pressure_peak_kpa'],pair_mechanical_work_J=work.tolist(),mean_ground_load_bw=float(ground.mean()/weight),ground_contact_fraction_at_100Hz_above_1N=contact_fraction,ground_load_range_N=[min_ground,max_ground],foot_mesh_y_range_m=proof['actual_mesh_min_y_range_m'],mesh_sampling_dt_s=.002,upward_rig_bw=summary['mean_upward_rig_bw'],lateral_rig_rms_bw=summary['lateral_rig_rms_bw'],knee_max_deg=np.rad2deg(summary['knee_max_rad']).tolist(),fresh_physics_matches_original=True,fresh_captured_state_matches=True,complete_restart_api_demonstrated=False,no_all_BPA_or_hardware_claim=True)
 (BASE/'RESULTS.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')
 plt.rcParams.update({'font.family':'Arial','font.size':10,'legend.fontsize':10,'pdf.fonttype':42,'svg.fonttype':'none'});fig,axes=plt.subplots(3,2,figsize=(9,8),sharex=True,constrained_layout=True)
 for side,color in enumerate(('#4477AA','#EE7733')):
  ax=axes[0,side];ax.plot(time,a['pressure'][:,side],color=color,label='Actual pressure');ax.plot(time,a['command'][:,side],color='#777777',ls='--',label='Reference');ax.set_ylim(0,640);ax.set_title(f'({chr(65+side)}) '+('Right' if side==0 else 'Left')+' vastus route, MN command 0.05');ax.legend(loc='upper left')
  ax=axes[1,side];ax.plot(time,force[:,side],color=color,label='Actual BPA tension');ax.plot(reference['t'],reference['force'][:,ids[side]],color='#777777',ls='--',label='Matched reference-force case');ax.set_title(f'({chr(67+side)}) Pressure-to-force response');ax.legend(loc='lower right')
  ax=axes[2,side];ax.plot(time,np.rad2deg(a['knee'][:,side]),color=color,label='Mixed actuation');ax.plot(reference['t'],np.rad2deg(reference['knee'][:,side]),color='#777777',ls='--',label='Reference actuation');ax.axhline(0,color='black',ls=':',lw=1);ax.set_title(f'({chr(69+side)}) Knee flexion, negative');ax.set_xlabel('Time (s)');ax.legend(loc='lower right')
 for row,label in enumerate(('Pressure (kPa)','Tension (N)','Knee angle (deg)')):
  axes[row,0].set_ylabel(label)
  for ax in axes[row]:ax.grid(alpha=.2);ax.set_xlim(0,2)
 fig.suptitle(f'CAD model: 2 assumed pressure BPAs and 52 reference tension actuators\nHybrid SNS and allocator, 5% coactivation, {100*summary["mean_upward_rig_bw"]:.1f}% BW harness support',fontsize=12)
 for ext in ('png','pdf','svg'):fig.savefig(BASE/('CAD_MIXED_PRESSURE_PAIR.'+ext),dpi=200)
 plt.close(fig)
 caption='Short heavily supported engineering trial, not autonomous neural control or an all-BPA robot. Two actual262cell SNS motor outputs, clamped by the original torqueallocator/MNcurrenttarget interface, feed an explicit static-force inverse and dynamic pressure bank. Single20mm sleeves, paired rest length,200kN/m series stiffness,40/60ms fill/vent and4ms delay are assumed, not identified. Other52routes retain reference tension with30ms activation filtering. Both vastus routes have identical5%neural recruitment floor; differing pressures/tensions follow differing geometry and dynamics. Original no-floor mixed trial produced0N in both BPAs and is preserved negative. Two independent repeats match original trajectories and each other’s captured states exactly. Actual visualfoot minimum−2.55mm; strong pelvis apparatus remains, intermittent contact is quantified separately. No stepping, independentbalance, endurance, hardware, or neuralcompensation claim.'
 (BASE/'CAD_MIXED_PRESSURE_PAIR_caption.txt').write_text(caption,encoding='utf-8');(BASE/'CAD_MIXED_PRESSURE_PAIR_alt.txt').write_text('Six bilateral panels show pressure lag, tension rise toward226.5N, and knees staying flexed in a2s heavily supported mixed-actuation model. Solid pressure-model trajectories are compared with dashed reference-force traces.',encoding='utf-8')
 command=f"& 'D:\\Anaconda\\envs\\myo\\python.exe' '{HERE / 'cad_mixed_pressure_followup.py'}'"
 record=f'''# CAD mixed pressure-pair integration

Machine: **easteregg2**. A fixed5% coactivation fixture completes2s with two pressure-actuated vastus routes and52reference-tension routes. Actual peak tensions are238.75N/252.82N; peakpressure499.43kPa. Two fresh processes reproduce the original physical traces and each other’s saved arrays and captured states exactly. This is a heavily supported hybrid-controller interface demonstration, not autonomous neural coordination, a complete BPA robot, or hardware identification.

## Controller and physics

The original262-cell,230-input SNS build in {J.C.ROBOT / 'run_neural_robot.py'} is preserved. Its online bounded torqueallocator and motor_target_current() background-synaptic-current cancellation remain explicit. There is no direct control bypass, but this engineering MN-current clamp prevents a whole-SNS-necessity claim. Neural normalized motor outputs are0.05 in this fixture. The existing5% activation-floor option applies to all54routes, not just the pressure pair; this is distinct from the no-floor baseline.

Pressure routes are vas_med_r/l; source master ID11 labels these Vastus Intermedius. The naming discrepancy is retained. MIF reference scaling is4530N per route, not measured BPA strength. A single20mm sleeve per side has source-derived rest length{cfg['actuators']['vas_med_r']['resting_length_m']:.9f}m. A constant remote-span offset{cfg['actuators']['vas_med_r']['constant_remote_path_offset_m']:.9f}m is added to original route length; its derivative is zero because the idealized remote sleeve/housing endpoints share the thigh body. Fixed nonmuscle length is chosen to give zero unloaded contraction at the source401-pose original-route maximum. This is an ideal lossless transmission assumption, not measured Bowden routing. Effective seriesstiffness200kN/m,40ms fill/60ms vent,4ms delay,25% maximum contraction, and620kPa ceiling are declared assumptions. No sleeve geometry was resized after a failure. Friction, hysteresis, valve identification, fittings, added mass, and continuous packaging remain unverified.

All native gains and biases are explicitly disabled in external-force modes. The complete moving-route force −JᵀT is added once; the native moving correction is not added again. Ia/II use complete dynamic length/velocity; Ib uses actual external bank/reference tensions, never native zero actuator_force. The other52routes use their existing30ms activation-filtered reference tensions. Existing support is preserved.

## Gates, failures, and comparisons

The original new0.1rad knee-flexed pose requires416.37N/295.32N merely to keep the assumed sleeves in domain at zero pressure. Its cold failure is preserved under mixed_pressure_pair/INITIAL_FLEXED_POSE_COLD_FAILURE.json. The successful setup uses existing source-certified frame1069, with knee coordinates−0.27790/−0.22834rad and zero initialvelocity. It is not a geometry or coordinate-sign repair. Earlier static PREREGISTRATION.json retained stale0.1rad text; LIVE_PREREGISTRATION.json explicitly corrects this metadata and identifies the actual source pose.

At the fixed sourcepose, static inverse/pressure equilibrium gives226.5N per route at371.575/363.741kPa. Reference-to-mixed reduced-torque change is roundoff-level. Complete54-route moving-Jacobian force parity has maximumerror{parity['full_moving_generalized_force_max_error']:.3g}Nm and zero virtual-power error; the disabled interface preserves original model/force arrays. This is a fixed-pose default gate, not a full-duration default-trajectory repeat.

The no-floor mixed2s case finishes but produces **0N through both pressure routes**. It is retained as a negative actuation result rather than counted as a BPA standing success. The single fixed5% fixture was registered before its tests; no sweep or uncertainty validation was run. Its matched reference-only baseline passes before the mixed test.

## Measured limits

Mean upward harnessforce is{100*summary['mean_upward_rig_bw']:.3f}%BW; downwardforce0. Lateral RMS is{100*summary['lateral_rig_rms_bw']:.3f}%BW. Mean absolute apparatus torques are{np.round(summary['mean_abs_rig_torque_Nm'],3).tolist()}Nm. Measured groundload is{100*stats['mean_ground_load_bw']:.3f}%BW; groundcontact above1N occurs in{100*contact_fraction:.1f}% of saved100Hz samples. Visualfoot mesh y ranges from{1000*stats['foot_mesh_y_range_m'][0]:.3f} to{1000*stats['foot_mesh_y_range_m'][1]:.3f}mm at2ms sampling. Positive minimum-height intervals mean continuous groundcontact is not claimed. Maximum knee angles are{np.round(stats['knee_max_deg'],3).tolist()}degrees; neither knee hyperextends.

Post-startup mean actual pair tension (t≥0.5s) is{np.round(mean_force,3).tolist()}N. Mechanical work −∫T·dL is{np.round(work,5).tolist()}J per route, not pneumatic or metabolic energy. The two independent repeats also match every saved neural/backend, pressure-history, and mechanical-state array. FINAL_COMPLETE_STATE.npz is the observer’s capture filename; a complete restart API including original runner-local scheduling variables has **not** been implemented or tested.

## Reproduction and files

{command}

The command reuses completed repeats and checks hashes/parity. The original coactivation fixture is cad_mixed_coactivation05.py; direct --mode commands rerun that mode and should not be used to overwrite retained evidence. Small reference-gate pointer fixtures in fresh_a/b explicitly identify the reused original baseline; they are not additional reference simulations. Rawoutputs occupy one external study folder; original model files and accepted Gait checkpoints are untouched.

RESULTS.json summarizes the measured scope. COACTIVATION_REGISTRATION.json, each LIVE_PREREGISTRATION.json, ASSUMED_PAIR_CONFIGURATION.json, MOVING_AND_DEFAULT_PARITY.json, and FRESH_AND_MESH_AUDIT.json provide provenance. CAD_MIXED_PRESSURE_PAIR.png/.pdf/.svg is a short bilateral pressure/tension/knee figure. Source snapshots include original runner, opt-in transformation, and experiment drivers. No manuscripts were edited.
'''
 record+='\nACTUAL_CORE_CONNECTOME.json exports 262 neurons, 650 synapses, and 230 input ports, with matrices verified against the saved fresh run. External allocator, MN current clamp, pressure conversion, and apparatus blocks are explicit; a core-only graph would omit essential controller behavior.\n'
 for old,new in {'fixed5%':'fixed 5%','2s':'2 s','20mm':'20 mm','52reference':'52 reference','54routes':'54 routes','262-cell,230-input':'262-cell, 230-input','single20':'single 20','200kN/m':'200 kN/m','40ms':'40 ms','60ms':'60 ms','4ms':'4 ms','620kPa':'620 kPa','0.1rad':'0.1 rad','0.05in':'0.05 in','sourcepose':'source pose','upward harnessforce':'upward harness force','groundload':'ground load','groundcontact':'ground contact','zero initialvelocity':'zero initial velocity','Actual visualfoot':'Actual visual foot','full-duration default-trajectory':'full-duration default trajectory','source401-pose':'source 401-pose','30ms':'30 ms','2ms':'2 ms','100Hz':'100 Hz','above1N':'above 1 N','not autonome':'not autonome','52routes':'52 routes','all54routes':'all 54 routes'}.items():record=record.replace(old,new)
 for old,new in {'peakpressure':'peak pressure','torqueallocator':'torque allocator','seriesstiffness':'series stiffness','maximumerror':'maximum error','downwardforce':'downward force','Visualfoot':'Visual foot','Rawoutputs':'Raw outputs','saved100':'saved 100','original262':'original 262','existing5%':'existing 5%','all54':'all 54','no-floor mixed2':'no-floor mixed 2','whole-SNS-necessity':'whole SNS necessity','%BW':'% BW','NFfull':'NFfull'}.items():record=record.replace(old,new)
 lines=[]
 for line in record.splitlines():
  if 'D:\\' not in line and not line.startswith('& '):
   line=re.sub(r'(?<=[A-Za-z])(?=\d)',' ',line)
   line=re.sub(r'(?<=\d)(?!(?:e|E)[+-]?\d)(?=[A-Za-z])',' ',line)
  lines.append(line)
 record=('\n'.join(lines)+'\n').replace('easteregg 2','easteregg2')
 (BASE/'EXPERIMENT_RECORD.md').write_text(record,encoding='utf-8');snap=BASE/'source_snapshot';sources=[Path(__file__),Path(J.__file__),Path(J.C.__file__),Path(J.L.__file__),HERE/'cad_mixed_pressure_followup.py',HERE/'export_cad_mixed_connectome.py',J.C.ROBOT/'run_neural_robot.py'];manifest={}
 for p in sources:
  target=snap/(p.stem+'_'+hashlib.sha256(p.read_bytes()).hexdigest()[:8]+p.suffix)
  if not target.exists():shutil.copyfile(p,target)
  manifest[str(p)]=dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),snapshot=str(target))
 (snap/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(stats,indent=2))
if __name__=='__main__':main()
