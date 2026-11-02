"""Save frozen reproducibility manifest and bounded-result handoff, no runs."""
import json,hashlib,sys,importlib.metadata
from pathlib import Path
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\source_domain_recruitment');HERE=Path(__file__).resolve().parent

def item(path):
 path=Path(path);return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())

def main():
 cfg=json.loads((ROOT/'PREREGISTRATION.json').read_text());selected=json.loads((ROOT/'pressure/SELECTED.json').read_text());validation=json.loads((ROOT/'pressure/VALIDATION_STATUS.json').read_text());proof=json.loads((ROOT/'STRONG_FRESH_AUDIT.json').read_text());core=ROOT.parent/'group_diagnosis/NORMALIZED_STANDING_CHECKPOINT.json';result=json.loads((ROOT/'pressure/fresh_training_strong/w2lvar_case.result.json').read_text());m=Path(r'D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml')
 manifest=dict(machine='easteregg2',software=dict(python=sys.version,executable=sys.executable,mujoco=importlib.metadata.version('mujoco'),sns_toolbox=importlib.metadata.version('sns-toolbox'),numpy=importlib.metadata.version('numpy'),numba=importlib.metadata.version('numba'),optimizer='Rustuna0.2.1 externalPython3.11 persistentmyo worker'),program=item(HERE/'mn_recruitment_experiment.py'),model=item(m),core_checkpoint=item(core),core_saved_configuration=json.loads(core.read_text()),force_law=item(result['force_law']),corrected_pose_checkpoint=item(result['checkpoint']),preparation=cfg['common_preparation'],experiment_registration=item(ROOT/'PREREGISTRATION.json'),source_manifest=item(ROOT/'source_snapshot/MANIFEST.json'),selected_parameters=selected['params'],additional_neuron_parameters=dict(tau_s=.006,capacitance_nF=6.,leak_uS=1.,input_cap_nA=10.,headroom_N=5.,error_span_N=5.,preview_s=.18,extra_transport_s=.006,excitatory_reversal_mV=40.,max_conductance_uS=selected['params']['mn_current_cap_nA']/40.,MN_full_command_voltage_mV=5.,native_physics_neural_step_s=.002),measurement_assumptions=cfg['measurement'],core_neurons=888,additional_valve_neurons=2,additional_MN_recruitment_neurons=2,physical_actuator_assumptions='Unchangedvirtual two20mmBPAs/route,series pairK400kN/m,25%stroke,50mmfixedlength,restderivedfromnominalrouteenvelope. Staticcharacterization reused; dynamic/geometry assumptions not mountedCAD identification.',physical_pressure_ceiling_kpa=620.,fixtures=cfg['training']+cfg['validation'],frozen_validation_passes=sum(c['accepted'] for c in validation['validation']),frozen_validation_total=len(validation['validation']),fresh_nominal_30s_exact=validation['fresh_exact'],fresh_known_strong_20s_exact=proof['accepted'] and proof['all_training_array_hashes_exact'],no_postvalidation_reselection=True,scope='Bounded gradedvirtualBPA guidedstanding only; failed0.78map remains. No fullcontrollerbasin,gait,spikingtransfer,CADhardware,independentbalance,or longstrongdurability claim.')
 frozen_path=ROOT/'FROZEN_CONFIGURATION.json'
 if frozen_path.exists():
  previous=json.loads(frozen_path.read_text())
  for key in ('program','model','core_checkpoint','force_law','corrected_pose_checkpoint','experiment_registration','source_manifest'):
   if previous[key]['sha256']!=manifest[key]['sha256']:raise ValueError('Frozen scientific dependency changed; do not rewrite checkpoint: '+key)
  if previous['selected_parameters']!=manifest['selected_parameters']:raise ValueError('Frozen selection changed')
 recorded={};spinal=HERE.parents[1]/'spinal';snap=ROOT/'source_snapshot'
 import shutil
 for name,expected in result['source_sha256'].items():
  p=HERE/name if (HERE/name).exists() else spinal/name
  if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise ValueError('Recorded numerical source changed: '+str(p))
  target=snap/(p.stem+'_'+expected[:8]+p.suffix)
  if not target.exists():shutil.copyfile(p,target)
  recorded[str(p)]=dict(sha256=expected,snapshot=str(target),scope='Matches actual fresh-result recorded hash')
 for name,expected in result['pressure_model_source_sha256'].items():
  p=HERE/name
  if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise ValueError('Recorded pressure source changed:'+str(p))
  target=snap/(p.stem+'_'+expected[:8]+p.suffix)
  if not target.exists():shutil.copyfile(p,target)
  recorded[str(p)]=dict(sha256=expected,snapshot=str(target),scope='Matches actual fresh-result recorded hash')
 (snap/'RECORDED_NUMERICAL_DEPENDENCIES.json').write_text(json.dumps(recorded,indent=2),encoding='utf-8');manifest['recorded_numerical_dependencies']=recorded
 frozen_path.write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 python=r'D:\Anaconda\envs\myo\python.exe';program=HERE/'mn_recruitment_experiment.py';text=f'''# Reproduce the frozen MN recruitment example

Machine: **easteregg2**. Program: Python, MuJoCo {manifest['software']['mujoco']}, SNS-Toolbox {manifest['software']['sns_toolbox']}, verified sparse graded/gated SNS backend, and Numba BPA force/pressure kernels. The four added relay equations are independently checked against SNS-Toolbox.

The precise configuration and dependency hashes are in `{ROOT / 'FROZEN_CONFIGURATION.json'}`. The program rebuilds the same 888-neuron core plus two valve and two MN relays. It uses a shared 5 nA MN current cap, a 10 kPa valve increment cap, capacity weight 0.75, 180 ms preview, 5 N headroom, 6 ms leak, and 6 ms additional sensor transport. All other neural and model settings come from the saved normalized-w2lvar core checkpoint and validated source files.

First check the frozen scientific files. This stops if their recorded hashes or the selected parameters differ:

```powershell
& '{python}' '{HERE / 'write_mn_recruitment_handoff.py'}'
```

To run a new nominal 30 s study from the original common preparation:

```powershell
& '{python}' '{program}' --case fresh_nominal_user
```

Outputs: `{ROOT / 'pressure/fresh_nominal_user'}`. An existing completed case is reused; choose another `fresh_nominal_...` suffix for a genuinely new repeat. These runs are headless and save body/neural/pressure traces. Full raw output is appropriate for deliberate repeats, not large searches.

Completed selected strong 20 s record: `{ROOT / 'pressure/fresh_training_strong/w2lvar_case.result.json'}`. Full body trajectory: `{ROOT / 'pressure/fresh_training_strong/w2lvar_case.npz'}` (`qfull`; original runner time includes 0.6 s pinned preparation). Pressure and motor telemetry use released-physics time. The neural/force figure is `{ROOT / 'MN_RECRUITMENT_STRONG.pdf'}`; its generator is `{HERE / 'plot_mn_recruitment_finalist.py'}`.

The frozen candidate passes all three known 20 s training conditions, two fresh 30 s nominal repeats exactly, and a separate fresh strong 20 s repeat exactly. Independent validation passes 3/4; the 0.78 force-map case fails at 18.072 s. The matched MN-off/valve-only control fails the strong case at 17.294 s. These results support a reproducible improvement for that finite challenge. They are guided standing tests with virtual BPA assumptions; walking, independent balance, calibrated CAD hardware, and general robustness remain unverified. Strong-case control becomes pulsatile under pressure delay. The plotted cycles are corrective feedback rather than gait CPG rhythms or action potentials. A complete-controller ±1% basin check remains outstanding.

Separate experiment record: `{ROOT / 'EXPERIMENT_RECORD.md'}`. Keep prior checkpoints intact. Rustuna BEST identifies the selected experiment, not an unconditional production winner. The global accepted checkpoint and manuscripts were unchanged.

Independent nominal-repeat archives remain physically independent. Automatic approval review rejected hard-link deduplication because it would merge independent evidence into shared files. Both originals were retained unchanged.
'''
 (ROOT/'HOW_TO_REPEAT.md').write_text(text,encoding='utf-8')
 final_report=f'''# Source-domain recruitment: frozen guided-standing result

Machine: **easteregg2**. The added MN pathway reproducibly completes the known strong 20 s condition (19.4 s released), while its matched valve-only control stops at 17.294 s. Independent validation passes **3/4**. The weaker 0.78 force-map case fails at 18.072 s, so general robustness is not established. Parameters remain frozen.

## Model and architecture

This is the corrected-pose Gait2392 human model with declared virtual BPA geometry and pressure dynamics. It retains strong lateral and rotational apparatus guidance. The 888-neuron normalized-w2lvar core, two existing valve relays, and two new MN relays represent 892 modeled cells, with additional gating and filter states. This is a standing-specific engineering variant, not a literature-controller reproduction or curriculum passer.

The new paired 6 ms relays compare measured tension against a predicted geometric minimum plus 5 N. The predictor uses measured route length, positive route velocity, assumed rest/fixed lengths and series stiffness, and a 180 ms preview. It never receives the case's true force multiplier. Positive conductance currents, capped at a shared 5 nA, excite the corresponding pectineus MNs. Original tonic inputs and synapses remain intact. The frozen valve observer retains its 10 kPa increment cap and capacity weight 0.75. Actual force, pressure, and membrane voltage are never assigned to a minimum or projected.

The static bilateral gravity fit sets the original neural standing bias once. It is not an online torque allocator. The common nearest-feasible nominal reference conversion remains explicit in both active and control conditions; consequently a low original neural force request can already produce a higher actual tension. The new signal uses actual tension to avoid confusing that latent command gap with missing physical force.

Sensor assumptions are an additional 6 ms transport delay, pressure noise SD 0.02 kPa, and tension noise SD 0.05 N, clipped at three SDs. These are declared model assumptions, not hardware identification. The graded relay equations and bounded conductance-current conversion match their independent checks exactly. Complete filters, queues, noise consumption, neural states, and actuator/mechanical preparation are audited.

## Source and load audit

All 86 configured routes were checked at 12 observed endpoints. Only pectineus targets were below the geometric minimum in those samples. Observed failing route lengths reconstruct within 1e−10 m; equality-aware pulling moment arms agree with local finite differences within 1e−6 m. Explicit offline 5 N headroom probes show overall restoring work and dissipative strong-case power, with an adverse lateral component. Loaded-contact COM margins were about 104–109 mm. This supported only a bounded guided-standing trial, with new lateral-force guards; it was not a proof of general balance.

Native coordinate signs require care: pectineus has the same hip-flexion torque sign as psoas and iliacus in this model. A positive native joint named hip_flexion_r moves the knee toward negative world x, opposite the measured toe-forward direction at this posture. No global sign or route repair was made, and original OpenSim coordinate parity was not newly tested. See REFERENCE_SIGN_CHECK.json.

## Training, fresh checks, and validation

Four Rustuna trials took 500.3 s. Two candidates passed all three known training conditions; the registered objective selected the shared 5 nA candidate before validation. There was no post-validation reselection.

| Check | Result |
| --- | --- |
| Nominal, moderate, and strong training | All pass with one frozen configuration |
| Two fresh nominal 30 s runs | Both pass; all arrays and states match exactly |
| Separate fresh strong 20 s run | Pass; every training array hash matches exactly |
| Independent interior uncertainty | Pass |
| Independent slower fill/delay | Pass |
| Independent 0.78 force multiplier | Fail at 18.072 s; unchanged pectineus domain guard |
| Independent forward pulse | Pass |
| Matched MN-off nominal and moderate | Pass |
| Matched MN-off strong | Fail at 17.294 s |

The active pathway is necessary for passing this particular finite strong challenge under the compared settings. This does not establish necessity of the whole SNS, universal actuator compensation, or long-term strong-case stability. No complete-controller ±1% basin check, gait result, spiking transfer, or calibrated CAD standing claim is made.

## Measured strong-case limits

Mean upward apparatus force is 1.137% BW; downward force is 0.0039% BW. Lateral RMS/peak forces are 0.278%/5.86% BW. Ground load is approximately 98.57% BW. Peak pressure is 338.66 kPa, below the unchanged 620 kPa limit. Added MN current peaks at 4.691 nA. Knee maxima remain approximately −1.93°/−1.98°, and the visual-foot clearance minimum is −4.057 mm. Strong lateral and rotational guidance remains present.

The new feedback becomes pulsatile under actuator delay. MN_RECRUITMENT_STRONG.png/.pdf/.svg shows a short right-pectineus window; these cycles are corrective feedback, not a gait rhythm or action potentials. The caption states the timing and model limitations. Smooth continuous force tracking is not demonstrated.

## Reproduction and records

Use HOW_TO_REPEAT.md and FROZEN_CONFIGURATION.json for full easteregg2 paths, software versions, controller settings, and dependency hashes. PREREGISTRATION.json records the pre-fit conditions, gates, and selection rule. BASELINE_GATE.json, INDEPENDENT_RECRUITMENT_CHECKS.json, STRONG_FRESH_AUDIT.json, and pressure/VALIDATION_STATUS.json hold verification evidence. TRAINING_COMPARISON.csv records per-case currents, pressure, apparatus forces, and gates. CIRCUIT.json describes the exact added pathway and log channels.

Search outputs are compact and windowed; full selected/fresh evidence and limited unique causal failures are retained outside the repository. The study occupies about 790 MB. Automatic approval review rejected hard-link deduplication; both independent repeat archives remain unchanged. Global accepted checkpoints and manuscripts are untouched.
'''
 (ROOT/'EXPERIMENT_RECORD.md').write_text(final_report,encoding='utf-8')
 (ROOT/'MN_RECRUITMENT_STRONG_alt.txt').write_text('Four panels show about four right-pectineus corrective feedback pulses during guided virtual-BPA standing. Error input drives a leaky neural relay and normalized motor-neuron command. Actual tension remains roughly70–100N while neural force commands pulse, and actual pressure is smoothed by the delayed actuator response. A vertical marker indicates the separate valve-only control stopping at17.294s. This is a standing feedback example, not gait rhythm or action potentials; the full frozen model passed only3of4 independent validation conditions.',encoding='utf-8')
 protocol={};snap=ROOT/'source_snapshot'
 import shutil
 for name in ('mn_recruitment_checks.py','mn_recruitment_followup.py','plot_mn_recruitment_finalist.py','write_mn_recruitment_handoff.py'):
  source=HERE/name;h=hashlib.sha256(source.read_bytes()).hexdigest();target=snap/(source.stem+'_'+h[:8]+source.suffix)
  if not target.exists():shutil.copyfile(source,target)
  protocol[str(source)]=dict(sha256=h,snapshot=str(target))
 (snap/'PROTOCOL_MANIFEST.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8');print(str(ROOT/'FROZEN_CONFIGURATION.json'))

if __name__=='__main__':main()
