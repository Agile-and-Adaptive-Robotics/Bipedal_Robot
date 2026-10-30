"""Curate the completed four-cell campaign without promoting failed stages."""
import os,json,hashlib,shutil,ast
from pathlib import Path
import numpy as np
H=Path(__file__).resolve().parent;names=['fourcell_hunt_20261029_v1','fourcell_multicondition_20261029_v1','fourcell_dynamic_20261029_v1'];records={}
for name in names:
    root=H/'runs'/name;done=json.loads((root/'TRAINING_COMPLETE.json').read_text());summary=json.loads((root/'curriculum/summary.json').read_text())
    if 'results' in summary:summary=summary['results']
    records[name]=dict(training=done,physical=summary,neural=json.loads((root/'weights/neural_validation.json').read_text()))
    if name!=names[0]:
        # The general convenience replay scored after 4 nominal periods.
        # Score the intrinsic-cadence condition after four of its own periods.
        z=np.load(root/'weights/step_C_deafferented.npz');keep=z['t']>=4*2.466
        score=float(np.sqrt(np.mean((z['activation'][keep]-z['target'][keep])**2)))
        records[name]['intrinsic_cadence_C_rmse_after_four_cycles']=score
        for block in ('pf','bi','mn'):
            a=json.loads((root/(block+'_afferented')/(block+'_weights.json')).read_text()) if (root/(block+'_afferented')/(block+'_weights.json')).exists() else None
            b=json.loads((root/(block+'_deafferented')/(block+'_weights.json')).read_text()) if (root/(block+'_deafferented')/(block+'_weights.json')).exists() else None
            if a is not None and b is not None:assert list(a['edges'])==list(b['edges'])
latest=H/'runs'/names[-1];old=H/'runs'/names[0]
benchmark=json.loads((old/'acceleration_check.json').read_text());dynamic=json.loads((latest/'DYNAMIC_TRAINING.json').read_text());repeat=json.loads((latest/'gate_verification/replay_comparison.json').read_text())
latest_bench_path=latest/'acceleration_check.json'
latest_bench=json.loads(latest_bench_path.read_text()) if latest_bench_path.exists() else None
if latest_bench and 'Concurrent physical probes were running.' in latest_bench['scope']:
    latest_bench['scope_original']=latest_bench['scope']
    latest_bench['scope']=latest_bench['scope'].replace('Concurrent physical probes were running.','This latest check ran after this campaign physical queue completed; numerical measurements are unchanged.')
    latest_bench['metadata_correction']='The concurrency clause was inherited from the earlier benchmark; only that context sentence is corrected.'
    latest_bench_path.write_text(json.dumps(latest_bench,indent=2))
ground=json.loads((latest/'gate_verification/stage4_replay.result.json').read_text());contacts=ground['ground_contacts']
inputs=[H/'runs/step_A/winner.json',H/'runs/step_B/winner.json',H/'runs/step_B_left/winner.json',H/'runs/targets/targets.npz',H/'runs/targets/metadata.json',H/'ground_contact_check.py']
meta=json.loads((H/'runs/targets/metadata.json').read_text());inputs += [Path(meta['source_ik']),Path(meta['source_grf'])]
model=Path(os.environ.get('AARL_MODEL',r'D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml'));inputs.append(model)
manifest={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
(latest/'INPUT_PROVENANCE.json').write_text(json.dumps(manifest,indent=2))
(latest/'CURATED_SUMMARY.json').write_text(json.dumps(dict(campaigns=records,benchmark=benchmark,latest_benchmark=latest_bench,dynamic_fit=dynamic,replay_differences=repeat,verified_ground=contacts),indent=2))
rows=[]
for i,name in enumerate(names):
    r=records[name];a=r['physical']['stage1']['assessment'];s2=r['physical']['stage2_plain' if i==0 else 'stage2']['assessment']
    rows.append(f"| {['Single cadence','Two cadences','Two cadences + native paired fit'][i]} | {r['neural']['activation_rmse']:.5f} | {a['cycle_phase_error']:.3f} | {min(a['knee_ranges_deg']):.2f} | {s2['cycle_phase_error']:.3f} |")
report='''# Four-cell RF/hamstring PF: training and physical validation

Session 1, easteregg2, 2026-10-29. Human Gait2392, 92 muscle actuators.
This is separate from CHATGPT_REPORT.md. No dissertation source was changed.

## Outcome

Three fits of the four-cell version followed Hunt's sequence and were physically tested.
No variant passed the complete corrected curriculum. Assisted standing passes;
ground stepping remains unsuccessful. Stage 5 was withheld because verified
alternating load-bearing steps were not obtained. Later-stage runs after an
earlier failure are diagnostics, not curriculum graduation.

The RF/HAM module has two graded PF integrators and two reciprocal inhibitory
relays per leg. It is not another intrinsic NaP oscillator. Sartorius and
gracilis retain separate function-specific modules and independent MN outputs.
The conditional census is 1,025 neurons in deafferented air, 1,027 with the
commissural candidate, and 1,035/424 inputs with the added PT branch.

## Hunt training

A/B were revalidated before fitting C: intrinsic period 2.466 s; isolated
sensory-driven right/left periods 1.236/1.212 s. The coupled expected-afferent
network gives 1.264/1.246 s. These are timing checks, not physical gait results.

C fits joint PF, RF/HAM integration, and MN output sequentially, without old D
weights contributing to the feedforward fit. D then fits Ia/Ib and optional
IaIN refinement with C frozen. All motor targets are inverse-model estimates,
not EMG. The 80% C target does not establish an additive physiological 80/20 split.

The first fit uses expected sensory-driven timing. The second fits C jointly
at that cadence and at the slower, deafferented intrinsic cadence. The third
also fits the native four-cell dynamics, including inhibitory conductances,
then refits MN output and D. Its small preallocated model reproduces the full
network's recorded module dynamics to 4.44e-16 mV before optimization.

| Fit | C+D activation RMSE | Stage-1 phase error | Stage-1 smaller knee ROM (deg) | Stage-2 phase error |
|---|---:|---:|---:|---:|
'''+ '\n'.join(rows)+'''

The table uses nominal, unscaled physical cases. Stage 1 in the first fit
fails coverage/frontal checks. The later fits improve deafferented phase and
frontal behavior, but the left knee remains below the 20-degree ROM gate.
Their afferented air motion still fails phase consistency. Smaller supervised
error alone is not a sufficient controller-selection objective.

A separate single-cadence calibration probe passes the air motion/frontal
screen: 4/4 knee cycles, phase error 0.1437, period 2.166 s, and no sampled
crossing. It changes commissural/reflex and motor gains after the fit, has no
all-synapse robustness certificate, and is not a curriculum winner.

## Actual ground contact and support

Ground diagnostics include static harnesses, gradual lowering by 40–120 mm,
explicit constant-speed towing, lateral restraint reduction, intact stance
output, and timing/motor-reflex ablations. Towing is mechanical assistance;
it is not self-propelled locomotion. None supplies a periodic stride clock.

The strongest loaded diagnostic from the native paired fit transfers
'''+f"{100*sum(contacts['mean_foot_bw']):.2f}% BW through the feet and {100*contacts['mean_harness_vertical_bw']:.2f}% through the vertical harness."+'''
Both feet remain in contact. This is a load-sharing result, not walking or a
successful support-weaning rung. Removing motor reflexes and timing inputs
also fails to produce sustained alternating ground contacts.

Recorded bilateral RG activity continues bursting on the ground. The initial
paired PF module is biased toward right RF and left HAM, and the feet fail
to unload. Thus continuous stance contact must not be described as proof that
the RG itself stopped. The motor/force transfer and phase coordination remain
unresolved; lower neural replay error has not solved them.

The corrected physical gate now uses once-counted GRFs, 0.05/0.025 BW
hysteresis, at least 80 ms preceding unloading and 100 ms loaded stance,
three contacts per leg, alternating contact phases, appropriate duty, physical
motion, frontal clearance, and pelvis tilt. Foot-body-origin heights are not
sole clearance; full foot-edge clearance/self-collision is not certified.
Legacy duplicated heel/toe contact totals remain diagnostics, not acceptance
evidence. Synthetic checks accept alternating loading and reject synchronous
loading, permanent double stance, and airborne motion. These synthetic cases
are tests of the gate, not successful walker trials.

A fresh replay with the corrected gate is bit-identical in qfull, control,
and recorded neural traces to the earlier loaded diagnostic. The new gate
does not change the model dynamics or rescue the failed result.

## Acceleration and records

The first complete fit takes 158.71 s; the two-cadence fit takes 383.13 s.
The native paired optimization performs 768 module trials in
'''+f"{dynamic['wall_s']:.2f} s, excluding MuJoCo physics and the subsequent MN/D refit."+'''
Sparse SNS updates, preallocated arrays, cached edge indices, and vectorized
afferents are used throughout. The scalar/vector afferent encoders agree exactly
in the sampled check. A 500-update fixed-conductance neural comparison takes
'''+f"{benchmark['fast_s']:.3f} s accelerated versus {benchmark['native_s']:.3f} s native ({benchmark['ratio']:.1f}x), with maximum voltage disagreement {benchmark['max_voltage_error']:.3g} mV."+'''
This excludes physics, changing conductances, setup, and logging, and was measured
during concurrent physical work. It is not a 100x full-simulation claim.
Earlier matched full-model measurements remain about 6.7–9.2x end to end.

Raw weights, traces, logs, source snapshots, hashes, and trial configurations
are under the ignored runs/fourcell_* directories. INPUT_PROVENANCE.json and
CURATED_SUMMARY.json accompany the native-pair campaign. Raw convenience
deafferented replay JSON used a nominal-period warmup; the curated summary
also recomputes that score after four intrinsic periods. No raw result was
silently overwritten to change its outcome. Earlier candidates are preserved.

The current viewer is network_figures/fourcell_dynamic/index.html; it includes
the actual four-cell subnetwork. The grounded diagnostic figure uses the
single-cadence full-feedback trial and is explicitly labeled with that scope.
No harness-free walking is demonstrated.
'''
if latest_bench:report+='\nLatest-conductance verification: '+f"{latest_bench['fast_s']:.4f} s sparse versus {latest_bench['native_s']:.4f} s native over 500 updates ({latest_bench['ratio']:.2f}x); maximum voltage disagreement {latest_bench['max_voltage_error']:.3g} mV. The physical queue had completed. The encoder comparison remains exact. This has the same neural-only scope; timings vary across runs.\n"
(H/'FOURCELL_HUNT_RESULTS.md').write_text(report,encoding='utf-8')
V=H/'network_figures';backup=latest/'gate_verification/previous_root_index.html'
if not backup.exists():shutil.copyfile(V/'index.html',backup)
redirect=latest/'gate_verification/root_redirect.tmp'
redirect.write_text('<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0; url=fourcell_dynamic/index.html"><title>Four-cell Gait2392 network</title><p><a href="fourcell_dynamic/index.html">Open the current dynamically trained four-cell network</a></p>',encoding='utf-8')
os.replace(redirect,V/'index.html')
shutil.copyfile(V/'fourcell_dynamic/NETWORK_CENSUS.json',V/'NETWORK_CENSUS.json')
viewer=V/'fourcell_dynamic/index.html';page=viewer.read_text(encoding='utf-8')
notice='<p id="training-status">Dynamically fitted four-cell candidate. Assisted standing passes; the completed alternating ground-step tests fail. <a href="../../FOURCELL_HUNT_RESULTS.md">Training results and limits</a>.</p><details><summary>Balance computation outside SNS</summary><img src="../balance_dataflow.svg"></details><p><a href="../fourcell_ground_diagnosis.pdf">Recorded grounded diagnostic from the single-cadence fit</a></p>'
if 'id="training-status"' not in page:
    page=page.replace('</h1>','</h1>'+notice,1);temp=latest/'gate_verification/viewer_index.tmp';temp.write_text(page,encoding='utf-8');os.replace(temp,viewer)
checked=[]
for p in H.glob('*.py'):
    if p.name.startswith(('fourcell_','train_fourcell','verify_','tune_fourcell','multicondition','summarize_fourcell','plot_fourcell')) or p.name in ('controller.py','motor_fit.py','physical_eval.py','ground_contact_check.py','draw_network.py'):
        ast.parse(p.read_text(encoding='utf-8'),filename=str(p));checked.append(p.name)
print(json.dumps(dict(report=str(H/'FOURCELL_HUNT_RESULTS.md'),physical_trials=sum(len(list((H/'runs'/name).rglob('*.result.json'))) for name in names),syntax_checked=len(checked),verified_ground_contacts=contacts['cycles']),indent=2))
