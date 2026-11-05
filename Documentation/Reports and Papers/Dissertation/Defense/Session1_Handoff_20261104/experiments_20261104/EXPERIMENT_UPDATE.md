# Session 1 calculation update: Gait2392 and postural control

Machine: **easteregg2**. Repository: `D:\GitHub\Bipedal_Robot`.
These experiments followed Ben's authorization to resume calculation on
4 November 2026. This update supersedes the older preparation-only status in
the handoff snapshots. Manuscript and Overleaf files were not edited.

## What works in the body simulation

**Controller identity correction:** `hunt` is an internal software selector for
Ben Bolen's **BiW2Lvar+ controller**, calibrated using Hunt's four-step
method and evaluated through Bolen's five-stage curriculum. It is not a published BiW2Lvar+ model. It extends W2Lvar through
the functional/axis PF and joint-integration implementations, including the
four-cell RF/hamstring PF module. The simpler `w2lvar` selector below is a
different, less expanded candidate. Historical filenames and selector names
remain unchanged so results can be reproduced.

The corrected four-command **graded biarticular multi-PF W2Lvar** candidate completed both suspended
air stages with the same learned weights and shared PF→MN expression scale of
0.1. All four descending command channels are explicit; stages differ in their
specified sensory/command conditions, not independently fitted left/right gains.

| Physical test | Principal knee cycles R/L | Period | Bilateral phase error | Maximum knee extension R/L |
|---|---:|---:|---:|---:|
| Stage 1, deafferented air | 4 / 4 | 2.466 s | 0.0337 cycles | −1.300° / −1.326° |
| Stage 2, afferented air | 7 / 7 | 1.300 s | 0.0933 cycles | −1.718° / −1.801° |

The stage-2 spring harness supplies **99.860% BW**, averaged from 5 s onward in
the 16 s protocol. Prescribed lift is zero, and ground contact is disabled.
This is suspended neural stepping. Frontal guides remain active. The measured
foot forces are zero, as expected in air.
The suspension-force measurement is bookkeeping for the air test, not evidence
that this controller requires that percentage of assistance for ground walking.
No minimum ground-support requirement was measured in this stage.

The first full-expression body probe triggered the predictive knee guard at
4.914 s, while the knee was still flexed at −0.767°. Higher air-expression
settings also stopped at the guard. The 0.1 candidate completed the full runs;
the guard was never disabled, and joint state was never clipped.

Stage 1 initially failed the legacy alternation counter because it counted a
secondary knee-flexion hump as another cycle. The separate principal-cycle
audit uses joint autocorrelation and separated flexion peaks. It passed known
antiphase double-hump data and rejected synchronous, frozen, and 10 Hz examples.
The original assessment is preserved alongside the revised audit. Stage 2 passes
both the legacy and principal-cycle checks. The secondary humps remain visible
in the traces and are a waveform-quality issue for further work.

Presentation assets copied beside this document:

- [Latest Gait2392 GIF](Gait2392_BiW2LvarX_pose_stage2_air.gif).
- [Still](Gait2392_BiW2LvarX_pose_stage2_air.png).
- [Aligned RG/PF/MN traces](Gait2392_BiW2LvarX_pose_stage2_neural.png).
- [Aligned joints, foot loads, and harness forces](Gait2392_BiW2LvarX_pose_stage2_mechanics.png).
- [Combined vector plots](Gait2392_BiW2LvarX_pose_stage2_verification.pdf).

## Ordered neural training across the four families

| Gait2392 family | Intrinsic rhythm A | Expected-afferent timing B | Feedforward peak C error | Reflex waveform D error | Outcome |
|---|---|---|---:|---:|---|
| BiW2Lvar+ | Pass | Pass | 0.2229 | 0.2148 | Offline A–D passed; physical air result above |
| s3k | Pass | Pass | 0.1221 | 0.3677 | D failed; no promotion to body curriculum |
| w2lvar | Pass | Pass | 0.1216 | 0.4027 | D failed; no promotion to body curriculum |
| syn6 | Pass | Pass | 0.1205 | 0.2471 | Offline A–D passed; slower physical air candidate needs cadence and ankle-waveform repair |

Errors are dimensionless normalized MN/activation errors, not joint-angle errors.
The declared peak/waveform gate is 0.25, with additional finite, rhythmic,
period, bilateral-phase, and dynamic-PF checks. Fresh-process repeats are saved.
The targets are engineering inverse-mechanics activation estimates from the
existing IK target file, rather than measured EMG. Their asymmetric and sometimes
saturated muscle peaks are enumerated in each TARGET_AUDIT.json.

These names identify different Gait2392 circuit architectures. BiW2Lvar+ includes
joint/functional PF channels and the four-cell RF/hamstring module. s3k uses
phase-window PF channels. w2lvar uses hip PF plus coupled knee/ankle PF pairs.
syn6 uses six fitted primitive/synergy channels. Innervating a biarticular muscle
does not by itself mean the controller has a separate biarticular PF module.
These runs are distinct from the 12-actuator AnimatLab/W2L virtual walker clips.

Initial attempts exposed three issues. First, the original intrinsic operating
points did not support the required free rhythm. The accepted BiW2Lvar+ RG design
was tested as an explicit seed for the other families, with PF intrinsic
parameters retained. Second, inherited commissural excitation could synchronize
the legs. Shared flexor/extensor commissural budgets were fitted while retaining
mirrored wiring. Third, one scale per broad muscle group was insufficient for
peak fitting. A bounded per-muscle scale, shared by the two sides, passed the
peak checks. s3k and w2lvar still require waveform/phase design work; their
remaining errors are not fixed by declaring peak success to be walking.

The final C/D fits separate feedforward peak fitting from local reflex refinement.
For Hunt, the expected replay also applies the pose-dependent recruitment gates
used by its live controller, at the same five-neural-step update interval. The
earlier ungated replay and body probes remain as diagnostics. All four selected
profiles passed strict graph fault injection for single inhibitory contact
entries and toe-to-RG entries. Heel excites ipsilateral InE and contralateral
RG-F, with extensor-side PF relays. There is no heel-to-both-IN fanout or baseline
toe-to-RG projection. Ground-triggered operation is a separate test.

Syn6's first body attempts exposed two historical-evaluator interface assumptions
(Hunt initialization-port names and a stock PF logger). Those were repaired in
the new body adapter, without adding neural populations or relabeling synergy
cells. Its corrected stage-1 run passes both physical counters: 3/4 cycles,
2.460 s period, phase error 0.00373, and 99.876% BW spring support. Stage 2
completes safely but contains fast joint/ankle activity; the legacy distance-based
peak counter aliases this into an apparent pass. The spectral guard rejects it.
Saved spectra show the RG near 0.44–0.56 Hz while ankle MNs/joints reach about
5.2 Hz in one reduced-central-feedback probe. This localizes a problem downstream
of RG timing. Central-feedback and local-II ablations are stored separately,
with all failed candidates retained.

## Postural premotor experiment

The tested 14-cell module routes descending task commands through antagonistic
premotor cells and inhibitory relays before MNs. Ia/II enter the premotor layer;
Ib polarity is explicitly selectable. It retains a small local Ia→MN pathway.
The task-error calculation uses a stated external PD law. The diagram and exact
edge ledger are [SVG](POSTURAL_SUBNETWORK.svg), [PNG](POSTURAL_SUBNETWORK.png),
and [JSON](POSTURAL_SUBNETWORK.json). [Response traces](POSTURAL_RESPONSE.png)
show premotor/MN voltages, joint angle, and tensile forces.

This is a fixed-base 10 kg ankle bench with synthetic tensile actuators. The
baseline and quiet-condition fitted circuits each passed 12/12 initial held-out
push/noise cases. The fit reduced quiet residual RMS angle from 0.302° to 0.036°.
However, with additional measurement delay, baseline/fitted passes were 9/12
versus 0/12 at 30 ms, and both were 0/12 at 70 ms. A new bounded robust fit includes
delay and noise in its objective and uses unseen seeds for final checks. It is
recorded separately. The quiet-condition fit is not deployed into the full
walkers on the strength of its nominal result.

The compact and normal SNS implementations gave identical bench trajectories.
For this tiny circuit the normal implementation was faster: approximately
1.99–2.40 s versus 2.66–2.73 s for 8 simulated seconds in the paired checks.
The larger Gait2392 models use the verified compact implementation. The latest
16 s BiW2Lvar+ body run took 32.60 s wall time, including the evaluator's overhead;
historical kernel-only acceleration ratios must not be applied to this full run.

## Full model and experiment paths on easteregg2

Exact compiled stage-2 body model:

`D:\Users\Ben Bolen\ControllerTraining_20261104\hunt_pose_air_validation_v3\stage2\capture\model.mjb`

Configuration and actual same-run state recording:

`D:\Users\Ben Bolen\ControllerTraining_20261104\hunt_pose_air_validation_v3\stage2\body.json`

`D:\Users\Ben Bolen\ControllerTraining_20261104\hunt_pose_air_validation_v3\stage2\capture\aligned_trace.npz`

Both-stage assessment and source bindings:

`D:\Users\Ben Bolen\ControllerTraining_20261104\hunt_pose_air_validation_v3\COMPLETE.json`

`D:\Users\Ben Bolen\ControllerTraining_20261104\hunt_pose_air_validation_v3\stage2\SOURCE.json`

Base converted model, loaded through the existing repair/support evaluator:

`D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml`

The exact compiled model includes the existing shared force-law patch named in
body.json. The native MuJoCo muscle law is used; no compliant-tendon, damping,
or BPA-physics variant was enabled by this campaign.
Its body mass is 78.4646 kg (769.7377 N weight); this is the saved human model's
mass, distinct from the CAD robot's requested mass target.

Final ordered neural output fits:

`D:\Users\Ben Bolen\ControllerTraining_20261104\ordered_output_fits_v5\COMPLETE.json`

Per family, `hunt`, `s3k`, `w2lvar`, and `syn6` contain C/D manifests, weights/scales,
costs, trials, source hashes, saved traces, and fresh-process repeat results.
Timing dependencies are identified in each C manifest. Strict contact tests:

`D:\Users\Ben Bolen\ControllerTraining_20261104\strict_contact_contract_v1\RESULTS.json`

Postural records and current robust fit:

`D:\Users\Ben Bolen\PosturalPilot_20261104\held_out_v1\SUMMARY.json`

`D:\Users\Ben Bolen\PosturalPilot_20261104\delay_noise_v2\SUMMARY.json`

`D:\Users\Ben Bolen\PosturalPilot_20261104\robust_delay_fit_v3\`

That first delay-aware fit passed 6/6 training cases but only 9/18 unseen cases
at the larger noise amplitude. It is not accepted. A follow-up includes that
noise level and shared premotor time constants in the fit:

`D:\Users\Ben Bolen\PosturalPilot_20261104\robust_delay_fit_v4\`

The follow-up passed 6/6 training cases, reproduced its training objective in a
fresh process, and passed 17/18 unseen cases. The remaining 70 ms-delay case had
0.694° final RMS angle and settled after 3.125 s, missing the predeclared 3 s gate.
The gate was retained. Thus it is an improved candidate, with a quantified
remaining failure, rather than a fully accepted postural controller.

Syn6 body records and targeted feedback tests:

`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_air_validation_v3\COMPLETE.json`

`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_body_feedback_v4\`

`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_localII_feedback_v5\`

Reducing local Type-II gain removed the fast dominant component at low gains,
but knee excursion remained about 13–16°, below the 20° movement gate. A bounded
12-trial Rustuna fit jointly calibrated that gain and PF expression in the actual
body, with a 1.233 s cadence target and fresh-process checks:

`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_body_fit_v6\`

The selected slow candidate uses expression 0.2913677891 and local-II scale
0.2114065546. It reproduced both principal air-motion gates in new processes,
with stage-1/2 periods 2.466/2.004 s and stage-2 support 99.867% BW:

`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_slow_candidate_repeat_v7\COMPLETE.json`

It fails the 1.233 s target-cadence criterion. Its saved-state spectral audit also
finds substantial secondary ankle power at 6.4–7 Hz, despite a slow dominant
component. Therefore, a principal knee-cycle pass does not make this a selected
presentation-quality gait. The audit is STAGE2_SPECTRA.json in the same folder.
A separate neural command map found that increasing equal F/E command from
0.15 to 1 shortened the period from 2.276 to 1.512 s but synchronized the legs
above 0.35; simply increasing drive did not solve the problem:

`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_drive_response_v2\RESULTS.json`

## Layer-specific current perturbations

The frozen BiW2Lvar+ candidate completed separate 100 ms, 2 nA pulses at RG-F,
hip-flexor PF, and iliacus MN, and a −2 nA PF pulse, with protected knees.
No voltage clamps were used. Both the 7.1 s silent-PF condition and 7.65 s
active-PF condition were tested. All comparisons have identical pre-pulse
neural trajectories and actually recorded current samples.

At the silent phase the PF pulse changes its voltage but remains below the
synaptic transmission threshold, producing essentially no MN response. At the
active phase it changes PF by up to 1.14 mV and iliacus MN by 0.16 mV. The
active-phase MN pulse changes iliacus MN by 1.79 mV and subsequently RG-F by
6.04 mV through the closed neuromechanical loop. These are paired maximum
differences over the displayed response window, rather than direct synaptic
gains. Sensory return means a PF/MN perturbation can eventually affect RG.

Raw trials and paired measurements:

`D:\Users\Ben Bolen\ControllerTraining_20261104\hunt_layer_stimulus_2nA_v2\COMPLETE.json`

`D:\Users\Ben Bolen\ControllerTraining_20261104\hunt_active_phase_stimulus_v3\COMPLETE.json`

`D:\Users\Ben Bolen\ControllerTraining_20261104\figures\Hunt_active_phase_stimulus\PAIRED_STIMULUS_METRICS.json`

Portable plots: [active-phase stimuli](BiW2LvarX_active_phase_stimulus.png),
[vector PDF](BiW2LvarX_active_phase_stimulus.pdf), and
[silent-phase stimuli](BiW2LvarX_silent_phase_stimulus.png). Each shows approximately
three cycles of RG, PF, MN, and hip/knee/ankle motion.

## Fixed-weight neural-noise verification

The BiW2Lvar+ stage-2 candidate completed all six ongoing-noise trials: three seeds
each at stationary RG-current standard deviations of 0.01 and 0.025 nA, with
20 ms correlation time. Learned parameters were fixed. Noise was independently
applied to the four RG cells and its actual samples were recorded. Every trial
retained the principal-motion, full-duration, frontal-separation, and knee gates.
The two additional cases were a zero-noise control and a repeated noisy seed;
both passed. The fresh-process repeated seed matched joint states, neural
voltages, and applied noise exactly (maximum absolute differences 0.0).

`D:\Users\Ben Bolen\ControllerTraining_20261104\hunt_rg_noise_v1\COMPLETE.json`

The compact [noise-result copy](BiW2LvarX_RG_NOISE_VERIFICATION.json) travels with
this handoff. ZERO_OPTION_PARITY.json beside the raw summary confirms that
adding the optional zero-noise feature leaves the saved sham trajectory exactly
unchanged in joint state, neural voltage, control, and activation.

This tests suspended-air coordination under continuous neural-current noise;
it does not test ground balance or change the connection parameters by ±1%.
No noisy trial was used to refit the candidate.

## Interpretation and next calculations

The current selected full-body result is BiW2Lvar+'s protected air motion. The four
families have completed ordered offline attempts, including recorded failures.
Standing and ground stepping with this revised controller require the separate
postural integration and physical load/clearance gates. The isolated postural
bench does not establish those full-body outcomes. CAD remains paused for
Ben's routing revision. No compliant-tendon or damping variant was started.

For subsequent work, s3k/w2lvar need PF waveform/phase refinement rather than
another uniform-gain sweep. Syn6 needs the ankle feedback loop resolved and
bilateral timing maintained as cadence increases. The postural candidate needs
its remaining delayed/noisy case resolved before selecting a full-model variant.
All jobs in this bounded batch have completed; no unattended continuation is
implied by the saved launcher commands.

Verification at handoff: 18 copied assets match their recorded SHA256 values;
all 38 Python source files in the two new modules parse; the original Hunt
controller retains SHA256
`00112d00d7e65c700c1c89e651f57354588d88335386209e5803b67f7154b57c`.
The six noise cases, two control/repeat cases, and both four-case stimulus batches
completed. MATLAB processes remain untouched. A compact source snapshot and
its manifest are at
`D:\Users\Ben Bolen\ControllerTraining_20261104\source_snapshot_final`.

Source, continuation commands, and separate training record are in:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\controller_profiles_20261104\TRAINING_RECORD.md`

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\postural_pilot_20261104\README.md`

The GIF, still, plots, and diagram in this folder travel with the repo to
EB475WS4. Compiled models and raw results remain on easteregg2. No agent commit or
push was made. MATLAB jobs and their modified files were left untouched. CAD
simulations remain paused for Ben's routing update. Full-model postural integration
and subsequent standing/ground stages remain separate from the completed air tests.

Canonical display name: **BiW2Lvar+**; plain identifier: **BiW2LvarX**.
Saved `hunt` strings and path names are deprecated compatibility identifiers,
not a model attribution. Hunt's four-step calibration and Bolen's five-stage
curriculum are distinct parts of the workflow; the revised candidate has so far
passed physical curriculum stages 1 and 2, not all five stages.



Latest W2L ground calculation: [32-trial result and exact paths](LiDeng_W2L108_RESULT.md), [GIF](LiDeng_W2L108_ground.gif), [still](LiDeng_W2L108_ground.png). Selected candidate repeated byte-exactly; 8/7 clear unloaded swings; 29.519% BW upward rail/guide assistance after 1 s.


Pelvic tilt/overlap follow-up: [full comparison record](PELVIS_OVERLAP_RECORD_20261104.md), [anterior/overlap tradeoff GIF](LiDeng_W2L108_pelvis_overlap_trial.gif), [lower-support posture GIF](LiDeng_W2L108_lower_support_posture_trial.gif). These are labeled diagnostic tradeoffs; no combined gait-quality finalist was accepted.


Short/slow shuffle follow-up: [record and sensory-gate design](SHORT_SHUFFLE_RECORD_20261104.md), [slow candidate with measured direction arrows](LiDeng_W2L108_slow_direction.gif). Fourteen alternating clear swings in 16 s, 10.1% BW upward assistance over 2–16 s; short-stride and sustained double-support targets remain unmet. Fresh numeric trajectories match exactly.

