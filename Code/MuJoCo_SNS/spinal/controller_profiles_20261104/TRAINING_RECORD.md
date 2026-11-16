# Session 1: four-command controller experiments

## Syn6 PF relay timing prototype, 2026-11-05

`syn6_s5_timing.py` first restores S5's two bursts, then applies the same
design to all six bilateral PF channels. RG rising transitions reset software
timers with Schmitt thresholds of 1.5 and 3.5 mV. Each burst uses its own onset
and shutoff thresholds, rise/fall slopes, and excitation gain. Delays scale by
the most recent period of the triggering RG channel. The relay conductances
enter the existing PF membrane equations as currents. This hybrid test
suppresses the original incoming PF connections, including direct sensory
inputs to PF, during offline D replay. Its relay timing is generated online
from RG voltage. Fitting uses reference cycles 2–6 and checks cycles 7–10 at
the same cadence with the same periodic reference.

S5 duty dropped from 100% to 34.0% left and 33.5% right, versus 33.0% and
32.5% targets. Mean measured burst windows were 12.5–25.0% and 73.0–94.5%
on the left, versus 13.0–25.0% and 73.5–94.5%. On the right they were
23.5–44.5% and 63.0–75.5%, versus 22.5–43.5% and 62.5–74.0%.

The simultaneous replay passed all 48 channel-stride checks across 12 PF
outputs and four strides. Worst onset/shutoff error was 2.0 percentage
points of a stride, worst fundamental phase error 13.27 degrees, and worst
samplewise activation error 0.3085. The gate checks matching burst counts,
onset/shutoff error at most 5 percentage points, fundamental lag at most
30 degrees, duty error at most 15 percentage points, and excursion at least
0.2. The combined run repeated exactly in a fresh process. RG traces and MN
conductance arrays matched baseline exactly. MN waveform RMSE was 0.15556
versus 0.15845 at baseline.

Results and source checkpoints are recorded under
`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_PF_relays_20261105`.
`S5\S5_TIMING.png` overlays the original, target, and relay-generated activity.
`combined\RESULT.json` contains every timing measurement and the relay settings.
`combined\FRESH_REPEAT.json` records the repeat. Native SNS implementation and
cadence/generalization and physical curriculum tests remain open. The
downstream optimization hold still applies to the original syn6 revision.

Naming correction: `hunt` below is the saved selector for Ben Bolen's
BiW2Lvar+ controller calibrated with Hunt's four-step method and Bolen's
five-stage curriculum,
not a published BiW2Lvar+ model. It differs from the simpler `w2lvar` selector.
Stage-2 spring forces describe suspension during afferented air stepping;
they do not measure a minimum ground-walking assistance requirement.

Machine: **easteregg2**. Ben authorized simulations and small training runs on
2026-11-04 after the handoff package was prepared. MATLAB jobs were inspected and
left running; this campaign uses one serial, single-threaded fitting lane. The
isolated postural test uses one additional worker. CAD work remains paused.

## Live records and reproducible source

- Final ordered A/B-dependent output fits: `D:\Users\Ben Bolen\ControllerTraining_20261104\ordered_output_fits_v5\COMPLETE.json`.
- Accepted BiW2Lvar+ physical air tests: `D:\Users\Ben Bolen\ControllerTraining_20261104\hunt_pose_air_validation_v3\COMPLETE.json`.
- Initial campaign status, including failed attempts: `D:\Users\Ben Bolen\ControllerTraining_20261104\recalibrated_v2\CAMPAIGN_STATUS.json`.
- Per family: `<campaign>\hunt`, `s3k`, `w2lvar`, or `syn6` contains the compiled
  static checks, exact profile, wiring changes, source hashes, and build settings.
- Per block: `training\A` through `D` contains MANIFEST, PREFLIGHT, TRIALS,
  BEST, FRESH_REPEAT, FINAL, and BEST_TRACE. Later blocks require every earlier
  FINAL to pass, including a new-process repeat. A failed block stops that family.
- Reusable runner: [campaign.py](campaign.py). It has an overall wall-time limit,
  reuses completed blocks, and preserves interrupted folders for inspection.

These are **graded neural replay experiments**. A calibrates the undriven rhythm
to twice the 1.233 s target period. B fits shared Ia/II/Ib input gains against
replayed expected afferents. The initial C used six broad output groups; the
final C fits one per-muscle scale shared across sides. D freezes those gains
and adjusts local reflex gains
against activation waveforms. Targets come from the existing inverse-mechanics
target file, not measured EMG. None of these offline blocks integrates the body.
Physical curriculum acceptance is recorded separately and is not inferred from
neural oscillation or activation error.

## Diagnosed and corrected during this campaign

1. The first launcher completed a neural replay but rejected the response because
   the worker protocol lacked its `ok` field. The protocol was corrected; the
   interrupted folder was retained.
2. BiW2Lvar+ A passed with 2.469/2.466 s periods and 0.000709-cycle bilateral phase
   error. A fresh-process replay reproduced the objective exactly.
3. The first B search, with gains 0.05–4, failed to accelerate the rhythm. Port
   names were checked against the builders and were correct. Individual Ia/II/Ib
   ablations showed that stronger Ia or II, unlike stronger Ib alone, entrained
   the rhythm. Those probes are preserved under
   `D:\Users\Ben Bolen\ControllerTraining_20261104\hunt_B_diagnosis_v2`.
4. The new B search includes those diagnostic candidates and logarithmic gains
   over 0.05–128. Its accepted BiW2Lvar+ candidate has 1.233/1.232 s periods, 9/10
   detected cycles, and 0.00669-cycle bilateral phase error. Its multiplier of 96
   scales the explicitly declared reference conductance budget; it is not a
   claim about biological firing rates or a measured human reflex gain. Source
   saturation, normalization, and the operating point must accompany that number.
5. Contact checking now flags any explicit inhibitory sensory entry, even to one
   channel, before traversing intrinsic RG/PF reciprocal inhibition. Toe-specific
   baseline paths, including the old local toe-dependent inhibition, are removed
   from these new candidates. Toe→ankle PF remains a separate proposed test.

## Isolated postural premotor test

Source: [postural pilot](../postural_pilot_20261104/README.md).
Records: `D:\Users\Ben Bolen\PosturalPilot_20261104`.

The descending task command passes through an antagonistic premotor pair and its
inhibitory relays before reaching MNs. Ia and II enter the premotor layer; Ib is
tested as off, inhibitory, or load-excitatory. A small local Ia→MN reflex remains.
The external task-error law is a declared PD law; the module tests the premotor
interface, not a complete neural implementation of task estimation.

The initial 10 kg fixed-base ankle bench recovered after an 8 N m, 50 ms push
with premotor Ib off or inhibitory. The load-excitatory Ib alternative failed the
final-angle gate. This result is specific to standing on this bench. It does not
select the Ib sign for locomotion. A small Rustuna fit of shared neural gains and
held-out tests are recorded separately. The circuit has not been deployed into
all four full controllers.

No existing tendon, damping, CAD, manuscript, or Overleaf configuration was
changed. The toy bench declares its own 0.2 N m s/rad hinge damping, 30 ms
actuator filter, and 5% excitation floor. It uses synthetic tensile actuators,
not a BPA pressure law. Training and physical support are separate quantities.

## Completed output fits and physical checks

The final ordered output study is
`D:\Users\Ben Bolen\ControllerTraining_20261104\ordered_output_fits_v5`.
All four families passed A and B. Their C peak errors / D waveform errors are
Hunt 0.222922 / 0.214819, s3k 0.122141 / 0.367688, w2lvar 0.121615 / 0.402723,
and syn6 0.120525 / 0.247142. BiW2Lvar+ and syn6 pass the offline D gate; s3k and
w2lvar fail and are not promoted. C now freezes timing, disables local sensory
readout for the peak fit, and fits one bounded gain per muscle base, shared
across sides. D restores local feedback and refines Ia/Ib with Rustuna.
BiW2Lvar+ replay includes reference-pose geometry-dependent recruitment; the initial
ungated replay was a transfer inconsistency and is preserved as a diagnostic.

`hunt_pose_air_validation_v3\COMPLETE.json` records protected 16 s physical runs
of the corrected BiW2Lvar+ fit. With the same weights and PF expression 0.1, stage 1
has 4/4 principal knee cycles, 2.466 s period, and phase error 0.03367; stage 2
has 7/7 cycles, 1.300 s period, and phase error 0.09328. All knees remain flexed.
Stage-2 measured spring support is 0.998602 BW, with ground disabled and no foot
load. It is suspended stepping. The first full-expression and higher-expression
tests stopped before predicted hyperextension and remain in their own folders.

The original stage-1 counter counted secondary flexion humps twice. A separate
principal-cycle audit, with antiphase/synchronous/fast/frozen self-checks, passes
the physical motion while retaining the legacy failed assessment. Stage 2 passes
both assessments. No production metric or dissertation checkpoint was overwritten.

Strict contact fault-injection checks pass all four families, including a single
inhibitory entry and a toe-to-RG entry. Exact results:
`D:\Users\Ben Bolen\ControllerTraining_20261104\strict_contact_contract_v1\RESULTS.json`.

The portable defense update, new GIF, traces, schematic, and full model paths are
in `Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104`.
The postural delay/noise test found that the initial quiet fit loses robustness
at 30 ms added delay. The first separate robust fit is stored in
`D:\Users\Ben Bolen\PosturalPilot_20261104\robust_delay_fit_v3`; that fit passed
6/6 training but only 9/18 unseen cases. A follow-up, robust_delay_fit_v4, passed
6/6 training and 17/18 unseen cases. Its remaining 70 ms-delay case settled at
3.125 s against the retained 3 s gate. Neither fit is promoted to all families.
This batch uses graded controllers; spiking
implementation work is not included in these new acceptance results.

Syn6 passed protected physical stage 1 with 3/4 cycles, 2.460 s period, and
0.00373 phase error, under 99.876% BW spring support. In stage 2, a legacy
distance-based knee counter passed fast multi-hump activity that the spectral
guard rejected. Recorded spectra place ankle MN/joint activity near 5.2 Hz in
one diagnostic, while RG stays near 0.44–0.56 Hz. Reducing RG sensory gain did
not solve it. Lower local II gain removed the fast component but left knee
excursion below the 20° gate. Actual-body tuning of local II and expression is
bounded to twelve Rustuna trials in syn6_body_fit_v6, with a fixed cadence
criterion and fresh-process finalist checks. It completed without meeting the
cadence gate. Its slow candidate (expression 0.2913677891, local II 0.2114065546)
passed principal air-motion checks again in syn6_slow_candidate_repeat_v7:
stage 1/2 periods 2.466/2.004 s, with stage-2 support 99.867% BW. The subsequent
STAGE2_SPECTRA.json detects substantial secondary ankle activity at 6.4–7 Hz.
That waveform-quality defect remains; this is not a fully accepted finalist.
The syn6_drive_response_v2 map shows that stronger equal F/E commands speed
the isolated network but synchronize the sides above command 0.35. No mechanical
property was changed.

## Layer stimuli and fixed-parameter noise

hunt_layer_stimulus_2nA_v2 and hunt_active_phase_stimulus_v3 record 100 ms
external-current pulses at RG-F, hip-flexor PF, and iliacus MN, plus an inhibitory
PF pulse. Pulses are at 7.1 s (PF silent) or 7.65 s (PF active). All cases
completed with protected knees. The recorded pre-pulse trajectories match
exactly. The silent-phase PF pulse stays below its outgoing synapse threshold;
the active-phase pulse produces a small MN/mechanical response. MN perturbations
also change RG later through mechanical/sensory feedback. Actual pulse samples
and paired differences accompany the few-cycle plots, not inferred markers.

hunt_rg_noise_v1 tests fixed learned parameters with independent 20 ms
Ornstein–Uhlenbeck current at the four RG cells, stationary sigma 0.01/0.025 nA,
three seeds per amplitude, and a fresh-process repeated seed. This evaluates air
coordination under ongoing noise; it does not change synaptic or intrinsic
parameters and does not establish standing robustness. Results are saved in
STATUS.json and COMPLETE.json under the external training root.

That fixed-weight noise batch completed: all six noise cases plus the zero
control and repeated noisy seed passed (8/8 total). The fresh-process noisy
repeat matched qpos, recorded neural voltages, and actual noise exactly, with
maximum absolute differences 0.0. Knees remained flexed in every case. This is
a small suspended-air robustness sample; no noisy case was used for fitting.

## Batch handoff

All bounded jobs launched in this calculation pass have completed. No ongoing
training, full-body postural result, new ground-curriculum pass, or new spiking
result is implied. Next work is PF waveform design for s3k/w2lvar, syn6 ankle
feedback/cadence repair, and a separately tested full-body postural integration.
The selected BiW2Lvar+ air GIF and aligned layer/joint/load plots are in the repo
handoff. Exact compiled models and raw results remain on easteregg2.

An external copy of this record is
`D:\Users\Ben Bolen\ControllerTraining_20261104\EXPERIMENT_RECORD.md`.
The final source snapshot and hash manifest are in
`D:\Users\Ben Bolen\ControllerTraining_20261104\source_snapshot_final`.

Canonical display name: **BiW2Lvar+**; plain identifier: **BiW2LvarX**.
Saved `hunt` strings and path names are deprecated compatibility identifiers,
not a model attribution. Hunt's four-step calibration and Bolen's five-stage
curriculum are distinct parts of the workflow; the revised candidate has so far
passed physical curriculum stages 1 and 2, not all five stages.


## s3k/syn6 PF readout and reflex follow-up (Session 1)

Completed: mirrored per-existing-edge PF-to-MN current-balance fitting,
32-trial local Ia/Ib/II refinements per new family fit, s3k relay-output causal
comparisons, syn6 live ankle-source comparisons, and bounded physical air probes.
Detailed record: S3K_SYN6_FOLLOWUP_20261104.md in this directory.

s3k C/D normalized waveform RMSE now 0.168362/0.173783; syn6 C/D
0.161525/0.158450. Both D finalists match in fresh processes. s3k needed
independent calibration of the Ib relay output, which is also driven by RG-E;
reducing raw afferent input alone did not suppress that output. No new PF cells
were added. Output weights and sensory gains remain mirrored across sides.

The separate syn6 ankle candidate retains its older D checkpoint, attenuates
local ankle Ia/II outputs by 0.15 while preserving knee-crossing sources, and
uses expression0.3. Both air stages pass at fixed settings (periods2.466/2.005s).
Stage2 repeats qpos, qvel, neural voltages, foot forces, and harness wrench
exactly. Absolute >3Hz ankle RMS improves4.52/4.84 to2.70/3.13deg; residual
ripple remains. Foot force is zero in air. This is not minimum BW support for
ground walking. The newly fitted syn6 readout has not passed physical stage2.

s3k physical expression0.2 passes stage1 but stage2 hip ranges127.6/142.4deg
exceed the80deg gate. Two hip-only gain reductions also reduce knee movement
and fail stage1. All failures retained; no ground advancement or physics edits.
All bounded follow-up jobs completed. Raw results and final source snapshot:
D:\Users\Ben Bolen\ControllerTraining_20261104\source_snapshot_s3k_syn6_followup.
Portable results/GIF/plots: Dissertation\Defense\Session1_Handoff_20261104\
experiments_20261104\s3k_syn6_followup. No commits or manuscript edits.

Visualization correction (2026-11-05): the air-stage GIF previously displayed
an inactive ground plane, so feet appeared to pass through a floor. The renderer
now hides that plane and states ground contact is disabled. Saved qpos, forces,
controller parameters, and acceptance results are unchanged. The earlier source
snapshot remains historical; the repository renderer contains this later fix.

## Physical IK response scoring (2026-11-05, Session 1)

kinematic_response.py now scores actual bilateral joint cadence, ROM gain,
mean bias, uncentered waveform error, and three harmonic gain/phase ratios.
body_probe.py automatically writes body.kinematics.json; new fit_body_feedback
studies select database_kinematics_v1 and freeze metric/reference hashes.
Old study objectives and saved results remain unchanged. The inherited
kine_ref.py already had several of these comparisons; the newer fit path had
not used their full content as its physical optimization objective.

Saved s3k and syn6 runs were scored against the original calibration IK/GRF
reference and all43 local reference files. Neither matches the primary reference
adequately: s3k hip ROM gain2.88/2.94; syn6 knee ROM gain0.29/0.32, ankle mean
bias-25.5/-26.2deg, and periods2.004/2.008s vs1.2334/1.2133s. Air-motion passes
are retained separately; no claim of database waveform validation follows.

Known amplitude2, lag45deg, offset, phase reversal, and half-cadence metric
checks pass. A fresh16s physical integration generated the metric and reproduced
joint positions/velocities, neural voltages, foot forces, and harness wrench
exactly. No new training/optimization campaign was run in this scoring pass.
Harmonic output/reference plots are not Bode input-output transfer functions.
Design/weights/limitations/full paths: KINEMATIC_OBJECTIVE_20261105.md here.
Portable figures/metrics: Dissertation\Defense\Session1_Handoff_20261104\
experiments_20261104\kinematic_response_20261105.

## Air-versus-ground scope correction (2026-11-05)

Ben correctly questioned grounded IK as an air-step target. Current default is
stage_specific_kinematics_v2: during air stages, database gain/phase/amplitude
and cadence are diagnostic only. Physical air ranking uses phase stability,
period CV, and cycle repeatability under existing motion/safety gates. Ground
reference differences cannot change air acceptance or ranking (self-check).
Hunt et al2017 section2.4.1 sets endogenous CPG PERIOD to twice intended stepping
period, not a universal2x physical air-step law. Our2.466s target was imposed
from2*1.233s, not independent confirmation. Compare intrinsic, afferented air,
and ground conditions separately; measure correlations in paired conditions.
The original v1 metrics remain historical; corrected saved-run reports are in
kinematic_response_20261105\air_scope_v2. Syn6 passes air gates; s3k still fails
its existing hip-motion gate. No training was run using the short-lived v1
objective, and no physics was changed. KINEMATIC_OBJECTIVE_20261105.md updated.
Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC5378996/ section2.4.

## Full goodness-of-fit statistics update (2026-11-05)

Primary result tables now show whole-window and per-cycle peak-to-peak amplitude,
min/max angles, sample SD and amplitude percent error; frequency/period mean+SD
and errors; per-cycle circular phase mean/SD and maximum absolute phase error;
and worst absolute waveform error across every retained cycle (plus per-cycle
maximum-error mean+SD). Joint spectral peaks and >3Hz power are also reported so
cadence cannot hide ankle ringing. RMSE is supplemental. Ground loss includes
normalized worst-cycle maximum error and per-cycle phase error. Ground IK
remains diagnostic for air; no air acceptance/ranking change. Known amplitude,
phase wrap, delay, and one-cycle100deg spike checks pass. No new training in
this update. Main portable RESULT.md now contains both models' full statistics,
not just a methods narrative. See air_scope_v2\{s3k,syn6}\RESULTS.md and JSON.

## Circuit beside neural traces (2026-11-05)

Created expanded s3k/syn6 offline-D figures with left RG/PF/MN schematic,
primary-group PF-to-MN matrix, and right three-cycle RG/all-PF/MN/target traces.
Separate exact individual-muscle maps, direct Ia/II/Ib-to-PF maps, and iliacus
feedback subnetworks are supplied as PNG/PDF/editable SVG. The extracted circuits
reproduce saved RG/PF/MN/target arrays exactly. All PF-to-MN homologues match
bilaterally with max difference0.0uS and no missing counterparts.

Flexor Ib is present locally: s3k directly inhibits the homonymous MN (0.035uS
for iliacus/bifemsh/tib_ant); syn6 uses Ib->IBIN excitatory0.0609156uS then
IBIN->MN inhibitory0.5uS. These checked flexors have no Ib->RG/PF pathway.
Central Ib sources are selected stance extensors; Ia/II centrally follow actual
muscle-carrying PF projections. No wiring was changed for the figures.

S3k uses four shared phase channels, not separate joint PF oscillators. Syn6
RG-E drivesS2/S3/S6, RG-F drivesS1/S4, both driveS5. Exact muscle allocations
are in PF_ASSIGNMENTS.md/CSV. Syn6 figure uses the new offline D fit, not the
older-readout air GIF candidate. Do not mix their weights.
Portable README/figures:
Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\
experiments_20261104\circuit_trace_20261105\README.md.
Full ledgers remain external under:
D:\Users\Ben Bolen\ControllerTraining_20261104\circuit_trace_audit_20261105.

## PF design failure confirmed; downstream fitting held (2026-11-05)

Ben identified hidden diagram routes and inadequate PF timing. Orthogonal A*
routing now validates every synapse against all neuron interiors: s3k20overview/
12flexor edges and syn627overview/20flexor edges, all zero interior crossings.
Neural data unchanged. PF timing audit confirms syn6 actual99%-variance rank3
vs intended6 on each side; S5 active100% vs intended~32.5/33%. RightS1/S2/S4
lags-84/-100/-61deg. All12 channels fail the provisional combined timing check.
S3k four channels have99%-variance rank2, with no declared numerical window
targets; its PF_SHAPE adaptation multipliers are discarded by _build_pf.

A bounded8-evaluation PF-only current-balance fit, with RG/MN weights fixed,
uses mirrored right primitive targets on the left with measured loading offset.
Loss2.51171->2.15643; fresh repeat exact, still0/12pass. Earlier diagnostic used
original independently factorized left targets and is retained separately.
Neither is installed into the body adapter or accepted as a controller repair.

pf_design_gate.py now blocks new C/D, readout, and body optimization for these
known failed PF revisions through the current fitting entry points. Historical
replay/figures and isolated PF-repair experiments remain available. No physics,
CAD, manuscript, or old checkpoint changes. Period/antiphase and MN-output
passes do not constitute PF timing validation.
Full source/report: controller_profiles_20261104\PF_DESIGN_AUDIT_20261105.md.
Portable report, corrected figures, timing overlays:
Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104\
circuit_trace_20261105\PF_DESIGN_AUDIT_20261105.md.
Raw repair results: D:\Users\Ben Bolen\ControllerTraining_20261104\
syn6_PF_input_repair_mirrored_20261105. All bounded experiments completed.
# Native syn6 onset and shutoff, 2026-11-05

The opt-in syn6 build adds bilateral graded SNS trigger, clock, onset,
shutoff, and quiet-shunt cells. The two clocks per side reset at the RG
extensor or flexor onset. All relay integration uses the compact SNS backend
at 2 ms. The new circuit replaces the legacy RG-to-PF projections and
PF-to-PF-interneuron drive. The sensory projections remain active, including
Ia and II, and Ib under its existing load gate. The suspended stage does not
deliver heel or toe events. PF-to-MN projections and all other fitted MN
input synapses have identical named endpoints, conductances, activation
ranges, and reversal potentials in the replay and live comparisons.

Final experiment on easteregg2 is in
`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_native_20261105_v7`.
See `timing\SPEC.json`, `timing\TIMING.json`, `final_verified\PF_TIMING.png`,
and `live_timing_knee-25\TIMING.json` and `live_timing_knee-35\TIMING.json`.
The legacy closed-loop comparison is in
`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_native_20261105_v5\baseline_body`.

| Measured replay period | Passed L | Passed R | Bilateral channels passing all four strides |
| --- | ---: | ---: | --- |
| 1.048 s | 12/24 | 4/24 | S1 |
| 1.233 s | 24/24 | 24/24 | S1–S6 |
| 1.418 s | 4/24 | 4/24 | S1 |

The cadence stress test rescales the legacy neuronal integration and RG
inactivation time constants and the afferent replay cadence together. The
native relay time constants, thresholds, and conductances stay fixed. Timing
is relative to each leg's RG extensor onset, with the target origin locked to
the nominal reference. The nominal worst absolute edge error was 4.5 percentage
points and the worst fundamental phase error was 14.81 degrees. All three
recorded neural traces repeated byte-for-byte in a fresh process.

The native Stage 2 Gait2392 probe with live Ia, II, and Ib signals hit the
predictive knee guard at 0.332 s with the initial 10-degree knee flexion.
At 25 and 35 degrees of initial flexion, the final native candidate hit the
guard at 1.254 and 1.686 s. Each longer capture contains one measurable left RG
interval with zero of six PF timing passes and no complete right interval.
The native candidate remains experimental and does not clear the live-body
or cadence-range gates.

## RG-cadence adaptation and live knee-guard diagnosis

The opt-in adaptive configuration estimates four E/F recurrence intervals
from the actual bilateral RG voltages. It scales each native clock's
integration rate inversely with its latest valid period. It starts with the
legacy PF timing inputs and hands over only after all four clocks have
measured a period. The scheduler changes PF timing edges and clock rates.
Native SNS neurons still perform the onset and shutoff integration. Frozen
MN synapse digests match the original replay and all live runs. Ia, II, and
load-gated Ib remain enabled, and the knee guard is unchanged.

Run directory on easteregg2:
`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_adaptive_20261105`.
`replay\TIMING.json` and `replay_live_range\TIMING.json` record 48/48
passes at each of six measured periods spanning 0.863–1.418 s. These
are cadence-controlled replays with sensory input. See `replay\PF_TIMING.png`.

In the original matched 10-degree flexion startup, right bifemsh MN
voltage was about 1.3 mV through the first extension. Fixed-clock native
timing dropped it below rest while the vastus intermedius MN rose to about
0.6 mV, and peak right knee extension speed rose from approximately
368 to 770 deg/s. After adaptive handover, the same pattern appeared on
the left with sustained S5 and rising S2 during the extension transient.
This is a measured loss of flexor recruitment relative to extensor
recruitment. It does not isolate individual muscle or passive joint torques.

| Adaptive live Stage 2 condition | Guard time | Side |
| --- | ---: | --- |
| One-cycle handover, 10-degree start | 3.582 s | Left |
| One-cycle handover, 25-degree start | 3.316 s | Right |
| Three-cycle handover, 10-degree start | 9.814 s | Left |
| Three-cycle handover, 25-degree start | 4.724 s | Left |
| Three-cycle handover with flexor phase correction, 10-degree start | 5.572 s | Left |

The longest run repeated byte-for-byte in a fresh process and spent 4.17 s
at full native timing before the guard. Its fully native capture passed
11/30 timing checks across five complete RG intervals
(`live_full_native\TIMING.json`). The measured left E-to-F transition near
the stop was about 45% of the cycle versus 57% in nominal replay. A separate
variant shifted F-triggered thresholds to compensate for this duty change.
It stopped during the handover at 98.1% native and is not selected.
No adaptive live Stage 2 trial passed.

## Syn6 knee-torque audit and flexor handover

Results on easteregg2 are in
`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_flexor_handover_20261105`.
The torque audit restores measured state, activation, and control into the
saved MuJoCo model and recomputes actuator force and equality-aware moment
arms. Positive moments extend the knee. The values are muscle contributions
and exclude inertial coupling, passive joint terms, and contact.

In the 118 ms before the original fixed-native right-knee stop, rectus
femoris averaged 11.71 N·m extension, vastus lateralis 4.27, vastus
intermedius 2.79, and vastus medialis 2.66. Right semimembranosus and both
biceps femoris heads fell from a combined 17.35 N·m flexion in the matched
legacy time window to 1.00 N·m with fixed-native PF. Net muscle extension
rose from 12.85 to 20.01 N·m despite lower rectus femoris extension.
In the adaptive left-knee failure, rectus femoris was again the largest
extensor at 10.58 N·m, with the three vasti totaling 6.58 N·m.

S3 strongly recruits semimembranosus and both biceps femoris heads.
S2 strongly recruits rectus femoris and the vasti. The first candidate
extended the S3 timing window across S2. The later candidates retained the
original adaptive windows and added bilateral joint-state premotor relays.
The relay input rises when the knee is extending toward neutral. It excites
S3, suppresses the S3 quiet shunt, and inhibits S2. Neither intervention
changes the PF-to-MN or other fitted MN synapses. All compared live motor
digests match. Muscle afferents remain enabled and grounded diagnostics
use separately measured hindfoot events and whole-foot load.

| PF intervention | Stage 2 | Stage 3 | Stage 4 |
| --- | --- | --- | --- |
| S3 overlap | Guard at 4.912 s | Completed, failed load gate | Guard at 7.928 s |
| Joint-state relay, 65 ms prediction | Guard at 5.416 s | Guard at 5.586 s | Guard at 6.358 s |
| Earlier joint-state relay, 120 ms | Guard at 12.002 s | Guard at 5.536 s | Guard at 13.126 s |
| Joint-state relay, 150 ms | Completed, failed rhythm gate | Completed, failed standing gate | Completed, failed load and movement gates |

The 120 ms candidate ran approximately 6.36 s fully native in Stage 2.
At the left guard, its short-head biceps femoris supplied about 3.00 N·m
flexion versus 0.11 N·m at the earlier adaptive stop, but the candidate still
hit the predictive limit. Stage 4 used 0.65 BW explicit upward support plus
0.173 BW mean vertical spring support in the 120 ms run. The feet contacted
the ground and averaged 0.177 BW combined vertical load over the captured
interval. The run stopped before meeting the supported-stepping criterion.

The 150 ms candidate completed all three evaluations with the predictive
knee guard enabled throughout. Stage 2 ran 16 s with 12 right and 11 left
knee cycles, a median cycle period of 0.796 s, and a period CV of 0.287.
Maximum knee extensions were −1.40° and −1.45°. After full handover the
PF timing gate passed 8/24 checks per side. The compact run used 42.87 s
wall time including setup and recording. Its aligned trace repeated
byte-for-byte in a fresh process.

Stage 3 ran 8 s. Feet carried approximately 0.24 BW, spring support
averaged 0.76 BW, torso tilt was 17.7°, and sway was 0.140 m. Stage 4 ran
16 s with 0.65 BW explicit support and approximately 0.18 BW mean spring
support in the scored window. The force gate counted six right and eight
left contact cycles. The mesh clearance gate counted three right and
five left swings, with peaks of 26.7 and 54.3 mm. Mean combined foot load
was 0.164 BW, phase error was 0.213 cycles, and the movement gate failed.
The grounded tests use the same frozen circuit with no added standing
or walking posture controller.

## Stage 2 PF regularity

Results on easteregg2 are in
`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_stage2_regularize_20261105`.
All trials retained the frozen MN synapses, live afferents, and predictive
knee guard. Stage 4 weaning was held. The reference 150 ms case showed one
or two S2 onsets and three or four S3 onsets per leg per complete RG cycle
in the 7–14 s comparison window. The extra S3 onsets coincided with
repeated knee-protection relay activation.

| Change | Stage 2 outcome | Knee period | Period CV | Bilateral phase error |
| --- | --- | ---: | ---: | ---: |
| Reference 150 ms | Completed, rejected | 0.796 s | 0.287 | 0.138 cycles |
| Shared same-leg clock, smoothing 0.35 | Guard at 5.374 s | — | — | — |
| Separate clocks, smoothing 0.65 | Completed, rejected | 0.678 s | 0.364 | 0.153 cycles |
| Five-cycle handover | Completed, rejected | 0.788 s | 0.344 | 0.278 cycles |
| 70 ms premotor release | Guard at 4.060 s | — | — | — |
| 120 ms premotor release | Completed, rejected | 0.520 s | 0.144 | 0.207 cycles |
| S2 latch until next ipsilateral RG extensor onset | Completed, rejected | 0.752 s | 0.235 | 0.203 cycles |

The S2 latch gave exactly one S2 onset in all 11 complete bilateral RG
intervals counted in the figure window. S3 still opened three or four
times per RG cycle, and the full PF timing score remained 16/48.
The compact S2 latch run took 28.31 s for 16 s of simulation. Figures
and per-joint range, phase, and maximum-error comparisons are in
`latch_figures`. The time trace spans three right RG cycles. The overlay
uses a shared right-knee maximum-flexion anchor for all six joint traces
and compares suspended motion to Subject01 and Ong 1.25 m/s ground
walking IK. Foot load was zero throughout the air trials.

## S3 single-burst Stage 2 diagnostics

Results are in
`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_s3_single_20261105`.
In the S2-latch recording, the knee-protection encoder repeatedly rose as
the knee re-entered extension. Each activation both excited S3 and
suppressed its quiet shunt. Representative right S3 onsets followed the
knee relay by about 6 ms. The reflex responds to the projected knee
angle and extension velocity every step, independently of RG phase.

After handover the single-window candidate routes S3 through a dedicated
RG-timed relay. It removes the repeated knee-relay drive into S3 while
retaining the knee relay's inhibitory S2 exit, the existing afferents,
the frozen motor digest, and the predictive guard.

| S3 drive window, fraction of ipsilateral RG cycle | Guard time |
| --- | ---: |
| 0.26–0.41 | 6.450 s |
| 0.20–0.65 | 8.050 s |
| 0.16–0.88 | 9.084 s |
| 0.40–0.88 | 6.890 s |

The widest window gave one S3 onset in each of five complete, fully native
RG intervals before the guard. Its measured onsets ranged from 12% to 24%
of the cycle, earlier than the nominal target near 30%. The target-length
window stopped before producing a complete fully native interval. No
candidate met the single-burst, timing, and knee-safety requirements together.

`SYN6_WIRING.html` lists S1–S6 separately, including left and right MN
destinations, frozen synapse values, other PF exits, and exact central
afferent routes. `phase_labeled_overlay\KINEMATIC_OVERLAY.png` adds separate
stance/swing bands for both benchmarks and both legs to the completed
S2-latch comparison. Both simulated feet remained unloaded in that air run.

The editable drawing is `SYN6_CONNECTOME_EDITOR.html`, with
`SYN6_MOTOR_AFFERENTS.json` for importing into the existing editor.
The navigation is leg → PF channel → muscle → off-page feedback targets.
Motor pools shared by several PF channels appear as references on each
relevant page. The editor's materializer preserved all nested edges.
`left_right_figures` has the corrected bilateral order and identifies
the Ong et al. benchmark as a 1.25 m/s walking simulation.
