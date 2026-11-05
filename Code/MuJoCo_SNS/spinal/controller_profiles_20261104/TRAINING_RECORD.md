# Session 1: four-command controller experiments

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
