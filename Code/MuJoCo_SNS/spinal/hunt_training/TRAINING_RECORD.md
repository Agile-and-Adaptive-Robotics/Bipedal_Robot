# Gait2392: Hunt four-step training and five-stage validation

Machine: easteregg2. Repository: `D:\GitHub\Bipedal_Robot`.
Python: `D:\Anaconda\envs\myo\python.exe`.

This is the working experiment record, separate from `CHATGPT_REPORT.md`.
Raw arrays, trial configurations, logs, and databases belong in the ignored
`runs/` directory. Existing production winners and the completed exploratory
campaign remain separate.

Current Session 1 checkpoint (2026-10-30): the tonic-abductor candidate passes
stages 1 and 2 with the same saved parameters in fresh-process repeats, and
passes assisted stage-3 standing. No fit passes the complete corrected
curriculum. Actual alternating ground stepping remains unresolved.
See the latest Session 1 entries below; `FOURCELL_HUNT_RESULTS.md` describes
an earlier comparison rather than the current candidate.
Sartorius and gracilis remain separate. Stage 5 is withheld until actual
alternating load-bearing ground steps pass.

## Current sequence

1. Audit full IK sources and model constraints for MTP, lumbar x/y/z, hip
   x/y, and calcaneal/subtalar x-related motion before selecting PF targets.
2. Prepare model-consistent kinematics, forces, expected afferents, and MN targets.
3. Fit undriven intrinsic rhythm, expected-afferent transition timing,
   PF-to-MN peaks, and Ia/Ib refinement as four separate functional steps.
4. Run deafferented air, afferented air, standing, harness-supported ground
   walking, and ground contact-PF refinement with a measured support-reduction
   ladder. Failed stage gates remain failures; later diagnostic tests do not
   certify progression. Stage 4 includes ground contact, as the stock curriculum
   specifies. The previous w2lvar-based air-pattern stage is not reused as
   though it were supported walking.

## Initial IK finding

The 43 cached gait references contain only bilateral hip flexion, knee, and
ankle curves plus summary statistics. Their intake code deliberately discards
the other coordinates. These caches cannot answer whether non-sagittal,
lumbar, or MTP patterns exist. Full original MOT/STO files and model locking
must be inspected; omission from a cache is not evidence of absent movement.

Kinematic independence can justify testing a separate temporal control pattern.
It cannot by itself identify a biological PF population or prove that one
dedicated oscillator per coordinate is necessary. The Gait2392 subtalar joint
is oblique; its coordinate and anatomical-x component will be distinguished.

## Completed intrinsic training: Hunt step A

The primary experimental reference has a 1.233 s right heel-strike period.
The target undriven RG period is therefore 2.466 s. The previous intrinsic
parameters latch without descending drive, so changing only the time constant
cannot satisfy this step. A scan of intrinsic activation/inactivation thresholds
and inhibitory coupling selected a self-sustained regime. Period fitting then
selected tau_h = 0.631371699 s, NaP activation threshold e_m = 0 mV,
inactivation threshold e_h = 2 mV, and laminated inhibitory gain 6.

The native SNS replay produces 2.466 s cycles, period CV 0.000313, and extensor
duty 0.5069. A finite initialization pulse ends at 100 ms; there is no tonic
input during the scored 10-40 s interval. The preallocated Numba training kernel
agrees with native SNS within 7.4e-13 mV and 2.8e-15 in the h gate. This is
equation agreement to a measured tolerance, not a bit-exact claim. MuJoCo
validation will retain the previously verified exact sparse backend and 2 ms step.

Evidence: ignored `runs/step_A/winner.json`, `regime_scan.json`, and
`validation.npz`. Other functional fits and the new physical curriculum are not
yet complete.

## Forward-kinematics sign audit

Compared a +0.1 rad coordinate perturbation in source OpenSim and converted
MuJoCo using the resulting body rotation matrices. All 17 axes align or oppose
with absolute dot product >0.99. Correct source-to-MuJoCo signs are -1 for
bilateral hip flexion and hip adduction, and +1 for hip axial rotation and the
other leg/lumbar coordinates. These differ from the old muscle-length-based
sign guesses. Future training and kinematic scoring use this explicit adapter;
production baselines and historical results have not been silently rewritten.

## Functional fitting and corrected physical tests

Hunt B independently calibrated both legs with intrinsic parameters frozen:
right period 1.236 s (target 1.2334), left 1.212 s (target 1.2133). Each passed
eight one-percent input-amplitude probes. Right-to-left transfer failed before
left-specific calibration; this is not independent-subject generalization.

Hunt C and D now have separate joint-PF, PF-to-MN, and Ia/Ib conductance fits.
Native expected-afferent replay has activation RMSE 0.127291 on a 0–1 scale.
Offline regression uses teacher-forced voltages; runtime and neural validation
do not clamp PF or MN voltages. Targets are inverse-force estimates, not EMG.
The experimental 1,021-neuron/420-input controller includes 17 joint-axis PF
pairs, causal phase filters, and biarticular integration. It is an engineering
adaptation, not a literal Hunt/Deng replication. No commissural coupling is
enabled initially; finite opposing initialization and body feedback coordinate
the legs. The fit does not certify successful closed-body control.

Moving the finite 100 ms startup pulse from 2 s to initialization, matching
the validated intrinsic fit, changes stage 1 from a symmetric latch to a passing
motion/frontal screen: 4/3 knee cycles, 44.4/29.2 degree excursions, period
2.464 s. This is an initial test, not a robustness certificate. Stage 2 remains
unsatisfactory. Ground stage 4 is genuinely contact-enabled. Initial height
tests load the feet but pin the left knee near its extension stop. Motor-block
adjustments restore bilateral cycles in several diagnostic tests, but loading
and other gates still fail. Do not claim curriculum completion or walking.

Full pelvis-restraint forces/torques, individual foot GRF vectors, and activation
effort are logged in ignored runs/corrected_curriculum_v*. Actual support is
measured from spring/damper force, never inferred from stiffness scale. Extra
out-of-plane joint restrictions remain declared assistance. The legacy runner
double-counts each foot contact as heel plus toe; the new evaluator independently
counts each contact once. Historical raw results and production code are intact.
See PHYSICAL_RESULTS.md for the current case table and measurement limitations.

Mean activation/excitation squared is a dimensionless secondary effort cost,
not metabolic watts. Valid motion, ground loading, and separation take priority.
Relevant constrained-effort precedent: https://doi.org/10.1371/journal.pone.0180320.
Speed-dependent pathway weighting requires separate contact/Ib/Ia/II ablations,
not a preset monotonic rule; sources and scope are in PHYSICAL_RESULTS.md.

## Standing-balance literature and instrumentation

BALANCE_DESIGN.md distinguishes quiet standing, postural recovery, and dynamic
foot-placement control. It reviews Peterka (2002), Ting and Macpherson (2005),
Welch and Ting (2008/2009), and Hof et al. (2007), with primary-source links.
The proposed delayed CoM task-level recruitment pathway is not implemented yet;
do not label the inherited standing pathways a validated Peterka/Ting controller.
At a 60 mm lower harness, measured foot loads sum to about 0.55 BW, nearly
symmetrically, but the positional balance check still fails. Upright posture
under assistance is not independent standing balance.

## Completed corrected-environment diagnostic batch

All five corrected environments have now been exercised, but only the initial
stage-1 screen passes. Stages 2–5 remain failed, so the curriculum has not
graduated. At the 80 mm lower stage-4 setting, both knees have four detected
cycles and the frontal ordering remains uncrossed, but phase/contact criteria
still fail. Contact-PF activation in stage 5 suppresses left-leg cycling in
these trials. Reducing out-of-plane restraint, quartering vertical stiffness,
and a 20 N/100 ms forward pulse were tested diagnostically, not promoted as
robust walkers. Quartering vertical stiffness barely changes measured support
(about 0.666 to 0.665 BW): stiffness reduction is not support-force reduction.

One half-restraint native process exited with code 3221225477 during setup;
its original log was retained, and a fresh-process retry completed. It is not
counted as a successful trial. Current per-case results are PHYSICAL_RESULTS.md.
HARNESS_DIAGNOSTICS.png shows simultaneous restraint forces, foot loads, and
knee trajectories from a failed stage-4 diagnostic, not accepted gait.

A 1,000-step comparison of this fitted 1,021-neuron network against native SNS
has maximum voltage difference 2.22e-15 mV; it is NOT bit-exact. Measured sparse
steps including conductance updates take 1.06 s versus 19.68 s for native
forward evaluation in that comparison. This is a neural-kernel observation,
not a whole-simulation speed ratio or a harness-free locomotion claim.

Next corrective work is the failed sensory/motor coupling and task-level
standing balance, not reducing assistance on an already-certified walker.
Keep the intrinsic fit fixed unless a dedicated revalidation justifies change;
any revised sensory timing must be followed by motor/reflex refitting.

## 2026-10-29: task-level balance implementation

Implemented the Peterka/Ting-inspired pathway in balance.py: four directional
task populations plus four inhibitory relays, bringing this opt-in controller
to 1,029 neurons and 424 inputs. Delayed CoM displacement/velocity and separately
optional filtered acceleration drive fixed, pose-calibrated muscle recruitment.
Somatic/orientation cue weighting is adjustable but fixed; this is not adaptive
reweighting or a literal receptor/circuit model. Native MN outputs remain free.

Initial replacement-mode gains did not improve standing. A neutral standing
configuration (hip -7 degrees, knee -10 degrees, ankle zero, in the native
coordinates) avoids using a mean walking posture as the standing initialization.
The inherited controller then passes the initial stage-3 screen: maximum
CoM-reference distance 0.04736 m. It fails when inherited balance feedback is
disabled. Additive task-feedback cases also pass, but their measured horizontal
harness reactions are slightly higher than the matched disabled-task control;
an incremental benefit has not yet been demonstrated. Do not attribute the
initialization benefit to the new controller or claim independent balance.

Campaigns balance_20261029_v1/v2 replace inherited BAL inputs; v3 explicitly
retains them. Source hashes distinguish revisions. A population-count setup
assertion was fixed by handling W2L's logging aliases separately from true new
neurons. Early setup failures are not successful trials.

The first matched 8 s full-simulation comparison measured 13.253 s accelerated
versus 122.120 s native dense SNS, approximately 9.21 times faster end to end.
All saved qfull, control, CoM, and neural arrays matched exactly in that pair.
Its neural-step sections measured 4.150 s versus 110.484 s. Timings include the
same model and neural topology, initialization, physics, sensing, and logging;
background trial load varied. These are observed times, not a controlled
hardware speed guarantee. Further active-pathway comparisons are being saved.

The full 3D foot-force vector is logged, but current neural contact/timing
pathways still mainly receive scalar contact/load. Fore-aft and lateral shear,
contact locations/COP, and ankle-referenced moments are distinct candidate
signals. Never treat a scalar normal force as the complete 3D contact state.

Lineage label proposed for this experimental branch: biw2lvar-Hunt-PT. It denotes
W2L-variant ancestry, biarticular integration, 3D PF additions, Hunt fitting, and
the task-balance pathway; it is not an unchanged historical W2L controller.
Current per-case results, failed gates, and matched timing data are in
BALANCE_RESULTS.md. CHATGPT_REPORT.md and manuscripts are not changed.

### Completed active-pathway comparison and continuation

For the active additive controller (kp 0.5, kd 0.2, ka 0, 100 ms delay), the
matched 8 s run took 19.478 s accelerated versus 131.266 s native, a 6.74-fold
end-to-end reduction. Neural-update sections were 6.594 versus 118.075 s,
approximately 17.9-fold. qfull, full control, CoM, and recorded neural traces
matched exactly. Three matched pairs observed 6.74–9.21-fold end-to-end ratios
under differing background loads; do not represent them as controlled hardware
benchmarks or faster-than-real-time full simulation.

The active stage-3 test passes with about 0.0474 m maximum CoM-reference
distance, 11.0 degree maximum tilt, and approximately 0.615 BW through the feet.
Mean vertical spring support is about 0.381 BW. The initial forward/lateral
30 N, 100 ms held-out pulses and half additional out-of-plane restriction also
pass. These are limited assisted tests, not a basin/generalization certificate.
The matched disabled-task control has slightly lower positional error and
horizontal harness reaction (0.10337 versus 0.10382 BW), so no incremental
benefit of the new pathway is established. The saved active winner is the best
nonzero candidate, not better than every control. Keep the extension opt-in.

Air stage 1 reproduces its passing screen from the calibrated air start;
stage 2 still fails. Neutral-standing-start air probes are recorded separately.
Ground stage 4 has genuine verified loading and uncrossed leg ordering, but
4/1 knee cycles and absent alternating contact fail the gait gate. Stage 5
also fails. Later tests remain diagnostic because the curriculum has not
graduated through every earlier gate.

Additional 60% and 40% prescribed vertical-support tests remove the vertical
spring/damper and apply explicit constant forces. Actual measured vertical
spring assistance is then zero. Both tested levels fail standing/stepping
criteria: the 60% supported standing case has excessive loading asymmetry,
while the 40% supported case exceeds the positional error limit. This is an
actual load-support test, not a reinterpretation of spring stiffness.

Reproduction on easteregg2, repository D:\GitHub\Bipedal_Robot:
run D:\Anaconda\envs\myo\python.exe with the full path to
Code\MuJoCo_SNS\spinal\hunt_training\balance_additive_campaign.py,
then support_force_campaign.py and summarize_balance.py in that folder.
Completed cases are reused; changed hypotheses require a new campaign folder.

## Network drawings, stage-2 repair, and PF clarification

Generated compiled-connectivity overview, SNS-Toolbox subnetwork panels, and a
visualization-only connectome-editor export. Current combined standing snapshot
has 1,031 neurons and 424 inputs, including two experimental commissural relays
and the eight task-balance cells. The viewer is network_figures/index.html;
full raw graph/export files remain ignored under runs/network_combined_diagram.
Enabled conductances and disabled outputs are distinguished; gait-expression
gating at a recorded pose is not a walking result. Generic drawing shapes do
not replace the model's native neuron/channel dynamics.

STAGE2_WIRING_AUDIT.md records the SADb/master-rule comparison and four-way
causal isolation. Timing feedback alone gives bilateral motion with poor phase
and crossing; motor reflexes alone suppress left knee motion; both together
pin both knees. The offline Ib encoder had omitted the physical runner's 0.5
factor. Only Hunt D was refitted into runs/encoder_fixed_v1; the original
weights are preserved. That correction alone does not fix the physical failure.

Restoring weak Ia reciprocal inhibition, reducing knee-extensor output gain,
and adding explicit axis/abduction feedback produces passing stage-2 movement
and frontal screens. The initial candidate has 4/5 knee cycles, about 1.844 s
period, phase error 0.08595, and no sampled crossing. Fresh replay and four
independent one-percent perturbations of its saved configuration gains pass.
This is not a perturbation of every fitted synaptic coefficient. The period
is slower than the 1.233 s reference, so reference-cadence transfer is unresolved.

Candidate ablations are important: removing reciprocal inhibition, removing
axis feedback, or restoring original motor gains causes failure. Removing the
weak commissural pathway still passes in this candidate. Thus its absence was
a design omission worth testing, not a proven necessary cause of this failure.
The engineered positional feedback and output-gain changes must remain explicit.

The revised controller passes the standing gate. Subsequent real-ground trials
still fail, usually with the left stance leg continuously loaded. Lowering or
raising the harness, a single finite walk-start pulse, tonic crossed contact,
contact-onset feedback, and lateral harness softening have not yielded a valid
alternating ground walker in these probes. No stage-5 progression was certified.

Ben clarified that his five stages were intended to emulate Hunt's four
functional training phases. The audit now records the mapping rather than
treating them as unrelated validation environments. Physical screen success
does not imply the corresponding nominal reference-cadence fit transfers.

The label HUNT_BANK was misleading: these are added engineering RG time filters,
not neurons prescribed by Hunt. Viewer labels now say RG phase filter. Old IDs
are retained only for saved-weight compatibility. The PT computation is also
explicitly documented as body-level CoM feedback computed outside SNS, distinct
from local myostatic Ia/II feedback; it overlaps the inherited CoM/trunk control.

Measured human-Gait2392 moment arms show RF hip flexion/knee extension, hamstring
hip extension/knee flexion, sartorius hip flexion/abduction with weak knee
flexion, and gracilis adduction/knee flexion with a small median hip-extension
contribution along these sampled gait poses. These are not full-RoM conclusions
or CAD-robot results. See BALANCE_AND_PF_CLARIFICATION.md and the raw function
sign file for ranges and coordinate conventions.

Added Ben's proposed four-cell RF/HAM module as a separate default-off variant:
the two existing integrators each excite an inhibitory relay targeting the
opposite integrator. Gracilis and sartorius remain function-specific, not
wholesale RF followers. Three initial stage-2 probes at pairing gains 0.1, 0.3,
and 0.6 pass the movement/frontal screen; periods remain about 1.85–1.87 s.
No full Hunt refit or ground-validation claim is made for this topology yet.
The four-cell candidate and actual balance dataflow have separate viewer panels.

Zotero: the requested Peterka (2002), Welch/Ting (2008), and Hof et al. (2007)
references were imported and verified in BolenDissertation (NLPW6UV5), item keys
8BBKLJ7S, TE2K8L4S, and 2TVQUUHV. No dissertation source or CHATGPT_REPORT.md edits.

## Session 1: fresh four-cell training, 2026-10-29

Ben authorized training the RF/hamstring four-cell variant while retaining
separate sartorius and gracilis modules. `train_fourcell_hunt.py` revalidated
Hunt A/B, then fitted C in three frozen-timing blocks (joint PF, RF/HAM
integration inputs, and MN drive), followed by D with C frozen. The pairing
gain is 0.3. D now optionally fits Ia reciprocal-inhibition outputs alongside
Ia/Ib refinement; the production/default and earlier candidates are unchanged.
The C fits start with no previously fitted motor-reflex weights, avoiding
contamination of feedforward identification by the older D fit.

Fresh isolated A: 2.466 s, 12 scored cycles. Isolated B: right 1.236 s and
left 1.212 s. With weak commissural coupling and the paired topology, expected-
afferent replay gives 1.264/1.246 s before fitting C. These timing checks pass.
The new C fit identifies 816 joint-PF, eight RF/HAM input, and 932 MN edges;
D fits 674 available Ia/Ib/IaIN edges, including zero solutions. C-only replay
RMSE against its 80% activation target is 0.10003; C+D against the full target
is 0.12700. These are model-derived activation targets, not EMG. They do not
establish an additive 80/20 physiological contribution or closed-loop walking.

The full fit completes in 158.71 s on easteregg2 using sparse/preallocated
SNS calculations, cached edge indices, and a vectorized afferent encoder.
Raw weights, teacher traces, source snapshots, and hashes are preserved under
`runs/fourcell_hunt_20261029_v1/`. Earlier weights were not overwritten.

Initial physical tests fail air-stage phase/frontal gates. Assisted standing
passes, but this is a diagnostic later-stage result while earlier gates remain
unresolved. Real-ground trials and sagittally restricted/reflex-ablation
diagnostics are in progress; no curriculum graduation or support reduction
has been certified for the fresh fit.

### Completed follow-up and current outcome

The in-progress checkpoint above is superseded by the completed comparison in
`FOURCELL_HUNT_RESULTS.md`. A second C fit combines expected-afferent 1.233 s
and deafferented 2.466 s conditions before refitting D. Its deafferented
physical phase error improves to 0.0739, but smaller knee ROM is only 16.85
degrees and fails the unchanged 20-degree gate. Its afferented air test still
fails coordination. C+D activation replay RMSE is 0.13502; it is not a walker.

A third fit identifies the native four-cell input and inhibitory conductances
before refitting MN output and D. A preallocated four-neuron simulation agrees
with the full network's recorded module to 4.44e-16 mV. The 768 module trials
take 10.91 s, excluding physics; subsequent MN/D fitting and replay checks take
175.30 s. Paired voltage MSE improves from 0.8103 to 0.7629 on the right and
0.7043 to 0.6883 on the left. The physical air stage still fails, despite
improved deafferented phase (0.0667), because left knee ROM is 17.33 degrees.

All three fits pass assisted standing diagnostics. Real-ground tests use
actual contact and include gradual lowering, constant-speed towing, reduced
lateral restraint, intact extensor output, and timing/reflex ablations. The
strongest loaded third-fit diagnostic places 55.81% BW through the feet and
44.05% through the vertical harness, but both feet remain planted. This is
load sharing, not a successful walking support-wean. There is no harness-free
walking and no stage-5 graduation.

Recorded RGs keep bursting while the feet remain loaded. The first paired PF
fit is biased toward right RF and left HAM. Calling the physical stance
failure an arrested RG would be unsupported. The direct paired refit improves
the supervised neural target only modestly and does not solve ground release.
The diagnostic plot shows measured neural and force traces explicitly.

The ground acceptance gate now requires independently counted GRF unloading
and re-contact, as well as motion, frontal, and tilt gates. Legacy duplicated
heel/toe totals no longer determine acceptance. Synthetic alternating,
synchronous, permanent-double-stance, and airborne checks pass their expected
outcomes. A fresh loaded replay is bit-identical in full coordinates, control,
and recorded neural traces; the corrected gate still rejects it.

The fixed-conductance, 500-update neural benchmark is 0.2449 s accelerated
versus 24.5338 s native (100.17x), with voltage disagreement 3.55e-15 mV. This
excludes physics, conductance updates, setup, and logging; concurrent physical
work was running. The vectorized/scalar encoder disagreement is zero in the
sampled check. Earlier full-model speed measurements remain 6.7–9.2x end to end.

Raw results, weights, source snapshots, input hashes, and the curated machine-
readable summary stay in ignored `runs/fourcell_*` directories. The current
viewer is `network_figures/fourcell_dynamic/index.html`; the original entry
page redirects there. Its snapshot has 1,035 neurons and 424 inputs. The
separate records retain failed gates and the unpromoted passing calibration
probe. No dissertation source or CHATGPT_REPORT.md changes were made.

The final latest-conductance check, after the physical queue completed, gives
0.2530 s sparse versus 14.3064 s native for 500 neural updates (56.55x), with
maximum voltage discrepancy 4.34e-19 mV and zero sampled encoder discrepancy.
The earlier approximately 100x timing is retained with its concurrent-work
context. Both are neural-only checks; neither is a full-simulation speed claim.

## Session 1: swing-release diagnosis, 2026-10-30 (active)

Force/current instrumentation confirms continued RG bursts but about 9.33 N m
right RF knee-extension torque during commanded swing, predominantly PF-driven
(about 1.52 model-current units versus 0.002 from motor reflexes). This is a
measured antagonist conflict, not proof of a stalled RG. A fresh net-active-
torque PF fit passes afferented air without manual motor scaling (3/4 knee
cycles, phase error 0.0915, no sampled crossing). Stage 1 still fails phase;
no later curriculum stage has been advanced.

Exact old activation-target replay also fails: 15/15 double-frequency cycles
and poor alternation. A low-gain sagittal torque-allocation teacher produces
alternating air motion at 2.466 and 1.233 s. These are explicitly clock/feedback
diagnostic teachers, not neural passes. Their successful physical traces are
being captured to replace inadequate motor targets and refit expected-afferent
timing, C output, and D refinement in order.

The v3.3 editor export is in Neuromechanical_Models/Mujoco_SNS_models as
Hunt_fourcell_session1_*. It includes the exact 1,035-neuron graph, 19 native
SUB modules, and right/left RF-HAM views. Four Hunt Session1 templates are now
registered in the original editor under Templates > SNS walker versions.
Its standalone preloaded copy uses an isolated storage key. The prior embedded
library, canonical 110/216 connectome, and original storage key are preserved.
Twenty-three contexts round-trip through the actual JavaScript materialize/
specOf functions. Both editors fix outer-model retention when saving inside a
SUB. No TEST_* entries were created. Live file:// UI access was rejected again
after Chrome reconnected following Ben's permission update. No server or
alternate-surface workaround was run. Live UI verification remains incomplete.

### Feasible-target fit and failed body recheck

`runs/feasible_hunt_20261030_v1` completed in 277.8 s. A was retained; B was
refitted against the successful diagnostic teachers' expected afferents, then
C joint PF, RF/HAM input, and motor weights were fitted at both cadences,
followed by D at the normal cadence. B periods were 1.233 s right and 1.232 s
left; the small timing kernel matched the native neurons within 1.16e-13.
The final expected-input motor replay RMSE was 0.02198 activation units.

That aggregate error concealed poor burst reproduction. Fresh body runs failed
stages 1 and 2: stage 1 had 0/4 counted knee cycles and a small sampled ankle
crossing; stage 2 had 0/0 counted knee cycles. Adding axis feedback did not
rescue stage 2. All runs were finite. No later stage was advanced.

Free expected-input PF traces already show a winner staying active while its
antagonist remains subthreshold. The largest afferented motor targets are
predicted mostly as near-constant means, with R-squared near zero. Therefore
the failure is present before coupling to the body. Whole-network RMSE is
insufficient as a Hunt C acceptance measure; the native validation now also
records per-muscle R-squared, active-muscle error, and predicted ranges.

Two isolated PF-only fits (`pf_subthreshold_1_20261030` and
`pf_subthreshold_2_20261030`) test inactive target voltages of -1 and -2 model
units, respectively. A/B and the measured teacher data are held fixed. Native,
unclamped hip/knee PF waveforms are checked before spending another fit on
motors/reflexes. These are preliminary fitting checks, not curriculum passes.
Raw traces and failed configurations remain under ignored runs folders.

## Session 1: native PF checks and live-editor verification, 2026-10-30

The initial file-URL restriction was followed by successful access to the HTTP
page Ben opened at http://127.0.0.1:8765/connectome_block_editor.html. The Hunt
module template and both RF/HAM detail templates were loaded through the
editor menu. Entering the right RF/HAM SUB exposed its four cells and boundary
ports. Right/left detail views contain 12 nodes and 10/11 enabled edges. The
two inhibitory projections were verified as red paths. Full numeric gains
remain in the data; labels display three significant digits. Focus canvas
hides the editing panels for a readable diagram. Live screenshots are saved
as Hunt_fourcell_session1_{right,left}_live.jpg in the export folder.

These drawings use the fourcell_dynamic_20261029_v1 snapshot, which failed
ground stepping. They are not a newly successful controller. Browser actions
stopped when Ben reported concurrent editing. Leave his current page and tabs
alone. The earlier statement that live verification was incomplete is superseded
by this bounded check; the entire flat graph was not loaded in the live test.

The -1/-2 inactive-voltage PF fits failed their free-waveform checks. Reducing
the inherited joint-PF inhibitory conductance to 0.15 of its prior value,
alone or combined with a -1 inactive target, also failed. These candidates
retained their four neurons and were not promoted to physical stages.

A free four-cell fitting kernel now handles each hip/knee PF pair, its two
inhibitory relays, and the recorded phase-filter inputs. It reproduces the
saved full-network PF trace within 6.67e-16 model-voltage units. The initial
joint fit took 62.56 s for four pairs, with at most 100 optimizer evaluations
per pair. This is a subset fitting time, not an end-to-end simulation speedup.
Fitting rectified native PF activity avoids forcing an inhibited cell to zero
voltage and then learning excitation to cancel that inhibition.

The two-cadence native fit improved slow-cadence target reproduction (7 of 8
hip/knee cell waveforms passed preliminary criteria) but failed all 8 normal-
cadence checks. A normal-cadence-only comparison improved the hip outputs but
still failed knee checks. Physical curriculum gates were unchanged. No new
motor/reflex fit or later physical stage was accepted from these results.

The full-network timing alignment was checked explicitly: afferented RG
periods are 1.234/1.232 s, phase errors 0.0101/0.00406 cycles, and stance MSE
0.0214/0.00687. Deafferented periods remain 2.466 s bilaterally. The upstream
timing check therefore does not explain the poor normal-cadence knee PF fit.

The next isolated fit separates joint motor primitives from biarticular torque:
hip targets use gluteus-maximus versus iliacus/psoas commands, and knee targets
use vasti versus biceps-femoris-short-head commands. RF and the three long
hamstrings retain their separate paired PF targets; sartorius and gracilis
remain separate. These are normalized simulation-derived command targets,
not EMG or a claimed literature replication. Raw failed fits remain in runs.

## Session 1 — mono-command checkpoint and extended run guide, 2026-10-30

The native joint-PF mono-command fit and its normal-cadence motor/reflex fitting
finished in `runs/native_joint_pf_mono_20261030_v1`. Its afferented stage-2 air
retry passes the current physical gate: 7/7 knee cycles, approximately 1.229 s
period, phase error 0.0537 cycles, and minimum ankle width 0.0323 m. The default
deafferented stage 1 still crosses. The gain-8 abductor release comparison passes
stage 1 (3/3 cycles, 2.466 s period, minimum width 0.0295 m) but fails stage 2's
frontal gate (minimum width -0.0663 m). These are different passing configurations;
the combined early-stage gate remains unmet. No later curriculum stage was
advanced. Next work is a bounded frontal sensory/torque diagnosis with fitted
weights frozen, followed by fresh-process repeats of both early stages.

Ben requested a detailed runnable guide and bounded jobs for other models.
Added `Code/MuJoCo_SNS/how_to_run_code_extended.md` and nine standalone jobs in
`Code/MuJoCo_SNS/ZCODE_SIMULATION_PROMPTS.md`. They name easteregg2 and exact
paths, distinguish the actual controller families and coordinate frames, separate
Hunt A/B/C/D fitting from corrected physical stages, and retain failed-stage
gates. The original `HowToRunCode.md`, CHATGPT_REPORT.md, dissertation, production
builders/parameters, saved winners, and the user's live editor were not changed.

New opt-in tools are `run_spinal_extended.py`, `run_w2l_bpa_extended.py`, and
`plot_recording_extended.py` in Code/MuJoCo_SNS, plus
`run_sns_model_extended.m` and `prepare_spiking_library_extended.m` in
Code/Matlab/SNS_Simscape. They select existing models, write separate results,
and provide named command/activation/pressure/neural/foot/harness plotting.
The Simulink spiking knee initially failed because the shared library lacked
SpikingNeuron/HybridSpikingSynapse. The preparation helper builds those blocks
in a private copy; the shared library remains unchanged.

Verification on easteregg2:

- s3k, w2lvar, and syn6: native versus preallocated fresh-process 15.3 s schedules
  have zero maximum difference in qfull, ctrl, and neuro. This is trajectory
  parity, not a walking pass or a controlled performance benchmark.
- W2L BPA graded: complete 12 s native/preallocated recordings match exactly in
  t, q, foot forces, height/progression, pressure, and RG-extensor trace. A fresh
  native repetition also matches. The spiking BPA version completes 12 s with
  finite state. These old drop-protocol metrics do not establish corrected
  alternating ground stepping; the graded trial has no right-foot load.
- Stock human hybrid spiking mirror: native 0.04 s standing smoke run completes.
  The graded fast kernels are explicitly rejected for spiking networks.
- MATLAB R2025a Update 1: graded knee and private-library spiking knee each
  complete five-second runs. Small fixed-point gain-quantization warnings remain
  in the spiking log. The Deng CPG also completes its earlier 0.02 s smoke run.
  No compiled MuJoCo bridge MEX was found in this checkout.
- Named joint, neural, command, native activation, pressure, contact, and harness
  selections were exercised. Command-only data are rejected as activation state.
  Generated plots were visually checked for labels and line/marker distinctions.

`render_air_preview.py` now uses MuJoCo 2.3.7 renderer destruction rather than
the unavailable close method. It exposes time range, frame density, and playback
speed, and derives GIF delays from the saved time axis. Cumulative centisecond
quantization prevents frame-by-frame rounding drift. A private 2.016 s source
clip encodes as 2.020 s at 1x playback; poses are replayed without correction.

Raw verification, plots, MAT/NPZ files, and the private library are ignored under
`tmp/extended_guide_checks`. `final_guide_verification.json` records parity,
plot checks, the activation-label rejection, and GIF timing. The guide also
flags the legacy BPA hip-RoM sign error and the evaluator's real JSON keys:
`diagnose_swing` enables mechanics logging, and the configuration filename stem
determines its output prefix. No claim that an unused JSON key redirects output.

The staged Simulink knee MAT-to-GIF/still recipe was also exercised. It exposed
an indexed-palette bug in `demos/sns_animate_demo_v2.m`: every frame was written
using the final frame's palette. The renderer now retains each frame's own
palette for knee and cup animations. No model dynamics or shared figure outputs
were changed. The legacy fixed frame delay and angle-derived decorative muscle
activity remain explicitly identified in the guide.

## Session 1 — shared motor fit at both cadences, 2026-10-30

Ben assigned worker jobs A/B to ZCode. Session 1 leaves their configuration
sweep and RF/HAM force comparison to those workers. Prompt B's mechanics path
was corrected to `spinal/functional_gait2392/mechanics.py`.

The current mono-command motor fit used the afferented cadence only, whereas
stage 1 evaluates the approximately twice-as-long free period. New isolated
experiment `runs/multicadence_mn_20261030_s1_v1` freezes timing, native four-cell
PF, and RF/HAM integration weights, fits one MN conductance set to both expected
cadences, and then fits residual feedback at the afferented cadence. No shared
winner or topology changes. The copied frozen weights are hash-checked.
Physical stages 1 and 2 will be rerun in fresh processes with one configuration;
no later physical stage is authorized by a regression score.

Stage-3 clarification: older `balance_20261029_v3` assisted-standing cases
passed, including the disabled-added-pathway control. They do not establish
a benefit from the added task-balance pathway. The current four-cell candidate
has not been promoted through stages 1/2 to stage 3. Harness-supported standing
and alternating ground stepping remain different acceptance targets.

The shared two-cadence MN fit finished in 181.9 s; its copied A/B/PF/BI hashes
were unchanged. Fresh accelerated body checks fail both stages: stage 1 has
1/4 counted knee cycles and stage 2 has 2/5, although both frontal screens pass.
The fit therefore trades motion quality for frontal separation and is retained
as a failure, not promoted.

Ben questioned unequal PF-IN return conductances. The current left-knee pair,
for example, has F-IN→E = 1.8537 and E-IN→F = 0.8649. These were independently
learned, not literature-prescribed conductances. In `runs/tied_pf_20261030_s1_v1`,
setting both opposing return conductances to their arithmetic or geometric mean
without retraining fails stages 1/2; the latter preserves each pair's product.
That ablation does not show the asymmetry is necessary, because the existing
motor mapping was fitted to the original PF dynamics. A new constrained fit
ties the two return gains during native PF optimization (49 rather than 50
free parameters per pair), then refits RF/HAM input, MN, and residual feedback.
The experiment is `runs/tied_native_pf_20261030_s1_v1`; four-cell topology,
intrinsic rhythm, and expected-afferent timing remain fixed.

Another isolated test, `runs/tonic_abduction_20261030_s1_v1`, adds small
conductances from the existing gait-gated tonic PF cell to hip-abductor MNs.
It separates constant support from scaling all phasic inputs to those muscles,
and is distinct from ZCode A's scale-8 positional-feedback sweep. Added
conductance is not described as a guaranteed activation floor. No body restraint,
clock, or imposed motor target is added.

The tonic-PF increment 0.02 passes both stages with one saved configuration:
stage 1 has 3/3 knee cycles, period 2.467 s, phase error 0.0773, and minimum
ankle width 0.0991 m; stage 2 has 7/7 cycles, period 1.214 s, phase error 0.0653,
and minimum width 0.1086 m. A separate fresh-process repetition passes both.
The increment is an added synaptic conductance, not a measured 2% activation
floor. Scaling the full phasic abductor mapping is not used. The 0.04 and 0.08
tonic increments fail at least one air gate; stronger support is not monotonic.

The six one-at-a-time +/-1% checks of the tonic increment, commissural gain,
and RF/HAM pair gain pass both stages. Two independent simultaneous jitters of
PF/BI/MN/reflex edge weights are being checked; intrinsic/expected-afferent fits
remain fixed, so this is not a basin test of every model parameter.

With these controller weights, stage 3 now passes under the declared assisted
standing environment: out-of-plane restriction 1.0, -0.02 m initial height offset,
hip/knee/ankle initialization -7/-10/0 degrees. Feet carry about 0.550 BW and
the vertical rig spring about 0.446 BW with the added task pathway enabled.
Horizontal harness RMS is 0.0770 BW, sway 0.0448 m, and maximum tilt 11.22 deg.
The identical-topology task-off control also passes with nearly identical
measurements. This establishes an assisted-standing screen, not a benefit of
the added pathway or independent balance. Restriction 1 adds 800 N m/rad and
40 N m s/rad on selected out-of-plane joints; it is substantial assistance.

Corrected stage-4 probes now begin only after those repeated air passes and
the standing pass. `ground_initial` retains the same weights and compares the
standing spring support with prescribed 70% BW support (vertical spring/damper
removed). Ground is enabled and acceptance uses actual alternating contact and
unloading gates. No stage-4 or later success is claimed before its result.

The constrained tied-inhibition PF fit completed in 50.9 s but failed its native
waveform screen; after separate RF/HAM/MN/reflex fitting, fresh physical checks
also fail (stage 1 0/2 cycles; stage 2 0/0). It is retained as a failed constrained
comparison. This does not establish unequal gains as biologically required.

The full basin check is 7/8: all six selected scalar probes and one independent
edge-jitter seed pass both air stages; seed 301031 fails stage 2 with 7/0 knee
cycles (stage 1 still passes). Nominal fresh repeats remain passing. This is a
broader synaptic sensitivity limitation, not a completed all-parameter basin.

The first stage-4 spring and constant-70%-support trials fail. Spring support
leaves only 0.176 BW on the feet and the left leg continuously loaded. Constant
70% support gives about 0.291 BW total foot load and large knee excursions,
but both feet remain continuously loaded. The left RG has extensor dominance
throughout the scored interval. The observed knee motion includes a 0.305 m
vertical CoM oscillation; it is not alternating ground stepping.

Onset initialization, weaker timing input, stronger commissural inhibition,
and a timing-input-off diagnostic all still fail the contact gate. With timing
input off, RG activity becomes largely antiphase, yet the feet remain loaded:
releasing the neural arrest alone does not solve the mechanical transfer of load.
The fitted left hip-extensor length input has gain 13.018 (right 0.00268) and
saturates under ground-posture signals. That is a marked distribution shift
from the expected-air afferent data; it is not evidence that the fitted gain
should be copied as a biological constant.

Contact wiring clarification: the passing air candidate has cross_contact_gain=0,
so no optional HEEL→contralateral RG-F excitatory edges are added. Recorded heel
cell activity is exactly zero in air and positive with real ground contact.
`contact_pathway` tests gains 1, 4, and 12 plus a spring-supported comparison,
and exports the compiled edges for verification. No pathway benefit is claimed
before its measured unloading/alternation results.

Zero-pathway clarification: zero conductance contributes no synaptic current.
The fixed-node sparse engine retains the dense pairwise summation tree. Adding
or removing relay cells changes that tree and may change sensitive trajectories.
Optional commissural, RF/HAM-inhibition, and crossed-contact motifs are added
conditionally; the learner also intentionally retains its fixed candidate slots.
The standing task-off comparator retained eight task cells with zero external
task drive. It was a same-topology control, not physical pathway removal.
A separate stage-3 comparison now sets balance=None, omitting those cells.

Physical task-module removal also passes assisted standing: 1,027 neurons
versus 1,035 in the same-topology task-off control. qfull, ctrl, and CoM traces
match exactly for this comparison. The harness and inherited posture/reflex
controller remain present. No independent-balance claim follows.

The direct crossed-contact tests at gains 1/4/12 and 40/100 do not pass ground
stepping. Positive-Ib-group output reinforcement by 0.3 does produce repeated
measured unloading/loading on both feet (5/9 accepted force onsets, duties
0.544/0.662) under constant 70% support, but contact coverage/phase and the
movement phase still fail. The same modified weight set passes both air stages.
These IBEXC group cells receive RG-E excitation as well as tendon-force input;
this change is not a pure force-sensor ablation or a direct Ib→opposite-F claim.

Extensor Ib commonly reinforces ipsilateral extension and delays flexion during
loaded stance; it is not interchangeable with the proposed crossed cutaneous
swing cue. Primary experimental reference checked: Pearson et al., 1992,
https://pubmed.ncbi.nlm.nih.gov/1426112/ . The code's loading and phase-specific
recruitment remain an engineering implementation.

An opt-in contact-event encoder now replaces the tonic whole-foot proxy with
a transient after measured 0.05/0.025 BW hysteresis and >=80 ms preceding
unloading. Its current is delivered to the existing heel neuron and the
explicit crossed RG-F projection; no stride clock or imposed joint trajectory.
`cross_contact_event_amplitude` defaults to zero. Three event gains and a lower
Ib-support variant produce additional bilateral contacts but still fail the
unmodified ground alternation criteria. Event currents are saved separately.
`harness_vertical_offset_m` is an explicit signed world-z spring-reference
offset (positive upward), default zero. Existing `harness_lowering_m` retains
its historical positive-world-z behavior; its name is misleading and should
not be interpreted as downward motion. Default-disabled new options reproduce
both original air trajectories exactly in qfull, ctrl, and neuro.

The constant-force harness permits large common vertical oscillation, while the
original 500 kN/m vertical spring with increased stance recruitment caused
early instability-guard termination. Correction to the initial classification:
the saved prefixes were finite (nan=False); these were incomplete trials,
not demonstrated nonfinite states. Zero-filled recording tails are invalid.
Neither result counts as stepping. A subsequent comparison combined 65–70%
prescribed support with 2–4% of that vertical stiffness and damping scaled by
its square root. All four compliant-height cases also ended early and fail
the unchanged completeness gate. Inertia/acceleration warnings require
separate diagnosis; they are not evidence of a postural fall.

ZCode A/B results were read from their job-local reports. A's three scale-8
axis-feedback cases all fail (left knee has zero counted cycles); no stage-1
promotion was performed by that worker. B confirms command/activation/torque
plumbing and much weaker, often braking-dominant motor activity in the newer
air model. Those are air diagnostics, not evidence of ground stance strength.

## Session 1 — symmetric PF reference and selective stance recruitment, 2026-10-30

`runs/direct_pf_reference_20261030_s1_v1` tests graded, directly RG-driven
four-cell hip/knee PF circuits: RG→PF gain 2.4, symmetric PF-IN return
conductances 2.749, and no phase-bank input to those hip/knee PF cells.
Other functional PF circuits remain unchanged. Hunt A/B are frozen, and
RF/HAM integration, MN, and residual-feedback weights are refitted.
This is a driven graded reference, not an intrinsic-NaP PF replication.
Stage 1 has 3/3 cycles and phase error 0.0942, but ankle-origin width reaches
−0.0166 m. Stage 2 has 6/8 cycles and phase error 0.2215. Both fail; neither
is promoted. Symmetric reference conductances alone do not solve transfer.

`tonic_abduction_20261030_s1_v1/0p02/mono_stance_support` restricts the added
positive IBEXC-group output to gluteus maximus, vasti, and soleus (14 existing
edges). It avoids additionally recruiting rectus femoris and gastrocnemius,
whose actions include flexion at another joint. The group neurons still
combine RG-E and force inputs; this is an engineering recruitment comparison,
not a pure Ib sensor experiment or proof of a physiological connection rule.

The 0.3 increment passes stages 1/2 with unchanged weights: stage 2 has
approximately 0.989 s period, phase error 0.1202, and minimum ankle-origin
width 0.1456 m. It completes the 8 s standing trial without a fall, with
mean once-counted foot load 0.888 BW and vertical spring support 0.111 BW.
It nevertheless fails stage 3: maximum tilt 22.1905° and sway 0.10417 m
exceed the unchanged 15°/0.06 m limits. Its ground trials are withheld.
The 0.6 increment fails stage 2 (phase error 0.1714 and crossing). No ground
or stage-5 promotion follows from greater loading alone.

The original 0.02 tonic-abductor candidate remains the accepted stages-1/2
configuration and assisted-standing checkpoint. After adding the default-off
direct-PF option, both stages were rerun in fresh processes. Full qpos,
excitation, and recorded neural-output arrays match the earlier repeat exactly
(maximum absolute differences 0.0). `CHECKPOINT.json` in that candidate folder
records weight hashes, source snapshots, accepted assessments, and limitations.
The source and raw trial outputs remain in ignored `runs/`; no production
winner, physics timestep, or acceptance threshold was changed.

An 8 s, 1× playback air-stepping preview and still are saved under its
`air_pair/preview/`, with a full-support caption. The anterior camera now views
from native world +x; an initial azimuth-0 render incorrectly faced posteriorly
and was replaced before delivery. Recorded poses were not corrected.

Ground transfer still has at least two separate problems: the fitted left
hip-extensor length timing input saturates in loaded poses, and the air-trained
motor mapping provides inadequate phase-specific stance/swing force. Timing
removal alone did not unload the feet; stronger stance recruitment alone did
not satisfy balance/alternation. Further fitting needs loaded-task targets
and explicit stance/swing recruitment checks, rather than accepting air motion
or increasing crossed-contact gain as sufficient evidence of walking.

## Session 1 — Hunt citation and regional foot sensing, 2026-10-30

Ben supplied the presentation citation: Hunt, A. J. (2019, February 15).
*Real-time control of dynamic walking in a dog-like robot*. Presented at the
2nd Intel Neuromorphic Research Community (INRC2) Meeting. This is a professional
meeting presentation, not a journal paper. Exact citation and bounded results
for the separate dissertation/defense session are saved in
`DISSERTATION_DEFENSE_UPDATE_20261030_SESSION1.md`; manuscript/deck sources
were not modified.

`regional_contacts.py` partitions actual ground contacts into hindfoot and
forefoot regions, using a configurable midpoint between calcaneus/toe origins.
Each contact force enters one region only. These are virtual regional sensors,
not calibrated anatomical pads. The experimental evaluator replaces HEEL/TOE
inputs, and subsequently also corrects LOAD to use their once-counted sum.
Default regional_contacts=False preserves the legacy production baseline.
Heel-event mode now receives hindfoot loading only; the forefoot signal is
separately logged. No new toe→RG/PF pathway has been promoted in stage 4.

Five `regional_ground_v1` diagnostics fail. Three `regional_ground_v2` tests
also fail after correcting LOAD; v1's evaluator snapshot preserves its earlier
LOAD behavior. Regional forces sum to independently measured foot forces to
within 9.1e-13 N. Without stance reinforcement, both feet remain continuously
loaded and the left RG remains extension-dominant (about 96% of samples after
5 s), consistent with the previously identified length-feedback problem.

With IBEXC output increment 0.3, measured contact onsets are 5/9, but phase
error 0.1974 and coverage 0.50/0.625 fail. All scored regional load in that
case is forefoot; hindfoot force is zero. Thus the heel-only event cue cannot
explain its repeated contacts. Removing timing feedback in this combination
ends before the scoring window, so it supplies no accepted gait evidence.
The provisional hypothesis is excessive plantarflexor recruitment/foot posture;
the proximal_stance_v1 ablation restores baseline ankle IBEXC outputs while
retaining hip/knee increments. These are diagnostic stage-4 trials, not
promoted replacements for the shared early-stage winner.

ZCode's clarification is accepted: the 0.02 tonic candidate already clears the
shared early-stage prerequisite for Jobs C/I. Reuse existing sensitivity tests;
do not repeat the failed gain-8 axis sweep. Stage-4 execution is authorized by
Ben's current instruction. Loaded-task Hunt fitting, rather than another air
pass, is the remaining training objective.

The proximal-only stance comparison completes but still fails stage 4.
Restoring baseline ankle IBEXC outputs brings hindfoot loading back. With
heel events, contact phase error improves to 0.1054, but coverage remains
0.333/0.600 (4/6 events), so this is not accepted alternating stepping.
The no-cross comparator has phase error 0.2779 and also fails.

Ben asked explicitly about dragging and tiptoe behavior. Recorded-state mesh
replay against the actual ground-plane normal confirms persistent left-foot
skimming in the heel-event trial: during force-defined unloading (<0.025 BW),
median minimum foot clearance is 0.887 mm; maximum left clearance in the
scored window is 8.615 mm. None of its sampled unloaded frames exceeds 10 mm.
Right clearance reaches 77.293 mm, but exceeds 10 mm in only 57.4% of unloaded
samples. Sampling is every 20 ms, and 10 mm is a diagnostic reference, not a
claimed literature-derived acceptance threshold. This is not a clean swing.
The left RG remains extension-dominant in about 79.5% of samples after 5 s,
versus 49.2% on the right. The earlier ankle boost produced forefoot-only
loading; removing it fixes that distribution but does not fix left clearance.
No backward-lean mechanism is inferred solely from pelvis tilt.

Raw clearance evidence: proximal_stance_v1/CLEARANCE_DIAGNOSIS.json, with
the replay driver preserved alongside it. Ground steps must demonstrate foot
clearance as well as force unloading before a walking claim. The next Hunt
fit should address loaded-posture timing (especially the left length-to-E
gain), swing knee flexion/dorsiflexion, and stance/swing force targets; simply
increasing heel-reset strength does not resolve this measured asymmetry.

### Harness-height and strength clarification

The latest proximal-stance stage-4 trials are force-controlled: vertical
spring/damping disabled, constant upward force 0.70 BW, and initial pelvis
height offset −0.02 m. There is no prescribed ongoing vertical harness height.
The heel-event run measures foot loads 0.09558/0.20266 BW and essentially zero
vertical spring support. Low left clearance is not evidence of suspension
above the floor; the left foot carries more load than the right.

A direct audit of the repaired model finds that bilateral maximum muscle
forces are not exactly equal. For example, soleus Fmax is 3395.50 N right
versus 3186.36 N left; tibialis anterior is 840.09 versus 852.155 N.
`proximal_stance_v1/SYMMETRY_AUDIT.json` records all 46 matched pairs.
Their provenance (including possible converter normalization) needs tracing
before any capacity correction; no symmetric-capacity change was silently made.
Neural gains are separately asymmetric: the fitted hip-extensor length→RG-E
gain is 0.0026777 right versus 13.0181 left. Thus neither mechanical capacity
nor learned neural drive should be described as enforced bilateral symmetry.
The earlier focus on neural asymmetry did not constitute a muscle-capacity
audit. The asymmetries must be isolated before attributing the whole failure
to a single pathway.

## Session 1 — strength provenance, shared parameters, and elastic physics

Ben requested correction of the bilateral feedback fit and a subagent for
compliant tendons/damping. The strength audit is in STRENGTH_ALIGNMENT_AUDIT.md:
source OpenSim soleus Fmax is 3549 N on both sides. Independent MyoConverter
Step-3 multipliers 0.9567481/0.89782008 create 3395.5/3186.36 N parameters.
Across 84 cached fits the formula reproduces XML within 0.00466 N. These are
surrogate curve parameters; equalizing Fmax alone does not align force laws.
Full-activation BFSH neutral forces differ (520.62/605.54 N) despite matched
sagittal routing. The requested physics subagent is auditing full force curves.

The independent left timing fit was inappropriate for the intended bilateral
controller. bilateral.py pools mirrored regression rows/columns into shared
parameters. HUNT_BILATERAL_TIED=1 skips independent solves;
require_bilateral_symmetry rejects unequal candidate timing/saved conductances.
Synthetic pooled-regression and mismatch-rejection checks pass. Historical
weights remain available but do not meet the new symmetry requirement.

symmetric_loadedstress_20261030_s1_v1 fits six common timing gains to both legs
and high-length-offset probes, with intrinsic A frozen. All four neural tests
pass. Shared extensor-length gain is 0.002815157, replacing 13.0181 on the left.
Native/reference errors are below 1.2e-13. These are stress probes, not loaded
walking ground truth. PF weights are explicitly projected right-to-left;
BI/MN/reflex weights are jointly fitted. The exact tonic variant passes stage 2
(7/7 cycles, phase error 0.06446) but fails stage 1 (3/4 cycles, coverage 0.5).
Common motor gains 1.25/1.5/2 do not repair its extra slow-stage knee pulse.
Two separately labelled stage-4 diagnostics also fail; no promotion is claimed.

symmetric_bothcadence_20261030_s1_v2 pools both cadences in the tied motor fit.
Its exact tonic-abductor variant passes stages 1 and 2 with one saved set.
Assisted standing and ground progression are being evaluated. An initial
two-cadence attempt stopped on duplicate solver-capture bookkeeping; that
duplicate was removed before the completed fit. No failure gate was relaxed.

The physics subagent implemented spinal/physics_fidelity/series_elastic.py plus
isolated tests and benchmarks. Matched OpenSim static equilibrium agrees within
2.16e-6 N. A warmed 92-unit compliant bank takes 2.47 s for 16 simulated s,
versus 0.0776 s for its rigid comparator. These are bank-only timings. Dynamic
FV remains approximate, and beta=.05 is an experimental regularizer. This is
not yet a validated production/Hyfydy replacement; see physics_fidelity/README.md.

The legacy nonlinear-contact option misindexes solimp width/midpoint and makes
impedance endpoints equal; it does not implement the claimed depth-dependent
profile. The old ~36 kg m^2 follower-inertia docstring is unsupported: measured
right-knee contributions are 0.03256 neutral and 0.02570 kg m^2 at swing,
6.19%/5.35% of reduced inertia. ID mode armature is .005, not zero. The current
feasible teacher uses live forward physics, so a direct ID-to-C mismatch is not
established. No production physics was silently changed.
