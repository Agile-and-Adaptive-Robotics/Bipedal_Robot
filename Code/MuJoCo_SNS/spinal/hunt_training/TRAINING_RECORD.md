# Gait2392: Hunt four-step training and five-stage validation

Machine: easteregg2. Repository: `D:\GitHub\Bipedal_Robot`.
Python: `D:\Anaconda\envs\myo\python.exe`.

This is the working experiment record, separate from `CHATGPT_REPORT.md`.
Raw arrays, trial configurations, logs, and databases belong in the ignored
`runs/` directory. Existing production winners and the completed exploratory
campaign remain separate.

Current Session 1 checkpoint: three RF/HAM four-cell variants trained and
78 fresh-process physical evaluations recorded, including replays and graph
snapshots. No fit passes the complete corrected curriculum.
See `FOURCELL_HUNT_RESULTS.md` for the current comparison and measured gates.
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
