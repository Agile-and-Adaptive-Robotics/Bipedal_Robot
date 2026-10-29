# Gait2392: Hunt four-step training and five-stage validation

Machine: easteregg2. Repository: `D:\GitHub\Bipedal_Robot`.
Python: `D:\Anaconda\envs\myo\python.exe`.

This is the working experiment record, separate from `CHATGPT_REPORT.md`.
Raw arrays, trial configurations, logs, and databases belong in the ignored
`runs/` directory. Existing production winners and the completed exploratory
campaign remain separate.

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
