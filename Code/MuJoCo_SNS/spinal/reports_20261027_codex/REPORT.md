# Simulation recovery results

> **UPDATE 2026-10-27 (later session): the Di Russo native-replay blocker
> below is CLEARED.** The author's SCONE fork was built on easteregg2 and
> `spinal.par` replays native walking (16 steps, 1.166 m/s, objective
> 17.2411, bit-exact fresh-process replay). See
> `native_baseline/NATIVE_BASELINE_20261027.md` for the run of record,
> the build manifest, and the controller audit.

The W2L candidate survives fresh-process verification with its effective gains
delivered to the compiled network. The AH001 CAD robot now compiles with its
four-bar knees, lumbar hinge, and 75 kg total mass. Di Russo walking remains
unreplicated. These results do not establish harness-free walking.

## Task status

| Task | Result |
|---|---|
| 1 W2L conductance tuning | Recovered and verified for supported stepping; physiological walking remains unproved |
| 2 Di Russo replication | Original source and parameters recovered; native controller baseline and walking replication incomplete |
| 4 Speed improvements | Exact W2L backend verified over a full trajectory; parallel fresh processes used; corrected Di Russo adaptation compiled |
| 5 Stimulus panels and actual robot | W2L panels regenerated with effective keys; CAD robot mechanics verified; actual robot neural locomotion incomplete |
| 3 Existing human and robot families | Not launched because tasks 2 and 5 have not met the requested completion condition |

## W2L recovery

The previous fresh-process verifier did not register the per-target gains before
the evaluator filtered its input. Registering them at the evaluator boundary
restored the saved candidates. Unknown settings now raise an error rather than
silently disappearing.

The compiled-matrix audit identifies 56 effective conductance dimensions and
three runtime sensor dimensions. Two saved settings, `ia_ank_rg` and
`ia_ank_pf`, have no corresponding edges in this template. The canonical
configuration omits these inert settings and reproduces the original trajectory
exactly. This is a search over conductance classes, not 56 independently tuned
individual anatomical synapses.

Trial 458, over the 4.5–30 s analysis interval, produced 28 left and 21 right
contact events, median phase 0.496063, forward displacement 2.246231 m, and fitted
forward speed 0.091786 m/s. Four independent replays match the saved trajectory
exactly. All 12 perturbations of every saved numeric setting by ±1%, excluding
the fixed 30 s horizon, pass the supported-stepping gate. The more balanced
candidate 360 also passes its separate 12-condition basin check.

The frozen trial 458 configuration passes all tested weight-support settings
down to `wsup=0`. That setting leaves the vertical PD harness active: its
average force is about 16.76% of body weight. At the baseline operating point,
average harness force is about 18.59% of body weight. Contact duties near 0.26
and unequal event counts further limit the interpretation. These are supported
stepping results, not established physiological or harness-free walking.

Evidence: `plumbing_audit.json`, `recovery_results.json`,
`validation_results.json`, `full_basin_results.json`, `verification.json`, and
`w2l_verified_config.json`. Run `run_verified_w2l.py` for a timestamped fresh
evaluation. Raw NPZ trajectories and logs remain local.

## Speed and numerical equivalence

The previous sparse backend fails full-trajectory parity: its maximum joint
angle deviation reaches 45.56° despite similar aggregate metrics. One-step
agreement does not establish closed-loop equivalence near these bifurcations.

`w2l_cpg/exact_fast_eval.py` preserves the reference reduction order and gated
channel equations. Its 30 s joint, pelvis, center-of-mass, force, neural, and
harness traces match the reference exactly. Set `AARL_W2L_BACKEND=exact_fast` to
select it; recovery scripts do this explicitly. The measured isolated neural
step improves from 145.76 to 102.29 μs, or 1.425×. This is not a whole-trial
speedup claim. Physics integration, initialization, and the support rig remain
in the authoritative runner. Eight parallel independent processes were used.

The corrected Di Russo-inspired adaptation has a compiled neural step matching
its Python equations within 3.34×10⁻¹⁶ over 2000 steps. Its isolated benchmark
is about 450× faster. This validates the implementation of that adaptation,
not equivalence to the original author's controller.

## Di Russo source and failed adaptation

The [original author repository](https://github.com/DiRussoAndrea/Spinal_controller)
and [custom SCONE fork](https://gitlab.com/simgait/SCONE) supersede the previous
report's attribution to `neumovelab/spinal-locomotion-control`. The latter is a
comparative reimplementation, not the original source. The fork checkout used
here is revision `58ea68762791e133b0b18f2f0db83e3eb2a18fe3`.

The original scenario uses `gait9dof18musc_Thelen_20170320.osim`: a reduced
planar Delp/Thelen model with nine active coordinates and 18 Thelen muscles.
Its tenth named coordinate, lumbar extension, is locked. It is related to the
Delp gait-model lineage but is not the full gait2392 or gait2354 model. The
source model loads and initializes in this machine's OpenSim 4.4.1.

The optimized `spinal.par` contains 256 parameters: one frequency, 16 state
offsets, two phase offsets, nine balance parameters, 45 pattern-formation
weights, 174 individual reflex connections, and nine length offsets. Original
source details include delayed state histories, muscle-specific delays, fiber
length normalization, and Renshaw inhibition of antagonist Ia interneurons.
The current adaptation does not implement all of these details.

Equation and sensor corrections in the adaptation address phase units,
synchronous phase updates, raised cosines, output saturation, motor-neuron
time constants, biarticular reversal, balance placement, sensor ordering,
contact-force normal, and initial-state/model repairs. A 960-trial exploratory
survival/speed search nevertheless selected a somersault: pelvis rotation
reached about 223°. The 6.98 s candidate at 90% support is **rejected** as
walking. Its plot is `dirusso_attempt.png`. Further optimization with this
surrogate was stopped; the script requires `--debug-surrogate` to reproduce it.

Native replay is blocked by the absence of the author's custom SegmentalController
runtime on this machine. The fork documents Ubuntu 18.04 as its supported build
route and says Windows builds require substantial work. The locally installed
OpenSim library alone does not supply that controller. Establish the native
original-model baseline before porting its full objective and controller to
gait2392 or spending further optimization effort. See `ZCODE_HANDOFF.md`.

## Perturbations and CAD robot

`w2l_perturbation_panels.png` and `feedback_comparison.png` use the verified
configuration at `wsup=0.25`. Mild and strong drive reductions, ±40 N pushes,
and Ia, Ib, Renshaw, crossed-heel, and heel-cascade deletions pass the existing
supported-stepping gate. Commissural deletion fails bilateral stepping.
Removing all feedback fails the contact-phase gate, although joint cycling and
forward motion persist. It does not demonstrate that feedback is necessary for
all oscillation. These are W2L adaptation experiments inspired by Deng's
protocol, not an original-model replication. Detailed outcomes are in
`perturbation_results.json`.

The intended Simulink import contains spherical hips and reduced revolute knees;
the old report's statement that it has no joints is incorrect. This machine's
Simscape Multibody license is unavailable, and its data-file path points to the
laptop. Earlier model-build logs also show rigidly connected knee frames.

The CAD fallback now provides nine rigid groups, current solid volumes and
inertias, matched imported joint anchors, both complete four-bar knees, and the
approved adjustable torso. The source CAD was opened read-only and not saved.
The complete model compiles at exactly 75 kg. See
`../../robot_torso/README.md`, `ah001_build_audit.json`, and
`ah001_verification.json` in that folder for mass assumptions and mechanical
checks. Six provisional passive trunk routes now have insertion points on the
spine. Their actuator capacity and neural coupling remain to be established.

The user subsequently selected adjustable routes based on the existing
robot-body OpenSim model. `../../robot_torso/ah001_with_opensim_routes.xml` now
contains 92 routes and 354 editable attachment/guide points. Source neutral
lengths reconstruct within 1.12 × 10⁻¹⁶ m. Twenty routes lose joint leverage
without foot articulation. Bilateral `add_mag2`, bilateral `grac`, and left
`tfl` reverse neutral sagittal action after transfer and need attachment review.
The 65-pose motion audit preserves the knee loops; tendon derivatives and
finite differences agree within 4.67 × 10⁻¹⁰ m. Moving points and conditional
membership remain frozen at source neutral, and all force interfaces have zero
gain. This is an adjustable routing scaffold, not a neural locomotion result.
The transferred erector routes require approximately 1354.84 N each for neutral
holding, compared with 1394.61 N for the earlier six-route template.

## Review and commit scope

Review this report directory, `robot_torso/`, `w2l_cpg/exact_fast_eval.py`, the
backend selection in `w2l_cpg/build_w2l_net.py`, gain registration in
`reports_20261016/w2l_wean/_w2l_wean_eval.py`, shared-XML handling in
`w2l_mujoco/test_w2l_ground.py`, and the Di Russo controller and plant changes.
The checkout already contains unrelated changes. No commits, pushes, or
dissertation edits were made in this recovery work. Generated meshes, raw
trajectories, helper binaries, and external source checkouts are ignored.
