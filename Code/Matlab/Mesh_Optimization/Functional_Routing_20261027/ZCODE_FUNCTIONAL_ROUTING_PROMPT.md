# Functional BPA routing and consolidation

Work on **easteregg2**, repository **D:\GitHub\Bipedal_Robot**. Check hostname
first. Python/OpenSim: **D:\Anaconda\envs\myo\python.exe**. MATLAB:
**D:\Program Files\MATLAB\R2025a\bin\matlab.exe**. Use MATLAB -batch and
machine-appropriate parallelism. Preserve the dirty tree; do not commit or push.

Develop a separate functional routing campaign for my master's robot. Optimize
the actuator system's ability to produce joint functions and combinations of
functions: hip flexion with external rotation, hip extension with external
rotation, their internal-rotation counterparts, abduction/adduction, knee
functions, multi-joint functions, and ankle, subtalar, and toe functions.
Retain the existing 27-actuator-per-side design as the baseline. Evaluate
consolidation as an explicit alternative with independently controlled,
pull-only routes, rather than silently changing the actuator count.

## Read and establish the baseline

1. D:\GitHub\Bipedal_Robot\AGENTS.md
2. Commits since 26 October, especially **34a6431d** and its survey3/capability-v3
   result files. Verify branch and source hashes before comparing results.
3. D:\GitHub\Bipedal_Robot\tmp\paper_bridge\capability_table_20261022\Opt_run_BiPulley_survey2_20261023.m
4. D:\GitHub\Bipedal_Robot\tmp\paper_bridge\capability_table_20261022\biPulleySpecsFromMasters_20261023.m
5. D:\GitHub\Bipedal_Robot\tmp\paper_bridge\capability_table_20261022\capability_v3_three_metric_20261027.csv
6. D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27_v2\actuator_map_27.json
7. D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\PULLEY_SWEEP_DIRECTIVE_20261021.md
8. D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\BiPam_pulley2_20261022.m
9. D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_robotbody.osim
10. D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_simbody.osim
11. D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\bpa_muscle.py
12. D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\README.md
13. D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\robotbody_routes.json
14. D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\dynamic_routes.py

The recent survey2 machinery holds master's attachment points fixed and
optimizes resting/tendon lengths and BPA/tackle configuration. The capability
comparison uses single-muscle graded surfaces and documented substitutions.
Reproduce one existing result first. Distinguish transmission improvement,
route relocation, and actuator consolidation in all comparisons.

## Represent mechanical function explicitly

At every admissible pose, calculate each route's complete signed moment-arm
vector across all independent coordinates. Use x anterior, y superior, z right.
Determine physical internal/external rotation signs from the model, including
left mirroring. Do not infer function from a muscle name or one neutral pose.
Record where moment arms change sign or become negligible. Distinguish a
multi-action muscle at one joint from a biarticular muscle crossing two joints.

Create a function-to-route map and target cases. Include pure functions and
combined torque vectors, such as hip flexion plus external rotation. Include
cross-joint functions for hamstring-, rectus-, and gastrocnemius-like routes.
Measure unwanted components: adduction, axial rotation, and lumbar loading
cannot disappear from the objective because the named primary torque passes.
Pure-flexion cases must allow cancellation of unwanted rotational torque by
other routes; combined-function cases must retain the requested rotation.

Human reference targets must come from stock human geometry in gait2392_simbody,
with the reference force assumptions documented. Robotbody and master's routes
provide robot geometry. Avoid double-counting muscles assigned to several
functional groups. Do not construct an impossible target by combining separate
per-axis maxima that the human muscles cannot generate simultaneously.

Use the coupled admissible pose envelope, including the existing graded-domain
mask. Do not extrapolate a sparse 5x5 surface into a verified full 3-D hip domain.
Extend sampling and validate interpolation where new rotation coordinates enter.
Keep the approved treatment of extreme combined poses explicit.

## Optimize attainable torque vectors

For each pose and design, use tau = R(q) F(P,L(q),design), with nonnegative BPA
tensions, pressure bounds, contraction limits, compliance, and slack behavior.
Test whether the collection of independently controlled routes can attain each
requested torque vector with bounded unwanted components. Report the full
residual vector. Maximizing only the projection onto the desired direction is
insufficient when large perpendicular torque is produced.

One pull-only channel has one torque direction at a given pose. It cannot
independently command flexion, extension, and rotation merely because its
moment-arm vector has several components. Preserve antagonist and independent
functional coverage when proposing consolidation. Check attainable sets and
controllability, including near-singular poses, rather than matrix rank alone.

Use a nested calculation: an inner pressure/tension allocation for each target
case, and an outer design search over feasible attachment patches, via points,
rest/tendon lengths, and BPA/transmission choices. Keep a fixed-route baseline
and compare it with route changes under the same task cases and force law.
Separate the objectives for unmet functional torque, useful RoM, unintended
torque, hardware count, and packaging. Report tradeoffs and infeasible cases;
do not hide them in one survival or coverage score.

Follow the standing nBPA/tackle logic: add parallel BPAs for torque shortfall,
add tackle travel for contraction/RoM shortfall, and re-solve lengths for each
configuration. Preserve force = nBPA*F/G and travel = G*BPA contraction. Use the
existing diameter rule and hardware assumptions as the baseline; any diameter
change is a separately labeled candidate. Load current Xi values from the
adopted sources rather than copying old README numbers. Include bracket and
pulley reactions, routing losses, clearance, and attachment feasibility.

## Validate and integrate in stages

Start with a hip demonstrator spanning pure flexion/extension, internal/external
rotation, and selected combinations, retaining sufficient independent routes.
Show signed torque-vector maps, feasible targets, pressure/strain limits, and
unwanted components before extending the campaign down the leg.

Check analytic length derivatives against constrained finite differences.
Preserve moving and conditional path points and closed four-bar knee coupling.
Verify coordinate/body transforms, source-model hashes, actuator energy/travel
consistency, and fresh-process replays. Add physically appropriate self-collision
and inter-leg clearance checks before a candidate can enter gait training.
For ordinary stance and stepping, reject crossed ankle ordering, overlapping
leg geometry, uncontrolled foot heading, and invalid knee motion.

The current AH001 replay has crossed legs and is not a valid gait target.
Its 65% support replay uses ideal reference tensions, not the optimized BPA
hardware. Import a selected design's actual dimensions, transmission, force
law, and pressure limits into a separately named MuJoCo candidate and compare
static torques with the design calculation before attempting neural training.
The current CAD feet are welded: assess ankle/subtalar/toe capability in an
explicitly articulated design model, and do not claim that capability for the
welded simulation.

Export a functional map suitable for PF modules containing two or four cells
per module. Keep the mechanical design and neural architecture decisions
separate. Do not silently reuse the current four-shared-phase-cell s3k wrapper
or call its MN current clamp an independently functioning PF controller.
Preserve the existing controller as a diagnostic reference.

Put new sources, configuration, and reports under
**D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Functional_Routing_20261027**.
Checkpoint each design, save all knobs and source hashes, and preserve the
existing campaign results. Use verified optimization/simulation speedups only.
First deliver the reproduced baseline, functional map, tested objective,
hip demonstrator, and a measured run-budget estimate. Proceed to the full
campaign only once those checks pass. Do not edit dissertation prose or send
messages to other sessions as part of this work.
