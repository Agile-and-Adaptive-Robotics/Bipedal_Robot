# Repeating the human-controller experiments

The four-muscle condition uses the source Gait2392 parameters for right and left rectus femoris and biceps femoris short head. At each coupled step the solver takes activation and path length from the live MuJoCo model, solves fiber–tendon equilibrium with active and passive tension, and applies the resulting tendon force through the same signed moment arm. MuJoCo continues to integrate the 88 other muscles, gravity, the joints, contact, and the recorded pelvic support. The matched native condition keeps all 92 MuJoCo muscles. Both conditions start with the same controller parameters and body-support settings.

Run the commands below from `Code/MuJoCo_SNS/spinal` in the project's `myo` Python environment. `OUT` denotes an empty results location. Each simulation records the muscle assignment, common starting state, and completed physics intervals.

The [joint-angle and stance study](pareto_reference/README.md) fits three independent objectives under the same mechanical and support conditions. It includes calibrated foot-force measurements, the complete set of evaluated controller settings, and fresh-process verification.

## Reference-controller tuning

`fit_reference_feedback.py` runs the original four rounds of feedback tuning. Each worker uses one BLAS thread. The four worker seeds in each round are recorded in `reference_feedback/trials.csv`. The 16 s evaluations score the trajectory beginning at 2 s. The 22 s production schedule uses the original score window beginning at 5 s. The shared descending drive remains 2.0773198435110425 nA while the four supplementary left/right drive channels can vary.

~~~text
python physics_fidelity/harness_sensitivity/fit_reference_feedback.py --output OUT/study --jobs 4
python physics_fidelity/harness_sensitivity/replay_feedback_trial.py 5 83 --output OUT/exact_trial
python physics_fidelity/harness_sensitivity/verify_reference_feedback.py physics_fidelity/harness_sensitivity/reference_feedback/selected.json --arm native --output OUT/matched_native
python physics_fidelity/harness_sensitivity/verify_reference_feedback.py physics_fidelity/harness_sensitivity/reference_feedback/selected.json --output OUT/seed3000 --basin-seed 3000
python physics_fidelity/harness_sensitivity/verify_reference_feedback.py physics_fidelity/harness_sensitivity/reference_feedback/selected.json --output OUT/zero_supplementary --scale 0
~~~

The search ranks finite candidates with at least three qualifying cycles on each side. The separate drive check requires at least six cycles on each side and a stance fraction of 0.40–0.75 on both sides. The eight perturbations use seeds 3000–3007 and the original three-cycle criterion. The verification output identifies which leg used a foot-loading onset and which required the extensor rhythm-generator fallback. A score under pelvic support measures the simulated trajectory under that support.

`reference_feedback/study.json` records the trial budgets, final weights, independent repetitions, drive checks, and eight-seed results. Each row of `trials.csv` preserves the real proposal and seed plus the resulting score and leg-specific metrics. Concurrent TPE proposals depend on the order the workers return. `replay_feedback_trial.py` uses the recorded proposal to reconstruct a specific trial directly. The historical search classified left-only trajectories as infeasible. The verification command reports their actual leg-specific cycle counts and computes the score for the active leg.

The two matched trajectory files in `reference_feedback` retain every 2 ms sample needed to recompute the historical kinematic score. The arrays contain time in seconds (`t`), joint angles in degrees (`q`), sampled neural channels (`neuro`), the historical right and left heel-plus-toe controller channels (`contact`), body center of mass in meters (`com`), and the lowest selected foot-body center height on each side in meters (`footz`). The runner assigned each foot contact normal to both its heel and toe channel. On this level-ground model the archived sum therefore equals twice the measured upward ground reaction. `key_joints` and `neuro_names` identify the angle and neural columns. The full and compact trace hashes appear in `study.json`.

`reference_angle_errors.py RECORDING --force-convention physical --output OUT/angles.json` uses an actual 50 N foot-loading threshold on the archived level-ground trajectories. The `historical` option reproduces the original threshold of 50 on the summed controller channel, equivalent to 25 N on this model. The experimental reference retains its physical 50 N threshold in both reductions. Each leg uses its own loading onset. Absolute hip, knee, and ankle trajectories are compared throughout the unshifted cycle. The pre-Pareto `reference_feedback/angle_comparison.json` records the historical reduction.

Set `AARL_FOOT_FORCE_NPZ` to a new output file to record calibrated `normal_N`, `upward_N`, and `contact_count` for the right and left feet on a new Gait2392 trial. The observer includes each foot–ground contact once, including simultaneous heel and toe contacts, and excludes contact between feet. It also records `t` on the runner's existing sampling clock. Evaluate the calibrated recording with `reference_angle_errors.py RECORDING --foot-loads CALIBRATED --output OUT/angles.json`. The feedback channels and the default trajectory format remain available for exact historical reproductions.

## BiW2Lvar+ and syn6 air and ground conditions

Supply the original local `ControllerTraining_20261104` directory as `STUDY`. The runner reuses the archived neural parameters, commanded currents, initial angles, feedback gains, support, and knee guard. A separate `--checkpoint` argument applies a newly fitted D checkpoint. The path names beginning with `hunt_` are the original filenames for the BiW2Lvar+ studies.

~~~text
python physics_fidelity/harness_sensitivity/repeat_body_condition.py STUDY hunt_pose_air_validation_v3/stage2 --output OUT/biw_air
python physics_fidelity/harness_sensitivity/repeat_body_condition.py STUDY syn6_flexor_handover_20261105/predict150_stage2 --output OUT/syn6_air
python physics_fidelity/harness_sensitivity/repeat_body_condition.py STUDY syn6_flexor_handover_20261105/predict150_stage4 --checkpoint STUDY/four_muscle_refit/syn6/waveform_D --output OUT/retrained_ground
~~~

`condition.json` identifies the sources and any knee-guard stop. The body folder contains the aligned joint, neural, actuator, and support measurements generated by the archive's own evaluator. For a new controller, `coupled_muscle_context.py` can also run an existing Gait2392 simulation script directly.

## Kinematics and recruitment

`gait_validate_all.py --output OUT/library` scores the four older reference configurations. The scorer accepts the joint angles already stored in degrees. The library contains 43 records. The left side of `ong_selfsel_Init200` is degenerate, leaving 42 bilateral comparisons. Only the leg that produced a cycle has a cycle-averaged duty and waveform.

`quasistatic_recruitment.py --input INPUT --output OUT/activation` recalculates estimated activation on the anatomically signed `bsolve_ik.py --fk-signs` recording. This comparison retains the original joint torques, moment arms, and held velocities and applies the source force calculation to the four selected muscles. The separate prescribed-trajectory controller training uses the walking history and its original 0.40-body-weight support protocol.
