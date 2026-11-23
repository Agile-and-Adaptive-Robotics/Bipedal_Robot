# Fitting joint angles and stance on the human model

The three objectives are the mean of six absolute joint-angle waveform RMSEs in degrees, the right-foot stance error, and the left-foot stance error. The experimental target is 61% stance on each side. Simulated stance is the fraction of all samples in the specified interval for which that foot's upward ground reaction exceeds 50 N. The waveform comparison uses each leg's own complete loading cycles, from one rising crossing to the next, and its recorded OpenSim hip, knee, and ankle angles on the original phase grid. Each angular score retains the mean posture and recorded phase.

The study uses the two-level nonspiking spinal feedback controller previously fitted to Gait2392. Its fitted stage-4 configuration compiles to 606 neurons, 386 scalar inputs, and 6,908 synaptic connections with the 92-muscle human model. The Gait2392 body, original support and ground-contact settings, 2 ms integration step, shared 2.0773198435 nA descending drive, and the four source-parameter compliant-tendon muscles are fixed throughout. The fit can change supplementary left/right motor drive, the existing neural and premotor weights, and timing. The experimental joint trajectory is used to score recorded movements after forward simulation.

The search requires a finite, complete run, a minimum whole-body center-of-mass height of 0.62 m, no more than 40° absolute pelvic tilt, and at least six **foot-loading** cycles per leg. A target hit requires mean error below 17.78° and both stance fractions between 0.58 and 0.64. The narrower check uses 0.59 through 0.63. The exact single-counted frozen comparison was 17.7779202°.

## Measured outcomes

The fit could match the right and left stance targets or improve the six-joint angle error. None of the 531 calibrated evaluations achieved both together. Of the 428 evaluations that met the physical and bilateral cycle requirements, 27 held both stance fractions within three percentage points of 61%, and 13 held them within two points. Three improved on the fixed controller’s exact 17.7779° angle error.

The entries give the mean error over six joint angles followed by the right and left stance fractions. Each schedule used a fresh process.

| Selected controller | 16 s | 22 s | 28 s |
| --- | --- | --- | --- |
| Smallest angle error (recruitment 38) | 16.98°, 80.8%/62.7% | 21.12°, 81.3%/56.7% | 21.53°, 77.4%/57.7% |
| Smallest error within ±3 percentage points (hip/drive bracket 26) | 18.48°, 62.3%/61.4% | 20.89°, 63.2%/53.7% | 21.20°, 69.3%/62.3% |
| Smallest combined shortfall (recruitment 17) | 18.16°, 64.3%/63.7% | 21.78°, 80.6%/57.5% | 21.07°, 71.5%/63.9% |

The same controller had the smallest angle error within both stance intervals.

For the stance-constrained selection, the constant-drive windows gave 18.52°, 77.5%/75.4% during 2–13 s and 21.87°, 74.8%/86.1% during 5–15 s. With independent ±1% changes to the fitted parameters, 8/8 trials retained at least six loading cycles per foot. 1/8 retained both stance targets, and 0/8 retained both stance targets and the angle target.

The parameter files are in `selected/`. The complete measurements, candidate settings, nondominated triples, and independent checks are in `trials.json.gz`, `front.csv`, `outcomes.json`, and `verification.json`. The two compact trajectory files include foot reaction forces and angle traces for the 16 s angle- and stance-selected controllers.

## Running a fit

Run from `Code/MuJoCo_SNS/spinal` in the `myo` environment. Use separate output directories for distinct studies.

The general calibrated fit used Optuna's constrained multivariate TPE sampler with three objectives, 24 startup trials, and seed 241005.

```text
python physics_fidelity/harness_sensitivity/pareto_feedback.py fit --force-convention physical_50_N --sampler tpe --jobs 7 --budget 160 --output OUT/physical_fit
python physics_fidelity/harness_sensitivity/refine_regional_feedback.py --force-convention physical_50_N --seed-search OUT/physical_fit --output OUT/regional --jobs 6 --budget 120
python physics_fidelity/harness_sensitivity/fit_swing_window.py --output OUT/swing --jobs 4 --budget 112
python physics_fidelity/harness_sensitivity/fit_premotor_split.py --parents physics_fidelity/harness_sensitivity/pareto_reference/split_parents.json --output OUT/split --jobs 6 --budget 96
python physics_fidelity/harness_sensitivity/bracket_hip_asymmetry.py --parent physics_fidelity/harness_sensitivity/pareto_reference/bracket_parent.json --output OUT/bracket --jobs 4
```

The general fit changes the existing pattern-formation drive to knee and ankle motor neurons, postural weights, descending currents, and premotor period. The regional search can change homonymous spindle excitation, force-dependent inhibition, and tonic input in specified muscle groups. The timing search uses the original contact-reset swing signal and separately varies its raised-cosine onset and offset, hip-extensor drive in stance, and knee-flexor drive in swing. The final search allows independent amplitudes for the existing premotor command to ankle dorsiflexors, plantarflexors, knee flexors, and hip flexors. All of these settings leave the mechanical model and externally applied support fixed.

The follow-up centered on the controller closest to the combined targets. It varied the right supplementary drive, left supplementary drive, and right hip abductor spindle excitation. Right hip adductor excitation varied reciprocally. All other weights stayed fixed. The design includes all 27 combinations at 98%, 100%, and 102% of each center value and 16 additional combinations within ±3% selected by a scrambled Sobol sequence.

`trials.json.gz` records the exact parameters and measured objectives of each evaluation, including the seeded and manually nominated settings. The adaptive optimizers incorporate completed evaluations as workers finish, so changing the worker count or completion order can change later proposals. To repeat an individual evaluation, extract its `parameters` object and run `pareto_feedback.py replay --parameters PARAMETERS.json --output OUT/repeat`.

Prior NSGA-II and premotor comparison trials used `--force-convention historical_stored_50`. That scorer thresholded the sum of the original heel and toe controller channels at 50. On this flat-ground model the physical force at that threshold was 25 N. Those trials document the sensitivity to the older convention. Their objectives cannot be pooled with the physical 50 N trials.

## Repeating a fixed candidate

```text
python physics_fidelity/harness_sensitivity/verify_pareto_feedback.py FIXED_PARAMETERS.json --output OUT/verification --jobs 3 --extended-holdout
```

The verifier uses a fresh process for each 16 s evaluation, 22 s production schedule, optional 28 s extension, and each of eight independent ±1% scalar perturbations. It saves the controller settings, compact trajectory, calibrated right/left foot forces, active muscle check, and absolute and historical-threshold results. The historical comparison windows are 2–16 s and 5–22 s. They include the original ramp-down and return to standing. Each validation also reports the constant-drive interval, 2–13 s for the evaluation, 5–15 s for production, and 5–21 s for the extended schedule. End times of the constant-drive windows are excluded.

`AARL_FOOT_FORCE_NPZ` enables calibrated observations for the spinal runner. `test_foot_loading.py` checks that two distinct contacts under the same foot contribute once each, that reversing MuJoCo's contact geometry order preserves upward force, and that foot-to-foot contacts are excluded. The original feedback still receives its historical sensor signals. A fixed trial replay confirmed that adding the observer preserved the entire earlier trajectory byte for byte.

The optional `pareto_feedback.py replay --foot-sensor-mode anatomical` condition feeds each calcaneus contact only to its heel sensor and each toe contact only to its toe sensor. Their sum supplies the physical load once. The saved raw `contact` field retains the archived scale for compatibility. `foot_ground.npz` contains the independently calibrated reaction. With the saved weights and support held fixed, the frozen controller produced 4 right and 0 left loading cycles with anatomically separated inputs, compared with 11 and 10 using its original inputs. The corresponding counts for recruitment trial 17 were 6 and 0, compared with 8 and 11. Both comparisons are in `foot_input_comparison.json`.
