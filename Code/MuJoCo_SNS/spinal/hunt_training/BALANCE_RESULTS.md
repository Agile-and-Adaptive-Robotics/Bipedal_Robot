# Stage 3 task-balance implementation and measurements

Machine: easteregg2. Eight added SNS cells: four directional task populations and four inhibitory relays. Controller total: 1,029 neurons, 424 inputs. No runtime MN target cancellation.

CoM displacement/velocity use 100 ms delay and a 40 ms velocity/acceleration filter by default. Acceleration and somatic/orientation cue weighting are separate experimental options. Ground-truth simulation signals are sensory surrogates, not a model of receptor transduction or adaptive sensory reweighting.

Recruitment is calibrated once at the initial pose from finite muscle-activation effects on CoM acceleration under contact. It is held fixed during each trial. Locomotor RG/PF fits are preserved.

| Campaign/case | Stage | Pass | Sway (m) | Horizontal harness RMS / BW | Foot load sum / BW | Wall (s) |
|---|---:|---|---:|---:|---:|---:|
| balance_20261029_force_support/s3_support0.4 | 3 | False | 0.0665 | 0.1417 | 0.599 | 20.99 |
| balance_20261029_force_support/s3_support0.6 | 3 | False | 0.0339 | 0.0634 | 0.399 | 20.53 |
| balance_20261029_force_support/s4_support0.4 | 4 | False | nan | 0.1072 | 0.599 | 31.67 |
| balance_20261029_force_support/s4_support0.6 | 4 | False | nan | 0.0724 | 0.399 | 31.91 |
| balance_20261029_v1/benchmark_native | 3 | False | 0.2598 | 0.0386 | 0.416 | 122.12 |
| balance_20261029_v1/benchmark_sparse | 3 | False | 0.2598 | 0.0386 | 0.416 | 13.25 |
| balance_20261029_v1/curriculum_stage_1 | 1 | True | nan | 0.0274 | 0.000 | 25.18 |
| balance_20261029_v1/curriculum_stage_2 | 2 | False | nan | 0.0418 | 0.000 | 30.42 |
| balance_20261029_v1/curriculum_stage_4 | 4 | False | nan | 0.0363 | 0.190 | 37.88 |
| balance_20261029_v1/curriculum_stage_5 | 5 | False | nan | 0.0980 | 0.314 | 34.50 |
| balance_20261029_v1/holdout_forward | 3 | False | 0.2598 | 0.0386 | 0.416 | 18.57 |
| balance_20261029_v1/holdout_lateral | 3 | False | 0.2598 | 0.0386 | 0.416 | 18.28 |
| balance_20261029_v1/less_restraint | 3 | False | 0.2598 | 0.0362 | 0.413 | 17.93 |
| balance_20261029_v1/less_vertical_stiffness | 3 | False | 0.2598 | 0.0392 | 0.414 | 23.97 |
| balance_20261029_v1/tune_00 | 3 | False | 0.2598 | 0.0386 | 0.416 | 14.03 |
| balance_20261029_v1/tune_01 | 3 | False | 0.2752 | 0.0118 | 0.377 | 12.33 |
| balance_20261029_v1/tune_02 | 3 | False | 0.2875 | 0.0294 | 0.309 | 13.13 |
| balance_20261029_v1/tune_03 | 3 | False | 0.2923 | 0.0603 | 0.236 | 12.79 |
| balance_20261029_v1/tune_04 | 3 | False | 0.2993 | 0.0571 | 0.247 | 13.61 |
| balance_20261029_v1/tune_05 | 3 | False | 0.2876 | 0.0294 | 0.309 | 14.98 |
| balance_20261029_v1/tune_06 | 3 | False | 0.2928 | 0.0601 | 0.237 | 13.05 |
| balance_20261029_v1/tune_07 | 3 | False | 0.2941 | 0.0587 | 0.237 | 12.77 |
| balance_20261029_v2/benchmark_native | 3 | False | 0.1404 | 0.0244 | 0.645 | 160.90 |
| balance_20261029_v2/benchmark_sparse | 3 | False | 0.1404 | 0.0244 | 0.645 | 17.86 |
| balance_20261029_v2/curriculum_stage_1 | 1 | False | nan | 0.0259 | 0.000 | 18.71 |
| balance_20261029_v2/curriculum_stage_2 | 2 | False | nan | 0.0410 | 0.000 | 26.28 |
| balance_20261029_v2/curriculum_stage_4 | 4 | False | nan | 0.1171 | 0.481 | 35.62 |
| balance_20261029_v2/curriculum_stage_5 | 5 | False | nan | 0.1184 | 0.572 | 38.24 |
| balance_20261029_v2/disabled | 3 | False | 0.1295 | 0.0096 | 0.521 | 16.14 |
| balance_20261029_v2/heldout_forward | 3 | False | 0.1404 | 0.0246 | 0.645 | 20.14 |
| balance_20261029_v2/heldout_lateral | 3 | False | 0.1404 | 0.0244 | 0.645 | 17.52 |
| balance_20261029_v2/inherited | 3 | True | 0.0474 | 0.1034 | 0.614 | 19.78 |
| balance_20261029_v2/less_restraint | 3 | False | 0.1405 | 0.0200 | 0.638 | 18.40 |
| balance_20261029_v2/pd_low | 3 | False | 0.1404 | 0.0244 | 0.645 | 18.14 |
| balance_20261029_v2/pda | 3 | False | 0.1405 | 0.0244 | 0.645 | 18.82 |
| balance_20261029_v3/acceleration | 3 | True | 0.0475 | 0.1062 | 0.621 | 19.38 |
| balance_20261029_v3/benchmark_native | 3 | True | 0.0474 | 0.1038 | 0.615 | 131.27 |
| balance_20261029_v3/benchmark_sparse | 3 | True | 0.0474 | 0.1038 | 0.615 | 19.48 |
| balance_20261029_v3/curriculum_air_stage_1 | 1 | True | nan | 0.0274 | 0.000 | 18.66 |
| balance_20261029_v3/curriculum_air_stage_2 | 2 | False | nan | 0.0418 | 0.000 | 24.00 |
| balance_20261029_v3/curriculum_stage_1 | 1 | False | nan | 0.0259 | 0.000 | 19.00 |
| balance_20261029_v3/curriculum_stage_2 | 2 | False | nan | 0.0410 | 0.000 | 26.72 |
| balance_20261029_v3/curriculum_stage_4 | 4 | False | nan | 0.1164 | 0.360 | 28.60 |
| balance_20261029_v3/curriculum_stage_5 | 5 | False | nan | 0.1128 | 0.428 | 27.15 |
| balance_20261029_v3/disabled | 3 | True | 0.0474 | 0.1034 | 0.614 | 23.84 |
| balance_20261029_v3/heldout_forward | 3 | True | 0.0474 | 0.1058 | 0.616 | 18.53 |
| balance_20261029_v3/heldout_lateral | 3 | True | 0.0474 | 0.1038 | 0.615 | 20.60 |
| balance_20261029_v3/high | 3 | True | 0.0479 | 0.1070 | 0.636 | 21.39 |
| balance_20261029_v3/less_restraint | 3 | True | 0.0474 | 0.0951 | 0.602 | 19.17 |
| balance_20261029_v3/low | 3 | True | 0.0474 | 0.1038 | 0.615 | 19.79 |
| balance_20261029_v3/medium | 3 | True | 0.0475 | 0.1059 | 0.622 | 19.43 |

## Matched backend comparisons

Same fitted model, network topology, initialization, assistance, duration, and logging. Native uses the dense SNS fixed-tau equations; accelerated uses the sparse preallocated Numba equations. Background trial load varies, so elapsed times are observations, not controlled hardware benchmarks. Whole-run times include model/network setup, standing solve, calibration, physics, sensing, and file output. Neural-step times also include live feature and conductance updates.

```json
{
  "campaign": "balance_20261029_v1",
  "simulated_s": 8.0,
  "sparse_wall_s": 13.252808100078255,
  "native_wall_s": 122.11963820003439,
  "end_to_end_ratio": 9.214623593569826,
  "sparse_neural_step_s": 4.149503400316462,
  "native_neural_step_s": 110.48374320252333,
  "maximum_array_differences": {
    "qfull": 0.0,
    "ctrl": 0.0,
    "com": 0.0,
    "neuro": 0.0
  }
}
```

```json
{
  "campaign": "balance_20261029_v2",
  "simulated_s": 8.0,
  "sparse_wall_s": 17.858282100060023,
  "native_wall_s": 160.89526799996383,
  "end_to_end_ratio": 9.009560219648622,
  "sparse_neural_step_s": 5.862119498196989,
  "native_neural_step_s": 146.82377060048748,
  "maximum_array_differences": {
    "qfull": 0.0,
    "ctrl": 0.0,
    "com": 0.0,
    "neuro": 0.0
  }
}
```

```json
{
  "campaign": "balance_20261029_v3",
  "simulated_s": 8.0,
  "sparse_wall_s": 19.477941900026053,
  "native_wall_s": 131.26625230000354,
  "end_to_end_ratio": 6.739225990802856,
  "sparse_neural_step_s": 6.593732902780175,
  "native_neural_step_s": 118.07478609634563,
  "maximum_array_differences": {
    "qfull": 0.0,
    "ctrl": 0.0,
    "com": 0.0,
    "neuro": 0.0
  }
}
```

v1 tests replaced all inherited BAL inputs. v2 used a neutral standing pose; the inherited controller passes its initial screen, while replacement task feedback does not. v3 retains the inherited controller and tests additive task feedback against an identical-topology disabled-pathway control. Distributed harness restrictions remain assistance. No independent balance or curriculum graduation follows merely from a stage-3 pass.

The v3 winner file selects the best **nonzero** task-feedback candidate for comparison, not an established improvement over the disabled-task control. The latter has slightly lower harness reaction and CoM displacement in these initial tests. Do not promote the added pathway as beneficial on this evidence.

Prescribed-support cases explicitly remove the vertical rig spring/damper and apply a constant 60% or 40% body-weight force. They differ from stiffness-reduction tests. Both stages 3 and 4 remain failed under those tested support settings; inspect symmetry, positional error, phase, and contact criteria rather than using force support alone.

Air stages named curriculum_air_stage_* retain their calibrated air initialization. Earlier curriculum_stage_1/2 neutral-standing-start probes are retained separately and do not replace the air-stage protocol.

Contact sensing: three forces per foot are already logged. Current neural contact/timing pathways mainly receive scalar contact/load, not the full 3D force vector. COP and ankle-referenced contact moments are useful next measurements, but are not currently neural inputs.

Suggested lineage label: biw2lvar-Hunt-PT experimental. This denotes W2L-variant ancestry plus biarticular integration, 3D joint PF channels, Hunt fitting, and the task-balance pathway. It is not an unchanged W2L model, a literature-established name, or a claim of validated walking.
