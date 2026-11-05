# Measured kinematic statistics

Machine: easteregg2. Source: `D:\Users\Ben Bolen\ControllerTraining_20261104\s3k_repaired_body_20261104\expression_0.2\stage2`.

Comparison role: **diagnostic_only_air**. Grounded IK differences do not reject air stepping.

Window: [5.0, 14.0] s; 6 complete right-knee-anchored cycles.

## Frequency and period

| Side | Intervals n | Frequency mean ± SD (Hz) | Reference (Hz) | Frequency error (%) | Period mean ± SD (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| R | 6 | 0.723 ± 0.006 | 0.811 | -10.769 | 1.382 ± 0.012 |
| L | 5 | 0.693 ± 0.010 | 0.824 | -15.893 | 1.443 ± 0.021 |

## Total amplitude and extrema

Amplitude means peak-to-peak excursion, not half-amplitude. Cycle SD is sample SD.

| Joint | Reference amplitude (deg) | Whole-window amplitude (deg) | Per-cycle amplitude mean ± SD (deg) | Mean amplitude error (%) | Window min / max (deg) |
| --- | ---: | ---: | ---: | ---: | ---: |
| hip_flexion_r | 43.281 | 127.628 | 124.621 ± 3.466 | 187.934 | -61.648 / 65.980 |
| knee_angle_r | 70.432 | 28.415 | 25.170 ± 1.769 | -64.263 | -30.152 / -1.737 |
| ankle_angle_r | 23.101 | 30.809 | 30.559 ± 0.271 | 32.286 | -42.101 / -11.292 |
| hip_flexion_l | 44.109 | 142.427 | 136.299 ± 7.390 | 209.005 | -67.749 / 74.678 |
| knee_angle_l | 70.618 | 25.953 | 21.035 ± 2.670 | -70.213 | -27.688 / -1.734 |
| ankle_angle_l | 23.937 | 34.516 | 33.954 ± 0.626 | 41.845 | -46.519 / -12.003 |

## Joint frequency content

Leg cadence alone can hide fast ankle oscillations. These spectra use the unwarped 5–14s window; bin width is approximately 0.111Hz. The 3Hz split is diagnostic.

| Joint | Dominant spectral frequency (Hz) | Strongest peak above 3Hz (Hz) | Power above 3Hz (%) | High-frequency amplitude RMS (deg) |
| --- | ---: | ---: | ---: | ---: |
| hip_flexion_r | 0.667 | 3.556 | 0.005 | 0.305 |
| knee_angle_r | 0.778 | 3.556 | 0.794 | 0.794 |
| ankle_angle_r | 0.667 | 4.333 | 0.949 | 1.067 |
| hip_flexion_l | 0.667 | 3.444 | 0.002 | 0.206 |
| knee_angle_l | 0.667 | 3.444 | 0.560 | 0.560 |
| ankle_angle_l | 0.667 | 3.444 | 1.183 | 1.324 |

## Phase error and maximum error

Phase is the fundamental lag versus reference, with one shared knee anchor. Positive means lag; circular statistics handle wraparound. Maximum error uses **every cycle**, not just the averaged waveform.

| Joint | Phase cycles n | Circular phase mean ± SD (deg) | Maximum absolute phase error (deg) | Worst absolute waveform error (deg) | Per-cycle maximum error mean ± SD (deg) | Mean bias (deg) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| hip_flexion_r | 6 | -158.732 ± 0.609 | 159.871 | 86.987 | 85.538 ± 1.597 | -0.136 |
| knee_angle_r | 6 | -10.247 ± 0.994 | 12.026 | 44.352 | 43.524 ± 1.130 | 16.814 |
| ankle_angle_r | 6 | 104.233 ± 0.748 | 104.945 | 51.116 | 51.098 ± 0.021 | -33.139 |
| hip_flexion_l | 6 | -129.380 ± 16.920 | 150.710 | 93.237 | 89.105 ± 3.748 | 3.646 |
| knee_angle_l | 6 | 22.693 ± 20.642 | 53.487 | 59.460 | 52.729 ± 6.124 | 15.557 |
| ankle_angle_l | 6 | 126.910 ± 17.684 | 154.644 | 60.033 | 56.716 ± 2.267 | -32.492 |

## Additional waveform diagnostics

| Joint | Mean absolute error (deg) | Mean-cycle RMSE (deg) |
| --- | ---: | ---: |
| hip_flexion_r | 51.963 | 57.927 |
| knee_angle_r | 17.010 | 21.346 |
| ankle_angle_r | 33.139 | 35.808 |
| hip_flexion_l | 52.224 | 56.832 |
| knee_angle_l | 17.590 | 23.755 |
| ankle_angle_l | 32.492 | 36.027 |

Descriptive within-run cycle statistics; SD is sample SD, not independent-subject uncertainty. Maximum waveform error includes every retained cycle, using one shared right-knee anchor; cadence error remains separate.

These few cycles are correlated observations from one run. They do not justify independent-sample significance tests or population confidence intervals.

Reference: subject01_walk1: existing calibration IK/GRF reference. See KINEMATIC_RESPONSE.json for all unrounded values, gains, per-side frequency errors, and provenance.
