# Measured kinematic statistics

Machine: easteregg2. Source: `D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_ankle_selected_repeat_20261104\stage2`.

Comparison role: **diagnostic_only_air**. Grounded IK differences do not reject air stepping.

Window: [5.0, 14.0] s; 3 complete right-knee-anchored cycles.

## Frequency and period

| Side | Intervals n | Frequency mean ± SD (Hz) | Reference (Hz) | Frequency error (%) | Period mean ± SD (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| R | 3 | 0.499 ± 0.000 | 0.811 | -38.453 | 2.004 ± 0.000 |
| L | 3 | 0.498 ± 0.000 | 0.824 | -39.557 | 2.007 ± 0.001 |

## Total amplitude and extrema

Amplitude means peak-to-peak excursion, not half-amplitude. Cycle SD is sample SD.

| Joint | Reference amplitude (deg) | Whole-window amplitude (deg) | Per-cycle amplitude mean ± SD (deg) | Mean amplitude error (%) | Window min / max (deg) |
| --- | ---: | ---: | ---: | ---: | ---: |
| hip_flexion_r | 43.281 | 25.431 | 25.420 ± 0.006 | -41.268 | -15.199 / 10.232 |
| knee_angle_r | 70.432 | 20.450 | 20.397 ± 0.018 | -71.040 | -22.383 / -1.933 |
| ankle_angle_r | 23.101 | 26.684 | 26.662 ± 0.004 | 15.413 | -31.371 / -4.687 |
| hip_flexion_l | 44.109 | 26.523 | 26.405 ± 0.005 | -40.137 | -16.266 / 10.257 |
| knee_angle_l | 70.618 | 22.840 | 22.756 ± 0.036 | -67.775 | -24.779 / -1.939 |
| ankle_angle_l | 23.937 | 27.991 | 27.910 ± 0.045 | 16.595 | -34.317 / -6.326 |

## Joint frequency content

Leg cadence alone can hide fast ankle oscillations. These spectra use the unwarped 5–14s window; bin width is approximately 0.111Hz. The 3Hz split is diagnostic.

| Joint | Dominant spectral frequency (Hz) | Strongest peak above 3Hz (Hz) | Power above 3Hz (%) | High-frequency amplitude RMS (deg) |
| --- | ---: | ---: | ---: | ---: |
| hip_flexion_r | 0.444 | 3.000 | 0.135 | 0.314 |
| knee_angle_r | 0.444 | 3.000 | 1.429 | 0.933 |
| ankle_angle_r | 0.444 | 6.000 | 9.118 | 2.721 |
| hip_flexion_l | 0.444 | 3.000 | 0.061 | 0.208 |
| knee_angle_l | 0.444 | 3.000 | 0.476 | 0.547 |
| ankle_angle_l | 0.444 | 6.000 | 12.423 | 3.151 |

## Phase error and maximum error

Phase is the fundamental lag versus reference, with one shared knee anchor. Positive means lag; circular statistics handle wraparound. Maximum error uses **every cycle**, not just the averaged waveform.

| Joint | Phase cycles n | Circular phase mean ± SD (deg) | Maximum absolute phase error (deg) | Worst absolute waveform error (deg) | Per-cycle maximum error mean ± SD (deg) | Mean bias (deg) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| hip_flexion_r | 3 | -88.914 ± 0.069 | 89.003 | 32.194 | 32.182 ± 0.015 | -6.489 |
| knee_angle_r | 3 | 9.895 ± 0.080 | 9.997 | 47.595 | 47.591 ± 0.005 | 15.171 |
| ankle_angle_r | 3 | 98.861 ± 0.071 | 98.938 | 40.070 | 40.065 ± 0.006 | -25.511 |
| hip_flexion_l | 3 | -99.876 ± 0.435 | 100.408 | 33.561 | 33.490 ± 0.070 | -6.319 |
| knee_angle_l | 3 | 4.656 ± 0.419 | 5.172 | 47.414 | 47.299 ± 0.116 | 9.840 |
| ankle_angle_l | 3 | 74.119 ± 0.409 | 74.620 | 39.644 | 39.587 ± 0.054 | -26.212 |

## Additional waveform diagnostics

| Joint | Mean absolute error (deg) | Mean-cycle RMSE (deg) |
| --- | ---: | ---: |
| hip_flexion_r | 16.391 | 18.590 |
| knee_angle_r | 17.470 | 22.749 |
| ankle_angle_r | 25.671 | 27.828 |
| hip_flexion_l | 17.507 | 19.785 |
| knee_angle_l | 14.696 | 20.166 |
| ankle_angle_l | 26.212 | 27.844 |

Descriptive within-run cycle statistics; SD is sample SD, not independent-subject uncertainty. Maximum waveform error includes every retained cycle, using one shared right-knee anchor; cadence error remains separate.

These few cycles are correlated observations from one run. They do not justify independent-sample significance tests or population confidence intervals.

Reference: subject01_walk1: existing calibration IK/GRF reference. See KINEMATIC_RESPONSE.json for all unrounded values, gains, per-side frequency errors, and provenance.
