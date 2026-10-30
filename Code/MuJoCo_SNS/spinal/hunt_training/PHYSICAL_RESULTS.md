# Corrected physical curriculum: measured outcomes

Machine: easteregg2. Native neural outputs, MuJoCo 2.3.7, implicitfast, 2 ms, sparse preallocated Numba SNS backend.

These are stage tests and diagnostic continuations. Later tests do not certify that earlier failed stages graduated.

| Campaign / case | Stage | Pass | Knee cycles R/L | Foot load R/L (BW) | Spring support (BW) | Lateral restraint RMS (BW) |
|---|---:|---|---|---|---:|---:|
| corrected_curriculum_v1 / s1 | 1 | False | [0, 0] | 0.000/0.000 | 0.999 | nan |
| corrected_curriculum_v1 / s2 | 2 | False | [0, 0] | 0.000/0.000 | 0.999 | nan |
| corrected_curriculum_v1 / s3 | 3 | False | standing | 0.002/0.059 | 0.938 | nan |
| corrected_curriculum_v1 / s4_lower20 | 4 | False | [4, 0] | 0.171/0.215 | 0.616 | nan |
| corrected_curriculum_v1 / s4_lower40 | 4 | False | [3, 0] | 0.218/0.207 | 0.598 | nan |
| corrected_curriculum_v1 / s4_restricted | 4 | False | [3, 0] | 0.300/0.106 | 0.539 | nan |
| corrected_curriculum_v2 / s1_init | 1 | True | [4, 3] | 0.000/0.000 | 0.999 | 0.023 |
| corrected_curriculum_v2 / s2_init | 2 | False | [0, 0] | 0.000/0.000 | 0.999 | 0.013 |
| corrected_curriculum_v2 / s4_ext0.1_ref0.0 | 4 | False | [4, 6] | 0.048/0.041 | 0.901 | 0.012 |
| corrected_curriculum_v2 / s4_ext0.1_ref1.0 | 4 | False | [4, 6] | 0.054/0.071 | 0.869 | 0.007 |
| corrected_curriculum_v2 / s4_ext0.25_ref0.5 | 4 | False | [4, 6] | 0.059/0.066 | 0.871 | 0.008 |
| corrected_curriculum_v2 / s4_ext0.5_ref0.0 | 4 | False | [4, 6] | 0.073/0.055 | 0.864 | 0.012 |
| corrected_curriculum_v2 / s4_ext0.5_ref1.0 | 4 | False | [4, 6] | 0.083/0.108 | 0.808 | 0.010 |
| corrected_curriculum_v3 / s3_height-0.04 | 3 | False | standing | 0.053/0.106 | 0.839 | 0.011 |
| corrected_curriculum_v3 / s3_height-0.06 | 3 | False | standing | 0.281/0.272 | 0.449 | 0.026 |
| corrected_curriculum_v3 / s3_height-0.08 | 3 | False | standing | 0.259/0.255 | 0.485 | 0.022 |
| corrected_curriculum_v3 / s4_height-0.04 | 4 | False | [4, 6] | 0.087/0.100 | 0.811 | 0.009 |
| corrected_curriculum_v3 / s4_height-0.06 | 4 | False | [4, 6] | 0.090/0.099 | 0.767 | 0.007 |
| corrected_curriculum_v3 / s4_height-0.08 | 4 | False | [4, 4] | 0.090/0.130 | 0.741 | 0.009 |
| corrected_curriculum_v3 / s5_baseline | 5 | False | [4, 0] | 0.089/0.244 | 0.666 | 0.015 |
| corrected_curriculum_v3 / s5_push20N | 5 | False | [4, 0] | 0.089/0.244 | 0.666 | 0.015 |
| corrected_curriculum_v3 / s5_restraint_half | 5 | False | [4, 0] | 0.088/0.242 | 0.663 | 0.017 |
| corrected_curriculum_v3 / s5_vertical_quarter | 5 | False | [4, 0] | 0.089/0.243 | 0.665 | 0.015 |

Vertical loads in current v2 tests count each foot-ground contact once. The legacy runner assigns each contact to both heel and toe, doubling its force sum. Its original contact metrics remain in raw files; current ground-load verification uses independent world-vertical forces. v1 s1/s2 were air tests; subsequent cases used the corrected force measurement.

The harness wrench is computed from the actual spring and damper forces on pelvis coordinates and projected into native world axes. It is a virtual instrumentation calculation, not a physical single-point six-axis load-cell attachment. Distributed hip/lumbar restrictions remain additional assistance and are explicitly configured.

Native world x is anterior, y is left, and z is superior; anatomical right force is negative native y. Record force and torque with contact timing; a large harness lateral reaction alone cannot identify which foot generated it.

Effort: mean activation squared and excitation squared are recorded. These are dimensionless regularizers, not calibrated metabolic energy, and never override movement, contact, or crossing gates.

Speed-dependent feedback evidence: Hurteau et al. 2017 reports nonlinear cutaneous-reflex modulation in spinal cats, strongest at intermediate speeds, not a universal monotonic decline (https://doi.org/10.1523/JNEUROSCI.3042-16.2017). Mayer et al. 2018 supports a role for spindle feedback in speed-dependent extensor recruitment in mice (https://doi.org/10.1152/jn.00250.2018). These do not establish that Ib becomes unnecessary or that Ia and II should both increase identically.

Separate pathway ablations across achieved speed, assistance, and phase are required before attributing causal importance. A change in signal amplitude or fitted conductance alone is insufficient.
