# Gait2392 OpenSim–MuJoCo strength provenance

Session 1, easteregg2, 30 October 2026. Raw audit:
`runs/strength_provenance_20261030_s1/AUDIT.json`.

## Origin of bilateral differences

The OpenSim input `Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim`
has identical bilateral maximum-isometric-force specifications for soleus
(3549 N), vastus medialis (1294 N), and tibialis anterior (905 N).
The corresponding subject01 inputs checked in this audit agree.

MyoConverter independently optimizes four force-curve parameters for each
muscle. `D:\GitHub\myoConverter\myoconverter\conversion_steps\O2MStep3.py`,
lines 335–347, assigns MuJoCo's force parameter to
`res_opt[3] * OpenSim_Fmax`, then copies the parameters to the passive law.
The other fitted values change normalized length limits and passive force.
The objective is RMS force-map error, normalized by source Fmax, over sampled
joint configurations. This is not a bilateral-constrained fit.

| Soleus parameter | Right | Left |
| --- | ---: | ---: |
| OpenSim Fmax, N | 3549 | 3549 |
| Fitted multiplier | 0.9567481 | 0.89782008 |
| MuJoCo force parameter, N | 3395.50 | 3186.36 |
| Stored normalized fit RMS | 0.03972 | 0.05838 |

Both soleus fits use seven sampled configurations, with independently fitted
normalized length/passive-force parameters. Across 84 cached optimized muscle
fits, the multiplier reproduces the saved XML force parameter within 0.00466 N
(XML rounding). The runtime trunk repairs are a separate intervention and
must not be conflated with these converter fits.

## Interpretation and remaining alignment work

The asymmetry was introduced during conversion, not specified as a biological
strength difference in these source muscles. However, these fitted MuJoCo
parameters compensate for a different muscle model; the force parameter alone
does not establish force/torque agreement at matched poses and activation.
Simply replacing both values with their mean or with OpenSim Fmax may worsen
the fitted force curve. No such correction has been silently applied.

Required comparisons include matched-pose active/passive force and moment arm,
length/velocity dependence, activation transients, tendon compliance, joint
passive moments, contact, and support forces. Converted moving/conditional
pathpoint followers and their stabilizing armature also affect dynamics.
The original Thelen source and MuJoCo rigid-tendon approximation require
explicit comparison; a new series-elastic variant is being investigated
separately. It is not claimed to replicate Hyfydy.

## Neural correction

The earlier independent left/right timing fit was inappropriate for the intended
symmetric controller. `train_symmetric_hunt.py` fits one six-value vector to
both legs. `bilateral.py` pools mirrored regression columns and rows so each
conductance pair has one solved value; it also ties midline mirrored inputs.
The opt-in motor fitter skips the former independent solves. Candidate option
`require_bilateral_symmetry` rejects unequal timing or saved edge conductances.
Historical configurations remain readable for reproducibility.

Sharing neural parameters does not force identical delivered torque when
geometry, muscle laws, sensor signals, or states differ. Those differences are
measured separately rather than hidden by independent controller tuning.
