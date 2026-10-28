# Quantitative verification

- Scope: only this new result folder is added under SNS_Simscape. No commits or pushes. Original-file SHA256 comparison passes for all eight preserved source, model, library, and result files (`original_integrity.csv`).
- Loaded model versus builder: exact functional dialog-parameter and input-edge equality (`builder_comparison.txt`, comparison flag 1).
- Existing outputs: original Simulink angle 14.8969 to 44.8709°; existing sequential Euler twin 14.8974 to 42.9548°, both in original positive-flexion coordinates. These are startup records, not transient recovery demonstrations.
- Original settling failure: 25–30 s angle range 0.7116°, velocity range 0.6973 rad/s, both MN outputs span 0–1. Original pulse trace retained as a failed operating-point diagnostic.
- Graded configuration passes explicit angle, velocity, six neural-state, and two activation-state settling gates. Baseline 25–30 s angle range 2.529e-10 rad, maximum velocity bounded well below 1e-5 rad/s, extensor activation range 1.711e-12. Full settling trace is provided in settling.pdf/.fig/.png.
- Protocol: 2 nA, 30 ms, at t=30 s, summed into SN_Ia_ext Iapp after the spindle transducer. No direct activation or mechanical injection.
- Final paper-coordinate baseline: -42.6392°. Extension excursion +1.0908°; onset above 0.01° at 19 ms; peak at 111 ms. Timing resolution is 1 ms.
- Recovery: 6.355 s after pulse onset under simultaneous 0.01° displacement and 1e-4 rad/s velocity thresholds. Late mean residual +9.3011e-7° over 45–50 s. Control angle range after 25 s: 1.4694e-8°.
- Numerical sensitivity: maximum full-trace angle difference 5.2027e-5° for 1 versus 0.5 ms maximum solver steps, unchanged tolerances. Fresh-process replay: resampled full pulse arrays exactly equal, maximum angle difference 0 rad (`verification.txt`).
- Raw MAT outputs preserve startup, all six membrane states, MN normalized outputs, stimulus, activations, forces, angle, velocity, and extensor Ia synaptic pairs. Final parameters are in replay_parameters.mat. The plot uses only a sign conversion for paper coordinates and no smoothing.
- Remaining limits: simple force law, arbitrary sensory-current amplitude, broadened saturation ramps, flexor inactive at baseline, no explicit tendon-impact mechanics, no physiological validation, and no separate pathway-disabled experiment. See README for equations and the proposed caption.
