# Isolated antagonistic postural test

Test before deployment into Hunt, s3k, w2lvar, or syn6. Machine: easteregg2.
External records: `D:\Users\Ben Bolen\PosturalPilot_20261104`.

- `initial_k80_d10`: six direct/premotor × Ib-policy tests, neural-backend and DC
  gain checks, saved mechanical and neural samples, exact model XML.
- `mirrored_neural_fit_v1`: 24 Rustuna trials fitting four shared neural scales.
  Both antagonists use identical parameters. Training uses mirrored ±0.04 rad
  starts and ±8 N m, 50 ms pushes.
- `held_out_v1`: frozen baseline and fitted candidates with ±0.08 rad starts,
  ±16 N m pushes, and seeded colored neural-current noise of 0, 0.025, or 0.05 nA.
  These are larger/different conditions than training; no parameter jitter is used.
- `delay_noise_v2`: 72 frozen-parameter tests, adding 0/30/70 ms measurement
  delay with multiple unseen noise seeds. The quiet-condition fit fails delayed
  tests despite its nominal success.
- `robust_delay_fit_v3`: 24 trials that include delay/noise during training;
  6/6 training cases passed, but only 9/18 held-out cases. Not accepted.
- `robust_delay_fit_v4`: follow-up with 0.05 nA training noise and shared
  premotor/relay time constants: 6/6 training and 17/18 unseen cases, with an exact
  fresh-process training repeat. The failed case settled at 3.125 s against a
  predeclared 3 s gate, so the candidate remains unaccepted. Do not infer full
  acceptance from BEST.json's training score.
- `figures`: short-window premotor, MN, tension, and joint traces, plus a saved-state
  replay GIF. This is an ankle bench, not a walking GIF.

The circuit contains 14 graded cells. Commands enter PP_POS/NEG; PP drives MN and
its inhibitory relay; each relay inhibits the opposite premotor and MN cells.
Ia/II excite the premotor channel of the corresponding tensile actuator. Ib
passes through a relay with an explicitly selected inhibitory or excitatory
connection. There are no contact inputs or RG/PF layers in this isolated bench.
POS/NEG refer to actual torque signs measured from MuJoCo tendon moment arms.

The comparison matches the isolated DC command-to-MN transfer. The task command
uses the same PD gains in both architectures, kp=80 N m/rad and
kd=10 N m s/rad, with torque normalized by the measured actuator moment arm and
300 N capacity. FSA gain is g=R Gm/(Erev−R), with R=5 mV, Gm=1 uS, and
Erev=8 mV. The extra premotor lag is retained; matching DC gain does not match
closed-loop dynamics. Local Ia→MN input remains in both alternatives.

The mechanism is a fixed-base, 10 kg single-ankle inverted pendulum with two
filtered tensile actuators. Neural/physics steps are 1 ms. No harness is applied;
the base is fixed. Passing requires a finite full 8 s trial, final RMS angle
below 1 degree, and settling within 3 s after the push. The held-out summary
reports successes and failures without extending those claims to Gait2392.

Scripts: model.py (circuit and mechanism), run.py (initial comparison), fit.py
(Python 3.11 Rustuna optimizer using a persistent myo scientific worker), verify.py
(fresh-process/held-out/reference checks), and figures.py (saved-state plots and
animation). Use full Python paths in AGENTS.md. No package/environment changes.

## Extension based on Ben's intern experiment

Ben described a physical leg-like inverted pendulum balanced on a ball-and-socket
joint with three BPAs, acting as soleus, tibialis anterior, and extensor digitorum
longus analogues. Interns used bang-bang control with IMU orientation/inertial
and angular-rate feedback. This is Ben's reported physical baseline; its geometry,
pressure settings, and measured performance have not been imported here.

The next comparison should reproduce that control problem before scaling up:

1. A ball-and-socket ankle pendulum with three tension-only BPA routes. Compute
   the moment-arm matrix and available torque directions before fitting controls.
   A ball joint permits three rotations; three actuator routes do not automatically
   provide independent control of all three, particularly axial yaw.
2. Compare bang-bang and antagonistic neural-premotor control on identical
   geometry, IMU signals, pressure bounds, and actuator dynamics. Keep the task
   error estimator and any actuator allocation explicit. Start with the existing
   pressure/length characterization; label assumed pressure dynamics separately.
3. Extend to linked double and triple inverted pendula. Distinguish extra movable
   joints from extra actuated joints; identify any underactuation. Evaluate segment
   angles/rates and whole-body COM, rather than only the top segment's orientation.

Use x anterior, y superior, z right. Matched tests should include startup tilt,
fore-aft/lateral pushes, IMU noise/bias, sensing/valve delay, pressure limits,
and held-out conditions. Report angle error, settling, falls, pressure activity,
and saturation. Do not infer neural superiority from different actuators or
different assistance. Preserve the existing 10 kg hinge test as the simpler
reference. These extensions are the next benchmark design, not completed results;
the 17/18 result continues to refer only to the original hinge bench.
