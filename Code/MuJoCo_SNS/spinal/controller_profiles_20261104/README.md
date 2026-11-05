# Four-command Gait2392 controller experiments

**Naming clarification:** the internal selector `hunt` identifies Ben Bolen's
**BiW2Lvar+ controller**, using Hunt's four-step calibration method and
Bolen's five-stage curriculum. Plain identifier: `BiW2LvarX`.
It is not a published controller supplied by Hunt. Its actual inheritance is
W2LVarnet → FunctionalNetwork → JointIntegrationNetwork → HuntNetwork. The
expanded candidate has joint/axis PF channels and a four-cell RF/hamstring
module; the separate `w2lvar` selector retains the simpler architecture.
Keep saved selector names for reproducibility, but use the descriptive name in
presentation captions. Stage 2 is afferented stepping in air. Its measured
suspension force is not a minimum walking-support requirement.

Machine: **easteregg2**, repository `D:\GitHub\Bipedal_Robot`.
The separate [TRAINING_RECORD.md](TRAINING_RECORD.md) is the execution record.
Raw configurations, studies, traces, and compiled models are under
`D:\Users\Ben Bolen\ControllerTraining_20261104`.

`profile.build()` composes opt-in graded BiW2Lvar+, s3k, w2lvar, and syn6 candidates
before compilation. Original builders and dissertation checkpoints remain intact.
Each candidate has explicit F/E commands for each side, mirrored connection
parameters, and declared Ia/Ib/II projections to central layers. Contact checks
flag any explicitly inhibitory sensory entry, heel fanout to both antagonistic
interneurons, and toe-to-RG entries. Intrinsic reciprocal inhibition is separate.
These new results do not validate a spiking implementation.

## Current results

All four families passed intrinsic-rhythm A and expected-afferent timing B.
Ordered per-muscle peak fitting C uses one bounded scale shared across sides;
reflex-waveform D freezes those scales while fitting local Ia/Ib gains. BiW2Lvar+'s
replay also applies its live pose-dependent recruitment gates. BiW2Lvar+ and syn6
passed the offline C/D criteria. s3k and w2lvar passed C but failed D.

BiW2Lvar+ completed protected physical stages 1 and 2 using the same learned weights
and 0.1 PF expression: 2.466 s and 1.300 s principal periods, respectively.
Stage 2 is suspended air motion with 99.860% BW spring support and zero foot
loading. Syn6 has a slower repeated air candidate, but target cadence and
ankle-waveform quality still need correction. Full-body postural integration,
standing, and ground stepping are subsequent gates for these revised candidates.
The isolated [postural pilot](../postural_pilot_20261104/README.md) is tested separately.

Portable figures, GIF, exact paths, and the concise defense update are in:

`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104\EXPERIMENT_UPDATE.md`

## Reusable workflow

- prepare.py: construct profiles and export wiring/source checks, without stepping.
- train.py, campaign.py: bounded Rustuna A/B timing studies through the persistent
  Python 3.10 scientific worker and isolated Python 3.11 optimizer.
- refit_motor_peaks.py, refine_reflex.py, finish_outputs.py: ordered C/D fitting,
  actual-network evaluation, manifests, and fresh-process checks.
- body_probe.py, air_validation.py: integrate frozen candidates using the existing
  physical evaluator, compact backend, measured assistance, and predictive knee
  protection. Higher-expression guard failures are retained.
- dominant_cycles.py, trace_spectrum.py: principal-cycle and spectral diagnostics.
  A knee-cycle pass alone is insufficient to exclude fast ankle modes.
- stimulus_trials.py, stimulus_figures.py: recorded external-current pulses
  separately at RG, PF, and MN levels, paired with an identical pre-pulse baseline.
- noise_trials.py: fixed-parameter, seeded colored RG-current perturbations with
  saved actual noise and a fresh-process reproducibility check.

All numerical outputs go to an explicitly supplied external directory.
The original preparation record remains at
`D:\Users\Ben Bolen\ControllerPreparation_20261104\EXPERIMENT_RECORD.md`.
No manuscript, CAD routing, tendon law, or damping settings are changed here.

Canonical display name: **BiW2Lvar+**; plain identifier: **BiW2LvarX**.
Saved `hunt` strings and path names are deprecated compatibility identifiers,
not a model attribution. Hunt's four-step calibration and Bolen's five-stage
curriculum are distinct parts of the workflow; the revised candidate has so far
passed physical curriculum stages 1 and 2, not all five stages.
