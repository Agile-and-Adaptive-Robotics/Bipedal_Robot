# CHANGES manifest — 2026-10-16 SNS_Simscape winner realization

Stage: Simscape realizer (overnight walk campaign). NO commits, NO pushes.
Everything below lives uncommitted in the easteregg2 working tree
(D:\GitHub\Bipedal_Robot, branch KneeTestSetup_BenBo_stw). The laptop
authoring copies are in tmp\sns_realizer_20261016\ (session scratch,
disposable).

## Modified TRACKED files (1)

- `Code\Matlab\SNS_Simscape\README_SNS_Simscape.md` — one short section
  appended at the end ("Winner circuit slice R8, 2026-10-16") pointing at
  the new demo + reports. Revert with
  `git checkout -- Code/Matlab/SNS_Simscape/README_SNS_Simscape.md`.

## New UNTRACKED files

Code\Matlab\SNS_Simscape\demos\:
- `SNS_WinnerSlice_R8.slx` — the deliverable model (34 neurons, 72
  synapses, ode1 @ 2 ms; built from SNS_Library_R2025a blocks + the local
  masked NapNeuron for the RG somata)
- `sns_build_winner_slice.m` — generator (reads
  spinal\reports_20261016\sns_simulink\slice_export.json)
- `sns_run_winner_slice.m` — headless 4-check verification (ALL PASS,
  2026-10-16 night log in the report)
- `README_winner_slice_20261016.md` — the README note for the model

Code\MuJoCo_SNS\spinal\ (probe scripts, evidence + reusable):
- `_slice_probe_20261016.py` (winner-load + slice export + zero-afferent
  refs), `_slice_probe2_20261016.py` (scripted-contact refs),
  `_slice_probe3_20261016.py` (modulation sweep),
  `_slice_probe4_20261016.py` (i0_ii control — negative result),
  `_slice_probe5_20261016.py` (slice-topology reference ->
  slice_ref_targets.json)

Code\MuJoCo_SNS\spinal\reports_20261016\sns_simulink\:
- `slice_export.json`, `slice_ref.npz`, `slice_ref2.npz`,
  `slice_ref_targets.json`, `winner_slice_verify.png`,
  `winner_slice_verify.mat`, `SNS_SIMULINK_WINNER_SLICE_20261016.md`
  (the stage report)

## NOT touched

- No spinal-stack file (params.py / build_network.py / runner.py / jsons)
  — the winner configuration is only READ, through
  wean_tune_20261016.apply_config.
- No existing .slx (SNS_Library*, demos, imports) modified; the new model
  links against the existing `SNS_Library_R2025a.slx` copy.
- No remote pushes, no commits, per campaign contract.
