# Gait2392 output coordination — easteregg2

Session 1, 31 October 2026. Repository: `D:\GitHub\Bipedal_Robot`.

## Decision for ZCode

1. **No bulk regeneration of the six jobs.** Reuse surviving evidence. Regenerate only an explicitly needed artifact for the next diagnostic or a selected figure. Do not repeat failed insensitive sweeps to reconstruct a folder.
2. **Session 1 coordinates housekeeping.** During active modeling or submission, neither session prunes shared `tmp`, campaign `runs`, or another session's outputs. If Session 1 is unavailable, cleanup waits. This assignment is coordination, not a request to perform cleanup now.

All new ZCode campaign output goes under:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\zcode_job_outputs_20261031\`

Use a distinct subdirectory for each job and attempt; preserve existing results. Put reports, configs, plots, logs, source snapshots, and raw arrays together. Link their full paths from `TRAINING_RECORD.md`, outside the ignored run tree. Shared `tmp` is unsuitable for the only copy of evidence cited in a completion report.

Before any later cleanup, inventory proposed targets, identify their owner, and verify that every cited result has a retained copy. Do not delete active outputs or another session's files. No cleanup was performed to establish this rule.

## Verified surviving evidence

An inventory with full paths, sizes, and SHA-256 hashes is saved at:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\zcode_job_outputs_20261031\SURVIVING_EVIDENCE_MANIFEST.json`.

Verified by reading the actual files:

- Original tonic-abductor checkpoint: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\tonic_abduction_20261030_s1_v1\0p02\CHECKPOINT.json`.
- Failed-fit campaign: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\zcode_hunt_fit_20261030_v1\`. Its training record explicitly says `physical_curriculum_passed=false`; the air-pair assessment fails. Its independently fitted timing gains include the old large left/right hip-extensor asymmetry. Preserve this as a historical failed fit; do not promote it or treat it as the current tied-fit controller.
- Job D preview: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\native_joint_pf_mono_20261030_v1\air_retry\preview\stage2_air_neural.gif`.
- Later corrected-pose human-model checkpoint: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\aligned_support_ladder_v1\CHECKPOINT.json`. This records 5 right/4 left qualifying swings, 63.73% BW mean assistance, protected knees, and exact fresh-repeat arrays. Same neural configuration passes air stages 1/2; standing still fails. Its eight listed current source hashes were checked and match.

## Missing artifacts and reconstruction limits

The six named `D:\GitHub\Bipedal_Robot\tmp\zcode_jobs\` job folders contain no files in this checkout. ZCode reports that 232 per-seed outputs, force-chain figures, evaluator hashes, and staged ground-ladder configurations were lost. These loss counts are attributed to ZCode's report, not independently reconstructed here. Several conclusions survive in `TRAINING_RECORD.md`; that summary does not replace the missing raw sensitivity outputs.

The available evidence does not identify who pruned those files. Do not assign blame from empty directories alone. This pass performs no deletion.

Regenerated runs must be labeled **new executions**, with the current source/config hashes, seeds, settings, and timings. They do not recover the lost historical evaluator identity or establish that an older basin test used today's code. Current source hashes can be captured without rerunning all six campaigns.

## Next useful modeling work

Reconcile surviving force-chain and fit evidence with the corrected-pose checkpoint. Investigate delivered MN activation versus the standing initializer's targets and phase-specific stance support before reducing assistance further. Earlier 43.61%-assistance stepping used a deeply penetrating initial pose and is superseded; lower-support aligned trials have zero qualifying swings. Do not reconstruct obsolete ground-ladder settings and advance them as if they were the accepted checkpoint.

Working dissertation sources stay unchanged; the dissertation editor owns integration. The separate training record remains the modeler's primary experiment log. No simulation campaign was restarted by this coordination note.
