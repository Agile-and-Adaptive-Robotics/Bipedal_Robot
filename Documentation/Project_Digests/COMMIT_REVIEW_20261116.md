# Pending commit review and cleanup

Reviewed the working tree on easteregg2 at commit `a6fa2fac5f36df1f1800e9f9d477fd0a3e187778`. The initial inventory contained 2,973 pending paths, with 34 modified files, 309 deleted files, and 2,630 untracked files. Nothing was staged.

## Cleanup completed

- Scoped `.gitignore` rules keep 2,522 generated routing and schematic files local. These include ten run directories, the erector sweep checkpoints and worker files, and six superseded schematic variants. Both selected erector JSON records remain visible because the figure generator reads them directly.
- The 309 tracked files removed from old defense workspaces were archived and CRC-checked. The source commit retains them in Git history. `Defense/current_deck.json` now identifies the existing local v34a presentation. The former v21 manifest remains in its historical section, and both older pickup notes identify their historical scope.
- The overwritten `Code/MuJoCo_SNS/spinal/spinal_run.npz` was copied to the local review folder before its committed contents were restored.
- A paused-campaign status note accompanies the routing source.

The resulting working tree has 456 pending paths. Of these, 309 are the deliberate old-deck deletions, leaving 147 additions or modifications. Source and experiment outputs remain available locally. No simulation or optimization was launched during this cleanup.

## Suggested commits

| Commit | Paths | Include |
|---|---:|---|
| Keep generated routing output local | 2 | Root `.gitignore` and this review |
| Preserve paused routing implementation | 101 | The 100 changed source/configuration files under `campaign_20261103`, plus `CAMPAIGN_STATUS.md` |
| Preserve neural-controller experiments and fixes | 30 | The 24 torso/controller/timing files and six postural/serial experiment files described below |
| Update the defense route schematic | 11 | Slide57 generator, final schematic/table, descriptions, rerun script, and the two selected erector JSON records |
| Retire superseded defense workspaces | 312 | The 309 deletions and three corrected defense pointer/handoff files |

The exact file list and proposed commit for every remaining path are in [COMMIT_PLAN.csv](D:/GitHub/Bipedal_Robot/tmp/commit_review_20261116/COMMIT_PLAN.csv). The initial full inventory is [BEFORE_PENDING_FILES.csv](D:/GitHub/Bipedal_Robot/tmp/commit_review_20261116/BEFORE_PENDING_FILES.csv).

## Routing source worth preserving

The 14 modified campaign files depend on new helpers among the 86 new source files. Examples include `fr_tibialLeadMargin.m`, `fr_parallelRows.m`, and `fr_robotKneeT.m`. Preserve the full source set in its paused experimental state.

The tracked [routing audit](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/audit_20261105/ROUTING_AUDIT_20261105.md) records eight of 54 cases passing its implemented constraints and zero passing those constraints plus all signed torque targets. Later local whole-leg reports recorded one of 27 rigid groups and two of 27 compliant groups clearing the modeled checks. The split-vasti redraw recorded zero of 29 rigid groups and one of 29 compliant groups clearing. CAD knee/ankle registration errors remained approximately 87/89 mm. Cam release continuity, physical reeving, and full assembled packaging still require verification.

The excluded run directories preserve failed candidates and checkpoints. This cleanup makes no claim that their route plots constitute a qualified robot design. Completing or resuming optimization is optional future work.

## Neural source worth preserving

The 30 files form several connected experiments. The exact paths appear in the CSV.

- Torso support requires new `Code/MuJoCo_SNS/torso_support.py` and `spinal/controller_profiles_20261104/lumbar_pf.py`, together with their modified runner, evaluator, body-probe, and CAD callers.
- The modified `profile.py`, `body_probe.py`, and `replay_adapter.py` also require new `native_pf.py`. The additional timing checks, wiring reports, and knee torque diagnostics belong with those changes.
- Capture-output handling and near-zero harmonic phase reporting fixes belong with the training changes.
- `postural_pilot_20261104/robust_fit.py`, the four `serial_*.py` files, and `gait2392_stage3_probe.py` can form a separate six-file experimental commit if preferred.

Successful fitted parameters and the campaign validation records remain in ignored `tmp/TorsoPF_20261110` or machine-local training directories. The current body-probe entry point still uses `D:/Users/Ben Bolen/NeuroAcceleration_20261104/body/compact_final.json`, and fit manifests refer to `D:/Users/Ben Bolen/ControllerTraining_20261104`. A fresh clone therefore needs those inputs to reproduce the run results. The Stage 1 principal-cycle recheck also remains in the local campaign folder, while the production evaluator retains the older knee-peak assessment. Preserve these changes as experiments with known limitations.

## Slide57 files to include

All paths below are under `Code/Matlab/Mesh_Optimization/Results/KompareXXVII/Slide57`.

- `README.md`, `route_visual.py`, and `rerun_erector.m`
- `route_schematic.png`, `route_schematic.pdf`, and `route_schematic_alttext.txt`
- `actuator_site_table.png`, `actuator_site_table.pdf`, and `actuator_site_table_alttext.txt`
- `erector_rerun/erector_compliant.json` and `erector_rerun/erector_rigid.json`

These are defense schematic assets. The README describes their mixed display poses and torque/packaging limitations. Six older render variants and the sweep MAT files remain local.

## Old defense workspaces

All 309 deletions belong to revision folders v10 through v21. They include scripts, handoffs, teaching scenes, and a submitted-source snapshot. Only 52 had identical surviving tracked copies elsewhere. The archive preserves every deleted tracked file.

The cleanup changes `Defense/current_deck.json`, `Defense/ZCODE_PICKUP_20261031.md`, and `Defense/ZCODE_PICKUP_20261031_v18.md` alongside those deletions. The v34a presentation itself stays at `C:/Users/Ben Bolen/Documents/Bolen_PhD_Defense_v34a.pptx`. This pass checked its existence, size, and SHA-256. Its slide content and playback were outside this cleanup review.

Local backup locations follow.

- [Retired workspace archive](<D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Archive/defense_workspaces_retired_20261116/retired_workspaces_HEAD.zip>)
- [Archive manifest](<D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Archive/defense_workspaces_retired_20261116/MANIFEST.json>)
- [Overwritten default-run backup](D:/GitHub/Bipedal_Robot/tmp/commit_review_20261116/spinal_run_working_20261116.npz)

The ZIP is 81,197,645 bytes, with SHA-256 `e14dbd4026fcb19f75f2d9ed9b508a9fe697a093a5b64e4d934001144e9e99b7`. It covers tracked files from the stated commit. Previously ignored generated deck files depend on their original local copies.

## Verification

The review checks exact Git path counts, ignore-rule boundaries, selected JSON visibility, archive CRC and file coverage, preservation of the existing source files, and absence of staged changes. [CLEANUP_VERIFICATION.json](D:/GitHub/Bipedal_Robot/tmp/commit_review_20261116/CLEANUP_VERIFICATION.json) records those checks. Ben retains control of committing and pushing through GitHub Desktop.
