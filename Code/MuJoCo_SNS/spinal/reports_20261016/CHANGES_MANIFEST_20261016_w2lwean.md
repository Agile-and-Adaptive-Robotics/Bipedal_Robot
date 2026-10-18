# CHANGES MANIFEST — feedback optimizer (robot) stage, 2026-10-16

Campaign stage 8 ("feedback optimizer (robot)"). Working-tree edits on
easteregg2 (`D:\GitHub\Bipedal_Robot`), NO commits, NO pushes. All other
modified files in `git status` belong to earlier stages (audit:
`reports_20261016\w2l_audit\CHANGES_MANIFEST_20261016_w2l.md`; bootstrap:
`reports_20261016\w2l_walk_sweep\CHANGES_MANIFEST_20261016_w2lwalk.md`;
human-track + drive-implementation stages have their own).

## Code\MuJoCo_SNS\spinal\w2l_mujoco\test_w2l_ground.py

5. **KV knobs `g_heel=0.5, g_toe_in=1.0, g_toe_df=5.0, g_toe_df_mn=2.749,
   g_ib=0.5, comm=1.0`** (feedback-gain passthrough for the wean study)
   + the `run_phase` build call now passes them:
   `B.build(comm=KV["comm"], alpha_l=…, alpha_r=…, gains=dict(aff_heel=…,
   aff_toe_in=…, aff_toe_df=…, aff_toe_df_mn=…, aff_ib=…))`.
   DEFAULTS equal `build_w2l_aff_net.AFF_GAINS` + comm 1.0 exactly, so
   the default build is unchanged BY CONSTRUCTION; gates re-run post-edit
   (both bit-exact vs pre-edit, this session):
   - default walk gate (`test_w2l_ground.py --phase=walk`): pelvis z min
     0.951, tilt 8.2, harness 49%, strikes 0/0, fwd +0.001, RG-E bursts
     11/10 period 1.539 s, r −0.776 — identical to the recorded baseline
     (bootstrap manifest lines 26–30).
   - bootstrap winner command (`--te=4 --tf=5 --legcoll=0 --sink=0.06
     --ctn_amp=2 --rig=0.25 --walker-yield=1.5 --dur-walk=30`): L33/R41
     strikes, +1.676 m (+0.066), 22/21 bursts, r −0.743, tilt 16.1,
     pelvis z min 0.905, harness 26% — identical to the pre-edit re-run
     and to the recorded w01 winner.
   No other source file touched.

## New files (all under reports_20261016\w2l_wean\ unless noted)

- `_w2l_wean_eval.py` — in-process evaluator (metrics = the bootstrap
  stage's definitions; appends the 402-row `wean_trials.csv`
  anti-lottery diagnostics log; CLI prints `RESULT {json}`)
- `_w2l_wean_campaign.py` — optuna ladder driver (TPE, studies
  `curr_s6_drive_robot_r100…r020` in `optuna_w2l_robot.db`)
- `_w2l_wean_deep.py`, `_w2l_wean_mid.py` — ladder extension 0.15–0.0 +
  0.175 midpoint (studies r015/r010/r005→not-run/r000→not-run, r018)
- `_w2l_wean_basin.py` — basin gate (12 seeded joint ±1% draws) +
  gain-sensitivity probe (`basin_results_r{100,060,020}.json`)
- `_w2l_cli_verify_rungs.py` — fresh-process CLI audit of every rung
  (`cli_verify_rungs.json`) — the audit that caught the in-process
  context artifact (§4 of the stage report)
- `_kv_diff_probe.py` — proves evaluator-vs-testCLI configs identical
  (only the unused `phase` KV key differs)
- `_w2l_finalize.py`, `_w2l_write_wean_params.py` — param-file writers
  (`w2l_wean_best_params.json`)
- `w2l_wean_run.cmd`, `campaign_launch.cmd`, `w2l_deep_run.cmd`,
  `campaign.log`, `deep.log` — WMI-detach launchers + logs (egg detach
  recipe: nohup/start/Start-Process all die with ssh; WMI Create works)
- `round_best_r*.json` (12), `ladder_summary.json` (with
  `cli_robust_view`), `wean_trials.csv`
- `..\W2L_ROBOT_WEANLADDER_20261016.md` — the stage report
- this manifest

Stray-file cleanup: `C:\Users\Ben Bolen\campaign.log`/`campaign.err`
(0-byte Start-Process redirect artifacts) deleted this session.
