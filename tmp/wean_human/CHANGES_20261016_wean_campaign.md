# CHANGES manifest — 2026-10-16 human feedback-tuning + drive-wean campaign

Stage: feedback optimizer (human), overnight walk campaign. NO commits,
NO pushes. All changes on easteregg2 (D:\GitHub\Bipedal_Robot, branch
KneeTestSetup_BenBo_stw) + the local authoring copies on the laptop
(C:\Users\Ben\Documents\GitHub\Bipedal_Robot\tmp\wean_human\).

## Modified tracked files

- `Code/MuJoCo_SNS/spinal/spinal_run.npz` — was ALREADY M from earlier
  campaign stages before this stage began; this stage additionally
  overwrote its working-tree content via UN-redirected eval runs (the
  sens/probe subprocess evals run without an AARL_NPZ spec). The
  committed v6-winner npz is recoverable from git history. NOT restored
  (the pre-stage working-tree content is not reconstructible; it was
  itself an intermediate campaign artifact).
- `Code/MuJoCo_SNS/spinal/spinal_run.png` — churned by full runs;
  RESTORED via `git checkout -- spinal_run.png` (verified clean).

NO source files were modified by this stage. In particular
runner.py / params.py / build_network*.py are untouched (their M state
is the earlier connectome-drive stage's, per
CHANGES_20261016_connectome_drive.md). No new params.G knobs were
introduced; all searched keys already had loaders (set_stage(4)/
set_params) or are the four hand-drive knobs loaded by the documented
post-set_stage recipe. The standing default-build contract is therefore
untouched, and the anchor-A baseline was re-verified bit-exact
(-161.56754173676563) twice this stage (before and after the driver
fix).

## New untracked files (spinal/)

wean_tune_20261016.py (campaign driver), wean_round_launch_20261016.py,
wean_probe_20261016.py, wean_full_20261016.py, wean_score_20261016.py,
wean_sens_20261016.py, wean_poll_20261016.py,
basin_gate_drive_20261016.py, basin_closed_20261016.py,
wean_{gate,r1,r2,r3,r4,r5,r6,r7,r8,probe,probe2,sens}.cmd (schtasks
launchers).

## New untracked files (reports_20261016/wean_human/)

WEAN_CAMPAIGN_20261016.md (the stage report), round{1,2,4,5,6,7,8}_w*.
jsonl (per-trial full kine diagnostics), round{1,5,6,7,8}_best*.json
(winner param files; rank suffixes = n-th best), probe_round1.jsonl +
probe.jsonl (drive-scale ladders), sens.jsonl (gain-sensitivity),
basin_round8_best3.json.json (network-only gate: inapplicable, baseline
decays) + basin_closed_round8_best3.json.json (closed-loop gate 8/8),
round{1,6,8}_full*.npz (production-schedule runs), finalist_spec.json +
finalist_eval.npz, task/round logs, CHANGES_20261016_wean_campaign.md
(this file).

## Optuna studies added to spinal/optuna_walk.db

curr_s6b_drive_human (100), curr_s7_wean07 (39, FB-inert — invalid as
a feedback re-tune, kept as record), curr_s6c_prodwindow (6, aborted),
curr_s6d_feedback (88), curr_s7b_wean07 (80), curr_s6e_prodwindow (64),
curr_s8b_wean04 (48).

## Scheduled tasks

All aarl_wean_* schtasks deleted after use (gate, r1..r8, probe, probe2,
sens). Verified none remain.

## Winner (for the supervisor gate)

`reports_20261016/wean_human/round8_best3.json` (curr_s8b_wean04
trial 35, drive_scale 0.3807379786697019). Reproduce:
`ssh -o BatchMode=yes easteregg2 "cd /d D:\GitHub\Bipedal_Robot\Code\
MuJoCo_SNS\spinal && set PYTHONIOENCODING=utf-8&& D:\Anaconda\envs\myo\
python.exe wean_tune_20261016.py one reports_20261016\wean_human\
finalist_spec.json"` — observed this session: kine_score
-145.07280220616852 (bit-exact vs the study value), 10+10 cycles, duty
0.71/0.65, no NaN.
