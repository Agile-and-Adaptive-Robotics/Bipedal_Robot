# MONITOR.md — run_0_weights_placeholder (launched 2026-11-04, easteregg2)

Run folder: `campaign_20261103/run_0/` (this folder).
Driver: `campaign_20261103/fr_launch_run0.m`, ONE detached MATLAB process
(`-batch`), `parpool(10)`, sequential contract layers, parfor groups within
each layer, per group rigid → compliant (Xi0–Xi3), refinement round 1 after
the deep-outward pass. Decision record: `LAUNCH_DECISION.md` (full-27
branch; weights PLACEHOLDER-PENDING-BEN).

## How to check progress (PowerShell or Git Bash, from the repo root)

- Master log (layer starts, per-group one-liners, refine summary):
  `campaign_20261103/run_0/logs/master.log` — `tail -f` it, or
  `Get-Content ... -Tail 30`.
- Checkpoints (one .mat per group per branch, written as each completes):
  `dir campaign_20261103\run_0\checkpoints` — expect up to
  27 `<group>_rigid.mat` + 27 `<group>_compliant.mat` + 27
  `<group>_refine1.mat`. Layer completion = all groups of the contract
  layer present (trunk/deep hip 5, thigh 10, lower leg 12).
- Live status JSON: `campaign_20261103/run_0/STATUS.json`
  (rewritten at layer/round boundaries and at the end).
- Per-group detail (solver trajectories, failures): `logs/group_<name>.log`.

## Expected completion

Launched 06:37:04 local 2026-11-04 (the running process; two earlier
launch attempts crashed on driver start-up bugs before any solve —
jsondecode column-vector layer ids, then a pool logging property — both
fixed; their crash lines remain at the top of master.log). Layer 1
(trunk/deep hip) completed ~06:45 with all five groups solved on both
branches. Expected total wall **≈ 2–4 h** (measured planning numbers:
seed-eval mean 0.354 s → 3.98 h serial for all 40,500 evals; parfor(10)
reduces wall below the serial bound; realistic worst bound ≈ 8 h). If
`STATUS.json` shows `stoppedCleanly: true`, the run was stopped via the
STOP procedure below and resumes from checkpoints.

## How to STOP cleanly

A stop must leave checkpoints intact and resumable.

1. **Graceful (preferred):** create the sentinel file
   `campaign_20261103/run_0/STOP` (empty is fine). The driver checks it at
   every layer boundary and before the refinement round; it finishes the
   in-flight layer (all its checkpoints land), writes `STATUS.json` with
   `stoppedCleanly: true`, and exits. Latency ≤ one layer (~15–60 min).
2. **Hard stop (if you cannot wait):** kill the MATLAB process
   (`taskkill /IM MATLAB.exe /F` — careful: kills ALL MATLAB sessions; or
   kill the specific PID from Task Manager, the `-batch` one whose command
   line contains `fr_launch_run0`). In-flight group solves are lost
   (no partial checkpoints are ever written); every completed branch's
   checkpoint survives. Do NOT delete anything in `checkpoints/`.

## How to resume / relaunch

- First DELETE `campaign_20261103/run_0/STOP` if a graceful stop wrote it
  (the driver refuses to start while STOP exists, writing STATUS and
  exiting without work).
- Relaunch the SAME detached command (safe any number of times):

```
"D:/Program Files/MATLAB/R2025a/bin/matlab.exe" -batch "cd('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027'); addpath('campaign_20261103'); fr_launch_run0"
```

The driver re-scans `checkpoints/` on startup (resume counts printed in
`master.log`), skips every completed (group, branch) pair and a completed
refinement, restores the envelope-of-record from the checkpoints, and
continues the layer/round sequence. The manifest is written once and kept
on resume. Psoas already carries full-budget rigid+compliant checkpoints
from the launch-plumbing smoke; the campaign skips them.

## After completion

`STATUS.json` carries per-group rigid/compliant f, violations, and the
refine verdict (accepted/rejected). Per-branch results for Ben's
weight-ruled re-solve live in `checkpoints/` (design, rep, info, envelope
per branch) — the refinement from checkpoints re-uses these, no rebuild.
