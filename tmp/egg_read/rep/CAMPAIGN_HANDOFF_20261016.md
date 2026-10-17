# Overnight walk campaign — handoff (written 2026-10-16)

> **SHUTDOWN ADDENDUM (2026-10-16 ~07:45 UTC / 12:45 AM PDT):** the orchestrating laptop
> was deliberately powered down mid-campaign (Ben moving to the lab). The live run is
> `dwfrun-5685cc77-2ff6-426f-a917-7ff04d8a964f` — it will show as stopped/interrupted;
> settled phases (recon + whatever audits finished) replay from its journal on resume, only
> in-flight asks re-run. To resume: open ZCode in `C:\Users\Ben\Documents\GitHub\Bipedal_Robot`
> on the laptop and resume that run id (a fresh session can find it via the runs list).
> No data loss beyond in-flight agent turns; nohup'd runs here finished on their own.

This file is the self-contained continuation guide for the campaign Ben launched tonight.
Read this if the orchestrating laptop (DESKTOP-5Q16KE9) goes down, or in the morning to
pick up wherever the campaign stopped.

## What is running, and from where

A ZCode dynamic workflow (`dwfrun-5685cc77-2ff6-426f-a917-7ff04d8a964f`, script on the
laptop at `.zcode\workflow-drafts\Walk-campaign-LR-drive-bootstrap--connectome-audit.dwf.ts`)
is orchestrating from the laptop over ssh. The laptop only orchestrates — every heavy run
is nohup'd ON easteregg2, so remote runs keep going even if the laptop disconnects.
Laptop is on AC with sleep disabled on AC (verified tonight), so a Windows lock (Win+L)
does NOT stop it. Only battery power (sleeps after 5 min) or a crash stops the laptop side.

Campaign outputs land in THIS folder (`spinal\reports_20261016\`) as stages complete.

## The plan (11 phases, in order)

1. Recon easteregg2 (+ AnimatLab reachability, egg→EB ssh test).
2. Connectome audit vs Ben's rules, both models.
3. Implement L/R alpha+gamma drive, default-off, gated bit-exact.
4–6. Human track: drive bootstrap sweep → (if not walking: connectome repair loop, max 3
   rounds) → optuna feedback tuning with drive-wean ladder → SNS-Simscape realization.
7–9. Robot track (W2L robot, w2l_mujoco, harness/virtual-walker protocol): same arc,
   realization in AnimatLab (headless AnimatSimulator on EB475WS4 if reachable).
10. Supervisor gate: re-runs every claimed winner, checks gates/sentinels/basin results.
11. Deliver the overnight report (also published as a card on the laptop).

## Ben's strategy encoded in this campaign (his words tonight)

- Separate LEFT and RIGHT hand drive, each split into alpha (tonic drive reaching MN
  pools / PF layer) and gamma (drive reaching the RG half-centers) — four default-0 knobs
  per model.
- Bootstrap walking with strong (possibly asymmetric) drive first.
- Then wean the drive down (ladder ~1.0 → ≤0.7 → ≤0.4) while re-optimizing feedback
  gains; report how far feedback alone sustains the gait.
- **If a model does not walk even under strong drive, treat the connectome as the
  suspect** — the workflow automatically loops into connectome repair (up to 3 rounds
  per model) before giving up on that model.
- Target: clearly better than the of-record −189.5 walking kine, ≥6 real cycles,
  duty 0.4–0.75, no NaN, stays up. No sentinel/frozen-model scores.

## Contracts that MUST hold for anyone continuing this work

- Default-off = bit-identical: default build counts (410, 376, 1186) unchanged; a recorded
  baseline eval reproduces exactly with new knobs at 0.
- Conditional topology: new neurons are built only when their gain > 0.
- JSON RULE: every params.G knob a study uses is saved in the winner json AND has an
  `if key in best` loader branch in the runner.
- NO git commits, NO pushes. All changes stay in the working tree + this reports folder.
- Frozen/NaN runs never count as walking; winners sitting at gate minimums get flagged.
- Finalist configs pass the basin gate (±1% param perturbation, rhythm persists).

## Environment facts (easteregg2 unless noted)

- Repo: `D:\GitHub\Bipedal_Robot`; spinal stack: `Code\MuJoCo_SNS\spinal\` (read DESIGN.md
  before touching runner.py / build_network.py / params.py).
- Python: `D:\Anaconda\envs\myo\python.exe` (py3.10.21, mujoco 2.3.7, sns-toolbox 1.5.2,
  optuna 5.0.0). MATLAB R2025a (`matlab -batch` works headless).
- Optuna dbs: `optuna_walk.db`, `optuna_w2lvar.db`, `optuna_syn6.db` (schema 12).
  New studies tonight use names starting `curr_s6_`.
- Gait references: `spinal\gait_refs\*.npz` (43 refs: Ong walking 050..200 + selfsel,
  Falisse, Arnold running).
- AnimatLab headless: EB475WS4 only —
  `D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe`.
  ssh aliases `eb475ws4` / `easteregg2` (from the laptop; laptop→EB verified tonight).
- W2L robot: `spinal\w2l_mujoco\` (M1 MJCF `w2l_mjcf_fixed.xml`) + `spinal\w2l_cpg\`
  (91-neuron contact-driven bilateral CPG; W2L_SMOKE gate = PASS baseline, period 1.000 s).

## How to check campaign state right now

- `dir Code\MuJoCo_SNS\spinal\reports_20261016\` — per-stage reports + changes manifest.
- `git status` in the repo — campaign edits, uncommitted by design.
- Optuna progress: load `curr_s6_*` studies from the dbs and list best trials.

## How to continue by hand (if the laptop side is dead)

- Reproduce winners: each stage report here records exact reproduce commands; the
  supervisor report re-ran them.
- Continue optimization: resume the `curr_s6_*` studies (optuna load-if-exists) with the
  same objective machinery (`optuna_walk.py` / `_curriculum*.py` family; robot side per
  the w2l drivers).
- The L/R drive knobs: search params.py / build_network.py / runner.py for the knob names
  recorded in the drive-implementation stage report (planned: alpha_l/gamma_l/alpha_r/
  gamma_r semantics — the implementer's report is authoritative).
- Workflow resume (laptop only): in a ZCode session, `ResumeWorkflowRun` on
  `dwfrun-5685cc77-2ff6-426f-a917-7ff04d8a964f` — finished phases replay free.

## Survival matrix (laptop dies mid-campaign)

- SURVIVES: nohup'd runs on easteregg2 (they run to completion), optuna dbs, everything
  already written here, all working-tree edits.
- STOPS: phases not yet started (later stages, realizations, supervisor gate), and the
  final summary card (which exists only on the laptop).
