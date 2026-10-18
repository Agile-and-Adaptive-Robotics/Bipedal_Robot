# HUMAN DRIVE SWEEP — phase-4 bootstrap (2026-10-16, overnight walk campaign)

Stage role (the ask): on the 92-muscle HUMAN body (spinal runner, s3k
winner config as the base), sweep the new L/R alpha+gamma hand-drive
knobs over a coarse grid — gamma in {1.5x, 2.5x, 4x} of the tuned drive
level, alpha in {0, 0.3, 1.0}, L/R imbalance in {0, +-20%} — enough
combos to answer: **does STRONG asymmetric drive make it walk (>= 6 real
cycles, duty 0.4-0.75, stays up 20 s)?** Report the best walking combo
and its param file, or exactly why nothing walked.

## 0. Headline answer

- **Strong ASYMMETRIC drive does NOT make it walk.** Every +-20%
  imbalance row either froze the left leg entirely (all dm20 rows:
  cycles_l = 0), left it permanently planted (dp20 rows at alpha 1.0:
  duty_l = 0.99 = the left foot "cycles" without ever really unloading),
  or froze BOTH legs (g1.5_a0.3_dp20, the -25.0 sentinel). Asymmetry
  never helped; it only removed the left leg.
- **Strong SYMMETRIC ALPHA does.** alpha_l = alpha_r = 1.0 (u[ALPHA_c_s]
  = drive*1.0 ~ 2.08 nA into each ALPHA cell, fanning to every
  ipsilateral MN at g = 0.5) recruits the frozen left leg at every gamma
  level; gamma alone (1.5x-4x, any imbalance) NEVER does (0 left-leg
  cycles in all 18 alpha<=0.3 rows).
- **Best walking combo: `g1.5_a1_d0` = gamma_l=gamma_r=1.5,
  alpha_l=alpha_r=1.0** (imbalance 0) on the anchor-A base at the s3k
  winner drive 2.0773198435110425.
  Param file: `reports_20261016/drive_sweep_human/winner_g1.5_a1_d0.json`.
  Honest gait description: bilateral in-place alternation — real
  contact-duty cycles on BOTH legs, no NaN, no fall, but no forward
  progress (dx -0.02 m) and hip/knee amplitudes far under reference.
  It is a bootstrap gait, not a tuned walker.

## 1. Method (exactly what ran)

- Base config = campaign **anchor A** (bit-exact baseline): v10
  multipliers + renshaw 0.5 + s3k trial-34 winner params
  (`reports_20260923/s3k_trial34_full_params.json`) + `full_rules=0`
  override, `set_stage(4, base)`, `runner --eval --drive 2.0773...`
  — the recipe in `spinal/_s3k_anchor_eval_20261016.py` flavor A.
- Knobs set exactly per `DRIVE_VERIFY_20261016.md` section 6: after
  `set_stage`, `P.G[k] = value` for gamma_l/gamma_r/alpha_l/alpha_r,
  before the build inside `R.main`. Port currents (runner.py step loop):
  `u[GAMMA_c_s] = drive*G[gamma_s]`, `u[ALPHA_c_s] = drive*G[alpha_s]`,
  so knob g means g x 2.0773 nA at the port.
- Grid: 3 gamma {1.5, 2.5, 4.0} x 3 alpha {0.0, 0.3, 1.0} x 3 imbalance
  {0, +0.2, -0.2} + a knobs-0 baseline row = **28 combos**. Imbalance
  delta scales the LEFT knobs x(1+delta) and RIGHT knobs x(1-delta)
  (whole-side hand-drive asymmetry; at alpha=0 rows only gamma is
  affected).
- Each combo ran in a FRESH subprocess (params.G is module-global), cwd
  `Code\MuJoCo_SNS\spinal`, `PYTHONIOENCODING=utf-8`.
  Scripts: `spinal/drive_sweep_human_20261016.py` (orchestrator,
  resume-safe via results.jsonl) -> `spinal/drive_sweep_one_20261016.py`
  (per-combo eval, appends JSON to
  `reports_20261016/drive_sweep_human/results.jsonl`).
- Launch pattern (egg has no nohup): plain foreground ssh died to a
  mid-run "Connection reset by peer" (the documented egg flakiness), and
  a PowerShell `Start-Process` child died with the session. The sweep
  ran to completion via **schtasks** (`aarl_drive_sweep`, deleted after):
  `schtasks /create /tn aarl_drive_sweep /tr "cmd /c
  ...\drive_sweep_launch.cmd" /sc once /st 23:59 /f && schtasks /run`.
  Log `drive_sweep_human/sweep.log`: "ORCH DONE", EXITCODE 0,
  28/28 combos, zero failed combos.

## 2. In-harness gate (this stage's own runs)

- Baseline row (all knobs 0) reproduced the anchor A eval **bit-exact**:
  kine_score = **-161.56754173676563** (expected
  -161.56754173676563). PASS.
- The anchor-A base itself is a ONE-LEGGED walker: 9 right-leg cycles
  (contact duty 0.47, T 1.15 s), left leg FROZEN (0 cycles, left foot
  planted 98% of the walk window) — the situation the bootstrap targets.

## 3. Results (16 s eval window, walk 2-13 s; full table in
`drive_sweep_human/summary.csv`, raw per-combo lines in `results.jsonl`)

Reading guide: cyc = real contact-onset cycles per leg; duty = per-leg
contact duty (reference duty_r/l ~ 0.61); kine = kine_ref score
(higher/less-negative = better; -25.0 = frozen sentinel).

| id | gamma l/r | alpha l/r | kine | cyc r/l | duty r/l | fell |
|---|---|---|---|---|---|---|
| base | 0/0 | 0/0 | -161.57 | 9/0 | 0.47/- | no |
| **g1.5_a1_d0** | **1.5/1.5** | **1.0/1.0** | **-173.10** | **6/6** | **0.67/0.60** | **no** |
| g2.5_a1_d0 | 2.5/2.5 | 1.0/1.0 | -173.30 | 4/7 | 0.73/0.58 | no |
| g4_a1_d0 | 4.0/4.0 | 1.0/1.0 | -179.23 | 3/4 | 0.68/0.65 | no |
| g1.5_a1_dp20 | 1.8/1.2 | 1.2/0.8 | -177.03 | 7/5 | 0.56/0.99 | no |
| g2.5_a1_dp20 | 3.0/2.0 | 1.2/0.8 | -177.94 | 7/6 | 0.56/0.99 | no |
| g4_a1_dp20 | 4.8/3.2 | 1.2/0.8 | -178.99 | 7/6 | 0.56/0.99 | no |
| g1.5_a1_dm20 | 1.2/1.8 | 0.8/1.2 | -255.32 | 10/2 | 0.47/0.69 | no |
| g2.5_a1_dm20 | 2.0/3.0 | 0.8/1.2 | -205.24 | 7/0 | 0.62/- | no |
| g4_a1_dm20 | 3.2/4.8 | 0.8/1.2 | -205.48 | 8/0 | 0.59/- | no |
| g1.5_a0.3_dp20 | 1.8/1.2 | 0.36/0.24 | -25.00 | 0/0 | -/- | no |
| (all 17 other a<=0.3 rows) | ... | ... | -161..-201 | 6-9/0 | 0.41-0.54/- | no |

Key patterns:
1. **No NaN, no falls anywhere in the grid** (bal_com_z_min 0.85-0.90,
   tilt 14-22 deg) — strong drive never destabilized the rig-held model.
2. **Alpha is the left-leg lever, gamma is not.** alpha>=0.3 raises
   right-leg duty toward reference; alpha = 1.0 symmetric is the only
   thing that unfroze the left leg (6/6, 4/7, 3/4 cycles at gamma
   1.5/2.5/4.0).
3. **Asymmetry is destructive in BOTH directions**: dm20 (drive biased
   right) re-freezes the left leg (0 cycles at gamma 2.5/4.0); dp20
   (biased left) leaves it at duty 0.99 (planted, micro-unloads).
4. **More gamma = fewer cycles**: bilateral cycle count falls
   12 -> 11 -> 7 as gamma goes 1.5 -> 2.5 -> 4.0; gamma 4x buys nothing.
5. Knee flexion in the eval window is near zero everywhere
   (knee_min -1.3..-1.8 deg vs ref -69.7) — the eval-window gait is
   stiff-legged; see section 4 for the full-schedule picture.

## 4. Winner verification — full 22 s run (the "stays up 20 s" check)

`spinal/drive_sweep_full_20261016.py` (non-eval default schedule:
stand 0-3, ramp 3-5, walk 5-15, ramp down 15-17, stand 17-22 s; npz
redirected via AARL_NPZ). Log: `drive_sweep_human/fullrun.log`.

- **winner g1.5_a1_d0** (`winner_g1.5_a1_d0_full.npz`):
  runner summary "stayed up"; **stays_up_20s=True** (com_z min 0.870 m
  over the full 22 s). Walk-window kine: **bilateral = true**, 4 right +
  3 left cycles (T_r 1.29 s, T_l 1.86 s), duty_r 0.77 / duty_l 0.43,
  **knee_min_r/l -60.2/-67.2 deg** (ref -69.7/-70.6 — near-reference
  swing flexion appears in the production schedule), interleg lag
  lag_rl 0.449 (ref 0.508), double support 0.243 (ref 0.232),
  kine_score **-169.21** (better than its eval-window -173.10).
- runner-up g2.5_a1_d0 (`runnerup_g2.5_a1_d0_full.npz`): stays up
  (com_z min 0.871) but the full schedule re-freezes its left leg
  (15 right cycles, T 0.75 s, duty 0.38, frozen_l = true, kine -203.14)
  — **schedule sensitivity is real**: the same knobs walk differently
  under the 16 s eval schedule (4/7 cycles) vs the 22 s production
  schedule. The winner stays bilateral in BOTH windows.

## 5. Verdict against the ask's bar

For the winner (numbers from the sweep's eval window unless stated):
- >= 6 real cycles: **YES** — 6 right + 6 left = 12 real contact cycles
  (full 22 s window: 4 + 3 = 7, bilateral).
- duty 0.4-0.75: **YES** in the eval window (0.665 / 0.597); full window
  0.77 / 0.43 (right slightly above 0.75).
- stays up 20 s: **YES** — verified by the full 22 s run (com_z min
  0.870 m, never below the 0.55 fall threshold).
- STRONG ASYMMETRIC drive specifically: **NO** (section 0) — the
  question's literal subject fails; symmetry is what works.

Honest caveats: no forward progression (dx -0.02 m; COM x -0.21..-0.06 m
— it alternates in place under the support rig); hips ride ~22 deg
flexed; hip range 15 deg vs ref 43; kine -173/-169 is far from the
reference (0 = perfect) and below the one-legged base (-161.57) because
the bilateral gait pays real kinematic penalties the frozen-leg base
dodges. This is a bootstrap gait for phase 5+ (feedback tuning + drive
wean), not a finished walker. The "-189.5 of-record" target quoted in
the handoff is a gait_lib score (different metric) — not comparable to
these kine_ref numbers.

## 6. Files

Added by this stage (all untracked, no commits, no pushes):
- `spinal/drive_sweep_human_20261016.py`, `drive_sweep_one_20261016.py`,
  `drive_sweep_full_20261016.py`, `drive_sweep_analyze_20261016.py`,
  `drive_sweep_write_winner_20261016.py` (+ two .cmd launchers).
- `reports_20261016/drive_sweep_human/`: results.jsonl (28 rows),
  summary.csv, sweep.log, fullrun.log, winner_g1.5_a1_d0.json,
  runnerup_g2.5_a1_d0.json, winner_g1.5_a1_d0_full.npz,
  runnerup_g2.5_a1_d0_full.npz, this report.

Side effects: the full runs' PNG save (not redirectable) overwrote
tracked `spinal/spinal_run.png`; restored via
`git checkout -- Code/MuJoCo_SNS/spinal/spinal_run.png` (it was unmodified
before this stage). npz outputs were redirected (AARL_NPZ), so
`spinal_run.npz` was untouched by this stage (its pre-existing M state
is from earlier campaign stages). Scheduled tasks aarl_drive_sweep /
aarl_drive_fullrun deleted after use.

## 7. Reproduce

    cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
    set PYTHONIOENCODING=utf-8
    D:\Anaconda\envs\myo\python.exe drive_sweep_human_20261016.py   # full grid (resume-safe)
    D:\Anaconda\envs\myo\python.exe drive_sweep_analyze_20261016.py # table + summary.csv
    D:\Anaconda\envs\myo\python.exe drive_sweep_full_20261016.py reports_20261016\drive_sweep_human\winner_g1.5_a1_d0.json winner_full.npz

Winner json is JSON-RULE ready: full param set incl. gamma_l/gamma_r/
alpha_l/alpha_r/alpha_to_mn; loader branch needed by any future study:
`P.G[k] = float(p.get(k, 0.0))` for the four knobs after set_stage.

## 8. Not run (stated plainly)

- Basin gate on the winner (spinal/basin_gate.py) — finalist-robustness
  check belongs to the phase-10 supervisor; not part of this sweep ask.
- Grid on the flavor-B base (full_rules=1 production topology): not
  swept; the verify stage already recorded flavor B = -25.0 (frozen) at
  the winner config under current physics
  (DRIVE_VERIFY_20261016.md section 2/3).
- Drive-wean ladder (gamma/alpha -> 0 while re-tuning feedback) — that
  is phase 6, downstream of this bootstrap.
