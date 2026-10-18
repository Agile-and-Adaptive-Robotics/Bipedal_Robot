# HUMAN FEEDBACK-TUNING + DRIVE-WEAN CAMPAIGN (2026-10-16, overnight walk campaign)

Stage role (the ask): run the optuna campaign Ben's strategy prescribes —
tune feedback at high drive first, then re-optimize with the hand-drive
forced down (<=0.7, then <=0.4), keeping the rhythm ONLY if feedback
sustains it; report how far the drive weaned before gait died. Anti-lottery
requirements baked in. PASS bar: kine clearly better than the -189.5
of-record walking score, >=6 real cycles, duty 0.4-0.75, no NaN, stays up.

## 0. Headline

- **WINNER (finalist): `round8_best3.json` = curr_s8b_wean04 trial 35** —
  feedback re-tuned at drive_scale **0.381** (the <=0.4 ladder rung).
  - Eval window (16 s --eval, the metric family of the -161.57 anchor /
    -173.10 bootstrap / -189.5 of-record comparisons):
    **kine_score -145.07280220616852**, **10 + 10 contact cycles**,
    duty **0.71 / 0.65** (both inside the 0.4-0.75 band), no NaN, kz 0.85.
    Reproduced bit-exact in a fresh solo process (the reproduceCmd below).
  - Production window (full 22 s schedule, `round8_full3.npz`): stays up
    the whole 22 s (com_z min 0.848 m), **bilateral 9 + 10 cycles**,
    knee_min -42.8/-60.5 deg (ref -69.7/-70.6), interleg lag 0.463
    (ref 0.508), kine -164.43.
  - Reference set (wean_score_20261016.py, kine_ref.compare): subject01
    of-record **-145.07** (vs the s3k walker's -189.5 on the same metric),
    Ong 0.75 m/s **-159.66**, Ong 1.0 m/s **-173.27**, Ong 1.5 m/s
    -197.39 — the gain set generalizes across independent references.
- **WEAN ANSWER (the structural test)**: with the FINALIST's feedback the
  scale ladder never kills the gait — bilateral walking persists at EVERY
  scale down to **0.0** (8+11 cycles, -157.07 at scale 0, still better
  than the full-drive bootstrap -173.10), with the optimum at scale ~0.3
  (-138.93). Full hand-drive (scale 1.0) is actively HARMFUL (-215.12,
  left leg frozen). For contrast, the bootstrap combo's own feedback dies
  between scale 1.0 and 0.85 (probe_round1.jsonl), and the 39-trial
  FB-inert round 2 produced ZERO bilateral trials at any scale <=0.7.
  **Feedback re-tuning is what sustains the gait — exactly Ben's
  hypothesis.**
- Higher-scoring alternate (flagged, not the headline winner because its
  duties sit above the band): round8 trial 31, **-138.42 @ scale 0.400**
  (17/12 cycles, duty 0.82/0.82 eval; production -138.00, 18+11 cycles,
  stays up) — `round8_best.json`.
- Honest caveats: still in-place alternation under the support rig (no
  forward progression), hips ride flexed (mean_hip ~19/12 deg vs ref
  4.3/3.5), hip range ~15 deg vs ref 43, knee range 25-28 vs ref 70,
  double support high (ds 0.41-0.67 vs ref 0.23), eval scores under
  concurrent load drift ~5-8 points (BLAS nondeterminism; all numbers
  quoted here are from solo runs unless noted), and one pair of lottery
  weights was found and is reported (section 5).

## 1. Method

- Base = campaign anchor A (v10 multipliers + renshaw 0.5 + s3k trial-34
  winner params + full_rules=0, `set_stage(4, ...)`), gated bit-exact at
  **-161.56754173676563** twice this stage (before and after the driver
  fix); bootstrap knobs from the drive-sweep winner g1.5_a1_d0
  (eval -173.0970129592704, reproduced exact in the same gates).
- Driver: `spinal/wean_tune_20261016.py` (+ `wean_round_launch_20261016.py`,
  4 workers via schtasks, the drive-sweep survival pattern). Searched
  space (20 keys): gamma_l/r, alpha_l/r, alpha_to_mn + 13 feedback keys
  (heel_rge, toe_rge, ib_rge, ia_in, ib_e/ii_e/ii_f/ia_f_central,
  phase_reset_e/f, contact_onset, f1_kneext_inh, f1_anklepf_inh) — all
  knobs the ask names. Objective = the established stage-4 eval score
  (kine_score, NaN -400 / no-cycles -320 sentinels, kz/tilt adjustments)
  + a **BILATERAL GATE** (-60 unless >=3 contact cycles on BOTH legs).
- Wean semantics: `drive_scale` multiplies ONLY the hand-drive knob
  currents (u[GAMMA_c_s]=drive*scale*G[gamma_s], u[ALPHA_c_s] likewise,
  runner.py:1366-1378); the s3k descending DRIVE stays at 2.0773.
- Anti-lottery item 1: every trial logs the FULL per-term kine diagnostics
  (rmse/mean/range per joint per leg, T, period_cv, lag_rl, ds, duty,
  knee_min, cycles, frozen flags) to `round{r}_w{worker}.jsonl`.
- Studies (optuna_walk.db, TPE, seeded trial 0 = incumbent):
  curr_s6b (100 trials), curr_s7 (39, aborted), curr_s6c (6, aborted),
  curr_s6d_feedback (88), curr_s7b_wean07 (80), curr_s6e_prodwindow (64),
  curr_s8b_wean04 (48).

## 2. The mid-campaign bug (honest record)

Rounds 1-2-4 ran with a driver bug: the sampled feedback keys were never
injected — `set_stage` re-loaded the s3k values every trial, so those
searches varied ONLY the drive knobs (+ scale). Symptom caught: round 2's
"feedback re-tune at <=0.7" produced 0/39 bilateral trials at scales
0.31-0.7 — exactly what a scale-only sweep predicts (probe_round1). The
bug was found by inspecting apply_config, fixed (fb merged into base
BEFORE set_stage), the gate re-passed bit-exact, and rounds 5-8 re-ran
the designs correctly:
- round 5 = redo of round 1 (eval window, 20 live keys)
- round 6 = redo of round 2 (wean <=0.7, seeded from round 5)
- round 7 = production-window objective (22 s schedule walk window
  scored; round-1 winners had walked in the eval window but froze the
  left leg in production — round1_best: eval -146.03, full run -231.92)
- round 8 = the <=0.4 rung (Ben's ladder), seeded from round 6.
The FB-inert rounds remain in the db + wean_human/round{1,2,4}_* as the
record; round-1's -146.03 is still reported as a (drive-only-tuned,
production-marginal) data point.

## 3. Round results (bilateral = >=3+3 cycles; scores = study objective)

| round | study | design | trials | bilateral | best bilateral |
|---|---|---|---|---|---|
| 1 | curr_s6b | eval, FB-inert (bug) | 100 | 51 | -146.03 (t87, drive-only tuning) |
| 2 | curr_s7 | wean<=0.7, FB-inert (bug) | 39 | 0 | none — invalid as re-tune |
| 4 | curr_s6c | prodwindow, bug | 6 | 2 | aborted |
| 5 | curr_s6d | eval, 20 keys live | 88 | 26 | -163.30 (t56) |
| 6 | curr_s7b | wean<=0.7, FB live | 80 | 48 | **-148.71 @ scale 0.452 (t34)** |
| 7 | curr_s6e | production window | 64 | 14 | -158.77 (t14; duty 0.94 high) |
| 8 | curr_s8b | wean<=0.4, FB live | 48 | 24 | **-138.42 @ 0.400 (t31); -145.07 @ 0.381 (t35 = finalist)** |

Production-schedule verification (wean_full_20261016.py, 22 s default
schedule, npz under wean_human/):
- finalist t35 @ 0.381: stays up (com_z min 0.848), bilateral 9+10,
  kine -164.43, knee -42.8/-60.5, lag 0.463.
- t35 gains @ 0.30: stays up (0.848), bilateral 8+10, kine -164.64.
- t31 @ 0.400: stays up (0.847), bilateral 18+11, kine -138.00 (duty
  0.89/0.75, ds 0.67 — draggy).
- round6 t34 @ 0.452: stays up (0.896), bilateral 9+7, kine -153.38,
  knee -55.0/-64.8, lag 0.475.
- round1 t87 @ 1.0 (FB-inert): stays up BUT left frozen (13/0, -231.92
  solo; -166.11 under load — marginal basin, see section 6).

## 4. Wean ladder answer

- Bootstrap feedback fixed (probe_round1.jsonl, round-1 winner gains):
  scale 1.0 -146.03 (5/6) -> 0.85 -221.57 (2/2) -> 0.7 -247.76 (6/2,
  left duty 0.99) -> <=0.55 one-legged -> 0.0 = anchor -161.57 exactly.
  **Death between 1.0 and 0.85.**
- FINALIST feedback fixed (probe.jsonl, round8_best3 gains):
  1.0 -215.12 (9/0!) / 0.7 -178.67 (10/4) / 0.55 -169.87 (11/6) /
  0.45 -157.88 (9/8) / 0.38 -152.82* (9/8) / 0.3 **-138.93 (8/9)** /
  0.2 -149.27 (6/8) / 0.1 -156.74 (10/8) / **0.0 -157.07 (8/11)**.
  (*concurrent-run number; solo = -145.07.) **No death point — bilateral
  gait at every scale including zero hand-drive.** At scale 0 the
  GAMMA/ALPHA cells are not even built (conditional topology): this is
  the anchor network + re-tuned feedback alone, and it beats the
  full-drive bootstrap (-157.07 vs -173.10).

## 5. Anti-lottery evidence (Ben's four requirements)

1. **Per-term diagnostics logged per trial** — every jsonl row carries
   rmse/mean/range for hip/knee/ankle x r/l, T, period_cv, lag_rl, ds,
   duty, knee_min, cycles, frozen flags (wean_human/round*_w*.jsonl).
2. **Gain-sensitivity probe** (sens.jsonl; baseline -145.07; each row =
   one key reset to its s3k/bootstrap seed value or zeroed, fresh eval):
   - STRUCTURAL (removal degrades): ia_in zero -> **-250.90 (5/0, left
     leg dies)**; ii_e_central seed -> **-25.0 (frozen sentinel)**, zero
     -> -173.95; f1_kneext_inh seed/zero -> -173.01/-183.26; heel_rge
     -> -156.38/-151.02; alpha_to_mn -> -149.41/-156.13. ia_in's value
     sits on a plateau (0.126 vs seed 0.825 score identically; 0 fails).
   - **LOTTERY WEIGHTS (removal does NOT degrade)**: toe_rge — removal
     IMPROVES the score (+13.0 seed / +14.2 zero); phase_reset_f (value
     2.29!) — zeroing changes nothing (-145.07, 10/10 cycles). Reported
     as such; neither is claimed as structure. phase_reset_e is ~0.06,
     effectively inert at this operating point.
3. **The wean ladder IS the structural test** — section 4: the winner
   walks at every drive scale including 0 (feedback-sustained), while the
   non-re-tuned bootstrap dies below 0.85. A full-drive-only winner would
   have failed exactly as round-1's t87 did (production re-freeze).
4. **Cross-reference sanity check** — finalist scored against Ong
   predicted walking at speeds different from the subject01 primary
   (wean_score_20261016.py): ong_speed_075 **-159.66**, ong_speed_100
   **-173.27**, ong_speed_150 -197.39 (eval-window npz; production npz:
   -164.43 / -189.94 subject01/ong_100). A lottery ticket trained on one
   reference would collapse on the others; this degrades gracefully and
   actually scores BEST against the slower 0.75 m/s gait.

## 6. Basin gate

- Network-only constant-DRIVE gate (basin_gate_drive_20261016.py, the
  basin_gate.py pattern + hand-drive ports): **inapplicable here** — the
  UNPERTURBED baseline itself does not oscillate (5 s-window swings decay
  8.3 mV -> 0, 0 peaks in the last 10 s; basin_round8_best3.json.json).
  This is the documented s3k-lineage property (AGENTS 2026-09-20:
  "tuned-config constant-DRIVE self-sustain is LOST (schedule-
  sustained)"). The rhythm is closed-loop: it needs the body + schedule.
- **CLOSED-LOOP BASIN GATE (basin_closed_20261016.py: basin_gate.perturb
  +/-1% on every scalar param, then the FULL ground eval; PASS = bilateral
  >=3+3 cycles survives): 8/8 PASS** (basin_closed_round8_best3.json.json).
  Perturbed scores -132.28..-158.99, cycles 6-12 per leg, no NaN, no
  falls; two perturbations IMPROVED the score (-132.28, -132.78) — the
  winner sits on a robust plateau, not a knife edge. Baseline
  (unperturbed) reproduced -145.07280220616852 yet again in this gate.
- Related honest observation: round-1's t87 production score flips
  between -231.92 (solo) and -166.11 (under load) — execution-level
  nondeterminism on marginal configs is real; the finalist's solo-vs-
  study numbers agree bit-exactly, and its ladder/probe rows are
  consistent under both conditions.

## 7. Reproduce (tested this session, bit-exact)

    ssh -o BatchMode=yes easteregg2 "cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal && set PYTHONIOENCODING=utf-8&& D:\Anaconda\envs\myo\python.exe wean_tune_20261016.py one reports_20261016\wean_human\finalist_spec.json"

(finalist_spec.json = {"winner": "round8_best3.json", "tag":
"finalist_eval", "npz": "finalist_eval.npz"}; observed RESULT kine_score
-145.07280220616852, identical to the study value.)
Production-schedule check: `wean_full_20261016.py reports_20261016\
wean_human\round8_best3.json <npz>` -> stays up, bilateral 9+10.

## 8. Files (all untracked; NO commits, NO pushes)

- spinal/: wean_tune_20261016.py, wean_round_launch_20261016.py,
  wean_probe_20261016.py, wean_full_20261016.py, wean_score_20261016.py,
  wean_sens_20261016.py, wean_poll_20261016.py, basin_gate_drive_
  20261016.py, basin_closed_20261016.py, wean_{gate,r1..r8,probe,probe2,
  sens}.cmd (launchers).
- reports_20261016/wean_human/: round*_w*.jsonl (per-trial diagnostics),
  round{1,5,6,7,8}_best*.json (winner param files), probe_round1.jsonl +
  probe.jsonl (ladders), sens.jsonl, basin_*.json, *_task.log/round*.log,
  finalist_spec.json, finalist_eval.npz, round{1,6,8}_full*.npz, this
  report.
- spinal_run.png: full runs churn it (not redirectable) — RESTORED via
  `git checkout -- spinal_run.png` (verified: no longer in git status).
  spinal_run.npz: tracked file was ALREADY M from earlier campaign
  stages; this stage's un-redirected probe/sens evals overwrote its
  working-tree content (round/prod runs WERE redirected via AARL_NPZ).
  The committed v6-winner npz remains recoverable from git history.

## 9. Not run (stated plainly)

- The <=0.4 rung used round-6's winner as seed; a <=0.2 rung was not
  run (the finalist ladder already probes 0.2/0.1/0.0 fixed-feedback,
  which answers the question the rung would ask).
- flavor-B (full_rules=1) grid — per the handoff, flavor B is recorded
  frozen at -25.0 under current physics; out of scope here.
- Forward-progress objective — the winner still alternates in place;
  no dx term was added to the objective (architecture-level, next study).
