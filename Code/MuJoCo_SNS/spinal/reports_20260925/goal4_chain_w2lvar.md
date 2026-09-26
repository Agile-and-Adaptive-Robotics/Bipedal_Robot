# goal4 chain log — `w2lvar` variant curriculum (5 stages)

Appended per stage as it is harvested. Study/db/npz per
`reports_20260925\goal4_wiring.md`; scripts `_curriculum_w2lvar.py`
(self-pins `AARL_NET=w2lvar`, `AARL_NPZ=spinal_run_w2lvar.npz`, db
`optuna_w2lvar.db`). Stage-1 launch: detached cmd
`reports_20260925\tmp\launch_curr_w2lvar_s1.cmd`, 2026-09-25 14:30:35,
18 trials on top of the 2 smoke trials.

## Stage 1 — `curr_w2lvar_s1_air_deaff` (deafferented air) — DONE

- **Run**: exit code 0 (`logs\curr_w2lvar_s1.done` line 1); log
  `logs\curr_w2lvar_s1.log` (122,440 B, last write 15:59:32); winner json
  `curriculum_w2lvar_stage1.json` written by the process itself.
- **Trials**: 20/20 COMPLETE in `optuna_w2lvar.db` (2 goal4 smoke + 18
  campaign, trial numbering continuous). Values (dump via
  `tmp\harvest_w2lvar_s1.py`): 0 NaN sentinels (<=−199.5), 9 rhythm-gate
  floors (−10.5..−3), 11 real-rhythm scores (>9).
- **Winner**: trial 17, **score 80.836** — stage-1 air objective
  `3·rises + 0.5·(−knee_min)` over t∈[5,17], higher is better,
  rhythm-gated (`_curriculum_w2lvar.py:213-216`). Params: drive 2.7179,
  rg_nap_h 0.2370, desc_e 1.6212, desc_f 2.0643, rg_to_pf 2.8221.
- **Exploit check: PASS (not exploited)** — re-ran the winner params
  through the exact stage-1 objective path (`tmp\harvest_w2lvar_s1.py`,
  log `logs\harvest_w2lvar_s1.log`): rises=14, knee_min=−77.67°,
  RG-E range 9.75 mV, recomputed score 80.835587 vs json 80.835587,
  **delta 0.000000** (bit-exact repro). Runner summary in the same run:
  RG_r 14 cycles, period 0.70 s (1.43 Hz), E-duty 0.45; knee_angle_r
  −78.3..+12.6 deg. Not on any sentinel floor.
- **Protection**: `spinal_run.npz` untouched (size 8,238,797, mtime
  2026-09-25 01:10:02, pre-session baseline); chain artifacts
  `spinal_run_w2lvar.npz` + db updated 15:59:32. The log's
  "saved spinal_run.npz" line is the runner's hardcoded print string;
  the write went to the chain npz (AARL_NPZ).
- **Diagnosis**: genuine deafferented air rhythm — 14 bursts/12 s window
  at 1.43 Hz with deep knee swing (−78°); winner beats the gate-c seed
  (56.002) by +44%; drive rose 2.5→2.72, nap_h dropped 0.35→0.24.
  Cadence is fast vs the 0.3 Hz Ivanenko air-stepping regime — the
  objective rewards burst count, so stage 2+ (afferents + interleg,
  then ground) is where duty/cadence realism gets shaped. Stage 2
  (`curr_w2lvar_s2_air_aff`) may proceed; it seeds KEYS1 from this json
  and the reflex keys at defaults 0.6/0.4/0.4/0.35
  (`_curriculum_w2lvar.py:289-304`).

## Stage 2 — `curr_w2lvar_s2_air_aff` (afferented air) — DONE

- **Run**: launched detached 16:10:43 (`tmp\launch_curr_w2lvar_s2.cmd`,
  18 trials = trial 0 seed + 17 TPE; NOTE: `optimize(n_trials=18)`
  INCLUDES the enqueued seed — the stage-1 launch note's "19 total" was
  wrong). The orchestrator's ~88 min poll budget expired pre-finish;
  process verified still progressing (15/18 runs at 17:37, log growing),
  waited a further bounded ~15 min; finished during that window. Exit
  code 0 (`logs\curr_w2lvar_s2.done` line 1); log
  `logs\curr_w2lvar_s2.log` (131,633 B, mtime 17:53:23); winner json
  `curriculum_w2lvar_stage2.json` written by the process itself.
- **Trials**: 18/18 COMPLETE (`tmp\harvest_w2lvar_s2.py` dump): 0 NaN
  sentinels, 13 at/below the rhythm-gate ceiling (<8), 5 real rhythms
  (>8). **Seed trial 0 = stage-1 winner + afferents/interleg ON at
  defaults scored −9.804 — LOST the rhythm**: the deafferented optimum
  is not drop-in compatible with the afferented stack; TPE had to
  re-find a rhythm from scratch.
- **Winner**: trial 13, **score 60.230** — same stage-1 air objective
  form under stage-2 conditions (`--no-ground`, afferents + interleg
  ON; `_curriculum_w2lvar.py:181-182`). Params: drive 2.2692,
  rg_nap_h 0.3179, desc_e 1.4390, desc_f 1.6949, rg_to_pf 1.8716,
  ia_to_mn 0.3564, ia_to_antagonist 0.3342, ii_to_mn 0.0480,
  ib_to_mn_inh 0.1511.
- **Exploit check: PASS (not exploited)** — exact objective-path repro
  (`tmp\harvest_w2lvar_s2.py`, log `logs\harvest_w2lvar_s2.log`):
  rises=10, knee_min=−60.46°, RG-E range 7.88 mV, recomputed
  60.230493 == json 60.230493, **delta 0.000000** (bit-exact).
- **Protection**: `spinal_run.npz` still untouched (size 8,238,797,
  mtime 2026-09-25 01:10:02) after stage 2 AND the harvest repro;
  chain npz updated 18:03:13 by the repro.
- **Diagnosis**: afferented air rhythm is genuine but SLOWER/shallower
  than deafferented (10 bursts @ knee −60° vs stage-1's 14 @ −78°;
  60.2 vs 80.8, −25%) — and TPE pushed the reflex gains DOWN
  (ii_to_mn 0.048, ib_to_mn_inh 0.151, ia_to_mn 0.356 vs defaults
  0.4/0.35/0.6): in AIR, length/load afferents mostly fight the
  rhythm; their value should appear on the ground (stages 3-5).
  Usable winner -> stage 3 (`curr_w2lvar_s3_balance`) may proceed;
  note stage 3 stands alone (keys vest_prop/rig_scale seed at
  0.0/1.0; the rhythm keys are NOT carried into its seed — script
  design at `_curriculum_w2lvar.py:62-64,86`).

## Stage 3 — `curr_w2lvar_s3_balance` (standing balance) — DONE

- **Run**: launched detached 18:05:30 (`tmp\launch_curr_w2lvar_s3.cmd`),
  18 trials; exited cleanly before harvest (exit code 0,
  `logs\curr_w2lvar_s3.done`; log 10,544 B, mtime 18:31:35 — standing
  evals are much faster than the 14 s air evals).
- **Trials**: 18/18 COMPLETE (`tmp\harvest_w2lvar_s3.py` dump): values
  25.08, −0.99, −150.0, −4.95, 14.46, 0.01, −17.77, −23.84, −4.60,
  6.71, −0.70, 11.41, 18.66, 18.28, 1.96, 13.96, 15.20, 21.86 →
  0 NaN sentinels, 1 fall (−150), 17 standers. NOTE: the stage-3 seed
  pins ONLY vest_prop/rig_scale; the objective suggests the 9
  rhythm/reflex keys regardless, so trial 0's other keys were TPE-
  sampled (drive 3.888, rg_nap_h 0.675, ...) — visible in the winner
  json's 11-key params dict.
- **Winner**: trial 0 (the first, random, sample), **score 25.077** —
  standing objective `100 − 400·sway − tilt_max − 40·|0.5−contact_sym|`
  over 8 s `--stand-eval` at rig_scale 1.0, higher is better, NaN
  −200 / fall −150 sentinels (`_curriculum_w2lvar.py:220-233`).
  vest_prop 0.0, rig_scale 1.0.
- **Exploit check: PASS (genuine stander)** — exact objective-path repro
  (`tmp\harvest_w2lvar_s3.py`, log `logs\harvest_w2lvar_s3.log`):
  bal_fell false, no NaN, sway 0.0996 m, tilt_max 25.32°, contact_sym
  0.2562, recomputed 25.077367 == json 25.077367, **delta 0.000000**
  (bit-exact); not on any sentinel floor. (The air-stage rises<3 check
  does not apply to a standing stage — standing still is the goal.)
- **Protection**: `spinal_run.npz` still size 8,238,797, mtime
  2026-09-25 01:10:02 after stage 3 + harvest repro.
- **Diagnosis**: standing WORKS but is weak — the best of 18 is the
  FIRST random sample and TPE could not beat it in 17 tries; the
  winner survives 8 s at full rig with ~10 cm sway, 25° peak tilt and
  poor weight symmetry (0.256 vs 0.5 ideal) ≈ 25/100 on the SCONE-3a
  analog scale. Rhythm/posture params for standing were sampled, not
  inherited (stage 3 stands alone by design), so the standing map is
  under-explored — a candidate for a later re-tune, not a chain
  blocker. Usable winner -> stage 4 (`curr_w2lvar_s4_walk_nocontact`,
  22 trials) may proceed.

## Stage 4 — `curr_w2lvar_s4_walk_nocontact` (walking, contact OFF) — DONE

- **Run**: launched detached 18:40:45 (`tmp\launch_curr_w2lvar_s4.cmd`),
  22 trials ADDED to the 1 goal-4 smoke trial (no seed re-enqueued —
  study was non-empty; confirmed by the absence of a "seeded" line in
  the log). Exited cleanly before harvest: exit code 0
  (`logs\curr_w2lvar_s4.done`); log 11,990 B, mtime 19:44:02 (~2.7
  min/trial — the `--eval` route prints far less than the air runs).
- **Trials**: 23/23 COMPLETE (`tmp\harvest_w2lvar_s4.py` dump): values
  −315.0 (smoke), then 16× the frozen sentinel −320.0 (no countable
  cycles), 5× the −315 clip, and exactly 2 real-kine scores: −271.959
  (trial 15) and −225.891 (trial 17 = winner). 0 NaN (−400).
- **Winner**: trial 17, **score −225.891** — stage-4 kine objective
  `max(kine_score, −315)` on `--no-ground --eval` (walking pattern-
  matched against the reference cycle in air; lower kine distance is
  better, study maximizes; NaN −400 / frozen −320 sentinels;
  `_curriculum_w2lvar.py:234-256`). 12-key params: drive 1.6783,
  rg_nap_h 0.8025, desc_e 1.5328, desc_f 2.1306, rg_to_pf 2.4605,
  ia_to_mn 0.5688, ia_to_antagonist 0.1621, ii_to_mn 0.5996,
  ib_to_mn_inh 0.2999, rg_mutual_inh 3.6819, rg_weak_exc 0.1150,
  pf_to_mn 0.8721.
- **Exploit check: PASS (real cycle pattern)** — exact objective-path
  repro (`tmp\harvest_w2lvar_s4.py`, log `logs\harvest_w2lvar_s4.log`):
  kine_score −225.891, tilt_max 1.07° (no tilt penalty), duty 0.6404,
  recomputed −225.890865 == json −225.890865, **delta 0.000000**
  (bit-exact); not on the −320/−400 floors nor the −315 clip. The
  kine score is only computed when the eval detects cycles, so the
  winner is a genuine gait pattern — and the 16 frozen trials prove
  motionlessness scores −320, far below the winner.
- **Protection**: `spinal_run.npz` still size 8,238,797, mtime
  2026-09-25 01:10:02; `optuna_walk.db` untouched (20 studies, mtime
  2026-09-25 01:10:02 vintage) after stage 4 + harvest repro.
- **Diagnosis**: air pattern-matching is HARD for this variant (only
  2/23 trials produced any countable gait), but the winner's pattern
  is genuinely good: kine −225.9 sits numerically in the stock s3k
  production walker's ground band (−222…−241 vs subject01/Falisse,
  AGENTS.md 2026-09-23) and its duty 0.640 nearly matches the
  reference 0.61 — achieved with NO contact forces (context, not an
  apples-to-apples ground comparison). TPE's search drove drive DOWN
  (1.68) and nap_h UP (0.80) from the air-rhythm optimum — a
  different regime than stages 1-2. Usable winner -> stage 5
  (`curr_w2lvar_s5_walk_contact`, 22 trials, the final stage) may
  proceed; stage 5 seeds KEYS5 from this json + contact keys at
  defaults.

## Stage 5 — `curr_w2lvar_s5_walk_contact` (walking WITH contact) — DONE — CHAIN COMPLETE

- **Run**: launched detached 19:53:28 (`tmp\launch_curr_w2lvar_s5.cmd`),
  22 trials; fresh study, full-dict seed = the 12 stage-4 winner keys +
  contact keys at defaults (contact_onset 0.0, contra_swing 0.0,
  ib_group_exc 0.5, ib_exc_to_mn 0.6). Exited cleanly before harvest:
  exit code 0 (`logs\curr_w2lvar_s5.done`); log 12,632 B, mtime
  20:43:26 (~50 min).
- **Trials**: 22/22 COMPLETE (`tmp\harvest_w2lvar_s5.py` dump): values
  −220.40, −320.0, −320.0, −320.0, **−165.04**, −275.39, −320.0,
  −246.08, −272.70, −276.81, −259.93, −259.74, −320.0, −172.11,
  −221.25, −274.80, −320.0, −320.0, −320.0, −248.25, −320.0, −320.0 →
  0 NaN (−400), 10 frozen (−320), 0 at the clip, **12/22 real ground
  kine scores** — vs 2/23 in air at stage 4: ground contact HELPS this
  variant's pattern matching.
- **Winner**: trial 4, **score −165.041** — stage-5 objective
  `max(kine_score, −315) − 20 if kz<0.62 − 10 if tilt>40` on the
  normal ground `--eval` (NaN −400 / frozen −320 sentinels;
  `_curriculum_w2lvar.py:238-256`). Full 16-key params: drive 3.6740,
  rg_nap_h 0.7275, desc_e 1.0396, desc_f 0.8557, rg_to_pf 2.0236,
  ia_to_mn 0.6774, ia_to_antagonist 0.5717, ii_to_mn 0.4743,
  ib_to_mn_inh 0.0933, rg_mutual_inh 4.0482, rg_weak_exc 0.0407,
  pf_to_mn 1.8081, **contact_onset 0.1339, contra_swing 0.2979,
  ib_group_exc 0.0866, ib_exc_to_mn 0.5695** — TPE adopted nonzero
  contact-pathway gains and pushed drive UP to 3.67 for ground
  walking (stage-4 air optimum used drive 1.68).
- **Exploit check: PASS (genuine upright ground walk)** — exact
  objective-path repro (`tmp\harvest_w2lvar_s5.py`, log
  `logs\harvest_w2lvar_s5.log`): kine_score −165.041 with ZERO
  penalties — kz 0.8297 (COM at 83% height, far above the 0.62 floor;
  no kz penalty), tilt_max 27.72° (under 40; no tilt penalty), duty
  0.4597; recomputed −165.040506 == json −165.040506, **delta
  0.000000** (bit-exact). Not on any sentinel floor; the kine score
  requires detected cycles, and 10 frozen trials prove motionlessness
  scores −320.
- **Protection**: `spinal_run.npz` still size 8,238,797, mtime
  2026-09-25 01:10:02; `optuna_walk.db` untouched (mtime 2026-09-25
  01:10:02 vintage) after stage 5 + harvest repro.
- **Diagnosis / chain outcome**: the w2lvar chain COMPLETED 5/5 stages.
  Its ground walker: kine −165.0 raw (no penalties), upright at kz
  0.83, duty 0.46 (reference 0.61), tilt 27.7° — the best pattern
  score of the whole chain and, for context only (different
  reference/scoring routes), numerically better than the stock s3k
  walker's recorded −189.5…−241 kine band; all five winners chain via
  `curriculum_w2lvar_stage{1..5}.json` in `optuna_w2lvar.db` studies
  `curr_w2lvar_s{1..5}_*`. Candidate caveats for Ben's review: duty
  still 0.46 vs human 0.61; tilt 27.7° is high; stage 3's standing
  map was under-explored (winner was a random sample); the winner
  config has NOT been basin-gated (basin_gate.py is tuned to the stock
  chain's dump state).
