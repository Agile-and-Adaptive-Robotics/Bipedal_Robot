# goal4 chain — syn6 variant curriculum harvest (2026-09-25)

One section per stage. Source studies live in `spinal\optuna_syn6.db`
(chain-private); winners in `spinal\curriculum_syn6_stage{N}.json`; run logs
in `reports_20260925\logs\`. Stage launch = the detached driver
`reports_20260925\tmp\run_curr_syn6_s1.cmd` pattern (18 trials stages 1-3 /
22 stages 4-5; wiring in `reports_20260925\goal4_wiring.md`).

## Stage 1 — `curr_syn6_s1_air_deaff` (deafferented air stepping) — DONE, REAL WINNER

- **Run**: `_curriculum_syn6.py 1 18` detached, 14:30:49–15:47:59; exit code 0
  (`logs\curr_syn6_s1.done` first line = 0); log
  `logs\curr_syn6_s1.log` (122,639 bytes). Study resumed over the 2 smoke
  trials → **20 trials COMPLETE** (sqlite census:
  `reports_20260925\tmp\harvest_check_syn6_s1.py` output), no new seed
  enqueued (study non-empty; `_curriculum_syn6.py:273`).
- **Score landscape**: 12/20 real rhythms (>25), 7 rhythm-gate (−5.3…−10.4),
  1 NaN (−200, trial 18), 0 frozen (−320). Winner 95.639 (trial 17) strictly
  above the seed 61.908 (trial 0) — no sentinel inversion.
- **Winner** (`curriculum_syn6_stage1.json`): score **95.639**, trial 17,
  params drive=3.3937, rg_nap_h=0.2350, desc_e=1.3654, desc_f=1.7006,
  rg_to_pf=2.5460 (+ pinned syn6=1.0, syn6_brainstem=0.0, renshaw=0).
- **Exploit check: PASS (REAL RHYTHM)** — winner config replayed bit-exactly
  (`tmp\replay_syn6_s1_winner.py` → `logs\curr_syn6_s1_exploit_check.log`):
  recomputed score 95.6388 vs json 95.6388 (**delta +0.0000**), RG_E_r
  rises=**16** (gate needs ≥3), span 9.74 mV, knee −95.3..+10.3° in t=5..17
  window; runner readout "cycles in window 15, period 0.64 s (1.57 Hz),
  E-duty 0.45". Not a static pose; not on any sentinel floor.
- **One-line diagnosis**: deafferented air rhythm converged at higher drive
  (3.39 vs seed 2.5) + faster Na-h recovery (nap_h 0.235): 1.57 Hz antiphase
  air-stepping with deep knee flexion (−95°) — healthy seed for stage 2.
- **Protection**: protected `spinal_run.npz` untouched (8,238,797 bytes,
  mtime 2026-09-25 01:10:02 AM, re-verified post-run); chain npz
  `spinal_run_syn6.npz` (15.1 MB, 15:47:59) is the real save target
  (AARL_NPZ; the runner's "saved spinal_run.npz" print at runner.py:1768 is
  a hardcoded literal). Replay scratch npz deleted after harvest.

**PROCEED → stage 2 (`curr_syn6_s2_air_aff`, afferented air, 18 trials).**

## Stage 2 — `curr_syn6_s2_air_aff` (afferented air, interleg ON) — DONE, REAL WINNER

- **Run**: `_curriculum_syn6.py 2 18` detached 16:00:23–17:28:34; exit code 0
  (`logs\curr_syn6_s2.done` first line = 0); log `logs\curr_syn6_s2.log`
  (122,924 bytes). Fresh study → **18/18 trials COMPLETE** (sqlite census:
  `tmp\harvest_check_syn6_s2.py` output); seed = stage-1 winner 5 keys +
  heel_rge/toe_rge 0.0 + syn6 1.0 (log "seeded" line — full dict, JSON-rule).
- **Score landscape**: 10/18 real rhythms (>25), 8 rhythm-gate
  (−7.1…−10.4), 0 NaN, 0 frozen. Winner 86.751 (trial 12) above the
  afferented seed 80.742 (trial 0) — no sentinel inversion. Note the seed
  re-scored 80.7 vs its stage-1 95.6: turning afferents+interleg ON costs
  rhythm quality (expected).
- **Winner** (`curriculum_syn6_stage2.json`): score **86.751**, trial 12,
  params drive=3.5551, rg_nap_h=0.2610, desc_e=1.3080, desc_f=1.8586,
  rg_to_pf=2.6804, **heel_rge=0.0110, toe_rge=0.0430** — both new gains
  tuned near-zero (same "afferents tolerated, not additive" outcome as the
  09-20 s2b stock chain).
- **Exploit check: PASS (REAL RHYTHM)** — winner replayed bit-exactly
  (`tmp\replay_syn6_s2_winner.py` → `logs\curr_syn6_s2_exploit_check.log`):
  recomputed score 86.7508 vs json 86.7508 (**delta +0.0000**), RG_E_r
  rises=**21** (gate needs ≥3), span 9.78 mV; runner readout "cycles in
  window 20, period 0.48 s (2.06 Hz), E-duty 0.36". Not a static pose; not
  on any sentinel floor.
- **One-line diagnosis**: afferented rhythm held (2.06 Hz, 21 bursts) but
  knee swing shallowed to −47.5..−12.3° (window) vs stage-1's −95..+10° —
  the winner trades flexion depth for burst count; heel/toe gains ≈0 say
  the contact ports add nothing in air, as in the stock s2b. Watch knee
  depth again at stages 4-5 (walking pattern-match).
- **Protection**: protected `spinal_run.npz` untouched (re-verified
  post-stage: 8,238,797 bytes, mtime 2026-09-25 01:10:02 AM); chain-private
  `spinal_run_syn6.npz` (15.7 MB, 17:28:33) is the real save target via
  AARL_NPZ (runner.py:1629); replay scratch npz deleted (Test-Path False).

**PROCEED → stage 3 (`curr_syn6_s3_balance`, standing balance, 18 trials).**

## Stage 3 — `curr_syn6_s3_balance` (standing balance, 8 s stand-eval) — DONE, GENUINE STAND

- **Run**: `_curriculum_syn6.py 3 18` detached 17:39:44–17:57:04; exit code 0
  (`logs\curr_syn6_s3.done` first line = 0); log `logs\curr_syn6_s3.log`
  (13,101 bytes — shorter than the air stages because 8 s standing evals
  print less and 5 trials ended in falls with MuJoCo instability warnings,
  which is what the done-marker tail shows). Fresh study → **18/18 trials
  COMPLETE** (sqlite census: `tmp\harvest_check_syn6_s3.py` output); seed =
  {vest_prop 0.0, rig_scale 1.0, syn6 1.0} (stage-2 winner has neither key —
  chain-merge correctly a no-op, log "seeded" line).
- **Score landscape**: 13/18 stood above the floors, 5 fell (−150 fall
  sentinel), 0 NaN-terminated. Winner 37.463 (trial 11) above the seed
  34.519 (trial 0) — no sentinel inversion, not on a floor.
- **Winner** (`curriculum_syn6_stage3.json`): score **37.463**, trial 11,
  params drive=3.6202, rg_nap_h=0.8356, desc_e=1.5719, desc_f=1.1327,
  rg_to_pf=1.8988, heel_rge=0.4165, toe_rge=0.4775, **vest_prop=0.3566,
  rig_scale=0.9459** (the objective samples the 5 core + heel/toe keys in
  every stage — `_curriculum_syn6.py:130-141` — so the 9-key dict is by
  design).
- **Validity check: PASS (GENUINE STAND)** — winner replayed bit-exactly
  (`tmp\replay_syn6_s3_winner.py` → `logs\curr_syn6_s3_exploit_check.log`):
  nan=False, bal_fell=False, bal_sway max 0.0878 rad / rms 0.0610,
  bal_tilt_max 18.76°, bal_contact_sym 0.284, kz 0.903; recomputed
  100−400·sway−tilt−40·|0.5−sym| = **37.4630 == json 37.4630 (delta
  +0.0000)**.
- **One-line diagnosis**: 8 s stand holds at rig_scale 0.946 (≈95% support —
  barely weaned; the search keeps rig high because lower rig directly lowers
  score) with side-biased loading (sym 0.28) and 18.8° max tilt; the
  mechanosensor gains that were dead in air (stage 2) won LARGE here
  (heel 0.42 / toe 0.48) — contact ports matter once feet load. CAVEAT for
  stage 4: per the documented chain (`_curriculum_syn6.py:281-295`) stage 4
  seeds its 7 shared keys from THIS stage-3 json (nap_h 0.836 — a
  standing-selected core, vs the walking-selected stage-2 nap_h 0.261); if
  the stage-4 seed underperforms, re-seeding from stage 2 is a legitimate
  later correction, not a bug in the wiring.
- **Protection**: protected `spinal_run.npz` untouched (re-verified:
  8,238,797 bytes, mtime 2026-09-25 01:10:02 AM); chain-private
  `spinal_run_syn6.npz` (4.0 MB, 17:57:04) is the real save target via
  AARL_NPZ (runner.py:1629); replay scratch npz deleted (Test-Path False).

**PROCEED → stage 4 (`curr_syn6_s4_walk_nocontact`, contact-free walking
kine, 22 trials).**

## Stage 4 — `curr_syn6_s4_walk_nocontact` (walking kine in air, contact OFF) — DONE, GENUINE AIR-GAIT (imperfect pattern)

- **Run**: `_curriculum_syn6.py 4 22` detached 18:03:19–18:54:57; exit code 0
  (`logs\curr_syn6_s4.done` first line = 0); log `logs\curr_syn6_s4.log`
  (11,904 bytes). RESUMED the study holding the wiring session's 1-trial
  route-proof smoke (trial 0 = seed defaults, −254.936) → **23/23 trials
  COMPLETE** (sqlite census: `tmp\harvest_check_syn6_s4.py` output); no new
  seed enqueued (study non-empty, `_curriculum_syn6.py:273`).
- **Score landscape**: 10/23 real kine scores, 13 frozen (−320, no cycles —
  sentinel-heavy walking landscape as in the stock chain), 0 NaN. Winner
  −237.253 (trial 14) beats the smoke seed −254.936 (trial 0) by 17.7 —
  above both the −315 clip and −320 frozen sentinel.
- **Winner** (`curriculum_syn6_stage4.json`): score **−237.253**, trial 14,
  params drive=2.7864, rg_nap_h=0.3526, desc_e=1.4796, desc_f=1.4494,
  rg_to_pf=2.0250, heel_rge=0.1231, toe_rge=0.1018, rg_mutual_inh=4.5943,
  **rg_weak_exc=0.1769** (>0 — the conditional RG↔RG weak-excitation
  topology is ENGAGED in the winner), pf_recip_inh=3.8038.
- **Validity check: PASS (GENUINE AIR-GAIT PATTERN)** — winner replayed
  bit-exactly (`tmp\replay_syn6_s4_winner.py` →
  `logs\curr_syn6_s4_exploit_check.log`): nan=False, kine detected
  (bilateral=True, n_cycles 3/leg, period_cv 0.0004); recomputed
  max(kine_score,−315)−tilt-penalty = **−237.2534 == json −237.2534 (delta
  +0.0000)**; RG_E_r rises=3 (t>5, span 10.1 mV), knee −49.6..+10.1°.
- **One-line diagnosis**: a real, very regular but SLOW in-air gait —
  T 2.45 s/cycle vs ref 1.23 (0.41 Hz), legs IN-PHASE (lag_rl 0.0 vs ref
  0.51), duty 0.41 vs 0.61, ankle PF-biased (mean −53.5° vs ref +4.2°),
  hip RMSE ~24°, knee_min −48.8 vs ref −69.7 — i.e. a genuine cyclic
  pattern with known architectural gaps (phase/cadence/ankle), NOT a
  static pose. Stage 5 starts contact tuning from this config; the
  in-phase-leg gap is the first thing contact feedback must attack.
- **Protection**: protected `spinal_run.npz` untouched (re-verified:
  8,238,797 bytes, mtime 2026-09-25 01:10:02 AM); chain-private
  `spinal_run_syn6.npz` (8.4 MB, 18:54:57) is the real save target via
  AARL_NPZ (runner.py:1629); replay scratch npz deleted (Test-Path False).

**PROCEED → stage 5 (`curr_syn6_s5_walk_contact`, ground walking, 22
trials).**

## Stage 5 — `curr_syn6_s5_walk_contact` (ground walking, contact ON) — DONE, GENUINE GROUND WALK; SEED WON — CHAIN COMPLETE

- **Run**: `_curriculum_syn6.py 5 22` detached 19:03:20–19:55:25; exit code 0
  (`logs\curr_syn6_s5.done` first line = 0); log `logs\curr_syn6_s5.log`
  (12,929 bytes). Fresh study → **22/22 trials COMPLETE** (sqlite census:
  `tmp\harvest_check_syn6_s5.py` output; trial 0 = the enqueued seed, 21
  searched). Seed = ALL 10 stage-4 winner keys + contact_onset/contra_swing/
  ib_rge at 0.0 + syn6 1.0 (log "seeded" line — full 13-key dict).
- **Score landscape**: 8/22 real kine scores, 14 frozen (−320, no cycles),
  0 NaN. **The SEED (trial 0, −197.215) won outright**; best searched trial
  only −261.95 (trial 17). A "seed wins" outcome: TPE could not improve the
  transplanted stage-4 winner on the ground in 21 trials, and the three
  stage-5 contact pathways (contact_onset, contra_swing, ib_rge) never
  earned nonzero values — the same "new pathways tolerated, not additive"
  pattern as stages 2 (heel/toe) and the stock s2b.
- **Winner** (`curriculum_syn6_stage5.json`): score **−197.215**, trial 0,
  params = the stage-4 winner verbatim + contact_onset=0.0, contra_swing=0.0,
  ib_rge=0.0 (13 keys).
- **Validity check: PASS (GENUINE GROUND WALK)** — winner replayed
  bit-exactly (`tmp\replay_syn6_s5_winner.py` →
  `logs\curr_syn6_s5_exploit_check.log`): nan=False, cfg=ground,
  **contact_frac_r 0.955 / contact_frac_l 1.0** (real foot loading; the
  stage-4 air walk had 0.0), right leg n_cycles=4 at T 1.88 s (period
  regular), lag_rl 0.395 (ref 0.51), kz 0.840 and tilt_max 28.07° both
  clear of the −20/−10 penalties; recomputed **−197.2146 == json −197.2146
  (delta +0.0000)**; RG_E_r rises=3 (span 10.1), knee −54.5..+10.0°.
- **One-line diagnosis**: a real right-leg-driven ground walk — but the
  left leg is planted ~100% of the time (no clean left cycle detected,
  bilateral=False), double-support is heavy (duty 0.745 vs ref 0.61), knee
  flexion shallower than ref (−49.7 vs −69.7). Score −197.2 lands in the
  same band the s3k STOCK walker scores against real gait references
  (−189.5…−241, goal-5 reports) — a respectable chain endpoint, with
  left-leg participation the obvious next architecture problem (same
  family as the stock chain's interleg/phase gaps).
- **Protection**: protected `spinal_run.npz` untouched (re-verified:
  8,238,797 bytes, mtime 2026-09-25 01:10:02 AM); chain-private
  `spinal_run_syn6.npz` (8.5 MB, 19:55:24) is the real save target via
  AARL_NPZ (runner.py:1629); replay scratch npz deleted (Test-Path False).

**CHAIN COMPLETE (5/5).** Final deliverables: winner jsons
`curriculum_syn6_stage{1..5}.json` (spinal cwd), studies
`curr_syn6_s1_air_deaff` (20 trials), `curr_syn6_s2_air_aff` (18),
`curr_syn6_s3_balance` (18), `curr_syn6_s4_walk_nocontact` (23), 
`curr_syn6_s5_walk_contact` (22) in `spinal\optuna_syn6.db`; validity
replays in `logs\curr_syn6_s{1..5}_exploit_check.log`. The syn6 production
walker = stage-5 winner params (replayable via the stage-5 command line).




