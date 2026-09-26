# Goal 3 — Data-driven test plan: how to TEST a robust balance + walking model

**Date:** 2026-09-25 · **Machine:** EB475WS4 · **Status:** plan (no protected file touched).
**Inputs read this session:** `reports_20260924\balance_data_catalog.md` (full),
`reports_20260923\goal5_gait_library_catalog.md` (full), `gait_refs\` (43 npz, `dir` + new
survey script `reports_20260925\tmp\ref_survey.py`, output `tmp\ref_survey_out.json`),
`gait_lib_score_all.py`, `gait_lib_loader.py`, `kine_ref.py` (key defs), `basin_gate.py`,
`wean_rig.py`, `runner.py` arg-parsing (manual argv, `runner.py:708-936`), and
`reports_20260923\goal5_allrefs_scoring.md` (existing s3k baseline). Staging dir
`D:\temp\gait_lib_staging\downloads\` listed read-only (17 items incl. `RawEMGData-latest`
→ `RAW_EMG_DATA`, 8 `Subject##-latest` zips, results-deficits, RunningSimulation,
Hamner2010, assistloadwalk, scone-setup-files-v2).

**Purpose (Ben's ask):** a concrete battery that uses the balance-study and gait data to
TEST — not train — a walker, so fragility is caught **before** any BPA bracket is 3D-printed.

---

## 0. Priority order (what catches a fragile walker first)

| P | test | cost | catches |
|---|---|---|---|
| **P1** | T1 multi-ref kine scoring + T2 Ong speed sweep | ~1 day tooling, minutes/eval | single-reference overfit (the s3k knee/hip misfit pattern), no speed headroom |
| **P2** | T4b in-sim perturbation battery (push / contact / strength derating) | ~1-2 days | basin-marginal tuning, zero torque margin before hardware exists |
| **P3** | T7 basin gate (exists) + T8 rig weaning (exists) on every finalist | minutes-hours, zero new code | bifurcation-edge winners that will not survive Simulink/co-sim or a real floor |
| **P4** | T4a Wang & van den Bogert standing-balance battery | download (264 MB processed only) + ~2 days | balance stage has no human reference at all today |
| **P5** | T5 synergy-timing vs raw EMG; T6 leave-one-family-out | ~2-3 days / compute-heavy | architecture-level (synergy structure), training-set leakage |

P1-P3 need **zero downloads** and are the pre-print gate. P4-P5 deepen the evidence for the
dissertation and the balance stage.

---

## 1. What each data family is good FOR (testing terms)

Verified family stats (this session, `ref_survey_out.json`):

| family | n | measured range (duty_r / T_r / knee_min_r) | testing role |
|---|---|---|---|
| `subject01` of-record cycle (kine_ref built-in, not an npz) | 1 | 0.61 / 1.23 s / −69.7° | **in-distribution regression anchor** — the one ref the current tuning was fitted against; any new model must reproduce its historical score (s3k: −189.5) bit-comparably |
| Falisse 2022 `Case_40_motion` (predicted walking) | 1 | 0.57 / 1.11 s / −59.2° | **in-distribution but model-free** — a *predictive* (non-tracking) walker, methodologically closest to what we build; score here says "matches a good optimizer gait", not "matches the trial we tuned to" |
| Ong predictive walking `ong_speed_050..200` + 3 selfsel | 10 | duty 0.53-0.62, T 2.10→1.02 s, knee_min −66.8..−73.3°, ds 0.25→0.08 (verified sweep) | **speed-generalization curve** — the ONLY local family spanning speeds; duty/ds/cadence all move monotonically with speed, so a scalar-tuned walker that cannot re-time will fail the far end (T 2.10 s = 0.48 Hz vs current sim 0.90 s ≈ 1.1 Hz — a predicted, quantifiable failure) |
| Arnold `subject{01,02,04,08,10,11,17,20}_Run_*` | 32 | duty 0.31-0.45, T 0.56-0.80 s, knee_min −80..−140° | **out-of-distribution stress only** — running is a different gait; quantify degradation, NEVER tune to it (explicitly excluded from curricula and LOFO training folds) |
| Wang & van den Bogert 2020 (Zenodo 3819630, CC BY 4.0) | 8 subj × 4 × 5-min | not local yet | **perturbation robustness reference** — random-pulse standing balance with markers + 6-DOF GRF + 9 EMG + **processed joint angles/torques** (no IK needed); feeds the stage-3 standing-balance gate (currently zero human references exist for it) |
| Hamner nmbl `RAW_EMG_DATA` (staged, read-only) | 10 subj × 4 trials | not integrated | **EMG/synergy architecture validation** — timing of muscle activations, the dimension `gait_refs` npzs lack entirely (no EMG inside any ref, verified) |
| results-deficits (weakness/contracture, staged) | 2 zips | not integrated | later: pathological-gait degradation references (pairs with the strength-derating test T4b-3) |

Honest structural limits of the local library, both verified: **no GRF waveforms and no EMG
are stored in any npz** (balance catalog §1), and there is **one reference cycle per ref**
(single cycle, phase-normalized) — so per-ref variance is unknowable and all scores are
point comparisons, not envelope comparisons.

---

## 2. The test battery

Notation: kine_score is **higher (less negative) = better, 0 = perfect** (`kine_ref.py:30,314`).
Note the report string at `gait_lib_score_all.py:106` says "lower kine_score = better" — the
opposite of kine_ref; ignore it (docstring bug); the recorded ranking confirms
higher-is-better (s3k subject01 −189.5 beats Falisse −200.5).

### T1 — Multi-reference kine scoring (43 refs)

- **How:** parameterize `gait_lib_score_all.py` (currently hardcoded to `spinal_run_s3k.npz`,
  `gait_lib_score_all.py:20-21`): take `(run_npz, out_prefix)` argv, keep the per-ref
  REF_CACHE swap (its `load_npz_ref` + `KR.compare(t,q,neuro,WALK_START,ref=ref,contact=...)`
  loop, lines 56-71). Output one CSV per run: ref, kine_score, **decomposed** terms
  (rmse_hip/knee/ankle, range/phase/period/duty penalties, guards) — the decomposition exists
  in the `compare()` output dict (`kine_ref.py:236-295`) and is needed to tell shape
  degradation from guard penalties (see T3).
- **Exclude from headline:** `ong_selfsel_Init200` (degenerate left cycle — flagged in the
  campaign brief; report it in a footnote row).
- **Acceptance (walking refs = subject01 + Falisse + 10 Ong; 12 refs):**
  - no-regression for any future candidate vs s3k: min over the 12 ≥ −235
    (s3k's measured worst is −224, Ong family; ~5% margin);
  - candidate *improves*: median(12 refs) > s3k median (compute once T1 runs; expected ≈ −212
    from the recorded −189.5..−224 range);
  - hard floor: no ref < −300 and no NaN/`frozen_*` flag on any walking ref.
- **Effort:** 2-4 h (tooling mostly exists).

### T2 — Ong 0.5-2.0 m/s speed-generalization curve

- **How:** same T1 output restricted to `ong_speed_*`, plotted vs the refs' own T (0.50→2.00
  m/s maps to T 2.10→1.02 s, verified). Two curves: (a) kine_score vs speed, (b)
  `|T_sim − T_ref|/T_ref` vs speed (period tracking = the DRIVE/cadence lever).
- **Acceptance:** (a) score-vs-speed curve is monotone-declining (any inversion means the
  model is gait-specific, not speed-specific); (b) period error ≤ 25% per ref point
  (currently T_sim 0.90 s fails every Ong ref except speed_200 by this bar — that is the
  point of the test); (c) total spread across 0.5-1.25 m/s ≤ 60 points.
- **Generalization claim for the dissertation:** model passes T2 ⇒ "re-times across the human
  walking speed range" — the single most valuable claim the Ong family buys.
- **Effort:** 2 h on top of T1.

### T3 — Arnold running OOD stress (quantify, don't tune)

- **How:** from T1's CSV, report median running score vs median walking score as the
  **degradation index**; decompose each running score into shape terms (rmse/phase) vs
  guard/penalty terms (foot-never-loads +20, no-cycles +12/+8, `kine_ref.py:217-226`) so the
  number is interpretable — at −332 (sprint, recorded) the score may be penalty-dominated,
  which says "not a runner", not "bad shapes".
- **Rule:** running refs never enter a curriculum, an objective, or a LOFO training fold.
  Report once per model generation.
- **Acceptance (printing gate):** degradation must exist (it will) but the model must still
  *complete* its walk schedule after exposure — no NaN, stays up. There is no score bar;
  running is stress, not target.
- **Effort:** free once T1 exists; interpretation 1 h.

### T4 — Perturbation battery (two halves)

**T4a — Wang & van den Bogert standing-balance references (download required).**

- **What to download:** Zenodo record 3819630, **Processed bundle only (264 MB)** — it ships
  processed joint angles + torques (inverse-dynamics, per the record description quoted in
  balance catalog §2) so no IK pipeline is needed; the 2.9 GB raw is optional (grab later
  only if CoP traces at native rate are required). Includes the MATLAB processing code +
  experiment PDF. 8 subjects × 4 trials (2 quiet + 2 perturbed), 5 min each.
- **Unverified at download time (catalog §2):** perturbation direction (treadmill-mounted
  device); confirm from the data files before building pulse-matched scenarios.
- **Metrics (sim vs human, per 60-s window):** 95% sway area (ellipse fit to CoM x-y path),
  sway velocity (path length/time), max excursion per pulse, recovery time back into ±2σ of
  the subject's own quiet band, fall = `com_z < 0.62 m` (same threshold as
  `wean_rig.py:31`) or schedule termination. Human-side thresholds are computed from their
  quiet-trial percentiles at integration time — do not guess numbers now.
- **Sim protocol:** stage-3 standing configuration (`--stand-eval` exists, `runner.py:936`;
  VEST knobs exist default-0), pulse train seeded to match Wang's pulse statistics once read.
- **Acceptance (provisional until calibrated):** quiet sway within 2× human median sway
  area; recovery from ≥80% of matched pulses; 0 falls at the human pulse amplitude.

**T4b — In-sim perturbation battery (no download; catches fragility NOW).**

New flags on runner.py, all **default-neutral** (protected default build must stay bit-exact:
standing gate −160.23425729850192, AGENTS 2026-09-24; rerun `_pf_layer_variant_test.py` +
`(410,376,1186)` counts after the edit):

1. **Push:** `--push "t,dur,Fx,Fy"` repeatable → `data.xfrc_applied` on the pelvis; grid:
   F ∈ {20, 40, 60} N × dur 150 ms × phase ∈ {mid-stance L, mid-stance R, double-support},
   9 cells × 3 seeds = 27 runs. Metrics: recovered-within-2-cycles / fell / completed;
   post-push cycle kine_score vs pre-push.
2. **Contact/stiffness jitter:** do NOT touch solref (direct k,b was proven UNSTABLE at
   2 ms, AGENTS 2026-09-23); instead (a) existing default-off `--contact-damp
   lessviscous|nonlinear` + `--adaptive-tol` variants, (b) floor friction ±20% and solimp
   jitter in a scratch copy of the XML (pattern: `make_simbridge_xml.py` generates patched
   copies; cvt3.xml itself is never touched), (c) the measured 2 ms→0.5 ms timestep ladder
   (`goal3_timestep_ladder.py`: 0.44° RMS ground) as the integration-error axis.
3. **Strength derating:** `--fmax-scale S` scaling `actuator_gainprm[:,2]` (= MuJoCo muscle
   Fmax, AGENTS) in a scratch XML, S ∈ {0.9, 0.8, 0.7, 0.6}. **Falsifiable prediction from
   the existing torque budget (hip 3.4×, knee 2.8×, trunk 2.1×, ankle 1.4× — AGENTS
   2026-09-13): ankle failure first, near S ≈ 0.7.** If hip or trunk fails first, the torque
   budget or the model is wrong — that alone is worth the test. Pairs with the staged
   results-deficits data later.
- **Acceptance:** recover ≥80% of 20 N pushes; complete walk schedule under friction ±20%
  and contact-damp variants; strength margin: completes at S = 0.8, degradation curve
  reported to S = 0.6.
- **Effort:** 1-2 days (runner flags + battery driver + seeded report).

### T5 — Synergy-timing validation vs raw EMG

- **How:** Hamner `RAW_EMG_DATA` (staged; 10 subj × 4 trials) → filter/high-pass, rectify,
  envelope (standard, matches their pipeline docs), map 11 lower-limb muscles → our 43-muscle
  leg set via `muscle_map` (`bsolve_ik.py:46`); compute per-synergy activation timing
  C = NNLS(W, A) as in `gait_lib_loader.py:186-197`, compare **phase of synergy peaks**
  (cycle-normalized) vs the sim's W_PF_MN-weighted group profiles.
- **Known-negative control (verified):** our subject01-derived six-synergy basis does NOT
  transfer across datasets — vs Falisse activations R²-vs-0 0.60 but VAF ≈ −0.02
  (`gait_lib_loader.py:169-202`, catalog §4). So: report the current basis's VAF on Hamner
  as a measurement; the **pass bar applies only to a re-trained multi-reference basis**
  (train on subject01+Falisse+Ong activations where available): held-out VAF ≥ 0.6 and
  peak-timing error ≤ 15% of cycle for the 3 dominant synergies.
- **Effort:** 2-3 days (EMG preprocessing dominates).

### T6 — Leave-one-reference-family-out (training-hygiene)

- **Folds:** walking folds only — {subject01}, {Falisse}, {Ong×9 (Init200 excluded)}; running
  always held out. A candidate trained on folds A,B is scored on held-out C with T1's CSV.
- **Acceptance:** held-out family median within 25 points of the in-family median (current
  in-sample spread is already 35 points: −189.5 vs −224 — so 25 points is ambitious but
  meaningful). LOFO cost is a re-*training* per fold (~hours-trial optuna); for existing
  winners it is report-only (they were s3k-tuned on subject01+library mixes).
- **Effort:** 1 day tooling; compute-heavy when actually run per fold.

### T7 — Basin gate (exists) — mandatory for every finalist

- `basin_gate.py` on the finalist's RUNNER_DUMP_STATE json, defaults (σ=1%, 12 trials, 20 s;
  PASS = last-window swing ≥ max(0.1 mV, ½ baseline first window), `basin_gate.py:22-23,162`).
- **Acceptance:** 12/12 pass, 0 exploded, margin ≥ 0.5×. History says this matters: v10-best
  passed 12/12 while the study config did not self-sustain at all (AGENTS 2026-09-16).
  **Known caveat (keep):** single borderline trials can be BLAS summation-order noise; one
  exploded trial → rerun seeded before declaring failure.
- **Effort:** 0 (exists); ~10-20 min per finalist.

### T8 — Rig-weaning ladder (exists) — mandatory for every finalist

- `wean_rig.py` (its PASS rule: full schedule, `kz > 0.62`, tilt < 35°, `wean_rig.py:49`).
- **Acceptance:** S = 0.8 minimum for any "walks" claim (the current known support boundary);
  S ≤ 0.6 pass = pre-print target; report the full ladder, not just the deepest stage.
- **Effort:** 0 (exists); hours per ladder (each stage = full stand→walk→stand).

---

## 3. Extend vs NEW (scripts, effort, gates)

| script | action | effort | notes |
|---|---|---|---|
| `gait_lib_score_all.py` | **extend** → argv-parameterized, per-term CSV columns, JSON summary | 2-4 h | keep read-only behavior for s3k reproduction (−189.5 anchor must still reproduce) |
| T2 plotting | **new** `reports_20260925\speed_sweep.py` (reads T1 CSV) | 2 h | |
| `runner.py` | **extend, default-neutral only**: `--push`, `--fmax-scale`, nothing else | 1 day incl. regression gate | every flag default = current bit-exact behavior; gate = standing −160.23425729850192 + (410,376,1186) + JSON-rule if any new G key ever appears |
| `reports_20260925\push_battery.py` | **new** driver (seeds, grid, report md+json) | 0.5-1 day | calls runner via its existing `main(argv)` API (`runner.py:708`, pattern proven by `wean_rig.py:47`) |
| `basin_gate.py`, `wean_rig.py` | **no change** — run as-is | 0 | |
| `kine_ref.py` | **no change** | 0 | fix nothing inside; the score_all report-string bug (:106) is corrected in the new caller's output header |
| Wang integration | **new** `wang_balance_loader.py` + `standing_battery.py` | 2-3 days + download | after stage-3 standing config is declared current |
| EMG pipeline | **new** `emg_timing.py` | 2-3 days | Hamner raw EMG is staged read-only; process into `reports_20260925\` derivatives only |

Compute reality check: EB475WS4 is not a heavy-runs box (AGENTS); 27-cell push grids and LOFO
retraining belong on easteregg2 or run overnight here capped.

---

## 4. Plugging into the curricula (which metric gates which stage)

Stage numbering per the 2026-09-24 renumber (AGENTS): **3 = standing balance** (curr_s3_balance),
**4 = ground walk** (curr_s4_nocross), **5 = per-PF-layer contact variant** (curr_s5_pfvariant,
stopped, resumable). Gates, in order a candidate meets them:

1. **Stage-3 exit gate (standing):** T4a quiet-sway + pulse-recovery metrics vs Wang —
   replaces today's purely sim-internal stand-eval with a human reference. Not gatable until
   the Wang download lands; until then stage 3 exits on the existing `--stand-eval` +
   rig/VEST tuning only (honest gap).
2. **Stage-4 exit gate (walking):** T1 acceptance (min ≥ −235, median beats s3k, no frozen
   legs) **plus** T2 monotone + period-error ≤ 25% **plus** T7 12/12 **plus** T8 S ≥ 0.8.
   T3 is reported, never gated. This 4-part gate is specifically designed to prevent a
   repeat of the 2026-09-20 static-pose exploit (a winner that scored well on one metric
   while not walking): kine-only winners fail T2's period-error term and T4b pushes.
3. **Stage-5 acceptance (variant):** the five contact-knob winners must beat s3k **on the
   full battery**, not on kine alone — same gates as stage 4, evaluated on the variant's
   winner. Any variant that wins kine but fails pushes/basin is a basin-marginal artifact.
4. **Finalist (pre-3D-print) gate = T1 ✓, T2 ✓, T7 ✓, T8(S≤0.6) ✓, T4b acceptance ✓.**
   T4a/T5/T6 strengthen the dissertation claim; they are not print blockers.

**Known gaps (stated plainly):**
- **No M/L GRF anywhere**; loader reads vertical only (`gait_lib_loader.py:27-30`) and npz
  refs store no GRF at all → all scoring is joint-kine + duty/T/lag/ds scalars. GaitRec /
  Gutenberg (both CC0, catalog §3) are the later fix; sim-side per-foot M/L contact-force
  extraction has its own work item (catalog §6).
- **One cycle per ref** → no human variance envelopes; a mean±SD envelope framing (as in the
  Dorn 2015 methodology, gait-library catalog §2) needs Schreiber 2019 / WBDS-style
  multi-trial data (catalog §3).
- **Wang direction unverified**; **ong_selfsel_Init200 degenerate left leg** (exclude);
  **kine_score guards** can dominate OOD scores (decompose before interpreting).
- **fitted-winner fragility:** v10-era `--fitted --bestN` winners no longer walk under
  current physics (AGENTS 2026-09-23) — the battery must always run against the *current*
  champion (s3k family), never against archived recorded numbers, except as historical
  anchors.

## 5. Licensing + provenance

| data | license (verified live, balance catalog §7 / gait catalog) | note |
|---|---|---|
| Wang & van den Bogert 2020, Zenodo 3819630 | CC BY 4.0 | cite Wang & van den Bogert; attribution in dissertation methods |
| Moore 2015 perturbed walking, Zenodo 13030 | CC0 | future walking-perturbation extension (AP/belt only) |
| GaitRec, Gutenberg | CC0 (all items) | M/L GRF norms when the loader gap is fixed |
| Schreiber & Moissenet 2019 | CC BY 4.0 | experimental multi-speed + EMG + 3D GRF companion |
| Hamner nmbl raw EMG (staged) | BSD-style Stanford 2013/2014 (per gait catalog §3) | SimTK login-held source; staged local copy is Ben's download |
| Arnold muscfib data | CC BY-NC 3.0 | **non-commercial flag** before any patent-adjacent reuse (gait catalog "Read this first") |
| Falisse predictsim_mtp | per its GitHub repo (verify at next use) | staged copy is outside the repo by design |

Provenance rules: downloads stay outside the repo (`D:\temp\gait_lib_staging\`, size-reduction
standing objective); only derived npz/CSV/refs are committed; every derived number in a
report cites its generating script + log under `reports_20260925\logs\`; nothing in this plan
modifies protected default behavior (runner/params/build defaults, `spinal_run.npz`,
existing optuna studies, reports_20260923/24).

---

## 6. What was and wasn't executed this session

- **Executed (read/verified):** catalogs read; `dir` of `gait_refs` (43 npz); new read-only
  survey `reports_20260925\tmp\ref_survey.py` → `tmp\ref_survey_out.json` (per-family
  duty/T/knee_min/hip_range/ds — numbers quoted in §1 come from it); `dir` of
  `D:\temp\gait_lib_staging{,\downloads}`; code reads cited by path:line above;
  existing s3k all-refs scoring report read.
- **NOT executed:** no walker runs, no battery runs, no downloads, no scoring rerun — this
  ask is the plan; every acceptance threshold above is either anchored to already-recorded
  results (s3k −189.5/−200.5/−224; −160.23425729850192 standing gate; basin 12/12 history;
  torque budget ratios) or explicitly marked provisional-pending-calibration (T4a sway bars).
