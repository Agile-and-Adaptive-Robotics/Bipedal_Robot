# Goal 4 — Optuna tuning review vs SCONE practice + Ben's ensemble question

**Date:** 2026-09-23 (EB475WS4). **Author:** optimization-review subagent (ZCode).
**Inputs read:** `Code/MuJoCo_SNS/spinal/{optuna_walk*.py, _curriculum.py, kine_ref.py, runner.py, basin_gate.py, optuna_walk.db}` (read-only; the db was **copied to %TEMP% and analyzed on the copy**), the local SCONE 2.4.4 install + writable tutorial copies, scone.software docs, and the papers cited inline.
**Companion files:** `goal4_db_analysis.json` (raw per-study output of the analysis I ran), `ensemble_objective_scaffold.py` (design sketch, smoke-tested, NOT wired into anything).

---

## 0. What changed / what Ben must decide

Nothing in the study stack was edited (per instructions). Three decisions for Ben:

1. **Adopt the ensemble objective?** My recommendation: **yes for random-seed/chaos robustness and small plant jitter (N=3), no for synthetic IK-reference perturbation** — the reference axis has a data problem, not a robustness problem (section 1). Cost is ~3 min/trial instead of ~1 (measured, section 4).
2. **Switch the next study to `constraints_func` + multivariate TPE + a ~10-param space.** The evidence: 22–70% of every stage-3 study's budget landed on a single sentinel value, and the curriculum search space grew to 32 params at 1.3–2.5 trials/param (section 3).
3. **Whether to add a runner-side Fmax/damping jitter hook** (the one piece the ensemble scaffold needs that doesn't exist yet). It must follow the default-off = bit-identical contract like every conditional knob.

---

## 1. DIRECT ANSWER — "slight variations in model parameters to get slightly different IK results, and use them all to train a solution?"

**Nuanced yes — but the three things bundled in that sentence are different levers with different evidence behind them, and one of them (IK perturbation) I recommend against as stated.**

### (a) Perturbing PLANT parameters (Fmax, damping, rig stiffness) — YES, small and inside the objective

This is standard neuromuscular-optimization practice, and SCONE does it in three ways simultaneously rather than as a separate ensemble:

- **Motor noise during optimization**: `NoiseController { base_noise = 0.02 proportional_noise = 0.15 random_seed = 0 }` inside the *optimized* controller stack — `C:\Users\Ben Bolen\Documents\SCONE\Tutorials\Tutorial 3b - Motor Noise Balance - OpenSim.scone:23-27`.
- **External perturbations during optimization**: backward and forward 100 N torso pushes every 4 s — `Tutorial 4c - Perturbed Gait - OpenSim.scone:26-45`.
- **Sampled initial state per candidate**: `initial_state_offset = 0~0.01<-0.5,0.5>` (a `~` distribution = an optimization parameter that CMA-ES samples) — same file, line 14, and every gait example.
- **Explicit plant variants**: `Gait2D - GeyerHerr2010 - Scaled - Hyfydy.scone:21-34` scales limb segments, weakens vasti/gastroc/soleus 10% (`MuscleModifier { max_isometric_force { factor = 0.9 } }`), and re-masses the model (`ModelModifier { mass = 65 }`) — then optimizes the *asymmetric* controller on that variant.

Literature: Koelewijn & van den Bogert 2022 (*PeerJ* 10:13085, listed on SCONE's publications page) show optimization **under uncertainty** finds co-contraction solutions that nominal optimization does not — i.e., perturbation during training doesn't just protect the solution, it *moves the optimum*. In sim-to-real robotics this is domain/dynamics randomization (Peng et al., ICRA 2018, "Sim-to-Real Transfer of Robotic Control with Dynamics Randomization"; Tobin et al., IROS 2017). For us the plant-side argument is even stronger than usual because our "plant" is a converted model with known fidelity gaps (rigid-tendon simplified Thelen, converter RoM normalization — AGENTS.md MuJoCo notes): perturbing Fmax/damping is a proxy for exactly that epistemic uncertainty. **Inference (mine, not a citation):** keep jitter small (±5–10%, not ±30%) — Koelewijn's result implies large uncertainty pressures toward stiff/co-contracted gaits, which we'd pay for in tracking score and (on hardware) air consumption.

We already have the seed of this: `basin_gate.py` perturbs every scalar param ±1% and requires the rhythm to persist 20 s (`basin_gate.py:8,115-118`; v10-best passed 12/12 per AGENTS.md). The SCONE lesson is to move that pressure *from a post-hoc gate into the training score* — or at least keep the gate but stop treating it as optional.

### (b) Perturbing the IK / reference data — NO as stated; the problem is we have too little reference data, not that it's over-precise

Measured facts (scripts in %TEMP%, outputs in section 6):

- `kine_ref.py:94-99` builds the reference from **one gait cycle per leg** (`t0, t1 = onsets[0], onsets[1]`) of `subject01_walk1_ik.mot`.
- That file spans only **0.50–2.50 s (121 samples)**; the GRF file yields **2 right onsets and 3 left onsets** — i.e., the right leg has exactly **one** cycle total, the left has two. There is nothing else to average.
- `subject01_walk1.mot` (the 15-s, 901-sample file that looks like it holds more) is **all zeros** in every coordinate column — a placeholder, not data.
- The repo holds only walk1 + static + the Tutorial1 crouch/normal motions (`dir /s /b *.mot *.trc` on `Solid_Models\OpenSim\Gait2392_Robotbody`).

So "slightly different IK results … use them all" cannot be fed by data we currently have. The two honest options: (i) **process more trials/subjects** (IK on additional walking trials — the raw data is not in this repo; this is a Ben decision about sourcing), or (ii) **synthetic reference jitter** (bootstrap the one cycle). I recommend (i) if more trials can be obtained, and *against* (ii) as a training signal: jittering a single cycle mostly regularizes against cycle-specific quirks of that one cycle (e.g., our reference `T_r = 1.233 s` exactly, `knee_min_r = −69.7°` — from my run of `kine_ref.ref_cached()`), which is a **bias** in the target, not noise around a true target. SCONE's analog, `MimicMeasure` (tracking a .sto), is explicitly "generally not recommended" for predictive sims (skill quick-reference, from the SCONE docs measures page) — the field's instinct is to avoid over-fitting one recorded motion. Where reference *averaging* is available (2nd left cycle), using it is free and worthwhile: `kine_ref.compare(..., ref=...)` already accepts a per-call reference (`kine_ref.py:193`).

### (c) Multiple random seeds — YES, and for us this is partly just measurement error control

Our simulator is chaotic with documented run-to-run nondeterminism even at identical seeds (BLAS summation order; AGENTS.md optimizer-insight (6) and the basin_gate caveat: an identical-seed trial exploded 1e164 in a 2 s smoke run and was stable in the 20 s run). Every optuna trial is **one 16-s rollout scored once** (`optuna_walk.py:114`, `optuna_walk_v8b.py:177`, `optuna_walk_v10.py:179`) — so the score TPE conditions on is a noisy observation of a chaotic system. SCONE hits the same issue (FAQ "Why do I get different results on different machines?" — floating-point accumulation) and its answer at optimizer level is **multi-start**: `CmaPoolOptimizer` "runs multiple CMA-ES optimizations in a prioritized fashion, based on their predicted fitness" (SCONE reference manual page; entry also in the local help keyword list `C:\Program Files\SCONE\resources\help\keywords.txt`). For us the cheapest version: N rollouts per trial with different perturbation draws (score = mean − λ·std) — this simultaneously (i) reduces score noise, (ii) selects basin-robust parameter sets, and (iii) is exactly Ben's ensemble proposal applied to the axis where it's cheapest.

**One-line answer:** yes to ensembling over *seeds* and *small plant perturbations* (SCONE-style, inside the objective); for *IK variants*, first get more reference data — perturbing our single-cycle reference synthetically would train robustness to a biased target.

---

## 2. What SCONE does that we do not — and vice versa

SCONE = CMA-ES-based shooting optimization (Geijtenbeek 2019, JOSS 4(38):1421, doi 10.21105/joss.01421; FAQ: "different flavors of Covariance Matrix Adaptation [Hansen 2006]"). Optimizer family in the reference manual: `CmaOptimizer`, `CmaOptimizerSpot`, `CmaPoolOptimizer`, `EsOptimizer` (fetched from scone.software/doku.php?id=reference; `cma_optimizer/cma_pool_optimizer/es_optimizer` in the local `resources\help\keywords.txt`).

| SCONE practice | Evidence | Do we have it? |
|---|---|---|
| **Robustness inside the objective** (motor noise, perturbation pushes, sampled initial state) | Tutorials 3b/4c (local .scone files, cited in §1a) | No — we gate *after* (basin_gate) or not at all during tuning |
| **Population optimizer with covariance adaptation** (learns parameter correlations; the drive↔pf_gain ridge we know exists) | FAQ (Hansen 2006); `CmaOptimizer { lambda mu sigma }` in every tutorial | Partially — our TPE is **univariate**: `TPESampler(seed=…, n_startup_trials=…)` with multivariate left at default None (verified in installed optuna 5.0.0; loaders at `optuna_walk.py:139`, `optuna_walk_v10.py:224`, `_curriculum.py:289-290`) |
| **Multi-start pool, budget prioritized by predicted fitness** | `CmaPoolOptimizer` doc page | No — our studies are single-chain with one enqueued seed (winner of previous stage) |
| **Stall detection** (`min_progress` + `window_size`) | skill §Optimizer + `min_progress = 1e-4` in every tutorial .scone | No — fixed `n_trials` (60/80/100/150) |
| **Warm-start as a distribution, not a point**: `.par` stores best/mean/**stdev**; `init { file std_factor = 10 }` re-widens | skill §Optimizer; `Uneven Terrain` example line 5 | No — `enqueue_trial(seed)` plants a point and TPE re-explores from scratch |
| **Tiered/lexicographic objectives** (`threshold`, `soft_threshold`, `use_first_non_zero_result` on CompositeMeasure; e.g. `Gait10.scone`: weight 100, threshold 0.05, termination_height 0.85) | local `Tutorials\measures\Gait10.scone`; skill §Measures | No — flat score + hand-tuned sentinel constants, which repeatedly collapsed (v4's −25 plateau; v8b's −65 collapse; s3c's −320 sweep — AGENTS.md + db numbers below) |
| **Population-parallel evaluation** (~15 threads/optimization) | FAQ hardware answer | No — serial `study.optimize` loop (1 eval ≈ 1 min, measured below) |

What **we** do that SCONE doesn't: staged curriculum with **conditional topology** (new pathways exist only when gain > 0, byte-identical at 0 — v5 notes); **bit-exact regression contract** (`repr(drive)`, seeded binaries); **post-hoc fANOVA param importance** (optuna) — SCONE has nothing built-in (sensitivity analysis is done in separate frameworks, e.g. Buchmann & Renjewski 2024 BioRob, SCONE publications page); **per-trial CSV + automatic full-22-s capture of new global bests**; and reference-tracking to human IK (SCONE's predictive style prefers cost+constraints; its tracking measure is discouraged). Our curriculum chaining (stage winner seeds the next stage's enqueue) is also exactly SCONE's warm-start pattern, just point-valued.

---

## 3. Concrete optuna config for the next study

Environment measured: **optuna 5.0.0** in the myo env (`python -c` via temp script → `optuna 5.0.0`). Sampler defaults verified by inspection: `multivariate=None`, `group=False`, **`constant_liar=True`** (changed from False in older optuna — our loaders already benefit), `n_startup_trials=10`, `n_ei_candidates=24`.

### 3.1 Replace sentinel scores with `constraints_func` — the single biggest win

Measured sentinel-plateau budget (fraction of COMPLETE trials at the single most common exact value, from the db copy):

| study | trials | at sentinel | value |
|---|---|---|---|
| curr_s3b | 60 | **41 (68.3%)** | −320 |
| curr_s3c | 80 | **56 (70.0%)** | −320 |
| curr_s3d | 40 | 24 (60%) | −320 |
| curr_s3e | 40 | 21 (52.5%) | −320 |
| curr_s3g | 40 | 14 (35%) | −320 |
| curr_s3f | 40 | 12 (30%) | −320 |
| curr_s3k | 80 | 21 (26.3%) | −320 |
| curr_s3h | 40 | 9 (22.5%) | −320 |
| curr_s3j | 40 | 6 (15%) | −320 |
| curr_s3i | 40 | 3 (7.5%) | −320 |
| v4b / v5 / v6 | 60/150/110 | 3.3% / 3.3% / 5.5% | −65 |

A quarter to seven-tenths of the stage-3 TPE budget sits on one flat value, which both wastes trials and poisons the TPE density model. Verified support: `TPESampler.__init__` accepts `constraints_func` in 5.0.0 (also NSGAIISampler, GPSampler; `CmaEsSampler` does not). Concretely:

```python
def feasibility(trial):           # >0 = violated; TPE samples feasibility separately
    return 0.0 if trial.user_attrs.get("real_gait") else 1.0
sampler = optuna.samplers.TPESampler(seed=…, multivariate=True,
                                     constraints_func=feasibility)
```
with the objective returning the raw kine_score (no −320 substitution) and raising/pruning on NaN. SCONE's equivalent is the `threshold` + `use_first_non_zero_result` machinery — same idea, tiered feasibility (§2 table).

### 3.2 Sampler options

- `multivariate=True` — our parameters are known to interact (v4b importance: `pf_gain` 30% *interacts with every weight multiplier* since it scales the same tables — `optuna_walk.py:76-91`). Univariate TPE cannot see ridges; CMA-ES (SCONE) would. This is the closest TPE gets.
- `group=False` for now — group-decomposed search needs the bigger spaces we're about to remove. Revisit only if a future study keeps >20 params with genuinely independent blocks.
- `constant_liar=True` — already the 5.0 default; **required** if we ever parallelize with `n_jobs>1` (don't, while MuJoCo/numpy threads contend — inference).
- `n_startup_trials=15–20` for a fresh study on the narrowed space (10 is thin when the first trials also have to map the new constraint).
- **No pruner.** Our objective is a single 16-s rollout — there are no intermediate reports to prune on; adding `MedianPruner` would do nothing. (Explicitly stating this because the ask asked.) The only "pruning" that helps is fail-fast on NaN, which the runner already does early (`t_end`-weighted NaN return, `optuna_walk.py:117`).
- Optional A/B: `GPSampler` (available in 5.0, constraints-capable) for a 40-trial pilot vs multivariate TPE. Keep TPE as default — it's what every recorded result in the db is calibrated against.

### 3.3 Shrink the space using the importance numbers I pulled

Trials-per-param collapsed as the curriculum accreted knobs (measured):

v5 10.0 (150/15) · v6 6.9 · s2b 3.85 · s3b 3.16 · s3c 3.81 · s3d 1.74 · s3e 1.67 · s3f 1.54 · s3g 1.43 · s3h 1.38 · s3i 1.33 · **s3j 1.29 (31 params/40 trials)** · s3k 2.50 (32/80).

fANOVA top-5 (full tables in `goal4_db_analysis.json`):

- **curr_s3k** (latest, 80 trials): `desc_e` .187, `ib_rge` .153, `ib_e_central` .131, `desc_f` .101, `contact_onset` .049 — then a tail of 27 params all ≤ .035, most ≤ .01.
- **curr_s3j**: `ia_f_central` .220, `desc_e` .142, `ii_e_central` .076, `ib_e_central` .060, `toe_rge` .060.
- **curr_s3i**: `ii_e_central` .122, `c1_gain` .108, `drive` .107, `pm_add` .086, `desc_f` .077.
- **ground_walk_v6**: `e2_pf` .241, `drive` .086, `rg_adapt` .082, `f1_kneext_inh` .069 (the study's own new knob only 4th), `rg_to_pf` .066.
- **ground_walk_v4b**: `pf_gain` .301, `desc_f` .142, `e2_pf` .130, `f1_df` .103.

Recommendation for the next stage-3-style study: **search ~8–10 params** — the consistent top tier (`desc_e`, `desc_f`, `ib_rge`, `ib_e_central`, `rg_to_pf`, `drive`, `contact_onset`, `ia_f_central`) **plus the one NEW knob under test** — and pin the rest at the chained winner's values. That restores ≥8–10 trials/param at an 80-trial budget. Two bound flags from the same data: s3k's best pins `pf_gain@lo` (0.3 — lower the bound to ~0.15 or re-center, log scale), and s2b's best (= its own seed, trial 0) pinned **all nine** new afferent gains at their lower bound 0 — those knobs were tolerated, not used; they are prime candidates to stop searching.

Also carry over unchanged: full-dict seeds (the JSON rule — `_curriculum.py:292-338`), and consider enqueueing **two** seeds (winner + runner-up) as a cheap multi-start nod to `CmaPoolOptimizer`.

---

## 4. Ensemble objective — design sketch (scaffold: `ensemble_objective_scaffold.py`)

Not wired in; pure-math parts smoke-tested (`python …\ensemble_objective_scaffold.py` → `ensemble smoke: -169.543`, `finalists smoke: [0, 1, 3]` — mean −167.67 minus 0.25×std 7.51).

**Score:** `score = mean_i(s_i) − λ·std_i(s_i)`, λ = 0.25 to start. On the s3k scale (best −159, spread of a few points between good trials) λ = 0.25 means a 10-point std across variants costs 2.5 points — robustness pressure without letting it dominate. **Inference (mine):** sweep λ ∈ {0.1, 0.25, 0.5} on a rerun of the existing s3k winner rather than inside the next study.

**Variants (N=3):**
1. **nominal** — exactly today's `--eval` (keeps the bit-exact regression anchor).
2. **plant jitter** — per-muscle-group Fmax × lognormal(σ=0.07), leg damping × lognormal(0.10), lateral rig k × lognormal(0.10); deterministic per trial number (same trial → same draws; different trials → different draws = domain randomization across the study, the SCONE NoiseController + sampled-initial-state pattern). **Needs one new runner-side Fmax hook** (gainprm[:,2] scale via env, default-off bit-identical) — doesn't exist yet; flagged in the scaffold.
3. **reference/pose variant** — left-leg target from its 2nd available cycle (measured: 3 left GRF onsets) via `kine_ref.compare(ref=…)`; pelvis_ty offset ±5 mm via the existing `AARL_PELVIS_TY` env (`_curriculum.py:104-110` already drives it).

Feasibility: ≥2 of 3 variants produce a real kine dict (and when 3.1 is adopted, this rides `constraints_func`, not the score).

**Budget:** measured eval cost = **median 1.0 min/trial (v10, 60 evals) and 1.1 min (v8b, 155 evals)** from the per-trial CSV timestamps → N=3 ≈ 3 min/trial ≈ 3 h for 60 trials. If that's too rich, use the scaffold's 2-stage mode: screen nominal-only, then ensemble-rescore trials within 5 points of the best (top 8) and pick the winner by ensemble score — basin_gate philosophy applied at study level.

**Per-variant scores go in `trial.set_user_attr("variant_scores", …)`** so post-hoc analysis (and the supervisor review) can see whether the ensemble changed the ranking or only re-ranked ties.

---

## 5. Notes and caveats

- **The v8b/v9/v10 studies are NOT in `optuna_walk.db`** — the db holds v1–v6 + curr_* only (19 studies; raw sqlite check shows study_ids run contiguously v6 → curr_s1, so v7–v10 were never in this file; that work ran on another machine per AGENTS.md). Their per-trial CSVs (`v8b/v9/v10_results.csv`) carry only 3–4 params + metrics, not full dicts, so param-importance can't be computed for them from here. The analysis above therefore grounds the recommendations on v4b/v5/v6 + the curr_s3 family, which is the current lineage anyway. If the laptop's db is reachable later, the analysis script (kept at `%TEMP%\goal4_importance.py`) reads it after one URL change.
- fANOVA importances at 1.3–2.5 trials/param are themselves poorly constrained — another reason to shrink the space before trusting any single ranking.
- SCONE citations are to the local install's files (2.4.4, verified on this machine by the scone skill) and the official docs pages fetched during this session; the `CmaPoolOptimizer` description quote comes from its doc page via search snippet — I could not open that page directly (wiki page-id not resolvable), so treat the wording as search-sourced.
- No existing study script, runner, or db row was modified. The db was analyzed on a %TEMP% copy.

## 6. Verified / not verified

**Verified (ran during this ask):**
- `cmd /c %TEMP%\goal4_run.cmd` (db copy listing) → 19 studies, counts/bests as tabled; `optuna 5.0.0`.
- `cmd /c %TEMP%\goal4_run2.cmd` (raw sqlite) → study_ids 1–19; no v7/v8b/v9/v10 rows.
- `cmd /c %TEMP%\goal4_run3.cmd` → full fANOVA + sentinel + bound-pinning analysis (archived: `goal4_db_analysis.json`).
- `cmd /c %TEMP%\goal4_run4.cmd` → optuna 5.0.0 sampler defaults + trials-per-param ratios.
- `cmd /c %TEMP%\goal4_run5.cmd` → `constraints_func` on TPESampler/NSGAIISampler/GPSampler; CmaEsSampler without; BoTorchSampler absent.
- `cmd /c %TEMP%\goal4_run6.cmd` → reference = 1 right cycle / 2 left cycles; `T_r 1.233 s, duty_r 0.61, knee_min_r −69.7°`.
- `cmd /c %TEMP%\goal4_run7.cmd` + `goal4_run8.cmd` → `subject01_walk1.mot` = 15 s but all-zero coordinates (placeholder).
- `cmd /c %TEMP%\goal4_run9.cmd` → median eval wall time 1.0 min (v10) / 1.1 min (v8b).
- `cmd /c %TEMP%\goal4_run10.cmd` → scaffold smoke: `-169.543`, `[0, 1, 3]`.
- Read (not run): loader/objective/kine_ref/runner code as cited by file:line; local SCONE tutorial/example .scone files; local scone skill; scone.software reference/FAQ/tutorials/publications pages; Geijtenbeek 2019 JOSS DOI page.

**Not verified / not done:**
- v8b/v9/v10 frozen-trial importance (studies absent from this db — see §5).
- No ensemble study was run; the scaffold's runner-hook (plant Fmax jitter) does not exist yet.
- `CmaOptimizerSpot`'s exact algorithm not characterized (page unreachable); only its existence in the reference manual is claimed.
- Song & Geyer 2015 (J Physiol 593(16):3495–3511) title/venue verified via search snippets; I did not open the full text in this session.
