# Goal 1 — SCONE tutorial campaign: full evaluation sweep + first optimizations on this machine

**Date:** 2026-09-25 (12:55–13:15 local) · **Machine:** EB475WS4 · **Status: PASSED** (all checks executed this session; no OpenSim fallbacks needed)
**Tooling:** `C:\Program Files\SCONE\bin\sconecmd.exe`, SCONE 2.4.4.3333, Hyfydy 1.12.6.1412 initialized (license active, per `sconecmd --version` output in every log below).
**Scenario source:** `C:\Users\Ben Bolen\Documents\SCONE\Tutorials3\` (untouched; optimizations ran on scratch copies under `%TEMP%`).

**What was never done before this session:** any optimization on this machine (`Documents\SCONE\results\` was empty at session start — verified) and any evaluation of Tutorials 1, 2a–2c, 3b, 4b–4d, 5a–5c, 6a–6d (3a + 4a were evaluated 2026-09-24, see `reports_20260924\scone_hyfydy_hands_on.md`).

---

## 1. Tutorial evaluations (all never-run scenarios; Hyfydy engine; `-l 2`; motion written per scenario)

Command pattern (one per row; full log = `logs\<tag>.log`, motion = `logs\motion_<tag>.sto`):

```
"C:\Program Files\SCONE\bin\sconecmd.exe" -e "C:\Users\Ben Bolen\Documents\SCONE\Tutorials3\<scenario> - Hyfydy.scone" -l 2 -r <...>\reports_20260925\logs\motion_<tag>
```

All 15 evaluations exited 0 on the Hyfydy engine — **no OpenSim fallback was needed anywhere** (the fallback branch in `tmp\run_evals.ps1` never fired). Scores are as printed by `-l 2` (lower = better for CMA-ES objectives; jump/straight-pose/6a/6d are maximization-style ScriptMeasure/JumpMeasure objectives — noted per row).

| Tutorial (Hyfydy variant) | result (default params) | key breakdown | sim time | x real-time | log / motion tag |
|---|---|---|---|---|---|
| 1 — Introduction (JumpMeasure, max_dur 2) | 42.3804 | jump_height 82.74 cm, early_jump_penalty 0 | 0.52 s | 92.7× | `tut1_intro` |
| 2a — High Jump (max_dur 2) | 42.3804 | jump_height 82.74 cm | 0.52 s | 83.8× | `tut2a_highjump` |
| 2b — High Jump Polynomial | 42.3804 | jump_height 82.74 cm (poly starts at constant 0.3 → identical controller to 1/2a at defaults; equality is a sanity check, not a copy-paste error) | 0.52 s | 63.4× | `tut2b_jump_poly` |
| 2c — Straight Pose Jump (init from `par/H0918v3.FC2.Jump.D2.par`, 27/45 imported) | 13.7741 (maximize) | JumpMeasure 101.27 cm, pelvis_tilt penalty −87.49 | 1.89 s | 110.6× | `tut2c_straight_pose` |
| 3b — Motor Noise Balance (max_dur 30) | 102.629 | BalanceMeasure 96.52, Effort 6.11, JointLimits 0 — **falls at 1.045 s** (same failure point as 3a baseline) | 1.045 s | 96.7× | `tut3b_motor_noise_bal` |
| 4b — Fast Gait (Gait15, max_dur 20) | 89.3969 | Gait 88.01: step_velocity 0.180 m/s (target 1.5), 11 steps; distance 4.65 m — falls | 3.73 s | 65.2× | `tut4b_fast_gait` |
| 4c — Perturbed Gait (±100 N pushes, max_dur 20) | 91.3812 | Gait 89.86: 0.113 m/s, 9 steps, 3.54 m — falls | 3.25 s | 53.1× | `tut4c_perturbed_gait` |
| 4d — Slippery Slope (`H0918v3_slope.hfd`, max_dur 20) | 97.7414 | Gait 96.53: 0.044 m/s, 3 steps, 1.60 m — falls fast on the slope | 1.195 s | 92.7× | `tut4d_slippery_slope` |
| 5a — Plantarflexor Weakness (Fmax×0.3 gastroc/soleus, H0914v3) | 79.4282 | Gait 78.68: 0.197 m/s, 7 steps, 2.56 m — best default-param walker of the pathology set | 2.415 s | 82.6× | `tut5a_pf_weakness` |
| 5b — Short Hamstrings (Lopt×0.8) | 96.933 | Gait 96.48: 0.020 m/s, 4 steps, 1.17 m — falls at 0.945 s | 0.945 s | 52.7× | `tut5b_short_hamstrings` |
| 5c — Hyper-reflexia (gastroc/soleus KF=1, C0=0.1 added) | 91.6098 | Gait 90.86: 0.097 m/s, 5 steps, 1.83 m — falls at 1.555 s | 1.555 s | 53.7× | `tut5c_hyper_reflexia` |
| 6a — Script Body Height (Lua ScriptController+ScriptMeasure) | 0.024403 (maximize) | runs headless; Lua measure on `calcn_r` height | 1.16 s | 67.3× | `tut6a_script_bodyheight` |
| 6b — Script Gyro Balance Gait (optional, tried) | 95.8342 | Gait 94.70: 7 steps, 1.94 m — falls at 1.755 s; Lua gyro messages fired (`gyro activated at t=1.035 torso_ori=-12.01 moment=14.30`) | 1.755 s | 77.1× | `tut6b_gyro_balance` |
| 6c — Script Reflex Modulation (optional, tried; loads `par/H0918GaitRS2Hfd4.par` 53/53) | 87.3623 | Gait 86.21: 10 steps, 4.45 m, sim 3.635 s — falls; Lua printed per-state reflex modulation | 3.635 s | 81.5× | `tut6c_reflex_mod` |
| 6d — Script Neural Delays (optional, tried) | 0.00946241 | per-muscle one-way delays 0.005–0.02 s applied via Lua; ScriptMeasure result | 0.615 s | 93.5× | `tut6d_neural_delays` |

No `CmaOptimizer.SimulationObjective.max_duration` override was needed: the longest scenario default is 30 s (3b), which costs ~0.35 s wall per single Hyfydy eval at ~90× real-time; shorter scenarios default to 2–20 s.

---

## 2. Optimizations — the first ever on this machine

Scratch copies were made with `robocopy` (`Tutorials3` → `%TEMP%\scone_opt_gait` and `%TEMP%\scone_opt_balance`); `Tutorials3` itself was never modified. The 4a scratch scenario got the warm-start block inserted after its `min_progress = 1e-4` line (the "init is an object child" gotcha — text edit, not an override):

```
	init { file = par/H0918GaitRS2Hfd4.par std_factor = 2 use_best_as_mean = 1 }
```

Optimization command template (full dot-path overrides; `lambda` = population; `random_seed` always set):

```
"C:\Program Files\SCONE\bin\sconecmd.exe" -o "<scratch>\Tutorial 4a - Gait - Hyfydy.scone" CmaOptimizer.SimulationObjective.max_duration=4 CmaOptimizer.max_generations=60 CmaOptimizer.lambda=16 CmaOptimizer.random_seed=42
```

Results land in `C:\Users\Ben Bolen\Documents\SCONE\results\<timestamp>.<signature>.R<seed>\` — each dir holds its own `config.scone` (includes inlined; duration baked in), model, init state, `history.txt` (columns: generation / best_fitness / median_fitness / predicted_fitness / fitness_progress), `optimization.log`, and `NNNN_<gen-fitness>_<best-so-far>.par` per new global best.

| run (log) | dims | gens / wall | best fitness | best par | results dir |
|---|---|---|---|---|---|
| **4a Gait, Hyfydy, warm-started, D=4** (`logs\opt_gait_4a.log`) | 53 | 60 / 4.35 s | **0.858** (gen 44; gen-0 warm start 0.901) | `...\260925.130052.H0918v3.RS2.S10WA3K1G14.D4.R42\0044_0.873_0.858.par` | `260925.130052.H0918v3.RS2.S10WA3K1G14.D4.R42` |
| **3a Balance, Hyfydy, from scenario seeds, D=12** (ask's ~40-gen budget) (`logs\opt_balance_3a.log`) | 36 | 40 / 0.97 s | **80.445** (gen 37; gen-0 96.56) | `...\260925.130211.H0918v3.R36.BW.D12.R42\0037_94.802_80.445.par` | `260925.130211.H0918v3.R36.BW.D12.R42` |
| **3a Balance EXTENDED, same seed/settings, D=12** (bonus: the 40-gen run found its biggest jump at gen 37/40 — clearly unconverged) (`logs\opt_balance_3a_ext300.log`) | 36 | 300 / 18.76 s | **2.207** (gen 293; 56 new-best generations, monotone) | `...\260925.130428.H0918v3.R36.BW.D12.R42\0293_37.931_2.207.par` | `260925.130428.H0918v3.R36.BW.D12.R42` |
| **2a High Jump, D=2 (scenario default)** (ask's "small" run) (`logs\opt_highjump_2a.log`) | 27 | 20 / 0.35 s | **93.606** (gen 19; gen-0 49.67) | `...\260925.130348.H0918v3.FC2.Jump.D2.R42\0019_82.742_93.606.par` | `260925.130348.H0918v3.FC2.Jump.D2.R42` |
| **4a Gait BONUS at D=20** (full scenario horizon; run to test the overfit finding below) (`logs\opt_gait_4a_d20.log`) | 53 | 60 / 17.00 s | **0.873** (gen 7; gen-0 0.883) | `...\260925.130924.H0918v3.RS2.S10WA3K1G14.D20.R42\0007_0.921_0.873.par` | `260925.130924.H0918v3.RS2.S10WA3K1G14.D20.R42` |

Convergence shapes (from the `B=` new-best markers in each log + `history.txt`):

- **4a Gait D4:** new bests at gens 0, 1, 2, 6, 12, 14, 44 (0.90→0.87→0.86…); `history.txt` min best_fitness = 0.858458; flat plateau gens 15–59 (population medians spiky: 0.88 to 31.9 — occasional falling samples).
- **3a Balance 40:** steady improvement 96.56 → 80.45 over 25 new-best generations, then the cap hit **one generation after the biggest jump (gen 37)** — unconverged, motivating the extension.
- **3a Balance 300:** monotone 96.56 → 80.44 (g37) → 17.72 (g172) → 2.95 (g202) → 2.207 (g293); population median collapsed from ~95 to ~2.4 (gens 283–299 medians 2.4–2.7 = whole population now stands).
- **2a High Jump:** 49.67 → 55.45 (g5) → 82.74 (g6) → 93.61 (g19); **still climbing steeply at the 20-gen cap** — treat 93.6 cm as a lower bound, not a converged optimum.
- **4a Gait D20:** new bests only at gens 0, 2, 7 (0.883 → 0.873); flat afterward — the warm start is already near-optimal on the full-horizon objective.

### Re-evaluations vs baselines (all exit 0; `-e <best.par> -l 2` inside each results dir; full logs `logs\reeval_*.log`)

Every optimizer best par reproduced its recorded fitness **bit-exactly** in re-evaluation (0.858458, 80.445, 2.20727, 93.6057, 0.873108) — same-machine determinism confirmed.

**Gait 4a (baseline: pretrained `par/H0918GaitRS2Hfd4.par` = 7 steps / 2.58 m at 3 s window, 09-24 session):**

| par | horizon | result | steps | distance | step_velocity | sim time (survived?) | log |
|---|---|---|---|---|---|---|---|
| pretrained (gen-0 warm start `0000_23.911_0.901.par`) | 4 s | 0.90058 | 8 | 4.38 m | 1.103 m/s | 4.0 s ✔ | `reeval_gait4a_init_d4` |
| **D4-optimized `0044_0.873_0.858.par`** | 4 s | 0.858458 | 8 | 4.29 m | 1.092 m/s | 4.0 s ✔ | `reeval_gait4a_best_d4` |
| pretrained | 20 s | 0.882806 | **36** | **21.48 m** | 1.076 m/s | **20 s ✔ walks the whole window** | `reeval_gait4a_init_d20` |
| D4-optimized | 20 s | 83.29 (fell) | 11 | 5.05 m | 0.194 m/s | **4.685 s ✖ FALLS** | `reeval_gait4a_best_d20` |
| **D20-optimized `0007_0.921_0.873.par`** | 20 s | 0.873108 | 36 | 21.08 m | 1.055 m/s | **20 s ✔** (GRF penalty 0.0057 vs pretrained 0.0314) | `reeval_gait4a_bestD20_d20` |

**Headline: the 4-s-window optimization (the ask's D4–6 guidance) IMPROVED the short-window fitness (0.901→0.858) but produced a gait that falls at 4.7 s — it does not survive the scenario's native 20-s horizon, while the pretrained par walks 20 s / 36 steps / 21.5 m.** This is the skill's "falling is efficient" trap in overfit form: the short window can't see the fall. Optimizing at the full 20-s horizon (bonus run) keeps robustness (20 s / 36 steps / 21.1 m) and still improves the true objective (0.883→0.873) with a 5.5× lower GRF penalty. **Recommendation: for any future gait campaign on this machine, set `max_duration` = full horizon (or ≥ 2× expected fall time), not a short window.**

**Balance 3a (baseline: default params fall at ~1.0 s; 09-24: 65.2/100 at 1.045 s):**

| par | horizon | result | BalanceMeasure | sim time | verdict | log |
|---|---|---|---|---|---|---|
| default/init (`0000_99.012_96.562.par`) | 30 s | 101.76 | 96.53 | 1.04 s | falls ~1 s — reproduces the 09-24 baseline ✔ | `reeval_bal_init_d30` |
| 40-gen best (`0037_94.802_80.445.par`) | 12 s (its own config) | 80.445 | 76.58 | 2.81 s | falls at 2.81 s — better than baseline, still falls | `reeval_bal40_best_d12` |
| 40-gen best | 30 s | 94.495 | 90.63 | 2.81 s | same fall, longer window | `reeval_bal40_best_d30` |
| **300-gen best (`0293_37.931_2.207.par`)** | 12 s | **2.20727** | **0** | **12 s ✔** | stands the full window, zero balance penalty, effort 220.7 | `reeval_bal300_best_d12` |
| **300-gen best** | 30 s (scenario default) | **2.19135** | **0** | **30 s ✔** | **stands the FULL 30-s default horizon without falling** | `reeval_bal300_best_d30` |

**Headline: CMA-ES found the vestibular gains.** The scenario's `BodyPointReflex` torso (vestibular-like) gains start at 0 (`$KP = 0~0.1`, `$KV = 0~0.1` in `controllers/ControllerReflexBalance.scone`); the 300-gen run drove them to clearly nonzero values (e.g. `Balance.iliopsoas-torso.KV = 0.727`, `Balance.bifemsh-torso.KP = -1.044`, `Balance.vasti-torso.KP = -0.369`, `Balance.hamstrings-torso.KP = -0.296`) and the model stands the full 30 s. The ask's 40-gen budget alone reaches only "falls at 2.8 s" — the improvement to "stands 30 s" happens in generations 54–293. The 40-gen and 300-gen runs are seed-identical, so the 300-gen run strictly contains the 40-gen trajectory.

**High Jump 2a (baseline: default params 42.38 / 82.74 cm):** best par re-eval = 93.6057, jump_height 93.90 cm (**+13.5 % over default**, `logs\reeval_jump2a_best.log`, sim 1.36 s, early_jump_penalty 0.299). Not converged at 20 gens (still improving at the cap) — a longer run is trivially cheap (~0.02 s/gen) if Ben wants the ceiling.

---

## 3. Where everything lives

- **This report:** `D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\goal1_scone_tutorials_and_optimization.md`
- **Campaign logs + motions (25 `.sto`, 56 files total):** `...\reports_20260925\logs\` — eval logs `tut*.log`, optimization logs `opt_*.log`, re-eval logs `reeval_*.log`, motions `motion_*.sto`. (`reports_20260924\scone_logs\` untouched.)
- **Optimization results (5 dirs, each self-contained with `config.scone`):** `C:\Users\Ben Bolen\Documents\SCONE\results\260925.{130052,130211,130348,130428,130924}.*.R42\` — listed in §2 with best-par filenames.
- **Batch scripts (reproducible):** `...\reports_20260925\tmp\run_evals.ps1` (evaluations + fallback logic), `tmp\setup_scratch.ps1` (robocopy + warm-start edit), `tmp\run_reevals.ps1` (in-place + full-horizon re-evals + 6b/6c/6d).
- **Scratch copies:** `%TEMP%\scone_opt_gait\`, `%TEMP%\scone_opt_balance\` (disposable; results dirs are the archive).

## 4. Exact commands actually executed

Evaluations (representative; all 15 in `tmp\run_evals.ps1`):
```
"C:\Program Files\SCONE\bin\sconecmd.exe" -e "C:\Users\Ben Bolen\Documents\SCONE\Tutorials3\Tutorial 4b - Fast Gait - Hyfydy.scone" -l 2 -r D:\...\reports_20260925\logs\motion_tut4b_fast_gait
```
Optimizations:
```
"C:\Program Files\SCONE\bin\sconecmd.exe" -o "%TEMP%\scone_opt_gait\Tutorial 4a - Gait - Hyfydy.scone" CmaOptimizer.SimulationObjective.max_duration=4 CmaOptimizer.max_generations=60 CmaOptimizer.lambda=16 CmaOptimizer.random_seed=42
"C:\Program Files\SCONE\bin\sconecmd.exe" -o "%TEMP%\scone_opt_balance\Tutorial 3a - Balance - Hyfydy.scone" CmaOptimizer.SimulationObjective.max_duration=12 CmaOptimizer.max_generations=40 CmaOptimizer.lambda=16 CmaOptimizer.random_seed=42
"C:\Program Files\SCONE\bin\sconecmd.exe" -o "%TEMP%\scone_opt_balance\Tutorial 3a - Balance - Hyfydy.scone" CmaOptimizer.SimulationObjective.max_duration=12 CmaOptimizer.max_generations=300 CmaOptimizer.lambda=16 CmaOptimizer.random_seed=42
"C:\Program Files\SCONE\bin\sconecmd.exe" -o "%TEMP%\scone_opt_gait\Tutorial 2a - High Jump - Hyfydy.scone" CmaOptimizer.max_generations=20 CmaOptimizer.lambda=16 CmaOptimizer.random_seed=42
"C:\Program Files\SCONE\bin\sconecmd.exe" -o "%TEMP%\scone_opt_gait\Tutorial 4a - Gait - Hyfydy.scone" CmaOptimizer.SimulationObjective.max_duration=20 CmaOptimizer.max_generations=60 CmaOptimizer.lambda=16 CmaOptimizer.random_seed=42
```
Re-evaluations (representative; run with cwd = results dir or scratch dir holding `config.scone` + renamed par):
```
cd C:\Users\Ben Bolen\Documents\SCONE\results\260925.130052.H0918v3.RS2.S10WA3K1G14.D4.R42
"C:\Program Files\SCONE\bin\sconecmd.exe" -e 0044_0.873_0.858.par -l 2 -r <...>\logs\motion_reeval_gait4a_best_d4
```
Full-horizon re-evals used a scratch `config.scone` (copy of the scenario, native 20/30-s duration) beside a renamed copy of the best par — the documented `-e` recipe (config next to par); no scenario-value override is honored reliably in `-e` mode, so none was attempted there.

## 5. Honesty notes / what was NOT done

- Tutorials 6b/6c/6d were OPTIONAL (ask: "try one"); all three ran headless, so all three are reported. Nothing was skipped for GUI dependence.
- No OpenSim-variant evaluations were run (nothing to fall back from). The ask's fallback path was armed but never fired.
- Wall-clock: whole campaign ≈ 20 min of the ~100-min budget; the two ask-mandated optimizations together took 5.3 s of solver time. Budget headroom went to the two bonus runs (balance ×300 gens, gait D=20) — both labeled above.
- Seed: all optimizations `random_seed=42` as asked. Same-seed reruns are deterministic on this machine (verified by bit-exact re-evals); cross-machine reproduction is NOT expected (floating-point accumulation, per skill FAQ).
- The 40-gen balance run's fitness at its cap (80.45) is real but NOT the capability ceiling — see the 300-gen result before quoting it anywhere.
- `history.txt` column semantics read from its header row (`generation best_fitness median_fitness predicted_fitness fitness_progress`); `best_fitness` is per-generation best, the global best is the minimum across generations (= the last `NNNN_*_<best>.par` filename).
