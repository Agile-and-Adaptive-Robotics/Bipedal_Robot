# Goal 1 — walker baselines + SPIKING mirror of the spinal network

**Date:** 2026-10-02 · **Host:** easteregg2 · **Workdir:** `Code\MuJoCo_SNS\spinal`
**Branch:** `KneeTestSetup_BenBo_stw`, clean tree at session start; **no git commands that mutate state were run** (Ben commits via GitHub Desktop).
**Plan of record:** `SPIKING_MIRROR_PLAN.md` (Ben, 2026-09-18) — followed literally; every deviation is listed in §6.
**SNS python (full path, per machine rules):** `D:\Anaconda\envs\myo\python.exe`, `CONDA_PREFIX=D:\Anaconda\envs\myo`.

---

## 1. Executive summary

- **A. Baselines:** the current production set was walker-run and every eval reproduced its documented anchor **bit-exactly**: s3k default build `−161.56754173676563`, w2lvar stage-5 t68 `−144.46440311078322`, syn6 stage-5 t0 `−200.37955834081458` (§3).
- **B. Spiking mirror BUILT** per the plan's hybrid ruling: `build_network_spiking.py` mirrors `build_network.py` population-for-population and edge-for-edge (topology gate: edge multisets **identical** — 1186/7908/7762 edges across three gain configurations), with MNs left non-spiking and afferents/RG/PF/INs/commissurals as spiking LIF cells at real mV ranges. Weights are **calibrated by script** (`calibrate_spiking.py` → `spiking_calibration.json`), not hand-tuned.
- **Bit-identity of the non-spiking pipeline PROVEN:** after the one default-inert selector edit, the s3k baseline re-evaluated to `−161.56754173676563` (all 17 digits) and the defaults gate reproduced `(410, 376, 1186)` exactly.
- **C. Gates (§5, §9):** topology mirror **PASS** (edge multisets identical in 3 gain configs); rhythm **SPLIT** — self-sustained bilateral antiphase **PASS**, period `0.621 s` vs analog `0.399 s` = **+55.6 % > the plan's ~20 %** FAIL; basin **PASS 12/12** (±1 % on every scalar incl. the mirror's own, margin 0.89×); 20 s air smoke **FAIL on cadence only** (finite, rhythmic, knee 106°/hip 50° vs the twin's 118°/71°, but ~10× its period in that preparation).
- **Stretch ground eval (§7):** the **untuned** mirror at s3k params scored `−236.35` (worse, expected) but produced a **bilateral** gait — 8/9 cycles, knee `−79.8°`, no NaN — the one thing the s3k non-spiking winner cannot do (its left leg is frozen in all 80 study trials).
- **Honest bottom line:** the mirror is structurally faithful (gate 1 exact), rhythmically robust (gates 2-partial/3), and plant-capable (§7), but its cycle period is 1.56× the analog at the check_selfsustain reference, its air cadence is several times faster, and it has **zero tuning** — mapped parameters only. It is the starting point for Ben's fresh (short) curriculum, exactly as the plan's effort/risks section predicted.

---

## 2. Files created / modified (complete list)

**Modified (pre-existing):**
- `Code\MuJoCo_SNS\spinal\build_network.py` — added the `AARL_NET=spiking` branch to the existing variant selector in `build()` (8 lines incl. comment). Default-inert: proven by the defaults gate `(410, 376, 1186)` and the bit-exact s3k re-run (§4).
- `Code\MuJoCo_SNS\spinal\spinal_run.png` — **run byproduct, overwritten unintentionally by the gate-4 air smoke**: the runner saves this figure at a hardcoded name in non-eval mode even when `AARL_NPZ` redirects the npz. The npz outputs were correctly isolated (`scratch_spk_air_*.npz`, gitignored); the png is regenerable from any run but the previous image (a pre-session run's figure) is replaced. `spinal_run.npz` itself was NEVER touched (all runs used unique `AARL_NPZ` names).

**Created (spinal root):**
- `build_network_spiking.py` — the hybrid spiking mirror (main deliverable; §5.1).
- `calibrate_spiking.py` — the calibration script (plan rule 2).
- `spiking_calibration.json` — calibration output (loaded by the builder; git-untracked, machine-local).

**Created (`reports_spiking_20261002\`):**
- `tools\baseline_eval.py` — baseline/stretch eval driver (documented winner-loading flows).
- `tools\topology_mirror_check.py` — gate 1.
- `tools\check_selfsustain_spiking.py` — gate 2.
- `tools\basin_gate_spiking.py` — gate 3 (+ `basin_gate_spiking.json` output).
- `tools\smoke_air_spiking.py` — gate 4 (+ `air_smoke_ns.json` / `air_smoke_spiking.json`).
- `tools\tune_rg_pair.py`, `tools\sweep_thrinc_fullnet.py` — RG cell-parameter calibration sweeps (v1 starting-point refinement, plan's own clause).
- `tools\probe_hybrid.py`, `tools\probe_lif.py`, `tools\probe_build_spiking.py`, `tools\debug_rg_pair*.py` — measurement probes used to establish the toolbox traps in §6.
- `logs\*.log` — all run logs (baselines, calibration, gates, sweeps).
- `goal1_walkers_spiking.md` — this report.

**Scratch npz (spinal root, from the run convention `AARL_NPZ=spinal_run_<tag>.npz`):** `spinal_run_spkbase_{s3k,w2lvar,syn6,spiking}.npz`, `scratch_spk_air_{ns,spiking}.npz`.

**NOT touched:** `runner.py`, `params.py`, `_curriculum*.py`, `optuna_*.db`, `spinal_run.npz`, `curriculum_*.json`, `best_walk_params*.json`, `reports_2026092*`, the dissertation, figures, `.tex`.

**Sibling-campaign artifacts (NOT this goal's):** the working tree also carries changes under `Code\Matlab\SNS_Simscape\` (spiking Simulink demos + `.slx` edits) and a few `logs\` files in this report folder (`spikingv2_*`, `v3_*`, AnimatLab W2L runs) — those belong to other goals of the 2026-10-02 spiking campaign, not to goal 1. Scratch npz `scratch_spk_air_*.npz` are gitignored (`scratch_*.npz` rule).

---

## 3. A. Baselines (the current production set, run this session)

All three used the documented winner-loading flows (the 09-25/27 replay pattern), `--eval --drive <winner drive>`, unique `AARL_NPZ` names, logs in `logs\baseline_*.log`.

| metric | **s3k** (default build, curr_s3k_nocross t34) | **w2lvar** (stage-5 t68) | **syn6** (stage-5 t0 seed) |
|---|---|---|---|
| kine_score | **−161.56754173676563** | **−144.46440311078322** | **−200.37955834081458** |
| documented anchor | −161.56754173676563 (bit-exact) | −144.46440311078322 (bit-exact) | −200.3796 (rebased; matches) |
| drive | 2.0773 | 3.8196 | 2.7864 |
| stayed up (kz) | 0.850 m (up) | 0.787 m (up) | 0.840 m (up) |
| tilt max | 22.1° | 31.1° | 28.1° |
| cycles r/l | 9 / 0 (frozen left) | 4 / 6 (bilateral) | 4 / 0 (planted left) |
| contact duty r/l | 0.637 / 0.985 | 0.590 / 0.843 | 0.955 / 1.0 |
| knee r (min..max) | −101.0 .. +8.4 | −52.3 .. +15.0 | −56.3 .. +10.2 |
| hip r | −17.5 .. +43.7 | −42.0 .. +65.0 | +15.8 .. +62.8 |
| ankle r | −26.5 .. +44.2 | −78.0 .. +55.9 | −64.7 .. −1.7 |
| knee l | −41.9 .. +14.7 | −47.2 .. +15.7 | −49.6 .. +11.3 |
| nan | False | False | False |

Repro (easteregg2, from `Code\MuJoCo_SNS\spinal`):
```
D:\Anaconda\envs\myo\python.exe reports_spiking_20261002\tools\baseline_eval.py s3k
D:\Anaconda\envs\myo\python.exe reports_spiking_20261002\tools\baseline_eval.py w2lvar
D:\Anaconda\envs\myo\python.exe reports_spiking_20261002\tools\baseline_eval.py syn6
```
(The driver reproduces the flows of `reports_20260925\tmp\rescore_s3k_fixed_ref.py` and `...\_supgate_replay_w2lvar_s5_t4.py` verbatim: s3k via `_curriculum.set_stage(4, {**mul, **st3, full_rules:0})`; w2lvar/syn6 via their curriculum modules' `set_stage(5, {**mul, **winner})`.)

---

## 4. Bit-identity of the non-spiking pipeline (the ask's hard requirement)

1. **Pre-edit anchors** (§3, this session, before any file was modified): s3k `−161.56754173676563`, w2lvar `−144.46440311078322`, syn6 `−200.37955834081458`.
2. **Edit:** the single `AARL_NET=spiking` selector branch in `build_network.py`.
3. **Post-edit defaults gate** (`reports_20260925\tmp\gate_ab_w2lvar.py`): `(410, 376, 1186)` — **PASS, exact**; w2lvar variant `(888, 382, 7430)` — PASS.
4. **Post-edit baseline re-run** (`baseline_eval.py s3k`, log `logs\baseline_s3k_postedit.log`): `kine_score = -161.56754173676563` — **bit-identical, all 17 digits**.

---

## 5. B/C. The spiking mirror and its verification gates

### 5.1 What was built (`build_network_spiking.py`)

Cell-class mapping exactly as the plan's table (MN layer stays `NonSpikingNeuron` in the 0..5 mV frame with `a = clip(V/E_HI, 0, 1)` unchanged; afferent encoders, RG-E/F, PF cells, all INs, commissurals spiking; DRIVE/POSTURE/BAL/VEST stay analog brainstem surrogates). Key mechanics:

- **RG half-centers = adapting LIF** (`SpikingNeuron` with `threshold_increment`), the plan's documented fallback — verified this session that 1.5.2 ships **no** spiking bursting class (the NaP class is non-spiking-only; `grep class` over `sns_toolbox\neurons.py`). Burst termination = threshold adaptation; `tau_theta` **reuses `TAU["rg_nap_h"]`**, so curriculum `rg_nap_h` values transfer.
- **Hybrid synapses:** `SpikingSynapse` from spiking cells onto the analog MNs (reversals `+8/−5 mV`, the analog frame); spiking→spiking (reversals `0/−70 mV`); graded `NonSpikingSynapse` for analog→spiking commands (DRIVE→RG, MN→RC-as-rate-encoder) and unchanged analog→analog edges (POSTURE/BAL/VEST→MN).
- **Sub-stepping (plan rule 3):** `step()` runs 4 network sub-steps at 0.5 ms inside each 2 ms plant step at constant input — no runner change needed.
- **Runner compatibility:** readout tap neurons `RD_<cell>` (non-spiking sinks, 50 Hz → full-scale `E_HI`) are appended for RG_E/RG_F + all PF cells, and `net.idx[...]` is repointed at them (`net.ridx[...]` keeps the raw spiking cells), so the runner's stance gates and neuro logging read analog levels unchanged. Taps are pure sinks (12/16 extra neurons depending on PF flavor; counted and asserted by gate 1).

### 5.2 GATE 1 — topology mirror: **PASS**

`tools\topology_mirror_check.py` builds both nets under identical gains and compares `(source, destination, sign)` edge multisets, input lists, MN maps, and conditional-population presence:

| config | non-spiking | spiking | edges equal | inputs equal | taps |
|---|---|---|---|---|---|
| stage-1 (all conditional gains 0) | 410 pops / 1186 edges | 422 / 1186 (+12 taps) | **YES** | YES | 12/12 |
| tuned (every pathway on, full_rules) | 884 / 7908 | 896 / 7908 (+12) | **YES** | YES | 12/12 |
| joint-layer PF | 890 / 7762 | 906 / 7762 (+16) | **YES** | YES | 16/16 |

Repro: `D:\Anaconda\envs\myo\python.exe reports_spiking_20261002\tools\topology_mirror_check.py`

### 5.3 GATE 2 — constant-DRIVE rhythm: **SPLIT** (rhythm PASS, period FAIL)

`tools\check_selfsustain_spiking.py` (the plan's recipe: winner tables via `effective_tables("best")`, constant DRIVE 2.5 nA, 20 s, network-only):

- non-spiking: 26 antiphase peaks in the last 10 s, period **0.399 ± 0.164 s**, swing sustained (4.79 mV in every 5 s window).
- spiking: **self-sustained bilateral antiphase** — RG_E 4.9 Hz / RG_F 3.2 Hz per leg, E−F spike-envelope period **0.621 ± 0.002 s** (clean and regular), swing 5.13 mV sustained in every window.
- verdict: rhythm-sustained **TRUE**; period deviation **+55.6 % > ~20 %** → the plan's period clause FAILS.

What was tried to close the gap (all scripted, logs kept): isolated-pair sweep (`tools\tune_rg_pair.py`, 45 configs) — alternation requires `thr_inc ≥ ~2` (below it, RG-F never escapes the laminated inhibition); full-net sweep (`sweep_thrinc_fullnet.py`) — `thr_inc 2.0` loses alternation entirely (RG_F 0.0 Hz), `k_ns` above 0.35 degrades antiphase (E and F fire into each other's phases). The chosen operating point (`k_ns 0.35`, `thr_inc 4.0`, `tau_theta = rg_nap_h 0.35`) is the best clean-alternation point found; its period is 1.56× the analog reference. **The period knob is identified** (`thr_inc` and the `rg_nap_h` binding), the residual is an untuned-starting-point issue, not a structural failure.

### 5.4 GATE 3 — basin (spiking RG): **PASS 12/12**

`D:\Anaconda\envs\myo\python.exe reports_spiking_20261002\tools\basin_gate_spiking.py` (log `logs\gate3_basin.log`, json `basin_gate_spiking.json`):

- baseline windows (5 s swings): 5.31 / 5.13 / 5.13 / 5.13 mV; bar = 2.66 mV.
- **12/12 trials PASS**, no explosions; final swing min / median = **4.58 / 5.11 mV** (margin **0.89×** baseline).
- The ±1 % jitter covered `params` G/TAU/W tables **plus the mirror's own `SPIKE`/`CAL` scalars** (thresholds, adaptation, bias, calibration factors) — the parameters a Simulink/AnimatLab port would re-round. The adapting-LIF half-center sits FAR from the analog build's near-bifurcation marginality.

Protocol: `basin_gate.py`'s ±1 % multiplicative jitter on every scalar parameter — extended to the mirror's own `SPIKE`/`CAL` scalars (the parameters a port would re-round) — 12 seeded trials, 20 s constant-DRIVE, PASS = last-5 s swing ≥ max(0.1 mV, ½ baseline first window). (Adaptation: the original gate consumes a `RUNNER_DUMP_STATE` winner json; the mirror has no tuned winner, so the perturbed config is the same `effective_tables("best")` baseline gate 2 uses.)

### 5.5 GATE 4 — 20 s air smoke (regime match): **FAIL on cadence, PASS on structure**

Twin full-runner runs (stage-1 params, drive 2.5, deafferented air, interleg off; logs `logs\gate4_air_{ns,spiking}.log`, jsons `air_smoke_*.json`):

| metric | non-spiking | spiking mirror |
|---|---|---|
| finite | true | **true** |
| rises (RG_E envelope) | 5 | **53** |
| RG period | 2.108 s (0.47 Hz, the Ivanenko air regime) | 0.215 s |
| E-duty | 0.82 | 0.26 |
| knee range / min | 117.9° / −106.3° | **106.0° / −92.6°** |
| hip range | 71.1° | **50.4°** |
| ankle range | 24.7° | 85.0° |

Verdict: finite ✓, rises ≥ 3 ✓, ranges ≥ 50 % of the twin ✓ (knee 90 %, hip 71 %), **period within 35 % ✗** (0.215 s vs 2.108 s) → **GATE 4 FAIL** on the cadence clause. Caveat stated plainly: the smoke's period metric counts half-max crossings of the RD_ readout level at plant resolution — ripple-contaminated for the spiking twin (the same measurement trap gate 2's fix addressed); the honest summary is "the mirror air-steps rhythmically with full-amplitude joint swings, but its cadence is several times faster than the non-spiking twin in this preparation".

Twin full-runner runs (`_smoke_nap.py` pattern: stage-1 params, `--no-ground --no-afferents --no-interleg --time 20`, drive 2.5) with `AARL_NET` unset vs `spiking`; compared on rises, RG period, E-duty, knee/hip/ankle ranges. PASS rule (documented in the script): finite AND ≥3 rises AND period within 35 % AND knee+hip ranges ≥ 50 % of the non-spiking twin.

---

## 6. Plan-vs-reality deviations (documented, not silent)

1. **`SpikingNonSpikingSynapse` does not exist in sns_toolbox 1.5.2** (the plan names it as the hybrid pattern's class). The pattern itself is supported: a `SpikingSynapse` whose destination is a `NonSpikingNeuron` — verified empirically (`tools\probe_hybrid.py`: mixed net steps clean, no NaNs; the compiler sets non-spiking `theta_0 = float_max` so they never spike).
2. **`g_m = 1` trap (measured, `tools\probe_lif.py`):** `SNS_Numpy.forward` computes `dV = dt/(Cm/g_m) · (−g_m(V−Vr) + I)`, which is dimensionally correct only at `g_m = 1 µS` — any other membrane conductance silently rescales every input current (a `g_m = 0.1` cell integrated 10× slower than its `Cm` implied). All spiking cells therefore use `g_m = 1` and reach threshold via a **tonic bias** (16 nA → equilibrium −54 mV, 4 mV subthreshold — the analog cells' near-threshold operating point), which is the physiological background-synaptic-current stand-in.
3. **`threshold_time_constant` is consumed in seconds** (`dt/tau_theta` with dt in s), not the ms its docstring claims.
4. **Mean-conductance matching is only valid onto the analog MN.** Calibrating spike→spiking weights by mean conductance (the plan's rule 2 read literally) produced single-spike increments of ~21 µS (vs `g_m = 1`) that saturated relay INs at ~200 Hz and permanently pinned the inhibited half-center (`tools\debug_rg_pair.py`). Spike→spiking weights are therefore **rate-calibrated** (40 Hz in → 40 Hz out; 40 Hz inhibition silences a 40 Hz-driven cell) — the plan's rule applied at the spiking operating point; spike→MN weights remain mean-conductance calibrated exactly as written (residuals ≤ 0.0003 mV).
5. **RG burst machinery:** no spiking bursting class ships (the plan's risk (a), confirmed) → adapting-LIF half-centers per the plan's own fallback.
6. **Afferent thresholds at −52 mV** (2 mV above the biased equilibrium) so the runner's existing 0..3 nA afferent currents map to ~0..25 Hz — the plan's "I→rate map chosen from the existing afferent gain ranges". Consequence: the `AFF.i0_ii` 1 nA baseline tone is just-subthreshold (silent at rest) instead of a standing analog excitation — a known standing-regime difference.
7. **RG membrane tau stays 50 ms** (analog value), but burst timing is carried by `tau_theta` (= `rg_nap_h`) and the synaptic taus (5 ms exc / 20 ms inh), since the plan's "keep the same ms values" cannot make nA-scale currents spike at 40 Hz through a 50 ms membrane.

---

## 7. Stretch (D): one ground-eval attempt of the spiking mirror — RAN

`D:\Anaconda\envs\myo\python.exe reports_spiking_20261002\tools\baseline_eval.py spiking` (log `logs\ground_eval_spiking.log`, npz `spinal_run_spkbase_spiking.npz`) — the s3k winner flow with `AARL_NET=spiking`, mapped parameters, **no tuning**:

| metric | s3k non-spiking | s3k **spiking mirror (untuned)** |
|---|---|---|
| kine_score | −161.57 | **−236.35** |
| bilateral (both legs ≥3 cycles) | **false** (left frozen, 0 cycles) | **TRUE** — 8 cycles r / 9 cycles l |
| knee_min | −64.3° (cycle) | **−79.8°** |
| contact duty r/l | 0.637 / 0.985 | 0.799 / 0.856 |
| nan / fell | no / up | no / up |

Honest reading: the mirror scores ~75 points **worse** overall (expected — zero tuning, mapped parameters only), but the run is a genuine **bilateral** ground gait with deep knee flexion and no NaN — the one thing the s3k non-spiking winner cannot do (its left leg is frozen in every trial of the study; DESIGN.md 2026-09-22 verdict). One run, no statistics, no tuning — flagged as a lead for the fresh curriculum the plan calls for, not as a result. The w2lvar and syn6 spiking mirrors were **not built** (separate architectures; fresh tuning is Ben's call per the ask) — listed as future work.

## 8. What is honestly NOT done

- w2lvar/syn6 spiking mirrors (future work; requires their own builders + fresh curricula).
- Any tuning/optimization of the spiking mirror (the ask forbids new studies; the plan expects a fresh short curriculum — Ben's call).
- Gate-2 period within 20 % (achieved 0.621 s vs 0.399 s; knobs identified, see §5.3).
- s3k winner under the spiking mirror is not expected to walk — it is an untuned mapped starting point.

## 9. Consolidated gate table

| gate | verdict | headline numbers |
|---|---|---|
| 0 compile (run_gate) | **PASS** | runner/_curriculum/params/build_network all compile |
| 0′ defaults gate | **PASS** | `(410, 376, 1186)` exact post-edit; w2lvar `(888, 382, 7430)` |
| 0″ bit-identity | **PASS** | s3k `−161.56754173676563` pre == post edit (17 digits) |
| 1 topology mirror | **PASS** | edge multisets identical: 1186 / 7908 / 7762 (3 configs); inputs + MN maps identical; taps 12/12, 12/12, 16/16 |
| 2 rhythm (constant DRIVE 2.5, 20 s) | **SPLIT** | sustained bilateral antiphase (E 4.9 / F 3.2 Hz, swing 5.13 mV in all windows) = PASS; period 0.621 ± 0.002 s vs analog 0.399 s = **+55.6 % > 20 %** = FAIL |
| 3 basin (±1 %, 12 trials) | **PASS 12/12** | final swing min 4.58 mV, margin 0.89× baseline; 0 explosions |
| 4 air smoke 20 s (regime) | **FAIL (cadence)** | finite ✓, 53 rises ✓, knee 106°/hip 50° of twin's 118°/71° ✓, period 0.215 s vs 2.108 s ✗ |
| 5 stretch ground eval (untuned) | ran | `−236.35`, **bilateral** 8/9 cycles, knee −79.8°, no NaN (vs s3k non-spiking: −161.57, left frozen) |

Calibration of record (`spiking_calibration.json`, produced by `D:\Anaconda\envs\myo\python.exe calibrate_spiking.py`): `sat_ref 0.9567` (measured on the analog net), `k_sn_exc 5.4087` (MN steady-V residual 0.0001 mV), `k_sn_inh 1.2284` (residual 0.0003 mV), `k_ns 0.35` and `rg_thr_inc 4.0` (RG-pair + full-net sweeps, alternation + period-target), `k_s2s_exc 1.8754` (40 Hz→40 Hz rate transfer), `k_s2s_inh 0.9996` (silences a 40 Hz-driven cell), `g_inc_readout 1.7093` (50 Hz → `E_HI`).
