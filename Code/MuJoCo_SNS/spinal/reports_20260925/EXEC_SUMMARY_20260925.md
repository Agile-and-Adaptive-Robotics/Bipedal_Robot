# EXECUTIVE SUMMARY — 2026-09-25 unattended campaign (goals 1–4)

**Written:** 2026-09-25 ~21:30, machine EB475WS4. This file indexes the whole campaign; every claim below is cited to the goal report it comes from (all in this folder), and the "what I verified myself" notes at the bottom say exactly which checks this summary re-ran vs. which are attributed to the session that wrote each report.

---

## TLDR

The campaign ran four goals in one day (12:55–21:15) and delivered three wins, one plan, and one honestly-diagnosed wall. **Goal 1 (SCONE):** first-ever SCONE evaluations and CMA-ES optimizations on this machine — all 15 tutorial scenarios ran on Hyfydy with zero OpenSim fallbacks, and optimization produced a balance controller that stands the full 30 s horizon (fitness 96.6 → 2.21, BalanceMeasure 0) plus the campaign's key methodological lesson: a short-window gait optimization *overfits* (improves 4-s fitness but falls at 4.7 s) while full-horizon optimization keeps a 20-s/36-step walk. **Goal 2 (AnimatLab walker → MuJoCo-SNS port):** the body transported exactly (pose dev 4.5e-7 m, mass dev 0.0 kg) and, after M3 found and fixed a real M1 defect (knee/ankle hinge axes were VERTICAL), the neural stack is proven end-to-end in air — including Ben's split-RG ask, with an ablation experiment proving each leg runs on its OWN RG (1.027 s free vs 1.357 s locked) — but ground walking does NOT walk, and the failure is now reduced to three measured causes, the decisive one being structural: **the rest pose's COM sits 3.09 cm outside the support polygon** (the source AnimatLab model always froze the pelvis). **Goal 3** delivered the data-driven test battery (T1–T8) that catches fragile walkers *before* any BPA bracket is 3D-printed — plan only, nothing executed. **Goal 4** built two new network architectures on the gait2392 body (`w2lvar` = W2L layout transplant, `syn6` = data-driven 6-synergy) behind a default-inert selector (defaults gate (410, 376, 1186) PASS throughout), wired Ben's 5-stage curricula for both, and ran BOTH chains to completion — 10/10 studies, 202 trials, every winner bit-exactly replayable, zero exploits — with both variants producing genuine ground walks (w2lvar kine −165.0 upright at kz 0.83; syn6 −197.2 with real foot loading but a ~100 %-planted left leg). Supervisor verdict on goal 4: **GO**. Protected artifacts (`spinal_run.npz`, `optuna_walk.db`, prior reports) untouched throughout — re-verified by stat during this write-up.

---

## Status table

| Goal | Outcome | Report(s) |
|---|---|---|
| 1 — SCONE tutorials + first optimizations | **PASSED** (15/15 evals Hyfydy exit 0; 5 optimizations; bit-exact re-evals) | `goal1_scone_tutorials_and_optimization.md` |
| 2 — M1 physical body port | **Complete** (12/12 gate PASS) | `goal2_m1_physical_model.md` |
| 2 — M2 Li neural architecture | **Partial** — closed loop runs; stepping gate "not yet" (body-side collapse) | `goal2_m2_li_architecture.md` |
| 2 — M3 2023-original air stepping | **Complete — GATE PASS** (+ fixed M1's vertical knee/ankle axes) | `goal2_m3_w2l_cpg_air.md` |
| 2 — M4 split-RG (Ben's ask) | **Complete — all 3 gates PASS** (incl. causality ablation) | `goal2_m4_rg_split.md` |
| 2 — M5 afferent feedback (Ben's rules) | **Partial** — gates a+b PASS, gate c (toe DF suppression) FAIL-as-specified | `goal2_m5_afferents.md` |
| 2 — M6 ground stand + walk attempt | **Partial** — rigged stand holds (27 % weight); unrigged fails structurally; walking = fall/march, 3 blockers named | `goal2_m6_stand_ground.md` |
| 3 — Data-driven test plan | **Plan delivered** (no runs — by design) | `goal3_data_testing_plan.md` |
| 4 — w2lvar build | **Complete** — gates a/b/c PASS | `goal4_build_w2lvar.md` |
| 4 — syn6 build | **Partial** — gates a/b + network-only rhythm PASS; full-runner air rhythm PARTIAL | `goal4_build_syn6.md` |
| 4 — curricula wiring | **Complete** — both 5-stage drivers smoked, defaults gate PASS | `goal4_wiring.md` |
| 4 — w2lvar chain | **Complete 5/5**, all exploit checks bit-exact | `goal4_chain_w2lvar.md` |
| 4 — syn6 chain | **Complete 5/5** (stage 5 won by the seed) | `goal4_chain_syn6.md` |
| 4 — final write-up + audit | **Complete — supervisor GO** | `goal4_variant_curricula.md` |

---

## Per-goal detail

### Goal 1 — SCONE (Hyfydy engine, `sconecmd`, seed 42 everywhere)

- **15 tutorial evaluations, all never-run scenarios** (Tutorials 1, 2a–2c, 3b, 4b–4d, 5a–5c, 6a–6d), all exit 0 on Hyfydy, **no OpenSim fallback fired anywhere**. Example numbers: 3b motor-noise balance falls at 1.045 s (same as the 3a baseline); 4b fast gait 0.180 m/s, 11 steps, falls; 5a PF-weakness is the best default-param walker of the pathology set (0.197 m/s, 2.56 m). Logs `logs\tut*.log`, motions `logs\motion_*.sto`.
- **First optimizations ever on this machine** (5 result dirs under `Documents\SCONE\results\260925.*.R42\`):
  - **Balance 3a, 300 gens, D=12:** fitness 96.56 → **2.207** (monotone, 56 new-best generations); the best par **stands the full 30-s default horizon** with BalanceMeasure 0 — CMA-ES found the vestibular-like torso gains (e.g. `Balance.iliopsoas-torso.KV = 0.727`). The ask's 40-gen budget alone only reaches "falls at 2.8 s" — the capability lives in gens 54–293.
  - **Gait 4a headline (the trap):** the 4-s-window D4 optimization improved short-window fitness (0.901 → 0.858) but its par **falls at 4.7 s**, while the pretrained par walks 20 s / 36 steps / 21.5 m. The bonus D20 full-horizon run keeps robustness (20 s / 36 steps / 21.1 m) and still improves the true objective (0.883 → 0.873) with 5.5× lower GRF penalty. **Recommendation banked: always optimize gait at the full horizon (or ≥ 2× expected fall time).**
  - **High Jump 2a:** +13.5 % over default (93.6 cm), **still climbing at the 20-gen cap** — a lower bound, not an optimum (resume item below).
  - All 5 optimizer best-pars reproduced their recorded fitness **bit-exactly** in re-evaluation (`logs\reeval_*.log`).

### Goal 2 — Walker_2_Layer_CPG → MuJoCo-SNS port (6 milestones)

- **M1 (body):** 13 bodies, 8 hinges, 12 muscle actuators + 2 toe springs; gate 12/12 PASS (`logs\validate_body_run6.log`). Proof numbers: **pose dev 4.5e-7 m, mass dev 0.0 kg (total 41.943 kg), joint-range dev ≤2e-10 rad**. Documented exclusions: Kse/Kpe/muscle-damping B not portable to MuJoCo 2.3.7 `<muscle>`; toe-spring damping capped 20000 → 172 (source value ~116× critical, un-integrable).
- **THE M3 FINDING (reframes M1/M2):** M1 shipped knee+ankle hinge axes VERTICAL (world −z; e.g. knee_L `(−0.052, 0.000, −0.999)`) — zero sagittal moment arm, knees could produce no flexion torque. Root cause: parent-relative vs world rotation in `make_w2l_mjcf.py`. `fix_joint_axes.py` → `w2l_mjcf_fixed.xml`, verification **8/8 lateral axes PASS**; post-fix moment arms are real (knee ext 0.96 mm/°, hip flx 1.73 mm/°). Standing lesson banked: **axis directions need a gate.** Also corrected the M1 record: the shipped Root has NO joint (welded pelvis), not a free joint — M2's "pelvis tunneling to −1.44 m" was actually the ankle angle.
- **M2 (Li architecture):** reference data re-generated from Li's own asim (period 1.305 s / 0.77 Hz, duty ≈ 0.5, antiphase, +0.64 m/s, never falls). The transcription runs end-to-end (56 neurons / 84 synapses / 12 MV relays) but the stepping gate is **"not yet"**: 4 runs, best still collapses (min pelvis −0.84 m, 0 stance episodes). Diagnosis is body-side: rest pose hovers 2.6–3.1 cm, no muscle damping B, rigid tendons — largely the M3 axis artifact plus the documented M1 exclusions.
- **M3 (2023-original air stepping): GATE PASS** — 20 s clean air stepping, 0.97 Hz (within 2× of the 0.77 Hz 2023 reference), antiphase r −0.623, hip range 41.8° vs ref 38°, zero ground contacts. Ankle is the least faithful joint (range +75 % over ref; transported range [−20°,−5°] excludes the rest pose — flagged).
- **M4 (split-RG, Ben's ask "disconnect the LH RG, mirror it, couple with commissurals"): all three gates PASS.** Coupled: period locks 1.357/1.358 s both sides, **L/R phase lag 0.507 cycle** (band-limited antiphase r −0.941). Ablated (`--comm=0`, commissurals not built — conditional topology): **both legs still oscillate at their intrinsic 1.027 s** — the causal proof each leg is powered by its own RG. M3 default path regression: PASS, identical numbers.
- **M5 (afferents per Ben's rules):** heel = stance reset at the PF layer (g 0.5), toe = DF inhibition (1.0/5.0/2.749), Ib load (0.5) — census 97n/188s (M4 exactly when off). **Gate (a) PASS** (affermilated 20 s air walk, period 1.389 s, Ib tracks extensor force r 0.97/0.99). **Gate (b) PASS** — causal heel reset: +54 ms first-onset shift, settling to a persistent +10–12 ms offset; sign flips with pulse phase (PPR signature; disclosed post-hoc pass bar, max|shift| over first 4 onsets ≥ 40 ms). **Gate (c) FAIL-as-specified:** the toe chain is built and correctly signed but the DF pool never excites in air, so ≥10 % suppression is **NOT DEMONSTRATED** — resolution belongs to M6/ground work (the dedicated DF half-center Ben drew doesn't exist in this 2023-template net).
- **M6 (ground):** unrigged stand **fails for a structural, measured reason** — **COM x = −3.4544 m sits 3.09 cm in front of the support polygon's front edge** (`tmp\m6_static_margin.py`; every free run topples forward at ~0.9 s, heel force decays 398→0 N). Rigged stand (S=1.0: kxy 2000 N/m, kz 2000 N/m, krot 400 N·m/rad) holds 11 s with 27 % weight carried — honestly a harness hold, feet unloaded (heel duty 0.00). **Ground walking does NOT walk:** free = fall at 0.89 s; supported = upright in-place march at the autonomous RG period (identical 14/14 bursts at 1.357 s, r −0.751 as the M4 air gate) with **zero heel loading in 19 s (heel SN max 0.00 mV)** — the contact-driven layer built in M5 stays silent on the ground. **Top 3 blockers, each with evidence:** (1) COM outside the polygon (structural; fix = standing-pose solve / balance layer); (2) ankle/foot cannot hold stance (range excludes rest pose → limit slam; PF latch lifts heels in 100 ms; muscle damping B unported); (3) the RG provably ignores the ground until 1+2 are fixed.

### Goal 3 — Data-driven test plan (plan only; explicitly no runs)

- Battery **T1–T8** with priority order P1–P5; P1–P3 need **zero downloads** and are the pre-print gate. Key acceptance anchors: T1 no-regression min ≥ −235 over the 12 walking refs (s3k measured worst −224); T2 period error ≤ 25 %/ref + monotone score-vs-speed (Ong 0.5–2.0 m/s family, T 2.10→1.02 s); T4b push battery (27 runs) with the falsifiable prediction **ankle fails first at Fmax-scale ≈ 0.7** (torque budget ankle 1.4× is tightest); T7 basin gate 12/12 and T8 rig weaning S ≥ 0.8 mandatory for every finalist. Finalist (pre-3D-print) gate = T1 ✓ T2 ✓ T7 ✓ T8(S≤0.6) ✓ T4b ✓.
- Verified library facts: 43 local refs (subject01 anchor −189.5, Falisse, 10 Ong speeds, 32 Arnold running = OOD stress only, never tuned to); **no GRF waveforms and no EMG in any npz**; one cycle per ref (point comparisons only). Staged-not-integrated: Hamner raw EMG (10 subj × 4 trials), results-deficits, SCONE setups.
- Honest gaps stated in the plan: Wang perturbation direction unverified until download; `ong_selfsel_Init200` degenerate (exclude); `gait_lib_score_all.py:106` report-string has an inverted "lower is better" docstring bug (ranking confirms higher-is-better); v10-era fitted winners no longer walk under current physics — always test the *current* champion (s3k family).

### Goal 4 — Two new architectures + Ben's 5-stage curricula

- **w2lvar** (`build_network_w2lvar.py`, W2L layout on the 92-actuator body): 888n / 382 in / 7430 syn; gates a/b/c all PASS (gate c after two documented fixes: per-muscle afferent division fixed a measured co-latch; Shevtsova rebalance fixed in-phase stepping → final E/F antiphase r −0.910/−0.916, period 1.80 s).
- **syn6** (`build_network_syn6.py`, 6-synergy layers from `synergy_basis.npz` via the audited Szczecinski Eq-18 conductance): 794n / 382 in / 3033 syn; gates a/b PASS, network-only rhythm PASS (period 2.133 s, E/F corr −0.989, S5 distinct after the mixed-channel mitigation), **full-runner air rhythm PARTIAL** (slow ~2–3 s events, interleave `FEFEFFE`, dies after ~2 cycles) — two E-latch mechanisms found and gated off (`syn6_brainstem` default 0; Ib→LBIN gated on `ib_rge>0`); residual slowness documented as a tuning surface.
- **Wiring** (`_curriculum_w2lvar.py` / `_curriculum_syn6.py`): self-pinning env, chain-private dbs (`optuna_w2lvar.db` / `optuna_syn6.db`) and npz; Ben's 5 stages (1 air-deaff, 2 air-aff, 3 balance, 4 walk-nocontact, 5 walk-contact); full-dict seeds, stock sentinel discipline; searches restricted to keys each variant actually consumes (regex-verified). Smokes: switch selects the right net (npz `neuro_names` proof), objective reads the right columns, stages 2/3/5 routes line-verified but **not separately executed** (labeled not run).
- **Chain results (both 5/5, 202 trials total, all COMPLETE, every winner bit-exact replay):**

| stage | w2lvar winner | syn6 winner |
|---|---|---|
| 1 air deaff | **80.836** (1.43 Hz, knee −78°) | **95.639** (1.57 Hz, knee −95°) |
| 2 air aff | **60.230** (reflex gains pushed DOWN in air) | **86.751** (heel/toe ≈ 0 — "tolerated, not additive") |
| 3 balance | **25.077** (winner = a random sample; under-explored) | **37.463** (heel 0.42 / toe 0.48 — contact ports win once feet load; rig 0.946 barely weaned) |
| 4 walk nocontact | **−225.891** (2/23 real cycles — air pattern-matching is hard) | **−237.253** (real but slow T 2.45 s and IN-PHASE legs, lag 0.0) |
| 5 walk contact | **−165.041** (zero penalties, kz 0.83, duty 0.46, tilt 27.7°; TPE EARNED 4 nonzero contact gains, drive 1.68→3.67) | **−197.215** (SEED won — 21 TPE trials couldn't improve; real loading contact_frac 0.955/1.0 but left leg planted ~100 %) |

- **Honest framing (per the 2026-09-25 correction in AGENTS.md):** the variants take Ben's rule files as the *dress* over two architecture families with documented deviations (D1–D6 / D1–D3) — they are **NOT** a rebuild of Ben's per-micro-layer connectome drawing and must not be described as "following his connectome rules." Variant kine scores (−165/−197) land numerically in the s3k band (−189.5…−241) but are **context only** — no goal-5-style multi-reference scoring has been run for the variants; no head-to-head claim is made. Neither winner is basin-gated; cadence/duty/phase gaps are the same architectural family the stock chain documented.

---

## What to resume tomorrow

1. **Goal-4 chains are COMPLETE (10/10 studies) — the resume commands are for extension, not completion** (each script self-pins env/db/npz; re-running a stage ADDS trials via `load_if_exists=True`):
   ```
   cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
   C:\Users\Ben Bolen\.conda\envs\myo\python.exe _curriculum_w2lvar.py 3 18   # re-tune the under-explored standing map
   C:\Users\Ben Bolen\.conda\envs\myo\python.exe _curriculum_syn6.py 5 22     # seed won; give TPE more shots (consider re-seeding stage 4 from stage 2)
   C:\Users\Ben Bolen\.conda\envs\myo\python.exe _curriculum_w2lvar.py 5 22   # stage-5 continuation (sentinel-heavy landscape)
   ```
   Then: multi-reference-score both variant winners against subject01/Falisse/Ong (goal-3 T1) for the first true head-to-head with s3k; adapt `basin_gate.py` to the variants' dump states and gate all four walking winners.
2. **SCONE follow-ups** (`goal1` §2/§5; scratch dirs `%TEMP%` are disposable — re-create with `tmp\setup_scratch.ps1` if gone): High Jump longer run (still climbing at the 20-gen cap, ~0.02 s/gen):
   ```
   "C:\Program Files\SCONE\bin\sconecmd.exe" -o "%TEMP%\scone_opt_gait\Tutorial 2a - High Jump - Hyfydy.scone" CmaOptimizer.max_generations=100 CmaOptimizer.lambda=16 CmaOptimizer.random_seed=42
   ```
   And per the banked recommendation: any future gait optimization at the **full horizon** (`CmaOptimizer.SimulationObjective.max_duration=20`).
3. **Goal-2 port, in the M6 report's own order:** (1) standing-pose solve (M1 open item #4) that puts the COM inside the support polygon; (2) ankle-range surgery + muscle-damping B in the generator (not runtime stand-ins); (3) re-run `test_w2l_ground.py` **unchanged** — its contact encoders will then have heel signals to reset on. After the body fixes, re-run the M2 Li gate (`test_li_stepping.py`, knobs `--tonic-scale/--mu-gain/--contact-gain-na/--hip-sign`) toward its still-unmet 20-s stepping gate; M5's toe-suppression gate stays deferred until real toe loading exists (the dedicated DF half-center is Ben's architecture call).
4. **Goal-3 battery, P1 first (zero downloads):** extend `gait_lib_score_all.py` to argv + per-term CSV (2–4 h), then T2 speed sweep; add default-neutral `--push`/`--fmax-scale` runner flags with the regression gate (standing −160.23425729850192 + counts (410, 376, 1186) + `_pf_layer_variant_test.py`); T4a needs the **Wang & van den Bogert processed-only bundle (264 MB, Zenodo 3819630, CC BY 4.0)**. Push grids and LOFO retraining belong on easteregg2 or run overnight here capped.
5. **Stock-chain standing item:** stage-3 (balance) tuning surface — the plan notes stage 3 has zero human balance references until the Wang download lands.

---

## Supervisor verdict

**GO** — "Audited the goal-4 variant curricula (w2lvar + syn6) before write-up; verdict GO — the winners are real and reproducible" (`goal4_variant_curricula.md` §4; audit scripts + logs in `audit\` and `logs\`). Findings: (1) read-only db census of all three dbs matches every claim (10 studies, exact COMPLETE counts 20/18/18/23/22 ×2, winners = argmax, jsons match db bit-exactly, `optuna_walk.db` untouched with zero name leaks); (2) no sentinel/exploit winners — the 2026-09-20 static-pose exploit pattern is absent everywhere; (3) defaults gate re-run PASS (410/376/1186) with `AARL_NET` unset and `syn6`/`syn6_brainstem` default 0.0; (4) spot replays of one winner per variant bit-exact (w2lvar s1 80.83558669459252, syn6 s1 95.6388239983187); (5) full per-number cross-check of both chain logs against db/logs, all 10 exit codes 0. Minor notes referenced by the verdict were not transmitted and none blocks the write-up. **Scope note:** this GO covers the goal-4 chains; the reports folder contains no separate campaign-wide supervisor verdict for goals 1–3 — their status rests on the per-report gates cited above.

---

## What this summary verified itself (vs. attributed)

Re-run during this write-up (~21:18–21:30): (a) `Get-ChildItem` stat of the campaign artifacts — `spinal_run.npz` 8,238,797 B and `optuna_walk.db` 4,403,200 B, **both still at the pre-campaign 2026-09-25 01:10:02 AM vintage** (protected-set claim confirmed first-hand); both variant dbs + all 10 winner jsons + both curriculum drivers + both builders + `w2l_mujoco\test_w2l_ground.py`/`fix_joint_axes.py` present with mtimes consistent with the chain timeline; (b) `audit\defaults_gate_output.log` read: `counts neurons/inputs/synapses = (410, 376, 1186) / GATE: PASS / syn6 = 0.0 / syn6_brainstem = 0.0`; (c) `logs\writer_replay_syn6_s1.log` read: `RECOMPUTED 95.6388239983187 vs json 95.6388239983187 delta +0.000000 (0.0000%) / REPLAY: PASS`, followed by the documented post-verdict WinError 32 housekeeping crash (PASS verdict unaffected). Everything else in this summary is attributed to the goal reports listed in the status table; I did not re-run any simulation, optimization, or study for this index.

### Artifact index (quick pointers)

- Campaign logs: `logs\` (goal-1 `tut*/opt_*/reeval_*`, goal-2 `test_*/m5_*/m6 probes in tmp\`, goal-4 `smoke_*/curr_*/harvest_*` + `*_exploit_check.log`), audit: `audit\` (db census, defaults gate, replay scripts + npz evidence).
- New code: `spinal\build_network_w2lvar.py`, `spinal\build_network_syn6.py`, `spinal\_curriculum_w2lvar.py`, `spinal\_curriculum_syn6.py`, `spinal\w2l_mujoco\` (M1–M6: parser, validators, builders, gates, `fix_joint_axes.py`), plan survey in `tmp\ref_survey.py`.
- Winner jsons: `spinal\curriculum_{w2lvar,syn6}_stage{1..5}.json`; dbs `spinal\optuna_{w2lvar,syn6}.db`; chain npz `spinal\spinal_run_{w2lvar,syn6}.npz`.
- SCONE results (5 self-contained dirs): `C:\Users\Ben Bolen\Documents\SCONE\results\260925.{130052,130211,130348,130428,130924}.*.R42\`.
