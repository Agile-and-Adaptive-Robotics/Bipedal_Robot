# Goal 7 — six-synergy neuromechanical model (first pass)

**2026-09-23, ZCode (EB475WS4).** Deliverables: `synergy_model.py`,
`synergy_basis.npz`, demo `synergy_replay_air_20260923_130556.npz` (all in
`Code/MuJoCo_SNS/spinal/`). Nothing in `runner.py` / `params.py` /
`build_network.py` / `_curriculum.py` was touched — `synergy_model.py`
imports `runner` read-only for the model path, patch stack, and harness.

## What changed, in one paragraph

The six synergies of record were **found, not recomputed**: they are the
per-leg rank-6 NMF basis `fsa_backsolve.py` produced on
`bsolve_out.npz['acts']` (the `bsolve_ik.py` back-solve of subject01_walk1
IK + measured GRF; 121 frames, 0.500–2.500 s). It is now exported verbatim
to the canonical layout `synergy_basis.npz` — `W_r`/`W_l` [43 muscles × 6],
`H_r`/`H_l` [6 × 121] — with provenance baked in. Replaying the reference
gait through **only** those six channels per leg reconstructs the reference
activations at **pooled uncentered VAF 0.945 / centered R² 0.920**
(and 0.900/0.860 against the *unsmoothed* raw back-solved acts). A suspended
no-ground MuJoCo open-loop replay of the synergy drive ran a **clean 6.00 s**
(3000 steps, no non-finite states) with leg-hinge damping 0.5; open-loop
joint tracking vs the reference IK is poor (RMSE 28–58°, correlations
negative) — expected, and the honest headline limit of a first pass.

## Step 1 — the basis (existing, verified)

- **Source of "the six synergies":** `fsa_backsolve.py` selects the
  smallest NMF rank with centered VAF ≥ 0.90 **per leg** — right needed 5,
  left 6 — and uses the shared count 6 (`fsa_backsolve.py:52` `VAF_TARGET`,
  `:93-105` rank_curve; result recorded in
  `fsa_results/fsa_backsolve_report.md:5` "Selected **6 NMF synergies per
  leg**"). DESIGN.md (2026-09-17 section, ~line 871) fixes the terminology:
  these are "the six NMF synergies … S1–S6", *not* the E1/E2/F1/F2 PF
  channels.
- **Underlying data:** `bsolve_out.npz['acts']` [121 × 92], 6 Hz low-pass
  smoothed + clipped to [0,1] into `{side}_target` [121 × 43] (fsa
  `smooth_matrix`/`fit_side`, `fsa_backsolve.py:76-82,199-201`); muscle
  filter = side + `Fmax > 5` + peak > 0.05 (`:190-196`) — which keeps all
  86 active actuators (43/leg) and drops only Ben's 6 prunes
  (quad_fem/gem/peri × 2 sides).
- **Robustness (prior session, read not rerun):** held-out-frame VAF at
  rank 6 was r 0.924 ± 0.003 / l 0.912 ± 0.003 (`fsa_results/
  fsa_rank_robustness.json`).
- **What I did:** extracted `{side}_gain` → `W` and `{side}_source` → `H`
  unchanged; round-trip check in this session reproduces the published
  reconstruction exactly: stored-H replay R² = 0.9261 (r) / 0.9126 (l)
  — the same numbers as `fsa_backsolve_report.md:13,21`. Basis coverage:
  86/92 actuators; the 6 outside are exactly the prunes (verified this
  session). Per-synergy muscle support (W > 1% of column max): 19–28
  muscles per synergy — these are broad, task-level modules, not single
  muscles.

## Step 2a — six-channel replay quality (NNLS through fixed W)

Per timestep t, solve `min_{c≥0} ||W c − a(t)||²` (`scipy.optimize.nnls`),
reconstruct `Â = C·Wᵀ`, and score vs the analysis target (all computed
this session, `synergy_model.py --stage replay`):

| metric | right (43) | left (43) | pooled (86) |
|---|---|---|---|
| VAF (uncentered, synergy convention) | 0.9501 | 0.9399 | **0.9452** |
| R² (centered) | 0.9261 | 0.9126 | **0.9196** |
| VAF / R² vs **raw** (unsmoothed) acts | 0.9060 / 0.8669 | 0.8938 / 0.8522 | 0.9001 / 0.8598 |

Per functional group (primary group per `muscle_map.py`; mean over the
group's muscles, NNLS replay):

| group | n/side | VAF r | R² r | VAF l | R² l |
|---|---|---|---|---|---|
| ankle_df | 3 | 0.979 | 0.966 | 0.970 | 0.947 |
| ankle_pf | 9 | 0.935 | 0.908 | 0.929 | 0.913 |
| hip_abd | 6 | 0.831 | 0.683 | 0.890 | 0.799 |
| hip_add | 6 | 0.740 | 0.663 | 0.758 | 0.660 |
| hip_ext | 4 | 0.852 | 0.767 | 0.776 | 0.676 |
| hip_flex | 4 | 0.877 | 0.803 | 0.877 | 0.802 |
| knee_ext | 4 | 0.928 | 0.904 | 0.934 | 0.908 |
| knee_flex | 4 | 0.957 | 0.926 | 0.932 | 0.891 |
| trunk_ext | 1 | 0.884 | 0.846 | 0.784 | 0.669 |
| trunk_flex | 2 | 0.952 | 0.866 | 0.940 | 0.824 |

**Finding (measured, my inference from it below):** the per-timestep NNLS
coefficients coincide with the stored NMF H to ~1e-6 of peak (max
|C−Hᵀ|/peak = 1.3e-6 r / 1.8e-7 l; coefficient correlations 1.000 for all
12 channels). *Inference:* the converged multiplicative-update NMF was
already the per-timestep optimum through its own W, so "NNLS back-solve"
adds no reconstruction power here — it only certifies that H is the best
six-channel drive for this W. Weakest groups are hip_abd/hip_add/hip_ext —
consistent with the earlier `_synergy_pf_study.py` finding (DESIGN.md
~line 1186) that abductor/adductor tone is where 6 modules run out.
Honest caveat on the data (from DESIGN.md ~line 884 and the fsa report):
the record is ~1.6 gait cycles, one complete stride per side — these are
cycle-normalized traces, not a multi-cycle statistical average.

## Step 2b — MuJoCo open-loop demo (suspended, no ground)

Configuration (`synergy_model.py --stage demo`): model = the runner's
default `gait2392_simbody_cvt3.xml` through `runner.apply_harness(
no_ground=True, pin_rot=True)` — contact section removed (suspended rig,
pelvis spring-pinned incl. rotation; `runner.py:360-361`), start pose =
normal.mot via `runner.apply_start_pose`, WARMUP 0.5 s pose-hold with ctrl
ramp (the runner's own warmup pattern, `runner.py:1409-1415`). ctrl =
NNLS-replayed activations interpolated from the 60 Hz analysis grid onto
the 2 ms sim grid, 2.0 s clip looped; 0 ctrl on the 6 prunes. **Standing
was not attempted** — the suspended/no-ground configuration was chosen
deliberately because open-loop activation replay has no balance closure.

- **Clean run (this session):** 6.00 s / 3000 steps, no non-finite states,
  leg-hinge damping 0.5 (runner repair-2d air-stability knob,
  `runner.py:370-388`; now the demo default). Legs alternate large
  sagittal sweeps (knee_r −105..+11°, knee_l −109..+11°, COM z
  0.96–1.00 m). Output: `synergy_replay_air_20260923_130556.npz` (t, q °,
  com, ctrl, reference q, zeroed actuator list, cfg).
- **Without damping it fails:** at default damping (`--leg-damping 0`)
  the same replay goes non-finite at t = 1.786 s ("Inertia matrix too
  close to singular at DOF 66" = `rect_fem_l` patella-follower slide)
  after knees are driven to ~−105° in air — the known follower-coupler
  stress family (DESIGN.md 2026-09-10 blocker notes; follower armature
  1.0 is already applied by `patch_xml`, `runner.py:299-316`). The failed
  npz was deleted after logging; the numbers above are from the logged
  run.
- **Open-loop tracking vs reference IK** (last 2.0 s clip loop): RMSE
  27.9/58.3/43.5° (r hip/knee/ankle R), 41.8/51.6/37.6° (L); all
  correlations negative (−0.07..−0.49). *My inference:* the activations
  were solved for loaded stance under measured GRF; replayed open-loop in
  a suspended unloaded rig the same drives over-flex the knees — the
  neuromechanical loop (reflexes + ground) is doing real work that this
  first pass deliberately omits.

## Step 3 — design note: wiring the six synergies into the walker

Route (from the s3k decision-tree verdict, DESIGN.md lines 143–174 and
5–54: gating levers exhausted → "the pattern layer itself must change";
recommended option "motor primitives with antagonist groups", phase source
= the pm_* contact-reset machine; spec LIT_CIRCUIT_AUDIT.md §10):

1. **Synergy → MN-pool mapping via W.** Each `W` column is the MN-pool
   weighting of one channel. In SNS terms: one PF source population per
   channel → MN pools with W > 0. The conductance mapping V_MN = E_HI·a
   with `g = k·R·Gm/(ΔE − k·R)` (Szczecinski et al. 2017 Eq. 18, DOI
   10.3389%2Ffnbot.2017.00037 — as already implemented and audited in
   `fsa_backsolve.py` `{side}_g_fit`/`{side}_g_analytic`) converts each
   W entry to a PF→MN synapse conductance; the fsa dynamic fit (VAF
   0.863 r / 0.835 l through the actual LIF dynamics) is the correction
   for overlapping inputs at our Eexc/R ratio.
2. **Six channels as phase-windowed drives.** Cycle-fold H (heel strike
   = 0, toe-off = 50, per `gait_phase.py`) and replace each channel with
   a raised-cosine window fitted to its phase profile (the Di Russo 2023
   eq.-8 shape Ben's decision tree already specifies). The pm_* machine
   supplies phase 0–1 per side. Weight matrix = measured W rather than
   searched — this is the data-driven variant of the recommended build.
3. **What this model establishes / does not.** Establishes: six channels
   per leg carry ≥ 0.90 VAF of the reference activation structure, and
   suffice to animate the suspended musculoskeletal model without any
   rhythm generator. Does **not** establish: ground gait, balance,
   closed-loop stability, or MN-membrane dynamics (ctrl bypasses them),
   and the s3k verdict's other missing piece — the crossed
   flexor-excitatory coupling (`contra_flex`, DESIGN.md lines 43–54) — is
   orthogonal to and still required for antiphase.

**Open items for Ben:** (a) 1:1 synergy↔PF-population mapping vs shared
half-centers (DESIGN.md's phase/PF-semantics warning stands: a synergy is
not a PF neuron); (b) S5/S6 are not bilaterally phase-stable (DESIGN.md
~line 878) — keep or drop per side; (c) phase-windowed H vs measured-H
replay; (d) whether to retune on the retargeted/longer records
(`bsolve_out_retarget.npz`) before porting into SNS topology.

## How to run

```
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
C:\Users\Ben Bolen\.conda\envs\myo\python.exe synergy_model.py                # all stages (demo = damped 6 s clean)
C:\Users\Ben Bolen\.conda\envs\myo\python.exe synergy_model.py --stage replay  # numbers only
C:\Users\Ben Bolen\.conda\envs\myo\python.exe synergy_model.py --stage demo --seconds 6 --leg-damping 0   # reproduces the 1.8 s failure
```
`--stage basis` regenerates `synergy_basis.npz` from
`fsa_results/fsa_backsolve.npz` (idempotent). The demo writes a fresh
timestamped npz each run and closes its file handle; no viewer is opened.

## Honest limits

- One stride per side in the source record; no between-cycle variance,
  and looping the 2 s clip puts a seam in the demo drive (visible as the
  2.0 s periodicity in the demo log).
- Reconstruction scores are against the 6 Hz-smoothed analysis target;
  against raw back-solved acts the ceiling is ~0.90 VAF (smoothing +
  six-channel truncation together).
- The demo is deliberately open-loop at the activation level: no MN
  membrane, no afferents, no ground — it is a connectivity/feasibility
  demonstration, not a gait result.
- `synergy_basis.npz` is a *copy* of the fsa basis in a new layout; if
  `bsolve_out.npz`/fsa is ever regenerated, rerun `--stage basis` (and
  the replay/demo) — no automatic staleness detection is wired.

## Verified / not verified

**Verified (run this session, myo env, cwd `Code\MuJoCo_SNS\spinal`):**
- `synergy_model.py` (bare, all stages) → exit 0; `synergy_basis.npz`
  written; stored-H replay R² 0.9261/0.9126 == fsa report values; clean
  6.00 s damped demo; `synergy_replay_air_20260923_130556.npz` written.
- `synergy_model.py --stage replay` → exit 0; numbers in the tables above.
- `synergy_model.py --stage demo --leg-damping 0` → non-finite at
  t=1.786 s, DOF 66 = rect_fem_l follower (run twice this session;
  output deleted after logging the numbers).
- NNLS-vs-H agreement 1.3e-6/1.8e-7 of peak; coverage 86/92 (probe script).

**Not verified / not run:** ground-contact or standing replay (not
attempted by design); MN-pool SNS circuit realization of W (design note
only — wiring would need Ben's connectome-spec route per AGENTS.md);
regeneration of `bsolve_out.npz` or the fsa basis (used as-is);
multi-cycle statistical robustness (impossible with this record).
