# Goal 2 — Standing-balance stage for the gait curriculum (2026-09-23)

**What changed:** a new curriculum stage 4 ("balance") plus the conditional
circuitry it exercises — a vestibular-analog VEST pathway and a stance-gated
proprioceptive II length-loop boost — all **default-off** (byte-identity at 0
verified). Modeled on SCONE Tutorial 3a "Balance". Nothing existing changes
until you turn a knob.

**What Ben must decide:** every new edge below is listed under
["Awaiting Ben's connectome-spec confirmation"](#awaiting-bens-connectome-spec-confirmation).
Nothing is wired into your connectome spec yet; the runner's
`connectome_gains.json` hook (runner.py:858-881) already honors `gain_key`
entries, so `vest_ext` / `vest_flex_inh` can be driven from your editor once
you add rules for them.

---

## 1. The SCONE reference (what the tutorial actually does)

All files verified locally on this machine (SCONE 2.4.4 tutorial copies):

| piece | what it is | source (verified) |
|---|---|---|
| sim setup | 30 s **standing** sim, `InitStateStand.zml`, `initial_load = 1` | `C:\Users\Ben Bolen\Documents\SCONE\Tutorials\Tutorial 3a - Balance - OpenSim.scone` lines 5-14 |
| proprioceptive feedback | autogenic **length** reflexes only: `MuscleReflex { L0 = ~0.65 KL = ~1 delay = 0.010-0.035 }` on 9 muscles (iliopsoas, glut_max, rect_fem, hamstrings, vasti, bifemsh, gastroc, soleus, tib_ant) | `Tutorials\controllers\ControllerReflexBalance.scone` lines 14-23 |
| vestibular feedback | `BodyPointReflex { source = torso, KP ~0.1, KV ~0.1, delay = 0.1, offset = [0 0.5 0], direction = [1 0 0] }` → all the same 9 muscles, `symmetric = 1` | same file, lines 25-43 |
| objective | `BalanceMeasure` weight 100 (`termination_height = 0.6` = fall at 60 % of initial COM height) + `EffortMeasure` 0.01 (Wang2012) + joint-force penalties 10 on knee/hip/ankle ("penalize locked joints") | `Tutorials\measures\MeasureBalance.scone` lines 1-33 |
| robustness variant | Tutorial 3b wraps the same controller in `NoiseController { base_noise = 0.02 proportional_noise = 0.15 }` | `Tutorial 3b - Motor Noise Balance - OpenSim.scone` lines 17-27 |

SCONE reference docs (fetched 2026-09-23, scone.software `ref:balance_measure`,
`ref:body_point_reflex`): `BalanceMeasure` "checks for balance, i.e. if a
specific vertical COP position is maintained"; `termination_height` = relative
COM-height factor (default 0.5). `BodyPointReflex` = PD(KA) on a **body point**
measured along a **global direction** (targets P0/V0/A0, gains KP/KV/KA,
neuromuscular `delay`, `min/max_control_value`).

So the tutorial's balance controller = *positive autogenic length feedback
(proprioceptive) + a torso-point PD with 0.1 s delay (vestibular)*, scored by
*stay-up + effort + don't-lock-joints*.

## 2. What our model already had (read, not changed)

- Runner-side balance surrogates: BAL_PF/BAL_DF ankle-strategy COM-x PD
  (runner.py:1163-1173 pre-edit), BAL_LAT stance-gated abductor strategy
  (1174-1187), BAL_TRK_EXT/FLX trunk-lean PD ("IMU/vestibular surrogate",
  1188-1203) wired in `build_network._wire_balance` (build_network.py:851-880
  pre-edit). All remain exactly as they were.
- Ia/II/Ib afferents per muscle with stance-gated, speed-modulated
  presynaptic gains (runner afferent loop; params.AFF / params.MOD).
- Heel/toe contact mechanosensors + LBIN (v11 P1a), all default-off knobs.

## 3. Design (what I added, and why)

### 3.1 Vestibular analog — VEST cells → extensor tone (conditional topology)

- **Cells:** `VEST_r` / `VEST_l` non-spiking neurons (tau = `TAU["descend"]`
  = 0.1 s — this membrane lag stands in for the SCONE 0.1 s vestibular delay),
  built **only when** `vest_ext > 0 or vest_flex_inh > 0`
  (build_network.py:271, 284-302). Input ports `VEST_c_r/l`.
- **Signal (runner-side):** rectified tilt deviation + tilt rate,
  `u_vest = clip(kp_trk·|lean − trk_ref| + kd_trk·|lean_dot|, 0, max_trk)`
  (runner.py:1240-1262) — reuses the Ben-tuned BAL_TRK gains and the existing
  torso up-vector lean. Bilateral symmetric (SCONE `symmetric = 1`).
  *Inference (mine):* rectified magnitude = an antigravity-stiffness pathway;
  the **directional** correction deliberately stays with the existing
  BAL_TRK/BAL_PF/BAL_DF controllers instead of duplicating them.
- **Edges:** `VEST_side → MN` direct, ipsilateral, same pattern as the
  existing BAL_* balance cells (build_network.py:901-927):
  - extensor groups `knee_ext, ankle_pf, hip_ext, trunk_ext` ← excitatory,
    gain `G["vest_ext"]` (21 muscles/side → 42 edges);
  - flexor groups `hip_flex, knee_flex, ankle_df, trunk_flex` ← inhibitory,
    gain `G["vest_flex_inh"]` (13 muscles/side → 26 edges; each edge gated by
    its own gain > 0).

Literature tracing for these edges:
1. **Ben's request text** (the ask): "vestibular analog = pelvis-tilt / COM
   sensors driving extensor tone + ankle strategy".
2. **SCONE Tutorial 3a** (verified): vestibular BodyPointReflex from torso
   drives all major leg muscles with KP/KV (lines 25-43).
3. **Di Russo, Ijspeert & Bouri 2023 JNE** (local PDF `spinal\_dirusso_2023.pdf`,
   text extracted this session): their trunk balance controller is the PD
   `ubalance = kp·(θ(t−tD) − θ0) + kv·θ̇(t−tD)` (**eq. 5**) on the forward-lean
   angle with neural delay, applied to hip MNs (ILPSO/GMAX/HAMS) during
   stance, summed outside the MN dynamics (eq. 4). Our VEST current has the
   same PD form (deviation + rate, reference angle `trk_ref`, delay via tau).
4. **Vestibulospinal → extensor (antigravity) tone:** *standard textbook
   physiology, NOT full-text verified this session* (canonical refs to check:
   Wilson & Yoshida 1969 J Neurophysiol; Brodal 1981). Labeled as inference;
   the edge authority is items 1-2.
5. **Flexor inhibition (vest_flex_inh):** LVST reciprocal flexor inhibition —
   same textbook status, flagged as the most speculative edge.

### 3.2 Proprioceptive balance — stance-gated II length loop (presynaptic)

`G["vest_prop"] > 0` scales the stance-gated **II length component** of the
afferent encode by `(1 + vest_prop·stance)` (runner.py:1170-1185, flag at
1132). This is the SCONE `KL·[L−L0]+` analog. It is a presynaptic gain, not a
new edge, matching the documented pattern "phase gating of afferent gains is
implemented presynaptically in the runner" (params.py:19-21).
*Decision (mine, flagged):* the ask said "Ia, possibly II" — I implemented II
only: Ia already carries stance-gated, speed-modulated homonymous + reciprocal
paths (params.MOD), and a second Ia boost would double-count the same loop.
Say the word and a `vest_ia` knob slots into the same four places.

### 3.3 Standing-balance eval — `--stand-eval` + bal_* metrics

- `runner --stand-eval T` (default 8 s): eval-mode run with a standing-only
  schedule (all windows collapsed; DRIVE 0 throughout) (runner.py:855-870).
- New eval metrics (runner.py:1591-1633), computed over the post-warmup
  prefix: `bal_sway` / `bal_sway_rms` (peak / RMS radial COM-x,y deviation
  from the ankle-axis reference), `bal_tilt_max` (|pelvis_tilt| envelope,
  deg), `bal_contact_sym` (mean load R/(R+L); 0.5 = even), `bal_com_z_min`,
  `bal_fell` (COM height < 0.55 m — the summary's existing fall criterion;
  SCONE's analog is `termination_height`).
- Runner flags, all default-off: `--vest X`, `--vest-flexinh X`,
  `--vest-prop X`, `--stand-eval T` (runner.py:871-889).

### 3.4 Curriculum stage 4 (balance)

`python _curriculum.py 4 [n_trials]` (study `curr_s4_balance`, sqlite db as
the other stages; **I did not launch it**):

- searched knobs (`KEYS4`, _curriculum.py:64): `vest_ext` [0, 0.5],
  `vest_flex_inh` [0, 0.3], `vest_prop` [0, 1.0] (small-entry ranges, the
  09-18 lesson) + `rig_scale` [0.05, 1.0].
- objective (_curriculum.py:262-278): 8 s standing eval at the trial's rig
  wean; `100 − 400·sway[m] − 1·tilt[deg] − 40·|0.5 − contact_sym|`
  (symmetry deviation from an even R/L split; 20 max cost);
  sentinels NaN −200 < fall −150 < any stander.
  **Correction (supervisor round 2, 2026-09-23):** the first version used
  `−40·(0.5 − sym)` — asymmetric: right-side overload (sym > 0.5) yielded a
  score BONUS up to +20. Now `abs(0.5 − sym)` — either side overloading costs
  the same. Defaults-off behavior unaffected (the term only executes inside
  stage-4 trials).
- seeding: vest gains start at 0, `rig_scale` starts at 1.0 → trial 0 is
  today's incumbent standing reproduced exactly.
- `set_stage(4, ·)` loader at _curriculum.py:125-130; stage-2/3 blocks are
  now guarded `2 <= stage <= 3` / `== 3` so stage 4 doesn't zero gait knobs.

**Why rig_scale is searched (my call, flagged):** at rig 1.0 the pelvis
springs do the balancing and sway is ~0 regardless of the controllers — the
stage would look converged while measuring nothing. Weaning the rig is the
point: AGENTS.md records the support boundary at S≈0.8-1.0 "until the
pelvis-balance piece exists" — this stage is that piece.

### 3.5 JSON RULE — all four places, per knob

| knob | params.G default | runner loader branch | winner-json save path | knob surface |
|---|---|---|---|---|
| `vest_ext` | params.py:285 (0.0) | runner.py:772-778 tuple | trial params → `curriculum_stage4.json` (whole-`best.params` save, _curriculum.py:394-398) | `--vest` flag; suggest at _curriculum.py:213 |
| `vest_flex_inh` | params.py:286 (0.0) | same tuple | same | `--vest-flexinh`; :214 |
| `vest_prop` | params.py:287 (0.0) | same tuple | same | `--vest-prop`; :216 |

(`rig_scale` is a runner arg, not a G knob; it rides the trial params the
same way `pelvis_ty`/`ky_scale` do.) `RUNNER_DUMP_STATE` dumps the whole `G`,
so the new keys appear in state dumps automatically.

## 4. How to run

```bat
:: from Code\MuJoCo_SNS\spinal, myo env; PROTECT spinal_run.npz:
set AARL_NPZ=goal2_stand.npz
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" runner.py --stand-eval 8 --vest 0.5 --vest-flexinh 0.2 --vest-prop 0.5

:: defaults regression (identical to pre-goal2 behavior — no flags):
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" runner.py --fitted --best10 --eval --drive 2.929

:: the balance stage (when you want it — I did NOT launch this):
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" _curriculum.py 4 30
```

## 5. What was verified (commands I ran, outputs I saw)

1. **py_compile** all four touched files → exit 0
   (`python -m py_compile params.py build_network.py runner.py _curriculum.py`).
2. **Defaults-off config diff (the ask's check):** dumped the whole params
   module (G/TAU/BAL/AFF/MOD/weight tables/SCHEDULE) before and after the
   edits; `fc` shows the ONLY difference is three added keys at 0.0:
   `"vest_ext": 0.0, "vest_flex_inh": 0.0, "vest_prop": 0.0`.
3. **Topology fingerprint at defaults:** built the full network before and
   after (92 muscles, interleg on): **identical** — 410 neurons / 376 inputs /
   1186 nonzero synapses / Σg = 668.17 both times.
4. **vest-ON build check:** with `vest_ext=0.5, vest_flex_inh=0.2`:
   412 neurons / 378 inputs / 1254 synapses / Σg = 694.37 — exactly
   +2 cells, +2 ports, +68 edges (42 extensor × 0.5 + 26 flexor × 0.2 =
   +26.20 conductance sum, arithmetic matches), ports `VEST_c_r/l` present.
5. **One short smoke (2 s standing, all knobs ON)** through the exact stage-4
   code path (`set_stage(4, p)` then `R.main(["--stand-eval","2",...])`):
   exit 0, no NaN, no fall. Metrics returned:
   `bal_sway 0.1205 m, bal_sway_rms 0.0938 m, bal_tilt_max 26.94°,
   bal_contact_sym 0.271, bal_com_z_min 0.757 m, bal_fell false`
   (tilt 26.9° is consistent with the standing runs on record, tilt 27-30°;
   asymmetry reflects the normal.mot asymmetric start pose: mean contact
   R 901 N vs L 2426 N).
6. **State dump from the smoke** (`RUNNER_DUMP_STATE`): `vest_ext 0.5,
   vest_flex_inh 0.2, vest_prop 0.5` — flags land in `G` before the build.
7. **Output hygiene:** run written to `spinal\goal2_smoke_balance.npz`
   (unique name; valid npz, 1000 steps, opened + closed cleanly);
   `spinal_run.npz` untouched (mtime 09/22/2026 17:15, before this session).
   No Dissertation tex/figures touched. No optuna study created.
8. **Round-2 fix (supervisor BLOCK):** stage-4 symmetry term corrected to
   `−40·|0.5 − bal_contact_sym|` (_curriculum.py:277-278, was
   `−40·(0.5 − sym)` which bonused sym > 0.5); comment updated to "20 max
   cost". py_compile re-run → exit 0; line-content check via findstr;
   arithmetic check run (see below). Defaults-off unaffected — the edit
   touches only the stage-4 objective branch, no params/G/network code.

### Verified / not verified

**Verified:** everything in the list above (compile, defaults-off params diff,
defaults-off topology fingerprint, vest-on topology arithmetic, 2 s smoke end
to end incl. metrics + state dump, output-name isolation).

**Not verified (honest gaps):**
- **No optuna study was launched** (forbidden by the task) — stage 4's
  *optimization* behavior (TPE over the 9-dim space, seed trial, sentinel
  ordering in a real study) is untested; only the single-trial code path was
  exercised.
- The 8 s full-duration stage eval was not run (smoke was 2 s, per "seconds
  of sim time").
- A bit-exact regression re-score of a recorded winner (e.g. `--best10`) was
  **not** run; the defaults-off proof is the config diff + topology
  fingerprint named by the task as the acceptable short check.
- The vestibulospinal-extensor citation is textbook-level, not full-text
  verified this session (§3.1 item 5).
- Knee/ankle-joint "locked-joint" force penalties (SCONE's third objective
  term) are not in the stage-4 objective — the RoM caps already enforce
  limits; add a term if you want the SCONE-equivalent.

## 6. Awaiting Ben's connectome-spec confirmation

All default-off; nothing below affects any existing run until you set a gain
(or add it to `connectome_gains.json` as a rule with the matching `gain_key`):

1. **VEST_r / VEST_l cells + `VEST_c_*` input ports** (brainstem-surrogate
   family, tau 0.1 s) — the vestibular analog itself.
2. **VEST → ipsilateral knee_ext / ankle_pf / hip_ext / trunk_ext MNs**
   (excitatory, `vest_ext`) — "extensor tone" per your request + SCONE T3a;
   direct-MN wiring following the BAL_* precedent.
3. **VEST → ipsilateral hip_flex / knee_flex / ankle_df / trunk_flex MNs**
   (inhibitory, `vest_flex_inh`) — LVST reciprocal flexor inhibition,
   textbook-labeled, most speculative.
4. **Rectified tilt-deviation + tilt-rate runner signal** (otolith/canal
   analog; no COM-x term — that stays with BAL_PF/DF). Inference (mine).
5. **`vest_prop` stance-gated II length-loop boost** (presynaptic; SCONE KL
   analog); and the **decision to leave Ia untouched** (§3.2).
6. **Stage-4 objective weights + 8 s duration + sentinels** (§3.4) and
   **searching `rig_scale` in the stage** (§3.4, my call).

— compiled 2026-09-23, EB475WS4, myo env; runner defaults regression-gated by
config diff + topology fingerprint this session.
