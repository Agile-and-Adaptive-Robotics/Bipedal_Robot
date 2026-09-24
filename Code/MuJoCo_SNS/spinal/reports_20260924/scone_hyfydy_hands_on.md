# SCONE 2.4.4 + Hyfydy hands-on — EB475WS4, 2026-09-24

All runs executed this session with `"C:\Program Files\SCONE\bin\sconecmd.exe"` (SCONE
2.4.4.3333; every invocation auto-initializes `OpenSim-3.3-2021-01-28` and `OpenSim-4.4`).
Scenario copies used: `C:\Users\Ben Bolen\Documents\SCONE\Tutorials3\` (the user-writable
copies made tonight; the Program Files set is read-only). sconestudio.exe was NOT touched
and the settings zml was NOT edited. No license keys are reproduced in this report.

Raw logs: `reports_20260924\scone_logs\` (this folder). Full `.sto` motion files + the
`.sto` parser stay OUTSIDE the repo at `D:\temp\scone_hands_on_20260924\`
(`sto_summary.py` = first/last/min/max of named columns; repo stays lean per the
size-reduction objective).

---

## 1. Tutorial 3a — Balance — OpenSim (headless, short duration)

Scenario: `Tutorial 3a - Balance - OpenSim.scone` — `CmaOptimizer > SimulationObjective`,
`max_duration = 30`, `ModelOpenSim3 { model_file = models/H0918v3.osim, state_init_file =
init/InitStateStand.zml, initial_load = 1, fixed_control_step_size = 0.005 }`, controller
include `controllers/ControllerReflexBalance.scone`, measure include
`measures/MeasureBalance.scone` (BalanceMeasure `termination_height = 0.6` weight 100 +
EffortMeasure Wang2012 weight 0.01 + locked-knee/hip/ankle DofMeasures).

Commands run (evaluate = one sim at the scenario's DEFAULT parameter values, i.e. the
`~mean` of every `~mean~std<min,max>` parameter):

```
sconecmd -e "...\Tutorial 3a - Balance - OpenSim.scone" SimulationObjective.max_duration=3 -r "D:\temp\scone_hands_on_20260924\balance_osim" -l 2
sconecmd -e "...\Tutorial 3a - Balance - OpenSim.scone" CmaOptimizer.SimulationObjective.max_duration=1 -r "D:\temp\scone_hands_on_20260924\balance_osim_md1" -l 2
```

Result (first command, log `scone_logs\balance_osim.log`):

```
Created model H0918v3; dofs=9 muscles=18 mass=74.5314
Terminating simulation at 1.050
result           = 102.92
  BalanceMeasure = 96.5  <- 100 * 0.965
  Effort         = 6.42043 <- 0.01 * 642.043
  JointLimits    = 0
simulation time           = 1.05
performance (x real-time) = 2.2357
```

**The model does NOT stay upright with default parameters — it falls at t = 1.05 s.**
From `balance_osim.sto` (parsed with `sto_summary.py`): COM height 0.9649 m → 0.5780 m =
**59.9% of initial**, which trips BalanceMeasure's `termination_height = 0.6` (relative COM
height at which the sim stops; per the SCONE reference wiki, default 0.5). Pelvis pitch
runs away −0.175 → +0.94 rad, COM drifts **backward** −0.057 → −0.814 m, and GRF collapses
to 0.08 BW by termination. This is the tutorial working as designed, not a defect: in
`ControllerReflexBalance.scone` the vestibular gains start at zero (`$KP = 0~0.1`,
`$KV = 0~0.1`), so only the fixed-gain muscle-length reflexes act; CMA-ES is supposed to
find the balance gains.

Second command (A/B for override syntax, log `scone_logs\balance_osim_md1.log`): with the
correct full dot-path the override applied cleanly (no "unused properties" warning), the
sim terminated at exactly `1.000` (max_duration, before the fall) and `BalanceMeasure = 0`
— i.e. **BalanceMeasure scores 0 when the model is still up when max_duration ends**, and
non-zero only for fall/drift.

### Gotcha banked: override dot-path needs the optimizer root

`SimulationObjective.max_duration=3` is **silently ignored** (just a
`Warning, unused properties:` line) — the scenario root is `CmaOptimizer`, so the working
form is:

```
CmaOptimizer.SimulationObjective.max_duration=<seconds>
```

Both forms were executed this session (first flagged unused, second applied — see logs).

## 2. Tutorial 4a — Gait — OpenSim (headless, short duration)

Scenario: `Tutorial 4a - Gait - OpenSim.scone` — same model (`ModelOpenSim3`, H0918v3,
`init/InitStateH0918Gait10ActA.zml`), controller include `controllers/H0918RS2v3.scone`,
measure = CompositeMeasure{ Gait10.scone (GaitMeasure `min_velocity = 1.0`,
`termination_height = 0.85`), EffortWangCubed2000.scone, DofKnee1.scone, Grf14.scone }.

```
sconecmd -e "...\Tutorial 4a - Gait - OpenSim.scone" CmaOptimizer.SimulationObjective.max_duration=4 -r "D:\temp\scone_hands_on_20260924\gait_osim" -l 2
```

Result (log `scone_logs\gait_osim.log`):

```
Terminating simulation at 2.780
result             = 51.7277
  Gait             = 49.5578 <- 100 * (0.495578 > 0.05)
    step_velocity  = 0.523106
    step_count     = 7
  Effort           = 0.74174 <- 0.1 * 7.4174      (effort 1382.62, distance 2.50099 m)
  MuscleActivation = 0.82625 <- 2000 * 0.000413125
  DofLimits        = 0.26133   (knee limit torques 7.42 / 5.65 exceed the 1.0 ceiling)
  GRF              = 0.340579 <- 10 * 0.0340579
simulation time           = 2.78
performance (x real-time) = 1.17147
```

**The model STEPS with default parameters: 7 steps covering 2.50 m in 2.78 s** (raw pelvis
speed ≈ 0.95 m/s; the measure's `step_velocity` = 0.523 m/s against the 1.0 m/s target),
then falls — termination at 2.78 s is GaitMeasure's `termination_height = 0.85` (parsed
`gait_osim.sto`: COM 0.9564 → 0.8124 m = 84.9% of initial). Final state at termination:
pelvis still at 0.979 m height, but both hips slammed into extension (−1.55 / −1.31 rad)
— a stumble, knees near extension (−0.13 / +0.05 rad); peak GRF 2.12 BW (right leg).
Again by design: these are unoptimized seed gains; the tutorial is a CMA-ES starting point
(`Tutorials3\par\H0918GaitRS2Hfd4.par` exists for warm-starting).

## 3. Tutorial 4a — Gait — Hyfydy: exact failure (diagnosis)

```
sconecmd -e "...\Tutorial 4a - Gait - Hyfydy.scone" CmaOptimizer.SimulationObjective.max_duration=3
```

Captured verbatim (log `scone_logs\gait_hyfydy.log`):

```
04:49:10 Evaluating C:\Users\Ben Bolen\Documents\SCONE\Tutorials3\Tutorial 4a - Gait - Hyfydy.scone
04:49:10 Error creating CmaOptimizer:
04:49:10   Error creating SimulationObjective:
04:49:10     This scenario uses a Hyfydy model, but no active license key was found.
04:49:10     Please check Tools -> Preferences -> Hyfydy, or visit https://hyfydy.com for more information.
```

Diagnosis: the scenario's `ModelHyfydy { model_file = models/H0918v3.hfd ... }` block
(scenario lines 9–20; the OpenSim twin uses the identical controller/measure includes)
fails at model creation because **no Hyfydy license is activated on this machine** — Ben
holds a trial key that has not been entered yet. Nothing else is wrong: the error is raised
before any simulation step, the rest of the scenario (controller, measures) never loads.

**One-minute activation for Ben** (Studio is currently open — do NOT close/kill it; do NOT
hand-edit `%LOCALAPPDATA%\scone\scone-settings.zml` while it runs):

1. In SCONE Studio: **Tools → Preferences → Hyfydy**.
2. Paste the provided license key into the key field → **Enable** (Studio writes the
   settings itself on OK).
3. CLI alternative (Studio may stay open): `"C:\Program Files\SCONE\bin\sconecmd.exe"
   --hyfydy <the provided license>`; `--hyfydy id` prints the hardware ID if the key
   request ever needs to be redone.

The Hyfydy scenarios are testable the moment he activates — `Tutorial 4a - Gait -
Hyfydy.scone` evaluates with exactly the command in section 3 (swap the model block for
the license; no other edits needed). Hyfydy is the 10–50× faster backend and the only one
with uneven-terrain support, so activating before any long optimization campaign is worth
the minute.

## 4. How SCONE expresses proprioceptive + vestibular control (controller reading)

Files read this session: `Tutorials3\controllers\ControllerReflexBalance.scone`,
`Tutorials3\controllers\ControllerGH2010v12.scone`, `Tutorials3\controllers\H0918RS2v3.scone`
(the one Tutorial 4a actually includes), plus the measure includes
`measures\MeasureBalance.scone` and `measures\Gait10.scone`.

### Proprioception: muscle sensors → delayed reflex arcs

- The atom is `MuscleReflex { target, source, delay, C0, KF, KL, KV }` inside a
  `ReflexController`. Stimulation law (SCONE reference / Geyer-Herr 2010 form):
  `U = C0 + KF·[F−F0]+ + KL·[L−L0]+ + KV·[V−V0]+` — F/L/V are normalized muscle force /
  contractile-element length / velocity; brackets are rectified unless `allow_neg_*`.
  `source` omitted = the muscle's own spindle/Golgi signal (autogenic mono reflex);
  `source = <muscle>` = cross-muscle coupling; negative gains implement reciprocal
  antagonist **inhibition** (e.g. `tib_ant source = soleus KF ~-1.0`,
  `iliopsoas source = hamstrings KL ~-3..-5` = late-swing hip-flexor switch-off).
- **Delays are explicit per arc** and tiered distally in `H0918RS2v3.scone`
  (`$hip_delay = 0.01, $ham_delay = 0.015, $knee_delay = 0.02, $ankle_delay = 0.035` s;
  `ControllerReflexBalance.scone`: 10 ms hip/proximal, 20 ms knee/biarticular, 35 ms
  ankle). GH2010v12 uses 5–20 ms with paper values quoted in comments.
- **Phase gating**: `GaitStateController` runs a per-leg 5-state FSM (`EarlyStance,
  LateStance, Liftoff, Swing, Landing`) driven by `leg_load` [body weights] vs
  `stance_load_threshold` (~0.2–0.3) and sagittal foot-position thresholds; each
  `ConditionalController { states = "..." ReflexController {...} }` switches a whole
  reflex set in/out per phase. `ConditionalMuscleReflex` adds a joint-angle window on top
  (vasti positive-force reflex active only while `knee_angle < −0.175 rad`, i.e. flexed).

### Vestibular / balance: two expression patterns

- **BodyPointReflex (standing balance, Tutorial 3a)**: `BodyPointReflex { target = <each
  of the 9 muscles> source = torso KP KV delay = 0.1 offset = [0 0.5 0] direction =
  [1 0 0] }` — PD feedback of a body point's position/velocity along a world direction,
  broadcast to every muscle with a **100 ms** (vestibular-timescale) delay; gains start at
  0 and are found by optimization.
- **DofReflex on pelvis_tilt (gait balance, both gait controllers)**: PD on a DOF signal —
  `DofReflex { target = hamstrings/glut_max source = pelvis_tilt KP ~1.91 KV ~0.25
  P0 = −0.105 C0 ~±0.05..0.1 delay = 0.005 }` plus mirrored negative-gain arcs to
  iliopsoas/rect_fem (push the trunk back / catch it forward). In H0918RS2v3 the same
  reflex exists **separately in stance and swing** with different gain sets and a 50 ms
  delay — i.e. balance feedback is phase-gated, not a fixed controller.
- Termination criteria live in the measures, not the controller: COM below
  `termination_height` (0.6 standing / 0.85 gait, relative to initial COM) ends the sim
  and scores the fall; GaitMeasure additionally requires `min_velocity` (1.0 m/s) and
  exempts the first `initiation_steps`.

### Mapping to our goal-2 balance stage (spinal/VEST work)

- SCONE's entire balance controller = trunk-pitch PD **distributed onto muscle
  stimulations with explicit delays** — exactly the role our VEST vestibular-analog cells
  (added 2026-09-23, `runner --vest/--vest-flexinh/--vest-prop`) play. Useful priors to
  steal: proprioceptive delays 10/15/20/35 ms per compartment (hip/hamstring/knee/ankle),
  vestibular-pathway delays 50 ms (phase-gated gait) to 100 ms (standing), and per-muscle
  tonic offsets (`C0`, `P0`) which correspond to our per-MN bias injection.
- Their phase gating (stance vs swing balance gain sets) argues for gating our VEST→MN
  gains by RG phase rather than one static gain; their antagonist-inhibition reflexes
  (negative cross-muscle gains) map onto our IaIN reciprocal-inhibition pathway.
- The `termination_height` convention (0.6 standing / 0.85 gait relative COM) is a
  ready-made acceptance metric for our `--stand-eval` / gait evals.

## 5. Operational notes (all verified this session)

- Evaluate (single sim, default parameters): `sconecmd -e <scenario.scone>` — works
  directly on a .scone; a `.par` instead needs `config.scone` beside it.
- Override dot-path: `CmaOptimizer.SimulationObjective.max_duration=<s>` (shorter paths
  are ignored with an "unused properties" warning — see §1).
- `-r out` writes `out.sto` (no extension in the argument); `-l 2` prints the full
  objective breakdown (step_count, step_velocity, distance, per-term weights);
  `-q` silences everything including the result line.
- Speed on this box, single-sim evaluation: standing balance 2.24× real-time, 4 s gait
  1.17× real-time (OpenSim3 backend; the pooled evaluator spawns 12 threads even for `-e`).
- Models: H0918v3 = 9 DOF / 18 muscles / 74.53 kg; both OpenSim backends initialize
  automatically, no external OpenSim install involved.
- Full `sconecmd --help` captured to `scone_logs\sconecmd_help.txt`.

## Artifacts

- `scone_logs\balance_osim.log`, `balance_osim_md1.log`, `gait_osim.log`,
  `gait_hyfydy.log`, `sconecmd_help.txt` (this folder).
- `D:\temp\scone_hands_on_20260924\`: `balance_osim.sto`, `gait_osim.sto` (motion
  replays, openable in Studio), `sto_summary.py` (column summarizer), same logs.
- Scenarios/controllers cited from `C:\Users\Ben Bolen\Documents\SCONE\Tutorials3\`.


## ADDENDUM (2026-09-24 morning): Hyfydy ACTIVATED by Ben - engine A/B done

Ben entered the trial license (settings now `hyfydy { enabled = 1 }`,
Hyfydy 1.12.6.1412 initializes). Both previously-blocked scenarios now
run headless via sconecmd:

| tutorial | OpenSim engine (overnight) | Hyfydy engine (now) |
|---|---|---|
| 4a Gait (max_duration=4) | 2.780 s sim, 7 steps, 2.50099 m, effort 1382.62, DofLimits 0.26133, 1.17x real-time | 2.33 s sim, 7 steps, 2.58029 m, effort 1212.87, DofLimits 0.03284, **90.5x real-time** |
| 3a Balance (max_duration=3) | terminated 1.050 s, BalanceMeasure 96.5/100, fell (com_y -> 59.9%) | 1.045 s sim, BalanceMeasure 65.2/100, fell, **99.1x real-time** |

Caveat: the two engines' scenarios carry DIFFERENT .par/.zml
initializations (Hfd4 vs Osim4 files), so treat this as a
both-engines-work comparison, not a same-controller A/B.

**Goal-3 headline**: the Hyfydy log literally shows the three features
MuJoCo 2.3.7 lacks: `muscle_force_m2012fast` (Millard-style elastic
tendon muscle), `contact_force_hunt_crossley_sb` (exact HC contact),
`predictive_integrator_psemtw` (adaptive/error-controlled integrator)
- plus HFD model format. The post-deadline engine port option from
goal3_mujoco_fidelity.md 3.2 is now CONCRETE and licensed until
2026-10-24. Logs: scone_logs/{gait_hfd,balance_hfd}.log + .sto results.
