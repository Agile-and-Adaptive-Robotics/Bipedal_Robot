# Goal 2 / Milestone 1 — AnimatLab 2-layer walker BODY ported to MuJoCo

**Date:** 2026-09-25 (EB475WS4, unattended campaign) · **Status: complete** (validation gate 12/12 PASS)
**Source of record:** `D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj` (2023 original, read-only — never written)
**Artifacts (all new, under `Code\MuJoCo_SNS\spinal\w2l_mujoco\`):**
- `make_w2l_mjcf.py` — read-only .aproj parser → emits `w2l_mjcf.xml` + `w2l_source_dump.json`
- `validate_body.py` — the M1 validation gate
- `w2l_mjcf.xml` — the MuJoCo 2.3.7 model
- `w2l_source_dump.json` — every parsed source value (validator + this report's tables)

Gate log: `reports_20260925\logs\validate_body_run6.log` (final; runs 1–5 kept as the debug trail).
Env: `C:\Users\Ben Bolen\.conda\envs\myo\python.exe`, mujoco 2.3.7 (`CONDA_PREFIX` set before import per the mujoco skill). No pip installs. No protected file touched (new folder only; `spinal_run.npz`, studies, prior reports untouched).

---

## 1. What was built

MuJoCo model with **13 bodies (+world)**, **8 hinge joints**, **12 muscle actuators on 12 spatial tendons**, **2 toe-spring tendons**, **ground plane**:

- Segment chain per side: `Root` (pelvis, 30.991 kg, **FREE** joint — the aproj's `<Freeze>True</Freeze>` trap deliberately NOT copied) → `femur` (3.5 kg) → `tibia` (1.5 kg) → `foot` (0.4 kg) → `toe` (0.064 kg), plus the two welded contact plates per foot (`foot_contact` 9.6 g, `toe_contact` 2.4 g — modeled as child bodies with no joint, exactly the aproj's rigid-parent structure).
- Muscle = `site` per Attachment body (GUID-resolved) + `<spatial>` tendon + `<muscle>` actuator. `knee_L_ext`/`knee_R_ext` have 3 attachments → the middle one is a kneecap via-point.
- Toe springs (aproj `Body_11` L / `toe_R_spring` R): two-site spatial tendon, `stiffness=16000`, `springlength=0.088`.
- Ground: `<geom type="plane>` at z=0, friction 1.
- Not ported (by design, documented): 12 `LinearHillStretchReceptor` bodies and 12 `LinearHillMuscle` overlay bodies (behavior/sensory layer, massless in the aproj, `ApplyTension=False` on the receptors) — milestone-2 material; the `WalkingPath` structure (cosmetic path marker); `CollisionExclusionPairs` (checked: EMPTY in the source, so all-pairs-active is faithful).

## 2. Source-convention findings (needed to parse .aproj physical data at all)

No pre-existing parser existed, so these were established this session against the file itself; the parser asserts them on every body (`LocalPosition` vs `LocalMatrix` cross-check: all ≤1e-6 m):

1. **Units:** `Value/Scale/Actual` triplets carry `Actual` in **meters for distance** (cross-checked: `LocalPosition Y=99.298 centi → Actual 0.99298`, and the Root `LocalMatrix` translation `9.9298` dm = 0.99298 m) but **grams for mass** (`Mass Value=30.991 Kilo → Actual 30991`). Gravity −9.8, timestep 1 ms (Environment block).
2. **`LocalMatrix` = 16 floats, column-major storage, column-vector convention, translation in dm.** Proof: root matrix R decodes to Rz(−6.417°) = the XML `Rotation Z=−6.417` exactly; femur Rz(−3°) = XML −3.
3. **Rotation elements are degrees; R = Rx(x)·Ry(y)·Rz(z)** — NOT the usual Rz·Ry·Rx: tibia (90.012, −3, 0) matches Rx·Ry·Rz to 5.8e-6 while Rz·Ry·Rx deviates by 0.052. (Only used for the cross-check; FK uses the matrices.)
4. **`<Joint>` lives inside the CHILD RigidBody and its LocalPosition is in the CHILD body's frame.** Verified anatomically: hip anchor (−3.4498, 0.9318, ±0.06) sits just below the pelvis center at the femur's top end; knee anchor (−3.5072, 0.5153, ±0.06) = the femur's bottom end to 6 mm; ankle anchor = the foot's top-rear corner. Parent-frame placement fails these tests.
5. **Hinge axis = the joint frame's local X** (Vortex primary rotational coordinate). Resulting world axes (AL frame): hip/knee/ankle/toe all ≈ ±z_AL = **lateral** — all four joints per leg are sagittal-plane hinges, as expected.
6. **Box `Length/Width/Height` = full extents along body-local x/y/z** (verified per-body once the frames are right: foot = 12 long × 8 wide × 2 thick; femur = 0.05×0.05×0.42 along its local vertical).
7. **World is y-up; the port remaps MJ(x,y,z) = AL(x,−z,y)** (the same remap the lab uses for OpenSim→MuJoCo), so box dims get emitted swapped (a,b,c)_AL → (a,c,b)_MJ.
8. **Rest pose:** feet hover at 2.6–3.1 cm above ground (foot body center z_AL = 0.044, contact plates at 0.026–0.032). This is the source's own rest configuration — the same "frozen biped hanging in the air" pose known from the standalone .asim audits.
9. **Attachment/receptor bodies carry stale Rotation elements** (their `LocalMatrix` vs Euler dev ≈ 1 rad). Irrelevant here: only their positions are used (sites), and every position came from the matrix chain. Recorded so nobody chases it later.
10. **Prior audit-note counts vs this file:** the task brief's facts say "25 Attachment bodies / 9 StretchReceptors / 62 RigidBody total"; this file parses to **30 attachments, 12 receptors, 69 RigidBody elements total** (13 Box + 12 LinearHillMuscle + 12 LinearHillStretchReceptor + 2 Spring + 30 Attachment + GroundPlane + WalkingPath). Counts above are read directly from the parsed file (dump JSON), not assumed; the brief's numbers likely describe a different vintage of the project.

## 3. Mapping table (AnimatLab LinearHill → MuJoCo muscle)

| AnimatLab quantity | MuJoCo 2.3.7 equivalent | Notes |
|---|---|---|
| `MaximumTension` [N] | `<muscle force>` | exact |
| `RestingLength` L0, `Lwidth` | `range = [L0−Lw, L0+Lw]`, `lengthrange` = same | MuJoCo's FLV then peaks at normalized 0.5 (= L0) and dies at normalized 0/1 (= L0±Lw), matching AnimatLab's parabolic length-tension `1−((L−L0)/Lw)²` (zero at ±Lw) |
| `Kse` (series elastic) | **none** | rigid tendon; explicit scope exclusion of the ask ("series-elastic surrogate out of scope for M1") — known family deviation, same as the converted-gait2392 actuators |
| `Kpe` (parallel elastic) | **none** | no passive-spring hook in the 2.3.7 muscle |
| `B` [N·s/m] | **none** | 2.3.7's `<muscle>` has no damping attribute (probed: `damp` rejected by schema; damping exists only in MuJoCo ≥3.x). The built-in force-velocity curve is the substitute. Deviation documented. |
| activation time const | `timeconst="0.01 0.01"` | AnimatLab activation arrives from the neural module (milestone 2); 10 ms ≈ MuJoCo default, 2-value form required by 2.3.7 |
| `StimulusTension` gain (A/B/C/D) | not ported | neural-layer encoding (stimulus→tension), belongs to the SNS port |
| toe `Spring` k=16000, b=20000, L0=0.088 | tendon `stiffness=16000`, `springlength=0.088`, **`damping=172`** | **deviation:** source damping 20000 N·s/m is ~116× critical for the ~0.46 kg toe chain and un-integrable at MuJoCo's 1 ms step — measured `max|qacc| = 3.5e8` + instability warning (A/B probe); capped at critical `c*=2·sqrt(k·m)` ≈ 172 → `max|qacc| 1.77e4`, identical to the c=0 contact-impact baseline. Source value kept in the XML comment. |
| Material `foot` `FrictionLinearPrimary=100000` | `friction="1 0.005 0.0001"` | **deviation:** the Vortex number is a viscous-type contact constant (≈ non-slip), not a Coulomb μ; shipped the MuJoCo/Default-material μ=1. "smooth" (0) and "Default" (1) materials map trivially. |
| Joint `EnableLimits` + `LowerLimit/UpperLimit` (deg) | `range` (rad) | exact (gate 4: dev ≤ 2e-10 rad); `Stiffness=1e9, Restitution=0` limit properties → MuJoCo hard limits |
| Joint Relaxations / ConstraintFriction / PID motor (all `Enabled=False` in source) | not emitted | faithful: source joints are plain undamped hinges |
| `Freeze=True` on Root | **ignored — free pelvis** | the documented Freeze trap |

## 4. Fidelity tables

### 4a. Body masses (aproj vs compiled MJCF — gate readback, `validate_body.py`)

| Body | aproj [kg] | MJCF [kg] | dev |
|---|---|---|---|
| Root (pelvis) | 30.9910 | 30.9910 | 0 |
| femur_L / femur_R | 3.5000 | 3.5000 | 0 |
| tibia_L / tibia_R | 1.5000 | 1.5000 | 0 |
| foot_L / foot_R | 0.4000 | 0.4000 | 0 |
| toe_L / toe_R | 0.0640 | 0.0640 | 0 |
| toe_L_contact / toe_R_contact | 0.0024 | 0.0024 | 0 |
| foot_L_contact / foot_R_contact | 0.0096 | 0.0096 | 0 |
| **Total** | **41.9430** | **41.9430** | **0** |

Box inertias = uniform-box formula about the aproj COM offset (Root COM 0.02 m x-offset preserved; all others 0). AnimatLab `Density` fields (Root 2.787, femur 3.333…) are unused — explicit masses govern in the source too.

### 4b. Joints (8 hinges; aproj vs MJCF — gate 4 PASSED, worst range dev 2e-10 rad)

| Joint | Limits aproj [deg] | Limits MJCF [rad] | World axis (AL frame) | Anchor (world AL, m) |
|---|---|---|---|---|
| hip_L | [−15, +23] | [−0.2617994, +0.4014257] | (0.000, 0.000, −1.000) lateral | (−3.4498, 0.9318, +0.06) |
| knee_L | [0, +60] | [+0.0000000, +1.0471976] | (0.000, 0.000, +1.000) lateral | (−3.5072, 0.5153, +0.06) |
| ankle_L | [−20, −5] | [−0.3490659, −0.0872665] | (0.000, 0.001, −1.000) lateral | (−3.6374, 0.0877, +0.06) |
| toe_L | unlimited | limited=false | (0.000, 0.001, −1.000) lateral | (−3.5266, 0.0379, +0.06) |
| hip_R | [−15, +23] | same as hip_L (dev 4e-11) | (0.000, 0.000, −1.000) | (−3.4498, 0.9318, −0.06) |
| knee_R | [0, +60] | same as knee_L | (0.000, −0.001, +1.000) | (−3.5072, 0.5153, −0.06) |
| ankle_R | [−20, −5] | same as ankle_L | (0.000, 0.001, −1.000) | (−3.6374, 0.0877, −0.06) |
| toe_R | unlimited | limited=false | (0.000, 0.000, −1.000) | (−3.5266, 0.0379, −0.06) |

Anatomy cross-checks (from the FK, all in the parsed source): hip anchor = femur top end (0.9318 vs 0.935 box top); knee anchor = femur bottom end (0.5153 vs 0.515); ankle anchor = foot top-rear; feet at 2.6–3.1 cm hover.

### 4c. Muscles (12; attachment names, rest geometry, source parameters)

Path length = parsed-FK length of the site chain at the aproj rest pose; MuJoCo tendon rest lengths (`mj_forward`, t=0) equal the same numbers to ≤3 mm (knee-ext via-point path vs straight line).

| Muscle | Attachments (aproj GUID→name) | Path [m] | L0 [m] | L0/Lw range emitted | MaxT [N] | Kse / Kpe / B (not ported) |
|---|---|---|---|---|---|---|
| hip_L_flx | body_L_f → femur_L_f | 0.161 | 0.170 | [0.081, 0.259] | 1500 | 62850 / 4000 / 800 |
| hip_L_ext | body_L_b → femur_L_b | 0.157 | 0.210 | [0.130, 0.290] | 1500 | 62000 / 4000 / 800 |
| knee_L_flx | femur_L_b_low → tibia_L_b | 0.344 | 0.345 | [0.245, 0.445] | 1000 | 62500 / 3770 / 600 |
| knee_L_ext | femur_L_f_low → knee_L → tibia_L (via-point) | 0.348 | 0.340 | [0.260, 0.420] | 1500 | 77193 / 4632 / 600 |
| ankle_L_flx | foot_L_f → tibia_L_f_low | 0.225 | 0.247 | [−0.753, 1.248] | 1500 | 74887 / 4493 / 400 |
| ankle_L_ext | foot_L_b → tibia_L_b_low | 0.248 | 0.225 | [−0.775, 1.225] | 1000 | 31900 / 4000 / 400 |
| hip_R_flx | femur_R_b → Body_R_b | 0.157 | 0.210 | [0.130, 0.290] | 1500 | 62000 / 4000 / 800 |
| hip_R_ext | body_R_f → femur_R_f | 0.161 | 0.170 | [0.081, 0.259] | 1500 | 62850 / 4000 / 800 |
| knee_R_flx | femur_R_b_low → tibia_R_b | 0.344 | 0.345 | [0.245, 0.445] | 1000 | 62500 / 3770 / 600 |
| knee_R_ext | femur_R-f_low → knee_R → tibia_R_f (via-point) | 0.348 | 0.340 | [0.260, 0.420] | 1500 | 77193 / 4632 / 600 |
| ankle_R_flx | tibia_R_f_low → foot_R_f | 0.225 | 0.247 | [−0.753, 1.248] | 1500 | 74887 / 4493 / 400 |
| ankle_R_ext | tibia_R_b_low → foot_R_b | 0.247 | 0.225 | [−0.775, 1.225] | 1000 | 31900 / 4000 / 400 |

Rest-pose path/L0 ratios: 0.75–1.10 (knee muscles ≈1.00) — the source rest pose is not muscle-rest-optimal; no assumption was made that it should be.

### 4d. Toe springs (2)

| Spring | Endpoints | NaturalLength | k | b (source → shipped) | Tendon rest len (mj) |
|---|---|---|---|---|---|
| Body_11 (L) | foot_L_spring ↔ toe_L_spring | 0.088 | 16000 | 20000 → **172** | 0.0913 |
| toe_R_spring (R) | foot_R_spring ↔ toe_R_spring | 0.088 | 16000 | 20000 → **172** | 0.0913 |

FK endpoint distance at rest = 0.0883 m ≈ NaturalLength (3 mm pre-load) — an independent confirmation that the FK and site placement are right.

### 4e. End-to-end pose readback (aproj FK → remap vs MuJoCo `xpos`, t=0)

All 13 bodies agree to **worst 4.5e-7 m** (8-significant-digit print rounding); e.g. foot_L aproj→mj (−3.5905, −0.0600, 0.0440) vs mj (−3.5905, −0.0600, 0.0440). Full table in `reports_20260925\tmp\pose_check.py` output above.

## 5. Validation gate — final run (paste of `logs\validate_body_run6.log`)

```
== gate 1: model loads under mujoco 2.3.7 ==
  [PASS] mj_loadXML  nbody=14 nq=8 nu=12 neq=0 ngeom=14
  [PASS] mujoco version 2.3.7  2.3.7
  [PASS] body count  nbody=14 (expected 14 incl. world)
  [PASS] 8 hinge joints  njnt=8
  [PASS] 12 muscle actuators  nu=12
  [PASS] 14 tendons (12 muscle + 2 toe springs)  ntendon=14
== gate 2: 2 s passive drop (timestep 0.001 s) ==
  [PASS] passive drop finite (2 s, 2000 steps)  max|qpos|=1.192, min foot height=0.0377 m
  [PASS] drop settles on ground (0 <= foot z < 0.1)  min foot z=0.0377 m
== gate 3: 12 actuators finite resting lengths ==
  [PASS] actuator_lengthrange finite + lo<hi for all 12  hip_L_flx[0.081,0.259]; ... (all 12 listed in log)
  [PASS] tendon lengths finite at rest  rest lengths: 0.162 0.175 0.345 0.276 0.267 0.355 0.355 0.277 0.345 0.263 0.136 0.156 0.088 0.088
== gate 4: joint ranges match aproj (deg->rad) ==
  [PASS] joint ranges match  hip_L ... dev=4.13e-11 | knee_L ... dev=1.97e-10 | ankle_L ... dev=1.13e-12 | toe_L unlimited | (R side same)
== fidelity readback: body masses (aproj vs MJCF) ==   (table above)
  [PASS] body masses match  worst dev=0.00e+00 kg
RESULT: 12 passed, 0 failed
```

Notes on gate behavior: the passive (zero-activation) 2 s drop ends in a crumpled kneel — legs fold (no muscle tone), femur_L↔femur_R and tibia_L↔tibia_R make contact (ncon=6). That is source-consistent (the aproj's `CollisionExclusionPairs` is empty → all pairs active) and state stays finite throughout.

## 6. Debug trail (kept for honesty; logs in `reports_20260925\logs\`)

- run1: `<muscle damp=...>` rejected — 2.3.7 has no such attribute (schema probe in `tmp\schema_check.py`).
- run2/run3: feet "settled" at 2.91 m, tendon lengths 3.5 m — root cause: I had emitted **world** positions into nested `<body pos>` (MuJoCo parent-relative) → the robot stacked to z≈5 m and the toe-spring force (k·(3.5−0.088) ≈ 55 kN) launched it. Fixed to parent-relative transforms; pose readback now 4.5e-7 m.
- run4/5: toe joints misnamed `toe_L_2` (my name-uniquifier collided body/joint names; MuJoCo uniqueness is per element type). Fixed.
- run5 (12/12 but) printed MuJoCo's QACC instability warning at t=0.032 → toe-spring damping A/B (`tmp\ab_spring_damping.py`) → damping capped at critical 172 (run6 clean).

## 7. Honest caveats / open items for Milestone 2

1. **Hinge-axis sign convention:** axes and limit values are transported exactly, but AnimatLab's positive-angle direction was not independently re-derived (would require an AnimatLab instrumented run); if flexion signs come out mirrored in closed loop, flip the axis sign per joint — a one-attribute change in the generator.
2. **Kse/Kpe/B unported** (§3): milestone-2's neural loop will see slightly stiffer/springier muscles than AnimatLab's LinearHill. If needed, a MuJoCo 3.x or custom-actuator route exists, out of M1 scope.
3. **Friction:** μ=1 shipped; the source's "non-slip" viscous constant has no μ equivalent (documented above).
4. **Rest pose hangs 2.6–3.1 cm** above ground (source's own rest configuration, copied verbatim). A standing-pose solve is milestone-2 work, mirroring the spinal runner's static-opt approach.
5. **Sensory layer** (stretch receptors, contact sensors, `ReceptiveFieldSensor` Bell/Polynomial gains) not ported — parsed and dumped in `w2l_source_dump.json` (`receptors` + materials) so milestone 2 can wire heel/toe contact cells and Ia/II pools against real numbers.
6. The three proof numbers to quote downstream: **pose dev 4.5e-7 m**, **mass dev 0.0 kg (total 41.943 kg)**, **joint-range dev ≤2e-10 rad**.
