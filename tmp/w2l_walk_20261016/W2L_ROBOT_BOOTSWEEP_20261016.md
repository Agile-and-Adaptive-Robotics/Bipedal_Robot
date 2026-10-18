# W2L ROBOT BOOTSTRAP DRIVE SWEEP — stage report (2026-10-16, egg)

Overnight-campaign stage: **"bootstrap driver (robot)"** (phase 7 of the
11-phase plan in `CAMPAIGN_HANDOFF_20261016.md`). All work ran ON
easteregg2 (repo `D:\GitHub\Bipedal_Robot`, branch KneeTestSetup_BenBo_stw)
via ssh; python `D:\Anaconda\envs\myo\python.exe`. NO commits, NO pushes.
Sweep outputs: `reports_20261016\w2l_walk_sweep\` (this report's data:
`results.csv` [128 runs], `logs\*.log`, probe npz's, driver scripts).

The ask: *sweep strong asymmetric drive under the harness/virtual-walker
protocol — does the ROBOT take alternating supported steps (>= 6 stride
cycles, antiphase L/R, forward progression > 0)? Report the best combo +
param file, or exactly why not.*

## 1. VERDICT

**YES — with the harness made faithful to Ben's protocol, the W2L robot
takes alternating supported steps with sustained forward progression.**
Winner (param file `w2l_walk_sweep\w2l_walk_best_params.json`):

| knob | value | | metric (25.5 s window) | value |
|---|---|---|---|---|
| te / tf | 4.0 / 5.0 nA | | primary stride cycles L/R | **22 / 21** |
| alpha_l/r, te_l/r, tf_l/r | 0 / unset | | period L / R | 1.199 / 1.201 s |
| legcoll | **0** (excl. L-R leg collision) | | inter-leg phase median | **0.498** (n=20) |
| sink | 0.06 m (lowering depth) | | RG-E L/R corr | −0.746 |
| ctn_amp | 2.0 nA (contact encoder) | | forward displacement | **+1.676 m** |
| rig | 0.25 (leash+tilt at 25%) | | sustained speed (fit, 8–30 s) | **0.063 m/s** |
| walker_yield | **1.5 s** (moving anchor) | | heel duty L/R | 0.060 / 0.069 |
| dur_walk | 30 s | | combined foot-loaded duty L/R | 0.128 / 0.171 |

Stays up the whole run: no fall (pelvis z >= 0.85 m vs 0.60 limit), tilt
max 16.1 deg (vs 45 limit), finite throughout, harness carries ~30% weight
by late run. Reproduce (one command, cwd `w2l_mujoco\`):

```
D:\Anaconda\envs\myo\python.exe test_w2l_ground.py --phase=walk --te=4 --tf=5 --legcoll=0 --sink=0.06 --ctn_amp=2 --rig=0.25 --walker-yield=1.5 --dur-walk=30
```

**Honest caveats** (this is a bootstrap, not a finished gait):
- The test's own strict verdict stays SHUFFLES-IN-PLACE: its WALKS gate
  additionally requires heel duty >= 0.3 (we have 0.05-0.07; even combined
  heel+toe loading is 0.13-0.17 vs the campaign's 0.4-0.75 target). The
  feet are loaded only a fraction of the cycle; much of the body weight
  rides the harness.
- Heels "double-tap" (strike + 0.3 s bounce); primary-cycle counting
  merges taps (0.45 s); cycles above are post-merge.
- The walk leans (tilt ~16 deg at rig 0.25). The more-upright runner-up
  (rig 0.5, walker_yield 3.0: tilt 7.3 deg, 24/37 strikes, +0.437 m,
  r −0.770) is in the param file too.
- Progression requires the yielding anchor (see §3.3): with the FIXED
  anchor the same combo shifts +0.10 m during the drop and then marches
  in place; with the leash fully off (rig 0) it topples (fwd +2.8…+17 m,
  tilt 97–180 deg — fall-slide, not walking).

## 2. HOW THE SWEEP WENT (128 runs, ~6.6 s per 20 s sim)

| stage | what | runs | outcome |
|---|---|---|---|
| A | symmetric te 2–6 × tf 3–7 | 25 | **zero heel strikes anywhere**; neural rhythm scales 4→21 bursts, antiphase holds (r −0.42…−0.88) |
| B | E-dominant + asymmetric gamma (±20 %, ±100 % splits) | 20 | best R-only strikes (a11/b14: 4-6); L heel NEVER loads |
| C | alpha_l/r 0.05–0.15 sym/one-sided | 12 | post-fix: R-only strikes 7-12, duty <= 0.05; no L |
| D | stance/harness probes (acap, wsup, sink, ctn_amp, cref) | 14 | nothing loads the L heel — prompted the geometry probe |
| E | legcoll=0 + sink sweep | 13 | **first bilateral strikes** (11/11 at sink 0.04-0.06) |
| F | refine (ctn_amp 2, drive, wsup, acap) | 18 | bilateral strikes everywhere; best 22/16 |
| G | rig/progression | 10 | rig 0 topples; rig 0.25 stays up, fwd +0.12 m |
| H | rig refinement (30 s runs) | 8 | h07 (te4/tf5, rig 0.25): 40/41 strikes, tilt 13.7 |
| W | yielding walker | 8 | **w01: +1.676 m sustained** |

Raw per-run rows: `w2l_walk_sweep\results.csv`; per-run logs in
`w2l_walk_sweep\logs\`.

## 3. THREE BLOCKERS FOUND (the "exactly why not" of the early stages)

### 3.1 [BUG — FIXED] alpha drive shadowed the physics loop index
`test_w2l_ground.py` line ~296 (audit-stage code): `for i, a_cur in
i_alpha:` INSIDE `for i in range(nphy):` rebound `i`. With alpha ports
built, every even-index log entry was stomped, the fall check read stale
tilt, and the rig anchor capture (`if i % NSUB == 0 ...`) never fired —
**the leash/tilt never engaged and the body wandered to x = −11.35 m,
y = −1.0 m** (measured in `probe_c05.npz`). The audit's "ALPHA 0.3 nA
regime change (harness 0 %, shuffle)" interpretation was read off this
corrupted logging. Fix: rename to `i_a` in BOTH
`test_w2l_ground.py` and `test_w2l_air_afferented.py` (line ~166).
alpha=0 → loop body never runs → bit-identical. Verified post-fix:
- default walk gate reproduces the recorded baseline EXACTLY (bursts
  11/10, period 1.539, r −0.776, tilt 8.2, harness 49 %, fwd +0.001);
- M5 canonical air gate = recorded state (a/b PASS, c FAIL-as-specified,
  suppression nan — identical signature to the audit's finding 5).

### 3.2 [BODY DEFECT — runtime knob] the two tibias interpenetrate at spawn
`_w2l_static_pose.py` (kept in this folder): the aproj is a PLANAR walker
— both legs share the sagittal plane (hips at y = ±0.06, near-straight
knees) — so in the 3D MJCF the tibia capsules overlap at the rest pose
and MuJoCo wedges them apart with **3368 N** (5 contact points, tibia_L /
tibia_R). The wedge splays the feet asymmetrically: the L heel plate
NEVER reached the ground in ANY of the first 71 runs (heelF L max = 0.0 N
measured), while the R foot did all the work. The source aproj has NO
leg-leg contact, so excluding L<->R pairs is a FIDELITY FIX toward the
source model. New default-off knob `--legcoll=0` (runtime
contype/conaffinity remask: L=1/1, R=2/2, Root+ground=3/3 — the ONLY
behavioral change is L<->R off). Default 1 = shipped behavior =
bit-identical (re-verified twice). With legcoll=0 the feet behave
symmetrically and, given deeper lowering (sink 0.04–0.06), BOTH heels
strike.

### 3.3 [PROTOCOL GAP — runtime knob] the fixed anchor makes progression impossible
The AnimatLab virtual walker is the MOVING WalkingPath body; the test
leashed the COM to a FIXED settled anchor (2000 N/m at rig 1). Winner
combo h07 with the fixed anchor: +0.105 m during the drop transient,
then FLAT (−3.33 ± 0.005 m for 22 s) — a stationary march. New
default-off knob `--walker-yield=TAU`: the anchor's x follows an EMA of
the COM (tau seconds) — the leash still damps fast deviations (trolley
inertia) but yields to slow step thrust. Default 0 = fixed anchor =
bit-identical (re-verified). With tau 1.5 s the same drive walks
+1.676 m at a steady 0.063 m/s (linear fit over 8–30 s).

Also relevant: the spawn pose is toe-biased (toe plates ~5 mm below heel
plates at qpos0; ankle slams from 0 deg into its PF range in ~20 ms) —
the sink 0.06 lowering is what finally seats the heels.

## 4. What the DRIVE sweep itself said (ask question, drive-only answer)

With the shipped body (legcoll=1) and fixed anchor, NO drive combination
walks: 57 drive combos (stages A–C) produced ZERO left-heel strikes; the
L heel plate carried exactly 0.0 N for entire runs. Strong/asymmetric
gamma/alpha changes cadence and stance tone (bursts 4→21, antiphase
r −0.42…−0.88, tilt 2.8–11 deg) but cannot put the left heel on the
ground. Alpha (>0) helps the R foot strike earlier (7–12 strikes at
0.05–0.15 nA, te4/tf3) but never the L. The winning combo ended up
SYMMETRIC drive (te 4 / tf 5); the asymmetry knobs were not needed once
the body and harness blockers were fixed — worth knowing before tuning
per-side drive further.

## 5. Files changed (working tree only; see CHANGES_MANIFEST for hunks)

- `w2l_mujoco\test_w2l_ground.py` — i-shadow fix; `legcoll` +
  `walker_yield` default-off knobs (both default bit-identical, gates
  re-run).
- `w2l_mujoco\test_w2l_air_afferented.py` — same i-shadow fix.
- Everything else this stage produced lives in
    `reports_20261016\w2l_walk_sweep\` (driver, probes, combos, CSV,
    logs, npz's, best-params json) — nothing else touched.

## 6. NOT RUN (honesty list)

- No optuna/curriculum tuning of the new knobs (te/tf were swept on a
  grid; tau/rig/sink/ctn_amp coarsely). The knobs are TEST-file KV knobs,
  not params.G — any study that tunes them must follow the JSON RULE.
- No AnimatLab-side realization (phase 9's job).
- The audit stage's M5/M6 recorded baselines were re-verified as stated
  in §3.1 but the M6 gate numbers with legcoll=0/walker_yield>0 are NEW
  configurations (not comparable to the recorded baselines).
- No forward-kinematic reference exists for the W2L robot, so no kine
  score is reported (kine_ref is the human-stack acceptance metric).

— bootstrap driver (robot), 2026-10-16
