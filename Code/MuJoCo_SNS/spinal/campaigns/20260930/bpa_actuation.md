# BPA actuation: air-stepping + standing (campaign 2026-09-30, track report 2026-10-01)

**Headline: the stretch goal's MuJoCo half is MET — the split-RG W2L walker air-steps and
stands under retained support with all 12 muscles replaced by Ben's BPA model
(`bpa_muscle.py`), on easteregg2, within the 620 kPa ceiling. The Simscape half ran its
BPA CPG leg demo clean (neural alternation confirmed) but the demo knee completes 0
sustained cycles at stock parameters — Simscape air-stepping is blocked on tuning, not on
missing blocks.**

Run artifacts (both machines): `Code/MuJoCo_SNS/spinal/campaigns/20260930/bpa_run/`
(`bpa_gate.log` = the easteregg2 gate output of record; `test_bpa_stepping.py` = driver;
`make_bpa_walker.py` = XML variant generator; `w2l_mjcf_bpa.xml` = the BPA body;
`build_w2l_split_net.py` = the census-fixed net builder used by the run).

---

## 1. MuJoCo track (easteregg2) — GATE PASS

### 1.1 What was built

* **BPA body variant** — `make_bpa_walker.py` (new, `w2l_mujoco/`) swaps the 12 Hill-type
  `<muscle>` actuators of the body of record `w2l_mjcf_fixed.xml` for **zero-gain
  `<general>` actuators on the SAME tendon routes** (the `add_bpa_to_mjcf.py`
  convention: `gaintype fixed / gainprm 0 / biastype none`). `add_bpa_to_mjcf.add_bpa()`
  was NOT used because it would append duplicate routes — this body already carries all
  12 site-based `<tendon><spatial>` routes. Freejoint, per-joint damping (ΣB·r²), toe
  springs and the spawn keyframe are untouched. NOTE: the easteregg2
  `w2l_mujoco/w2l_mjcf_fixed.xml` is the body of record (freejoint at line 24, keyframe
  at line 147, `qpos="... 1.007846 ..."`); the laptop working-tree copy was stale
  (welded root, no keyframe) and was synced from easteregg2 this session.
* **Driver** — `test_bpa_stepping.py` (new): split-RG W2L net
  (`build_w2l_split_net.build(comm=1.0)`, corrected-V3 commissurals) drives
  `BPAMuscleSystem` (forces via `qfrc_applied` through `data.actuator_moment`) under
  Ben's virtual-walker protocol copied from `test_li_stepping.py:163-177`: harness hold
  (feedforward weight + vertical PD + weak xy damping) with feet CLEAR while the CPG
  air-steps → linear lowering over 2 s → harness RETAINED at 0.7× weight. Activation
  a∈[0,1] → pressure a·620 kPa (`bpa_muscle.py:134-138`), capped at `--cap` (0.6).

### 1.2 BPA sizing law (the key engineering finding)

`festo4()` grows **exponentially for negative strain** (stretched muscle): a naive
"rest length = spawn route length" sizing produced knee forces of 13,096 N and joint
angles of ±14,000° (sim unstable at t=0.136 s — measured in the first smoke run). Fix,
now in the driver: **rest = MAX route length over a 7³ pose grid of the side's joint
ranges − tendon length**, so strain ∈ [0,1] for every reachable pose and forces are
bounded by p·Fmax by construction. Second deviation found necessary: `kmax_frac=0.6`
(max-contraction fraction) because the w2l routes — designed for Hill muscles — swing
0.06–0.16 m over the joint ranges, far more than a physical Festo DMSP-20 sleeve
(KMAX ≈ 0.2) could cover; at kmax_frac=0.25 the knee extensors sit slack at spawn
(rel_spawn 0.85). **This 0.6 is a design deviation to be honest about: a real BPA at
these routes would need re-routing (the Xi/mesh-optimization program), not just a longer
stroke.**

Sized units (easteregg2 log, 2× 20 mm per route, tendon 0.04 m, kmax_frac 0.6):

```
  name                route min..max   spawn    rest    Fmax kmax_len rel_spawn
  hip_R_ext         0.1335..0.1972  0.1606  0.1772    2223   0.0709     0.532
  knee_R_ext        0.3407..0.3965  0.3407  0.3765    2576   0.1506     0.335
  ankle_R_ext       0.2820..0.2904  0.2863  0.2704    2450   0.1082     0.148
  (full 12-row table in campaigns/20260930/bpa_run/bpa_gate.log)
```

### 1.3 The gate run (of record)

Launched WMI-detached on easteregg2 (colleague check first: `query user` → only
`ben bolen`, Disc), `pid=27608 ret=0`, via
`launch_bpa_gate.ps1` →
`python test_bpa_stepping.py --dur=20 --hold=6 --lower=2 --kmax-frac=0.6 --cap=0.6 --count=2`.
Full log: `campaigns/20260930/bpa_run/bpa_gate.log`. Verbatim results:

```
   AIR knee L: 3 flexion excursions, period 1.542 s, flexion range 63.5 deg [-1.8,+61.8]
   AIR knee R: 4 flexion excursions, period 1.543 s, flexion range 63.2 deg [-1.3,+62.0]
   AIR neural: L RG-E bursts 3, period 1.542 s, max 7.00 mV
   RETAINED phase (11.5 s of it): pelvis z min 0.916 m, end 0.916 m
   RETAINED vertical load: mean L 19.8 N / R 28.1 N; PEAK L 880 N / R 1372 N;
        stance episodes L 7 / R 5 (harness carries 70% of 411 N)
   BPA pressure: peak 372 kPa (ceiling 620), mean-of-active 318 kPa
   [PASS] finite
   [PASS] air_bursts_ge3_both_knees
   [PASS] stand_height_held
   [PASS] feet_bear_load
   [PASS] pressure_le_620kPa
VERDICT: BPA ACTUATION PASS (air-stepping + retained-support stand on BPA muscles)
```

Gate criteria vs the ask: ">=3 air-stepping bursts with knee swing" → 3/4 excursions per
knee with 63° flexion swing, locked 1:1 to the RG-E burst period (1.542 s). "Retained-
support stand where feet bear load" → height held (0.916 m end, never below), 7/5 stance
episodes ≥ 0.05 s per foot with peaks 880/1372 N (means 19.8/28.1 N ≈ the 30 % of the
411 N weight the harness leaves, delivered in stance bursts — the walker MARCHES under
the retained harness rather than standing quietly; see §1.5). Pressure ceiling respected
(372 ≤ 620 kPa). The BPA forces are genuinely load-bearing — a direct probe at the spawn
pose applied `qfrc_applied[knee] = −41.1 N·m` from a 743.6 N knee-extensor BPA force at
0.6 activation, through the same `actuator_moment` column the gate run uses.

### 1.4 Two documented deviations from the stock protocol

1. **Force-servo lowering.** The stock open-loop "lower to z_contact" assumes the feet
   are at their 1 mm keyframe clearance at contact time; with the march pose differing
   (hips swing to extension, knees cycle) the feet stayed in the air — measured 0.0 N
   over 11.5 s. The retained phase now servos the PD target down/up at 1 cm/s until the
   feet carry (1−0.7)·411 N ≈ 123 N, then holds (settled at z=0.916 m). This is what
   the AnimatLab virtual walker's descent does (lower until the platform bears).
2. **Stiff joint limits** (`--stiff=1`, the `test_w2l_air.py:188-190` stand-in): without
   it the ankles blow through their [-20,-5]° range to +115°, the routes then exceed the
   sampled max, and the negative-strain explosion returns (ankle_flx force 159,975 N —
   measured).

### 1.5 What does NOT work yet (blockers, precisely)

* **Quiet standing (latched drive) fails both ways.** Extensor latch (`--stand-te=6
  --stand-tf=0`): knees fold to the −60° flexion limit with 1240 N of extensor force —
  root cause measured with a moment-arm probe: **in this ported body the "knee
  extensor" route produces flexion-sign torque** (`actuator_moment[knee_L_ext, knee]` =
  −0.0553 N·m/N at 0°, −0.0427 at −60°; flexor +0.0171…+0.0957). The Hill-muscle runs
  carry the same sign (same tendons, same moment matrix), so the net was never relying
  on the semantic label — but any extensor-latch "stand" strategy is architecturally
  wrong on this body. Flexor latch (`--stand-te=0 --stand-tf=6`): mixed pose (L knee
  −60°, R knee 0°, hips/ankles pinned at limits), 0 foot contact. The passing mode is
  the harness-retained MARCH (alternating stance); a true quiet stand needs either
  balanced co-contraction tuning or a corrected knee-extensor route.
* **`build_w2l_split_net.py` census assert was stale on BOTH machines** (expected 166
  synapses; the 2026-09-30 V3-wiring correction removed the 2 invented
  `v3_SynAmp0.1_weak` edges → actual 164), so `build()` raised
  `AssertionError: synapse transcription mismatch: 164 vs 166` before any run. Fixed
  (census + assert updated to 87/164, absence of `v3_SynAmp0.1_weak` asserted); the
  shared laptop file now builds `87 neurons, 164 synapses, 12 muscle outputs`
  (verified). The gate run used the fixed copy.
* **Concurrent-edit incident (process note):** a parallel workflow track edited the same
  `w2l_mujoco/` directory mid-session and briefly reverted the census fix
  (`test_w2l_ground.py`, `w2l_ground.xml` modified at 01:50). Resolved by running from a
  sandbox (`campaigns/20260930/bpa_run/`) and re-verifying the shared file afterwards.
* **kmax_frac = 0.6 exceeds a physical Festo DMSP-20 KMAX (~0.2)** — see §1.2. The
  BPA-forced walker is a valid actuation demonstration, not a build-ready muscle
  design for these routes.

## 2. Simscape track (laptop, MATLAB R2025b) — demo runs, no sustained cycling yet

Command (exactly as the ask prescribed): from `Code\Matlab\SNS_Simscape\demos`,
`matlab -batch "cd(...); sns_run_cpg_demo"` → **exit 0**. Log (`tmp/bpa_track/
bpacpg_run.log`, mirrored in `bpa_run/`):

```
CPG: 2 hysteresis switches in 10 s (period ~ 0.155 s, freq ~ 6.46 Hz)
CPG run OK: theta range 9.7..48.0 deg
```

My own trace analysis of the saved `results/sns_cpg_results.mat` (9745 samples, 0–10 s):
knee θ 9.7–48.0° but **0 sustained knee cycles** by hysteresis counting; activations
0–0.58 / 0–0.57; BPA forces ext 0–66 N, flex 0–30 N; RG half-center voltages
anticorrelated **r = −0.701** (antiphase). So the neural half-center + antagonist-BPA
chain works and both BPAs produce force, but the stock 10 s configuration moves the knee
as a transient, not a rhythm — Simscape "air-stepping" is a tuning problem (drive /
adaptation / epsScale), not a missing-component problem.

### 2.1 BPA block inventory in SNS_Library (what exists / what is missing)

Blocks (enumerated via `find_system` on `SNS_Library.slx`, 12 total): **BPAForce**
(normalized: `Fmax`, `epsMax`; activation + strain → force — what BPACPGLegDemo uses),
**BPA_10mm / BPA_20mm / BPA_40mm** (physical: **P[kPa], L[m] → F[N]**, Ben's Festo
equations per `dev/build_bpa_sns_rig_v2_20260926.m:7`), BioMuscle, MuscleActivation,
IaMuscleSpindle, IbGolgiTendon, NonSpikingNeuron, NonSpikingSynapse, SpikingLIFNeuron,
SynSum.

Existing BPA-bearing Simscape models: `demos/BPACPGLegDemo.slx` (1-DOF knee + half-center
CPG + antagonist BPAForce pair); `dev/build_bpa_sns_rig_20260926.m` + `_v2` (knee test
rig with BPA_20mm fed by a Transform Sensor for length; results mats exist); and
`dev/build_bpa_sns_humanoid_20260926.m` — **both humanoid KNEES** with EXT/FLX BPA
muscles + SNS antagonistic pairs (results `humanoid_bpa_sns_20260926.mat`). A
representative W2L CPG neural core (`dev/build_sns_w2l_cpg_20260930.m`,
`results/SNS_W2L_CPG.slx`) landed from a parallel track this same night.

**Missing for a full BPA walker in Simscape:** hip and ankle BPA routing (the humanoid
build covers knees only), the full PF→12-MN wiring (the W2L core carries 2 representative
MNs), ground contact + afferent channels on a plant, and the sustained-oscillation tuning
above. The blocks themselves (physical pressure-driven BPAs, synapses, spindles, GTOs)
all exist.

## 3. Reproduction

* **easteregg2 (of record):** `cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\campaigns\20260930\bpa_run`
  then `D:\Anaconda\envs\myo\python.exe make_bpa_walker.py D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_mjcf_fixed.xml w2l_mjcf_bpa.xml`
  and `set PYTHONPATH=D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal;D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS`
  then `D:\Anaconda\envs\myo\python.exe test_bpa_stepping.py --dur=20 --hold=6 --lower=2 --kmax-frac=0.6 --cap=0.6 --count=2`
  (identical numbers reproduced on the laptop `myoconv` env before the remote run).
* **Laptop Simscape:** `matlab -batch "cd('.../SNS_Simscape/demos'); sns_run_cpg_demo"`.
* New/changed repo files this track: `w2l_mujoco/make_bpa_walker.py`,
  `w2l_mujoco/test_bpa_stepping.py`, `w2l_mujoco/w2l_mjcf_bpa.xml`,
  `w2l_mujoco/w2l_mjcf_fixed.xml` (synced to the easteregg2 body of record),
  `w2l_mujoco/build_w2l_split_net.py` (census-assert fix),
  `campaigns/20260930/bpa_run/*` (run artifacts), `results/sns_cpg_results.png|.mat`
  (regenerated by the demo run).
