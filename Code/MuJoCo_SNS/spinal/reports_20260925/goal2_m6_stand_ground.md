# Goal 2 / Milestone 6 — GROUND: standing balance + contact-driven walking attempt

**Date:** 2026-09-25 (EB475WS4, unattended campaign) · **Status: partial — as the ask anticipated:
(a) unrigged stand FAILS for a measured, structural reason (the COM is outside the support
polygon — no open-loop controller can stand this pose); a harness-supported stand holds
(27% weight carried, reported); (b) ground walking does NOT walk — fully characterized
(fall at 0.89 s free; supported in-place march at the autonomous RG period with ZERO heel
loading), top 3 blockers named with evidence.**

**The ask:** (a) afferented walker on the ground, real heel/toe contact, pelvis released,
COM sway + tilt over ≥ 10 s; light support rig allowed but stiffness must be reported
(rigged stand = partial, unrigged = full). (b) contact-driven ground walking attempt; full
honesty — if it does not walk, characterize what it does and name the top 3 blockers with
evidence.

**Artifacts (all under `Code\MuJoCo_SNS\spinal\w2l_mujoco\` unless noted; no protected file
touched):**
- `test_w2l_ground.py` — **the M6 gate** (both phases). Writes `w2l_ground.xml`
  (NEW on disk) = `w2l_mjcf_fixed.xml` (M3 axis-fixed body) + a FREE joint on the Root
  pelvis — the ask's "release the pelvis". Source XMLs never modified.
- Probe scripts + outputs: `reports_20260925\tmp\m6_probe_fall.py` (passive first-0.6 s
  trace), `tmp\m6_probe_closed.py` (closed-loop latch trace with ctrl/forces),
  `tmp\m6_static_margin.py` (**the decisive static-balance margin measurement**).
- Gate logs (`reports_20260925\logs\`): `test_w2l_ground_stand_final.log` (**canonical
  gate a**), `test_w2l_ground_walk_free.log` (**canonical gate b, free**),
  `test_w2l_ground_walk_rig1.log` (**gate b, harness-supported**), plus the smoke trail
  `test_w2l_ground_smoke1.log` and un-logged foreground smokes 2–5 (their numbers are
  quoted in §3 where relevant).

Env: `C:\Users\Ben Bolen\.conda\envs\myo\python.exe`, cwd `w2l_mujoco\` (mujoco 2.3.7,
sns-toolbox 1.5.2; CONDA_PREFIX set before import per the mujoco skill). No pip installs.
Protected set untouched: `runner.py`/`build_network.py`/`params.py`, `spinal_run.npz`,
optuna studies, `reports_20260923/24`. M3/M4/M5 artifacts untouched (new files only).

---

## 1. Setup (what "on the ground, pelvis free" means here)

- **Body:** `w2l_ground.xml` = M3's axis-fixed `w2l_mjcf_fixed.xml` with `<freejoint
  name="root"/>` added to the Root. The aproj rest pose is kept verbatim EXCEPT two
  documented, runtime/parameter-level modifications:
  1. **`--drop` 0.026 m** (default): the free-joint z is lowered so the rest pose STARTS
     grounded (the aproj plates hover 2.6–3.1 cm, M1 §2.8; Li's model also starts
     feet-on-ground, M2 §1). Measured t=0 plate heights after the drop
     (`m6_static_margin.py`): heel plates +0.0038/+0.0052 m; toe plates **stab 4 cm below
     ground** (rotated plate corners — the source pose is itself inconsistent with its
     ground plane; MuJoCo pushes them out in the first 20 ms).
  2. **Ankle limits runtime-widened to [−20°, +5°]**: the transported range [−20°, −5°]
     EXCLUDES the aproj's own rest pose (ankle = 0°) — the M3 §5 transport inconsistency.
     Measured consequence if left (`m6_probe_fall.py`): the stiff limit solver slams both
     ankles to −4.8° within 20 ms of t=0 (max|qacc| 5.6e3) and the run never gets past
     the transient. The −20° PF bound is kept; widening is runtime-only, same family as
     the M3 joint-damping stand-in.
- **Contact encoders (the ask's "heel/toe from MuJoCo contact forces per the rules"):**
  the aproj's four dedicated contact plates are the sensors (heel S = `foot_S_contact`,
  toe S = `toe_S_contact`). Per step the plate normal forces are accumulated via
  `mj_contactForce` and encoded with a saturating linear map
  `I_port = ctn_amp · clip(F/cref, 0, 1)` (defaults ctn_amp 4.0 nA = the M5
  gate-b saturating/strike calibration; cref 50 N ≈ 12% body weight). Same PORT wiring
  as M5 — only the scripted scheduler is replaced by real forces. Ib group per M5
  (`ibnA`, run at 0 — see §3).
- **Support rig (stiffness reported, as the ask requires):** world-frame PD on the Root
  via `xfrc_applied`: horizontal leash (kxy 2000 N/m) + **soft vertical harness
  (kz 2000 N/m — ~25% of the 411 N weight at a 5 cm sink, ~0 near upright, so feet stay
  loaded by design)** + tilt assist (krot 400 N·m/rad), all × S (the rig-scale number),
  critical-rate damping, anchored at the pose captured at t = settle = 0.4 s, S blended
  1.0 → target over 0.6 s. **S=0 = no rig at all.** The gate prints the MEASURED
  weight-bearing ("harness carries X% weight") per run.
- **Drive regimes:** stand = `latch` (tonic te=3 on BOTH RGs, no F tonic, no antiphase
  kickoff → steady extension posture; the first smoke drove neither side — kickoff alone
  latches only L in the split net — fixed); walk = the M5 config (te=3/tf=4, antiphase
  kickoff) with real contact replacing the scripted scheduler. Body stand-ins carried
  from M3/M4: joint damping 3.0 N·m·s/rad, stiffer limits, ctrl cap 0.5, plus one new
  documented knob `--acap` (ankle-specific ctrl cap; see §3).

## 2. GATE (a) — standing balance (final run, `logs\test_w2l_ground_stand_final.log`)

Command: `C:\Users\Ben Bolen\.conda\envs\myo\python.exe test_w2l_ground.py --phase=stand
--dur_stand=12 --acap=0.15 --ibnA=0` (knobs `--acap 0.15` = ankle PF stance trim, see
§3; `--ibnA 0` = the M5 air-latch positive loop disabled — with Ib on, extensor force
feeds RG-E/PF-E and latches the net rigid even in air, M5 §3.1).

```
== M6 gate (a): GROUND STANDING BALANCE, afferented walker, pelvis FREE, 12 s ==
   knobs: sregime=latch ctn_amp=4.0 cref=50.0 ibnA=0.0 damp=3.0 settle=0.4 (+0.6 s blend) drop=0.026 m
   (rig at S=1: kxy=2000 N/m, kz=2000 N/m [~25% weight at 5 cm sink], krot=400 N m/rad)
   fall: pelvis z<0.6 m or tilt>45.0 deg
   rig S=0.00: finite=True fall=t=0.92 s | sway_x std/range 2.3/24.0 cm | sway_y std/range 0.3/2.7 cm
     | COM z mean/min 0.103/0.091 m | tilt mean/max 72.5/81.8 deg | heel duty L/R 0.00/0.00 | harness carries +0% weight
   rig S=0.25: finite=True fall=t=1.02 s | sway_x std/range 3.6/12.3 cm | sway_y std/range 1.4/6.4 cm
     | COM z mean/min 0.247/0.224 m | tilt mean/max 41.2/43.7 deg | heel duty L/R 0.00/0.00 | harness carries +25% weight
   rig S=0.50: finite=True fall=t=1.23 s | sway_x std/range 1.8/6.8 cm | sway_y std/range 0.4/1.5 cm
     | COM z mean/min 0.420/0.410 m | tilt mean/max 22.3/23.5 deg | heel duty L/R 0.00/0.00 | harness carries +28% weight
   rig S=1.00: finite=True fall=no | sway_x std/range 0.8/3.3 cm | sway_y std/range 0.2/0.8 cm
     | COM z mean/min 0.533/0.530 m | tilt mean/max 10.5/11.0 deg | heel duty L/R 0.00/0.00 | harness carries +27% weight
   [PARTIAL] unrigged_stand
M6 VERDICT: PARTIAL (rigged stand only)
```

Analysis window = 11 s (t ∈ [1.0, 12.0]) — meets the ask's ≥ 10 s.

**Verdict: PARTIAL, as the ask's own rubric anticipated** ("a rigged stand is a partial
result, an unrigged stand is a full one"). The unrigged stand FAILS — and the failure is
STRUCTURAL, not a tuning gap:

**THE COM IS OUTSIDE THE SUPPORT POLYGON.** `tmp\m6_static_margin.py` (measured on the
grounded rest pose, all support geoms = 4 plates + feet + toes): support x-range
[−3.6555, −3.4852] m (17.0 cm feet), **COM x = −3.4544 m = 3.09 cm IN FRONT of the front
edge** (margin to rear edge +20.1 cm; y margin +10.8 cm). A rigid statue with this COM
topples forward; no open-loop muscle pattern can stand it — the CoP is bounded by the
polygon, so there is no torque solution. This is consistent with the source's own
history: the aproj ships the Root with `<Freeze>True</Freeze>` — the AnimatLab walker was
NEVER a free-standing model (AGENTS/2026-09-15 note; M3 §6), and Li's variant stood on his
own proportions + LinearHill muscle damping B (M2 §5). Measured dynamics agree: free
runs topple forward over the toes at t ≈ 0.9 s in every regime tried (heel plate force
decays 398 → 0 N while toe force grows — the exact forward-pivot signature,
`m6_probe_fall.py`). Fix belongs to M1 open item #4 (a standing-pose solve that puts the
COM over the polygon) and/or a balance layer (vestibular — Ben's stage-4 work in the main
spinal project), not to more knob-turning here.

**What the S=1.0 hold IS:** an upright, harness-supported squat (tilt 10.5° mean, sway
3.3 cm range over 11 s, no fall) with the harness carrying a measured **27% of body
weight**; COM z 0.53 m vs 0.85 standing, and heel duty 0.00 — at that sink depth the
feet are unloaded, so honestly it is a harness hold, not a foot-loaded stance. Reported
as-is.

## 3. GATE (b) — contact-driven ground walking attempt (20 s each)

Canonical commands (exit 1 = does not walk):
`... test_w2l_ground.py --phase=walk --dur_walk=20 --acap=0.15 --ibnA=0 --rig=0`
(log `test_w2l_ground_walk_free.log`) and the same with `--rig=1.0`
(log `test_w2l_ground_walk_rig1.log`).

**FREE (S=0) — paste from `test_w2l_ground_walk_free.log`:**
```
   finite=True | fall=t=0.89 s | pelvis z min 0.092 m (start 0.993) | tilt max 136.7 deg | harness carries +0% weight
   cycles (heel strikes >20 N, window 19.0 s): L 0 (interval nan s) | R 0 (interval nan s)
   heel duty L 0.00 R 0.00 | forward disp -3.311 m (-0.174 m/s)
   flexion ranges (deg): hip L 37.2 R 33.9 | knee L 68.9 R 68.4 | ankle L 33.8 R 33.8
   flexion min..max: hip L [-27.8,+9.4] knee L [-4.5,+64.4] ankle L [-11.0,+22.9]
   neural: L RG-E bursts 14 (period 1.357 s), R RG-E bursts 14; L/R RG-E r -0.751
   contact encoder evidence: heel SN max 0.00 mV, toe SN max 3.76 mV, heel port current max 0.00 nA
   walk verdict: FALLS/NO-GAIT
```
**What it does:** topples forward at t = 0.89 s (blocker 1 below), ends prone; the RG
keeps firing and the legs keep flexing at full amplitude while lying (the ranges above
are prone-flailing, not stepping) and the body slides −3.31 m. Zero heel strikes, zero
duty. NOT a shuffle, NOT a gait — a fall, honestly.

**HARNESS-SUPPORTED (S=1.0) — paste from `test_w2l_ground_walk_rig1.log`:**
```
   finite=True | fall=no | pelvis z min 0.660 m (start 0.993) | tilt max 11.7 deg | harness carries +39% weight
   cycles (heel strikes >20 N, window 19.0 s): L 0 (interval nan s) | R 0 (interval nan s)
   heel duty L 0.00 R 0.00 | forward disp -0.002 m (-0.000 m/s)
   flexion ranges (deg): hip L 34.1 R 28.2 | knee L 61.2 R 55.5 | ankle L 33.0 R 15.2
   flexion min..max: hip L [-25.7,+8.4] knee L [+1.9,+63.1] ankle L [-10.7,+22.3]
   neural: L RG-E bursts 14 (period 1.357 s), R RG-E bursts 14; L/R RG-E r -0.751
   contact encoder evidence: heel SN max 0.00 mV, toe SN max 4.00 mV, heel port current max 0.00 nA
   walk verdict: FALLS/NO-GAIT   (fall=no on the metrics line — the label is coarse; see below)
```
**What it does:** stays upright (no fall, tilt 11.7° max) and MARCHES IN PLACE — knee
61/56°, hip 34/28° flexion excursions at the neural period — but zero forward motion
(−0.002 m over 19 s) and, the key measurement, **the heels never load: heel SN max
0.00 mV, heel port current max 0.00 nA, duty 0.00 in 19 s.** The toe plate is the only
plate that ever touches (toe SN max 4.00 mV). So the contact-driven stance-reset
architecture is OPEN at the heel: what remains on the ground is the AUTONOMOUS central
rhythm, not a contact-driven gait. (The verdict label "FALLS/NO-GAIT" is the script's
coarse 3-way bucket; the run itself is a supported march, `fall=no` printed on the same
line.)

## 4. Top 3 blockers (evidence per the ask)

1. **The rest pose is statically unbalanceable — the COM sits 3.09 cm outside the front
   edge of the support polygon** (measured: `tmp\m6_static_margin.py`, support x
   [−3.6555, −3.4852] vs COM −3.4544). Consequence measured: every free ground run
   topples forward at 0.89–0.92 s (gate a S=0; gate b free) before any stepping can
   develop; heel force decays 398→0 N while toe force grows (the forward-pivot
   signature, `m6_probe_fall.py`). The aproj's own `Freeze=True` confirms the source was
   never meant to free-stand. → Fixes: M1 open item #4 standing-pose solve (COM over the
   polygon), or a balance layer; nothing in the current M6 scope can fix it.
2. **The ankle/foot system cannot hold or exploit stance.** (i) The transported ankle
   range [−20°,−5°] excludes the aproj's own rest pose (0°) — measured limit slam at
   t=0, max|qacc| 5.6e3, ankles pinned at −4.8° in 20 ms (`m6_probe_fall.py`; M3 §5
   flagged the inconsistency); (ii) with the widened range, the open-loop extension
   latch drives both plantarflexors at the air-walk cap (ctrl 0.50, ankle PF force
   1483 N) and lifts the heels within 100 ms — plate forces heelL 40→0 N, toe R 423 N,
   pelvis 0.967→0.10 m, knees fold to the −60° limit under gravity
   (`m6_probe_closed.py`) — hence the new `--acap` ankle-stance-trim knob; (iii) the
   feet are 17 cm plates whose toe corners already stab 4 cm below the plane at the
   grounded rest pose. Missing muscle damping B (M1 §3, LinearHill 400–800 N·s/m
   unported) removes the element that damps exactly these transients in AnimatLab.
3. **The RG is open-loop and provably ignores the ground, so "contact-driven" never
   engages.** In BOTH ground runs the neural signature is IDENTICAL to the M4 air gate:
   14/14 RG-E bursts at period 1.357 s, L/R r −0.751 (M4 report: 13/13 at
   1.357/1.358 s, r −0.750 — same rhythm within window-length resolution). The heel
   reset (heel IN → InE/InF + PF-layer INs, g 0.5 — Ben's rules, built in M5) requires
   heel loading, and heel loading NEVER occurs on the ground as built (heel SN max
   0.00 mV in 19 s, both runs). The causal chain is: blocker 1 tips the body before
   stance establishes; blocker 2 keeps heels unloaded even when supported; the M5
   heel/toe gains (0.5/1.0) are therefore untestable on the ground until 1+2 are fixed.

## 5. Honest notes / deviations

1. **`--acap` (new knob, documented):** ankle-specific ctrl cap. The uniform 0.5 cap is
   an M3 air stand-in; on the ground it made the stance posture impossible (blocker
   2-ii). The gate runs above use acap 0.15 for BOTH phases, printed in every log.
   Even with it, nothing free-stands (blocker 1 is structural).
2. **`--ibnA 0` in all final runs:** the M5 air-latch positive loop (force→Ib→RG-E) is
   destabilizing on the ground too; the M5 report itself required ground recalibration
   of the encoder scale (M5 §3.1). With it at 0, the Ib cells are built (conditional
   topology intact, 97/188 census unchanged) but receive no current.
3. **The stand "latch" mode is a run regime, not a new net:** te on both RGs, tf=0,
   no antiphase kickoff — the steady-extension posture. The first smoke accidentally
   drove neither side (kickoff alone latches only L in the split net) — caught and
   fixed in smoke 2; the shipped gate drives both.
4. **Settle device** (0.4 s rig at S=1 + re-anchor + 0.6 s blend) is a runtime gate
   device, reported here; it cannot mask falls because every metric window starts
   AFTER the blend (t ≥ 1.0 s) — and the unrigged runs still fall at ~0.9 s inside the
   window.
5. **Not run (no budget):** a posture-solved start pose (M1 #4), ankle-range surgery in
   the generator itself (runtime widening shipped instead), muscle-damping B port,
   vestibular balance layer. None of these are M6 scope; each is listed where it bites
   (§2–§4).
6. Determinism: same-binary reruns reproduce these numbers; the standing cross-platform
   chaos caveat applies (AGENTS 2026-09-13).

## 6. PORT-LEVEL VERDICT — the six milestones in one paragraph

**The AnimatLab 2-layer walker port to MuJoCo-SNS is: body faithful, air-neural PROVEN,
ground NOT yet walking — with the ground failure now reduced to three measured, named
causes rather than a mystery.** M1 transported the physical model exactly (pose dev
4.5e-7 m, mass dev 0, joint ranges ≤2e-10 rad) but shipped vertical knee/ankle axes —
M3 caught and fixed them (`fix_joint_axes.py`, 8/8 lateral) after M2's closed loop
failed for that hidden reason; the lesson "axis directions need a gate" is banked. M3/M4
proved the neural stack end-to-end: the 2023 original and Ben's split-RG (his
"disconnect the LH RG, mirror it, couple with commissurals" ask) both air-step on the
axis-fixed body, with the ablation proving each leg runs on its OWN RG (1.027 s free vs
1.357 s locked, causal). M5 wired Ben's contact rules (heel = stance reset at the PF
layer, toe = DF inhibition, Ib load) with a causal heel-reset gate PASS in air. M6
shows the remaining gap is the GROUND INTERFACE, in order: the source rest pose cannot
statically stand (COM 3.09 cm outside the support polygon — the aproj always froze the
pelvis), the ankle/foot interface cannot hold stance (transported range excludes the
rest pose; PF latch lifts the heels in 100 ms; muscle damping B unported), and until
those two are fixed the heel contact encoder never fires, so the contact-driven layer
built in M5 stays silent and the ground behavior reduces to the autonomous M4 rhythm
(identical 1.357 s / r −0.751 signature on the ground). The port's next session should
be, in order: (1) a standing-pose solve (M1 #4) that puts the COM inside the polygon,
(2) the ankle-range + muscle-damping fidelity fixes in the generator, (3) THEN re-run
this gate — `test_w2l_ground.py` is written to be re-run unchanged, and its contact
encoders will finally have heel signals to reset on. All six reports:
`goal2_m1_physical_model.md`, `goal2_m2_li_architecture.md`, `goal2_m3_w2l_cpg_air.md`,
`goal2_m4_rg_split.md`, `goal2_m5_afferents.md`, this file.

## 7. File manifest (this milestone)

- NEW `w2l_mujoco\test_w2l_ground.py` (the M6 gate; writes `w2l_ground.xml`)
- NEW `w2l_mujoco\w2l_ground.xml` (axis-fixed body + free pelvis; generated)
- LOGS `reports_20260925\logs\test_w2l_ground_{stand_final,walk_free,walk_rig1}.log`
  (canonical) + `test_w2l_ground_smoke1.log` + the foreground smoke trail
- tmp `reports_20260925\tmp\m6_probe_fall.py`, `m6_probe_closed.py`,
  `m6_static_margin.py` (+ their inline outputs quoted above)
- No protected file touched; M3/M4/M5 artifacts unchanged.
