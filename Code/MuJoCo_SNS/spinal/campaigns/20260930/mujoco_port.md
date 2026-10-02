# MuJoCo port of the AnimatLab walkers — Li + W2L drop-protocol gates (2026-10-01)

Track: land the AnimatLab walkers in MuJoCo on easteregg2 under Ben's protocol
(suspended in air -> lowered onto the platform WITH the harness retained,
Ben 2026-09-30). All runs on easteregg2 (`D:\Anaconda\envs\myo\python.exe`,
repo `D:\GitHub\Bipedal_Robot`), logs in `Code\MuJoCo_SNS\spinal\w2l_mujoco\`;
this report written on the laptop with the remote logs fetched verbatim.

## 1. Colleague check (shared-box rule)

`query user` over ssh before any launch — only Ben's own disconnected session:

```
 USERNAME              SESSIONNAME        ID  STATE   IDLE TIME  LOGON TIME
 ben bolen                                 2  Disc      1+00:53  9/16/2026 9:33 AM
```

No other user on the box -> launches authorized (Ben's full-tilt rule). All
jobs WMI-detached (`Invoke-CimMethod Win32_Process.Create`) at BelowNormal.

## 2. Pre-existing breakages found and fixed (the gates could not run without these)

Three breakages, all exposed by simply trying to run the gates on the current
easteregg2 working tree. Each was reproduced on the UNMODIFIED tree before
any edit of mine.

### 2a. Stale synapse censuses after the 2026-09-30 V3-wiring correction

`B.build(comm=1.0)` failed on the unmodified easteregg2 tree:

```
AssertionError: synapse transcription mismatch: 164 vs 166
  (build_w2l_split_net.py line 338; run 2026-10-01 via ssh)
```

Cause: the V3 correction (2026-09-30, Ben-approved — V3 -> contra RG ext IN
only; the two `v3_SynAmp0.1_weak` HC edges were a generator invention) removed
2 commissural synapses from the split template, but the hardcoded count asserts
had not been updated with it (builder mtime 09/30 12:35 AM; the aff census in
`build_w2l_aff_net.py:236` hardcoded the same stale 166 baseline). A census
probe (`tmp\li_fetch\probe_aff_census.py`) confirmed the per-tag counts are
otherwise exactly as asserted — only the totals were stale.

Fix (commented in both files): split coupled baseline 87/**164** (was 87/166),
aff baseline **164**+22 (was 166+22), and the split tag assert now requires
`comm_v3 == 2 and v3_to_contra_InE == 2` (the `v3_SynAmp0.1_weak == 2`
assert is gone with the edge). Verified after the fix (laptop, same builders
hash-synced to easteregg2):

```
SPLIT OK 87 neu 164 syn
```

### 2b. `write_ground_xml()` duplicate root joint after the MJCF rebuild

The 2026-09-30 MJCF rebuild made `w2l_mjcf_fixed.xml` ship the root freejoint
at the source (`w2l_mjcf_fixed.xml:24` = `<freejoint name="root"/>`); the M6
ground test still INJECTED a second one:

```
ValueError: Error: repeated name in joint array, position 1
Object name = root, id = 1
```

Fix: `write_ground_xml()` is now a read-only copy + header comment with an
assert that the source freejoint exists (test_w2l_ground.py).

### 2c. `test_w2l_air.py` FAILS `airborne` on the freejoint-ed body (NOT fixed — reported)

Re-ran `test_w2l_air.py` on easteregg2 after the census fix (log
`w2l_air_rerun.log`): the neural rhythm reproduces exactly the numbers on
record (18 RG-E bursts, period 1.027 s, hip antiphase r = -0.414), but the
gate now reads:

```
   ground contacts: 28890 (must be 0) | leg-leg self contacts: 22982
   ...
   [FAIL] airborne
VERDICT: FAIL airborne
```

The rebuilt body is no longer welded to the world, so the air test's
suspended-body assumption broke with the same rebuild. Out of scope for this
track (the ask covers the Li sweep + the ground-gate conversion); flagged as
the next fix if the air gate is to be re-cited. Note the inherited STATE's
"test_w2l_air PASS (18 bursts, 1.027 s, -0.414)" predates the freejoint
rebuild of `w2l_mjcf_fixed.xml`.

### 2d. Machine drift synced (laptop <- easteregg2)

The laptop's `w2l_mujoco` tree was behind the easteregg2 working tree on two
inputs (MD5s, `certutil`/`certutil`-equivalent both sides):

| file | laptop (before) | easteregg2 | action |
|---|---|---|---|
| `build_w2l_split_net.py` | 8a86161c… | f0841935… (V3-corrected) | pulled east->laptop |
| `w2l_source_dump.json` | 1fbb7d80… | 2f14f866… | pulled east->laptop |

All other `.py`/`.xml` inputs in `w2l_mujoco\` were already hash-identical.
After my fixes, the three changed files were pushed back east with
`scp -O` (legacy SCP protocol — plain `scp`/SFTP rejects the spaced username
"ben bolen"; `copy con`/stdin pipes hang remotely) and hash-verified:
`build_w2l_split_net.py` 08fdc01e…, `build_w2l_aff_net.py` 4a4ccc7d…,
`test_w2l_ground.py` 8ed4e4f0… — identical on both machines.

## 3. Li contact-driven CPG under the harness-retained drop protocol — verdict table

Protocol per `test_li_stepping.py` (air-hold 3 s at keyframe z + clear ->
lower 1 s -> retained support at SUPPORT fraction + PD at contact height;
knobs parsed as `--key=value`). Success criterion (the ask): stance episodes
> 0 on BOTH feet and a real period. Seven runs total: the three earlier
blocks in `li_drop.log` + the six-knob sweep (cap six, all launched via WMI,
all completed; each knob exactly once, otherwise defaults). Every verdict
line, verbatim from the logs:

| run | pelvis height min/end (m) | L / R stance episodes | period (s) | verdict |
|---|---|---|---|---|
| li_drop.log block 1 (earlier drop attempt) | 0.099 / 0.100 | 0 / 0 | nan | not yet |
| li_drop.log block 2 (earlier drop attempt) | 0.138 / 0.140 | 0 / 0 | nan | not yet |
| default (support=0.7) — li_drop.log block 3 | 1.003 / 1.004 | 0 / 0 | nan | not yet |
| `--support=0.85` | 1.005 / 1.007 | 0 / 0 | nan | not yet |
| `--support=0.5`  | 1.001 / 1.004 | 0 / 0 | nan | not yet |
| `--support=1.0`  | 1.006 / 1.007 | 0 / 0 | nan | not yet |
| `--clear=0.02`   | 1.005 / 1.006 | 0 / 0 | nan | not yet |
| `--clear=0.06`   | 1.002 / 1.003 | 0 / 0 | nan | not yet |
| `--hold=2.0`     | 1.003 / 1.005 | 0 / 0 | nan | not yet |

Representative quoted lines (`li_sw_support05.log`):

```
== M2 gate: Li CPG on the M1 MuJoCo body, 20 s ground run ==  knobs: defaults
  pelvis height: min 1.001 m, end 1.004 m
  mean forward speed: +0.132 m/s   (Li ref +0.638)
  L_foot: 0 stance episodes, period nan s (Li ref 1.3), stance nan s, duty 0.00
  R_foot: 0 stance episodes, period nan s (Li ref 1.3), stance nan s, duty 0.00
  hip L: -3.2..+14.1 deg  knee L: -66.2..-11.0 deg  ankle L: -5.4..+70.6 deg
  hip R: -17.9..-9.5 deg  knee R: -33.6..+6.9 deg  ankle R: -22.0..+69.7 deg
VERDICT: not yet
```

(The `knobs: defaults` print is a known logging gap: the script pops the
protocol knobs before printing, so each sweep run wrote to its own
`li_sw_<knob>.log` for attribution.)

**Diagnosis (from the sweep, no extra runs needed).** The pelvis NEVER
descends to the platform under any support fraction 0.5–1.0 or clearance
0.02–0.06 m: it hangs on the retained PD at `z_contact − (1−SUPPORT)·W/K_H`
(the support=0.5 run's min, 1.001 m, is exactly keyframe z 1.0078 − 7.5 mm =
the predicted zero-contact sag). The spawn keyframe z (1.007846,
`w2l_mjcf_fixed.xml:148`) has the feet 1 mm clear **in the keyframe pose**,
but the CPG runtime pose runs knees at −64…−67° flexion — the feet stay well
clear of the platform even with the pelvis below the keyframe plane. The
harness carries the body; the heel/toe contact SNs never fire; the
contact-driven CPG never engages. Meanwhile the two earlier release-style
runs (0.099 / 0.138 m) show the opposite failure: without retained support
the body buckles through the soft joint limits (the known no-Kse/Kpe
body-fidelity blocker, probe_hold.py).

**Next levers for Li (in order):** (1) lower the retained PD target BELOW
keyframe z by the knee-flexion margin (a `--z-offset`-style knob; the sweep's
knob set had no such knob) or make the lowering contact-gated (descend until
the heel SN fires, then hold); (2) the standing body-fidelity fix from the
09-30 probe (stiff solimp on the 8 hinges + Kpe as tendon stiffness), which
is what lets the stance leg bear load without the harness; (3) only then
re-judge antiphase stepping against Li's 1.3 s reference.

## 4. W2L ground gate converted to Ben's drop protocol

`test_w2l_ground.py` changes (laptop copy edited, pushed east, hash-verified):

- The spawn-time grounding (`d.qpos[jadr+2] -= drop`, feet-down at t=0) is
  REMOVED. The rest pose spawns verbatim (plates hover 2.6–3.1 cm per the M1
  report), as AnimatLab does — suspended in air.
- New virtual-walker vertical harness on the root free-joint z (the
  `test_li_stepping.py` runner-rig idiom: feedforward weight + stiff PD,
  K_V = 4·W/drop, not scaled by the rig S): hold at the AIR height for
  `--hold=3.0` s while the split-net air-steps -> linear descent over
  `--lower=1.0` s to `z_air − drop − sink` (0.99298 − 0.026 − 0.01 =
  0.95698 m) -> feedforward retained at `--wsup=0.7` of body weight + PD at
  contact height for the rest of the run (drop WITH the harness).
- The old soft COM-z kz spring (2000 N/m — could never hold the body in air;
  it would need a 20 cm sag to carry 411 N) is replaced by that channel. The
  horizontal COM leash + tilt assist stay, scaled by `--rig` (walk-phase
  default now 1.0 = retained).
- Contact encoders (heel/toe plate force -> saturating current; Ib from
  extensor force) and all gate criteria unchanged. Metric windows now start
  at hold+lower+0.5 s (contact is impossible earlier by construction).

Pre-flight: `py_compile` + a 6 s `--phase=walk` smoke on the laptop (myoconv,
same mujoco 2.3.7) ran the full protocol path clean:

```
   finite=True | fall=no | pelvis z min 0.951 m (start 0.993) | tilt max 7.7 deg | harness carries +49% weight
   cycles (heel strikes >20 N, window 1.5 s): L 0 (interval nan s) | R 0 (interval nan s)
   contact encoder evidence: heel SN max 0.00 mV, toe SN max 2.32 mV, heel port current max 0.00 nA
```

i.e. the body descends, the platform is reached (toe SN fires; harness <
wsup because the feet take part of the load) — the W2L conversion does NOT
have Li's never-touches failure, because drop+sink lower the target ~5 cm
below the keyframe contact plane.

**Full gate on easteregg2** (`test_w2l_ground.py`, defaults: --phase=both,
4 stand rigs + 20 s walk, WMI-detached): see the verdict block appended in
section 4a below (filled from `w2l_ground_protocol.log` when the run
completed).

## 4a. W2L ground gate — verdict

Full gate on easteregg2 (defaults `--phase=both`: 4 stand rigs x 12 s +
20 s walk; `w2l_ground_protocol.log`; first stand row and the walk block
verbatim, the three middle stand rows abridged to their differing fields —
full log in `tmp\li_fetch\w2l_ground_protocol.log` on the laptop):

```
   rig S=0.00: finite=True fall=no | sway_x std/range 0.5/1.9 cm | sway_y std/range 0.0/0.1 cm | COM z mean/min 0.851/0.851 m | tilt mean/max 2.0/2.1 deg | heel duty L/R 0.00/0.00 | harness carries +50% weight
   rig S=0.25: finite=True fall=no | sway_x std/range 0.0/0.0 cm | ... | tilt mean/max 0.6/0.6 deg | heel duty L/R 0.00/0.00 | harness carries +50% weight
   rig S=0.50: finite=True fall=no | ... | tilt mean/max 0.4/0.4 deg | heel duty L/R 0.00/0.00 | harness carries +50% weight
   rig S=1.00: finite=True fall=no | ... | tilt mean/max 0.2/0.2 deg | heel duty L/R 0.00/0.00 | harness carries +50% weight
   balance metrics table: 4 rig scales
   [PASS-FULL] unrigged_stand

== M6 gate (b): CONTACT-DRIVEN GROUND WALKING attempt, real heel/toe contact forces, 20 s, rig S=1.0 ==
   finite=True | fall=no | pelvis z min 0.951 m (start 0.993) | tilt max 8.2 deg | harness carries +49% weight
   cycles (heel strikes >20 N, window 15.5 s): L 0 (interval nan s) | R 0 (interval nan s)
   heel duty L 0.00 R 0.00 | forward disp +0.001 m (+0.000 m/s)
   flexion ranges (deg): hip L 33.9 R 29.3 | knee L 63.3 R 66.8 | ankle L 38.3 R 36.9
   neural: L RG-E bursts 11 (period 1.539 s), R RG-E bursts 10; L/R RG-E r -0.776
   contact encoder evidence: heel SN max 0.00 mV, toe SN max 2.57 mV, heel port current max 0.00 nA
   walk verdict: FALLS/NO-GAIT

M6 VERDICT: stand=FULL; walk=FALLS/NO-GAIT
```

Reading (honest): the protocol conversion WORKS mechanically — the walker
air-steps in suspension, is lowered, stays UP (no fall anywhere, tilt
<= 8.2 deg), the split net keeps a clean antiphase rhythm through the drop
(L/R RG-E 11/10 bursts, period 1.539 s, r = -0.776), and the feet DO reach
the platform (toe SN 2.57 mV; harness carries only ~49% so the feet bear
the rest). The stand gate passes its printed criteria at every rig scale
including S=0 — under the retained protocol harness (~50% weight), which is
exactly Ben's "drop WITH the harness" semantics. What fails is locomotion:
heel plates NEVER load (heel SN 0.00 mV, heel duty 0.00, 0 heel strikes) —
all load goes through the toes — so there are no heel-driven contact
cycles and no forward progression (+0.001 m over 15.5 s). Same terminal
symptom family as Li (afferents starved at the foot), but a different
mechanism: W2L reaches the platform, it just bears on the toe plates only,
with the knee cycling 63-67 deg flexion range instead of a stance-extended
leg.

## 5. Next blockers (honest list)

1. **Li feet never reach the platform under the retained harness** — the PD
   target must go below keyframe z by the knee-flexion margin (z-offset or
   contact-gated lowering); no knob in the swept set does this.
2. **Body fidelity (both walkers): no elastic load path.** Li's LinearHill
   Kse/Kpe is absent (rigid tendons, MuJoCo 2.3.7 muscles), and soft joint
   limits yield under the 411 N body weight — probe_hold.py showed
   extensors-only hold collapses through knee/ankle limits 0.5 s after
   release. Stiff solimp on the 8 hinges + Kpe-as-tendon-stiffness (goal-3
   §3.2 recipe: >= 1 kg virtual tendon mass or 0.5 ms timestep) is the
   structural fix.
3. **test_w2l_air.py airborne gate broken by the 09-30 freejoint rebuild**
   (28890 ground contacts; §2c) — needs its root pinned/lifted for air runs
   before its PASS can be cited again.
4. **W2L loads through the toes only** (heel SN 0.00 mV, heel duty 0.00,
   toe SN 2.57 mV): no heel strikes -> no contact-driven cycles, no forward
   progression. Levers: larger `--sink` (press the heel plates past their
   3.1 cm hover), check which rebuilt foot geom actually bears in the
   latch posture, and/or the ankle PF trim (the M6 `acap` knob exists for
   exactly this) to stop the toe-first landing.
5. Census/log hygiene: the Li script's `knobs:` print drops the protocol
   knobs (sweep attribution relied on per-run log files); the laptop tree
   was 2 files behind easteregg2 (§2d) — worth committing the pair together
   per the standing UNCOMMITTED-on-both-machines note.

## 6. Artifact trail

- Logs (easteregg2, `Code\MuJoCo_SNS\spinal\w2l_mujoco\`): `li_drop.log`
  (3 blocks), `li_sw_{support05,support085,support100,clear002,clear006,
  hold20}.log`, `w2l_air_rerun.log`, `w2l_ground_protocol.log`. All copied
  to the laptop at `tmp\li_fetch\` (fetched via ssh `type`, byte-faithful).
- Changed files (hash-identical on both machines, uncommitted like the rest
  of the campaign): `test_w2l_ground.py` (protocol conversion),
  `build_w2l_split_net.py` + `build_w2l_aff_net.py` (census fixes),
  `tmp\li_fetch\probe_aff_census.py` (census probe),
  `tmp\li_fetch\run_li_sweep.ps1` / `run_w2l_gate.ps1` / `run_ground_only.ps1`
  (launchers).
- Ops note: `scp -O` (legacy SCP) is the working laptop->easteregg2 push;
  plain scp/SFTP rejects the spaced username, and `copy con`/PowerShell
  stdin pipes hang remotely. Remote->laptop pulls: ssh `type` with
  single-quoted paths (Git Bash mangles `\\$var` in double quotes).
