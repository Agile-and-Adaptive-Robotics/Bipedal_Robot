# Session 1: task-5 supported-control checkpoint

**Superseding geometry audit:** the recorded ankles reverse left/right order
in 91.8% of post-ramp samples, with approximately 28°/32° median toe-out.
Leg self-collision is disabled. The replay is a failed crossed-leg diagnostic;
its numerical equivalence and support-force measurements remain valid, but
it is not an acceptable physical stance/stepping result. The current ideal
tension actuators also do not implement the optimized BPA hardware.

Machine: **easteregg2**. Repository: `D:\GitHub\Bipedal_Robot`.
Python: `D:\Anaconda\envs\myo\python.exe`. MuJoCo 2.3.7; 2 ms step.

The 75 kg CAD robot now completes 40 s with actual vertical rig force capped
at 65% body weight after 15 s. The final 25 s include intermittent bilateral
foot lifting, but physical phase coverage fails. Regular supported stepping
and forward walking remain unverified. No five-stage curriculum has been run.

## Current construction and controller

- Coordinates: x anterior, y superior, z right; gravity acts downward along y.
- Lumbar ball joint, rigid 0.25 m spine, payload split equally between torso
  and ilium; total mass remains 75 kg. Welded feet use the supplied foot CAD.
- Robotbody right-side routes mirrored left; 27 actuators per side. Moving and
  conditional points update at runtime with their complete force Jacobians.
  Construction and dynamic-route verification pass; zero reversed neutral actions.
- Existing layered SNS architecture, rebuilt as 262 neurons and 230 inputs for
  54 actuators. RG events drive timing. An engineering posture/stance-load
  allocator supplies compensated MN target currents. This is hybrid control.
- The checkpoint uses a 5% allocated activation floor. Recorded activation
  never falls below 0.05 beyond floating-point error in the assessment window.
  The general runner default remains zero for explicit comparison. An activation
  floor is distinct from a measured coactivation index or intrinsic PAM stiffness.
- Forces use reference MIF caps and a 30 ms activation filter, not measured
  pneumatic force/pressure dynamics. Twenty ankle/toe actuators have no independent
  joint leverage because the feet are welded.

## Verified checkpoint

Evidence lives in `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso`.

| Check | Result |
|---|---|
| Fresh reference and exact-fast replay | Both complete 40 s |
| Complete logged trajectory comparison | Bit-identical in every saved array |
| Final support-cap interval | 15–40 s; actual vertical force at 65% BW |
| Activation floor | 5% |
| Peak right/left mesh clearance | 55.2 / 21.8 mm |
| Clearance-qualified swing bouts | Six per side |
| Neural RG bursts in final interval | 17 per side |
| Net forward movement in final interval | +8.2 mm |
| Physical alternating-cycle coverage | Fails; many neural cycles are skipped |
| Maximum pelvis / lumbar rotation | 8.32° / 17.12° |
| Mean absolute horizontal rig force, x/z | 13.20 / 31.53 N |
| Mean absolute rig torque, x/y/z | 13.93 / 3.50 / 17.24 N·m |

The rig includes angular stabilization, horizontal centering, and a 30 mm
RG-timed lateral pelvis target. The 65% number describes total actual vertical
rig force, not all forms of assistance. The raw model also permits excursions
past the CAD knee reference; physical stops and admissible ranges need explicit
verification before claiming hardware-equivalent gait. Poor hip-yaw tracking
causes the observed duck-footing. No trip perturbation was applied, so the
replay does not establish trip recovery.

Primary files: `supported65_backend_verification.json`,
`runs\verified65_reference.json/.npz`, `runs\verified65_exact_fast.json/.npz`,
and `runs\verified65_exact_fast_assessment.json`. Each verified run saves its
options, source hashes, and an exact runner snapshot captured at run start.
Exploratory trials precede this provenance improvement; the final fresh-process
verification is the reproducible checkpoint.

## What changed during support reduction

The original posture test needed about 95% vertical support. Lowering the rig
force alone caused collapse. Explicit ground-load compensation through the
muscle allocator allowed lower support. Separating knee lift from hip motion,
unloading the RG swing leg, and adding disclosed lateral support improved
clearance and survival. In one matched 75% support comparison, the zero-floor
case failed lumbar posture at 19.408 s while the 5% case completed 30 s.
Other comparisons showed that 5% can reduce clearance; its benefit is conditional.

Higher explicit gains, activation-tracking gain 3, and initial implicit-posture
experiments failed and are not promoted. Slower intrinsic RG time constants
were also checked: 0.7 s completed 40 s but failed physical phase coverage;
0.525 and 1.05 s lost lumbar posture. These are bounded diagnostics rather than
a completed curriculum or a successful optimization.

## View or reproduce on easteregg2

Recorded post-ramp animation: `robot_support65.gif`. Evidence plot:
`verified65_exact_fast_evidence.png/.pdf`. For the interactive recorded replay:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\replay_neural_robot.py' --run 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\runs\verified65_exact_fast.npz' --view
```

For a new simulation and assessment, with timestamped output:

```powershell
& 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\run_support65.ps1'
```

Next work is curriculum-based yaw/heading control, admissible knee motion,
consistent swing clearance and alternation, posture, and actual support-force
gates. Fresh perturbation basins remain necessary before promoting a gait
winner. No commits, pushes, or dissertation source changes were made here.
Only Session 1's bounded block in the shared `CHATGPT_REPORT.md` was updated;
bytes outside that block were verified unchanged.
