# Foot-mount sensitivity and w2lvar timing comparison

Session 1, easteregg2, 28 October 2026. Repository `D:\GitHub\Bipedal_Robot`.
Environment `D:\Anaconda\envs\myo\python.exe`. All four conditions completed
40 simulated seconds. Values below describe the 15–40 s interval, after
vertical assistance reaches 80% body weight. Lateral and angular assistance,
provisional knee stops, and all other engineered gains remain active.

| Measure | s3k hybrid, original mount | s3k hybrid, proximal 60 mm | w2lvar timing, original mount | w2lvar timing, proximal 60 mm |
|---|---:|---:|---:|---:|
| Median ankle separation (mm) | 159.0 | 57.3 | 63.8 | 99.0 |
| Minimum ankle separation (mm) | 142.4 | 43.2 | 53.6 | 55.2 |
| Reversed ankle ordering | 0% | 0% | 0% | 0% |
| Maximum sampled interleg penetration (mm) | 0 | 1.316 | 0 | 0 |
| Maximum sampled knee-loop position error (mm) | 0.387 | 0.361 | 0.516 | 0.784 |
| Mean right TFL activation | 81.8% | 17.2% | 37.9% | 22.5% |
| Mean left TFL activation | 99.8% | 72.5% | 99.9% | 95.1% |
| RG burst counts, R/L | 17/17 | 17/17 | 15/14 | 15/14 |
| Bilateral sustained unloading gate | Fail | Pass | Fail | Pass |
| Bilateral repeated lifting gate | Fail | Fail | Fail | Pass |
| Repeated alternating foot-apex gate | Fail | Fail | Fail | Fail |
| Forward progression gate | Fail | Fail | Fail | Fail |

The attachment change materially affects loading, separation, and foot motion.
No condition establishes regular stepping. The shifted w2lvar timing condition
is a useful follow-up candidate because it passes bilateral unloading and
lifting counts, but only 25% of right-foot apex intervals contain a qualifying
left-foot apex. Neural antiphase alone does not establish alternating foot motion.

## What was changed

The baseline foot is an entire welded CAD assembly placed at the distal shank,
with a top-of-assembly datum approximately 103.7 mm above the sole. No verified
mechanical mating datum is available. The hypothesis moves each whole foot
60 mm toward the proximal origin of its shank body, along the original
origin-to-mount vector. This also moves foot-attached route sites. It does not
shorten tibial meshes, remove material, change tibial inertia, or redesign the
interface. The full foot's local mass and inertia are retained; total model mass
remains 75 kg. Identical initial joint angles are used, with soles placed at
ground again. The original model is preserved.

## What the w2lvar comparison means

The variant is built from `spinal\build_network_w2lvar.py` with the actual 54
actuators: **544 neurons and 230 inputs**. Four additional named ankle indices
alias the knee PF cells; the inherited summary's 548 index count must not be
quoted as a neuron census. There are two antagonistic PF pairs per leg,
HIP E/F and shared KNEE/ANKLE E/F, each with inhibitory interneurons. This is
the existing w2lvar mapping, not a new independently optimized functional map.

The inherited engineering MN target-current clamp remains active and cancels
incoming PF/reflex MN synaptic currents. These runs therefore compare RG
timing and afferent-central interactions within the same postural force
allocator. They do **not** compare full neural PF-generated muscle patterns,
and they do not show that w2lvar is superior as a locomotor controller.
The next neural-control experiment must retain PF/reflex synaptic current and
bound the engineering postural contribution separately, then use the staged
air/ground curriculum. All actuators here remain ideal reference tensions,
not calibrated BPA pressure/strain/compliance actuators.

The exact-fast backend is used for every MuJoCo run. The baseline has a fresh
full-trajectory reference replay. The 544-neuron variant has zero voltage
difference from the reference backend across 400 seeded-input steps; no
full-trajectory replay, parameter sensitivity, or robustness claim is made
for these exploratory variant runs.

## Virtual attachments

Route points belong to prescribed model bodies. A point can float outside a
surface and still transmit tendon force through that body's Jacobian. This
assumes an ideal rigid bracket or guide at that location. No nearest-solid
snapping, bracket mass, bracket deflection, physical mounting design, or
collision-aware wrapping is present. Moving and conditional path functions
remain active, with their virtual-work derivatives.

## Reproduce on easteregg2

Run from `D:\GitHub\Bipedal_Robot`:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\probe_foot_mount_w2lvar.py' --prepare
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\probe_foot_mount_w2lvar.py' --network s3k --shifted
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\probe_foot_mount_w2lvar.py' --network w2lvar
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\probe_foot_mount_w2lvar.py' --network w2lvar --shifted
```

These deterministic probe names overwrite their own exploratory outputs.
Evidence is in `runs\abd_defense_exact.*`, `runs\s3k_foot60.*`,
`runs\w2lvar_foot0.*`, and `runs\w2lvar_foot60.*`, including the associated
abduction, geometry, and stepping assessment JSONs. Source-model hashes and
the translation vectors are in `foot_mount_sensitivity_spec.json`; variant
gains, source hash, true census, and adapter hash are saved in each variant
summary. `w2lvar_probe_backend_verification.json` records backend parity.
