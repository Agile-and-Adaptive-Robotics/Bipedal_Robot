# Registration report — physical packaging obstacle model (task 1)

JSON of record: `registration_report.json` (same folder). Built by
`tmp/functional_routing_20261103/packaging_v2_build/build_sdf.py`
(easteregg2, 2026-11-04).

## Registration chain

- CAD pieces: ah001_cad_meshes split meshes (build_ah001.py localize_mesh), vertices already in robot body frames; body frames ISB (x anterior, y superior, z right) at the pelvis datum (AGENTS.md 2026-10-27 ruling; build_ah001.py cad_to_isb normalizes through the pelvis part 02_01_PE_ datum)
- robot joint anchors = CAD anchors (build_ah001.py uses the audit joint anchor_a; MuJoCo body offsets = the anchor deltas)
- OpenSim neutral transforms parsed from gait2392_robotbody.osim and validated EXACTLY (max|dT| < 1e-9) against geometry_20261103/extraction_robotbody.json bodyTransformsInGround.neutral
- registration anchor: align the HIP centers (translate(hipOS - hipR)); pelvis frames both ISB so rotations are identity; knee/ankle mismatches are then measured, not absorbed
- robot knee state vs OpenSim knee_angle_r: four-bar solved for (theta_back, theta_front) at each sampled theta_shank; flexion mapping theta_shank(kappa) = sFlex*(kappa*-1 - flex0) verified by reconstructing the femur/shank axis angle

## Registration error table (mm)

| anchor | error (mm) | definition |
|---|---|---|
| hip | 26.66 | robot hip ball (CAD knob-pin anchor, MuJoCo `thigh_r` offset) vs OpenSim `hip_r` pelvis_offset; registration translation of record `treg = hipOS − hipR` |
| knee | 86.97 | robot KB pin (four-bar knee-bottom pivot = `shank_r` origin) carried on the hip-aligned femur segment vs OpenSim knee center of record (the TIBIA BODY ORIGIN at neutral — gait2392 knee: rotation axis at the femur origin, knee-angle splines carry the tibia). Dominated by the segment-split difference (robot femur link longer, shank shorter); per-segment anchoring keeps each SDF piece correct |
| knee station over sweep | mean 53.5, max 87.0 (53.5 @ −60°) | pose-dependent residual: robot four-bar KB pin vs OpenSim moving-pivot tibia origin at the same κ — the two knee KINEMATICS differ (splines vs four-bar). Inherent to CAD-mechanism + OpenSim-routing; quoted, not hidden |
| ankle station | 88.80 | robot `foot_r` origin = mesh-derived distal tibia end center (`foot_geometry.py distal_center`), welded-foot station vs OpenSim `ankle_r` location; the robot has NO ankle pin in this build (welded foot) — a station error, not a joint error |

## Cross-checks

- Knee flexion axis: OpenSim knee flexion axis vs robot knee pin axis misalignment (abs): **0.000 deg**.
- Shank length: OpenSim **430.0 mm** (knee→ankle) vs robot KB→foot **324.3 mm**.
- Four-bar: max closure residual **4.44e-16 m**; knee coordinate of record κ = −(q_back + q_shank) (`robot_kinematics.py knee_pose`, verified to 1.38e-9 m in dynamic_route_verification); chain composition cross-validated against MuJoCo `mj_kinematics` in the build (max KB-pin deviation 5.6e-17 m).
- Physical flexion at robot kneeAngle = 0: mesh-PCA **22.4 deg** (pin-based 27.0 deg; CORRECTED_SCENE quoted 17.20 deg by its own axis definition) — the robot model zero carries physical flexion; the obstacle table uses the model convention κ_robot = κ_OpenSim (δ = 0).
- OpenSim neutral composition validated against `geometry_20261103/extraction_robotbody.json` bodyTransformsInGround.neutral: max |dT| ≤ 2.1e-06 m.

### Sparse-cloud cross-checks (OpenSim joint centers vs bone clouds, owning body frames)

| check | distance (mm) |
|---|---|
| hipCenter to FemurCloud mm | 19.61 |
| kneeCenter to FemurCloud mm | 3.41 |
| kneeCenter to TibiaCloud mm | 30.31 |
| ankleCenter to TibiaCloud mm | 18.75 |
| ankleCenter to TalusCloud mm | 4.86 |

These distances include cloud granularity (sparse surface samples) — ordering context, not clearance tolerances.

## Cloud-only bodies and uncertainty bounds

| body | frame | points | NN mean (mm) | NN median (mm) | NN max (mm) |
|---|---|---|---|---|---|
| Spine | torso | 5434 | 3.2 | 3.1 | 11.83 |
| Pelvis_R | pelvis | 125 | 11.6 | 11.1 | 28.13 |
| Sacrum | pelvis | 192 | 6.6 | 6.6 | 14.55 |
| Talus | talus_r | 99 | 4.9 | 4.9 | 7.77 |
| Calcaneus | calcn_r | 1000 | 3.4 | 3.2 | 9.18 |
| Toes | toes_r | 958 | 1.9 | 1.7 | 5.91 |

Nearest-cloud distance OVERestimates the true surface distance by up to the local NN spacing → conservative for clearance. Cloud-only bodies are: Spine, Pelvis_R, Sacrum, Talus, Calcaneus, Toes (no robot CAD part in `ah001_cad_meshes`).

## SDF coverage per body

| body | triangles | grid | voxel (mm) | inside % | vertex |sdf| max (mm) | near-surface parity agreement |
|---|---|---|---|---|---|---|
| pelvis | 447300 | 81×91×110 | 3.0 | 9.5 | 3.764 | 0.592 |
| thigh_r | 193300 | 54×225×64 | 2.5 | 6.6 | 3.085 | 0.723 |
| back_r | 6928 | 47×40×37 | 2.5 | 10.6 | 1.921 | 0.906 |
| front_r | 6920 | 37×51×37 | 2.5 | 9.4 | 2.880 | 0.886 |
| shank_r | 61414 | 87×161×55 | 2.5 | 7.6 | 3.184 | 0.822 |
| foot_r | 108002 | 107×55×56 | 3.0 | 10.1 | 3.598 | 0.768 |

SDF convention: meters, negative inside, clamped ±30 mm; queries outside a grid return +clamp (conservative). Trilinear interpolation in `fr_obstacles_v2.m`.

## Usage rule

Attachment searches may use this model for clearance, but a large attachment excursion justified by uncertain registration is FORBIDDEN: the campaign patch bounds stay anatomically tight (surface patch radius + mount standoff in mm, fr_patch/fr_setup conventions).

## Limitations

- cloud-only bodies (Spine, Pelvis_R, Sacrum, Talus, Calcaneus, Toes): sparse-sample distance, conservative; NN bounds recorded per cloud
- robot hip is a BALL joint mapped by the OpenSim hip axes sequence (exact for the OpenSim joint model; the physical robot hip has the same axes at its zero pose, both ISB)
- the SDF sign uses the closest-triangle normal with a ray-parity consistency check on the near-surface shell (agreement recorded per body)
- pelvis/foot grids at 3 mm voxels, others 2.5 mm; SDF clamped at +/-30 mm; queries outside the grid return +clamp (conservative for clearance)

## Files

- `sdf/sdf_<body>.mat` — pelvis, thigh_r, back_r, front_r, shank_r, foot_r (grid + chain + provenance + stats)
- `sdf/cloud_<body>.mat` — Spine, Pelvis_R, Sacrum, Talus, Calcaneus, Toes (point sets + uncertainty)
- `sdf/knee_table.mat` — κ → (q_back, q_shank, q_front), KB-pin track, flexion reconstruction, per-κ knee-station residual
- loader/query: `campaign_20261103/fr_obstacles_v2.m`
