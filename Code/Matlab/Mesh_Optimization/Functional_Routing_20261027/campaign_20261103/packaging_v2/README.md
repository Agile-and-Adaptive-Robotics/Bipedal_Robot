# packaging_v2 — physical packaging deliverables (tasks 1–5)

Built 2026-11-04 on easteregg2 for the Functional Routing campaign
(directive items 3/5/6/7 as corrected by HARDWARE_AND_ZCODE_HANDOFF.md).
Nothing outside this folder and `tmp/functional_routing_20261103/` was
modified; `run_0` was not touched (read-only; the driver owns it).

## Contents

| file | task | what it is |
|---|---|---|
| `sdf/sdf_<body>.mat` | 1 | per-body signed-distance grids from the ACTUAL robot CAD meshes (`Code/MuJoCo_SNS/robot_torso` split meshes; build_ah001.py registration): pelvis, thigh_r, back_r, front_r, shank_r, foot_r. Grid in the piece's own frame; the frame-to-pelvis transform chain + provenance + validation stats are INSIDE each .mat |
| `sdf/cloud_<body>.mat` | 1 | cloud-only bodies (no robot CAD part): Spine, Pelvis_R, Sacrum, Talus, Calcaneus, Toes — Bone_Mesh_Plots point sets, conservative nearest-point distance, uncertainty bound inside |
| `sdf/knee_table.mat` | 1 | AH001 knee table of record: κ → (q_back, q_shank, q_front) per the model's own knee coordinate κ = −(q_back+q_shank) (`robot_kinematics.py knee_pose`); KB-pin track, mesh-PCA flexion reconstruction, per-κ knee-station residual |
| `registration_report.json` | 1 | registration errors (mm) + uncertainty bounds + method (JSON of record) |
| `REGISTRATION.md` | 1 | human-readable registration report |
| `attachment_patches.json` | 2 | per-group bounded anatomical surface patches (center + radius mm + patch frame + mount standoff) over all 27 contract groups, extracted from the `fr_cases27` builders; z=0 centerline rule citations inside |
| `fr_obstacles_v2.m` (campaign root) | 1 | loader + query: ground-frame points at pose q → per-body distance + nearest-body id; chains composed per pose with the recorded constants + knee table |
| `fr_camVasInt.m` (campaign root) | 3 | vas_int cam-assisted transmission candidate (polar cam on a femoral pivot + tibial patellar-tendon attachment patch); `README_camVasInt.md` here documents variables/bounds/units |
| `fr_kneeAssembly.m` (campaign root) | 4 | knee-crossing package model: pinned proximal + distal BPA sections, central Bowden/tendon span, brackets/ferrules/stops/exits, bend allowance, articulation clearance over knee 0..−120° |
| `fr_hardwareLoads.m` (campaign root) | 5 | rated-load check: T, R = 2·T·sin(β/2) vs the S/M/L catalog capacities (S 400 N / M 667 N / L 3003 N) with the UNRESOLVED-convention caveat attached to every use; overload = FLAGGED, never resized |
| `logs/check_*.log` | all | verification logs of record (MATLAB −batch, serial) |

## Hardware status (read before using any load number)

The S/M/L families are CATALOG-DIMENSIONAL CANDIDATES ONLY
(`hardware_candidates.json`, retrieved 2026-11-03 Pacific). Their
capacities are published weight capacities whose dynamic radial
interpretation is UNRESOLVED; no assembly WLL exists. **Final hardware
approval is NOT claimed by any file in this folder** — `fr_hardwareLoads`
returns `approved=false` with the caveat attached to every result.

Bend radius 75 mm is a BICYCLE PREVIEW ASSUMPTION, not a vendor rating;
ferrule length 12 mm is the Blender library PLACEHOLDER.

## Verification logs

- `logs/check_obstacles_v2.log` — SDF/chain/registration checks
- `logs/check_camVasInt.log` — cam candidate: per-pose geometry, arm
  curve, closure, 14/14 per-variable sensitivity
- `logs/check_kneeAssembly.log` — per-pose sweep tables + FLAG demo
- `logs/check_hardwareLoads.log` — point checks + built-mech walk

All checks run serially via `"D:/Program Files/MATLAB/R2025a/bin/
matlab.exe" -batch` with NO parpool (the run_0 driver owns the workers).
