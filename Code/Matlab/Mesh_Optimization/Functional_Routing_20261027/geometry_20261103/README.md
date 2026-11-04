# geometry_20261103 — Functional Routing campaign, directive item 1 (geometry before optimization)

Generated 2026-11-03 on easteregg2, OpenSim 4.4.1 Python (`D:/Anaconda/envs/myo/python.exe`),
from:

- `Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_robotbody.osim` (robot seeds)
- `Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim` (stock human reference)
- `Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json` (27-group map, seed preferences, human groups)
- `Code/Matlab/Bone_Mesh_Plots/Open_Sim_Bone_Geometry/*_Mesh_Points.txt` (bone surface clouds)

Generator scripts (scratch, reusable): `D:/GitHub/Bipedal_Robot/tmp/functional_routing_20261103/geometry/`
— `survey_models.py`, `extract_geometry.py`, `analyze_attention_bones.py`
(+ `probe_api*.py` API checks). Rerun: `D:/Anaconda/envs/myo/python.exe extract_geometry.py`
then `analyze_attention_bones.py`. No model, mesh, or existing code file was modified.

## Files in this folder

| file | content |
|---|---|
| `extraction_robotbody.json` | per-route exact XML declaration (owning bodies, exact location strings, conditional ranges in rad+deg, moving-point function arrays), per-point frame chains to ground, per-pose evaluated points vs OpenSim's own computed path (body-frame + ground-frame errors), body-in-ground transforms for the full right-side chain at all six poses, joint/coordinate defaults, winners of record |
| `extraction_simbody.json` | same, for the stock human reference routes (human_group targets) |
| `extraction_meta.json` | map-derived route requests + full census of every Conditional/MovingPathPoint in both models |
| `comparison_maxErrorMm.csv` | 462 rows (77 routes × 6 poses): model, actuator ids, role, route, seed status, pose, n_declared, n_matched_active, maxErrorMm, maxGroundErrorMm, membershipMatch, flags |
| `bone_distance_layering.csv` | 126 rows: group, attachment, role (origin/via/insertion), owning body, neutral bone-surface distance (mm) + which cloud |
| `bone_distance_layering_detail.json` | method, cloud point counts, nearest cloud point per attachment, erector-spinae seed-of-record distances |
| `attention_items.json` | quantified vasti / TFL / foot / erector-spinae findings + pelvis-frame z=0 audit |

## Verification protocol (what "maxErrorMm" means)

For each route and pose (neutral; hip flexion +30°; hip flexion −15°; knee 0°/−60°/−120°;
all other coordinates at model defaults):

1. The extracted representation is evaluated per pose: static points at their declared
   owning-body coordinates; ConditionalPathPoints by declared-range membership at the
   realized coordinate value; MovingPathPoints by evaluating the model's own declared
   location functions (`get_x/y/z_location().calcValue`) at the realized coordinate value.
   **Nothing is frozen at the default pose** — moving-point locations demonstrably track
   the pose (e.g. `vas_med_r-P5` x: 0.0560 m @ 0° → 0.0308 m @ −60° → 0.0098 m @ −120°,
   in `tibia_r`), and conditional membership demonstrably switches (vas_med_r P3/P4
   inactive at 0°/−60°, active at −120°, matching `isActive`).
2. OpenSim's own computed path is read via `GeometryPath.getCurrentPath(state)`
   (per-point stateful `getLocation(state)`, each point's owning body frame, plus
   membership) and a ground-frame cross-check is composed from
   `Body.getTransformInGround(state)`.
3. `maxErrorMm` = max over the active (optimizer hand-off) points of
   ‖my evaluation − OpenSim evaluation‖ in the point's owning body frame, ×1000.
   **Result: 0.0 for all 462 route–pose rows** (bit-identical: the same evaluator on the
   same declared inputs). `maxGroundErrorMm` is likewise 0. `membershipMatch` true
   everywhere; no flags.

This is exact per-point coordinate agreement, not bounding-box agreement.

Pose-setting respected model ranges and clamping: no clamping occurred except at the
knee −120° request, which sits on the model range endpoint (`knee_angle_r` range
[−2.09439510, +0.174533] rad); realized −119.99999986° (1.4e-7° range-precision
rounding, recorded per pose in the JSONs). No coordinate couplers exist in either
model (no `CoordinateCouplerConstraint`); all joints are `CustomJoint`s.

## Wrapping

Neither model contains any wrap object: 13 `WrapObjectSet`s per model, 0 objects;
92 `PathWrapSet`s per model, 0 `PathWrap`s. All "wrapping" behaviour in gait2392 is
realized by ConditionalPathPoints with knee-coordinate ranges (captured per point in
the JSONs, e.g. vas_med_r-P3 [−150°, −69.3°], med_gas_r stock P2 [−45°, +10°]).
The ask's "path wrap objects (regions and active ranges)" therefore resolves to an
empty set in both models; the functional equivalent with its active ranges is the
conditional-point census.

## Conditional / moving point census (both models, all 92 muscles)

- simbody: 26 conditional + 8 moving. robotbody: 18 conditional + 8 moving.
- All 8 moving points per model are the quadriceps tibial insertions
  (`rect_fem`, `vas_med`, `vas_int`, `vas_lat`, right and left), SimmSpline × 3
  (x,y,z) driven by the knee coordinate, owned by `tibia_r`/`tibia_l`.
- **Robotbody repair artifacts (Connor-era), preserved exactly and flagged:** on the
  right side, stock ConditionalPathPoints were replaced by plain `PathPoint`s
  literally named `"default"` (default-pose coordinates), or dropped outright:
  `rect_fem_r` stock 3 pts (P2 conditional on femur [−150°,−83.7°]) → robotbody 5 pts
  (P2, P3, 'default' plain, P3_0 moving); `med_gas_r`/`lat_gas_r` stock 3 pts (P2
  conditional on femur [−45°,+10°]) → robotbody 3 pts (P1, 'default', P3); `vas_int_r`
  stock 4 pts (P3 conditional on femur [−150°,−81.4°]) → robotbody 4 pts (P1, P2,
  'default', P4 moving); `semimem_r` stock 3 pts (P2 conditional on **tibia**
  [−32°,+10°]) → robotbody 5 pts (P1 pelvis, 'default' tibia, P3/P4/P5 tibia — the
  conditional frozen plus two extra tibia points); `grac_r` stock 4 pts (P2 conditional
  on tibia [−25°,+10°]) → robotbody 6 pts (adds P2/P2_0 on femur, conditional →
  'default' tibia); `iliacus_r` and `psoas_r` stock 5 pts (P3 conditional on pelvis,
  hip range [−90°,+45°]) → robotbody 3 pts with the conditional dropped entirely
  (psoas_r also starts on torso instead of pelvis). The stock model keeps the
  conditional forms with their ranges. The optimizer must implement these behaviours
  (or explicitly replace them with optimized mechanical guides) — the robotbody
  "default"-frozen forms must not silently stand in for them.

## Attention items (see `attention_items.json` for exact numbers)

1. **Vastus intermedius / combined vasti (actuator 11).** Seed preference
   `[vas_med_r, vas_int_r, vas_lat_r]`; winner `vas_med_r` (first well-formed entry;
   the map notes it carries the patella-tracking route). The "insertion" in both
   models is a **MovingPathPoint in the tibia frame** (Delp-1990 quadriceps
   mechanism; gait2392 has no patella body), identical in robotbody and stock
   (neutral x = 0.0560/0.0555/0.0600 m for med/int/lat; x → 0.0098 m at −120°).
   Deep-flexion routing is carried by femur-frame conditional points activating
   between −69.3° and −150°. The previous render's "too far anterior" insertion is
   contradicted by the model data; the optimizer hand-off uses exactly these
   function-evaluated tibia-frame values, and the actuator-map seed preference is
   confirmed (no arbitrary anterior substitute).
2. **TFL (actuator 7).** Robotbody `tfl_r` has 5 points (stock: 4). Robotbody P2
   `(0.037, −0.037, 0.071)` m sits 63.97 mm from stock P2 `(0.0294, −0.0995, 0.0597)`
   in the `femur_r` frame and 11.3 mm off at P3; lateral z 71.0/47.0 mm vs stock
   59.7/35.7 mm. At neutral the femur frame is aligned with ground (+z lateral for
   the right leg), so the robotbody route indeed runs higher and further lateral on
   the femur before the femur→tibia frame change at P3→P4 — the quantitative basis
   for the previous render's "too far lateral" artefact. Both routes are extracted
   exactly; the seed of record remains the map's `tfl_r` (robotbody), with the stock
   route available as anatomical bounds context.
3. **Foot routes (CalcnA / calcnSub / toes).** Master's-joint mapping: CalcnA →
   `ankle_r`/`ankle_angle_r` (tibia→talus), calcnSub → `subtalar_r`/`subtalar_angle_r`
   (talus→calcn), toes → `mtp_r`/`mtp_angle_r` (calcn→toes). All foot routes are
   plain PathPoint chains `tibia_r×2 → calcn_r×2–3 → toes_r×2–3`; no attachment sits
   on `talus_r`, so the tibia→calcn segments span **two** joints (ankle + subtalar)
   and their geometry depends on both transforms. Full transform chains and per-pose
   ground transforms for tibia/talus/calcn/toes are in `extraction_robotbody.json`.
   Verified 0.0 mm at all six poses.
4. **Erector spinae (actuator 1).** Stock `ercspn_r` = 2 points (P1 pelvis
   (−0.14, 0.0439, 0.0436) m; P2 torso (−0.055, 0.11, 0.0241) m). Robotbody
   `ercspn_r` = 9-point zigzag loop (torso y alternating 0.171/0.178/0.027/0.027/
   0.175/0.182/0.048 m) — the legacy manual wrap that predates the formal
   transmission scheme; its pelvis P1 also differs from stock P1 (16.05 mm).
   **Seed selected: the stock 2-point route** (original OpenSim coordinates,
   expressed in the robot pelvis/torso frames; identical joint definitions verified),
   per the directive's "initialize from the original OpenSim coordinates; do not
   preserve the manually wrapped loop". Loop coordinates remain documented in
   `extraction_robotbody.json`. Seed-route bone distances: P1 5.36 mm (sacrum cloud),
   P2 18.63 mm (spine cloud).

## Bone-distance layering table

`bone_distance_layering.csv` (126 attachments over the 27 winner routes): Euclidean
distance from each neutral attachment (declared location; moving points at their
neutral evaluated location, marked `evaluated_at_neutral(moving_point)` in
`location_source`) to the nearest point of the assigned cloud — Femur/Tibia/Calcaneus/
Talus/Toes/Pelvis_R+Sacrum (nearest of the two)/Spine. These are **ordering metrics**,
not clearance tolerances: the clouds are sparse surface samples (femur 456 pts, tibia
144, pelvis_R 124, talus 98, sacrum 191, calcaneus 999, toes 957, spine 5433), so
distances carry cloud granularity. Empirical alignment check (attachment NN-distance
means): calcn 3.7 mm, toes 2.1 mm, sacrum 7.5 mm, pelvis_R 18.6 mm, tibia 14.5 mm,
femur 26.8 mm, spine 35.8 mm — dense clouds agree tightly; femur/pelvis/spine values
are coarser. Right side only (campaign convention), neutral pose.

## Pelvis-frame centerline audit (handoff item 8 precondition)

All winner-route attachments transformed into the pelvis frame at neutral:
minimum z = **5.0 mm**, zero attachments on or across z=0. No seed needs centerline
disambiguation.

## API notes (opensim 4.4.1 bindings gotchas, for reuse)

- `GeometryPath.getCurrentPathPoints` does not exist; use `getCurrentPath(state)`
  (returns `ArrayPathPoint` with names/bodies/stateful locations).
- `isinstance(p, osim.MovingPathPoint)` is unreliable on `PathPointSet.get(i)`;
  use `osim.MovingPathPoint.safeDownCast(p)` (works), then `get_x_location()` /
  `get_y_location()` / `get_z_location()` (NOT `getXFunction`) and
  `calcValue(osim.Vector(1, 0.0))`.
- `ConditionalPathPoint.getRangeMin/getRangeMax` are not exposed; declared ranges
  were taken from the XML (authoritative source) and cross-checked against
  `isActive(state)` and computed-path membership.
- Vec3 iteration is broken in this binding; index explicitly `v[0] v[1] v[2]`.
- `Transform * Vec3` is not wrapped; compose `R@x + p` from `T.R().get(i,j)` /
  `T.p().get(i)`.

## Scope statement

This deliverable is geometry only (directive item 1): extraction, per-pose
verification, and the layering metric. No optimization, no Xi factors, no packaging —
those are later directive items and consume these files as inputs.
