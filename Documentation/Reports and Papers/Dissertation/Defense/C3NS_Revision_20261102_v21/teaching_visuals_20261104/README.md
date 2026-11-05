# Defense teaching visuals — 2026-11-04

**Superseded rendering instructions:** use [REPAIRS_20261105.md](REPAIRS_20261105.md)
and `build_fixed_20261105/` for the current corrected visuals. The cam now uses
Ben's authoritative `09_BA_003.SLDASM`, its native component placements, and
its named Theta1 datum. The older inferred cam registration described below
is historical and must not be used to rebuild the current deck.

Three scientific teaching animations for the defense deck. Nothing in the
dissertation, controllers, optimization winners, or experimental data was
modified: every rendered quantity is read from the sources of record listed
below. Rendered outputs (GIFs, stills, frames) are under `build/`
(git-ignored, stays on this machine — the defense computer embeds the GIF
files directly, no external links). Authoring sources (`scripts/`,
`data/`, `scenes/`) are meant to be committed.

## Deliverables

| Visual | Looping GIF (embedded in deck) | Labeled still | Unlabeled still | Source scene |
|---|---|---|---|---|
| 1. Cam-assisted extensor | `build/visual1/cam_extensor_loop.gif` (8 s) | `build/visual1/still_labeled.png` (1080 px) | `build/visual1/still_unlabeled.png` | `scenes/visual1_cam_extensor.blend` |
| 2. Reverse pulley | `build/visual2/reverse_pulley_loop.gif` (7 s) | `build/visual2/still_labeled.png` (1080 px) | `build/visual2/still_unlabeled.png` | `scenes/visual2_reverse_pulley.blend` |
| 3. Knee coordinate frames | `build/visual3/knee_frames_loop.gif` (9 s) | `build/visual3/still_side.png` + `still_oblique.png` (1080 px) | — | `scenes/visual3_knee_frames.blend` |

All: white background, right-side view (anatomical +z camera; visual 1 adds a
short oblique intro/outro), Arial Bold labels, no raw coordinates or
model-code names on screen.

## Rebuild

```cmd
:: geometry extraction of record (MATLAB R2025b; read-only on the mat)
matlab -batch "run('.../scripts/dump_cam_visuals_20261104.m')"
:: derived per-pose geometry + validation
python scripts/prep_geometry.py
python scripts/check_cam_geom_20261104.py        :: moment-arm replica check
python scripts/verify_gait2392_fk_20261104.py    :: Gait2392 FK vs MuJoCo (myoconv env)
:: scenes + stills + frames + GIFs
blender --background --factory-startup --python scripts/build_visual1_cam.py -- [--stills-only]
blender --background --factory-startup --python scripts/build_visual2_pulley.py -- [--stills-only]
blender --background --factory-startup --python scripts/build_visual3_frames.py -- [--stills-only]
python scripts/make_gifs.py
```

Blender: portable 5.2.1 at
`C:\Users\Ben\.zcode\cli\plugins\data\blender@zcode-plugins-official\engine\blender-5.2.1-windows-x64\blender.exe`
(engine = EEVEE, 32 samples, 720 px frames / 1080 px stills, 15 fps).

---

## Visual 1 — cam-assisted extensor

**Source configuration.** Adopted cam design of record
`Code/Matlab/Mesh_Optimization/Results/Vas_Pam_20mm_Result_cam_20261014_0416.mat`
(`ctxCamBest` + `predBest`, minted by `save_cam_result_20261014.m`), the same
mat behind the dissertation figures `extensor_cam_AF_bones.pdf`,
`extensor_cam_GH_bones.pdf`, `cam_pose_A.pdf`, `cam_pose_H.pdf` (generator
`Figures/94-AppendixC/Components/registration_20261028/diss_cam_c_registered.m`,
ported from `Export_AppendixRouteGeo.m` / `Opt_run_Ext.m`). Rendered per pose
(100 solver poses, −120°…+10°):

- **Route p1…p9** — `predBest.routeInfo.raw`, rows 1–5 native femur frame,
  rows 6–9 t1 frame mapped to femur exactly as the figure generator does
  (`RowVecTrans(ctx.T_ICR_t1)`, then `RowVecTrans(ctx.T_Pam)`). Drawn as one
  blue tube: the dissertation figure of record draws the whole optimized
  route as the blue "BPA path", and the sources do not delineate a BPA/cable
  transition point, so the label reads "BPA / cable".
- **Shaped cam** — the actual `TibiaClip_002` hardware (both halves, machined
  steel tone), registered to the t1 frame as described in the missing-assets
  section below (wrap circle + p8 clip-edge rule). The wrap circle equals the
  circumcircle of the full-flexion contacts p6/p7/p8
  (`geo.seedP6/P7/P8XY`; circumradius 59.8 mm; center t1 (71.8, −24.9) mm),
  and the black contact markers sit at those adopted seeds on the real
  surface. The **96.9 mm** value of record is the *loss-model drum radius*
  `routeInfo.wrap.radius(2,:)` under Ben's 2026-10-14 p7x−p9x rule — the
  published figure shows it as an annotation, not drawn geometry.
- **Knee pivot** — the migrating instant center, `T_Pam(1:3,4)` per pose
  (robot tables `knee_x_Pam/knee_y_Pam`): migrates 4.1→23.4 mm in x,
  −417.4→−399.1 mm in y across the sweep.
- **Moment arm** — perpendicular from the pivot to the line of action through
  the crossing point p6 along the incoming femoral segment: this is exactly
  the evaluator's `predBest.momentArm` construction
  (`MonoPamDataExplicit_balanceX3.Mom`, cross-point 6). Validated: replica vs
  evaluator mean difference 0.18 mm (max 1.06 mm) over the sweep; range of
  record 59.4–102.8 mm. Drawn gold with a force arrow along the same line of
  action (tension pulling the tibia toward the femoral side).
- **Femoral pulley** — `geo.femurCylCenter` (53.9, −274.8) mm,
  radius `geo.femurCylRadius` = 10 mm (the femoral cylinder the route wraps
  at p2; p2 releases at −21.5°, the contact-release schedule comes from
  `routeInfo.active`/`eliminatedAngleD`).
- **Bones** — OpenSim `femur_r.stl`/`tibia_r.stl` from the MyoConverter
  conversion (`Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/Geometry/`,
  vertex bounds identical to `Bone_Mesh_Plots/Open_Sim_Bone_Geometry` clouds —
  verified), registered exactly like the dissertation figures
  (`appendix_c_bone_background.m` offsets: femur condyle
  [6.52, −408.00] mm → `geo.femurProfileCenter`, tibial groove
  [6.00, −31.255] mm → `geo.tibiaLowerCenter`). Bones fade to 22% opacity
  during the oblique intro so the pivot/cam/force path are readable.

**Animation.** 8 s, seamless loop, 15 fps: oblique overview of the physical
assembly with opaque bones (1.2 s) → bones fade while the camera settles into
the side view (1.2 s) → slow flexion→extension sweep −120°→+10° (2.8 s) →
extension→flexion (2.8 s) → camera returns to oblique and bones unfade (1.2 s).
The labeled still is at −30° (all contacts active).

**Teaching point.** As the knee extends, the pivot-to-line-of-action
perpendicular (moment arm) grows 59→103 mm: the same cable tension produces
more extensor torque (τ = F·r), so the same torque needs less pulling force
when the cable runs farther from the pivot.

## Visual 2 — reverse tackle, two parallel BPAs

**Source configuration.** The exact schematic geometry of
`Figures/30-torque/reverse_pulley_principle.pdf` and
`reverse_pulley_biomimetic.pdf` — generators
`Figures/30-torque/Components/reverse_pulley_principle/draw_reverse_pulley_publication.py`
and `draw_biomimetic_reverse_pulley_20261029.py` (femur/tibia polygons,
crossed-linkage knee, BPA spans y 5.55/6.45 x 1.15–3.8, moving sheave
(4.8, 6.0) r 0.45, fixed tackle mount x 7.35, fixed redirect (7.35, 5.1),
tendon tangent x 7.8, insertion (7.8, 2.35) — all coordinates copied from
those scripts). Ideal kinematics only: each BPA contracts by Δℓ, the moving
sheave translates Δℓ, and the output tendon past the fixed redirect feeds out
2Δℓ (two line parts between moving sheave and fixed mount). Labels: fixed
pulley, moving pulley, output tendon, plus the ideal force relation
"output force = half the combined BPA force" (n = 2 BPAs, G = 2 →
F_t/F_b = n/G = 1 per actuator = half of the combined 2F_b).

**Kept separate by design (as instructed):** measured efficiency and any
dynamic BPA validation result are NOT shown here.

**Animation.** 7 s seamless loop: BPAs contract (arrows track), moving
pulley translates, output tendon + insertion diamond feed out twice the
distance, then reverse. The insertion diamond rides the tendon end to
visualize the travel relation; in the hardware the tibia rotates under this
pull (the schematic source figure states "Geometry is schematic, not to
scale").

## Visual 3 — Gait2392 knee coordinate frames

**Source configuration.** Stock OpenSim Gait2392 (`gait2392_simbody.osim`,
"3DGaitModel2392"): right femur and tibia bodies with all muscles hidden.
Bone meshes are the conversion's `femur_r.stl`/`tibia_r.stl` in OpenSim body
frames. The knee is the model's `knee_r` CustomJoint — parent frame
`femur_r_offset` (= femur_r body frame, origin at the hip), child frame
`tibia_r_offset` (= tibia_r body frame, origin at the knee), single
coordinate `knee_angle_r` (range −2.0944…+0.1745 rad, flexion NEGATIVE):

- tibia origin in the femur frame: (tx(q), ty(q), 0) from the model's
  SimmSpline tables (natural-cubic-spline evaluation, same tables the robot
  context builder uses as its human reference);
- tibia orientation: Rz(q) about that displaced origin (migrating
  tibiofemoral contact point, not a pure hinge).

This composition was verified against the repo's converted MuJoCo model
(`mjc/gait2392_simbody/gait2392_simbody_cvt3.xml`, MyoConverter): rotation
error 0.000°, position ≤ 0.87 mm over the range (quartic-fit difference) —
`scripts/verify_gait2392_fk_20261104.py`.

**Rendered.** Parent (femur) and child (tibia) frame triads — RGB arrows
x/y/z with letters; the child triad is rigid on the tibia and shows the
orientation change through flexion. Labels "Femur (parent frame)" /
"Tibia (child frame)", live readout "Knee angle: −N deg".

**Animation.** 9 s seamless loop, 15 fps: prescribed flexion 0°→−120°→0°
(eased turnarounds), side view; oblique companion still at −92°. This is a
kinematic reference illustration only — no loads, no control.

## Missing / approximated source assets (nothing silently invented)

1. **Shaped-cam hardware: REAL, updated 2026-11-04.** An earlier revision of
   this README wrongly said no shaped-cam CAD existed (the `CAM/` folder holds
   CNC programs for the rectangular plate, and the first sweep of
   `Knee assembly/` missed it). Ben pointed to
   `Knee assembly/TibiaClip_002.SLDPRT`; its same-day STL exports
   (`STL Files/TibiaClip_002_lateral.STL` / `_medial.STL`) are rendered as
   the shaped cam in visual 1. Registration to the evaluator t1 frame:
   (a) wrap circle — fixed-radius Hough (59.83 mm, equal to the evaluator
   circumcircle of the adopted contacts) places the part-frame wrap center at
   (119.1, 52.9) mm, mapped to the evaluator center (71.79, −24.92) mm;
   (b) rotation — the evaluator documents p8 as the "clip-edge" contact, so
   the registered arc ENDPOINT lands on p8's angular position (Δ = +145.95°),
   cross-validated by p6/p7 landing strictly inside the arc surface;
   (c) z — the halves straddle z 0..92.2 mm with a 0.5 mm cable slot at the
   mid-plane, so scene z = part z − 46.1 mm and the cable runs in the slot.
   The rotation is a single-DOF determination from the p8-endpoint rule, not
   an assembly-mate readout; if Ben's CAD mounting differs, only that rotation
   changes and the wrap kinematics stay exact.
2. **Origin bracket / insertion bracket 3D models:** not extracted anywhere;
   rendered as point markers at the adopted design points p1 (femur frame)
   and pEnd/p9 (tibia frame). (`Tendon_Extensor.SLDPRT`, `Path_guide.SLDPRT`
   and the `Tibia.SLDPRT` bone CAD also exist in the assembly folder and were
   not needed for these visuals.)
3. **Display-only thicknesses:** the evaluator is a sagittal-plane model.
   Femoral pulley width 20 mm, cable radius 5 mm, marker sizes — display
   choices, no mechanical meaning. The clip halves render at their true CAD
   thickness.
4. **Bone vtp geometry** (`femur_r.vtp` etc.) is not installed on this
   machine; the conversion-extracted STLs (same vertex geometry, verified by
   identical bounds) are used instead.
5. **Clearance envelopes** (gold/pink circles in the dissertation figures)
   are analysis overlays, not hardware — deliberately omitted from the
   animation to keep the pivot/cam/force path readable.
6. **Visual 2 output travel** is illustrated as the output tendon feeding out
   past the fixed pulley with the insertion diamond riding its end (ideal
   tackle relation); the physical tibia rotation under load is out of scope
   for this schematic, matching the source figure's own framing.

## Frame-convention fix (2026-11-04, clip integration)

Integrating the real clip exposed a missing translation in the tibia object's
motion: the figure-of-record chain is
`femur = Rz(phi) * (p_t1 - t1toICR(phi)) + pivot(phi)`, and the builder had
dropped the `- t1toICR(phi)` term (a 24-30 mm pose-dependent offset on the
tibia bone, clip, and tibia-frame markers; the route/cable was always correct
since it lives in the femur frame). Fixed by keyframing the tibia object
location as `pivot - Rz(phi)*t1toICR(phi)`; contact markers, insertion,
clip, and the cam label use raw t1 coordinates as tibia children.

## Validation summary

| Check | Script | Result |
|---|---|---|
| Moment-arm replica vs evaluator of record | `check_cam_geom_20261104.py` | mean −0.18 mm, max 1.06 mm over 100 poses |
| Cam loss radius assertion | `prep_geometry.py` | 96.90 mm (p7x−p9x rule) |
| Gait2392 knee FK vs converted MuJoCo model | `verify_gait2392_fk_20261104.py` | 0.000° rotation, ≤ 0.87 mm position |
| Bone STL vs registered point clouds | `inspect_bone_stl_20261104.py` | identical vertex bounds |

Run outputs (`build/`) are git-ignored per the generated-defense-files rule;
Ben commits `scripts/`, `data/`, `scenes/`, and this README.
