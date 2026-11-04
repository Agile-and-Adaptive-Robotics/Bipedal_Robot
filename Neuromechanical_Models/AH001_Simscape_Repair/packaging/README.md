# AH001 packaging inspection

Current corrected Blender scene: `output/AH001_Packaging_Corrected.blend`. See [CORRECTED_SCENE.md](CORRECTED_SCENE.md) for the 20 mm counts, pinned/Bowden sections, foot registration, anatomical knee limits, white background, and camera orbit. The files described below are the earlier preview.

Open `output/AH001_Packaging.html` in Chrome or Edge. It runs offline. The selected KompareXXVII side has 27 routes and 69 individual BPAs. The opposite-side checkbox mirrors the routes for visual inspection; it is not a second independently optimized design.

## Files

- `output/AH001_Packaging.html`: orbit, pan, isolate a route, inspect its saved state, sweep both closed-loop CAD knees, show maximum diameter envelopes, save a PNG, or export the visible geometry as OBJ.
- `output/AH001_Packaging.glb`: CAD, spine, BPA sleeves, fittings, and tendon tails at the source-grid pose nearest neutral. Import into Blender with File > Import > glTF 2.0. Components and BPA segments remain separate objects. This is a static snapshot.
- `output/AH001_Packaging.blend`: native Blender 5.3 scene, now built and saved through the running Blender Python console. It includes 69 animated sleeves, 19 CAD poses, and a 217-frame forward/back sweep, with separate route collections, fittings, and tendons. The original imported scene is preserved as a separate scene in this file. `Build_Blender.cmd` is the command-line rebuild entry point; launching the desktop executable from the restricted shell still returns 0xC0000022, while the running GUI works.
- `output/Defense_Packaging.png` and `output/Defense_Erector_Spinae.png`: clean stills.
- `output/Defense_Knee_Sweep.gif`: the vastus-intermedius envelope over the CAD knee sweep. Caption it as a packaging preview.
- `output/bundle_overlap_candidates.csv`: conservative capsule gaps between individual BPA sleeves, including members of the same bundle. Negative values are overlap candidates under this preview placement.
- `output/route_packaging_screen.csv`: actuator dimensions, sleeve length, route length, and sharp-corner flags.
- `output/mesh_registration.json`, `verification.json`, and `browser_verification.json`: placement and UI checks.

## What changed

The inherited STL files use millimeters and independently translate each part into positive coordinate space. Scaling alone does not restore the assembly. The builder recovers each translation by matching STL vertices to the original STEP geometry. The foot STEP files use inches; the leg STEP files use millimeters. All output geometry uses meters.

The restored local part geometry is placed with the audited SOLIDWORKS component transforms, normalized through the pelvis datum. It does not use `cad_to_mujoco_fit.json`. The rigid groups and joint frames come from the existing AH001 audit and body tree. Both four-bar loop closures are solved at every CAD sweep pose. The maximum positional closure residual across the 19 poses is 7.66e-12 m. This verifies the joint-frame loop calculation, not assembly clearance or the survey attachment registration.

The new view uses recorded `sContraction` for BPA length: `L = RestingL - sContraction`. It does not treat whole-route length as actuator length. The expanding diameter follows the existing calibrated braid law. At the recorded calibration endpoints, 10/20/40 mm uninflated BPAs expand to 22.5/38/75 mm at 17/25/32% contraction. Values beyond these endpoints are flagged rather than represented as additional calibrated expansion. Saved fitting lengths and parallel-bundle offsets are retained.

## Scope of the packaging preview

1. Survey attachment registration remains provisional. The saved source transforms are recovered directly from the MAT files, folded at zero degrees, and attached to corresponding CAD bodies for this preview. This is not a measured reconciliation between the master's survey geometry and the newer AH001 CAD.
2. Each BPA sleeve occupies an arc-length interval of the saved polyline beginning after the origin fitting. Its remaining path is drawn as tendon. This is an explicit provisional installation rule, not a physical mount layout recovered from the optimizer. A sleeve crossing sharp waypoint corners is a redesign warning.
3. G is a travel multiplier. The saved designs do not specify sheave positions, diameters, bearings, brackets, or cable reeving. No invented sheaves are presented as saved hardware.
4. The knee sweep uses the CAD four-bar geometry. With Follow saved knee contraction enabled, actuator contraction is interpolated from the saved source knee grid, with the other source angle at its sample nearest zero. This does not re-solve equilibrium on the registered CAD geometry. Other saved-state sliders inspect actuator envelope extremes without moving all CAD joints.
5. Feet inherit the existing model's mounting assumption. The separate foot joints are not articulated in this view. The rigid 0.25 m spine is included; its approximate lumbar mount comes from the existing model.

The foot visualization was corrected after Ben's import inspection: assembly coordinates are mapped through the recorded distal mount and width/up/anterior-to-anatomical basis, matching `foot_geometry.py`. The prior visualization incorrectly treated the Ankle Exterior part datum as the foot-body datum and mirrored the wrong side. `AH001_Packaging_FeetFixed.glb` is the corrected static import; the native `.blend` also includes this correction. The model's underlying mount assumption remains unchanged.
6. Capsule screening covers BPA-to-BPA sleeve envelopes on the selected side at the saved grid sample nearest neutral. It excludes CAD surface tests, fittings, tendon contacts, opposite-side contacts, and the full movement domain. It is not a manufactured-clearance or contact-force calculation.

## Erector spinae and the next route solve

The saved erector-spinae route retains its loop and nine waypoints. Its selected actuator is one 40 mm BPA, resting length 538.96 mm, with 50 mm saved tendon and 25 mm fittings at each end. At the saved sample nearest neutral, its solved sleeve length is 471.85 mm and the calibrated diameter is 57.87 mm. The sleeve crosses three polyline corners under the provisional installation rule. The loop therefore needs an explicit physical interpretation before a clearance-constrained redesign.

There is also a saved-field discrepancy to reconcile before using this route as a physical installation: its saved closure contraction is positive, while its saved geometric `g.strain` is negative at every sample (about -55.8% nearest neutral). These fields have not been silently reconciled or recomputed. The displayed diameter uses `sContraction`; the geometric strain is shown separately. `saved_state_field_audit.csv` records this comparison for all sites. A taut flag alone is not evidence that the displayed sleeve, saved tendon length, and route geometry form a physically consistent installation.

Ben clarified that this older erector-spinae loop was intended as a reverse-pulley concept and will be reoptimized separately. Retain that intent. For the next optimization, resolve which segments are contractile and which are reeved tendon, with explicit fixed and moving blocks and their mounts. Re-solve resting length and tendon length when changing that layout. Evaluate the expanding radius at every movement sample; do not use nominal diameter for packaging. Keep parallel offsets as design variables or enforce spacing using the largest required diameter. Add fitted endcaps, brackets, and real sheaves before accepting a route.

The capsule screen supplies a reusable separation calculation. A clearance constraint has the form `required_clearance + radius_a + radius_b - segment_distance <= 0`. Its current CSV is screening evidence, not an installed constraint in the existing MATLAB optimizer. CAD constraints require registered geometry and surface-distance tests. Retain the existing torque, contraction, tension, and compliance requirements when adding these packaging constraints. No reoptimization, refit, or recorded-result replacement was performed in this pass.

## Rebuild

Run from the repository root. MATLAB exports fields from saved objects without constructing or updating them:

```powershell
matlab -wait -batch "run('Neuromechanical_Models/AH001_Simscape_Repair/scripts/export_packaging.m')"
& C:\Users\Ben\.anaconda3\envs\myoconv\python.exe Neuromechanical_Models\AH001_Simscape_Repair\scripts\build_packaging_view.py
& C:\Users\Ben\.anaconda3\envs\myoconv\python.exe Neuromechanical_Models\AH001_Simscape_Repair\scripts\screen_packaging.py
```

The exporter follows `KompareXXVII/route_configurations.csv` to the selected survey MAT filenames. Replace the selected results and regenerate the CSV roster before rebuilding. CSVs alone do not contain the saved contraction grids. Browser verification and clean media use the two accompanying Node scripts and the locally installed Playwright package; `make_packaging_gif.py` assembles the rendered frames.

Authoring scripts and this README are the reviewable source changes. `build/` and `output/` are ignored generated artifacts. Existing ZCode files, the dissertation sources, and optimizer results were not overwritten. No commits or pushes.
