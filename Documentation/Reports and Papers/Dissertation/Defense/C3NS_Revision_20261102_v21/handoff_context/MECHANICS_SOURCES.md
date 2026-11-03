# v21 mechanics assets and provenance

All outputs are in assets/mechanics/ relative to this revision. Dissertation, paper, model, measurement, and saved design sources were read only. No parameter identification, route optimization, manuscript edits, or simulation campaign was run.

## Assets for the main talk

| Source slide/use | Asset | Source and scope |
|---|---|---|
| 4: muscle force produces joint torque | Sutton2023_fig1B_swing.png / .pdf | Exact crop of Fig.1B in Sutton et al. (2023), obtained from the full publisher PDF in NLM's public archive. It is a **generic leg joint in swing**, not a knee-specific anatomical model. Final crop 867×839 removes the neighboring-panel fragment. |
| 8 and 31: historical BFLH comparison | historical_bflh_common_scales.png / .pdf / .fig | BFLH-only restoration of the archived 2019 source calculations. Human and robot share z and color limits of −60 to 0 N·m, camera, and displayed hip/knee axes. Original sampled hip domains remain distinct. Parent adds native Human/Robot labels above the panels. |
| 29: authentic model mechanics | opensim_knee_BFSH_overlay_landscape.png / .pdf | Actual Gait2392 femur/tibia meshes and BFSH path at −55° knee angle. The instantaneous knee axis comes from OpenSim body spatial velocity, including coupled translations. The exact signed force cross product agrees with computeMomentArm to floating-point precision. Portrait and JSON geometry are backups. |
| 33: single pinned example | pinned48_preoptimization_torque.png / .pdf / .fig | The 48 cm pinned-flexor configuration; 90 recorded points from its constituent repeated runs. Same measurements as the next length plot. Blue measured, gold hybrid, gray Hunt, peach BPA prediction. Blue circles are drawn last. |
| 34: length mechanism | pinned48_measured_vs_predicted_length.png / .pdf / .fig | **Actual recorded active BPA lengths**, not inferred from torque or a corrected model. Blue points are InflatedLength; peach line is saved rigid-geometry MuscleLength after fitting/tendon subtraction. Full curve fits within 380–500 mm. |
| 45: single-BPA 20 mm shortfall | biomimetic_flexor_publication.png / .pdf / .fig | Current dissertation native Components/figure_polish_20261030/FlxGrp.fig; use its 20 mm panel with the full legend. Parent owns the crop. Underlying values, colors, and axes are unchanged. |
| 48: redesigned flexor validation | redesigned_flexor_publication.png / .pdf / .fig | Dissertation native Components/directives_20261029/Redesigned_flexor_measured_last.fig. Blue measured circles are opaque and above yellow. Yellow is **prediction at measured pressure**, not the pinned-test hybrid calculation. |
| After 50: cam comparison | cam_force_and_pressure_comparison.png / .pdf / .fig | Left: **total force at 620 kPa**, directly from saved baseline/cam designs. Right: saved pressure to match selected **vastus medialis** torque. Do not call the left panel “required force.” |
| 51: physical cam validation | cam_extensor_physical_publication.png / .pdf / .fig | Saved nine-point publication audit plus saved final cam curve. Blue measured circles are opaque and drawn last. Gold is model prediction at measured pressure. Saved target is vas_med_r; legend says “Human model (VM).” |
| Before 52: reverse pulley | reverse_pulley_biomimetic.png / .pdf | Exact active dissertation Figures/30-torque/reverse_pulley_biomimetic.pdf. Two parallel BPAs and two cable parts: output force equals one BPA force and output travel doubles. No new optimization or waived individual muscle example. |

## Measurement and design provenance

- Testing_Data/2022_02_Festo/Plot_KneeFlxPin10mm_48cm.m names source load-cell and fish-scale sheets and explicitly records InflatedLength, Angle, Torque, ICRtoMuscle, and pressure data. Its saved MAT contains 90 measured length records, 0.398–0.467 m. The export asserts sorted angle/torque/length triples match the current evaluator's 48 cm record. Matching CSVs: pinned48_measurements.csv and pinned48_predictions.csv.
- The zero-correction pinned curve is a read-only call to current minimizeFlxPin(0,Inf,Inf,1). This evaluates the model once; it does not fit parameters. The Hunt curve follows the published companion generator's direct HuntEq call with its own 0.171 maximum-contraction basis. No legacy lookup-table pressure-axis error was reintroduced.
- Baseline/cam sources: Code/Matlab/Mesh_Optimization/Results/Vas_Pam_20mm_Result.mat and Vas_Pam_20mm_Result_cam_20261014_0416.mat. The exporter asserts baseline xBest equals the cam file's xRef8, and adopted χ values agree.
- Pressure source: Code/Matlab/Mesh_Optimization/Results/Ext_20mm_pressure_for_human_torque.csv, generated 2026-10-28. The source script explicitly uses final cam contact/tendon values. This matches the unmultiplied VM target; it is not the 5% overshoot target. Table copied to cam_saved_pressure_comparison.csv.
- Cam physical points: tmp/paper_bridge/figscripts/Ext20mm_series_audit_20261026.csv, nine ExtTest20mm_2 measurements stored to 0.001 N·m. This is the declared dataset used by the final publication figure, not the superseded six-point reduction. Curves come from predBest.TorqueZ and ctxCamBest.humanTorque. Original CSV copied alongside plot. Yellow points do not use measured BPA length or measured moment arms; those were not recorded for this test.
- ctxCamBest.targetName was directly checked as vas_med_r, consistent with column 5 of Testing_Data/2022_02_Festo/OpenSim_Vasti_Results.txt. The target is not a sum of vasti muscles.

## Historical BFLH restoration

build_legacy_bflh.m reads archived tmp/old_code_copy_20261022/ files: AttachPoints.mat, AttachPoints_Robot.mat, MovePoints.mat, MuscProp.mat, PAM_Prop.mat, Size.m, festo.m, and forz.m. It restores the two-coordinate BFLH calculation documented in OpenSim_kinematics_bi_hip.m and bipedal_kinematics_bi_hip.m, whose BFLH branches are currently commented out. Their equations and historical sizing rule are retained.

This is a regenerated historical comparison, not recovered raw test data or a newly optimized robot. The force-balance solver in the human muscle function evaluates the historical muscle law. Human sampled hip domain: −20 to 138.4°; robot: −15 to 108.75°. Both sampled knee domains: −120 to 8.7°. Result ranges: human −50.8089 to −3.12232 N·m; robot −32.6546 to 0 N·m. These reproduce the old plots' magnitude ranges. Numeric surfaces are preserved in historical_bflh_surfaces.mat.

## OpenSim overlay

- Model: Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim.
- Authentic meshes: mjc/gait2392_simbody/Geometry/femur_r.stl and tibia_r.stl.
- BFSH path comes from OpenSim getCurrentPath at the displayed pose. Force is applied to the first tibial path point after the femoral point; the downstream tibial segment does not define the joint-crossing force direction.
- The knee is a CustomJoint. Instantaneous center is computed from tibial angular/linear velocity at unit knee speed. Signed geometric cross product and OpenSim moment arm both equal −0.033654416847766454 m. Graphic shows M_z = (r_k × F)_z, with **vector** r_k. It does not label a projected perpendicular distance as the exact generalized moment arm.
- Missing VTP display-file warnings during OpenSim loading do not remove displayed geometry: the renderer explicitly reads the corresponding STL meshes above.

## Sutton citation and retrieval

G. P. Sutton, N. S. Szczecinski, R. D. Quinn, and H. J. Chiel (2023). *Phase shift between joint rotation and actuation reflects dominant forces and predicts muscle activation patterns*. PNAS Nexus 2(10):pgad298. [DOI](https://doi.org/10.1093/pnasnexus/pgad298). Figure 1B. CC BY 4.0.

The figure was cropped from the actual PDF, with no redrawing or inpainting. Full source PDF, source metadata, and crop script are retained. [NLM archive PDF](https://pmc-oa-opendata.s3.amazonaws.com/PMC10563792.1/PMC10563792.1.pdf). Current archive metadata verified PDF checksum and open-access status. The older 739 px preview is only a source-discovery artifact; use the final PDF crop / 867 px PNG.

## Reproduction and visual checks

- build_mechanics.m: pinned pair, cam force/pressure, and native publication exports.
- build_legacy_bflh.m: common-scale historical surfaces.
- build_cam_physical.m: saved nine-point cam validation.
- render_opensim_knee.py: actual model projection, exact moment-arm check, and wide closeup.
- crop_sutton.py: source-PDF crop. Ghostscript renders final PNG at 400 dpi.
- Every delivered plot and illustration was opened and inspected individually. Blue/yellow overlap, axis clipping, source-panel residue, and the wide OpenSim frame were corrected after inspection. These are asset-level checks; parent must still inspect actual slide renders after composition.

## Native render follow-up

- `single_bpa_shortfall_clean.png/.pdf/.fig` copies the native 20 mm axis from `biomimetic_flexor_publication.fig` into its own figure. It removes the panel letter and internal title, retains all three curves and both measured point groups, and renames the model legend entries to `Rigid-body model` and `Corrected model`. This supersedes cropping the two-panel raster for the single-BPA shortfall slide.
- `assets/force/maximum_force_pressure_10.png/.pdf/.fig` now has its legend below the axes. Both labels have explicit normalized positions to avoid the old custom label positions colliding with the legend or output edge.
- Reproduction: `assets/mechanics/clean_review_panels.m`. Assertions compare every X/Y/Z/color array and curve/marker style before and after copying/layout. Both final raster exports were inspected individually. No deck modules or manuscript files changed.
