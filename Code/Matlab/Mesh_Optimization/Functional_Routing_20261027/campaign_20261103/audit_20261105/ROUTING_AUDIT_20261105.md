# Routing audit and repaired candidates — 2026-11-05

The previous campaign implemented mathematical BiPulley gains and parallel actuator forces. It did **not** demonstrate a complete, interference-free, target-matching mechanical layout. The signed-results report already acknowledged substantial failures; the additional geometry and bookkeeping errors below invalidate treating its feasibility labels as packaging approval.

The separate repaired pass covers all 27 groups and both rigid and Xi0–Xi3 compliant branches. Of 54 branch results, **8 pass the implemented global hard constraints**, and **0 both pass those constraints and meet every sampled signed torque target**. These are bounded search results, not a proof that a better design cannot exist. The solver searches sparse joint grids; full simultaneous multi-DoF qualification is not established.

## What the latest ZCode results actually contained

- Source of record: `../../SIGNED_RESULTS_20261104.md` in the Functional_Routing_20261027 directory, superseding the older RUN1 report; actual `run_1/checkpoints` were also read and exported. Old results were preserved.
- Parallel counts and transmission gains were implemented separately. Examples in saved checkpoints: vasti 3 parallel BPAs/G2, soleus 2/G2, gluteal aggregate 4/G2. These were not automatic four-actuator substitutions.
- The sign-blind objective had already been corrected in the signed campaign. Earlier magnitude-only success claims were superseded by that report.
- The report explicitly acknowledged that robot-CAD distance fields and the vasti cam were not connected to the optimization equilibrium. Mathematical gain is not a documented physical reeving design.
- The in-layer word `solved` meant constraint status under an incomplete geometry model. It did not mean full torque coverage or all-group packaging success.

## Verified defects repaired for this pass

| Defect | Consequence | Repair/evidence |
|---|---|---|
| Moving-path helper populated `out` without returning it | NaN insertion seed could silently pick the first cloud vertex | Return the updated rows; reject non-finite seeds |
| Surface patch used the second SVD output (singular values) as surface directions | Search axes and normal were geometrically wrong | Use the third SVD output; rotated-plane fixture |
| Patch displacement was added twice | Attachments could leave the intended local search region | Apply only the tangent residual once and clamp the disk radius |
| Fixed Z projection used for every joint coordinate | Hip adduction, rotation and translating knee torques scored incorrectly | Coordinate wrench projection; independently agrees with native OpenSim moment arms for seven leg coordinates, maximum error 1.36e-11 m |
| Diameter used whole-route strain rather than solved vessel contraction | Reverse-transmission shortening produced misleading swelling | Use solved `sContraction/L0`, with explicit contraction ceiling |
| Parallel envelope capped at three channels | Fourth and further actuators disappeared from packaging checks | Emit every selected channel; apply the declared bundle spread consistently |
| Two-point routes could have zero-length sleeve envelopes | Large pressure vessels disappeared from collision checking | Emit full solved vessel lengths, including candidate pinned end sections |
| Local bone clouds compared directly with world routes | Clearance changed with inconsistent frames | Transform obstacles with their owning bodies at each pose |
| Retained groups used their own pose indices | Collision checks compared different configurations | Re-evaluate retained envelopes at the current candidate's poses |
| Segment-pair clamping/sign errors | Skew or overlapping parallel segments could be assigned the wrong distance | Exact convex minimum over the interior and all four edges; crossing, skew, parallel, endpoint and point tests |
| Fixed seed tendon procurement ignored optimized tendon length | Tendon-length search could not resolve required spans | Compare required spans with the chosen tendon length |
| Foot rule inspected only a guide | A vessel could occupy the foot despite passing the check | Test full sleeve/fitting extent against a conservative ankle plane |
| End fitting envelopes did not cover their full chosen length | A clear pressure vessel could hide fitting intrusion | Extend candidate collision segments and use enclosing fitting spheres in final refinement/checks |

Interrupted exploratory repair passes are retained as `reoptimization_pre_pca_fix` and `reoptimization_pre_full_sizing`; neither is the exported repaired result. Source hashes and the original source snapshots are retained. The 35 mm maximum bracket-pin reach is a **proposed packaging assumption**, not an approved bracket design; bone-root patch radii were not enlarged.

## Repaired search scope and results

The run uses three local MATLAB workers. It preserves each anatomical region's deep-to-superficial order, then performs one coordinated all-group refinement round. The preliminary seed/configuration search retained the inherited segment-distance and short fitting-envelope conventions; the final coordinated round restarts all workers and uniformly uses the verified segment-distance routine and full fitting extents. Each final branch re-screens both diameters, counts near its selected count, gains 1/2/3 and the three tackle families (up to 30 seconds), then searches continuous geometry (up to 200 pattern evaluations/20 seconds). The combined-vasti output is then constrained to 20 mm and reoptimized (up to 30 seconds for count/gain/family screening and 300 pattern evaluations/40 seconds for continuous geometry); all groups are checked again against that layout. Final checks use matching rigid layouts for the rigid branch and matching compliant layouts for the compliant branch, at simultaneous poses. The preliminary pass screens both actuator diameters, parallel counts, gains and all three tackle families at seed and old warm-start designs, then searches continuous attachments, lengths and split variables for the winning screen configuration. It does not fully optimize every discrete combination: up to 60 surrogate evaluations/12 seconds and 300 pattern evaluations/40 seconds per branch are explicit first-pass budgets. Heavy aggregate groups additionally screen 8, 12 and 16 parallel actuators. The final re-screen is local and time bounded; a wider joint discrete/continuous search is still needed before claiming an optimal selection. The search is not a global optimum certificate.

| Group | Rigid target poses met | Compliant target poses met | Globally constraint-feasible (R/C) |
|---|---:|---:|---|
| ercspn | 0% | 0% | yes/yes |
| intobl | 0% | 0% | yes/yes |
| extobl | 0% | 0% | no/no |
| bifemlh | 0% | 0% | no/no |
| semimem | 0% | 11% | no/no |
| sar | 100% | 44% | no/no |
| tfl | 56% | 89% | no/no |
| rect_fem | 0% | 0% | no/no |
| grac | 56% | 67% | no/no |
| bifemsh | 0% | 0% | no/no |
| vas_int | 0% | 0% | no/no |
| soleus | 0% | 0% | no/no |
| med_gas | 0% | 0% | no/no |
| lat_gas | 0% | 0% | no/no |
| tib_post | 0% | 0% | no/yes |
| tib_ant | 0% | 0% | no/no |
| per_brev | 22% | 0% | no/no |
| per_long | 0% | 0% | no/no |
| per_tert | 0% | 0% | no/no |
| flex_dig | 0% | 0% | no/yes |
| flex_hal | 0% | 0% | no/no |
| ext_dig | 0% | 100% | no/no |
| ext_hal | 78% | 89% | no/no |
| glut_max3 | 11% | 0% | no/yes |
| add_mag3 | 0% | 0% | no/no |
| iliacus | 0% | 22% | no/no |
| psoas | 0% | 56% | no/yes |

`repaired_result_summary.csv` gives the worst constraint for each branch. `repaired_torque_curves.csv` gives every scored target/robot pair. `repaired_contraction_and_packaging.csv` separates nominal size, parallel count, mathematical gain, physical contraction and rendered vessel dimensions.

`Repaired_Torque_Comparison.pdf` and `.png` plot all 27 groups against their signed human targets. `repaired_sampled_hardware_loads.csv` reports the sampled all-pass tendon load, bracket reaction and reversing-sheave 2T bound; these are model loads, not an approved assembly rating. The inherited objective's catalog load constraint is primarily scored on the primary pass, so the all-pass export must also be reviewed.

Independent 25-cell checks are in `repaired_dense_verification.csv`. Human targets there are explicitly resampled from the original target grid, rather than recomputed human-muscle equilibria. Pose-to-pose torque and sleeve-length jumps are reported without dismissing them based on a guessed conditional-point location.

## What the visuals show

1. `gait2392_robotbody_latest_compliant_blue.osim` and its rigid counterpart reproduce the prior saved primary-pass routes. Original paths are gray; changed candidate paths are blue. They are display-only PathActuators, so they do not add muscle force to the model.
2. `gait2392_robotbody_repaired_compliant_blue.osim` and its rigid counterpart show the repaired candidates. Each serialized frame and blue color is reopened and checked through OpenSim.
3. `Repaired_Routes_OpenSim_Matched.blend` contains compliant and rigid scenes, the same OpenSim meshes and body transforms, gray reference paths, blue candidate routes, and all solved pressure-vessel channels/sections. Colored vessels change length and diameter from the solved state; fitting solids are envelopes. Tackle wireframes are envelopes with unresolved reeving and are hidden in final renders. The white world/compositor, 180-degree camera orbit and 0 to -120 degree knee cycle are saved. The original CAD foot-registration transform is not reused.

The knee cycle moves only the knees. Hip-only, ankle-only and trunk muscles need not contract during that cycle. Primary-joint contraction ranges are separately reported in the CSV; lack of visible movement is not replaced with fabricated contraction. Physical vessel placement is an explicit route-tangent candidate, so protrusions and collisions remain visible when the optimizer fails.

The saved Blender file was reopened headlessly. Evaluated pressure-vessel length/radius was checked at three knee poses for every channel/section (429 mesh checks), together with the native bone frames, white world and 180-degree camera action. `blender_reopen_verification.json` records the numerical errors. `Repaired_Packaging_Closeup.png` gives a closer knee/foot view; an optional packaging-closeup camera is also saved in each scene.

The stock and robotbody JointSets are identical. A direct MATLAB/native body-frame comparison is saved in `matlab_vs_native_body_frames.json`; the inherited natural-cubic joint interpolation differs from native OpenSim by up to **0.2349 mm** on the displayed knee grid. The native OpenSim bone/path transforms are used for rendering; world-space pressure-vessel candidates retain the MATLAB mechanics frames. The moment-arm verification checks the generalized-wrench formula independently; it does not certify every inherited route approximation.

`requirement_audit.csv` records each requested feature separately, and `target_constituents.csv` lists the exact muscle aggregates used. The dissertation was checked directly: `Final/chapters/30-torque.tex` and `94-AppendixC.tex` document directional stiffness, including the biomimetic-flexor origin assignment [Xi1,Xi2,Xi2] and insertion [Xi1,Xi2,Xi1]. The campaign transfers those prototype assignments to the new routes; this is not calibration of every new bracket.

`initial_sizing_screen.csv` evaluates your approximate half-Fmax10 rule using the longest constituent path and the repository's 620 kPa resting-force equation. It is an initial screen, not a substituted torque/packaging feasibility test; the repaired search still compares both allowed diameters.

## Remaining blockers — these results are not defense-ready design validation

- **Vasti cam:** no torque-coupled cam equilibrium or final patellotibial cam/bracket design. Vasti still uses the campaign's direct-route approximation. This is a missing requirement, not a solved feature.
- **Series knee sections:** geometric end sections honor the selected total sleeve length and split. Independent force/strain equilibrium of unequal sections has not been established. A central path is not a qualified Bowden assembly. For knee-crossing routes ending on the foot, placing an active end section at the anatomical insertion conflicts with the no-foot-BPA rule; a separate tibial bracket/stop and distal tendon route are needed. The current fixed topology does not supply those independent mount variables, so its failures are retained rather than relaxing the foot ban.
- **Pulleys:** gain, parallel count and series sections are distinct. The inherited `sheaveCountPhysical` metadata incorrectly inferred a count from gain; the wrapper now explicitly marks that value unresolved (NaN), and neither force class consumes it. Preserved old audit records retain their original metadata. Exact reeving, pin reactions and catalog assembly ratings remain unqualified.
- **Robot CAD and bones:** newly sampled signed fields cover closed displayed femur, tibia and talus meshes, on a 3 mm grid with a conservative cell-diagonal allowance of approximately 5.2 mm plus segment sampling allowance. Reported near-contact violations include those numerical allowances; they are not exact measured physical penetrations. Open pelvis/foot/torso surfaces retain sparse-cloud checks, which can miss penetration. The robot-CAD registration has large knee/ankle mismatch and is not silently accepted. Own tackle versus bone, bracket geometry and complete rigid-to-rigid collision qualification are also incomplete.
- **Attachment kinematics:** the NaN-return bug is fixed, but some moving/conditional path points remain frozen or dropped by the inherited pass builder. A blue static point in OpenSim does not reproduce a missing dynamic wrapping model.
- **Bundle and fitting placement:** parallel channels use a fixed ring/spread convention around the common route, not independently optimized multi-layer branch placements. Fitting spheres enclose the chosen end-fitting extents; they are not detailed catalog fitting, pin and bracket geometry, and the inherited case builder shortens fitting assumptions on short routes. Bracket deflections are coupled to the force solve but the preview does not fully reconstruct every displaced hardware surface.
- **Load components:** the existing three candidate tackle families and bicycle-style Bowden placeholders are not vendor-qualified for the computed loads. In particular housing compression, tendon strength, fatigue, minimum bend radius and assembly WLL need verified components. No invented McMaster part approval is claimed.
- **Aggregation:** the campaign's target contract is reused (combined vasti, hamstring/adductor/gluteal targets). This run does not invent a three-way gluteal repartition or anatomically qualified psoas-to-iliacus joining bracket. Inspect the exact target contract before treating that topology as implemented.
- **10 mm diameter:** the inherited 10 mm external-diameter law is 13.6 mm at rest to 17.837 mm at 16.5% contraction; this is different from a nominal 10 mm bore. It needs confirmation from the actual actuator specification. The 20 mm user calibration is 20 to 38 mm at 24% contraction.
- **Computing:** SSH to easteregg2 was unavailable in this session; the repaired search ran locally. The bounded search cannot be represented as an overnight remote exhaustive run.

The next valid campaign needs the cam and physical transmission/series topology connected to equilibrium, verified full robot geometry registration, and component load qualification before spending a larger optimization budget. Longer runs of an incomplete packaging model cannot establish those requirements.
