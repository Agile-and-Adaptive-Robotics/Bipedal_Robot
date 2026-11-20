# Folder: Mesh_Optimization
This folder contains scripts that optimize the PAM placement by using points contained in the mesh files that define the structure of the OpenSim bone models.
Note: Add folder \Bipedal_Robot\Code\Matlab\Mesh_Optimization\ and subfolders to Matlab path.  At a minimum add subfolders \Functions, \Mesh_Optimization, \Mesh_Optimization\Results, and Robot_Data. Older scripts contain file paths specific to my Connor’s PC. These can be updated to be more general, as they are only using scripts and functions contained within the GitHub repository.
# Files and Descriptions:
## Newer optimization files (Bolen et al. 2026)
### For the Flexor:
Opt_run.m 		% run file for optimizer
Opt_sanity.m 		% Sanity check before running full optimizer
buildKneeFlexorContext20mm.m    % build/load constants once
predictKneeFlexor20mm.m         % x -> BPA torque/length/strain
objective_KneeFlexor20mm.m      % x -> scalar optimizer cost
buildGeoExclusion.m		%Geometry for p2 exclusion
nonlconExclusion.m		%Nonlinear exclusion constraint for p2
objconstrExclusion.m       % Objective/constraint wrapper for surrogateopt
buildKneeFlexorRoute20mm.m    % build the wrap/release route using the selected BPA radii
predictOriginalKneeFlexor20mm.m    % evaluate the original route used for comparison
signedDistanceToTibia20mm.m    % measure clearance from the simplified tibial surface
### Related flexor files
bpaR.m    % calculate the BPA radius from physical contraction
minimizeFlxPin10_results_20260730_2transforms_Z2.mat    % load identified Xi0, Xi1 and Xi2 values
minimizeExtPin10_results_20260819_2transforms_Z2.mat    % load the identified Xi3 value
Bifemsh_20mm_Result.mat    % optional saved optimizer output for subsequent flexor analysis
### For the Extensor:
Opt_run_Ext.m 		% run file for optimizer
Opt_sanity_Ext.m 		% Sanity check before running full optimizer
buildKneeExtContext20mm.m    % build/load constants once
predictKneeExt20mm.m         % x -> BPA torque/length/strain
objective_KneeExt20mm.m      % x -> scalar optimizer cost
objconstrExt20mm.m    % package objective and constraints for surrogateopt
nonlconExt20mm.m    % check extensor route and length feasibility
buildDistalRingLocation20mm.m    % build the changing distal-ring routing points
### Related extensor files:
minimizeExtPin10_results_20260819_2transforms_Z2.mat    % load the identified Xi0-Xi3 values
Vas_Pam_20mm_Result.mat    % optional saved optimizer output for subsequent extensor analysis
###  Shared reference files
MonoPamDataExplicit_balanceX3.m    % calculate stiffness-aware BPA geometry, force and torque
MonoPam_mult.m % two-route equal-force/shared-bracket BPA calculations
TestMonoPam_mult.m % mirrored-route compliance and torque test
festo4.m    % evaluate normalized BPA force at contraction and pressure
maxBPAforce.m    % calculate the maximum force used to scale BPA force
RowVecTrans.m    % transform routing points between reference frames
RpToTrans.m    % assemble rotation and translation into a transform
### Implementation
### How the newer optimizer is organized
ctx means context. It is an ordinary MATLAB struct, not a special MATLAB object. A context groups values that stay fixed during one run: angle grids, transforms, geometry, BPA and stiffness constants, target data, penalty settings, initial values, and bounds. The eight values that the optimizer changes remain in x = [p1(1:3), pEnd(1:3), rest, tendon].
Related containers – geo is the geometry portion of the context; pred is the prediction returned for one candidate x; routeCtx is the smaller saved flexor context needed to rebuild the accepted route later. Saving ctx or routeCtx with xBest prevents a later analysis from silently using changed constants or geometry.
### Flexor optimizer call flow
- Opt_run.m calls buildKneeFlexorContext20mm once to create ctx and ctx.geo.
- surrogateopt calls objconstrExclusion. That wrapper evaluates predictKneeFlexor20mm once and sends the same radius-coupled prediction to objective_KneeFlexor20mm and nonlconExclusion.
- predictKneeFlexor20mm builds the route with buildKneeFlexorRoute20mm, constructs MonoPamDataExplicit_balanceX3, and returns torque, length, strain, wrap, and bpaR-radius diagnostics.
- patternsearch then refines the global result using objective_KneeFlexor20mm and nonlconExclusion with the same ctx.
- After feasibility checks, Opt_run.m packages routeCtx. Save xBest, routeCtx, and Xi3 in Bifemsh_20mm_Result.mat for Knee_Flexor_Data_20mm.m.
### Extensor optimizer call flow
- Opt_run_Ext.m calls buildKneeExtContext20mm once to create ctx, including the nine-row route definition, extensor CAD geometry, OpenSim target, Xi0-Xi3, and BPAcount.
- surrogateopt calls objconstrExt20mm, which packages objective_KneeExt20mm and nonlconExt20mm in the structure required by surrogateopt.
- predictKneeExt20mm builds the changing route with buildDistalRingLocation20mm, evaluates MonoPamDataExplicit_balanceX3, and returns the same categories of prediction diagnostics.
- patternsearch refines the selected starting point. The adjusted-seed block in the supplied Opt_run_Ext.m is active and overwrites xBest unless it is commented out.
- Save ctx and xBest together in Vas_Pam_20mm_Result.mat. Knee_Extensor_20mm.m intentionally reloads the saved ctx rather than rebuilding it.
### Post-optimization, experiment, and plotting workflow
- Flexor: run Opt_run.m; accept and save Bifemsh_20mm_Result.mat; run Knee_Flexor_Data_20mm.m (called Knee_Flexor_20mm.m in some working folders); then save the resulting theoretical BPA objects and geometry in the test-specific MAT file.
- Extensor: run Opt_run_Ext.m; accept and save Vas_Pam_20mm_Result.mat; run Knee_Extensor_20mm.m; then save the resulting theoretical BPA objects and geometry for the planned test.
- Perform the physical experiment and add the measured knee angle, torque, pressure, inflated length, and moment-arm data to the test-specific result file.
- Run a comparison script such as Plot_KneeFlx_20mm_42cm.m. These plotting scripts combine the saved class predictions, experimental measurements, hybrid force-times-moment-arm calculations, and OpenSim reference data; they do not rerun the route optimizer.
### Relationship to the Xi0-Xi3 minimizers
- minimizeFlx, minimizeFlxPin, minimizeExt, and minimizeExtX3 identify or validate Xi stiffness/length-loss parameters against saved experimental and hybrid data. Their kf and ke structures are data containers whose core geometry, BPA inputs, and predicted fields correspond closely to the MonoPamDataExplicit class properties. Experimental fields such as Aexp/Mexp and hybrid fields such as A_h/Lm_h/mA_h/M_h remain outside the class because they describe a particular test, not the BPA model itself.
- The bracket frame is configuration-specific. minimizeExt, minimizeExtX3, and minimizeFlx update a femur/hip-side point through a bracket frame expressed relative to that body. minimizeFlxPin updates a tibia/knee-side point through a tibial bracket frame. This is why bracket origin, rotation, parent body, and deformed Location row must be explicit before one common stiffness routine can safely cover every experiment.
### Consolidation without creating many more files
- Good consolidation candidate – the duplicated robot-knee transform and angle-grid setup in buildKneeFlexorContext20mm and buildKneeExtContext20mm. Consolidate it only after a regression check confirms identical T_Pam, T_ICR_t1, T_t1_ICR, T_t1_f, phiD, and pos outputs.
- Good consolidation candidate – common prediction fields and plotting style. predictKneeFlexor20mm.m and predictKneeExt20mm.m currently create their own pred structures; no shared prediction-schema function currently exists. plotStyle() and loadColors() are local functions in Opt_run_Ext.m, while Opt_run.m defines comparable colors, line widths, and axes settings directly in its plotting section. If this duplication becomes difficult to maintain, proposed shared functions could be named buildKneePredictionResult.m for the common result fields and kneePlotStyle.m for colors and formatting. Those two names are proposals, not existing files.
- Keep separate – buildKneeFlexorRoute20mm versus buildDistalRingLocation20mm, and nonlconExclusion versus nonlconExt20mm. Their geometry, activation logic, and collision constraints are genuinely different.
- Do not consolidate the femur-side and tibia-side bracket calculations until the bracket parent frame, bracket transform, stiffness principal axes, and deformed route row are explicit inputs. Hiding those differences inside a generic helper would make the result harder to verify.
- Recommended boundary – keep Opt_run.m and Opt_run_Ext.m as readable orchestration scripts; keep physics in the MonoPamDataExplicit_balance class; keep route and collision geometry in the side-specific builders/constraints; and add a shared helper only when it replaces a substantial duplicated block rather than adding another thin wrapper.
## Used by Morrow et al. (2020)
Muscle_Name_Mesh_Opt.m – These scripts use the skeletal mesh from OpenSim to determine the best placement for PAMs in order to replace one or more muscles. Each script follows the general structure:
- Add file paths to the human muscle class structure, robot PAM class structure, functions folder, and OpenSim data folder
- Creates transformation matrices for the joints that the muscles will cross over. These are made by creating a linearly spaced vector of angles between the joints minimum and maximum orientations. Those angles are then used to create a rotation matrix and then a transformation matrix by using the distance between the first point in a muscle’s reference frame and the following joint frames that the muscle crosses over into.
- Runs the class structure for the human muscles.
- Runs the class structure for the PAM based on an initial starting point. This initial starting point is typically the same as the human muscle’s locations, but may have one or more points removed to make a more linear path.
- If the joint has more than one degree of freedom, the torque tensor is unstacked into multiple torque matrices that only capture information about one axis.
- If the PAM is replacing a group of muscles, those human muscles have their torques added together for comparison later against the PAM torque.
- The cost function runs to get a starting value of the cost.
- The mesh optimization begins, by importing location data from the skeletal locations that create the model mesh. For each run through the optimization loop, one location is updated to be one of the locations from the mesh. Every permutation between two bones that the muscle connects to is used, with the cost value being updated for ever run, looking for the lowest value.
- Once the optimization is finished, the best location is updated. Torque profiles are generated for the human and optimal pam location. A plot of the angle between the human torque vector and the pam torque vector is generated. A plot of the human skeleton with the human muscle in red and the pam in green is created.
These scripts can be run without any changes on the users part. These are a starting point for optimization. The results from these scripts should always be the same, since the locations they are pulling from are constant, so the solution should be a global minimum for the data set.
The intention was for these scripts to serve as a starting point, for both the optimization code and for generating decent locations for physical PAM attachments. This code can be improved by adding onto it a more traditional machine learning algorithm, in which after the best location is found from the mesh data points, that would serve as a starting point to do a first order gradient descent search around that point. The implementation of this is in “Bifemsh_Mesh_Refine_Opt.m”.
Bifemsh_Mesh_Refine_Opt.m – This script is an example extension of the Mesh_Opt code that was discussed previously. It contains all of the previous code, but adds a component that uses a machine learning algorithm to search the space around the best attachment location based on the skeletal data points. The algorithm looks at 8 points around the muscle locations before and after a point and calculates the cost function between all 64 connections using those points. Once a new minimum is found between all of those points, it uses the minimum cost locations to repeat the process. It outputs the same form of results as the previous code.
costFunctionKnee.m – This cost function is an adjustment from the generic one found in the “Functions” folder. Because the knee only actuates in one direction, this cost function more heavily weights the torque and torque angle around the z direction of the knee, and reduces the value of those parameters in the other directions.
Inputs: This function takes the human muscle torque and the PAM replacement torque as inputs.
Outputs: The output of this function is a C or cost value. The value is unitless and can be scaled to any arbitrary amount through the use of the gains Gt and Ga.

---

> Agent-readable copy of the .docx readme (the .docx is the canonical human version).
> After editing this file, also update the .docx.
