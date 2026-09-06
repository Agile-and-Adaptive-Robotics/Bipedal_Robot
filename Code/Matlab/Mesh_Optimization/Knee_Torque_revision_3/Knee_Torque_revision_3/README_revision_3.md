# Revision 3 file guide

## Fixed versus calculated BPA radii

Edit `ctx.bpaRadiusMode` in `buildKneeFlexorContext20mm.m`.

```matlab
ctx.bpaRadiusMode = "scalar"; % use geo.bpaRb and geo.bpaRs directly
ctx.bpaRadiusMode = "bpaR";   % calculate one radius per knee angle
```

For fixed radii, select `"scalar"` and edit these lines in `buildGeoExclusion.m`:

```matlab
geo.bpaRb = 0.016; % wrap-point standoff, scalar or N-by-1 array
geo.bpaRs = 0.013; % collision radius, scalar or N-by-1 array
```

Do not replace the call to `bpaR` with zero. `bpaR(0,Dia,KMAX)` means the physical BPA radius at zero contraction. Scalar mode bypasses `bpaR` and holds the selected values over the motion.

In `"bpaR"` mode, `predictKneeFlexor20mm.m` calculates `predRadius = bpaR(bpa.strain_p,Dia,KMAX)`. It then uses `predRadius + bpaRbOffset` for wrap placement and `predRadius + bpaRsOffset` for collision. The route and contraction are iterated because each affects the other.

## Moving wrap point and bend length

`buildKneeFlexorRoute20mm.m` calculates `pWrap`, an N-by-3 array in the t1 frame. Its y-coordinate is fixed, its z-coordinate follows the p1-pEnd line in the t1 YZ projection, and its x-coordinate solves `signedDistanceToTibia20mm(pWrap) = bpaRb`. The array is returned as `routeInfo.pWrapT1`. At active-contact frames it is transformed and stored in `Location(2,:,i)`.

The current `MuscleLength` is the sum of the straight p1-pWrap-pEnd segments. `bendMeasure = wRap*abs(turnWrapped)` is an arc-length estimate used only by the Xi3 bend-loss calculation. There is no torus-contact centerline arc in the route length.

## Class corrections

`MonoPamDataExplicit_balance.m` and `MonoPamDataExplicit_balanceX3.m` in this folder are the corrected versions supplied after the earlier revision. They retain the `_i` result-variable pattern, `Force_p` for the predicted force-vector calculation, and `fortz` for the force-equilibrium solver. X3 does not contain `delta_L_final`.

`get.Force` is MATLAB's getter for the dependent property named `Force`; `Force_p` and `fortz` are separate local functions with distinct jobs.

## MonoPam_mult

`MonoPam_mult.m` follows the section order, calculation pipeline, and function naming of `MonoPamDataExplicit_balanceX3.m`, with cell-aware changes for two actual routes. Pass `Location = {Location1;Location2}` and `BPAcount = 2`. Resting length, contracted length, pressure, tendon length, tendon wraps, and bend measure are common values. There is one physical tendon on each route, but both use the same scalar tendon specification; tendon stiffness is not multiplied by two.

`predictOriginalKneeFlexor20mm.m` is intentionally different: it constructs the original two-point route with one BPA, passes `BPAcount = 1`, and passes `Xi3 = 0`. Thus the plotted original curve is the original one-BPA/no-Xi3 baseline, whereas the optimized curve is the summed torque from two actual mirrored BPA routes.

For route directions `u1` and `u2`, expressed in the same bracket frame, the common force magnitude produces:

```text
d_br = Fmag*(u1 + u2)*C_bracket
c_bracket(i) = ui*C_bracket*(u1 + u2)'
c_eff(i) = c_bracket(i) + 1/kSpr(i)
k_eff(i) = 1/c_eff(i)
```

Each tendon carries the common per-BPA force and therefore has the same calculated stretch `Fmag/kSpr`. Only the independent, per-angle `fzero` equilibrium solves inside `fortz` use `parfor`; the rest of the stiffness-aware pipeline remains in the same serial order as X3. The class reports `forceMismatch`, the difference between the two force-law predictions after enforcing one force magnitude. A nonzero value shows that the supplied routes/lengths are not mechanically symmetric enough for the equal-force assumption to be exact.

## Running the test

Place the revision files, `Vas_Pam_20mm_Result.mat`, `Bifemsh_20mm_Result.mat`, and the normal BPA dependencies on the MATLAB path. `FestoLookup.mat`, `festo4.m`, `maxBPAforce.m`, and `RowVecTrans.m` must be available. Then run:

```matlab
clear classes
TestMonoPam_mult
```

The Vastus test shifts every saved route point by +30 mm and -30 mm in z and compares the result with the saved parallel-BPA X3 configuration. The Bifemsh test mirrors the saved route about z=0 using its saved transformation matrix. It reports high-flexion paired/single torque, force, and torque-per-force geometry ratios so a torque deficit can be separated into force-law and moment-arm effects.

`TestMonoPam_multi.m` is included as a compatibility entry point for the earlier test name. It redirects to `TestMonoPam_mult`. The earlier `MonoPam_multi.m` wrapper is obsolete and should be removed from the MATLAB path.

## Named consolidation candidates

`predictKneeFlexor20mm.m` and `predictKneeExt20mm.m` currently create their own `pred` structures; there is no shared prediction-schema function. `plotStyle()` and `loadColors()` are local functions in `Opt_run_Ext.m`, while `Opt_run.m` defines comparable colors, line widths and axes settings directly in its plotting section.

If that duplication becomes difficult to maintain, `buildKneePredictionResult.m` and `kneePlotStyle.m` are reasonable names for future shared functions. They are proposals, not files included in this revision.
