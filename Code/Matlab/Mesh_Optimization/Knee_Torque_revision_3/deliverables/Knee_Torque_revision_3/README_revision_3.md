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

`MonoPamDataExplicit_balance.m` and `MonoPamDataExplicit_balanceX3.m` use your `_i` result-variable pattern. Lowercase `fortz` is the BPA force-vector function. The six-input deformation solver is named `equilibrium`, so it cannot collide with `fortz` or the dependent `Force` property. X3 no longer contains `defer` or `sharedGeometry`, and no `pAold` variable remains.

`get.Force` is MATLAB's getter for the dependent property named `Force`; `fortz` is a separate helper function. There are no two local functions named `Force`.

## MonoPam_mult

`MonoPam_mult.m` is a standalone two-route class rather than a wrapper around X3. Pass `Location = {Location1;Location2}` and `BPAcount = 2`. Shared values are used for resting length, contracted length, pressure and stiffness. Tendon length, tendon wraps and bend measure may be shared numeric values or two-element cells.

`predictOriginalKneeFlexor20mm.m` is intentionally different: it constructs the original two-point route with one BPA, passes `BPAcount = 1`, and passes `Xi3 = 0`. Thus the plotted original curve is the original one-BPA/no-Xi3 baseline, whereas the optimized curve is the summed torque from two actual mirrored BPA routes.

For route directions `u1` and `u2`, expressed in the same bracket frame, the common force magnitude produces:

```text
d_br = Fmag*(u1 + u2)*C_bracket
c_bracket(i) = ui*C_bracket*(u1 + u2)'
c_eff(i) = c_bracket(i) + 1/kSpr(i)
k_eff(i) = 1/c_eff(i)
```

Each tendon stiffness remains independent. Knee angles are solved with `parfor`; both BPAs at one angle remain inside the same iteration because their forces jointly deform the bracket. The class reports `forceMismatch`, the difference between the two force-law predictions after enforcing one force magnitude. A nonzero value shows that the supplied routes/lengths are not mechanically symmetric enough for the equal-force assumption to be exact.

## Running the test

Place the revision files and normal BPA dependencies on the MATLAB path. `FestoLookup.mat`, `festo4.m`, `maxBPAforce.m`, and `RowVecTrans.m` must be available. Then run:

```matlab
clear classes
TestMonoPam_mult
```

The rigid-bracket part checks force agreement and cancellation of the two off-axis torque components. The finite-bracket part checks `c_eff = c_bracket + 1/kSpr` and `k_eff = 1/c_eff`. It does not require the two projected bracket compliances to be equal because the bracket principal-axis frame can be tilted relative to the femur z=0 symmetry plane.

`TestMonoPam_multi.m` is included as a compatibility entry point for the earlier test name. It redirects to `TestMonoPam_mult`. The earlier `MonoPam_multi.m` wrapper is obsolete and should be removed from the MATLAB path.

## Named consolidation candidates

`predictKneeFlexor20mm.m` and `predictKneeExt20mm.m` currently create their own `pred` structures; there is no shared prediction-schema function. `plotStyle()` and `loadColors()` are local functions in `Opt_run_Ext.m`, while `Opt_run.m` defines comparable colors, line widths and axes settings directly in its plotting section.

If that duplication becomes difficult to maintain, `buildKneePredictionResult.m` and `kneePlotStyle.m` are reasonable names for future shared functions. They are proposals, not files included in this revision.
