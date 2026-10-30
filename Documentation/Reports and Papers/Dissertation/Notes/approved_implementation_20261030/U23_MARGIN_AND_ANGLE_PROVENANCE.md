# U23: direct-route margins and worst-case angle resolved

30 October 2026. Read-only extraction from the adopted direct-route MAT files. No force model, objective, optimizer, or prediction was called. The only postprocessing was locating the minimum in the already saved flexor margin array, looking up saved angles, extracting extrema of saved human-reference arrays, and formatting the existing fractions/angles for display. Chapter 4 was not edited by this reviewer.

Detailed raw values and MAT SHA256s are in `U23_SAVED_DIAGNOSTICS_EXTRACTION.json`; the extraction script is retained beside it. The earlier `DIRECT_ROUTE_SAVED_FIELDS.json` lists the same saved variables but replaces nonscalar arrays with dimension descriptions, so the original MAT arrays were read to recover their actual entries.

## Values for Table 4.10

| Route | Stored minimum margin fraction | Percent display | Worst-angle source | Angle for consistent one-decimal display |
|---|---|---|---|---|
| Flexor |0.0500687314632795 |5.01% | `Bifemsh_20mm_Result.mat`: argmin of stored `torqueMarginFraction` over stored `validHumanTorque` gives MATLAB index 11; `ctx.phiD(11)` is −106.86868671763432 degrees. |−106.9° |
| Extensor |0.05024928158006919 |5.02% | `Vas_Pam_20mm_Result.mat`: stored `iWorst=77`, stored `minMargin`; `ctx.phiD(77)` is −20.20201995901284 degrees. |−20.2° |

Both saved contexts contain `ctx.requiredTorqueMargin=0.05`. Their stored `ctx.phiD` arrays each contain 100 knee angles spanning approximately −120° to +10°. All 100 flexor entries are included by the stored `validHumanTorque` mask. The flexor file lacks a dedicated stored worst-angle/index scalar; its reported angle is explicitly the argmin/index lookup from the precomputed saved diagnostic array, not a new prediction.

The direct-route diagnostic is the pointwise margin over the **plain human-reference torque magnitude**:

`margin_j = abs(robot_torque_j) / abs(human_reference_j) - 1`.

For the positive extensor curves, the saved code writes the equivalent `(robot-human)/human` form. The minimum is taken over the stored evaluated knee-angle grid. It is neither a margin divided by the 1.05-times-human envelope nor a margin divided by `ctx.torqueScale`.

Therefore comparing achieved minima 5.02%/5.01% with a required 5% design allowance is correct when **both margins are stated relative to the human reference**. The design requirement is 1.05times that reference magnitude. Suitable row labels are “Minimum predicted margin above human reference” and “Required design allowance above human reference,” with percent indicated once. No new allowance or experimental uncertainty is inferred.

## The torque-range row contains plain human reference values

The extrema below are entries of the saved `ctx.humanTorque` arrays:

| Route | Exact saved minimum / maximum (N m) | Existing table display |
|---|---|---|
| Extensor |20.24783502 /50.32136527 |20.25 to 50.32 |
| Flexor |−23.15874514 /−2.46121431 |−23.16 to −2.461 |

The separately saved interpolated human-reference arrays agree at the table's displayed precision. These are the unscaled human references, so **“Human reference torque”** is supported for that row; the existing values should be retained. Calling them the 1.05 requirement envelope would misidentify them.

## Code and field provenance

- `Code/Matlab/Mesh_Optimization/Opt_run.m:719–725` constructs the flexor diagnostic by dividing absolute predicted torque by the interpolated plain `ctx.humanTorqueAbs`, then subtracting 1. Lines 728–731 report its minimum. The extraction above uses the stored resulting array, without rerunning that construction.
- `Opt_run_Ext.m:244–249` interpolates plain `ctx.humanTorque` into `humanBest`;307–314 define the direct extensor margin and save/use `minMargin`, `iWorst`, and the corresponding `ctx.phiD` entry.
- `buildKneeFlexorContext20mm.m:128–135` reads the biceps-femoris-short-head reference column directly into `ctx.humanTorque`;144 sets the separate 0.05 allowance. `buildKneeExtContext20mm.m:120–139` reads the selected human vastus column directly;152 sets the separate allowance.
- `objective_KneeFlexor20mm.m:25` and `objective_KneeExt20mm.m:37` form the required magnitude as `(1+ctx.requiredTorqueMargin)*humanAbs`. That requirement is separate from the plain-human denominator of the saved diagnostic margin.
- The cam objective's `marginVsRequired` at `objective_KneeExtCam20mm.m:49–51` instead subtracts the required envelope and divides by `ctx.torqueScale`. It is a different quantity. This U23 resolution does not settle or change the separate cam-margin issue T01.

MAT files read: `Code/Matlab/Mesh_Optimization/Results/Bifemsh_20mm_Result.mat` (SHA256 `178107C1BB1D0721A9E0DFC89061646F7BCDABE0A6A25450F2DEC313FD6B73B6`) and `Vas_Pam_20mm_Result.mat` (SHA256 `057F787FBF8C328C937D4B9A0BEBC6727922F99D5C1E789E0E7D80E37B1DC21F`). No stored value or MAT file was changed.

