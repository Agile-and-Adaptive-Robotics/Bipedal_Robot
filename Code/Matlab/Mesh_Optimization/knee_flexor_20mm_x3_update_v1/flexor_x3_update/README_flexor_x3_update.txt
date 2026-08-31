FLEXOR 20 MM X3 UPDATE

Run order
1. Put these files with the rest of the existing flexor project files.
2. Run: clear classes; clear functions; rehash
3. Run Opt_sanity.m.
4. Inspect the reported wrap state, route-transition figure, three strain
   curves, off-axis baseline, and collision constraints.
5. Run Opt_run.m only after the sanity output is correct.
6. Run Knee_Flexor_Data_20mm.m after Opt_run saves Bifemsh_20mm_Result.mat.

Route definition
- p1 is in the femur frame.
- pWrap and pEnd are in the t1 frame.
- pWrap x-y starts at [-0.0026, -0.0105] m.
- pWrap z is the midpoint between candidate p1.z and pEnd.z.
- The fully extended sample is always active.
- Beginning with the next sample while flexing:
    aDirect = atan2d(p1T1.y-pEnd.y, p1T1.x-pEnd.x)
    aWrap   = atan2d(pWrap.y-pEnd.y, pWrap.x-pEnd.x)
  pWrap is removed when aDirect > aWrap.
- The route uses no dot product, cross product, or vector-angle helper for
  this release calculation.

Objective changes
- Maintains the pointwise 2 percent z-axis torque target.
- Penalizes pointwise off-axis resultant torque above the x0 baseline.
- Penalizes Contraction/KMAX above 1.
- Requires the wrap to release within the modeled range.

X3
- Uses MonoPamDataExplicit_balanceX3.
- strain_f includes Xi3 and drives force.
- strain_p excludes Xi3 and remains the stretch check.
- Contraction is reported separately.
- The geometric bend measure is Rwrap*abs(aOut-aIn) while pWrap is active,
  with Rwrap = 0.022 m in buildKneeFlexorContext20mm.m.

Existing project dependencies that are still required
- minimizeFlxPin10_results_20260730_2transforms_Z2.mat
- minimizeExtPin10_results_20260819_2transforms_Z2.mat if the flexor result
  does not itself contain Xi3
- RowVecTrans.m, RpToTrans.m, MonoMuscleData.m, and the plotting helpers
  already used by Knee_Flexor_Data_20mm.m
