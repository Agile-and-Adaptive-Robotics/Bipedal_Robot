# Routing27 campaign (2026-10-03)

The corrected 27-actuator BPA routing campaign. The human reference is
always STOCK Gait2392 (muscles on their stock OpenSim routes); the robot
side uses the post-09-26 robotbody routes. The prior attempt's bug
(human muscles routed on robot paths used as the human reference) is
structurally impossible here: targets come only from
`gait2392_simbody.osim` via `gen_gait2392_torque_targets.py`.

## Files

- `actuator_map_27.json` — the documented actuator map: master actuator
  to stock human muscle/group and to robotbody seed route. Every group
  decision is grounded in the master's own MIF values (exact single or
  exact group sums; the two judgment calls, the gluteal group and the
  duplicate Tibialis Posterior rows, are recorded under
  `mif_grounding.judgment_calls`).
- `gen_gait2392_torque_targets.py` — target generator (OpenSim 4.6,
  `D:\Anaconda\envs\opensim\python.exe` on easteregg2). Sweeps each
  actuator's primary DOF over the campaign RoM at 101 poses, computes
  tauA = moment arm x MIF per stock muscle, validates moment arms
  against central-difference dl/dtheta at 8 spot angles per muscle (the
  add_mag3 8/8 pattern; a mismatch above 5 percent fails the run), sums
  groups elementwise, and writes `Human_Torques_27\targets27_*.csv` plus
  `targets27_manifest.json` (model and map hashes, RoM sources, FD
  results). tauB (equilibrium force) is intentionally omitted; the
  add_mag3 verdict section 6 documents why.
- `routingSpecsFromRobotbody.m` — robot-route spec builder (adapted from
  `biPulleySpecsFromOpenSim.m`). Handles the full right chain including
  torso (back joint on lumbar_extension) and toes (mtp), freezes
  MovingPathPoints at their coordinate default, resolves
  ConditionalPathPoints at the default pose, folds psoas's torso rows
  into the pelvis frame through the back joint at default, inserts the
  true talus middle frame on tibia-to-calcn spans so subtalar torque is
  expressible, and takes each actuator's first well-formed seed from the
  map's preference list (logged in `Results\routing27_spec_report.txt`).
- `buildRoutingContext.m` — per-actuator context: transform grid on the
  primary DOF (the other crossing holds its default, matching the human
  target generation), target interpolation (pchip), Xi of record
  (flexor 2brk pick 77: Xi0 3.94 mm, Xi1 3.998e4, Xi2 1.4734e4 N/m),
  per-cell bone clouds in the proximal frame, footprint cap, seeds and
  bounds. All maps are resolved here so the context serializes to
  parallel workers.
- `routingChainLib.m` — hinge/rolling-knee edge transforms, route
  segment mapping, Ericson segment-segment distance, cloud clearance.
- `Opt_run_Routing27.m` — the campaign driver. Modes via
  `ROUTING27_MODE`: `plan` (dry-run), `validate` (actuators 11 vasti
  single, 24 gluteal single, 12 soleus pulley), `full` (all 27 x 5
  configs), `smoke`. Budgets `ROUTING27_SURROGATE`/`ROUTING27_PATTERN`
  (defaults 400/2000 full, 60/200 validate). Checkpoint/resume per
  actuator x config; summary CSV; dated full-workspace mats behind
  `liveRun`; rolling diary log.
- `run_routing27.cmd` — easteregg2 launcher (full path to MATLAB
  R2025a; `where matlab` fails there, the full executable path is
  mandatory).

## Modes (Ben's three actuation modes)

| config     | nBPA | tackle gain G | meaning                              |
|------------|------|---------------|--------------------------------------|
| single     | 1    | 1             | one 20 mm BPA on the route           |
| par2/par3  | 2/3  | 1             | parallel BPAs, force scales with n   |
| pulley2    | 1    | 2             | 2:1 tackle: more RoM, force to the body divided by G, pulley-axis reaction penalized |
| pulley2par2| 2    | 2             | 2-BPA bundle through the 2:1 tackle  |

Best-of selection per actuator is feasibility-first, then objective.

## Objective and constraints

- Meet-or-beat (ported from `objective_KneeExt20mm.m`): worst-angle
  deficit vs 1.05 x the stock group torque dominates (1e5), mean-square
  deficit next (1e3), small overshoot penalty, mild design-change
  penalties keeping routes near Connor's seed.
- Pulley penalty: 5e-2 x (max pulley-axis reaction / (nBPA x single-BPA
  max force))^2, zero when no tackle (Ben: the objective must penalize
  the internal reaction at the pulley axis).
- Hard constraints (nonlcon): strain floor (class convention), zero
  pulley-infeasible cells, bone clearance >= inflated BPA radius + 5 mm
  at every primary-DOF pose, lateral footprint <= cap, spacing >= BPA
  diameter + 5 mm vs every route finalized earlier in the campaign.

## Documented simplifications

- Primary-DOF-only sweep; the secondary crossing holds its default
  (both sides of the comparison use the same convention, so the
  meet-or-beat test is consistent).
- Cross-actuator spacing and footprint are envelope metrics in the
  global pelvis frame (each route swept over its own thinned pose set);
  same-joint pairs are near-exact, cross-joint pairs are conservative.
- Bone clouds are subsampled surfaces; clearance is exact only to the
  sample spacing (the 5 mm margin covers it).
- Design variables are bounded deltas on the origin, insertion, and the
  first two via points (max 30 mm each), plus resting length and tendon
  length; routes stay near the seed design by construction.

## Running on easteregg2

```
D:\Anaconda\envs\opensim\python.exe gen_gait2392_torque_targets.py
set ROUTING27_MODE=validate && run_routing27.cmd validate
run_routing27.cmd full        (detached; see README lines in the cmd)
```

Resume: rerun the same mode; the checkpoint skips done configs. Restart
from scratch: delete `Results\routing27_checkpoint_<mode>.mat`.
