# Optimizer improvement plan — pulley sweep directive execution

Plan for the next optimizer campaign per `Code/Matlab/Mesh_Optimization/PULLEY_SWEEP_DIRECTIVE_20261021.md` (sections 1-6, read this session). Audit date: 2026-10-22, machine easteregg2. Every "exists today" claim below was read from the cited file:line in this session; every "searched and absent" claim was searched repo-wide with the named command. No Ben file is edited by this plan; every new file carries the `_20261022` suffix and lives in `Code/Matlab/Mesh_Optimization/` or `tmp/paper_bridge/opt_improve_20261022/` (doctrine-copy pattern: `tmp/paper_bridge/Opt_run_BiPulley_doctrine_*_20261021.m` lines 1-9 establish the copy-and-mark-edits convention).

**Gain doctrine of record (directive §1, lines 14-16):** insertion force gain = nBPA/G, travel gain = G, `TackleLineParts` = G. Per-crossing gains in `BiPam_pulley` must be written `[G 1]` or `[1 G]`; both > 1 is rejected (`BiPam_pulley.m:280-288`, `TwoTackles`).

---

## Item 1 — Configuration sweep grid nBPA {1,2,3} × G {1,2,3}, sizes re-solved per configuration

### What exists today
- **Flexor (mono) side:** `Opt_run_pulley.m:145-147` `ctx.pulleyConfigSweep` behind env `OPT_PULLEY_CONFIGS=1`, with candidate lists `ctx.pulleyBPACountList = [1, 2]`, `ctx.tackleLinePartsList = [1, 2, 4]`. The discrete config loop (`Opt_run_pulley.m:219-374`) re-runs the FULL surrogateopt (1000 evals) + patternsearch (15000 evals) per config, so the whole 8-variable design including `rest` (sizes) is re-solved per configuration; one dated mat per config (`Opt_run_pulley.m:352-362`, `Bifemsh_20mm_Result_pulley_nB%d_G%d_<stamp>.mat`).
- **Biarticular side:** `Opt_run_BiPulley.m:303-455` `runMuscle` optimizes `[rest; tendon]` per manifest spec (bounds at `Opt_run_BiPulley.m:358-361`), i.e. sizes re-solved per configuration; grid N1×N2 = 5×5 (`Opt_run_BiPulley.m:114-115`).
- **Doctrine round coverage (2026-10-21):** only (2, [2 1]), (2, [4 1]), (4, [4 1]). Manifest injections at `tmp/paper_bridge/Opt_run_BiPulley_doctrine_2BPA_G2_20261021.m:268-269`, `Opt_run_BiPulley_doctrine_2BPA_G4_20261021.m:258`, `Opt_run_BiPulley_doctrine_4BPA_G4_20261021.m:268` (`nPulleyBPA`, 4, `tackleLineParts` `[4 1]`). Dated mats + CSVs exist in `Results/` (`BiPulley_opensim_*_{2BPA_G2,2BPA_G4,4BPA_G4}_20261021_*.mat`, `BiPulley_batch_summary_doctrine_*_20261021.csv`; verified by directory listing this session).
- **Slack mechanism (why sizes must be re-solved):** `MonoPam_pulley.m:63-67` documents that required contraction `(DeltaL + delta_t)/G` DECREASES with G, so a direct-drive-sized actuator runs slack; the solve clamps slack cells to s = 0 (`MonoPam_pulley.m:864-875`), `Ftendon(slack) = 0` (`MonoPam_pulley.m:552-553`), and `F_ins` is zeroed (`MonoPam_pulley.m:595`). The directive records this as proven 2026-10-21 (`PULLEY_SWEEP_DIRECTIVE_20261021.md:25-26`).
- **Sanity gates:** `Opt_sanity_pulley.m` / `Opt_sanity_BiPulley.m` assert the G = 1 regression identity (1e-8) and G = 2 closure (1e-6); the smoke self-check is embedded at `Opt_run_pulley.m:1669-1781` (`smokeSanityCheck`).

### What is missing
1. nBPA = 3 has never been run anywhere (flexor lists cap at 2; doctrine used 2 and 4).
2. G = 3 has never been run anywhere (lists use 1, 2, 4).
3. The directive's 3×3 grid {1,2,3}×{1,2,3} is not covered by either driver's lists.
4. Flexor pipeline hard limit: `predictKneeFlexorPulley` errors for BPAcount other than 1 or 2 (`Opt_run_pulley.m:1395-1397`); nBPA = 3 on the flexor side must ride the `BPAcount == 1` + `nPulleyBPA = 3` fold-into-transmission path (`Opt_run_pulley.m:1344-1355`) with the symmetric-bundle convention from `biPulleySpecsFromOpenSim.m:381-422`.
5. The biarticular batch torque target is the heuristic `fmax * 5 cm` arm (`Opt_run_BiPulley.m:353-358`), not per-joint human targets (the header says so); per-site human targets exist only in Routing27_v2's generated CSVs (`Routing27_v2/buildRoutingContext.m:131-140`).
6. `TackleLineParts` validation accepts any value ≥ 1, integer not enforced (`MonoPam_pulley.m:240-247`; `BiPam_pulley.m:252-260`) — harmless for {1,2,3} but worth an assert in the new drivers.

### New-file design (never edits Ben's files)
- `tmp/paper_bridge/opt_improve_20261022/Opt_run_pulley_sweep33_20261022.m` — verbatim copy of `Opt_run_pulley.m` with marked edits (doctrine pattern): (a) `ctx.pulleyBPACountList = [1 2 3]`, `ctx.tackleLinePartsList = [1 2 3]`; (b) nBPA = 3 flexor runs force `ctx.BPAcount = 1` + `ctx.pulleyBPACount = 3` (transmission-fold path, with a printed banner explaining the single-corridor modeling); (c) new dated-mat prefix `Bifemsh_20mm_Result_pulley33_nB%d_G%d_<stamp>.mat`; (d) assert integrality of both lists; (e) sweep forced on (`ctx.pulleyConfigSweep = true`) with the smoke path intact.
- `tmp/paper_bridge/opt_improve_20261022/Opt_run_BiPulley_sweep33_20261022.m` — copy of `Opt_run_BiPulley.m` with an outer `for nB = [1 2 3], for G = [1 2 3]` loop that injects `struct('nPulleyBPA', nB, 'tackleLineParts', [G 1], 'routingMode', 'moving_via')` into the single-builder manifest (the exact doctrine-copy edit, generalized), new checkpoint `tmp/paper_bridge/opt_improve_20261022/BiPulley_checkpoint_sweep33_20261022.mat`, new CSV `Results/BiPulley_batch_summary_sweep33_20261022.csv` with two added columns `nBPA`, `G`, and dated mats `BiPulley_<source>_<name>_nB%d_G%d_<stamp>.mat`. One copy instead of nine keeps the grid in one log; per-config chapters mimic the doctrine structure.
- Optional follow-on (flagged, not blocking): retarget the biarticular objective from the 5 cm heuristic arm to the Routing27 per-site CSVs — separate copy, only if Ben asks.

### Verification gate
- Pre-flight: `OPT_PULLEY_SMOKE=1` and `OPT_BIPULLEY_SMOKE=1` passes on the new copies (the G = 1 identity at 1e-8 and G = 2 closure at 1e-6, `Opt_run_pulley.m:1753-1779`).
- Real run (easteregg2, `parpool` ≤ 10, budgets 1000/15000 flexor and 1000/5000 biarticular, launched in background with a log per the ~8 min rule): all 9 configs complete; every winner mat records `sum(PulleySlack(:)) == 0` over the task envelope and `worstConstraint <= 1e-6` (the `Opt_run_BiPulley.m:382-384` constraint pair), which is the per-configuration size re-solve proof.
- Regression anchor: the (1,1) and (2,2)-adjacent cells must reproduce the 2026-09-26 G = 1 baseline behavior (`Results/BiPulley_batch_summary.csv`, gains `1/1`, per `tmp/paper_bridge/PULLEY_FACTS.md` section 1) within run-to-run tolerance, and the (2,2)/doctrine cells must agree with `BiPulley_batch_summary_doctrine_2BPA_G2_20261021.csv` in sign and scale.
- Records of record untouched: `Results/BiPulley_batch_summary.csv`, `BiPulley_checkpoint.mat`, all `*_20260926_*`/`*_20261021_*` mats are only ever read. Deliverables per directive §4: margin surfaces over the 3×3 grid, quotable RoM/peak tables, dated mats/CSVs, `tmp/paper_bridge/SURFACES_INDEX.md` updated with the sweep33 artifacts.

---

## Item 2 — Fuzzy escalation: torque-short → +BPA, RoM-short → +tackleLineParts (directive §5, lines 48-57)

### What exists today
- Nothing. The only "selection" machinery is the blind discrete sweep (`Opt_run_pulley.m:189-224`) and the fixed doctrine configs. `Routing27_v2/buildRoutingContext.m` has per-DOF task envelopes but no escalation logic. No code classifies a site as torque-short versus RoM-short.

### What is missing
The per-site, per-pose-region diagnosis that chooses the next rung of the configuration ladder. Inputs it needs all exist as class outputs: insertion torque `Torque_ins` (`MonoPam_pulley.m:606-611`; `BiPam_pulley.m:688-690`), infeasibility/slack flags `PulleyInfeasible`/`PulleySlack` (`MonoPam_pulley.m:622-623`), travel budget terms `deltaL`, `RefIndex`, `PulleyTravel` (`MonoPam_pulley.m:149-153`), and the closure `G*(s - s0ref - cb*F) = DeltaL + nBPA*F/(G*kSpr)` (`MonoPam_pulley.m:850-857`; `BiPam_pulley.m:1014-1022`).

### New-file design
- `tmp/paper_bridge/opt_improve_20261022/escalation_advisor_20261022.m` — pure function `[advice, report] = escalation_advisor_20261022(finalObjOrMat, siteSpec, targetTau)`. Read-only; loads a solved result mat (e.g. a doctrine `*_2BPA_G2_20261021_*.mat`) or accepts a constructed `MonoPam_pulley`/`BiPam_pulley`. Per site and per pose region (grid cell):
  - **Torque-short** at cells with a target: `min(|Torque_ins|/tauTarget - 1) < 0` → advice `nBPA + 1` (insertion force scales nBPA/G, `MonoPam_pulley.m:548`).
  - **RoM-short**: any `PulleyInfeasible` cells, or a travel deficit `max(DeltaL) + nBPA*F/(G*kSpr) > G*(KMAX*Rest - s0ref)` from the closure terms → advice `G + 1` (travel gain = G).
  - Both → both; neither → converged. Emits an ordered ladder per site, printed and saved as `tmp/paper_bridge/opt_improve_20261022/escalation_advice_20261022.csv`, consumed as the PRUNED candidate list for the Item 1 sweep33 drivers (the sweep stays as the completeness check; the advisor orders/prunes, never replaces).
- Anti-gaming: the advisor's report includes both the raw worst margin and the per-region breakdown, so a ladder step cannot be justified by an artifact-corner cell — it consumes Item 6's corner weights as a mask.

### Verification gate
- Replay on the 2026-10-21 doctrine mats must reproduce the known verdicts: e.g. `doctrine_rerun_20261021.log:1160` records `margin = -1.0000` for med_gas_r at (2, [2 1]) → torque-short → advise nBPA 3, matching Ben's ladder.
- Unit polarity test `tmp/paper_bridge/opt_improve_20261022/test_escalation_advisor_20261022.m`: (a) a deliberately direct-drive-sized `rest` evaluated at G = 2 must classify as slack/RoM-side (the proven slack result, Item 1); (b) a scaled-down target-pressure object must classify torque-short; (c) deterministic, zero optimizer calls (runs in seconds, no parpool).

---

## Item 3 — Moving-via deflections from the BPA origin with the tackle-increased internal force (directive §2)

### What exists today
Most of §2 is already implemented in the 2026-09-24 classes:
- Bracket compliance references the BPA origin (row 1, proximal frame): `MonoPam_pulley.m:739-741` (`pA = L(1,:,iRef)`), per-orientation force direction from the origin `computeForceVector` (`MonoPam_pulley.m:645-687`); the biarticular generalizations `computeForceVectorBi` (`BiPam_pulley.m:853-885`) and `LokBi` with `pA = L(1,:)` (`BiPam_pulley.m:911`).
- The equilibrium carries the tackle terms through the solve: `psiFun = G*s - G*cb*F(s) - nBPA*F(s)*cSpr/G - DeltaL - G*s0ref` (`MonoPam_pulley.m:856-857`; per-cell reference `s0(i) - DeltaL(i)` at `BiPam_pulley.m:1014-1020`).
- moving_via honors the deformed exit: `u_t` recomputed per orientation from the DEFORMED exit point (`MonoPam_pulley.m:568-578`; `BiPam_pulley.m:639-665`).
- Load split is carried: BPA-side segments at `nBPA*F_single`, tendon at `nBPA*F/G` (`BiPam_pulley.m:593-601`); mount reaction at the tackle block `R = nBPA*F*u_bpa + F_t*u_t` (`MonoPam_pulley.m:613-620`; `BiPam_pulley.m:692-724`).
- `Opt_run_pulley` documents the real-flexor-route caveat (moving wrap segment inside `MuscleLength` but outside the pulley span, ~3.3 N·m of ~59 N·m deviation, `Opt_run_pulley.m:123-131`).

### What is missing
The tackle-block mount reaction is an OUTPUT only. No compliance anywhere consumes `ReactionF`: only row 1 (the BPA origin bracket, Xi1/Xi2) deforms, and the moving-via pulley block's own mount sees the tackle-amplified internal force (`nBPA*F*u_bpa + F_t*u_t`) with no stiffness model. That is the concrete gap in "the increased internal force accounted for" for the deflection path at the moving via.

### New-file design
- `tmp/paper_bridge/opt_improve_20261022/mountDeflectionAudit_20261022.m` — post-solve function `[delta, report] = mountDeflectionAudit_20261022(obj, Cmount)`. Takes a solved `MonoPam_pulley`/`BiPam_pulley` (or loads a doctrine mat's `finalObj`) plus a 3×3 compliance matrix; extracts `ReactionF`/`ReactionFmag`, `u_bpa`, `u_t` per grid cell; computes the block-mount deflection `delta = Cmount * R` and its first-order effect on the exit point and moment arm (`mA_ins` shift). Reports worst-per-pose values and the induced `u_t` rotation. Additive analyzer: the classes and their G = 1 regression are untouched.
- `tmp/paper_bridge/opt_improve_20261022/pulleyDeflectRescore_20261022.m` — one deflection iteration at the driver level: solve → deflect the exit row by `delta` → rebuild the route (`buildKneeFlexorRoute20mm` for the flexor; manifest points shifted for the biarticular) → resolve once; reports both the undeflected and deflected margins so the sweep33 drivers can quote the compliance-honest margin. (A fixed-point inside the class would require editing Ben's classes; forbidden. One iteration bounds the effect and is quoted as such.)

### Verification gate
- Rigid-consistency: with `Cmount = 0` the audit returns exactly zero and the rescore reproduces the loaded mat's recorded margins bit-for-bit (same cells).
- Force scaling: across the doctrine G2 vs G4 mats at the same site, worst `|delta|` must scale with `ReactionFmag` approximately linearly (print the ratio; large deviation flags a bug).
- Monotonicity: deflection magnitude increases with nBPA at fixed G (bigger bundle pull) and with G at fixed nBPA on the BPA side where `ReactionF` grows.
- Log everything to `tmp/paper_bridge/opt_improve_20261022/mount_deflection_audit_20261022.log`; no writes to any Results mat.

---

## Item 4 — Bowden branch: packaging enforcement + bracket deflection at aluminum/CroMoly stiffness (directive §3, lines 37-46)

### What exists today
- routingMode `'bowden'` with `u_t` pinned at the theta0 cell in the child frame: `MonoPam_pulley.m:575-577`; `BiPam_pulley.m:649-664` (with the documented `sub2ind` stride fix).
- Envelope fields with Shimano-type constructor defaults: `BowdenBossDia = 0.007` m, `BowdenRunClearance = 0.005` m (`MonoPam_pulley.m:229-232`; `BiPam_pulley.m:238-241`), carried in `ctx` knobs (`Opt_run_pulley.m:136-137`, display backfill `Opt_run_pulley.m:456-460`).
- Driver echo only: "verify the housing path against both" (`Opt_run_BiPulley.m:473-477`).

### What is missing (searched)
1. **Enforcement:** the envelope is printed, never checked. No constraint or objective term uses `BowdenBossDia`/`BowdenRunClearance`.
2. **Bowden bracket deflection:** the housing anchor bracket on the distal body has no stiffness term anywhere; the class deforms only row 1 under the BPA-side force (`MonoPam_pulley.m:783-790`).
3. **Materials:** repo-wide grep for `6061|4130|chromoly|cro-moly|alumin` over `*.m`/`*.py` in `Code/` and `tmp/paper_bridge/` (run this session) returns only coincidental digit matches; no bicycle-standard material constants exist. `Xi1`/`Xi2` are fitted effective system stiffnesses (AGENTS.md, Xi1/Xi2 semantics), NOT material-parametrized bracket models.
4. Real hardware dimensions: the directive itself says "revisit with real hardware dimensions" (`PULLEY_SWEEP_DIRECTIVE_20261021.md:41-42`); today's values are Shimano-type stand-ins.

### New-file design
- `tmp/paper_bridge/opt_improve_20261022/bowdenBracketStiffness_20261022.m` — material + section model returning `Cmount` (consumed by Item 3's audit): table with `aluminum6061T6` (E = 68.9 GPa) and `chromoly4130` (E = 205 GPa), cantilever bracket of user-set section (b, h, L): axial `k_a = E*A/L`, bending `k_b = 3*E*I/L^3` per axis; `C = diag([1/k_a, 1/k_b, 1/k_b])`. Section dims are inputs with PENDING-BEN defaults (5 mm plate section guess), never silent.
- `tmp/paper_bridge/opt_improve_20261022/nonlconBowdenPackaging_20261022.m` — composable constraint `c(9)`: worst over poses of `(BowdenBossDia/2 + BowdenRunClearance) - min signed distance from the bowden housing line (exit-to-anchor line at the theta0 direction, per `MonoPam_pulley.m:575-577`) to the tibia/femur SDFs`. Reuses `signedDistanceToTibia20mm` and the femur SDF helpers exactly as `nonlconExclusion.m:278-303` does, by CALLING them, not editing `nonlconExclusion.m`. Appended by the sweep33 drivers' local nonlcon wrappers (the `nonlconExclusionPulley` pattern, `Opt_run_pulley.m:1519-1550`, is the precedent for extending the vector without touching the base file).
- Bowden-mode runs in the sweep33 grid: at least the (1, G3) column with `routingMode = 'bowden'` so the branch is exercised per the directive.

### Verification gate
- Material check: swapping aluminum ↔ CroMoly must change `k` by the E ratio (205/68.9 ≈ 2.97) and worst deflection by its inverse — asserted analytically in the audit log.
- Polarity test: with the envelope radius shrunk well below actual bone clearance the packaging constraint returns ≤ 0; inflated beyond it, > 0 (synthetic route, no optimizer).
- Bowden-mode object check: `all(obj.u_t == repmat(obj.u_t(obj.RefIndex,:), N, 1))` on a constructed bowden object (direct assert on class output, `MonoPam_pulley.m:575-577` behavior); smoke-budget end-to-end run of the (1,3) bowden config completes with `PulleySlack == 0`.
- Real-dimension swap: when Ben supplies hardware dims, only the constants in the two new files change; gate re-runs.

---

## Item 5 — Hip adduction/abduction coupling limits (directive §5, lines 59-64): SEARCH RESULT

### Search performed (this session, repo-wide)
- `grep -rniE "adduct"` over `Code/**.m`; `grep -rniE "(adduct|abduct).{0,40}(limit|cap|bound|clamp|narrow)|(limit|cap|bound|clamp|narrow).{0,40}(adduct|abduct)"` over `*.m *.py *.md *.json` repo-wide; `grep -rniE "coupling.*(hip|adduct)|(hip|adduct).*coupl"`; flexion-dependent-range patterns; `CoordinateCouplerConstraint` in the OpenSim models; plus reading every hit's context.

### Verdict: the coupling limit code is ABSENT from the repo
Nearest artifacts, none of which are it:
1. `Code/Matlab/Mesh_Optimization/Routing27_v2/buildRoutingContext.m:29-31` — header note: "Full hip flexion combined with full abduction is excluded as impractical" (a documented exclusion, not code), and `buildRoutingContext.m:70-83` — per-primary-DOF task-RoM table (`'hip_adduction_r', [-25, 20]`, `'hip_flexion_r', [-45, 120]`) intersected with the campaign range. One DOF swept at a time, other crossing frozen at default (`buildRoutingContext.m:60-66`); no flexion-dependent adduction narrowing. `TASK_ROM_TABLE.md` does not exist yet (`find` returned nothing; the header says "Refined values pending").
2. `Code/Matlab/Mesh_Optimization/Add_Mag_Mesh_Opt.m:29-45` — pure rectangular tensor-product sweep: adduction [-45°, 20°] × rotation {-45, 0, 40} × flexion {-25, 0, 85}, no coupling. Same grid in `Code/Matlab/testBiMeshCalculation.m:27-45` (three-point version).
3. `biPulleySpecsFromOpenSim.m:29` (header note 2) — secondary coordinates (hip adduction/rotation, subtalar, mtp) pinned at ZERO; the bi-pulley builder never sweeps adduction at all.
4. The OpenSim models carry no `CoordinateCouplerConstraint` (grep count 0 in `Solid_Models/OpenSim/Gait2392_Robotbody/*.osim`).
5. `Code/MuJoCo_SNS/spinal/runner.py:514` sets capsule STIFFNESS on `hip_adduction_{side}` (passive stiffness in the MuJoCo SNS stack, not a kinematic limit, and not part of the optimizer).

### New-file design (the directive's "or write it")
- `tmp/paper_bridge/opt_improve_20261022/hipAdductionEnvelope_20261022.m` — pure function `mask = hipAdductionEnvelope_20261022(hipFlexDeg, adductionDeg)` implementing `adducLimDeg = f(hipFlexDeg)` as piecewise-linear interpolation over a control-point table in `tmp/paper_bridge/opt_improve_20261022/hipAdductionEnvelopeTable_20261022.csv` (columns: `hipFlexDeg, adducLimDeg, source`). Control points ship as `PENDING-BEN` placeholders with citation slots (the envelope tightens near full flexion and full extension); the code is complete, the physiology numbers are not invented. Returns a soft weight (1 inside, taper outside — shared taper kernel with Item 6) so it feeds the corner machinery directly.
- Until Ben supplies and approves the control points, every hip-site run in the sweep33 drivers must print `HIP ENVELOPE UNTRUSTED — coupling table pending Ben` and treat hip-site extreme corners with the Item 6 floor weight. The directive's own bar ("before the hip sites' requirement envelopes are trusted", `PULLEY_SWEEP_DIRECTIVE_20261021.md:63-64`) is honored by refusing to trust silently.

### Verification gate
- `tmp/paper_bridge/opt_improve_20261022/test_hipAdductionEnvelope_20261022.m`: (a) at the control points the interpolated limit equals the table exactly; (b) the full-flexion + full-abduction corner is masked (consistent with `buildRoutingContext.m:29-31`); (c) with the PENDING table the untrusted-banner rule fires (asserted by calling the driver's guard helper); (d) zero dependence on any Ben file.

---

## Item 6 — Soft extreme-pose penalties (directive §5, lines 66-68)

### What exists today
- The pattern exists ONLY in Routing27_v2: task-RoM envelope `ctx.isTaskPose`/`envCells`/`outerCells` with hard constraints inside the envelope and soft penalties outside, weights `ctx.romW = 1e3`, `ctx.outerW = 1e2`, `ctx.anatW` (`Routing27_v2/buildRoutingContext.m:79-99`). Single-DOF harness.
- The pulley campaigns have NO corner softening: `Opt_run_BiPulley.m:367-385` takes `margin = min(validTau)/tauTarget - 1` over ALL 5×5 cells; `objective_KneeFlexorPulley` prices the worst-angle shortfall over the whole sweep (`Opt_run_pulley.m:1196-1199`).

### What is missing
Combined-pose (two-joint-grid) corner weighting: the biarticular grid corners (e.g. full knee flexion combined with extreme ankle, or hip flexion × adduction corners once Item 5 lands) are enforced at full weight today, so the solver is over-penalized for extreme combined poses the directive says humans may not need.

### New-file design
- `tmp/paper_bridge/opt_improve_20261022/softCornerWeights_20261022.m` — pure function `W = softCornerWeights_20261022(theta1Deg, theta2Deg, env1, env2, opts)`: 1 inside both per-joint task envelopes, cosine/quadratic taper to a floor `wMin` (default 0.1, PENDING-BEN) across a corner band (default 10°, PENDING-BEN), optionally AND-ed with the Item 5 hip-coupling mask. No side effects; unit-testable.
- Consumption (driver copies only, never Ben's objectives): in `Opt_run_BiPulley_sweep33_20261022.m` the local `batchObjective` becomes `J = -Wweighted-margin + ...` with `shortfallWeighted = W .* shortfall` (min over cells of the weighted shortfall replaces the raw min); in `Opt_run_pulley_sweep33_20261022.m` the same weighting wraps the local `objective_KneeFlexorPulley` call site (the local function is already sweep-owned). Anti-gaming stays intact: mid-envelope cells keep W = 1 exactly, so the softening can never hide a functional-RoM failure; `JoffAxis`, `JwrapRelease`, `JrestLength`, `Jstrain*` terms are untouched (`Opt_run_pulley.m:1206-1242`).
- Reporting rule: both raw and softened worst margins are printed and stored in the CSV so no softened number can be mistaken for a raw one.

### Verification gate
- Identity: with `W = ones` the weighted objective equals the unweighted one bit-for-bit on a loaded doctrine mat's cells (recompute and compare).
- Polarity: a synthetic design failing ONLY in the tapered corner scores strictly better than one failing mid-envelope; a mid-envelope failure is never softened (unit test with hand-built shortfall grids).
- Record: sweep33 CSVs carry `worstMarginRaw` and `worstMarginSoft` columns; the campaign report quotes both.

---

## Item 7 — Gait2392 via-point/wrap-flip path-artifact audit helper (directive §6)

### What exists today
- `Routing27_v2/gen_gait2392_torque_targets.py` validates moment arms against a central difference of muscle length at 8 spot angles per muscle, 5% tolerance (`FD_SPOTS = 8`, `FD_TOL = 0.05`, lines 25-33), and records `opensim_version` + hashes in the manifest (lines 155, 194). That is a numerical cross-check on the targets, not an artifact audit.
- `tmp/paper_bridge/audit_osim_paths_20261020.py:1-16` dumps the 5-muscle list's path points (read-only print).
- **Model facts verified this session** (scripted scan of `gait2392_robotbody.osim` with `D:/Anaconda/envs/myo/python.exe`): the model carries 52 `ConditionalPathPoint` elements; `vas_med_r` and `vas_lat_r` each carry **4 ConditionalPathPoints + 2 MovingPathPoints** (the vasti are the extensor requirement source, `Testing_Data/2022_02_Festo/OpenSim_Vasti_Results.txt`, loaded by `buildKneeExtContext20mm.m:120`); `semiten_r` carries 2 ConditionalPathPoints; **no muscle carries an actual PathWrap object** (all 92 `PathWrapSet`s are empty).
- **Silent-drop bug class found:** `biPulleySpecsFromOpenSim.osimPathPoints` reads only exact-tag `PathPoint` (`biPulleySpecsFromOpenSim.m:193-212`), so ConditionalPathPoint/MovingPathPoint rows are silently omitted from the MATLAB routes — the doctrine log shows `semiten_r` would enter at pts=4 while the model carries 6 path rows (the current 5-muscle defaults: `med_gas_r`, `soleus_r`, `tib_ant_r`, `bifemsh_r` carry none; `rect_fem_r` carries none). The `insertWrapMidpoints` refusal (`biPulleySpecsFromOpenSim.m:215-235`) can never fire on this model (no PathWrap objects exist).

### What is missing
The per-site audit the directive orders: visually/numerically verify each muscle's path across the full requirement range — no via-point activation jumps (moment-arm steps at the graded angles), no wrap-side flips, no bone penetration — and flag sites whose requirement curve inherits an artifact. Nothing in the repo detects conditional-point activation or arm discontinuities.

### New-file design
- `tmp/paper_bridge/opt_improve_20261022/audit_osim_path_artifacts_20261022.py` — read-only audit, pinned to opensim 4.4.1 (`D:\Anaconda\envs\myo\python.exe`), on `gait2392_simbody.osim` (the model `gen_gait2392_torque_targets.py:44-46` used) with the robotbody twin scanned for path-point inventory. Per site muscle and per DOF: sweep finely (0.5°), record moment arm and path length vs angle from the opensim API; flag (1) any angle where the active path-point set changes (conditional activation/deactivation) or the arm steps more than 2% of peak between adjacent samples; (2) wrap-side flips (API wrap contact changes; expected zero on this model — the check stays for the simbody twin and future wrapped models); (3) bone-penetration proxy via point-set changes plus a TODO-marked optional STL proximity check against the MyoConverter STLs (never silent). Emit `tmp/paper_bridge/opt_improve_20261022/path_artifact_audit_<model>_<stamp>.csv` with `opensim_version` in the header (mirroring `gen_gait2392_torque_targets.py:194`) and a per-site `artifactFlag` column.
- Consumers: the Item 6 corner weights down-weight flagged activation-adjacent angle bands, and the Item 2 advisor ignores shortfall cells inside flagged bands when advising ladder steps; the sweep33 campaign report carries the flag table so no flagged site's envelope is quoted without the caveat (directive §6's "flag any site whose requirement curve inherits an artifact").
- The `vas_med_r`/`vas_lat_r` finding above also needs a one-line caveat recorded against `OpenSim_Vasti_Results.txt` consumers in the campaign report (the extensor requirement of record inherits conditional-point paths by construction).

### Verification gate
- The audit MUST flag `vas_med_r`, `vas_lat_r` (4 conditional + 2 moving points each), and `semiten_r` (verified present this session), and MUST report zero flags for clean two-point muscles (e.g. `glut_med1_r`, staticPP=2).
- Version pin asserted: recorded `opensim_version` must equal the pinned 4.4.1 env string; the run reads only (no `.osim` writes — verified by checking file mtimes before/after in the gate log).
- Gate command of record: `"D:/Anaconda/envs/myo/python.exe" tmp/paper_bridge/opt_improve_20261022/audit_osim_path_artifacts_20261022.py --model gait2392_simbody --out tmp/paper_bridge/opt_improve_20261022/` run once before any sweep33 real budget is launched.

---

## Execution order and campaign rules

1. **Order:** Item 7 audit → Item 5 + 6 masks (they gate envelope trust) → Item 2 advisor replay on doctrine mats → Items 3/4 compliance analyzers → Item 1 sweep33 real runs (easteregg2, background + log, ~8 min class, `parpool` ≤ 10) → deliverables per directive §4 (surfaces, margin tables, dated mats/CSVs, `SURFACES_INDEX.md` update).
2. **Never edit in place (records of record):** `Opt_run*.m`, `BiPam_pulley.m`, `MonoPam_pulley.m`, `biPulleySpecsFromOpenSim.m`, `nonlconExclusion.m`, `objective_/nonlcon_/predict*` evaluators, all recorded mats/CSVs in `Results/` and `Testing_Data/`. All new code is new files, suffix `_20261022`, in `Mesh_Optimization/` or `tmp/paper_bridge/opt_improve_20261022/`.
3. **Real budgets only** for the campaign runs (no token runs); smoke gates first (`OPT_PULLEY_SMOKE=1`, `OPT_BIPULLEY_SMOKE=1`); anti-gaming objective terms stay in.
4. **No git commits**; the working tree stays clean for Ben, with a commit list at the end.
5. PENDING-BEN values (hip coupling control points, corner band and floor, bracket section dims, real bowden hardware dims) are carried as named constants in the new files and printed in every run's preamble — never guessed silently. If any PENDING value blocks a directive gate, escalate rather than substitute.
