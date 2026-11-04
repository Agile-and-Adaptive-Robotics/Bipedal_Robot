# Functional routing stage 1 report (2026-11-04)

Machine: easteregg2. Contract of record: `MASTER_DIRECTIVE_20261103.md` (base) as corrected by
`HARDWARE_AND_ZCODE_HANDOFF.md` (10/20 mm sizing, contraction limits, target substitutions,
z=0 centerline, Xi0–Xi3, hardware families). Workflow run: `dwfrun-64c78ae0` (66/66 steps settled, 0 failed).

## Directive-item status (short form)

| Contract area | Status | Evidence |
|---|---|---|
| Geometry extraction + per-pose verification vs OpenSim (Master 1) | DONE | 462 route-pose rows (77 routes × 6 poses), zero nonzero maxErrorMm, owning-body frames (`geometry_20261103/comparison_maxErrorMm.csv`) |
| Vasti / TFL / foot / erector-spinae attention items (Master 1) | DONE | `geometry_20261103/attention_items.json` |
| Conditional/moving points not frozen (Master 1) | DONE | `fr_evalRoute.m:29-39`; med_gas P2 and psoas P3 flip active/inactive in-gate (`output/gate_r4.log:52-54`) |
| Transmission closure / force balance / work; no fictitious joint (Master 2) | DONE | identities at 1e-16 level (`mechanics_20261103/mech4_transmission_20261103.log`; `output/gate_r4.log:30-42`) |
| nBPA / series / sheaves / reeving kept distinct (Master 2) | DONE | `fr_mechBuild.m:55-60`; two-tackle rejected by class error |
| 20 mm diameter endpoints 20→38; 10 mm own law (Master 2, Handoff 2) | DONE, with a convention caveat (see findings) | `fr_diameter.m:17-37`; `mech2_diameters` log |
| Contraction limits 16.5%/24% wired; no universal 25% (Handoff 2) | DONE in campaign code; legacy drivers still carry 25% (see findings) | `fr_setup.m:53-54`, `fr_mechBuild.m:53`; grep: zero 25% hits in campaign code |
| Target corrections vas_int / semimem / add_mag3 / glut_max3 / psoas-iliacus (Handoff 3–7) | DONE | `targets_20261103/target_contract.json` (27 groups, both count conventions recorded) |
| z=0 centerline, sideSign rule (Handoff 8) | DONE | `objective_spec.json` hard_constraints.centerline_z0 |
| Xi0–Xi3 only, adopted rows 77/32 (Master 4, Handoff 9) | DONE | mechanics mech6 Xi log; Xi4 absent from campaign code |
| Objective structure: no off-axis penalty, soft overlap, separate shape term (Master 8) | DONE, weights PLACEHOLDER-PENDING-BEN | `campaign_20261103/objective_spec.json` |
| Layered deep-outward packaging + refinement (Master 5, Ben's ruling) | DONE and PROVEN in-gate | layering proof below |
| Foot/knee tendon spans (Master 6) | DONE as computed constraints | `objective_spec.json` required_tendon_spans |
| Sensitivity gate (Master 9) | DONE, 3 fresh-process rounds exit 0 | `output/gate_r3.log`, `gate_r4.log` (gate of record) |
| Baseline reproduction (prompt "reproduce one existing result") | DONE | tib_ant_r M1 = 100.0000 bit-identical to the 2026-10-23 record, 33 s, scratch copy, originals untouched |
| Full 27-group campaign | LAUNCHED, RUNNING | `campaign_20261103/run_0/` |

Campaign-stage items deliberately NOT-DONE (staged): signed per-pose human target torque curves;
exact per-pose Xi3 bend-loss coupling inside BiPam_pulley2 (currently max-envelope through the Xi0
channel, profiles stored); CAD-level SolidWorks obstacle registration (documented hook
`fr_obstacles:cadHook`); variable-row-count mechanics grids (conditional rows documented fixed);
final tackle-family load ratings; full packaging refinement across all layers.

## Gate of record (output/gate_r4.log, round 4; rounds 1–3 also exit 0)

- Per-variable sensitivity: every attachment/routing variable perturbation changes geometry AND
  objective/constraint outputs; bounds have nonzero extent; seed-restore check passes.
- Representative cases (hip flex/ext, hip rot, flex+ext-rot, knee/Bowden, foot span, trunk erector
  spinae): rigid and compliant both solve with viol 0 (iliacus 4.1e-2, intobl 1.2e-2 worst).
- Rigid-limit recovery: |dtau| ~ 1e-5 to 1e-6 N·m.
- Both bracket contributions: bracket deflection up to 38.26 mm, tendon stretch gama, dstrain vs rigid 1.6e-2.
- Parallel load sharing exact; transmission work consistency 3e-16; compliant closure 1e-15.
- Diameter endpoints: 20 mm → 20.000/38.000 mm exact; 10 mm own law 13.600 → 17.837 mm.
- Conditional points: med_gas condyle wrap active(-44°)=1 / active(-46°)=0 (0.99 mm route jump,
  quantified); psoas P3 flips at +44/+46°.
- **Layering proof (Ben's deep-outward ruling, demonstrated):** superficial route seeded inside the
  retained deep corridor violated the envelope by 6.77 mm; constrained solve deflected 61.11 mm to
  feasibility; refinement (budget 300 vs 150) improved the superficial objective 10135 → 1135 while
  the deep solution degraded by exactly 0.000 (tolerance 1e-6).
- Timing: median 0.128 s / mean 0.177 s per evaluation; full campaign 40,500 evals ≈ 2 h serial.

## Baseline reproduction of record

tib_ant_r capability_v3 M1 = 100.0000 (meetCells 25/25), bit-identical to the 2026-10-23 record
(rest 0.268868 m, tendon 0.0890625 m, ladder (1,1)[14/20]→(2,1)[25/25]→(3,1)[25/25] CONV, 180 evals,
33 s), run from a scratch copy under `tmp/functional_routing_20261103/mech_repro/` (only mechanical
edit: resultsRoot relocated so no new mats land in the repo). Caveat of record: the survey2 machinery
sizes with the pre-handoff kmax = 0.75·rest; this reproduction validates the recorded number, not the
corrected contraction limits.

## Launch (run_0_weights_placeholder)

- Escape criterion written BEFORE launch (`run_0/LAUNCH_DECISION.md`): ≤36 h wall → full branch.
  Arithmetic: 40,500 evals × 0.177 s ≈ 2 h serial; ×10 conservative envelope-inflation bound +
  startup ≈ 20.1 h; pre-flight (`fr_preflight27`, all 27 groups built + seed-evaluated on both
  branches) refined to ~4 h direct / ~8 h worst bound → FULL-27 branch.
- Driver: one detached MATLAB process, parpool(10), layers SEQUENTIAL in contract layerOrder
  (trunk 5 / thigh 10 / lower leg 12), groups parallel within a layer; rigid then compliant per
  group; refinement round 1 after the first deep-outward pass.
- Live at report time: layer 1 complete (psoas 121.011/102.183, iliacus 165.870/156.704 (viol 4.1e-2),
  ercspn 115.499/124.714, intobl, extobl — rigid/compliant), 15 checkpoints on disk, layer 2 running
  with 5 obstacle envelopes carried.
- Monitor: `run_0/logs/master.log` (tail -f), `run_0/checkpoints/<group>_<branch>.mat`,
  `run_0/STATUS.json` (rewritten at layer/round boundaries), per-group `logs/group_<name>.log`.
- CLEAN STOP: create `run_0/STOP`; the driver finishes the in-flight layer, writes
  STATUS.json stoppedCleanly:true, exits; checkpoints intact. RESUME: delete STOP and relaunch
  `"D:/Program Files/MATLAB/R2025a/bin/matlab.exe" -batch "cd('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027'); addpath('campaign_20261103'); fr_launch_run0"`
  — it re-scans checkpoints, skips completed (group, branch) pairs, restores the envelope-of-record.
  Hard stop = kill the MATLAB PID; no partial checkpoints are ever written.

## Findings requiring Ben's ruling or input

1. **Legacy drivers still carry the 25% contraction limit (verified, high).**
   `Opt_run_BiPulley.m:363` and `Opt_run_BiPulley_survey2_20261023.m:439` size with
   `kmax = 0.75*rest`. The campaign code uses the corrected limits. Fix is one argument
   (Kmax = the class's 7th constructor argument, demonstrated in the mech4 fixtures) before any
   sizing run reuses those drivers.
2. **The 20/38 mm diameter pair is not one physical relation (verified, high).**
   The braid law (`bpaR.m`) gives 23.6 mm outer diameter at rest (nominal 20 mm is the Festo inner
   tube; Dout adds 2×1.8 mm membrane) and 34.8 mm at 24%, reaching 38 mm only at 35.1% contraction.
   The 38.5 mm value (`buildKneeExtContext20mm.m:195-197`) and the 38.1 mm circle
   (`buildGeoExclusion.m:69`) are separate measured constant fully-inflated envelopes. The campaign
   declares linear interpolation between calibrated endpoints (`fr_setup.m:62-66`) and prints the
   discrepancy in every gate run. RULING NEEDED: braid law, or measured inflated envelope, as the
   packaging convention.
3. **Objective weights are PLACEHOLDER-PENDING-BEN** (`objective_spec.json`): capacity deficit 100,
   usable-RoM shortfall 40, torque-shape 0.5, soft BPA–BPA overlap 2.0. Stage-1 seeds still show a
   large meet-or-beat deficit (psoas demo 2–7 N·m vs 34–49 N·m target curve). Your weights trigger a
   checkpoint-based re-solve, not a rebuild.
4. **Compliance stand-ins are 6061-T6 [25 5 40] mm** (PENDING-BEN section in the campaign headers).
5. **10 mm calibration file not located**: the 10 mm endpoints (13.600 → 17.837 mm) were reproduced
   from the repo's own bpaR 10 mm branch. If you have a dedicated 10 mm calibration, point to it.
6. **BPA-count convention conflict, documented not resolved**: contract uses ceil(MIF/Fmax20)
   (vas_int 4, soleus 3, glut_max3 4, add_mag3 2, semimem 2, ercspn 3, others 1) while the located
   `bpa_force_screen.csv` used ceil(MIF/(0.5·Fmax20)) (vas_int 8, soleus 6, glut_max3 8, ...). Both
   recorded per group; counts are initial estimates either way.

## Judged-only findings (evidence cited, not independently confirmed)

- Bone point clouds are sparse surface samples (femur 456 pts, tibia 144, pelvis_R 124, sacrum 191,
  calcaneus 999, toes 957, spine 5433): layering distances are ORDER metrics, not clearance
  tolerances (attachment NN-distance means: pelvis 18.6 mm, femur 26.8 mm, spine 35.8 mm).
- Cloud-to-body-frame registration validated only empirically; no authoritative registration
  document found in the repo.
- Robotbody seeds carry Connor-era frozen conditional points (rect_fem_r, med_gas_r, lat_gas_r,
  vas_int_r, semimem_r, grac_r; iliacus_r/psoas_r dropped theirs; psoas_r starts on torso) —
  preserved exactly and documented so they cannot silently stand in; the optimization stage must
  implement the stock behaviors or explicitly replace them with optimized guides.
- Neither model contains any WrapObject/PathWrap (13 and 92 empty sets per model); functional wrap
  behavior lives in ConditionalPathPoint ranges, all captured.
- knee −120° sits on the model's stored range endpoint (realized −119.99999986°, clamping recorded).
- Left-side conditional points censused, not per-pose verified (campaign optimizes the right side;
  left is its mirror).
- opensim 4.4.1 Python binding workarounds documented in `geometry_20261103/README.md`.
- L0 seeds are sampled lower bounds: fresh product-grid recomputation matched the screen CSV within
  0.05 mm on 20/27 groups, differed up to 20.1 mm (gracilis) on 7; no 10 mm screen outcome changed.
- Hardware ratings remain open per the handoff (sheave capacity convention, cable fatigue/bend,
  termination ratings, housing compression, ferrule/stop geometry); the trunk overload case
  (3×20 mm at G=2 → ~3942 N radial) is FLAGGED, not silently resized.
- layerOrder carries the provisional order + the psoas/iliacus coupled-junction adjustment; the
  geometry layering table was not on easteregg2 when the contract was built — apply bone-distance
  adjustments when it lands.

## File inventory

- Geometry: `geometry_20261103/` (extraction JSONs per model, `comparison_maxErrorMm.csv`,
  `attention_items.json`, `bone_distance_layering.csv`, README with the verification protocol)
- Mechanics: `mechanics_20261103/` (mech1–mech6 scripts and logs: force law, diameters, contraction,
  transmission, Xi, reproduction)
- Targets: `targets_20261103/target_contract.json` + provenance README
- Machinery + gate: `campaign_20261103/` (`fr_*.m`, `objective_spec.json`, `output/gate_r*.log`)
- Campaign run: `campaign_20261103/run_0/` (LAUNCH_DECISION.md, MONITOR.md, manifest.json,
  checkpoints/, logs/, STATUS.json)
- Workflow status table: 55 directive-item rows (artifact "Directive item status")
