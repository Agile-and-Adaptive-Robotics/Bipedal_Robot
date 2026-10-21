# SURFACES INDEX — pulley re-optimization outputs (EDIT_SPEC 20261020 §5)

Machine easteregg2, session 2026-10-20/21. Every number below was produced
this session by the named script; records of record were only READ (the
five `BiPulley_opensim_*_20260926_*.mat`, `Bifemsh_20mm_Result_pulley_
20260926_1652.mat`, `Vas_Pam_20mm_Result_20260925.mat`, `Vas_Pam_20mm_
Result_cam_20261014_0416.mat`, `Bifemsh_20mm_Result_20260921_2037.mat`).
New outputs all carry `_20261020` (gain-2 re-solve mats carry `_gain2_`
plus a 20261021 stamp).

## Bowden or moving force point (state verbatim in prose)

The 2026-09-26 batch used a **MOVING FORCE POINT, not bowden**: every batch
mat quotes `RoutingMode={moving_via, moving_via}` (tendon exit body-fixed
to the proximal body, u_t recomputed per grid cell; `BiPam_pulley.m:90-93`),
and the mono bifemsh mat has `ctx.pulleyRoutingMode=moving_via`. Bowden is
implemented but was **never exercised on record**: `BiPam_pulley.m:268-269`
validates `{'moving_via','bowden'}` with the Shimano-type envelope
(`BowdenBossDia` 0.007 m, `BowdenRunClearance` 0.005 m — present in all mats
as unused constructor defaults per `BiPam_pulley.m:238-241`),
`MonoPam_pulley.m:575` actively branches on bowden, and
`Opt_run_BiPulley.m:473-477` has a bowden-only echo branch that never fired.
Verified executable in the audit session (in-memory, nothing written):
soleus at gain [2 1] gives bowden max|tau| = 13.249 N·m vs moving_via
8.967 N·m — live machinery, zero recorded use. All sweeps and the gain-2
re-solve here stay `moving_via`.

## Files produced (this section)

| artifact | path |
|---|---|
| Sweep mats (5) | `Code/Matlab/Mesh_Optimization/Results/PulleySurfaces_<site>_20261020.mat` (soleus_r, tib_ant_r, bifemsh_r, med_gas_r, rect_fem_r) |
| Sweep CSV | `Code/Matlab/Mesh_Optimization/Results/PulleySurfaces_summary_20261020.csv` (45 rows: 5 sites x gains {1/1, 2/1, 3/1} x pressures {300, 450, 620}) |
| Gain-2 re-solve mats (5) | `Code/Matlab/Mesh_Optimization/Results/BiPulley_opensim_<site>_gain2_20261021_*.mat` |
| Gain-2 CSV | `Code/Matlab/Mesh_Optimization/Results/BiPulley_batch_summary_gain2_20261020.csv` |
| Required-torque curves (new §5 input) | `tmp/paper_bridge/reqtorque_20261020/OpenSim_ReqTorque_<site>_<joint>_20261020.txt` + `reqtorque_all_20261020.csv` |
| Figure 1 | `ProofFinal/Figures/PulleyTorqueSurfaces_20261020.eps` (QC: `figscripts/QC_PulleyTorqueSurfaces_20261020.png`) |
| Figure 2 | `ProofFinal/Figures/PulleyMarginSurfaces_20261020.eps` (QC: `figscripts/QC_PulleyMarginSurfaces_20261020.png`) |
| Figure 3 | `ProofFinal/Figures/RouteVariants_20261020.eps` (QC: `figscripts/QC_RouteVariants_20261020.png`) |
| Scripts | `tmp/paper_bridge/figscripts/pulley_sweeps_20261020.m`, `gen_PulleySurfaces_20261020.m`, `gen_RouteVariants_20261020.m`, `gain2_summary_20261020.m`; `tmp/paper_bridge/gen_required_torque_20261020.py`; driver copy `tmp/paper_bridge/Opt_run_BiPulley_gain2_20261020.m` (log `gain2_rerun_20261020.log`) |

## Method + regression gates (all PASSED, run this session)

Evaluation-only fixed-gain sweeps of the recorded 2026-09-26 batch:
`BiPam_pulley(name, points, cross, 20, T, rest, kmax, tendon, 0.025, p,
5e-3, 5e5, 1e4, 6, struct('nPulleyBPA',1,'tackleLineParts',gpair,
'routingMode','moving_via'))` from each mat's saved `spec`/`points`/`T`/
`finalObj` sizes; gains {[1 1], [2 1], [3 1]} x pressures {300, 450, 620}
kPa; dense grid 25 primary-joint x 5 secondary angles (T rebuilt inline
with the driver's own two-hinge logic; the class requires >= 2 joint slabs
on dim 4, so the secondary-pinned margin slice uses a duplicated theta2 = 0
page). Gates, per site (console of `matlab -batch pulley_sweeps_20261020`):
recorded-grid rebuild at gain 1 / 620 kPa reproduces `finalObj.Torque_ins`
with **max |delta| = 0.000e+00 on all five sites** (bit-exact), and the
rebuilt dense T matches the recorded shared cells with max |dT| = 0.000e+00.
Torque reported = Torque_ins projected on the recorded hinge axis
(off-axis at the peak cells <= ~10 % of the vector norm; e.g. bifemsh
max |tau_ax| 17.18 N·m vs vector-norm 19.78 on the recorded grid).
Margins compare MAGNITUDES (the batch objective itself used elementwise
|·|); RoM width = trapezoid measure of the angle set where margin >= 0 on
the 25-point primary-joint slice (secondary joint pinned at 0 deg).

## Required-torque inputs (the section-5 missing input, produced this session)

Stock human model `Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim`
via opensim 4.4.1 (`D:/Anaconda/envs/myo/python.exe`), `computeMomentArm x
max_isometric_force`, 100 points, other coordinates at default; ranges =
each site's recorded `thetaRangesDeg`:

| site / joint | tau_req range (N·m) | fmax (N) |
|---|---|---|
| soleus_r / ankle [-50, 20] deg | -171.35 … -134.08 | 3549 |
| tib_ant_r / ankle [-50, 20] deg | +24.60 … +39.35 | 905 |
| med_gas_r / knee [-120, 10] deg | -36.61 … -14.23 | 1558 |
| med_gas_r / ankle [-50, 20] deg | -74.13 … -58.56 | 1558 |
| rect_fem_r / hip_flexion [-15, 110] deg | +13.31 … +63.16 | 1169 |
| bifemsh_r / knee (cross-check) | -27.23 … -3.57 | 804 |

All five fmax values are identical to the batch mats' recorded fmax.
bifemsh margins use the RECORDED requirement, the mono mat's
`ctx.humanTorque` = -2.46 … -23.16 N·m (per the auditor's plan); the stock
cross-check above differs because the legacy txt's provenance differs
(2026-10-03 probe). The batch's own stored target (fmax x 5 cm heuristic:
soleus 177.4, tib_ant 45.3, med_gas 77.9, rect_fem 58.5, bifemsh 40.2 N·m)
remains the objective of the re-solve below.

## Quotable numbers per site (fixed-gain sweeps, recorded sizes)

Ankle/knee/hip angle grids; peak |tau| in N·m at 620 kPa unless stated.

**soleus_r (ankle, mono-pulley/uniarticular).** Taut only at 620 kPa
(fully slack at 300/450). G=1/1 peak 25.3; G=2/1 peak 11.9; G=3/1 8.0.
Never meets the human soleus requirement (-134 … -171 N·m; ~6.8x short at
best); RoM 0 deg at every (G, p).

**tib_ant_r (ankle, mono-pulley/uniarticular).** The best site. G=1/1:
RoM 0 deg at 300 kPa, 14.6 deg at 450 kPa, **48.1 deg (46.7 contiguous) at
620 kPa**, peak 52.0, min margin over the meeting span -45.6 %. G=2/1 at
620 kPa: peak 26.5, RoM 4.4 deg. G=3/1 at 620 kPa: peak 18.0, RoM 0.
Pressure quote: at G=1/1 the target-meeting RoM grows 0 -> 14.6 -> 48.1 deg
from 300 -> 450 -> 620 kPa.

**bifemsh_r (knee, mono-pulley/uniarticular).** G=1/1: 8 of 25 grid cells
infeasible at every pressure; peak 25.5 at 620 kPa vs the recorded human
requirement peak 23.16 N·m — meets the human target only at isolated
points (RoM 5.4 deg total-width, 0 contiguous). **G=2/1 and G=3/1: ZERO
infeasible cells — the full 130 deg knee span is feasible** (the doctrine's
range gain), at the cost of peak torque 25.5 -> 15.7 (G=2) / 9.8 N·m (G=3).

**med_gas_r (knee+ankle, bi-pulley/biarticular, duplicate crossing [3 3]).**
The recorded 5x5 grid is dead at gain 1 (10 infeasible + 15 slack cells —
reproduced exactly). On the ankle = 0 deg slice a taut window exists for
knee in [-44, -28] deg: G=1/1 peaks 199.6 / 243.6 / 402.5 N·m at 300 / 450
/ 620 kPa with RoM 10.8 / 16.2 / 21.7 deg (7 infeasible cells remain,
knee >= -22.5 deg). At G >= 2 with the recorded gain-1 sizes the route goes
slack everywhere (peak 0). The gain-2 RE-SIZE (below) revives it.

**rect_fem_r (hip crossing of the recorded route, bi-pulley classified).**
Recorded route = femur-ending, crossing 2 rigid, hip live only; margins are
over HIP angle. G=1/1: peak 91.6 / 128.8 / 192.3 at 300 / 450 / 620 kPa,
RoM 15.6 / 31.2 / 46.9 deg (41.7 contiguous). G=2/1: peak 119.7 at 620 kPa,
RoM 13.0 deg. G=3/1: at most 2.6 deg.

**Doctrine quote (torque <-> range exchange, 620 kPa):** raising the
tackle from 1/1 to 2/1 at fixed sizes removes every infeasible cell on the
knee sites (bifemsh 8/25 -> 0/25 over the full 130 deg span; med_gas 7/25
infeasible on its slice at 1/1) while cutting peak insertion torque by
~38-50 % (bifemsh 25.5 -> 15.7 N·m; tib_ant 52.0 -> 26.5 N·m); pressure
buys both (tib_ant RoM 0 -> 14.6 -> 48.1 deg from 300 -> 450 -> 620 kPa).
With the gain-1-optimal sizes, the 2/1 tackle cannot re-tension med_gas
(all-slack), which is what the re-solve fixes.

## Gain-2 full re-optimization (ran; every solve < 10 min)

Driver copy `tmp/paper_bridge/Opt_run_BiPulley_gain2_20261020.m` (5
marked edits; the driver of record untouched) with
`biPulleySpecsFromOpenSim([], [], struct('tackleLineParts',[2 1],
'routingMode','moving_via'))`, FULL budgets (surrogateopt 1000 +
patternsearch 5000, parallel). Wall clock per muscle from the CSV stamps:
med_gas ~7 min including parpool start, soleus 113 s, tib_ant 102 s,
rect_fem 99 s, bifemsh 90 s — all under the 10-minute rule, so the full
re-run is INCLUDED (fullRerunPending = false). Objective = the batch's own
fmax x 5 cm heuristic target (not the human curves).

| site | rest / tendon (m) | margin vs heuristic | peak abs(tau_ax) (N·m) | infeas / slack |
|---|---|---|---|---|
| med_gas_r | 0.5131 / 0.0200 | -1 (exitflagP -2, worstC +0.43) | **184.6** | **0** / 10 |
| soleus_r | 0.2117 / 0.0220 | -0.9982 | 12.5 | 0 / 0 |
| tib_ant_r | 0.3310 / 0.0293 | -0.9645 | 25.4 | 0 / 0 |
| rect_fem_r | 0.4365 / 0.0200 | -0.9957 | 61.3 | 0 / 0 |
| bifemsh_r | 0.0599 / 0.0200 | -1 (exitflagP -2, worstC +0.32) | 12.4 | **0** / 10 |

Re-sizing at gain 2 vs the fixed-gain [2 1] sweep at 620 kPa: med_gas
0 -> 184.6 N·m peak with zero infeasible cells (vs 10 infeasible in the
recorded gain-1 run); bifemsh keeps 0 infeasible (recorded gain-1: 10) at
a lower peak (15.7 -> 12.4 N·m); rect_fem's worst-cell-oriented re-solve
trades peak down (119.7 -> 61.3 N·m); soleus/tib_ant change little. Against
the batch's own heuristic target every margin stays ~ -1, confirming the
09-26 screening conclusion (single 20 mm BPA per reeved route cannot carry
the fmax x 5 cm target) — the human-target margins in the sweep section are
the operative comparison for the doctrine.

## Route-variant figure (Table 5 as curves; R3's replacement input)

`RouteVariants_20261020.eps`, two panels, curves verified from the mats
this session: (A) extensor — direct route `Vas_Pam_20mm_Result_20260925`
TorqueZ 35.45-109.36 N·m; cam design of record `..._cam_20261014_0416`
TorqueZ 40.20-105.91 N·m; vasti requirement ctx.humanTorque 20.25-50.32
N·m + 5 % envelope. (B) flexor — direct `Bifemsh_20mm_Result_20260921_2037`
TorqueZ -81.23 … -3.86 N·m; mono-pulley insertion-side
`Bifemsh_20mm_Result_pulley_20260926_1652` TorqueInsZ -85.56 … -4.02 N·m
(off-axis max 0.803 N·m); bifemsh requirement ctx.humanTorque -2.46 …
-23.16 N·m + 5 % envelope. All four curve ranges match the Table 5
(tab:Redesign) values already in prose — no number changes.

## Disclosures and flags (for the prose agent; do not silently absorb)

1. **rect_fem route choice:** the RECORDED 4-point femur-ending route
   (crossing 2 rigid, hip live) is used everywhere here. The current
   `gait2392_robotbody.osim` + builder reproduce the SAME 4-point spec
   (live batch plan of the gain-2 run: pts=4, cross=[2 3], frames
   pelvis/femur/femur-rigid), so the audit's 5-point live-knee expectation
   does not reproduce on today's file; there was no route divergence to
   choose around. rect_fem margins are over HIP angle only.
2. **med_gas slice:** margins evaluated at ankle = 0 deg (the recorded
   5-angle ankle grid is entirely dead at gain 1); the taut window is knee
   in [-44, -28] deg and collapses to all-slack at G >= 2 with gain-1 sizes.
3. **FLAG for Ben (batch builder, not touched here):** the 09-26 builder
   applies the subtalar fold only when the middle body is not calcn
   (`biPulleySpecsFromOpenSim.m:465-467`), so the tibia->calcn sites
   (soleus_r, tib_ant_r) carry CALCANEUS-body coordinates in an
   ankle(talus)-centered frame — offset by the talus->calcn vector
   (-0.04877, -0.04195, 0.00792) m. Evidence: osim soleus_r-P2 on calcn_r
   is [0, 0.031, -0.0053]; the class receives it unshifted while the hinge
   pivot is the ankle. Consequences visible this session: the reconstructed
   soleus/tib_ant crossing-1 torques are SAME-SIGNED (+) although the
   robotbody model's own moment arms are anatomically opposite (verified
   with opensim on gait2392_robotbody: soleus -170.9 … -136.4, tib_ant
   +23.8 … +40.3 N·m), and soleus's moment arm collapses (peak 25 N·m
   where the osim route implies ~0.05 m of arm). med_gas (femur->calcn) got
   the fold and keeps a large arm. The recorded mats ARE this geometry (the
   sweeps reproduce them bit-exactly), so all soleus/tib_ant numbers above
   inherit it; a corrected fold rule would change those two sites.
4. **Sign conventions:** the bi-class insertion torque is positive for
   knee/ankle flexors (opposite to the mono-class negative-flexion
   convention); all comparisons here are on magnitudes, matching the batch
   objective's elementwise absolute value.
5. **Fixed-gain caveat:** the sweep's G > 1 curves hold the recorded
   gain-1-optimal [rest, tendon]; they show the transmission trade at
   fixed sizes, not a re-sized design (the gain-2 re-solve covers that).
6. **bifemsh requirement provenance:** recorded mono-ctx humanTorque used
   (per plan); stock-model cross-check differs by up to ~4 N·m mid-range.
7. **Figures:** campaign style per §0 (Arial 12 bold, LineWidth 2, palette
   c1/c3/c5/c7 + gray dashed, legend 8 pt, panel letters (A)-(E) on the
   title line). Margin curves are clipped at +300 % for display; -105 %
   marks infeasible/zero-torque cells.

---
## SUPERSEDED — 2026-10-21 DOCTRINE RE-SOLVES (Ben's gain ruling)

ALL gain-ladder quotables above (1/1, 2/1, 3/1 single-BPA fixed-size
re-evaluations + gain-2 re-solves) are SUPERSEDED and must NOT be quoted.
Direction error: TackleLineParts is TRAVEL gain G (insertion force =
nBPA*F/G, travel = G*BPA contraction); force gain = nBPA/G.

Results of record (PulleyDoctrine_summary_20261021.csv, mats
BiPulley_opensim_<site>_{2BPA_G2,2BPA_G4,4BPA_G4}_20261021_*.mat, driver
copies tmp/paper_bridge/Opt_run_BiPulley_doctrine_*_20261021.m, post script
figscripts/pulley_doctrine_post_20261021.m; regression 0.000e+00 everywhere):

Meeting range vs human requirement, 620 kPa (contig):
- tib_ant (ankle): baseline 48.1 (46.7); 2BPA G=2 56.9 (55.4) peak 53.0;
  4BPA G=4 59.8 (58.3) peak 54.5. Baseline pressure ladder 0/14.6/48.1 at 300/450/620.
- bifemsh (knee): baseline 5.4 isolated + 8/25 infeasible; 2BPA G=2 27.1 (21.7)
  peak 34.1 zero infeasible; 4BPA G=4 56.9 (54.2) peak 25.9 zero infeasible;
  2BPA G=4 (force gain 0.5) 24.4 peak 13.3.
- med_gas (knee): baseline 21.7 (16.2) + 7/25 infeasible; 2BPA G=2 65.0 (59.6);
  2BPA G=4 105.6 (102.9) peak 118.0 ZERO infeasible; 4BPA G=4 119.2 (113.7) peak 236.1.
- rect_fem (hip): baseline 46.9 (41.7) peak 192.3; 2BPA G=2 99.0 (93.8) peak 128.0;
  4BPA G=4 96.4 (93.8) peak 160.3.
- soleus: NEVER meets at any transmission (peak 25.3 vs req 134.08-171.35);
  force shortfall, tackle cannot add force.

Fixed gain-1 sizes at G>=2 go slack (proven, psiLo>0) - sizes MUST be re-solved
per configuration. Paper figures PulleyTorqueSurfaces_20261020.eps (now the
2BPA G=2 re-solved configs) and PulleyMarginSurfaces_20261020.eps (now the
matched-force ladder 1BPA/G1, 2BPA/G2, 4BPA/G4 at 620 kPa) regenerated
2026-10-21 by pulley_doctrine_post_20261021.m; paper prose rewritten same day.
