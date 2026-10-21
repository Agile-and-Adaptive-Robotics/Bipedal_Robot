# PULLEY_FACTS — Pulley machinery audit (EDIT_SPEC 20261020 §4)

Read-only audit of the 2026-09-26 bi-pulley batch and the mono-pulley
bifemsh run. Machine: easteregg2 (EASTEREGG2, MATLAB R2025a U1, verified in
session). Every claim below was read from the named file or produced by the
named command in this session. Nothing was edited; no records touched.
Audit scripts: `tmp/paper_bridge/audit_pulley_mats_20261020.m` (full log
`audit_mats_output_20261020.log`), `tmp/paper_bridge/audit_pulley_demo_20261020.m`,
`tmp/paper_bridge/audit_osim_paths_20261020.py`.

---

## 1. Gains of record per site

**CSV of record** `Code/Matlab/Mesh_Optimization/Results/BiPulley_batch_summary.csv`
(6 lines, read directly): the `pulleyGain` column is **`1/1` on all five
rows** — `med_gas_r`, `soleus_r`, `tib_ant_r`, `rect_fem_r`, `bifemsh_r`,
all `source=opensim`. Row verbatim:
`med_gas_r,opensim,0,-2,101,10,-1,1/1,20.000,Inf,20260926_1655:32` (the
other four rows likewise carry `1/1`).

**Mat evidence** (loaded per site in MATLAB; identical in all five
`BiPulley_opensim_*_20260926_*.mat`):

| field (finalObj) | value, every site |
|---|---|
| `NPulleyBPA` | `[1 1]` (1 BPA per crossing) |
| `TackleLineParts` | `[1 1]` (G = 1 at both crossings) |
| `PulleyGain` | `[1 1]` |
| `PulleyActive` | `[false false]` — **no tackle was active anywhere** |
| `RoutingMode` | `{moving_via, moving_via}` |
| `Pressure` | 620 kPa (single value) |

`spec.pulley` (builder input) in every mat: `nPulleyBPA=1`,
`tackleLineParts=1`, `routingMode='moving_via'`. This is the DEFAULT path:
`Opt_run_BiPulley.m:112` `defaultManifest()` calls
`biPulleySpecsFromMasters()` / `biPulleySpecsFromOpenSim()` with no
pulleyConfig (`Opt_run_BiPulley.m:261,266`), and both builders default
`1 / 1 / 'moving_via'` (`biPulleySpecsFromOpenSim.m:376-379`,
`biPulleySpecsFromMasters.m:208-211`). The class default is identical
(`BiPam_pulley.m:229-235`: `[1,1] / [1,1] / gain=TackleLineParts /
{'moving_via','moving_via'}`).

Gain notation: the CSV `1/1` is the PER-CROSSING gain pair `[g1 g2]`
(`Opt_run_BiPulley.m:247-253` `gainString` joins the array with `/`). A
scalar `tackleLineParts>1` broadcast to both crossings is REJECTED
(`BiPam_pulley.m:280-288` `TwoTackles` — hit live in this session's demo),
so a gain-2 solve must be written `[2 1]` or `[1 2]`.

Mono-pulley bifemsh run (`Bifemsh_20mm_Result_pulley_20260926_1652.mat`,
product of `Opt_run_pulley.m` at 16:52, 3 min before the batch):
`ctx.pulleyBPACount=2` (the rig's two parallel BPAs), `ctx.tackleLineParts=1`,
`predBest.pulleyConfig` = `nPulleyBPA=1, tackleLineParts=1, gain=1,
routingMode=moving_via`, `predBest.pulleyGain=1`.

**Verdict: every recorded pulley run of 2026-09-26 — batch and mono — ran
at gain 1 (pulley disabled / straight tendon). No gain > 1 was ever solved
or evaluated on record.**

## 2. Bowden or moving force point?

**The batch used a MOVING FORCE POINT: `routingMode = 'moving_via'` in every
recorded mat** (quoted above). In that mode the tendon exit is body-fixed to
the proximal body and the tendon direction `u_t` is recomputed per grid cell
as the joint rotates (`BiPam_pulley.m:90-93`); the moving exit point can
increase the moment arm with rotation.

**Bowden exists and works but was NEVER exercised on record.**
- Supported: `BiPam_pulley.m:268-269` validates `{'moving_via','bowden'}`;
  header `:93-99` defines it (housing anchor on the distal body, `u_t` FIXED
  in the child frame; Shimano-type envelope M7×1.0 mm boss, 5 mm housing +
  5 mm bracket plate). `MonoPam_pulley.m:575` actively branches on
  `strcmp(obj.RoutingMode,'bowden')`. `Opt_run_BiPulley.m:473-477` prints
  the Bowden envelope only when a routing mode equals `'bowden'` — that
  branch never fired in the batch (every mat is moving_via).
- Never exercised: all five batch mats carry `RoutingMode={moving_via,
  moving_via}` and `BowdenBossDia=0.007 m / BowdenRunClearance=0.005 m`
  (constructor defaults, `BiPam_pulley.m:238-241`); the mono mat carries
  `ctx.pulleyRoutingMode='moving_via'`. The CSV has no routing-mode column.
- Verified executable in this session (in-memory, nothing saved):
  reconstructing soleus with `'bowden'` runs and gives DIFFERENT torque than
  `moving_via` at the same gain (at gain [2 1]: bowden max|τ| = 13.249 N·m
  vs moving_via 8.967 N·m), so the mode is live machinery, not dead code.

## 3. BFSH verdict

**bifemsh_r is UNIARTICULAR (knee flexor only) — confirmed from both the
model file and the recorded spec.**
- `gait2392_robotbody.osim` (extracted with `audit_osim_paths_20261020.py`):
  `bifemsh_r` PathPoints = `[femur_r, tibia_r, tibia_r]` — origin on the
  femur, insertion on the tibia, NO pelvis point → crosses the knee joint
  only (the biarticular hamstring is the long head; `bifeml_r` does not
  exist under that name in this model).
- Batch mat `BiPulley_opensim_bifemsh_r_20260926_1702.mat`: `spec.frames =
  {femur_r, tibia_r, tibia_r (rigid)}`, `cross=[2 3]`, `kin.types =
  {hinge, fixed}`, `kneeSlot=1` ("crossing 1 = knee_r hinge … crossing 2
  rigid"), `thetaRangesDeg=[-120 10; 0 0]` — one live crossing, the knee.

**What actually happened in the batch** (row 6 of the CSV, quoted from the
mat): at gain `[1 1]` the 2-variable `[rest; tendon]` search over the short
femur→tibia route (RestingL 0.0487 m) could not meet the heuristic target:
`exitflagS=0, exitflagP=-2`, **10 of 25 grid cells PulleyInfeasible, 10 of
25 slack**, |Torque_ins| ≤ 17.18 N·m vs `tauTarget = 40.2 N·m`
(fmax 804 N × 5 cm nominal arm), margin −1 (zero torque at the worst
cells), penalized objective 101 (= −(−1) + 10×10). It is a gain-1 SIZING
failure on a short uniarticular route — **not a biarticular-routing
infeasibility, and not evidence that "bifem short head failed" as a muscle
class** (R4's deletion of that narrative is correct).

The mono-pulley world already covers it: `Bifemsh_20mm_Result_pulley_
20260926_1652.mat` (Opt_run_pulley, 2-BPA rig, gain 1, 620 kPa) has
`infeasibleCount=0, slackCount=0` and a full torque curve
`TorqueInsZ(φ)` = −85.56…−4.023 N·m over the 100-point knee grid
φ ∈ [−120°, +10°]. BFSH belongs to the mono-pulley / Opt_run
redesign lineage (its historical seeds `Bifemsh_20mm_Result_2026*.mat` are
`Opt_run.m` flexor products; `Opt_run_Ext.m` is the extensor sibling of
that same 8-variable lineage).

## 4. Surface-data availability

**What exists in the mats:**
- Five batch mats: per site, `Torque_ins` is `[5×3×5×2]` — τ(θ1, θ2) on a
  5×5 angle grid per crossing — at ONE pressure (620 kPa), ONE gain
  ([1 1]), at the recorded `[rest, tendon]`. Grids quoted: med_gas θ1 =
  knee [−120…10]°, θ2 = ankle [−50…20]°; soleus/tib_ant θ1 = ankle
  [−50…20]°; rect_fem θ1 = hip [−15…110]°; bifemsh θ1 = knee. Full class
  grids (Ftendon, PulleyTravel, strain, PulleyInfeasible/PulleySlack) and
  the transforms T are saved, so re-evaluation at other (G, p) is direct.
  Peak |τ| per site (active crossing): soleus 27.2 N·m, tib_ant 54.4 N·m,
  rect_fem 139 N·m, bifemsh 17.2 N·m, med_gas 0 N·m (10 infeasible + 15
  slack cells).
- Mono bifemsh mat: `predBest.TorqueInsZ` [100×1] τ(φ) on the 100-point
  knee grid at 620 kPa, gain 1; AND the human reference curve
  `ctx.humanTorque` (−23.16…−2.46 N·m on `ctx.humanAngleD` [−120…10]°) —
  a margin curve τ_robot(φ)/τ_human(φ) is computable from this one mat at
  that single pressure/gain. `ctx.Pbins = 620` (single bin).

**Missing for (a) torque surfaces τ(θ, p):** the pressure axis — no
recorded pulley run stored more than one pressure (all 620 kPa); the batch
angle grids are only 5×5. Missing for (b) margin surfaces over (θ, gain):
the gain axis — no G>1 solve or evaluation exists anywhere on record (CSV
all 1/1; all mats PulleyGain [1 1], PulleyActive false); and required-
torque curves τ_req(θ) per joint exist only in the mono flexor context
(humanTorque), not for the batch sites, whose stored target is the
fmax×5 cm heuristic (`Opt_run_BiPulley.m:355`), not an OpenSim curve.

**Route discrepancy flag (rect_fem_r):** the recorded mat has a 4-point
route ending on the FEMUR with a rigid crossing 2 (frames
{pelvis, femur_r, femur_r (rigid)}, note "slot 3 = rigid clone"), but the
osim file — in git since 397a990a (2026-09-08, before the batch) and
unchanged since c2d95c08 (2026-09-27, which moved only an add_mag3 point)
— carries a 5-point route `[pelvis, femur_r, femur_r, femur_r, tibia_r]`
with a LIVE knee crossing. The other four mats match today's osim exactly.
So the 09-26 batch consumed a rect_fem route variant that no longer matches
the file; any re-run must pick one and disclose it.

## 5. Re-run plan (minimal; records of record untouched; outputs `*_20261020`)

**Recommended: fixed-gain EVALUATION sweeps (no optimizer).** New script
(e.g. `tmp/paper_bridge/figscripts/` → outputs to `Mesh_Optimization/Results/`
with `_20261020` suffix). Per batch mat: load `spec`/`points`/`T`/`xBest`
(or `finalObj.RestingL/TendonL/Kmax`); loop gain pairs {[1 1], [2 1],
[3 1]} × pressure {300, 450, 620} kPa × `'moving_via'`; construct
`BiPam_pulley(name, points, cross, 20, T, rest, kmax, tendon, 0.025, p,
5e-3, 5e5, 1e4, 6, struct('nPulleyBPA',1,'tackleLineParts',gpair,
'routingMode','moving_via'))`; store τ(θ1, θ2) grids + margins.
Measured in this session: **0.06–0.72 s per construction** on the 5×5 grid
(soleus 0.70 s first call, 0.06–0.18 s warm) → 5 sites × 3 gains ×
3 pressures ≈ **45 builds < 1 min** of compute (plus ~30 s MATLAB start per
invocation). A denser angle grid (e.g. 25×5) needs T rebuilt via the same
two-hinge construction (`Opt_run_BiPulley.m:585-634` logic, inline copy) —
≈25× cells per build, still tens of minutes total.
Margins vs τ_req: use the mono mat's `ctx.humanTorque` for bifemsh; for the
other sites τ_req(θ) must come from the gait2392 moment arms × fmax (the
`Human_Data/OpenSim_*_Results.txt` style) — §5's input to produce; without
it only the fmax×5 cm heuristic margin is available.
RoM improvement per site = angular width where margin ≥ 0 at each (G, p).

**Full re-optimization at gain > 1 (only if §5 opts in):** measured 09-26
wall clock was 100–125 s per muscle at FULL budgets (budgetS=1000,
budgetP=5000, parallel — verified in every mat; CSV stamps
16:55:32 → 17:02:34 across five muscles), so per-solve cost is inside §5's
<10 min rule. Caveat: `Opt_run_BiPulley` hardcodes `defaultManifest()`
(gain-1 defaults) and `BiPulley_checkpoint.mat` would skip done muscles —
a gain>1 re-solve needs a small driver COPY in tmp with the pulleyConfig
injected (`biPulleySpecsFromOpenSim({}, [], struct('tackleLineParts',[2 1]))`),
writing only NEW dated mats. Do not edit the driver of record.

**Site categorization** (spec §4.5): mono-pulley / uniarticular —
`soleus_r` (ankle, tibia→calcn, 1 live crossing), `tib_ant_r` (ankle),
`bifemsh_r` (knee, femur→tibia). Bi-pulley / biarticular — `med_gas_r`
(knee + ankle, duplicate crossing [3 3], both hinges live);
`rect_fem_r` is anatomically biarticular (hip + knee) BUT its RECORDED
batch route has only the hip live (femur-ending, rigid crossing 2) —
today's osim route reaches the tibia, so a live-knee rect_fem is available
from the current builder; the re-optimizer must choose and disclose which
route it uses (see §4 flag).

---

### Evidence appendix (commands run this session)

- `matlab -batch audit_pulley_mats_20261020.m` → `audit_mats_output_20261020.log`
  (whos inventory of all 6 mats; verbatim pulley fields; checkpoint
  `doneKeys` = the five opensim keys; mono ctx/predBest pulleyConfig).
- `matlab -batch audit_pulley_demo_20261020.m` (corrected per-crossing
  torque stats; budgets; timed gain/bowden demo; mono φ grid + humanTorque).
- `D:/Anaconda/envs/myo/python.exe audit_osim_paths_20261020.py` + bounded
  `git show 397a990a:...osim` extraction (PathPoint frames, pre- and
  post-batch).
- `ls "D:/Bipedal humanoid"` → absent on easteregg2, so the masters builder
  was skipped by its try/catch (`Opt_run_BiPulley.m:260-264`) — the batch
  is the 5-muscle OpenSim set only.
- File reads: `BiPulley_batch_summary.csv`, `Opt_run_BiPulley.m`,
  `Opt_run_pulley.m` (header + ctx block), `biPulleySpecsFromOpenSim.m`,
  `biPulleySpecsFromMasters.m`, `BiPam_pulley.m`, `MonoPam_pulley.m`
  (greps + targeted reads; line numbers cited above).
