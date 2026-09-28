# add_mag3_r P2 verdict — performance ruling (2026-09-26)

**WINNER: the full-precision thumb-drive point.** `add_mag3_r-P2` (pelvis frame)
= **(-0.1357, -0.0929, +0.0591) m** — APPLIED and regenerated. The rounded repo
point (-0.059, -0.108, -0.03) is disqualified as a transcription artifact.

> **FOR BEN — read before accepting:** your rule ("keep whichever path gives
> human-value-or-greater torque magnitude for hip ADDUCTION/ABDUCTION **and**
> hip FLEXION/EXTENSION") returned **"neither path passes both"** — the two
> candidates trade torque between the axes (table below). The tie-break applied
> on your behalf: (1) the rounded repo point is not a design but a corruption
> artifact (it numerically duplicates the femur-side P2_0, creating a
> zero-length via segment at the default pose — this is what produces the
> mid-range flexion null and the +/− sign flip the human reference doesn't
> have), so it is not a legitimate candidate; (2) the adduction miss of the
> winner (0.82x human) is a **muscle-sizing gap, routed to the
> Opt_run/Mesh_Optimization lever** (MIF/route sizing — that program's literal
> purpose), not a path choice; no path point fixes both axes anyway. **You can
> overrule any of this** — the one-line diff in `gait2392_robotbody.osim` is the
> only repo model change (git-clean before, single `<location>` line changed).

## 1. What was compared

- **variant_repo** — P2 = (-0.059, -0.108, -0.03): as shipped in
  `gait2392_robotbody.osim` (and _hip, robot, gait2327). Numerically identical
  to `add_mag3_r-P2_0` (femur_r) → zero-length P2→P2_0 segment at q=0.
- **variant_thumb** — P2 = (-0.1357, -0.0929, +0.0591): the master's
  thumb-drive original (per AGENTS.md 2026-09-22 packaging block).
- **human_ref** — stock gait2392 route (`ConnorBipedal.osim`'s add_mag3_r:
  P1 pelvis (-0.11108, -0.11413, 0.04882) → P2 femur_r (0.007, -0.3837,
  -0.0266)) with the project's own human data: **Adductor Magnus 3, MIF 488 N,
  OFL 0.131 m, TSL 0.249 m, pennation 0.08726646 rad**
  (`Code\Matlab\Mesh_Optimization\Add_Mag_Mesh_Opt.m:93-99`). This is the human
  reference the Mesh_Optimization program uses for this muscle. Thelen params
  are byte-identical across ConnorBipedal / robotbody / stock
  `gait2392_simbody.osim`, so **P2 was the only variable** (verified: 338 path
  points scanned, exactly one differs).
- Thumb drive NOT mounted on easteregg2 (`D:\Bipedal humanoid` absent); AGENTS.md
  + full git history (`git log --all -S"-0.1357"` → only a coincidental
  `-0.13572558` in `ResultsBSolve\zz_bsolve_MuscleAnalysis_TendonPower.sto`,
  commit 6e3776fb) confirm the two coordinate sets are the complete difference.

## 2. Method

OpenSim 4.6 python API (`D:\Anaconda\envs\opensim\python.exe`), scripts +
logs + CSVs in `add_mag3_compare_20260926\`:

- Grid = the project's hip RoM for this muscle (Add_Mag_Mesh_Opt.m:41-52):
  **flexion -25..+85 deg, adduction -45..+20 deg**, 23x14 = 322 poses
  (model clamped ranges are ±120 deg — anatomically impossible poses kept out
  of the ruling; supplementary `arm_tau_*_full.csv`).
- Moment arm = `computeMomentArm` about `hip_flexion_r` / `hip_adduction_r`;
  **validated** vs central-difference dl/dq (eps 1e-5 rad): exact magnitude
  match, 8/8 probes (`add_mag3_validate_and_profile.py`, `validate_log.txt`).
- Torque conventions, both reported:
  - **tauA = |arm| x 488 N** (project convention, MonoMuscleData) — primary.
  - **tauB = arm x F_tendon**, Thelen2003 isometric equilibrium at
    activation = 1 (add_mag3_r only; full-model equilibrium fails on Ben's
    rerouted semimem_r, irrelevant here).

## 3. Results (primary grid, peaks; full curves in `arm_tau_*_primary.csv`)

| hip torque envelope | repo (rounded) | thumb (full-precision) | human (Add Mag 3) |
|---|---|---|---|
| **flexion** peak (tauA) | 19.14 N·m = **0.63x human**, arm crosses **ZERO at ~45 deg flexion** (0.30 N·m) and **flips sign** (+3.5 cm at -25 deg → -2.8 cm at +85 deg) | 34.28 N·m = **1.12x human**, above human at 252/322 poses (min ratio 0.91 only at the extreme 85-deg corner), sign-consistent like human | 30.61 N·m (peak arm 6.27 cm) |
| **adduction** peak (tauA) | 49.96 N·m = **1.49x human**, above human at 319/322 poses | 27.62 N·m = **0.82x human**, below human at 310/322 poses | 33.60 N·m (peak arm 6.88 cm) |
| mean tauA flex / add | 8.80 / 41.23 N·m | 26.62 / 11.22 N·m | 22.80 / 17.35 N·m |

tauB (equilibrium, act=1): same split directionally — repo passes adduction at
every pose (min ratio 10.4) but fails flexion at 53/322; thumb passes flexion
at every pose (min ratio 1.92) but fails adduction at 46/322.

## 4. Ruling

1. **Rounded repo point DISQUALIFIED** as a corruption artifact (P2 == P2_0
   numerically → degenerate via → flexion null + sign flip; its 1.49x adduction
   is a property of broken geometry nobody chose).
2. **Full-precision thumb point KEPT**: passes flexion outright (1.12x),
   sign-consistent on both axes like the human.
3. **Adduction shortfall (0.82x) recorded as a MUSCLE SIZING gap** → handle in
   Opt_run/Mesh_Optimization (MIF/route sizing), not by path-point choice.

## 5. Where the winner landed (and what was NOT touched)

- **APPLIED**: `gait2392_robotbody.osim` — one line, `add_mag3_r-P2`
  `<location>` → `-0.1357 -0.0929 0.0591` (git diff = 1 file, 1 insertion,
  1 deletion; file was git-clean before the edit — commit via GitHub Desktop).
- **REGENERATED + VERIFIED**: `repair_robotbody_muscles.py` rerun (myo env
  python) → `gait2392_robot.osim`: "mirrored 46 right->left GeometryPaths",
  "non-mirrored pairs: 0", "VALIDATION OK", left paths equal to official
  simbody left = 21/46 (21 expected, unchanged). `add_mag3_l-P2` =
  (-0.1357, -0.0929, **-0.0591**) (exact z-mirror). OpenSim probe: regenerated
  robot.osim add_mag3_r **== thumb variant bit-exactly (max dev 0.00e+00)**
  at 4 poses spanning the grid (`verify_regen.py`, `verify_log.txt`).
- **NOT touched** (correct, or needs your call):
  - `ConnorBipedal*.osim` + thumb drive master's data — never modified.
  - `gait2392_robotbody_hip.osim` — carries the old rounded P2; side copy, not
    in the repair chain. Sync it with the same one-line change if you use it.
  - `gait2327.osim` (generated) — still carries the old P2; regenerate via the
    laptop-side XML-surgery build (`D:\GitHub\myoconverter\build_gait2327.py`).
  - **Downstream MuJoCo/spinal**: `mjc\gait2392_robot\` conversion + any
    spinal references/pre-tuned curricula were built from the OLD path — they
    are now one path-point stale for add_mag3(_r/_l). Regeneration (MyoConverter
    route) + any gait re-checks are yours to schedule (spinal npz/db files left
    untouched).

## 6. Separate model-quality finding (independent of this ruling)

**Thelen-equilibrium forces are inflated ~10x for BOTH robot variants by
passive overstretch of the robot P1/route** (robot P1 (-0.163, -0.013, 0.005)
makes the route longer than the human's; at the default pose lce = 0.253 m =
1.93x OFL → 4509 N tendon force vs 488 N MIF; full-range peaks reach ~14 kN at
±120 deg clamped poses). tauB peaks therefore do not discriminate between
variants and are not usable as "capability" numbers without addressing the
route length / resting-length mismatch. Worth its own look (P1 choice, TSL, or
route rework) — it affects every pose of the current add_mag3(_r/_l) routes in
the converted models.

---
*Artifacts: `add_mag3_compare_20260926\` — variant .osim copies, analysis
scripts, `run_log.txt`, `validate_log.txt`, `verify_log.txt`,
`arm_tau_*_{primary,full}.csv`, `add_mag3_torque_summary.json`.*

---

## Ben's framing correction (2026-09-26 late) — read before using the table above

Ben clarified the comparison doctrine: the meet-or-beat test is **human muscle on its OpenSim route
vs a BPA on a MODIFIED route over the same RoM** — NOT route-variant vs route-variant. The peak-ratio
table above therefore settles **PROVENANCE ONLY** (rounded repo P2 = transcription corruption — it
duplicates P2_0, giving the zero-length via segment and the flexion null/sign flip; full-precision
thumb P2 = the original; Ben confirmed such transcription errors are likely in his old brute-force
pre-objects/classes code). The performance question belongs to Opt_run/Mesh_Optimization: a BPA on
an OPTIMIZED route vs the human target, judged on GROUP-TOTAL torque profiles over fuller RoM where
needed (one 40 mm BPA vs the adductor group; sometimes 2 BPAs in parallel for one human muscle),
under physical-packaging constraints. The 0.82x adduction number is NOT a failure verdict — it is
the starting point for that sizing optimization.

Context facts Ben supplied: Festo unrestrained contraction at 620 kPa ~ 14-17% (10 mm), ~25%
(20 mm), ~28% (40 mm), vs >= 33% for biological muscle; Festo BPAs produce maximum force at
RESTING length, whereas human muscle peaks active force at OPTIMUM FIBER LENGTH and produces force
on either side of it.
