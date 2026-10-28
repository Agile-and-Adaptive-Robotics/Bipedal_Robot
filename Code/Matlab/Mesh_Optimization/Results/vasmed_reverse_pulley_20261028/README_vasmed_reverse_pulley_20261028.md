# Reverse-pulley vastus medialis placement optimization — result folder (2026-10-28, easteregg2)

Folder: `Code\Matlab\Mesh_Optimization\Results\vasmed_reverse_pulley_20261028\`
Purpose: verified result fragments + figures for the torque-paper Section 4.4
opening example (reverse-pulley mesh optimization of vastus medialis BPA
placement, before extending to the master's-routing survey). No manuscript
sources were edited; Ben copies what he wants back to the editor.

## Result verdict (updated at the end of the run — see VERDICT section)

See the bottom of this file. The numbers below in [BRACKETS] are filled by
the run logs (`vasmed_rp4_final_20261028.txt`, `vasmed_rp6_verify_20261028.txt`).

## What was actually wrong with the 10-27 study (audit result)

1. `study1_vasmed_canon_20261027.m` optimized ONLY `[rest, tendon]` on the
   FIXED route of record — a length-sizing sweep. Attachment placement
   (origin p1, insertion pEnd) was never optimized. Its "4/5 nodes, best
   n=3/G=2" claim is therefore not evidence of placement optimization.
2. Its grading used `min(margin)` over nodes, and MATLAB `min()` silently
   omits NaN, so infeasible nodes could not lower the reported minimum
   margin — the "misleading minimum-margin" risk.
3. **Structural finding (the important one):** `MonoPam_pulley`'s closure
   assumes `mL = const + crossed-segment span` (its (1,1) regression
   identity is asserted only for "a geometry whose crossed segment is the
   tendon line", `Opt_sanity_pulley.m`). On the distal-ring extensor route
   the femur-side contacts RELEASE toward extension (probe
   `vasmed_rp0_probe_20261028.m`: at −120° the crossed segment is 0.107 m
   with a 0.607 m femur chain; at −19° the chain is collapsed and the
   crossed segment is 0.634 m), so the identity fails by up to ~0.6 m of
   route bookkeeping through the flexion region. Both the 10-27 "success"
   and a naive re-grade are artifacts of that violation.
   `BiPam_pulley2_20261022` (the canonical BiPulley class of the capability
   survey) exists precisely for this case: its tendon span is the whole
   exit-to-insertion routed path and its closure references each cell to
   the per-cell `s0(i) − DeltaL(i)` (`fortzBi2`).
4. Xi3/wrapping loss: `MonoPam_pulley` carries no Xi3 term at all. Per
   Ben's 2026-10-28 ruling the BPA itself bends on this route (femur
   cylinder p2 + condyle/clip p3–p5 are BPA-side contacts; only the tibia
   ring p6–p8 and crosspoint are cable-side), so Xi3 MUST be incorporated.
   It is carried INSIDE the equilibrium through the evaluator-of-record
   formula `delta_L = Xi3·bendMeasure·comp²` (`MonoPamDataExplicit_balanceX3
   Contraction_k`), applied per pose.

## Evaluator of this study (no class edits, no whole-class 1-D port)

The `BiPam_pulley2_20261022` closure equation (`fortzBi2`, cAdd = 0)
applied per pose:

```
psi(s) = G·s − G·cb_i·F(s) − nB·F(s)·cSpr/G − DL(i) − G·(s0x3(i) − DL(i))
```

with `s0x3(i) = rest − mL(i) + Xi0 + tendon + 2·fitting + delta_L(i)`.
Every per-pose input is machinery of record: route/mL/bend from
`buildDistalRingLocation20mm`; DL = the class's exit(row 5)→insertion
(row 9) multi-segment span change; cb = the verbatim `LokBi2` bracket
projection (dia-20 Pbr = [9.48, −36.21, 30.86]/1000, pulley-family table);
F = festo4/maxBPAforce; cSpr = wraps(3)·Aeff·E/TendonL (pulley-family Spr).
Two gates verify the implementation:

- **GATE 1** (nB=1, G=1) vs `MonoPamDataExplicit_balanceX3` (BPAcount 1,
  per-pose Location, Xi3 native) on the route of record: [NUMBER].
  Residual sources disclosed: the X3 family's kSpr = 18·Aeff·E/L vs the
  pulley family's 3·Aeff·E/L and Pbr y = −33.38 mm vs −36.21 mm (both
  quantified below); both are negligible next to the bracket compliance
  (cSpr ≤ 3.5e-7 vs cb ≈ 2.5e-5 m/N).
- **GATE 2** (nB=2, G=2) vs `BiPam_pulley2_20261022` itself on a FIXED-ROW
  control route (uncollapsed seed20 rows, 100-pose grid): [NUMBER].

## Conventions and disclosures

- Pressure 620 kPa, 20 mm BPAs, KMAX 0.255, fitting 0.021 m, wraps 3,
  adopted extensor chi pick 32 (Xi0 −0.00637802 m, Xi1 39978.7 N/m,
  Xi2 14733.8 N/m, Xi3 0.158289).
- Sign/domain: knee angle −120…+10 deg, flexion NEGATIVE; human vastus
  medialis (vas_med_r, single muscle, pinned OpenSim Gait2392 simbody,
  Fmax 1294 N) extension torque POSITIVE (30.2–59.9 N·m over the domain).
- Arm methods kept separate: (a) OpenSim ds/dθ arm; (b) routed cable arm
  about the measured moving ICR (taut-cable construction of study1,
  verbatim). The optimization grades on (b); (a) is reported alongside.
- Grading: taut = ~PulleySlack & ~PulleyInfeasible & Ftendon > 0 (actual
  positive delivered force; never reconstructed strain_p or nonzero
  torque). Slack/infeasible poses are FULL shortfalls. Worst margin is
  over the ENTIRE domain, never over successful nodes only. Both plain
  and 1.05× (route-design 5% margin) targets tracked.
- Tendon floor 0.075 m = pulley packaging constraint (2026-10-28
  directive); it is the ROUTE-MODEL tendon parameter only — the physical
  line looped through a G-part tackle and terminated at the insertion is
  longer; the route model does not carry that internal tackle line.
- Attachments: p1 = femur-frame BPA ORIGIN (route row 1), pEnd = t1-frame
  tibial INSERTION (row p9) — the tibia-route solver may optimize BOTH
  (Ben 2026-10-28). The patella-insertion solver family (insertion given
  on the femur as a function of knee angle) must keep its insertion fixed
  and optimize only its origin — not used here.
- Deviation penalty vs the anatomical seed rows (routeSeed 1/9): weight
  0.5·dev/0.10 m in the attachment objective (magnitude Ben left open;
  the weight is small enough that a 0.01 worst-domain meet gain dominates
  any plausible deviation change). Rest-length term 0.5·rest/0.83 per
  Ben's ruling that 0.8277 m is too long.
- kSpr convention audit (the "unverified strand multiplier" item):
  identification evaluator `minimizeExtX3` uses 2·Aeff·E/L (no wraps
  term); `MonoPam_pulley`/`BiPam_pulley2` use wraps·Aeff·E/L (3×);
  `MonoPamDataExplicit_balanceX3` uses 6·wraps·Aeff·E/L (18×). All three
  are effectively rigid next to the bracket compliance (force effect
  < 0.1%); the study uses the pulley family's 3×. FLAGGED for Ben: the
  strand-count interpretation needs a hardware call.
- Mount/bowden deflection (BiPam_pulley2 variant-2 fields) is OFF,
  matching the recorded default configuration of the survey batch. The
  mount reaction R = nBPA·F·u_bpa + F_t·u_t (the tackle-appropriate
  internal force, NOT the reduced tendon force) is the class's reported
  quantity; bracket deflection per BPA uses the per-BPA force (each BPA
  its own bracket) per the pulley-family convention.
- The poster file named in the task ("Bolen et al. - 2024 - Self-balancing
  Robot Leg.pdf") does not exist in this repo (checked 2026-10-28; nearest
  on disk: `4008_Bolen_Ben_PDF.pdf`). The schematic follows the
  MonoPam_pulley/BiPam_pulley2 hardware contract and the actual vas-med
  route of record instead.

## Files

- `vasmed_rp0_probe_20261028.m` — closure-identity probe (the audit
  evidence); `vasmed_rp0_gate_diag_20261028.m` — first-gate diagnostics
  (metric-artifact analysis; superseded by rp4's gates).
- `vasmed_rp1_gate_regrade_20261028.m` + outputs — the SUPERSEDED
  MonoPam_pulley regrade (kept as the record of why it is wrong).
- `vasmed_rp4_final_20261028.m` — the study of record (gates, 15-config
  sizes sweep on the fixed route, attachment optimization).
  Outputs: `vasmed_rp4_final_20261028.mat`, `vasmed_rp4_sizesweep.csv`,
  `vasmed_rp4_attachopt.csv`, `vasmed_rp4_final_20261028.txt`.
- `vasmed_rp6_verify_dense_20261028.m` — fresh-process 131-pose dense
  verification. Outputs: `vasmed_rp6_dense_20261028.mat`,
  `vasmed_rp6_dense_curves.csv`, `vasmed_rp6_verify_20261028.txt`.
- `vasmed_rp5_fig_schematic_20261028.m` / `vasmed_rp7_fig_torque_20261028.m`
  — figure sources; exports
  `vasmed_reverse_pulley_schematic_20261028.{fig,pdf,eps,png}` and
  `vasmed_reverse_pulley_torque_20261028.{fig,pdf,eps,png}`.
  Unique `vasmed_reverse_pulley_*` names (not coordinated with the
  Knee_Torque_Test/Dissertation figure folders — nothing was copied there).
- `Components\` — verbatim input copies (human surface CSVs + manifest,
  the 10-27 study1 summary, the cam design mat of record).
- `SHA1SUMS_20261028.txt` — hashes of inputs, machinery of record, and
  all scripts.

## Solver settings

- Sizes sweep: `patternsearch`, MeshTolerance 1e-4, StepTolerance 1e-4,
  MaxFunctionEvaluations 150, five seeds; rest bounds [0.30, 0.90]·Lchar,
  tendon [0.075, 0.30] m.
- Attachment optimization: `surrogateopt` (150 evals, UseParallel,
  MinSampleDistance 1e-3, ConstraintTolerance 1e-6, two seeds) +
  `patternsearch` polish (150 evals, tolerances 1e-4), constraints
  `nonlconExt20mm` (canonical: tendon vs geometry-dependent max,
  minimum-possible tendon, pEnd envelope, p8 x-limit, zero-strain rest
  feasibility, Contraction/KMAX ≤ 1 — evaluated through
  `MonoPamDataExplicit_balanceX3` with the adopted Xi3).
- Objective (arm b): −1000·worstFrac − 100·meanFrac + 0.5·rest/0.83 +
  0.5·dev/0.10.


## VERDICT (session close, 2026-10-28 evening, easteregg2)

1. AUDIT (fully verified, evidence in this folder): the 10-27 study
   (`study1_vasmed_canon_20261027.m`, output `VasMed_minBPA_study1_*`)
   (a) optimized only [rest, tendon] on the fixed route of record -- a
   length-sizing sweep, NOT attachment placement; (b) graded through
   MonoPam_pulley, whose closure identity (mL = const + crossed span)
   is VIOLATED on this route through the flexion region by up to ~0.6 m
   of route bookkeeping (probe: vasmed_rp0_probe_20261028.m -- at
   -120 deg the crossed segment is 0.107 m with a 0.607 m femur chain;
   at -19 deg the chain is collapsed and the crossed segment is
   0.634 m); (c) reported min-margins via MATLAB min(), which silently
   omits NaN nodes; (d) carried no Xi3 although the BPA itself bends on
   this route (Ben's 2026-10-28 ruling). Its "4/5 nodes, n=3/G=2"
   numbers are WITHDRAWN as evidence for anything about placement
   optimization; do not copy them into the manuscript.

2. CORRECTED EVALUATION (this study, rp4): the BiPam_pulley2_20261022
   closure (whole exit-to-insertion span, per-cell reference
   s0(i) - DL(i), fortzBi2 lines 1515-1517, cAdd = 0) applied per pose,
   with Xi3 inside the equilibrium (delta_L = Xi3*bend*comp^2, the X3
   evaluator's own formula; Xi3 = 0.158289, pick 32). Gate status:
   GATE 1 vs MonoPamDataExplicit_balanceX3 at (n=1, G=1): max rel
   force err 8.56% (mean 3.42%) -- NOT the 1e-3 target; kSpr
   convention ruled out as the source (aligned variant: 8.52%);
   remaining suspects: the X3 family's Pbr (y -33.38 vs -36.21 mm) and
   its r-bracket edge behavior; UNISOLATED at session close. GATE 2
   (fixed-row G=2 control vs BiPam_pulley2 itself): NOT RUN -- the
   class's TransformationMat 4th dim is the joint index; constructing
   the campaign-driver layout ran past the deadline. So the G>1 psi
   path is verified by verbatim transcription only, not bit-level test.

3. VERIFIED SWEEP RESULT (fixed route of record, tendon >= 0.075 m,
   sizes re-solved per config, full 100-pose domain, slack/infeasible
   poses = full shortfalls; log rows in vasmed_rp4_final_20261028.txt):
     - nB=2, G=1 (direct drive, two parallel 20 mm BPAs, rest 0.7027 m,
       tendon 0.0925 m): meets the human single-muscle vastus medialis
       target at ALL 5 nodes AND the +5% design-margin target, worst
       full-domain margin +6.77 N*m (+11.4%), 0 slack, 0 infeasible,
       100/100 taut. Arm-method (a) worst +5.05 N*m.
     - nB=3, G=1: worst +39.59 N*m (+67.0%), 5/5 (+5% 5/5).
     - nB=4, G=1: worst +71.99 N*m (+122.5%), 5/5 (+5% 5/5).
     - EVERY tackle configuration (G = 2 or 3, all nB): FAILS -- 25-50
       slack poses, worst -59.91 N*m (-100%) at the +10 deg node.
     - nB=1, G=1: worst -15.10 N*m (-25.2%), 3/5 nodes.
   Classification of this result: a REGRADED SAVED CONFIGURATION (fixed
   attachments of the cam design of record, sizes re-solved). NOT a
   newly optimized attachment route, NOT physically validated.

   DENSE VERIFICATION (rp6, fresh process, 131 poses): the nB=2/G=1
   winner HOLDS -- worst +6.73 N*m (+11.4%), 5/5 nodes, 5/5 at +5%,
   0 slack, 0 infeasible; reproduces the 100-pose search value (+6.77)
   to 0.04 N*m. All other designs reproduce to <= 0.05 N*m. The saved
   design of record as-is (n=1, direct drive) verifies at -47.4 N*m
   worst (2/5) -- the direct-drive baseline for comparison.

   STAGE C (attachment optimization, REDUCED BUDGET: surrogateopt 150 +
   patternsearch 150 evals on 8 dimensions): did NOT improve on the
   fixed-route designs -- best attachment-optimized result 4/5 nodes,
   worst -36.1 N*m (nB=4), versus the fixed route's 5/5 at +6.8 N*m
   (nB=2). VERDICT: budget-starved search failure, not a placement
   finding -- the objective values prove the fixed-route design scores
   far better than what the 8-D search found; the campaign-of-record
   budgets (1000+ evals) were not affordable in this session. Do not
   read Stage C as placement optimization does not help; read it as
   not demonstrated tonight.

   Physical reading (for the Discussion): on this extensor route the
   tackle never pays -- the per-cell reference shows the BPA-side route
   shortens toward extension, so the tackle's span-demand credit is
   outweighed by the 1/G force division; the tendon goes slack at the
   flexion end. The enabling mechanisms are PARALLEL BPAs and (pending
   Stage C) attachment placement, not the reverse-pulley transmission.

4. UNVERIFIED / INCOMPLETE at session close:
   - GATE 1 residual isolation (8.6% max vs the X3 class; kSpr ruled
     out) and the GATE 2 fixed-row control (never constructed).
   - Attachment placement optimization at campaign-of-record budget
     (Stage C ran at 150+150 evals and underperformed the fixed route).
   - Figures RENDERED (schematic + torque-angle, FIG/PDF/EPS/PNG,
     unique vasmed_reverse_pulley_* names; nothing copied into the
     manuscript figure folders). The torque figure plots all six
     rp6-verified designs; for the paper prefer the saved-design
     baseline + the nB=2 winner + the best attach-opt row.
   - No physical validation anywhere in this folder.
   - Nothing was copied to Knee_Torque_Test\Figures or
     Dissertation\Figures-torque (unique vasmed_reverse_pulley_*
     names only, no coordinated names).

5. Ready-to-copy manuscript fragments: LaTeX_FRAGMENTS_20261028.tex in
   this folder (Methods / Results / Discussion), written so every number
   matches the logs above; bracketed items mark what changes if Stage C
   or rp6 land.
