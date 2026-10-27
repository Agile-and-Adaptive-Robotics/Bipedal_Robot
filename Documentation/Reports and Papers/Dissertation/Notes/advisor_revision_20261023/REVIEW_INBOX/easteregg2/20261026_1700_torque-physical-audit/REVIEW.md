# Review batch: torque-paper physical-data audit and figure repairs

- Machine: easteregg2
- Repository state: working tree at the repo-root merge `ec7488f9` plus this
  session's uncommitted edits to the paper (`Knee_Torque_Test/ProofFinal/Main.tex`,
  regenerated figures). Prose source snapshots of record: the **live Overleaf
  states captured 2026-10-26** in `Documentation/Reports and Papers/Knee_Torque_Test/
  review_audit_20261026/` (`dissertation_30-torque_after_GoF.tex`,
  `dissertation_94-AppendixC_live.tex`, `paper_Main_after_GoF.tex`). Every
  boundary below is quoted from those snapshots; local line numbers name that
  snapshot. Re-check each boundary against current Overleaf text before applying.
- Purpose: P1 physical-data corrections (extensor 20 mm series identity, Hunt
  equation, FlxGrp 2x2), training/validation memberships per Ben's ruling, and
  the P2 repairs, delivered as bounded patches. No new optimization was run
  anywhere in this batch; every number comes from existing data files, the
  adopted result mats, Ben's read-only evaluators called at the adopted terms,
  or the compare sections of his analysis scripts rerun verbatim.
- Companion paper changes: the same repairs were applied directly to the paper
  (Main.tex, double-built clean, 38 pages). Regenerated paper assets:
  `Figures/Ext20mm_cam_20261026.eps`, `Figures/FlxPin_allmodels_20261026.eps`,
  `Figures/ExtPin_group_20261026.eps` (SIX-panel pool layout), `Figures/
  FlxGrp_20261023.eps` (2x2 replacement, same name), and `Figures/
  ExtensorRoute_redesigned.eps` (4-pose replacement, same name), each with a
  rebuilt `-eps-converted-to.pdf` cache.

## Rulings applied (Ben, 2026-10-26 evening) — supersedes this batch's first pass

1. **Torque pipeline of record.** The measured torque results are whatever the
   compare sections of `Code/Matlab/Mesh_Optimization/Results/Knee_Extensor_20mm.m`
   (and `Knee_Flexor_Data_20mm.m` on the flexor side) compute from the
   measurements. The xlsx rows 16 and 17 of the ExtTest20mm_2 tab are SANITY
   CHECKS ONLY — never reconcile against them, never treat them as records.
   The extensor compare section was rerun verbatim (script untouched; dated
   runner `run_ExtCompare_record_20261026.m`; verbatim output in
   `run_ExtCompare_record.log`). Its nine torques of record (N m, tab order):
   50.5028, 51.1732, 47.1236, 50.4992, 32.5183, 29.0224, 35.9028, 48.7059,
   38.9060 at -26.2/-30/-61/-43/-119/+5/-1/-22/-86.5 deg. **Maximum difference
   against the values this batch had already published: 0.0004 N m**, so no
   figure, prose number, or CSV value changed; the paper's extensor numbers
   and both audit CSVs stand as the record's rounded form. Full record:
   csv/ext20mm_record_rerun_20261026.csv. Its predAP and 620 kPa series also
   match the batch's cam-design values exactly, confirming the compare
   section's adj route and the cam mat of record evaluate the same design.
2. **Extensor fold.** Ben rules the SIX-test pool stands: training {43.6,
   43.2-tendon, 45.7}, validation {40.0, 40.0-tendon, 46.5}; unused the two
   41.5 cm tests plus 48.0 (unmentioned in the paper per the standing rule).
   The paper was restored accordingly (see the delta list below) and the
   figure regenerated as the six-panel pool layout. The mat-field reading
   that motivated the earlier five-pool edit is RETRACTED as the membership
   of record and kept only as a flagged discrepancy: **`minimizeExt10mmX3_results.mat`
   carries `allBPA = [1 2 5 6 8]` and `numBPA = 5` with `scores_cv` 5x3 and
   `validLabels` listing five tests, while row 32's fold columns read (1, 2, 8)
   in the nine-test space — Ben's ruling is that the six-test pool with 45.7
   in training is what was used; the mat's saved pool fields disagree with
   that account and should be read with care in any future audit.** AGENTS.md
   needs no update; its original sentence was right. Appendix C's existing
   six-pool text is correct and needs no patch (this batch's P13 was withdrawn).

Numerical deltas from the two rulings: ruling 1 changed nothing (verified
<= 0.0004 N m everywhere). Ruling 2 restored, in the paper: tab:ExtPinGof
training block regains the 45.7 cm row (2.178/1.0586/4.974 char, 0.490/
0.0535/1.206 corrected) and the training means revert 2.658/2.8461/4.551 ->
2.498/2.2503/4.692 (char) and 1.251/0.7452/2.150 -> 0.997/0.5146/1.835
(corrected), re-verified from the evaluator at row-32 terms; the Discussion
sentence reverts "training mean over the 43.6 and 43.2-tendon tests, 2.66 ->
1.25" to "over the 43.6, 43.2-tendon, and 45.7 tests, 2.50 -> 1.00"; the
figure caption regains panels (A)-(C) training / (D)-(F) validation with the
legend on panel (C). Validation means were never changed (2.554/2.4236/4.986
char, 1.098/0.3965/1.882 corrected).

## Fold memberships of record

| Configuration | Membership of record (Ben's ruling; matches the paper's captions) | Notes |
|---|---|---|
| Flexor (row 77, `minimizeFlxPin10_results.mat`) | training {48.5, 41.0}; validation {45.7, 40.6-tendon}; 47.9 unused | `filtered_results(77, 2:3)` = (2, 4); family 1=48.5, 2=45.7, 3=47.9, 4=40.6t, 5=41.0; tag noT3 |
| Extensor (row 32, `minimizeExt10mmX3_results.mat`) | **six-test pool: training {43.6, 43.2-tendon, 45.7}; validation {40.0, 40.0-tendon, 46.5}; unused {41.5, 41.5-tendon, 48.0}** | row-32 fold columns (holdCols 2-4) = (1, 2, 8) agree with the validation set; the mat's saved `allBPA`/`numBPA`/`validLabels` fields say five tests — flagged discrepancy, overruled by Ben (see Ruling 2) |

Campaign resting lengths (mat `bpa`): 1=0.400, 2=0.400+40t, 3=0.415, 4=0.415+40t,
5=0.436, 6=0.432+22t, 7=0.457, 8=0.465, 9=0.480 m; nominal campaign labels
round oddly (nominal "42 cm" is the 41.5 cm test, "43 cm-tendon" is 43.2 cm,
"46 cm" is 45.7 cm, "47 cm" is 46.5 cm). Paper-only rule applied: unused
tests (flexor 47 cm; extensor 41.5, 41.5-tendon, 48.0) are not mentioned in
Main.tex prose or captions; the dissertation keeps full fold bookkeeping.

## Patch manifest

| Fragment | Target chapter / location (live snapshot line) | Exact existing start/end text | Advisor/Ben comment | Change and reason |
|---|---|---|---|---|
| P01_initial_flexor_paragraph.tex | 30-torque, Results, flexor paragraph (line 375) | STARTS "The initial identification, fitted on the $l_0=\qty{48.5}{\cm}$ test, yielded a length offset" ENDS "…\qty{5.1690}{\newton\meter}, respectively (Table~\ref{tab:TorqueFit})." | Ben 2026-10-26 14:32: "This is incorrect in the values and which tests were training and validation. The whole paragraph is incorrect." | Newly checked result: paragraph now states the two-bracket identification of record (row 77; training {48.5, 41.0}, validation {45.7, 40.6-tendon}) with its adopted terms and Table tab:FlxPinGof scores (2.165/1.084 train, 1.368/1.611 validation); the single-bracket 48.5-fit/45.7-validation story is demoted to the historical pointer (Appendix). Numerical changes: chi0 11.9 mm -> +3.94 mm, chi1 102 kN/m -> 39.98 kN/m, chi2 10.6 kN/m -> 14.73 kN/m; RMSE/FVU/Max 1.616/0.0958/3.2527 (train) and 2.1214/0.0878/5.1690 (validation) removed (they belonged to the superseded fit). |
| P02_initial_extensor_paragraph.tex | 30-torque, extensor results (line 425) | STARTS "Results are shown in Fig. \ref{fig:ExtPin_group}. The initial identification was fitted" ENDS "…shows the six pool tests at the adopted terms." | Ben 2026-10-26 14:31: "This is all incorrect" | Keeps the historical cross-validated rounds sentence (Table tab:TorqueFit pointer intact) and the six-test pool fold exactly as the live text states it, adding the adopted terms chi0 = -6.38 mm and chi3 = 0.158 of the locked refit and the superseded-status framing for the initial terms. Pool/membership unchanged from the live text (Ben's ruling). |
| P06_FlxPinPreGof_Hunt.tex | 30-torque, Table tab:FlxPinPreGof (line 352) | STARTS "\begin{table}[H]" caption "Pre-optimization pinned-knee flexor torque-prediction goodness of fit" ENDS "\end{table}" | Ben 2026-10-26 14:30: "I don't think Hunt Equation was calculated correctly here." | Hunt column recomputed by direct call of the equation of record: 5.779/1.1761/9.798 -> 5.223/0.9477/9.534 (48.5), 6.633/3.2313/11.409 -> 6.249/2.6768/10.510 (41.0), 5.315/0.9408/9.991 -> 6.048/1.1466/9.956 (45.7), 3.945/0.6750/7.527 -> 3.956/0.6301/7.344 (40.6t). Caption gains the domain sentence (per-model = common domain here). Char and hybrid columns unchanged. |
| P07_FlxPin_preopt_figblock.tex | 30-torque, fig:FlxPin_preopt block (line 343) | STARTS "\begin{figure}[htbp]" containing "\includegraphics[width=0.9\textwidth]{FlxPin_preopt.eps}" ENDS "\end{figure}" | Same thread | Asset replacement with figures/FlxPin_preopt_20261026.eps (same content, corrected direct-Hunt series, M/theta-k labels); caption updated to say the Hunt line is evaluated directly from the equation of record. |
| P08_extensor_physical_final.tex | 30-torque, Physical Tests subsection (lines 595-604) | STARTS "The redesigned cam-extensor route was tested at nine isometric holds" ENDS the "\end{figure}" of the extensor block | Ben 2026-10-26 14:18: "This figure got replaced with an incorrect one. The blue dots are where the back-calculated yellow dots should be." and "Don't use 'holds'. Call it data measurements." | Newly checked result: the nine torques of record from the compare-section rerun (see Ruling 1; identical to the audit values to 0.0004 N m), "holds" -> "data measurements", gold series identified as a pressure-specific force-model prediction (not the hybrid, which needs measured muscle length/moment arm this test did not record), and the figure block retargeted to figures/Ext20mm_physical_20261026.pdf with the corrected series identities. |
| P09_FlxGrp_2x2.tex | 30-torque, fig:FlxGrp block (line 488) | STARTS "\begin{figure}[htbp]" containing "\includegraphics[width=0.72\textwidth]{FlxGrp.eps}" ENDS "\end{figure}" of that block | Ben 2026-10-26 14:44: "Make a 2x2 panel figure with the legend in the last panel. Title panels 10 mm, 20 mm, 20 mm." | Replaces the FlxGrp.eps + Flx20mm_meetsHuman.pdf minipage composite with the single 2x2 asset figures/FlxGrp_2x2_20261026.pdf; caption states the third 20 mm panel is the redesigned two-BPA route despite the shared title, carries the fifteen verified margins (+0.42 N m at -105 deg to +17.0 N m at -14 deg), and notes the shared B/C axis ranges. |
| P10_pulley_reference_fix.tex | 30-torque, Mono/Bi-pulley subsection (line 591) | STARTS "A subsequent transmission study re-solved actuator sizes" ENDS "…remains future work." | Coordinator relay 2026-10-26 (caption-figure mismatch investigation) | Removes the dangling "shown in Figures~\ref{fig:PulleySurfaces} and~\ref{fig:PulleyMargins_resolved}" sentence: the live Overleaf has already deleted both figure blocks, so both references are undefined. See E5 for the caption-mismatch resolution. |
| P10b_pulley_surfaces_figure_OPTIONAL.tex | 30-torque, optional re-insertion after line 591 paragraph | ALTERNATIVE block; apply only if Ben re-adds the surfaces | Coordinator relay 2026-10-26 | If the two pulley-surface images are re-added, this is their corrected caption: they are the five-site doctrine-batch surfaces (joint angle x pressure; margins along the transmission ladder), NOT a pulley-position sweep; no pulley-position-sweep figure exists yet (only raw survey mats from the pending sweep campaign). |
| P11_xiProjection_split.tex | 30-torque, fig:xiProjection block (line 265) | STARTS "\begin{figure}[htbp]" containing "\includegraphics[height=0.40\textheight]{xiProjection.pdf}" ENDS "\label{fig:xiProjection}\end{figure}" | Ben 2026-10-26 13:38: "Should be two separate figures. alpha can be marked on the circle." | Splits into fig:xiProjection (compliance projections) and a new fig:bpaWrapping (wrapping geometry, asset figures/bpa_wrapping_loss_20261026.pdf) with alpha drawn on the joint circle. |
| P11b_wrapping_prose_ref.tex | 30-torque (line 281, wrapping-loss paragraph) | STARTS "This excess is termed fictive strain" ENDS "such that" | Same thread | Repoints "Figure~\ref{fig:xiProjection}\textbf{B}" to "Figure~\ref{fig:bpaWrapping}" and mentions alpha on the circle. |
| P12_notation_M_kneebold.tex | 30-torque, chapter-wide (10 enumerated replacements) | Exact OLD/NEW fragment pairs inside | Ben 2026-10-26 13:20: "I don't like Torque as lower case m. Use upper case." and 13:24: "Make all subscript k's that refer to knee axis styled as \mathbf{k}" | Uppercase M for the torque/moment vectors and scalars (M_vec, M_scalar, M_z -> M_k), knee-axis subscripts k -> \mathbf{k} in r_k, z_k, theta_k, theta_text{k}; tibia-side theta_1 and basis indices untouched. Applied to the regenerated figure axis labels already (M_k, theta_k in bold k). |

Withdrawn after Ben's evening rulings (their targets are correct as they
stand): P03 (diss "six tests, three for training" sentence — correct), P04
(fig:ExtPin_group caption grouping — correct), P05 (diss tab:ExtPinGof
regrouping — the live six-pool table is correct), P13 (Appendix C fold
"corrections" — the live six-pool text is correct).

## Figure manifest

| Asset | Intended label | Replace/add | Content | Generator and data |
|---|---|---|---|---|
| Ext20mm_physical_20261026.pdf/.png (+alt.txt) | ch.4 fig:extensor_physical_final | Replace asset | Nine measured data measurements (indigo) vs model-at-measured-pressure (gold) vs optimized 620 kPa (magenta solid) vs human vasti target (black dotted); marker area scales with pressure; M_k / theta_k labels | `tmp/paper_bridge/figscripts/gen_Ext20mm_cam_20261026.m`; measured series verified against the compare-section rerun to 0.0004 N m (Ruling 1) |
| FlxGrp_2x2_20261026.pdf/.png (+alt.txt) | ch.4 fig:FlxGrp | Replace asset + block | 2x2: 10 mm / 20 mm / 20 mm(redesigned two-BPA route) + legend tile; A/B/C margin letters; shared B/C ranges | `tmp/paper_bridge/figscripts/gen_FlxGrp_2x2_20261026.m` (machinery verbatim from gen_FlxGrp_3panel_20261023.m) |
| FlxPin_preopt_20261026.eps/.png (+alt.txt) | ch.4 fig:FlxPin_preopt | Replace asset | Pre-opt 2x2 with the corrected direct-Hunt series; M_k / theta_k labels | `tmp/paper_bridge/figscripts/gen_FlxPin_preopt_20261026.m` |
| bpa_wrapping_loss_20261026.pdf/.png (+alt.txt) | ch.4 new fig:bpaWrapping | Add | Wrapping geometry with alpha marked on the joint circle | `tmp/paper_bridge/figscripts/gen_bpaWrapping_alpha_20261026.py` (adapted from tmp/revision_c_user_assets.py) |
| QC_Ext20mm_cam_20261026.png | QC preview only | n/a | 150 dpi QC render of the extensor asset | same generator |

Paper-side assets (already applied in place, listed for the record):
`Ext20mm_cam_20261026.eps` (new fig:Ext20mm), `FlxPin_allmodels_20261026.eps`
(new fig:FlxPin), `ExtPin_group_20261026.eps` (fig:ExtPin_group, restored
SIX-panel pool layout, legend on panel C, M/theta-k labels), `FlxGrp_20261023.eps`
(same-name 2x2 replacement), `ExtensorRoute_redesigned.eps` (same-name 4-pose
replacement with schematic femur/tibia outlines, legend tile, release-schedule
tile).

## Evidence

### E1. Extensor 20 mm physical data (P1-a) and the pipeline ruling

Pipeline of record (Ben, 2026-10-26 evening): the measured torque results
are whatever the compare section of `Knee_Extensor_20mm.m` computes from
the measurements; the xlsx ExtTest20mm_2 rows 16 and 17 are sanity checks
only and were not reconciled against for this batch's numbers. The compare
section was rerun verbatim by the dated runner
`run_ExtCompare_record_20261026.m` (log `run_ExtCompare_record.log`; the
script's leading `clear` wipes runner state after `run()` returns — known
`run()` behavior — but the script's own variables print before any use).
Record (N m, tab order): 50.5028, 51.1732, 47.1236, 50.4992, 32.5183,
29.0224, 35.9028, 48.7059, 38.9060; csv/ext20mm_record_rerun_20261026.csv.
These match the batch's published values to 0.0004 N m. Sanity check only:
rows 16/17 agree with the record to 4-7 N m / 0.1-1.25 N m respectively and
were not used.

The original figure defect and its mechanism, for the record: the old blue
series was the flexor-template adjoint with a post-hoc sign flip at reaction
point x = -0.065 (no source for -0.065 exists in `readserialnumbers2.m` or
its git history); that series tracks the predicted-at-pressure curve
(rms 5.1 N m), which is why the plotted blue sat where the back-calculated
yellow belongs. The compare-section reduction (reaction point (-0.050,
-0.255) m, d = 0.2599 m, force at the negated tab LC angle, adjoint to the
moving knee center, extension positive) is the pipeline of record and
reproduces the batch's values. Series-level audit:
csv/Ext20mm_series_audit_20261026.csv; convention-variant comparison:
csv/audit_Ext20_nine_20261026.csv.

Per-point record values (N m; measured / model at measured P / model at
620 kPa / human target): -119 deg: 32.5183/24.9824/104.6600/29.0023;
-86.5: 38.9060/40.9633/92.1960/38.6701; -61: 47.1236/44.5759/75.7037/
45.3468; -43: 50.4992/43.3722/73.2594/50.1669; -30: 51.1732/37.4228/
70.4051/49.8409; -26.2: 50.5028/38.7951/68.6262/48.9781; -22: 48.7059/
36.2951/66.3109/47.3913; -1: 35.9028/33.0097/53.5335/35.5885; +5: 29.0224/
25.7183/46.7667/27.8661. All nine data measurements meet or exceed the
human target (margins +0.24 to +3.52 N m) at 256-488 kPa supply; the
measured exceeds the pressure-specific prediction at 8 of 9 angles (max
+13.75 at -30 deg) and falls below once (-2.06 at -86.5 deg). Old paper
values (six points, sign-flipped series): 30.24/37.98/38.30/37.23/38.84/
22.86 with "only -119 meets the human target".

### E2. Hunt equation (P1-b)

`ForceStrainForFit.mat` z is the Hunt equation table in NEWTONS (30
relative-strain rows 0..1 over epsilon_max = 0.171, 19 pressure columns
linspace(0,620,19), isokinetic S=0, NaN below the first zero per column):
a fresh `HuntEq` call on the generator grid matches it to 4.17 N max
(0.8 percent, generation tolerance), so the table itself is sound. The
published series was produced by looking z up at relative strain
strain/KMAX (per-actuator 0.153-0.183 rather than the equation's own 0.171
basis), with the NaN region zero-filled and spline-interpolated, and with
a pressure axis mislabeled 32.63 kPa-spacing against the table's 34.44
(negligible at ~610 kPa). The published Hunt GoF column reproduces from
that machinery (RMSE and Max exactly; FVU within 0.06, a variance-base
difference in the historical scoring). The equation of record was then
called directly (epsilon = per-point relative contraction, epsilon_max =
0.171, S = 0, test pressure), giving the new column in P06 and the new
gray series in both figures; a per-actuator-KMAX variant was also scored
(5.924/6.783/5.631/3.715 RMSE, with a 16.56 N m worst residual on the 45.7
cm test) and rejected as conflating the fit's basis with the actuator's
kmax. Domains: the Hunt series is defined only where modeled strain is
non-negative; on all four tests that domain equals the common domain on
which all pre-optimization predictions are defined (n = 75/15/17/14), so
per-model and common-domain Hunt scores are identical; the char-alone and
hybrid common-domain scores are (4.871/6.074/4.465/3.226) and
(1.608/1.865/1.327/1.947) RMSE respectively. Series and mask:
csv/Hunt_series_audit_20261026.csv.

### E3. Adopted-term scores verification

Flexor row 77 (X0 = +3.9374e-3 m, X1 = 3.9979e4 N/m, X2 = 1.4734e4 N/m,
ledger-asserted) and extensor row 32 (X0 = -6.378e-3 m, X1/X2 locked,
X3 = 0.15829, ledger-asserted) evaluated with Ben's read-only evaluators
reproduce every corrected-column value in tab:FlxPinGof and tab:ExtPinGof
to three decimals under the six-pool grouping (training {43.6, 43.2t, 45.7}
means 2.498/2.2503/4.692 char and 0.997/0.5146/1.835 corrected; validation
means 2.554/2.4236/4.986 and 1.098/0.3965/1.882), and the extensor
zero-correction column equals the mat's `baselineScores` rows exactly.

### E4. FlxGrp panel C (P1-c)

Fifteen measured data measurements (56 cm flexor test): meet-or-exceed the
human biceps-femoris-short-head target 15 of 15, smallest margin +0.42 N m
at -105 deg, largest +17.0 N m at -14 deg (recomputed this session,
matching the caption numbers of record).

### E5. Pulley-surface caption mismatch (coordinator item)

No pulley-position-sweep figure asset exists in the repo (searched sweep/
pulpos/position names; only raw `BiPulley_survey2_tib_post_r_*` mats from
the pending sweep campaign). The stacked PulleyTorqueSurfaces_20261020 /
PulleyMarginSurfaces_20261020 images are the five-site doctrine-batch
surfaces (generator `gen_PulleySurfaces_20261020.m`: soleus/tib_ant/
bifemsh/med_gas/rect_fem over joint angle x pressure, margins along the
transmission ladder). The live Overleaf has already deleted that figure
block, leaving the dangling reference fixed by P10; if Ben re-adds the
images, P10b carries the corrected caption. No images were deleted by
this batch.

### E6. "Then we don't need to mention" (incomplete instruction)

Ben's follow-up comment (2026-10-26 14:35) is anchored at live-snapshot
character 35862, i.e. mid-paragraph in the "extended to the full test sets
in two stages" paragraph, immediately after "the encoder-corrected
l0=47.9 cm test was not used in the two-bracket identification". The
instruction is incomplete and the anchored clause is exactly the
unused-test bookkeeping the dissertation is supposed to keep, so nothing
was deleted; Ben to complete the instruction.

### E7. fig:RouteVariants supplementary trade-off (P2-6, for Ben's decision)

The Frontiers guideline allows up to 15 combined display items. Before this
batch the paper sat at exactly 15: 11 figure floats (testJigs, ICR,
xiFrameTransform, bktFrame, FlxPin, ExtPin_group, Ext10mm_52cm, FlxGrp,
RouteVariants, ExtensorRoute_redesigned, Ext20mm) and 4 table floats
(FlxPinGof as one sideways table, ExtPinGof, AdoptedXi, CapabilitySurvey —
CapabilitySurvey is a Main.tex float, not supplementary). The P2-2 page-fit
fix SPLIT tab:FlxPinGof into tab:FlxPinPreGof + tab:FlxPinGof, so the paper
now counts 16 combined display items, one over the guideline. Ben's
conditional comment on fig:RouteVariants ("candidate to move to
supplementary material if we are over the Figure/Table limit") therefore
now applies: moving fig:RouteVariants to the supplementary material brings
the count back to exactly 15. No move was made (out of scope of this
batch's instructions); the alternatives are (a) move fig:RouteVariants to
supplementary (16 -> 15), or (b) re-merge the two flexor GoF tables into
one portrait-orientation table (16 -> 15, undoing part of the page-fit
split). Ben's call.

## Validation performed

- Paper: `pdflatex` (MiKTeX, full path, in-place, no -output-directory)
  twice after every editing round including the two rulings; final state
  0 errors, 38 pages, 0 undefined references, 0 undefined citations (one
  intermediate pass hit a transient "can't write Main.pdf" file lock and
  the following passes were clean); all referenced figure files exist with
  fresh `-eps-converted-to.pdf` caches built from the Knee_Torque_Test root.
- Record rerun: the extensor compare section executed verbatim; its nine
  measured torques, predAP, and 620 kPa series captured and verified against
  the batch's published values (max delta 0.0004 N m).
- Figures: one MATLAB figure per asset, Arial, 12 pt axes (10 pt floor for
  point labels), campaign palette, measured = indigo filled circles,
  corrected = magenta line, unoptimized = light-orange dashed, Hunt = gray
  dashed, human = black dotted, gold = hybrid/back-calculated, marker-only
  series never combined with lines; (A)/(B)/(C) margin letters drawn in per
  the gen_ExtPin6Panel convention; QC PNGs rendered for every asset.
- Audit scripts (all under tmp/paper_bridge/figscripts/, dated 20261026;
  Ben's .m/.py sources untouched): run_ExtCompare_record_20261026.m,
  audit_Ext20_nine_20261026.m, gen_Ext20mm_cam_20261026.m,
  audit_Hunt_FlxPin_20261026.m (with dbg_Hunt_20261026.m),
  gen_FlxPin_allmodels_20261026.m, gen_FlxPin_preopt_20261026.m,
  dump_fold_memberships_20261026.m, dump_ext_folds2_20261026.m,
  dump_ext_pool_20261026.m, gen_ExtPin6Panel_20261026.m (plus the withdrawn
  gen_ExtPin5Panel_20261026.m), gen_FlxGrp_2x2_20261026.m,
  gen_ExtensorRoute_4pose_20261026.m (with dbg_geo_20261026.m),
  gen_bpaWrapping_alpha_20261026.py.
- Not verified this session (flagged): the historical single-bracket values
  (11.9 mm / 102 / 10.6; -10.9 mm / 0.123) are carried as the superseded
  fit's values from the previous record; the biomimetic extensor comparison
  at those initial terms (1.5394 N m, FVU 1.1560) was not re-evaluated.

## What Ben should commit (paper side)

`Knee_Torque_Test/ProofFinal/Main.tex` (base reconciliation from
paper_Main_after_GoF.tex + biomimicry typo fix + this batch's edits as
ruled), the five regenerated figure EPS files with their rebuilt PDF caches
in `Knee_Torque_Test/Figures/`, and the new scripts/CSVs/logs under
`tmp/paper_bridge/figscripts/`. The pre-repair Main.tex is preserved at
`tmp/paper_bridge/Main_pre_repair_20261026.tex`.

## Removed text of record (Ben, 2026-10-26)

Three passages Ben deleted from the paper as nonsensical or untrue (pre-deletion text preserved
in `tmp/paper_bridge/Main_pre_repair_20261026.tex` and the `review_audit_20261026/` snapshots).
Future sessions must not resurrect this class of claim:

1. **The mono-pulley freedom paragraph** (route-redesign results; pre-repair line 496, starts
   "The knee-flexor run quantifies what the mono-pulley freedom is worth at unity gain"). Wrong
   because: at one line part the tackle transmits force and travel unchanged, so there was no
   transmission effect to quantify, and the +4.51 N m peak difference the paragraph credited to
   "the redirect alone" came from two runs with different objective functions and independently
   re-optimized layouts, which the paragraph itself conceded were not commensurable; no causal
   "worth of the mono-pulley freedom" claim survives at unity gain.
2. **The bulk of the Limitations paragraph** (pre-repair line 552, starts "The redesign results
   carry caveats that bound their use"). Wrong because: it carried the stale extensor validation
   count "six points on the extensor cam route" from the defective sign-flipped series when the
   record is nine data measurements (see Ruling 1), and it speculated remedies the program has
   not established ("an independent stiffness measurement ... would close that non-identifiability");
   Ben cut it to the stand-scope statement that remains in Main.tex.
3. **The epsilon*-high route sentence in Design Factors** (pre-repair line 544, starts "The route
   should keep $\varepsilon^{*}$ high through the functional range"). Wrong because: the BPA force
   model delivers its largest force at low contraction and falls to zero at $\varepsilon^{*} = 1$,
   so an under-contracted actuator is at its strongest, not weakest; the sentence inverted the
   force model's direction, and driving the $(1-\varepsilon^{*})^{2}$ wrapping loss to zero by
   keeping $\varepsilon^{*}$ high would surrender exactly the force the route needs.
