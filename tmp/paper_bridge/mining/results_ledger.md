# Results Ledger — Isometric Torque Paper (Knee_Torque_Test)

Numbers ledger for the manuscript "Design Factors for Isometric Robotic Knee Torque Using
Braided Pneumatic Actuators". Built from a read-only inventory of `Testing_Data/2022_02_Festo/`
(minimize* family + `Dig_out/`), the summary tabs of
`Testing_Data/2026_06_Festo/Results_table_10mm_pinned.xlsx` and `Results_table_20mm.xlsx`
(read with `D:/Anaconda/python.exe` + openpyxl 3.0.10, READ-ONLY), and
`Code/Matlab/Mesh_Optimization/Results/` (route-redesign outcomes).

**Values of record = the committed dissertation chapter
`Documentation/.../Dissertation/ProofFinal/chapters/30-torque.tex` (committed state 0d23da03)
and its Appendix C (`94-AppendixC.tex`).** Disagreements are flagged, never averaged.

Confidence codes:
- **DR-verified** — dissertation-record value, and I reproduced/confirmed it against the named
  artifact (mat/log/CSV/xlsx) or by running the producing code during this session.
- **DR** — dissertation-record value; artifact exists but was not individually re-verified here.
- **AO** — artifact-only (no dissertation record yet); needs a chapter/paper home before use.
- **FLAG** — disagreement or provenance caveat (see §Flags).

Verification tooling for this session lives in `tmp/paper_bridge/mining/dump_*.py`,
`read_xlsx.py`, `csv_ranges.py` (all read-only).

---

## 1. Adopted correction terms (two-stage identification of record)

| Value | Unit | Source (file / row / field) | Paper claim it supports | Confidence |
|---|---|---|---|---|
| Flexor χ0 = **+3.94** (3.9374) | mm | `30-torque.tex:502` + `94-AppendixC.tex:201` (Table `tab:app_xi_values`); artifact: `Testing_Data/2022_02_Festo/minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat`, `results_sort_actual` MATLAB row 77 (scipy idx 76) = `[…, 3.937370e-03, 3.997873e+04, 1.473375e+04, …]` | Constant length offset correcting measured-vs-predicted muscle length (Discussion additional results) | DR-verified (mat row 77 read this session) |
| Flexor χ1 = **3.998e4** (39,978.7) | N/m | same mat row 77 (col 5) + appendix Table `tab:app_xi_values` | Effective axial stiffness (bracket-frame x,z) | DR-verified |
| Flexor χ2 = **1.473e4** (14,733.8) | N/m | same mat row 77 (col 6) | Effective bending stiffness (bracket-frame y) | DR-verified |
| Extensor χ0 = **−6.38** (−6.3780) | mm | `30-torque.tex:502`, appendix Table `tab:app_xi_values`; artifact: `minimizeExt10mmX3_results_20260920_flxr77.mat`, `x_actual[32] = [-6.378022e-03, 3.997873e+04, 1.473375e+04, 1.582887e-01]`; also `Dig_out/Opt_flxr77_20260920.log:9` "EXT ctx Xi = [-0.00637802 …]" | Extensor-specific length offset (opposite sign to flexor) | DR-verified (mat row 32 read this session) |
| Extensor χ3 = **0.158** (0.158289) | unitless | same mat/log as above | Wrapping-loss (fictive-strain) scale factor | DR-verified |
| Lock provenance: extensor refit locks (χ1, χ2) to flexor front row 77; `xi1lock=39978.73`, `xi2lock=14733.75`, `flxRow=77`, `LOCKSRC` documented | — | `minimizeExt10mmX3_results_20260920_flxr77.mat` fields | "single stiffness pair shared across flexor and extensor configurations" (advisor requirement) | DR-verified |
| Adoption criterion: both redesigns meet the 5% overshoot simultaneously — flexor redesign min margin **+5.00%**, locked extensor **+5.043% at −20.202°** | % | `Dig_out/frommat_flx_20260920.log:118` (+0.050000) and `Dig_out/Opt_flxr77_20260920.log:15721` ("+5.043% at -20.202 deg"); narrative in `94-AppendixC.tex:192` | How the adopted pair was selected from the flat Pareto front | DR-verified |
| Search bounds (flexor): χ0 ∈ [0, 20] mm, χ1 ∈ [3e4, 1e6], χ2 ∈ [5e3, 2e4] N/m; (ext refit): χ0 ∈ [−20, +5] mm, χ3 ∈ [0, 1] | — | `94-AppendixC.tex:182-183`; code `minimizeFlxPin10mm_2brk.m:82-83` / `minimizeExt10mmX3.m:67-68` (bounds seen: lb Xi0 = −2 cm, locked pass Xi0 [−2, 0.5] cm) | Methods: optimizer parameterization | DR-verified |
| Initial population (flexor): χ0 ∈ [5,15] mm, χ1 ∈ [5e4, 5e5], χ2 ∈ [7e3, 1.5e4] N/m | — | `94-AppendixC.tex:176` caption | Methods: seeding of the identification | DR |

Superseded intermediate sets (do **not** quote in the paper, but know they exist when reading
archives): pk107lock adopted extensor `Xi0=−0.0101197, Xi1=43535.7, Xi2=17014.3, Xi3=0.620932`
(`Dig_out/xi_audit_tmp/confirmer_reeval_bio.log:4`); flexor row-107 alternative
`[0.00885741, 56237.8, 18542.3]` (confirmer log). Both were re-evaluated and rejected in favor
of the row-77 pair.

## 2. Initial single-bracket identification (reported with the results figures, Table `tab:OptCoeff`)

| Value | Unit | Source | Paper claim | Confidence |
|---|---|---|---|---|
| Flexor X0 = 0.0119 m, X1 = 1.02E+05 N/m, X2 = 1.06E+04 N/m | m, N/m | `30-torque.tex:357` (Table `tab:OptCoeff`) | First-generation correction terms (Results narrative) | DR |
| Extensor X0 = −0.0109 m, X3 = 0.123 (stiffness pair locked to flexor values) | m, unitless | `30-torque.tex:358` | Extensor offset ≈ flexor magnitude, opposite sign; wrapping factor | DR |
| GoF: flexor pin 48.5 cm (fit) RMSE 1.616, FVU 0.0958, Max 3.2527; 45.7 cm (validation) RMSE 2.1214, FVU 0.0878, Max 5.1690 | N·m, –, N·m | `30-torque.tex:348` + Table `tab:TorqueFit` rows 1-2 | Initial-identification prediction quality | DR (FLAG-4: artifacts overwritten, see §Flags) |
| Extensor bio 51.8 cm (initial set): RMSE 0.8920, FVU 0.3396, Max 1.4188 | N·m, –, N·m | `30-torque.tex:430` + Table `tab:TorqueFit` row 3 | Validation of extensor terms on the biomimetic knee | DR |
| Flexor bio 10 mm 41.5 cm: RMSE 0.3007, FVU 0.0102, Max 0.4647 (at 604 kPa) | N·m | `30-torque.tex:451`, Table `tab:TorqueFit` row 4 | Validation of flexor terms on the biomimetic knee | DR |
| Flexor bio 20 mm (620 kPa): RMSE 1.8024, FVU 0.0499, Max 3.5678 | N·m | `30-torque.tex:451`, Table `tab:TorqueFit` row 5 | 20 mm two-BPA biomimetic flexor prediction | DR (FLAG-3: length label inconsistent) |

## 3. Per-test goodness of fit at the adopted terms (two-stage of record)

Flexor pinned (corrector: `minimizeFlxPin`, mixed convention; re-evaluated 2026-09-26 in
`Dig_out/xi_audit_tmp/confirmer_reeval_gof.log` / `supgate_reeval_20260927.log`):

| Test | RMSE (N·m) | FVU | Max (N·m) | Zero-correction RMSE/FVU/Max | Source | Claim | Confidence |
|---|---|---|---|---|---|---|---|
| 48 cm (training) | 2.165 | 0.1718 | 5.467 | 7.617 / 2.1278 / 25.531 | chapter Table `tab:FlxPinGof` (30-torque.tex:373); reeval log test 1 (2.16465/0.171825/5.4668); baseline col from `Dig_out/encoder_campaign_noT3newXi_20260908.log` | corrected ≫ improved-BPA-only prediction | DR-verified |
| 46 cm (training) | 1.368 | 0.0468 | 3.614 | 8.803 / 1.9380 / 26.264 | same, test 2 | " | DR-verified |
| 47 cm (held-out validation) | 2.287 | 0.1641 | 7.818 | 7.670 / 1.8463 / 21.650 | same, test 3 | pure held-out validation | DR-verified |
| 40 cm tendon (training) | 1.611 | 0.0427 | 5.891 | 6.373 / 0.5935 / 16.859 | same, test 4 | tendon-in-series case | DR-verified |
| 41 cm (training) | 1.084 | 0.0277 | 1.766 | 9.694 / 2.1149 / 21.076 | same, test 5 | " | DR-verified |
| Mean (training, 4 tests) | 1.557 (RMSE) | 0.0723 | 4.184 | 8.122 / 1.6935 / 22.433 | chapter `tab:FlxPinGof:379` | "training-mean RMSE from 8.12 to 1.56 N·m" (30-torque.tex:362) | DR-verified (mean recomputes from the 4 training rows) |

Extensor pinned pool (corrector: `minimizeExtX3`, 2trans default; same reeval logs, tests
{1,2,5,6,7,8} per `Dig_crossPredict.m:46`; labels use true l0):

| Test | RMSE | FVU | Max | Rigid RMSE baseline | Source | Claim | Confidence |
|---|---|---|---|---|---|---|---|
| 40.0 cm (fitting) | 1.244 | 0.1303 | 2.196 | 3.092 | chapter Table `tab:ExtPinGof:409`; reeval test 1 (1.24371/0.130342/2.19576); rigid col = appendix Table `tab:app_gof` + `baselineScores` | RMSE improves on every pool test | DR-verified |
| 40.0 cm tendon (fitting) | 1.350 | 0.8878 | 1.890 | 3.622 | reeval test 2 (1.34986/0.887751/1.88953) | " | DR-verified |
| 43.6 cm (fitting) | 0.978 | 0.1634 | 1.818 | 2.437 | reeval test 5 (0.977636/0.163404/1.81786) | " | DR-verified |
| 43.2 cm tendon (fitting) | 1.525 | 1.3270 | 2.481 | 2.878 | reeval test 6 (1.52455/1.32704/2.48131) | only test with FVU > 1 | DR-verified |
| 45.7 cm (validation) | 0.490 | 0.0535 | 1.206 | 2.178 | reeval test 7 (0.48955/0.0534948/1.20601) | held-out validation | DR-verified |
| 46.5 cm (fitting) | 0.699 | 0.1714 | 1.559 | 0.949 | reeval test 8 (0.699429/0.171366/1.55922) | " | DR-verified |
| Mean (fitting, 5 tests) | 1.159 | 0.5360 | 1.989 | 2.596 | chapter `tab:ExtPinGof:416` | "2.60 → 1.16 N·m, roughly halved" (30-torque.tex:502) | DR-verified |
| Excluded tests 3, 4, 9 (42 cm, 42 cm-tendon, 48 cm) GoF 1.299/1.289/0.784 RMSE | N·m | reeval tests 3,4,9 | not part of the fit of record | AO (context only) |

Biomimetic extensor 51.8 cm at the adopted terms:

| Value | Source | Claim | Confidence |
|---|---|---|---|
| corrected RMSE 1.17 N·m, FVU 0.67, Max 2.21 N·m; rigid baseline 3.03 / 4.49 / 5.59 | `30-torque.tex:502` + `94-AppendixC.tex:209`; artifact `Dig_out/xi_audit_tmp/confirmer_reeval_bio.log` (pick32 pure 1.1702/0.6681/2.2088; baseline 3.0322/4.4856/5.5870) | shared-pair model explains bulk of biomimetic-extensor variance ("improves RMSE from 3.03 to 1.17 and FVU 4.49 → 0.67") | DR-verified |

## 4. Test-configuration constants (values the model consumed)

From Appendix C Tables `tab:app_flxpin_tests` / `tab:app_extpin_tests` (read from
`FlxPinBPASet.mat` kf(1..5), `ExtPinBPASet.mat` ke(1..9); I did not re-open the BPASet mats —
appendix is the record):

| Set | Values | Confidence |
|---|---|---|
| Pinned flexor (⌀10 mm): kf(1) 48 cm l0=0.485 P=605.7 Fm=443.4; kf(2) 46 cm 0.457/614.2/441.4; kf(3) 47 cm 0.479/623.2/443.0 (held out, +5.3° encoder correction); kf(4) 40 cm-tendon 0.406/623.4/436.9, tendon 0.035; kf(5) 41 cm 0.410/615.4/437.3 | m, kPa, N | DR |
| Pinned extensor (⌀10 mm): ke(1) 40 cm 0.400/624.9/436.3; ke(2) 40 cm-tendon 0.400/618.4/436.3 t=0.040; ke(3) 42 cm 0.415/605.2 (excl); ke(4) 42 cm-tendon 0.415/605.2 t=0.040 (excl); ke(5) 43 cm 0.436/621.8/439.7; ke(6) 43 cm-tendon 0.432/622.9/439.3 t=0.022; ke(7) 46 cm 0.457/605.9/441.4; ke(8) 47 cm 0.465/605.9/442.0; ke(9) 48 cm 0.480/605.9 (excl) | m, kPa, N | DR |
| Pool = {1,2,5,6,7,8}; excluded = {3,4,9} | — | DR (also in `Dig_crossPredict.m:46` per appendix) |
| Flexor route: origin p1 = [−75.0, 100.0, 32.8] mm (femur), insertion p2 = [−50.1, −45.0, 32.3] mm (tibia); sweep −120°…+10° | mm, deg | DR |
| Extensor route of test 1 at −125° (7 points, p1 [30.0,−50.0,0] … p7 [42.5,−75.9,0] mm); wrap released at +35° | mm, deg | DR |

## 5. Raw measured data (xlsx summary tabs; read this session, READ-ONLY)

`Testing_Data/2026_06_Festo/Results_table_10mm_pinned.xlsx` — 13 extensor tabs + 7 flexor tabs
(⌀10 mm pinned-knee raw lab record: loads, knee angles, load-cell angles, muscle lengths, hand
moment arms, torque formulas, failure notes). Key observations:

| Item | Value / location | Confidence |
|---|---|---|
| Measured pinned-flexor torque range across tests | e.g. `FlxTest10mm_1` torque −0.11…−14.88 N·m; `FlxTest10mm_3` (46 cm) down to −17.62 N·m | AO (raw; chapter reports fits, not these point lists) |
| Measured pinned-extensor torque range | e.g. `ExtTest10mm_1` (48.5 cm RML) 4.5…6.1 N·m; `ExtTest10mm_9` up to 10.99 N·m | AO |
| Hardware-failure context for excluded tests | `ExtTest10mm_4` r18 "bracket broke, failure via material pull out"; `ExtTest10mm_7` r18 same | AO (context for exclusions) |
| kmax per actuator (10 mm) | tabs record e.g. contraction fraction 0.1475–0.1823 (r3: "0.1475", "0.1667", "0.1811") | AO |

`Results_table_20mm.xlsx` — ⌀40 mm extensor (2 tabs), ⌀20 mm extensor (2 tabs), ⌀20 mm flexor
(5 tabs incl. blank `FlxTest20mm_51cm`):

| Item | Value | Tab / cell rows | Paper claim | Confidence |
|---|---|---|---|---|
| **First physical torque points, redesigned 20 mm extensor route** — loads 3.7 / 25.2 / 76 N at knee 10° / 11° / 10°, gauge 388.5 / 388 / 592 kPa, LC angle (tibia) 72.5 / 70 / 68° | raw inputs | `ExtTest20mm_1` rows 6-10 (cols 3-5) | hardware validation of the redesigned extensor route has begun | DR-verified (chapter 30-torque.tex:544 quotes these) |
| Measured torque at those points = **1.163 / 7.849 / 23.437 N·m** | N·m | chapter 30-torque.tex:544; **verified live this session**: `matlab -batch` run of `Code/Matlab/Mesh_Optimization/Results/Knee_Extensor_20mm.m` printed `ExtTest20mm_1 Torque actual: 1.1630 7.8494 23.4367` | "twentyfold range of held torque near extension" | DR-verified (see FLAG-1 for the constants caveat) |
| Tab's own row-16 "Torque at t1" = 1.0382 / 6.9718 / 20.7581 N·m (F·sin(LC+0.83°)·d, d=292.9 mm) | N·m | `ExtTest20mm_1` r16 | — | AO — **disagrees ~12% with the chapter values** (FLAG-1); row 17 (Torque actual) is still empty in the tab |
| 20 mm flexor raw test of record | RML 415 mm, tendon 14 mm, 22 points incl. 620 kPa set (torque-actual −5.69…−29.90 N·m) and reduced-pressure set (385/451/296.6/421/281/324.8/325.58/500/560 kPa) | `FlxTest20mm_42cm (3)` rows 6-17 | source of the FlxGrp-B measured curve; caption's pink dots at 560/421/325/281 kPa match these tab pressures | AO (FLAG-3: RML/tendon labels differ across tabs and chapter) |
| 40 mm extensor context | RML 557 mm, tendons 30/90 mm; measured 19.91 / 29.42 N·m at 604 / 530 kPa vs "Human Vas_int Torque" 30.557 / 48 N·m; "Leg assembly started breaking when load cell was reading 124 N" | `ExtTest40mm_1`, `ExtTest40mm_2` | 40 mm scale-up context (probably out of paper scope) | AO |
| Pending test | `FlxTest20mm_51cm` tab is blank | — | "the full sweep, and the corresponding flexor measurements, are pending" (30-torque.tex:544) | DR-verified (blank confirmed this session) |

## 6. Route redesign (Mesh_Optimization outcomes; Discussion additional results)

Shared actuator plan (chapter `sec:route_redesign_results`): two ⌀20 mm BPAs per route,
620 kPa ceiling, maximum usable contraction KMAX = 0.255 — KMAX 0.255 verified in
`Dig_out/frommat_flx_20260920.log:66` and redesign mats.

### 6.1 Direct extensor route (superseded by cam, still quoted)

| Value | Unit | Source | Claim | Confidence |
|---|---|---|---|---|
| Predicted two-BPA torque span **35.45…109.4 N·m** | N·m | `Vas_Pam_20mm_Result_20260925.mat` `predBest.TorqueZ` min 35.4525 / max 109.3595 (read this session); chapter 30-torque.tex:510 | exceeds the human envelope at every angle | DR-verified |
| Human-scaled requirement span **20.25…50.32 N·m** | N·m | `Ext_20mm_pressure_for_human_torque.csv` `human_target_Nm` min 20.248 / max 50.321 (computed this session); chapter :510 | " | DR-verified |
| Minimum margin **+5.02%** at worst angle **−20.2°** | % | mat scalar `minMargin = 0.0502493`; `iWorst = 77`; release schedule `routeInfo.eliminatedAngleD` (p7 at −20.202°); chapter :510 | overshoot consumed exactly as designed | DR-verified |
| Optimized endpoints p1 = [40.0, 137.0, 0] mm (femur), pEnd = [34.4, −197.1, 0] mm (tibia), rest 827.7 mm, tendon 36.8 mm | mm | `Vas_Pam_20mm_Result_20260925.mat` `xBest = [0.04, 0.137, 0, 0.03442, −0.1971, 0, 0.82774, 0.03684]`; appendix `94-AppendixC.tex:236` | design definition of the extensor route of record | DR-verified |
| Release schedule p5 −89.8°, p4 −55.7°, p6 −47.8°, p7 −20.2°, p3 +6.1°, p8 +7.4° | deg | mat `routeInfo.eliminatedAngleD = […, −55.66, −89.80, −47.78, −20.20, 7.37]`; appendix Table `tab:app_ext_seed` caption | wrap-release schedule of the route-elimination rule | DR-verified |
| Direct-route moment arm span 44.9…61.2 mm | mm | mat `predBest.momentArm` (min 0.0449 / max 0.0612) | not currently quoted in the chapter (AO) — distinct from the CAM values below | AO |
| Objective fBest 0.109 (exitflagG 0 surrogate, exitflagP 1) | — | mat scalars `fBest = 0.1090998` | internal consistency (not for paper text) | DR-verified |

### 6.2 Cam redesign (current extensor design of record)

| Value | Unit | Source | Claim | Confidence |
|---|---|---|---|---|
| Min torque margin **+30.83%** vs required envelope; **+39.77%** vs plain human target | % | chapter `30-torque.tex:512`; artifact `Vas_Pam_20mm_Result_cam_20261014_0416.mat` (verified by the supervisor MATLAB recompute recorded in `Documentation/.../Notes/revision_campaign_20261015/CAMPAIGN_REPORT.md:173-176`; I read `XiUsed = [−6.37802e-03, 3.99787e+04, 1.47338e+04, 0.158289]` and `xFullBest` in the mat, but did not recompute the margins myself) | cam converts razor-thin margin into design headroom | DR (margins) / DR-verified (Xi + geometry) |
| Pre-adjustment winner margins +31.30% / +42.34% | % | chapter :512; campaign report attributes to run-8 winner log (tmp, not re-checked) | honest reporting of the adjustment step | DR |
| Pressure to meet envelope **269–507 kPa** | kPa | chapter :512 (campaign-report-verified vs mat); VTUG floor 150 kPa | headroom vs regulated supply floor | DR |
| Tibial wrap bend radius **96.9 mm**; condyle clearance **43.1 mm at 10°** extension; moment arm **59.4–102.8 mm, mean 72.4 mm** | mm | chapter :512 (same verification path) | headroom not bought with routing risk | DR |
| Cam geometry: xFixed = Vas0925 winner; cam radius param 0.09249; refPts 3; tendonRef 36.84 mm | m | mat `xFullBest`, `xBestCam`, `refPts`, `tendonRef` (read this session) | design definition | DR-verified |

### 6.3 Flexor route (direct, of record)

| Value | Unit | Source | Claim | Confidence |
|---|---|---|---|---|
| Predicted knee-axis torque span **−81.23…−3.864 N·m** | N·m | `Bifemsh_20mm_Result_20260921_2037.mat` `predBest.TorqueZ` min −81.2343 / max −3.8643 (read this session); chapter :514 | exceeds human envelope at every angle | DR-verified |
| Human-scaled requirement **−23.16…−2.461 N·m** (peak magnitude 23.16) | N·m | same mat `humanTorque` min −23.1587 / max −2.4612; also `Flx_20mm_pressure_for_human_torque.csv` target span; chapter :514 | " | DR-verified |
| Minimum margin fraction **0.0501** vs required 0.05 | – | same mat `torqueMarginFraction` min 0.0501; chapter :514 | meets the 5% overshoot requirement at its thinnest | DR-verified |
| Applied Xi (incl. extensor χ3 0.158289); objective fBest 74.47; rest length 561.5 mm; tendon 25 mm; endpoints p1 [−125, 135, 72.08] mm, p2 [−18.72, −32.69, 38.96] mm | –, –, mm | mat `XiUsed = [3.93737e-03, 3.99787e+04, 1.47338e+04, 0.158289]`, `fBest = 74.4733`, `xRefined[…, 0.56149, 0.025]` | flexor redesign applies the adopted terms incl. wrapping factor; direct-route objective quoted in pulley comparison (chapter :542 "74.47") | DR-verified |
| Model assumption: both actuators share one attachment route (mirroring pending) | — | chapter :504/:514 | stated modeling assumption | DR |

### 6.4 Mono-pulley flexor variant

| Value | Unit | Source | Claim | Confidence |
|---|---|---|---|---|
| Peak predicted flexion torque **85.74 N·m** (knee-axis) / **85.56 N·m** (insertion-side) vs 81.23 direct | N·m | `Bifemsh_20mm_Result_pulley_20260926_1652.mat` `predBest.TorqueZ` min −85.7402, `TorqueInsZ` min −85.5621, `fRefined = 0.0846509`; chapter :516/:542 | redirect alone deepens peak flexion torque | DR-verified |
| Gain **4.51 N·m (5.6%)** | N·m, % | 85.7402 − 81.2343 = 4.506; 4.506/81.234 = 5.55% (recomputed this session from the two mats); chapter :542 | what the added routing freedom buys | DR-verified |
| Off-axis torque ≤ **0.803 N·m** (< 1% of peak) | N·m | pulley mat `predBest.offAxisTorque` max 0.8034; chapter :542 | screening quality of the redirect | DR-verified |
| Pulley gain fixed 1:1 (redirect, not multiplier); collision-feasible, zero bound violations, pattern-search exit flag 1; objective 0.08465 not commensurable with direct 74.47 | – | mat `fRefined`, `outputP.function`; chapter :516/:542 | honest cross-objective caveat | DR-verified |

### 6.5 Bi-pulley batch (five OpenSim muscle sites)

All values verified this session against
`Code/Matlab/Mesh_Optimization/Results/BiPulley_batch_summary.csv` (also chapter Table
`tab:bipulley_batch`):

| Site | Objective | Outcome | Clearance (mm) |
|---|---|---|---|
| soleus | 0.9968 (CSV 0.996794) | converged flag 1, nonpositive margin | 19.142 |
| tibialis anterior | 0.9288 (0.928794) | converged flag 1, nonpositive margin | 15.558 |
| rectus femoris | 3.423 (3.42299), worst constraint 0.2423 | flag −2, nonpositive margin | 7.024 |
| medial gastrocnemius | 101 sentinel, worst constraint 10 | flag −2 | Inf (none recorded) |
| biceps femoris short head | 101 sentinel, worst constraint 10 | flag −2 | 4.260 |

Supporting numbers: clearance requirement 25 mm = 2×10 mm actuator radii + 5 mm bundle margin
(chapter :518; CSV `footprintWidth_mm` 20.000 / `minRouteClearance_mm`); BFsh target 40.2 N·m at
MIF 804 N (chapter :518; DR). All five runs closed nonpositive — the reeved routes do not yet
fit these corridors.

### 6.6 Whole-body 27-actuator extension (status, not a completed result)

Chapter `30-torque.tex:538` only: first pass 19/27 actuators, all-but-one chosen routes
infeasible (BF long head worst-angle margin +0.03%); restructured pass: 12 terminated on a
route-builder concatenation error, 15 finalized, 8 constraint-feasible (semimembranosus, BF
long head, TFL, soleus, lat gastrocnemius, tibialis anterior, iliacus, psoas), 6 of those meet
torque at every angle with worst-angle margins +0.25% (BF long head) to +250.7% (TFL).
**Confidence: DR; no result artifact located in `Mesh_Optimization/Results/` this session —
if the paper quotes any of these, ask Ben for the whole-body run artifacts.**

### 6.7 Adductor-magnus restoration (robot-body model update)

Chapter `30-torque.tex:536`: P2 restored to (−0.1357, −0.0929, 0.0591) m; flexion torque 34.28
N·m = 1.12× human target 30.61 N·m; adduction 27.62 N·m = 0.82× human 33.60 N·m; mirrored model
validated with zero non-mirrored pairs; MuJoCo conversion predates the restoration (one stale
path point). **Confidence: DR** (chapter-only; related table
`Solid_Models/OpenSim/Gait2392_Robotbody/bpa_actuators_27.json` exists and holds the
27-actuator MIF table — BFLH mif_n = 896 N confirmed there; the 804 N short-head entry was not
isolated this session).

## 7. Optimizer settings (Methods)

| Setting | Value | Source (verified this session) | Claim |
|---|---|---|---|
| Flexor gamultiobj: population 150, ≤600 generations, spread tolerance 4e-3, adaptive-feasible mutation, scattered crossover 0.8 | – | `minimizeFlxPin10mm_2brk.m:39-40` (`POP=150; MAXGEN=600`) and `:102-105` (MutationFcn/CrossoverFraction 0.8/FunctionTolerance 4e-3); chapter :319 | identification solver settings |
| Extensor refit: population 50, ≤150 generations, tolerance 1e-3 | – | `minimizeExt10mmX3.m:111-116` (`PopulationSize 50`, `MaxGenerations 150`, `FunctionTolerance 1e-3`) | " |
| Route redesign: surrogateopt budgets extensor 3,000 / flexor 1,000 evaluations, then patternsearch 15,000 | – | chapter :328; campaign report verified vs `Opt_run_Ext.m:29` / `Opt_run.m:79` + mat funccounts (not re-run here) | redesign pipeline budgets |
| Cam variant: 1,000 surrogate + 1,500 refinement evaluations as written | – | pulley mat `maxEvalsG = 1000`, `maxEvalsP = 15000`; chapter :328 — FLAG-5: the winning cam run actually used 2,000 surrogate evaluations via env override (campaign report §30-torque round-1) | " |
| Cost weights: w_w = 1e5 (worst-angle shortfall), w_t = 1e3 (mean-square), 5% torque margin requirement | – | chapter Eq. `eq:opt_cost` + :326-330 | redesign cost function |

## 8. Fictive-strain framing support

| Item | Source | Confidence |
|---|---|---|
| ε_620 ≡ (l0 − l620)/l0, unloaded contraction at 620 kPa; ε* = ε/ε_620 | `20-muscles.tex:22-33` (Eqs. `max_strain`, `rel_strain`) — definition of record; the force characterization values themselves live in chapter 20 / the Actuators paper (bolen_2026), outside this ledger's three locations | DR (definition located; numeric ε_620 values not inventoried here) |
| ζ = 1 − ε* (additive complement) enters χ3 via Δl = χ3 R |θ_wrap−θ_k| ζ² | `30-torque.tex:278-297` (Eqs. `deltaL`, `comp`, `x3`) | DR |
| Wrap geometry: pinned extensor arcs R1 = 12 mm (27°…−23°), R2 = 40 mm; biomimetic w_rap = 25 mm | `94-AppendixC.tex:173` | DR |

## 9. Pressure back-calculation CSVs (direct-route; `Mesh_Optimization/Results/`)

| File | Span (computed this session) | Confidence |
|---|---|---|
| `Ext_20mm_pressure_for_human_torque.csv` (100 rows, −120…+10°) | human target 20.248…50.321 N·m; required pressure 229.4…596.6 kPa; adjusted variant 246.8…504.5 kPa; all rows status "ok" | AO (single-run artifact; not in the chapter. The chapter's 269–507 kPa is the CAM design — do not mix) |
| `Flx_20mm_pressure_for_human_torque.csv` (100 rows) | human target −23.159…−2.461 N·m; required pressure 114.4…609.4 kPa | AO |

---

## Flags (disagreements / provenance caveats — for Ben, never harmonized silently)

1. **FLAG-1 — First measured 20 mm extensor torques, two variants in circulation.** The xlsx
   tab `ExtTest20mm_1` row 16 carries 1.0382 / 6.9718 / 20.7581 N·m (direct formula
   F·sin(LC+0.83°)·d, d = 292.9 mm), while the chapter quotes 1.163 / 7.849 / 23.437 N·m, which
   I verified live as the print of the Adjoint-transform section of `Knee_Extensor_20mm.m`
   (lines 581-651; `matlab -batch` printed `1.1630 7.8494 23.4367`). The section's rig
   constants `dExt1 = 292.9 mm`, `angExt1 = −90.83°` are explicitly the 40 mm test-1 defaults
   with an in-code "UPDATE after the 20 mm rig setup is measured" comment, and the tab's row 17
   ("Torque actual") is still empty. If the paper quotes the measured points, confirm the 20 mm
   rig was actually measured against those constants (or update and re-quote), and populate
   row 17.
2. **FLAG-2 — Undated mats are provenance traps.** `minimizeFlxPin10_results.mat` and
   `minimizeExt10mmX3_results.mat` (no date suffix) have been overwritten by the 2026-09-08/20
   campaigns (`TAG='_noT3newXi'`, `pick=77`; extensor `sol_actual` = the flxr77 pick). The
   2025-era initial-identification artifacts live only in the dated mats
   (`minimizeFlxPin10_results_2025*.mat`, `minimizeExtPin10_results_2025*.mat`), which I did
   not individually open; the chapter's Table `tab:OptCoeff` numbers (0.0119/1.02e5/1.06e4;
   −0.0109/0.123) remain DR (not re-verified against a specific dated mat this session).
3. **FLAG-3 — 20 mm flexor biomimetic length labels disagree.** Chapter Fig. `FlxGrp` caption
   says l0 = 42.0 cm with a 15 mm tendon; Table `tab:TorqueFit` row 5 says l0 = 41.5 cm; the
   xlsx tabs say RML 415 mm/tendon 14 (`FlxTest20mm_42cm (3)`), 420/15 (`(2)`), 423/17
   (`FlxTest20mm_42cm`). Pick one label of record for the paper (the tab actually carrying the
   plotted data appears to be `FlxTest20mm_42cm (3)`, RML 415 mm).
4. **FLAG-4 — Chapter mean-RMSE uses the 4-test training mean (8.122), the run log prints the
   5-test mean (8.0314).** Both are correct for their own set (`encoder_campaign_noT3newXi_
   20260908.log`); the paper should quote the 4-test training mean (8.12 → 1.56 N·m) and say
   which set it averages. Also note that log's "Performance with sol_actual (pick=1)" block is
   the pick-1 solution, NOT the adopted row-77 values — use the §3 reeval numbers above for
   the paper.
5. **FLAG-5 — Cam surrogate budget.** Chapter says the cam variant used 1,000 surrogate
   evaluations; the winning run used 2,000 via env override (campaign report §30-torque
   round-1). If the paper states budgets, state 2,000 for the cam run or avoid the number.
6. **FLAG-6 — Superseded side-chain numbers exist in Dig_out.** `Dig_FlxBio_refine_20260908.log`
   (grid values like 0.91/1.27 RMSE at Xi ≈ [0.012, 5e5, 1e4]) and the pk107lock-era extensor
   set (Xi0 −0.0101 / X3 0.621, GoF 2.195/2.351/4.727 on the bio test) are exploration-era;
   the paper must use only §1/§3 values. Cited here so nobody mines them by accident.
7. **FLAG-7 — 10 mm raw-tab lengths vs model constants.** Raw tabs record RML values like
   485/457/479/406/417 mm for the flexor tests; the appendix mat-of-record values are
   0.485/0.457/0.479/0.406/0.410 m — the 41 cm entry differs (tab 417 vs mat 410). Appendix
   (mat) is the record; do not quote tab lengths without checking against
   `FlxPinBPASet.mat`/`ExtPinBPASet.mat`.

## Session integrity note

The verification run of `Knee_Extensor_20mm.m` (needed for the live check of the
ExtTest20mm_1 measured torques) silently regenerated two tracked artifacts in
`Mesh_Optimization/Results/` — `Ext_20mm_pressure_for_human_torque.csv` (numeric-format churn,
100 lines) and `Knee_Extensor_20mm.gif` (regenerated animation). Both were restored from HEAD
immediately (`git checkout --`, verified clean afterwards). The CSV value ranges in §9 were
computed from the committed version before that run, so the ledger values are unaffected.
Lesson for the figure agents: that script writes files despite having no `save`/`saveas`/
`exportgraphics` calls — its writers are `writetable`/`imwrite`-style; grep for those before
assuming it is read-only.

## Not done / not found

- `Knee_Flexor_NegativeX_Sweep_Raw.csv` / `_ReachCorrected.csv` and
  `Knee_Flexor_Design_Sweep_Results.mat` were inventoried but their contents were not
  extracted (no chapter claim currently cites them).
- The whole-body 27-actuator redesign artifacts were not located in
  `Mesh_Optimization/Results/` (§6.6).
- The 2025-era dated identification mats behind Table `tab:OptCoeff` were not opened (FLAG-2).
- ε_620 numeric values (force characterization) belong to chapter 20 / bolen_2026 and were not
  inventoried here (definition located, §8).
