# Bridge Digest — Dissertation Chapter 30-torque "Isometric Knee Torque"

Produced for the Frontiers paper bridge run (BRIEF.md contract). Everything below was
read directly from the committed dissertation sources in this session. Line references
are to the file named in each section header.

**Provenance / state check (ran this session):**
`git status --porcelain` on `chapters/30-torque.tex`, `92-AppendixA.tex`,
`93-AppendixB.tex`, `94-AppendixC.tex` returned no output (all clean);
`git log --oneline -1 0d23da03` → `0d23da03 dissertation edits`, and the last commit
touching `30-torque.tex` is that same `0d23da03`. This matches the BRIEF's "committed
state = 0d23da03" claim. Working tree = committed state for all four files.

**Files read in full:**
- `Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex` (547 lines)
- `Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/94-AppendixC.tex` (255 lines, `app:model_geometry` — the settled-Xi appendix)
- `Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/92-AppendixA.tex` (244 lines, `app:muscle_calcs`)
- `Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/93-AppendixB.tex` (scanned by grep; `app:neuron_equations` — **no torque-relevant content**: its five tables are `tab:app_neuron`, `tab:app_spk_cells`, `tab:app_spk_cal`, `tab:app_spk_gates`, `tab:app_spk_animatlab_full`, all neuron/spiking/connectome; its "knee" mentions are joint-angle traces, not torque)
- `Documentation/Reports and Papers/Dissertation/Figures/30-torque/` (inventory; note: figures live at the **Dissertation** level, one above ProofFinal — `\graphicspath{{../Figures/30-torque/}}` at 30-torque.tex:2 resolves `..` from `ProofFinal/`)

---

## 1. Section map (with publication-note overlap marking)

Publication note (30-torque.tex:5): Sections **`sec:robot_model`–`sec:improved_model`**,
**`sec:results_torque`**, and **`sec:disc_error`** overlap the manuscript
`bolen_isometric_2025` (Bolen, Burgess, Morrow, Elzein, Pang, Scharzenberger, Hunt —
7 contributors; "unpublished manuscript in preparation", not submitted). CRediT split
is spelled out there. **Dissertation-only sections** (the additional material the paper's
LONG Discussion must absorb per the BRIEF): `sec:optimization` (partially), `sec:ongoing`,
`sec:route_redesign_results`, `sec:disc_pulley`.

| # | Section (label) | Line | Overlap w/ bolen_isometric_2025 |
|---|---|---|---|
| — | `\chapter{Isometric Knee Torque}` (`ch:torque`, also `ch:results`) | 1 | — |
| 3.1 | **Methods** (no label) | 7 | container |
| 3.1.1 | Robot Leg Model and Knee Test Stands (`sec:robot_model`) | 9 | **YES** |
| 3.1.2 | Knee Torque Data Collection and Equipment (`sec:torque_equipment`) | 52 | **YES** (inside robot_model→improved_model range) |
| 3.1.3 | Hybrid Torque Calculation Method (`sec:hybrid_method`) | 63 | **YES** |
| 3.1.4 | Improved Torque Model (`sec:improved_model`) | 80 | **YES** |
| 3.1.4a | — Length Offset and System Stiffness (`sec:offset_stiffness`) | 86 | **YES** (subsubsection of improved_model) |
| 3.1.4b | — BPA Fictive Contraction (`sec:fictive_contraction`) | 276 | **YES** |
| 3.1.4c | — Multiobjective Identification Algorithm (`sec:identification`) | 299 | **YES** |
| 3.1.5 | Optimization Algorithms and Cost Functions (`sec:optimization`) | 313 | NO (not named in the note; identifies χ + defines route-redesign cost/pipeline) |
| 3.2 | **Results** (`sec:results_torque`) | 332 | **YES** |
| 3.2.1 | Flexor on Pin Joint | 336 | **YES** |
| 3.2.2 | Extensor on Pin Joint | 385 | **YES** |
| 3.2.3 | Extensor on Biomimetic Knee | 423 | **YES** |
| 3.2.4 | Flexor on Biomimetic Knee | 448 | **YES** |
| 3.3 | **Discussion** (no label) | 478 | container |
| 3.3.1 | Sources of Error in Isometric Torque Generation (`sec:disc_error`) | 480 | **YES** |
| 3.3.2 | Follow-Up Identification and Route Redesign (`sec:ongoing`) | 500 | NO — dissertation-only, carries the adopted Xi0–Xi3 + redesign |
| 3.3.2a | — Redesigned-Route Results on the Gait2392 Robot-Body Model (`sec:route_redesign_results`) | 506 | NO |
| 3.3.3 | Mono- and Bi-Pulley Tendon Routing (`sec:disc_pulley`) | 540 | NO |

External cross-refs the chapter leans on: original rigid-body torque model = `sec:bg_knee`
(background chapter 15, 30-torque.tex:334); improved BPA force characterization =
`ch:muscles` (chapter 20); force-model equations `eq:Fstar`/`eq:Fmax_3D` live in the
muscles chapter (referenced from 92-AppendixA.tex:190, 202).

---

## 2. Equations (27 labeled + 1 unlabeled inline)

All in 30-torque.tex unless noted.

**Kinematics / moment arm (Methods 3.1.1, 3.1.3):**
- `eq:momentarm_def` (:37) — scalar moment arm r = ‖d × û‖ (perpendicular distance joint-axis→line of action).
- `eq:project_ma` (:44) — Young-method projected scalar moment arm r_k about ẑ_k from the free (frame-crossing) muscle segment.
- `eq:M_vec` (:68) — vector moment m = d × F (classical cross product).
- `eq:M_scalar` (:74) — scalar knee torque M_z = r·F; basis of the hybrid method.

**Length offset (3.1.4a):**
- `eq:Lm` (:90) — corrected muscle length l_m = l_LMT − 2·l_fitting − χ0 (constant length offset).

**Bracket-frame compliance chain (3.1.4a):**
- `eq:bracket_frame_transform` (:103) — two-rotation bracket frame R_sb = R_z(θz)R_y(θy)^T and homogeneous T_sb at p_br (first rotation toward sagittal projection, second pitch onto the full 3-D attachment point).
- `eq:bktstiffness_general` (:134) — general diagonal bracket stiffness K_br = diag(k_x, k_y, k_z).
- `eq:bktstiffness` (:143) — identified y-symmetric form diag(χ1, χ2, χ1); the k_x=k_z symmetry is what makes the frame convention invariant (argument at :139).
- `eq:bktcompliance` (:155) — C_br = K_br⁻¹.
- `eq:bktforcevector` (:162) — F_br = F·û (force in bracket frame by magnitude × unit direction).
- `eq:bktvectorhooke` (:169) — Hooke's law F_br = K_br·ε_br.
- `eq:bktdeflection` (:175) — ε_br = C_br·F_br.
- `eq:bktdeflectionforceunit` (:182) — ε_br = F·K_br⁻¹·û.
- `eq:bktcomponentdeflection` (:189) — component form [F·u_x/χ1, F·u_y/χ2, F·u_z/χ1].
- `eq:bktprojecteddeflection` (:196) — scalar deflection along line of action δ_br = ûᵀ·ε_br.
- `eq:bktprojecteddeflection2` (:203) — δ_br = ûᵀ·K_br⁻¹·F_br.
- `eq:bktprojectedcompliance` (:210) — δ_br = F·ûᵀ·K_br⁻¹·û.
- `eq:keff` (:217) — scalar stiffness along the force path k_û,br = (ûᵀK_br⁻¹û)⁻¹.
- `eq:hookereduced` (:224) — reduced Hooke F = k_eq·δ.
- `eq:equivalentstiffness` (:230) — series k_eq = (1/k_û,br + 1/k_ten)⁻¹ (bracket + artificial tendon).
- `eq:activelengthcorrected` (:238) — compliance-corrected active length l_m(δ) = l_m(0) − δ.
- `eq:complianceforcebalance` (:246) — implicit force balance F_max·F_BPA(ε*(δ), P*) − k_eq·δ = 0; solved by bracketed `fzero` to machine-precision bracket collapse, with two δ=0 safeguards (:250).
- `eq:bracket_series_projection` (:256) — two-bracket generalization: per-bracket 3-D deflection η_br,j = Ση_i,j ê_i, its projection δ_br,j = û_jᵀη_br,j (j ∈ {1,2}), and δ_series = δ_br,1 + δ_tendon + δ_br,2.

**Fictive strain / wrapping loss (3.1.4b) — the ε_620 framing:**
- `eq:deltaL` (:279) — wrapping-length loss Δl = χ3·R·|θ_wrap − θ_k|·ζ² for θ_k ≤ θ_wrap (0 otherwise); R = wrapping radius, θ_wrap from solid modeling.
- `eq:comp` (:288) — ζ = 1 − ε*, the additive complement to **relative strain ε***; loss grows quadratically toward resting length, vanishes at full contraction.
- `eq:x3` (:295) — combined corrected length l_m = l_LMT − 2·l_fitting − χ0 − Δl.
- The chapter's prose definition of fictive strain is at :278 ("We call this excess fictive strain: the additional contraction the force model must see to reproduce the observed force, beyond the actuator's physical contraction" — beam-like behavior of the BPA bent over the joint, wrapped arc stops contributing usable contraction).
- **Where ε_620 itself is defined:** 92-AppendixA.tex:187 — "The maximum free-load contraction at 620 kPa is ε_620 = (l_0 − l_620)/l_0, where l_620 is the BPA length at 620 kPa with no external load (code variable `KMAX`; `MonoPamDataExplicit.m` line 278). The relative strain is ε\* = ε/ε_620." So the chapter's ε\* (hence ζ = 1 − ε\*) is normalized by the unloaded 620 kPa contraction — that is the BRIEF's "fictive-strain via epsilon_620" framing. Pressure is likewise normalized at 620 kPa: P\* = P/620 kPa (92-AppendixA.tex:197). Raw strain definition including χ0: `eq:app_strain` (92-AppendixA.tex:182).

**Optimization cost (3.1.5):**
- `eq:opt_cost` (:322) — route-redesign cost J = w_w·J_worst + w_t·J_torque + Σ_a J_a + J_geom + J_len (w_w = 1e5, w_t = 1e3).

**Unlabeled inline relation (Discussion, :493–498):** m·ε̈\* + c(ε̇\*)·ε̇\* + k(ε\*)·ε\* = F with
k(ε\*, P\*) = ∂F/∂ε\*|_{P\*} = −c0·c1·e^{−c1·ε\*} − 2·c2·ε\*·P\*·e^{−c2(ε\*)²}; setting m and F to
zero gives c(ε̇\*) = k(ε\*)·ε\*/ε̇\*. Presented as an aside on damping interpretation.

**Torque-relevant appendix equations (for the paper's Methods/appendix):**
- 92-AppendixA.tex: `eq:app_strain` (:182, strain with χ0), ε_620/ε\* (:187), `eq:app_fstar` (:191, closed-form F\* = c0(e^{−c1·ε\*} − 1) + P\*e^{−c2(ε\*)²}, zero-clamped), `eq:app_perpfoot` (:209, perpendicular-foot moment arm — equivalent to eq:project_ma), `eq:app_kbr` (:218, K_br = diag(χ1, χ2, χ1) + k_û,br), `eq:app_tendon` (:227, k_ten = n·A_c·E/L_ten), `eq:app_forcebalance` (:234, fzero balance on [0, l_m − l_620]).
- 94-AppendixC.tex: `eq:app_forcedir` (:10, route force direction from cross-point row), `eq:app_Tk` (:96, per-frame knee transform), `eq:app_bracketangles` (:111, yaw/pitch from bracket geometry), `eq:app_bracketframe` (:120, bracket pose), `eq:app_compliance` (:148, scalar path compliance + series sum).

---

## 3. Figures

### 3.1 Figures in the chapter (10, all files verified present in Figures/30-torque/)

| Label | File | Line | Caption gist / evidences |
|---|---|---|---|
| `fig:testJigs` | `testJigs_reviewed.pdf` | 15–20 | (A) CAD assembly of revolute knee in test stand; (B) pinned robot knee in apparatus, θk = −30°, components labeled; (C) photo of extensor isometric test at θk = −120° showing compressed BPA shape. → hardware/setup evidence. |
| `fig:ICR` | `KneeICR.eps` | 26–31 | ICR migration of the biomimetic knee vs knee angle (measured points −150°→+10°), overlaid on knee model at 0/−60/−120°. → biomimetic joint kinematics evidence. |
| `fig:xiFrameTransform` | `xiFrameTransform.pdf` | 114–124 | Two-rotation construction of the bracket frame (θz toward sagittal projection, then −θy onto full point vector). Has authored ALT TEXT (:117–121). → methods figure for the compliance frame. |
| `fig:bktFrame` | `bktFrame.pdf` | 126–132 | Bracket frame orientations for (A) pinned flexor, (B) pinned extensor, (C) biomimetic flexor+extensor 10 mm, (D) biomimetic flexor 20 mm. → per-configuration geometry. |
| `fig:xiProjection` | `xiProjection.pdf` | 265–270 | Two-bracket 3-D deflections vs their scalar projections along the BPA force path (F2 = −F1); CVD-safe palette noted in caption. → two-bracket series-compliance figure. |
| `fig:FlxPin_group` | `FlxPin_group.eps` (sideways) | 342–347 | Pinned-knee flexor: measured (blue) vs original model (black dashed) vs hybrid (yellow) vs optimized (magenta); (A/B) pre-optimization 48.5 & 45.7 cm, (C/D) with optimization results. → the headline flexor identification figure. |
| `fig:ExtPin_group` | `ExtPin_group.eps` | 389–394 | Pinned-knee extensor, same color code; (A) 48.0 cm, (B) 45.7 cm, (C) 41.5 cm. → extensor identification incl. the χ3 wrapping-loss effect. |
| `fig:Ext10mm_52cm` | `Ext10mm_52cm.eps` | 432–437 | Biomimetic knee, extensor 10 mm l0 = 51.8 cm: measured vs original vs optimized vs hybrid. → cross-configuration validation of extensor terms. |
| `fig:Ext_passthru` | `KneeExt_10mm_passthru.png` | 439–444 | Photo-like view of the hard-to-model configuration: BPA passing between femoral condyles close to the ICR. → modeling-difficulty evidence (used in Discussion :486). |
| `fig:FlxGrp` | `FlxGrp.eps` | 453–458 | Biomimetic flexor: (A) ⌀10 mm l0 = 41.5 cm + 21 mm tendon (measured/optimized/char.-only); (B) ⌀20 mm l0 = 42.0 cm + 15 mm tendon at 620 kPa with human max-isometric torque as black dotted line; sub-620 pressures (560/421/325/281 kPa) as pink dots (femur ultimate strength limited P). → human-torque comparison evidence. |

### 3.2 Inventory of `Figures/30-torque/` (committed == on-disk, 22 files; `git ls-files` + `find` both run)

Referenced by the chapter (10): `testJigs_reviewed.pdf`, `KneeICR.eps`, `xiFrameTransform.pdf`,
`bktFrame.pdf`, `xiProjection.pdf`, `FlxPin_group.eps`, `ExtPin_group.eps`, `Ext10mm_52cm.eps`,
`KneeExt_10mm_passthru.png`, `FlxGrp.eps`.

Not referenced (12): epstopdf byproducts `FlxPin_group-eps-converted-to.pdf`,
`ExtPin_group-eps-converted-to.pdf`, `Ext10mm_52cm-eps-converted-to.pdf`, `FlxGrp-eps-converted-to.pdf`,
`KneeICR-eps-converted-to.pdf`; source/editable variants `xiFrameTransform.png`, `xiProjection.png`,
`xiProjection.svg`, `xiProjection.pptx`; `Deprecated/` holds `testJigs.pdf`,
`testJigs_alternate1.pdf`, `paired_coefficient_tables_2026-09-29.jpg`. No `Components/`
subfolder exists for this chapter (only `Deprecated/`).

### 3.3 Appendix-C figures (route redesign of record)

- `fig:app_route_geometry` — `route_geometry_ext20mm_20260925.pdf` (94-AppendixC.tex:240–246): eight-panel femur-frame view of the optimized 9-point extensor route at the elimination-event poses.
- `fig:app_bone_mesh` — `bone_mesh_static_extensor20mm.png` (94-AppendixC.tex:248–254): optimized route on OpenSim bone point clouds at θk = 0°, total path length 0.784 m.

---

## 4. Tables

### In chapter 30-torque.tex (5)
1. `tab:OptCoeff` (:350, sidewaystable) — **initial single-bracket identification** coefficients: Flexor X0 = 0.0119 m, X1 = 1.02E+05 N/m, X2 = 1.06E+04 N/m; Extensor X0 = −0.0109 m, X1/X2 same pair, X3 = 0.123.
2. `tab:FlxPinGof` (:364) — pinned-flexor per-test GoF, improved-BPA-char-alone vs corrected (adopted two-stage terms); 47.9 cm held out. Training-mean RMSE 8.122 → 1.557 N·m; validation 2.287 N·m.
3. `tab:ExtPinGof` (:400) — pinned-extensor per-test GoF at adopted terms (stiffness pair locked to flexor); fitting-mean RMSE 2.596 → 1.159 N·m; validation (45.7 cm) 2.178 → 0.490 N·m.
4. `tab:TorqueFit` (:462) — cross-configuration summary (5 rows, at the *initial* identification values).
5. `tab:bipulley_batch` (:520) — bi-pulley batch on five muscle sites: soleus 0.9968 / tibialis anterior 0.9288 (both flag 1, nonpositive margin), rectus femoris 3.423 (worst constraint 0.2423), medial gastrocnemius + biceps femoris short head = 101 infeasibility sentinel (flag −2).

### Appendix C (94-AppendixC.tex) — the settled-Xi tables (9)
1. `tab:app_flxpin_tests` (:17) — 5 pinned-flexor tests: kf(1) 48 cm l0 0.485 / kf(2) 46 cm 0.457 / kf(3) 47 cm 0.479 (held-out validation, +5.3° encoder correction) / kf(4) 40 cm-tendon 0.406+0.035 / kf(5) 41 cm 0.410; pressures 605.7–623.4 kPa; Fm ≈ 436.9–443.4 N; all ⌀10 mm.
2. `tab:app_ext_seed` (:33, sideways) — biomimetic extensor 20 mm 9-point seed route (p1–p9 femur/tibia frames) + release schedule from `Vas_Pam_20mm_Result_20260925.mat`: p5 −89.8°, p4 −55.7°, p6 −47.8°, p7 −20.2°, p3 +6.1°, p8 +7.4°.
3. `tab:app_extpin_route` (:53) — pinned-extensor test-1 route at θk = −125° (7 points, mm).
4. `tab:app_extpin_tests` (:71) — 9 extensor tests ke(1)–ke(9) with l0/tendon/P/Fm; validation pool = {1,2,5,6,7,8}; tests 3, 4, 9 excluded. (Labels round: "42 cm" is l0 = 0.415 m, "43 cm-tendon" is 0.432 m.)
5. `tab:app_bracket_points` (:128) — bracket reference points (mm): Pbri (−48.11, −107.81, 13.80) pinned-flexor insertion; Pbr2 (−52.61, 0, 75.06) pinned-flexor origin-side; Pbr (−3.84, −46.44, 62.50) pinned extensor rib midpoint; (8.38, 20.75, 25.10) biomimetic ext 10 mm cantilever centroid; (9.48, −33.38, 33.90) biomimetic flx 20 mm bolt-pattern centroid; (−19.00, 22.00, 27.60) biomimetic flx 10 mm cantilever centroid.
6. `tab:app_stiffness_arrays` (:156) — as-executed stiffness arrays by configuration: pinned-flexor insertion [X1,X2,X1]; pinned-flexor origin [X1,X1,X2]; pinned extensor 2trans [X1,X2,X2]; biomimetic extensor [X1,X2,X1]; biomimetic flexor [X1,X2,X2]; legacy single-bracket [X1,X2,X1]. With code-line provenance.
7. `tab:app_opt_bounds` (:175) — adopted-identification bounds: flexor X0 ∈ [0,20] mm, X1 ∈ [3e4, 1e6] N/m, X2 ∈ [5e3, 2e4] N/m; extensor refit X0 ∈ [−20,+5] mm, X3 ∈ [0,1], X1/X2 locked to flexor front. Initial population X0 ∈ [5,15] mm, X1 ∈ [5e4,5e5], X2 ∈ [7e3,1.5e4].
8. **`tab:app_xi_values` (:194) — the settled values of record:** Flexor χ0 = +3.94 mm, χ1 = 3.998e4 N/m, χ2 = 1.473e4 N/m; Extensor χ0 = −6.38 mm, χ3 = 0.158 (same locked χ1/χ2). Provenance (:192, :207): flexor = pooled-front row 77 of `minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat` (212 members); extensor = row 32 of `minimizeExt10mmX3_results_20260920_flxr77.mat`; row-77 pair adopted because both redesigns then meet the 5% overshoot simultaneously (flexor +5.00%, extensor +5.043% at −20.2° in the adoption logs; refined redesigns +5.01%/+5.02%).
9. `tab:app_gof` (:211) — per-test GoF **at the adopted values** (flexor 5 tests: RMSE 1.08–2.29, FVU 0.03–0.17; extensor pool 6 tests with rigid-RMSE baseline: corrected 0.49–1.52 vs rigid 0.95–3.62).

### Appendix A (92-AppendixA.tex), torque-relevant (1 table + key text)
- `tab:force-r2` (:98) — MuJoCo vs OpenSim muscle-force R² along subject01 walk (12 muscles; 8 of 12 R² = 0.87–0.98). Cross-platform muscle-model support, not a torque table per se.
- Text of record: `sec:app_opensim_targets` (:29–52) — human torque targets are OpenSim tendon-excursion moment arms r(θk) = ∂l_MT/∂θk × max-isometric force; exported tables `OpenSim_Vasti_MomentArm.txt`, `OpenSim_Vasti_Results.txt`, `OpenSim_Bifem_MomentArm.txt`, `OpenSim_Bifem_Results.txt`; optimization loads **vastus medialis** (extensor) and **biceps femoris short head** (flexor) columns.
- Constants: l_fitting = 25.4 mm (constructor default, :185); tendon spring n = 2 strands (×6 for 20 mm assembly), A_c = 1.51e-6 m², E = 193 GPa (:231); torque values discarded for ε < −0.02 (:213); 20 mm Fmax coefficient CIs a1 ∈ [1.4745, 1.5011] N/kPa, a2 ∈ [0.0238, 0.0258] /kPa/m (:202); unused 40 mm fit c0 = 0.1224, c1 = 10.47, c2 = 2.023 (:197).

### Appendix B (93-AppendixB.tex)
No torque-relevant tables (verified by grep; all five tables are neuron/spiking/connectome).

---

## 5. Numbers of record (with chapter locations)

### 5.1 Identification parameters
- **Initial single-bracket set** (results figures + `tab:OptCoeff`, 30-torque.tex:348, :357–358, :396):
  flexor χ0 = 11.9 mm, χ1 = 102 kN/m, χ2 = 10.6 kN/m; extensor χ0 = −10.9 mm, χ3 = 0.123
  (stiffness pair shared). Text notes extensor χ0 ≈ equal magnitude, opposite sign to flexor (:396).
- **Adopted two-stage set of record** (30-torque.tex:502; `tab:app_xi_values` 94-AppendixC.tex:201–202):
  flexor χ0 = +3.94 mm, χ1 = 3.998e4 N/m, χ2 = 1.473e4 N/m; extensor χ0 = −6.38 mm, χ3 = 0.158 on the same locked pair.
- Stiffness pair is **not uniquely determined** by the pinned-flexor data (family of near-equal pairs);
  adopted pair chosen by the external requirement that both redesign configurations meet the 5% overshoot simultaneously (30-torque.tex:502; 94-AppendixC.tex:192).

### 5.2 Flexor pinned-knee GoF
- At initial values (30-torque.tex:348, `tab:TorqueFit`:469–470): 48.5 cm fit RMSE 1.616 / FVU 0.0958 / Max 3.2527; 45.7 cm validation RMSE 2.1214 / FVU 0.0878 / Max 5.1690.
- At adopted values (`tab:FlxPinGof`:373–380): per-test corrected RMSE 1.084–2.287, FVU 0.0277–0.1718; training-mean 8.122 → 1.557 N·m; held-out 47.9 cm 2.287 N·m / FVU 0.1641.

### 5.3 Extensor pinned-knee GoF (adopted values, `tab:ExtPinGof`:409–417)
- Corrected RMSE 0.490–1.525 per test; fitting-mean 2.596 → 1.159 N·m; validation 45.7 cm 0.490 / FVU 0.0535 / Max 1.206; FVU < 1 on all but the 43.2 cm tendon test (1.3270).

### 5.4 Biomimetic-knee GoF
- Extensor 10 mm, l0 = 51.8 cm — **two of-record pairs, different parameter sets (do not merge):**
  - at the *initial* extensor terms (30-torque.tex:430, `tab:TorqueFit`:471): RMSE 0.8920 / FVU 0.3396 / Max 1.4188;
  - at the *adopted* terms (30-torque.tex:502; 94-AppendixC.tex:209): RMSE 3.03 → 1.17 N·m, FVU 4.49 → 0.67, Max 2.21 (rigid baseline 3.03/4.49/5.59).
- Flexor 10 mm l0 = 41.5 cm + 21 mm tendon @ 604 kPa (30-torque.tex:451): RMSE 0.3007 / FVU 0.0102 / Max 0.4647 N·m.
- Flexor 20 mm + 15 mm tendon @ 620 kPa (30-torque.tex:451): RMSE 1.8024 / FVU 0.0499 / Max 3.5678 N·m; single-BPA torque does not exceed human magnitude ca. 60–100° flexion; within ±20% of human mean, insufficient 78–103°.
- **Label discrepancy to flag:** text/figure say l0 = 42.0 cm (30-torque.tex:451, :457) while `tab:TorqueFit`:473 lists the same GoF triplet as "l0 = 41.5 cm". Ben's call which is right for the paper.

### 5.5 Test-set constants (Appendix C)
- Pinned flexor: 5 tests, ⌀10 mm, l0 0.406–0.485 m, P 605.7–623.4 kPa, Fm 436.9–443.4 N (`tab:app_flxpin_tests`).
- Pinned extensor: 9 tests, pool {1,2,5,6,7,8}, excluded 3/4/9 (`tab:app_extpin_tests`).
- Flexor pinned route (test 1): p1 = [−75.0, 100.0, 32.8] mm femur, p2 = [−50.1, −45.0, 32.3] mm tibia (:15).
- Wrap radii (pinned extensor): R1 = 12 mm (27°→−23°), R2 = 40 mm; normalizers 27.75/(ang1−ang2), 92.46/(ang2+120°); biomimetic bend measure w_rap = 25 mm (94-AppendixC.tex:173).
- Knee transforms: ICR spline fits x 23.30→10.47 mm, y −416.65→−395.66 mm (θk 0.17→−2.62 rad); archived flexor set constant translation [4.5, −395.8, 0] mm, sweep −120°→+10° (:95–106).

### 5.6 Optimization settings (30-torque.tex:319, :328, :326, :330; :508)
- Identification: `gamultiobj` (NSGA-II), adaptive-feasible mutation, scattered crossover 0.8, Pareto fraction 0.35; flexor two-bracket pop 150 / ≤600 gens / spread tol 4e-3; extensor refit pop 50 / ≤150 gens / tol 1e-3; normalized space (cm + log10 stiffness).
- Route redesign: scalar penalty cost `eq:opt_cost` with w_w = 1e5 (worst-angle shortfall), w_t = 1e3 (mean-squared shortfall); `surrogateopt` budgets 3,000 (extensor) / 1,000 (flexor) / 1,000+1,500 (cam variant); `patternsearch` refinement 15,000 evaluations; hard constraints via prohibitive penalty; human targets from Gait2392 (vastus medialis extensor, biceps femoris short head flexor); required torque margin 5% above the human curve; two BPAs per joint (flexor mirrored pair summed; extensor single shared route × count); redesign actuator plan: two ⌀20 mm BPAs/route, 620 kPa ceiling, max usable contraction 0.255.

### 5.7 Route-redesign results (Discussion, dissertation-only)
- Direct extensor: predicted two-BPA torque 35.45–109.4 N·m vs human-scaled requirement 20.25–50.32 N·m; zero shortfall angles; thinnest margin +5.02% at −20.2° (:510).
- Direct flexor: −81.23 to −3.864 N·m vs requirement −23.16 to −2.461 N·m (peak magnitude 23.16 N·m); min margin fraction 0.0501 vs required 0.05 (+5.01% at −107°, :504, :514); assumption: both actuators share one attachment route (mirroring = next refinement, :504).
- **Cam redesign (current extensor design)** (:512): min margin +30.83% vs required envelope, +39.77% vs plain human target (optimizer winner pre-adjustment +31.30%/+42.34%); pressure to meet envelope 269–507 kPa (VTUG floor 150 kPa); tibial wrap bend radius 96.9 mm; condyle clearance 43.1 mm at 10° extension; moment arm 59.4–102.8 mm (mean 72.4 mm).
- Mono-pulley flexor variant (:516, :542): converged collision-feasible, flag 1, objective 0.08465 (direct-route objective 74.47 — different functions, not commensurable); single tackle-line part at fixed 1:1 gain (redirect, not multiplier); peak flexion torque 81.23 → 85.74 N·m (+4.51 N·m, 5.6%); 85.56 N·m insertion-side measure; off-axis torque ≤ 0.803 N·m (<1% of peak); peaks exceed the human flexor scale 23.16 N·m by factors 3.51/3.70.
- Bi-pulley batch (5 sites, `tab:bipulley_batch` :520–534, :518, :542): all five nonpositive torque margin; soleus 0.9968, tibialis anterior 0.9288 (flag 1); rectus femoris 3.423 (worst constraint 0.2423); medial gastrocnemius + biceps femoris short head at the 101 infeasibility sentinel (worst constraint 10, flag −2); finite min pair clearances 4.260–19.14 mm vs 25 mm required (two ⌀10 mm radii + 5 mm bundle margin); BFSH site target 40.2 N·m at 804 N max isometric force.
- Adductor magnus restoration (:536): P2 restored to (−0.1357, −0.0929, 0.0591) m; flexion torque 34.28 N·m = 1.12× human 30.61 N·m; adduction 27.62 N·m = 0.82× human 33.60 N·m; regenerated mirror model passed with zero non-mirrored pairs; MuJoCo conversion still carries the stale rounded P2.
- Whole-body 27-actuator extension (:538, honest status): first pass 19/27 actuators, all-but-one chosen routes violated feasibility (exception BFLH +0.03%); restructured pass: 12/27 terminated on a route-builder concatenation error, 15 finalized, 8 constraint-feasible (semimembranosus, BFLH, TFL, soleus, lat. gastrocnemius, tibialis anterior, iliacus, psoas), of which 6 met torque at every angle (worst-angle margins +0.25% BFLH to +250.7% TFL); semimembranosus + soleus geometry-feasible but torque-short.
- First physical measurements of the redesigned extensor route (:544): ⌀20 mm, near full extension — 1.163 N·m @ 10° (load 3.7 N, 388.5 kPa), 7.849 N·m @ 11° (25.2 N, 388 kPa), 23.437 N·m @ 10° (76 N, 592 kPa); full sweep + flexor measurements pending.

### 5.8 Hardware/equipment numbers (`sec:torque_equipment`, :52–61)
80/20 1515 frame; MODERN STEP 300 kg crane scale; CALT DYLY-103 100 kg S-beam load cell + HX711; Freescale MPX5700 GP 5 V pressure sensor; Parker 20R113GC 0–120 psi + Husky HDA72200 regulators in series; Festo VTUG-10 manifold; Sparkfun BlackBoard C (Arduino Uno style); 3/16 inch Kevlar winch rope; leg geometry from bone scans of a 6 ft (1.83 m) person (:22); biomimetic knee per steele_biomimetic_2018 refs (:24).

---

## 6. How the Discussion carries the follow-up identification (Xi0–Xi3) and route-redesign results

Structure (all dissertation-only content, i.e., exactly the "ADDITIONAL RESULTS" the paper's
LONG Discussion must absorb per the BRIEF):

1. **`sec:disc_error` (:480, overlap section)** — error sources narrative: compliance in
   winch cable/tendons/brackets; χ0/χ1/χ2 introduction recap; wrapped-path approximation
   error (short moment arms amplify, citing young_analyzing_2019); BPA slip ±20 mm in Z;
   the load-dependent divergence in `fig:ExtPin_group` motivating χ3; biomimetic extensor
   validation statement; closes with the "progressing from simple to complex" framing and
   the unlabeled damping aside (:493–498).
2. **`sec:ongoing` (:500)** — the identification follow-up: declares `tab:OptCoeff` values
   to be the *initial single-bracket* identification; states the adopted two-stage values
   (+3.94 mm / 3.998e4 / 1.473e4; extensor −6.38 mm / 0.158); states the non-identifiability
   of the stiffness pair and the external selection rule (both redesigns must meet the 5%
   overshoot simultaneously); quantifies what the corrections buy (tables FlxPinGof/ExtPinGof;
   extensor fitting-mean halved 2.60 → 1.16 N·m; biomimetic extensor 3.03 → 1.17 N·m,
   FVU 4.49 → 0.67); assigns solvers (gamultiobj for identification because the Pareto front
   is needed for holdout validation; surrogate+patternsearch for redesign).
3. **`sec:route_redesign_results` (:506)** — the redesign as *additional results* on the
   Gait2392 robot-body model: shared actuator plan (two ⌀20 mm BPAs, 620 kPa, 0.255 max
   contraction); direct extensor (+5.02% thinnest margin), cam redesign (+30.83%/+39.77%,
   pressure/headroom/geometry numbers), direct flexor (+5.01%, −81.23 to −3.864 N·m vs
   −23.16 to −2.461 N·m requirement), mono-pulley (+4.51 N·m), bi-pulley batch (all failed,
   geometry-bound), adductor-magnus model restoration, whole-body 27-actuator honest status,
   and the caveat that these are optimizer-stage results not yet propagated into the
   OpenSim/MuJoCo trees.
4. **`sec:disc_pulley` (:540)** — what the pulley freedom buys (reach and screening),
   mono-pulley numbers in full, bi-pulley screening conclusion (redirects design effort to
   clearance-relaxed sites), and the first physical torque measurements of the rebuilt
   extensor route (:544).

Narrative devices the paper can reuse: "initial → adopted" two-stage arc; stiffness pair
selected by an external design requirement, not residual minimization; every redesign family
tied back to the same adopted χ set; honest-status reporting for unfinished extensions.

---

## 7. Flags / notes for the paper bridge (discrepancies — recorded, not harmonized)

1. **Biomimetic flexor 20 mm resting length:** 42.0 cm in text+figure (30-torque.tex:451, :457)
   vs 41.5 cm in `tab:TorqueFit`:473 for the same GoF triplet (1.8024/0.0499/3.5678). Paper
   must pick one — Ben's call.
2. **Two biomimetic-extensor (51.8 cm) GoF pairs** belong to two different χ sets:
   0.8920/0.3396/1.4188 at the *initial* extensor terms (:430) vs 1.17/0.67/2.21 at the
   *adopted* terms (:502; Appendix C :209). The paper should state which set every number
   uses; the BRIEF prefers dissertation-chapter values of record (adopted set for the
   Discussion's additional results).
3. **`sec:identification` array wording** (30-torque.tex:303: "χ2 in the X_br direction,
   χ1 in Y_br and Z_br" for the origins bracket) describes the older single-bracket
   identification and does not match the as-executed arrays in
   `tab:app_stiffness_arrays` (94-AppendixC.tex:163–165: insertion [X1,X2,X1], origin
   [X1,X1,X2], legacy single-bracket [X1,X2,X1]). Appendix C is the as-executed record
   (with code-line provenance); the paper should follow Appendix C.
4. **ε_620 location:** the chapter uses ε\*/ζ = 1−ε\* throughout; the literal ε_620 definition
   is 92-AppendixA.tex:187 (ε_620 = (l_0−l_620)/l_0, `KMAX`; ε\* = ε/ε_620). The paper's
   fictive-strain framing should define ε_620 once and use it, per BRIEF.
5. **Figure color conventions:** the chapter's existing result figures use blue (measured),
   yellow (hybrid), black dashed (original), magenta (optimized/uncorrected-vs-corrected
   variants). The BRIEF/PlotPrompt spec for *regenerated* figures is indigo measured /
   gold corrected / gray dashed original / c{5} corrected prediction. Regenerated figures
   will therefore look different from the dissertation versions by design.
6. **Figure QC rule** (BRIEF, from figure_qc.md): every caption states "negative values
   denote knee flexion" — the current chapter captions do NOT carry that sentence; new
   figures/captions for the paper must.
7. **Author-count mismatch** already flagged by the BRIEF: Main.tex lists 4 authors vs 7
   contributors in the publication note (30-torque.tex:5) — Ben's call, unchanged.
8. tab:OptCoeff uses X0–X3 symbols; the chapter body uses χ0–χ3; Appendix C mixes X_i (arrays,
   code-facing) and χ_i (values). The paper should standardize on χ.
9. The chapter's `tab:TorqueFit` rows 1–2 (pin-joint flexor at initial values) and figure
   panels C/D of `fig:FlxPin_group` are *initial-identification* artifacts; the per-test
   tables of record at the adopted values are `tab:FlxPinGof`/`tab:ExtPinGof`/`tab:app_gof`.
