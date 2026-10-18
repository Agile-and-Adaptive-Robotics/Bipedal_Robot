# evidence_d — Appendix batch (slides_d.py)

Sources are `Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/<file>.tex:<line>`
unless an artifact path is given. Chapter/appendix numbering verified against
`ProofFinal/build_r5/main.toc` (Ch. 3 muscles, Ch. 4 torque, Ch. 6 SNS; Apps. A/B/C =
92/93/94-Appendix*.tex) and equation numbers verified in `build_r5/main.pdf` via Ghostscript
text extraction (Ch. 3: 3.1–3.6 on printed pp. 47–48; Ch. 4: 4.18/4.20/4.22/4.24–4.26 on
printed pp. 85–91). Deck-internal slide numbers are assigned at build time by
`build_deck.py` (`start` chaining); rows below use the appendix letter of each slide.

| # | Slide (assertion) | Numbers on the slide | Source |
|---|---|---|---|
| A | "The complete actuator force model" — eqs (3.1), (3.4), (3.5), (3.6) shown | Eq. numbers 3.1, 3.4, 3.5, 3.6 | chapters/20-muscles.tex:17 (eq:strain=3.1), :39-50 (eq:Fstar=3.4), :54-58 (eq:Fmax_3D=3.5), :60-64 (eq:force=3.6); implementation files: chapters/92-AppendixA.tex:190-202 (§A.3) |
| B | "Coefficient values with 95% confidence intervals" | F₆₂₀ arctan: 303.5 (300, 308), 19.03 (17.48, 20.57), R² 0.9854, RMSE 14.72; 922.4 (914.2, 930.7), 15.37 (14.75, 15.98), R² 0.9945, RMSE 23.83. F*: c₀ 0.5682, c₁ 4.254, c₂ 0.5597, R² 0.9998/0.9994; c₀ 0.2579, c₁ 6.477, c₂ 1.321, R² 0.992/0.9943. Festo-tool RMSE 189.9 N (10 mm), 668.4 N (20 mm). 40 mm fit unused: 0.1224/10.47/2.023 | chapters/20-muscles.tex:167-199 (Table 3.1 rows), :221-248 (Table 3.2 rows), :173 and :191 (Festo-tool RMSE columns); 40 mm fit: chapters/92-AppendixA.tex:197 |
| C | "Full goodness-of-fit against published BPA models" | Bolen/Sárosi/M&B RMSE-FVU-max for 4 tests: 22.9/0.04/46.0; 16.8/0.02/30.3; 35.6/0.01/114.0; 67.6/0.03/166.6; 143.8/1.39/235.8; 132.7/1.17/222.6; 132.1/0.20/313.4; 177.7/0.18/294.8; 73.8/0.37/122.9; 62.6/0.26/101.2; 72.8/0.06/137.8; 862.1/4.20/1321.1 | chapters/20-muscles.tex:369-372 (Table 3.3; FVU rounded to 2 dp from 0.0353/0.0187/0.0144/0.0258, 1.3903, 1.1655, 0.1992, 0.1785, 0.3668, 0.2590, 0.0605, 4.2022) |
| D | "From bracket stiffness matrix to solved deflection" | Eq. numbers 4.18, 4.20, 4.22; fzero bracket [0, l_m − l₆₂₀]; δ ≡ 0 if uncorrected ε* > 1; two-BPA coupled balance on [0, F₆₂₀] | chapters/30-torque.tex:215-220 (eq:keff=4.18), :228-232 (eq:equivalentstiffness=4.20), :244-248 (eq:complianceforcebalance=4.22), :250 (fzero/termination); implementation detail: chapters/92-AppendixA.tex:233-238 (§A.5) |
| E | "Wrapping-loss implementation, per contact arc" | R₁ = 12 mm, 27° to −23°; R₂ = 40 mm; normalizers 27.75/(ang₁−ang₂), 92.46/(ang₂+120°); Δl = χ₃·b(θₖ)·(1−ε*)²; wRap = 25 mm; eqs. (4.24)–(4.26) | chapters/30-torque.tex:279-297 (eq:deltaL=4.24, eq:comp=4.25, eq:x3=4.26); chapters/92-AppendixA.tex:242 (§A.6, radii/normalizers/wRap/force-strain-only); chapters/94-AppendixC.tex:173 (§C.4 constants) |
| F | "Identification campaign details and data decisions" | +5.3° encoder correction on the 47.9 cm test; 4 training tests + 47.9 cm held out; extensor pool {1,2,5,6,7,8}, tests 3/4/9 excluded; pooled-front row 77 (212 members), locked refit row 32; +5.01% / +5.02% redesign margins; flexor pop 150, ≤600 gens, spread tol 4×10⁻³; extensor pop 50, ≤150 gens, 10⁻³; NSGA-II | chapters/94-AppendixC.tex:18 (encoder, Table C-1 caption), :72 (extensor pool/exclusions, Table C-4 caption), :192 (row-77 adoption + 5% requirement + margins), :201-207 (row 77 of 212 members; flxr77 row 32), :176-188 (bounds table); chapters/30-torque.tex:301-305 (holdout + two-stage), :319 (solver settings, NSGA-II), :502-504 (adopted values, +5.01/+5.02) |
| G | "Frame convention is irrelevant for y-symmetric stiffness arrays" | invariance iff k_x = k_z; identical torque predictions to numerical precision; arch bracket 75 mm in Z, χ₂ in X_br, χ₁ in Y_br/Z_br | chapters/30-torque.tex:139 (invariance statement + R_y K R_yᵀ identity), :303 (arch-bracket axis assignment, 75 mm in Z); frames/arrays: chapters/94-AppendixC.tex:126 (§C.3 pitch-only vs two-rotation), :156-171 (Table stiffness arrays) |
| H | "Bracket reference points and adopted correction terms — recorded together" | Pbri [−48.11, −107.81, 13.80]; Pbr2 [−52.61, 0, 75.06]; Pbr [−3.84, −46.44, 62.50]; Pbr [8.38, 20.75, 25.10]; Pbr [9.48, −33.38, 33.90]; Pbr [−19.00, 22.00, 27.60] (all mm); adopted: χ₀ +3.94 mm, χ₁ 3.998e4 N/m, χ₂ 1.473e4 N/m (flexor); χ₀ −6.38 mm, χ₃ 0.158 (extensor, pair locked); row 77 / row 32 provenance | chapters/94-AppendixC.tex:128-143 (Table C-5, all six points verbatim), :194-205 (Table C-8, adopted values verbatim), :207 (row 77 / flxr77 row 32 provenance), :4 (point-recorded-with-values rationale); flexor/extensor values also at chapters/30-torque.tex:502 |
| I | "Moment arm and torque conventions" | perpendicular foot = cross-product projection of Eq. (4.2); r_k = hypot(in-plane) in tibial frame about z; discard ε < −0.02; human targets r(θₖ) = ∂l_MT/∂θₖ | chapters/30-torque.tex:37-46 (eq:momentarm_def=4.1, eq:project_ma=4.2), :48 (tibial body frame); chapters/92-AppendixA.tex:213 (§A.4: hypot, ε < −0.02, equivalence statement); human targets: chapters/92-AppendixA.tex:32-36 (§A.1 derivative definition), chapters/30-torque.tex:326 |
| J | "Synthetic nervous system: Simulink library and spiking twin" | 7 core blocks → 12 → 15 with spiking additions; base Simulink; twin transient 15.79° at 0.05 s (both), settle +0.9°; single-spike conductance 6.2×10⁻⁷ µS; 62/62 spikes; regression 6.07×10⁻⁴ mV; SolidWorks → URDF → Simscape verified, tendon modeled in CAD not yet inserted; Figs. 6.5, 6.11 | block counts: chapters/60-sns.tex:91 (seven core + twelve), :123 (fifteen); twin numbers: chapters/60-sns.tex:277, Table 6.4 rows at :286-288 (15.79/15.79, settle 43.5→44.4 = +0.9); unit checks: chapters/60-sns.tex:123 (6.2e-07 µS, 11.75 µs), :275 (62/62, 6.07e-04 mV); CAD coupling: chapters/60-sns.tex:93; figure files: Figures/60-sns/simulink_sns_library_grid.png (Fig. 6.5), simulink_knee_spiking_twin.png (Fig. 6.11) — copied as deck_assets/figs/d_simulink_library_grid.png, d_simulink_knee_spiking_twin.png |
| K | "Walker configurations of record: scores, library sweep, and pathway deletions" | tuning-ref / library-mean: s3k −186.4/−8776; s3kpruned −176.1/−6924; syn6 −200.4/−9961; w2lvar −144.5/−5808; worst refs −8841 (Case_40_motion), −6986, −10040, −5896 (subject08_Run_20002); s3k deletions: −186.4 full, Ia −25 (no cycles), all afferents −267.8, Renshaw −267.3, II −255.2, Ib −180.3, contact −183.9, commissural −174.9; drive: −201.9/−177.0/−171.2 for +0.01/−0.01/−0.02, freeze at +0.02; 15 push evals, max sway 0.123 m; syn6/w2lvar fall under ±0.01/±0.02 (tilt 39–40°); s3kpruned −0.02 → −168.2 | chapters/60-sns.tex:179-182 (Table 6.1 rows), :187 (worst-reference sentence), :202-217 (Table 6.2 rows verbatim), :224 (drive-perturbation and push-test paragraph); notes-only extras: walking-only −8816 vs running-only −8761 (:187), syn6 all-afferents gait −190.5, 10+10 cycles, knee −71.2°, 0.684 s, w2lvar Ia −273.5, weak-RG-exc −279.9 (:193) |

Notes on provenance decisions:
- A6/F: the old deck's "fold study ≈80% cross-fold duplicates" and "fresh GA seed reproduced the 4-member front" cards had no dissertation anchor and were replaced with the chapter-anchored non-uniqueness/adoption story and the recorded solver settings.
- A7/G: the old "agree to ≈1×10⁻¹⁷" claim is not in the dissertation; the chapter states "identical torque predictions to numerical precision" (30-torque.tex:139), which is the wording now used.
- A8/H: the old rib-midpoint vs lower-bolt-hole comparison numbers (χ₃ 0.62/0.60, bio-ext RMSE 2.195/2.204, filter 72/90 vs 83/90) are not in any chapter; the slide was reshaped to the recorded bracket points (App. C Table C-5) and adopted values (Table C-8). The pre-freeze comparison is acknowledged in the speaker notes only, with no numbers.
- A9/I: the old "coupled pathpoints stay equality-constrained" line is MuJoCo-conversion lore not present in the chapters; replaced by the human-target derivative-definition line (App. A §A.1).
- A10/J: bottom-left CAD card retained — now anchored at 60-sns.tex:93; the "insertion pending" phrasing follows the chapter ("tendon hardware is modeled in CAD but not yet inserted").
- No post-freeze/"since submission" numbers appear in this batch; every number traces to a chapter line.

## QC round 1 fixes (2026-10-18, this file only; no slide numbers or sources changed)

- Appendix K (deck slide 54, a11) HIGH: the bottom banner card covered the last deletions row
  ("interleg commissurals", −174.9). Fix: right table compressed to start y 2.35 in with 0.375 in
  spacing and 11.5 pt rows (slides_d.py a11 loop). Verified three ways: shape math (lowest row box
  bottom 5.73 in < card top 5.85 in), rendered-pixel check of deck_render/tmp_d/slide_11.png
  (banner top edge 702 px = 5.85 in; lowest right-table text pixel 662 px = 5.52 in; 40 px clear
  gap), and zero out-of-bounds shapes on rebuild.
- Appendix J (deck slide 53, a10) LOW: twin-figure axis text borderline. Figure boxes enlarged
  3.25 → 3.45 in tall (captions 10.5 → 11 pt; banner cards moved 5.85 → 5.90). PIL bbox audit of
  the source PNGs showed only ~2% croppable margin (41-42 px of 2250/1560), so autocrop was
  rejected as a no-gain change; the slide keeps the uncropped dissertation Fig. 6.11 asset.
- Appendix I (deck slide 52, a9) LOW: ~150 px dead band under the title. Equation card moved
  2.15 → 1.95 in and grown 1.5 → 1.55 in; conventions card moved 3.95 → 3.72 in and grown
  2.4 → 2.75 in (bottom 6.47 in).
- Cross-batch routing note (NOT this batch's slides): slide_index.json shows batch c's deck
  slides 35/36 citing the gait-validation sweep as "Table 6.2" and pathway removal as
  "Table 6.3"; the compiled dissertation numbers them Table 6.1 and Table 6.2
  (ProofFinal/build_r5/main.lot, entries at pp. 159 and 162). Batch d's own cites (Table 6.1/6.2
  on slide 54; Table 6.4 and Figs. 6.5/6.11 on slide 53) verify correct against the same .lot/.lof.
  slides_c.py belongs to the parallel batch and was not touched.
