# Paper state map — Knee_Torque_Test/ProofFinal (as of 2026-10-17)

Read-only audit for the paper_bridge session. Everything below was verified in this
session by reading the file or running the named check; nothing is inferred from the
brief alone except where marked "(from BRIEF.md)".

Checks run (commands, verbatim purpose):
- `ls -la` on `ProofFinal/`, `ProofFinal/Figures/`, `ProofFinal/PDF examples/`,
  `Knee_Torque_Test/` (parent), `Knee_Torque_Test/Bibliography/`, `Knee_Torque_Test/Tables/`.
- Full read of `Main.tex` (553 lines), `frontiers_SupplementaryMaterial.tex` (347 lines),
  `latexmkrc` (85 lines); head read of `FrontiersinHarvard.cls` via grep (class package
  loads L32–40, title-page/logo machinery L484–504, sectioning L335–354).
- Cross-check script `tmp/paper_bridge/mining/check_paper_state.py` run with
  `D:/Anaconda/python.exe`: cite-keys-vs-bib, includegraphics-vs-disk, ref-vs-label,
  TODO-marker scan, section outline. (Its printed line numbers for the outline are in
  comment-stripped text; all line numbers cited below are REAL file lines from the full
  read of Main.tex.)
- Bib-mirror comparison (Python, both copies normalized) and byte/CRLF census.
- `where pdflatex / latexmk / bibtex` and `ls D:\MiKTeX\miktex\bin\x64\*.exe` —
  all NOT FOUND (see §9).
- Word counts per section via Python (active, non-commented lines only).

All paths below are relative to
`Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/` unless prefixed.

---

## 1. Inventory

ProofFinal/ contains: `Main.tex` (54,821 B), `frontiers_SupplementaryMaterial.tex`
(18,835 B), `BolenFrontiers22.bib` (106,903 B, LF), `Frontiers-Harvard.bst`,
`FrontiersinHarvard.cls` (v2.6, 2016/11/10, cls L2–3), `frontiers_suppmat.cls`,
`stfloats.sty`, `latexmkrc`, `test.bib` (1,721 B, unused Frontiers sample — the live
bibliography is `BolenFrontiers22` per Main.tex:515), logo files `logo1.eps`
(required by the class title page, cls L486), `logo2.eps`, `logos.eps`, `YM-logo.eps`,
and `PDF examples/` (frontiers.pdf + frontiers_SupplementaryMaterial.pdf, template
samples). `Figures/` holds 14 files (list in §4).

Parent `Knee_Torque_Test/` additionally holds (untouched per folder rules): Ben's
Overleaf export `Torque_Values_About_Robotic_Knee_Using_Braided_Pneumatic_Actuators.pdf`
+ `.zip` (both 2026-10-17 21:39), the Connor-comments PDF (2025-09-29),
`Bending_bracket_math.jpg`, `Paper Gantt chart.xlsx`, `Bibliography/BolenFrontiers22.bib`
(107,993 B, CRLF), `Tables/CoefficientTable.xlsx` (2026-05-27), and local
`Figures/` + `Figure components/`.

**Bibliography mirror status:** the two `BolenFrontiers22.bib` copies have identical
key sets (62 = 62; Python check output "identical after newline normalization: True").
The only difference is line endings: `Bibliography/` = CRLF, `ProofFinal/` = LF.
Content is currently in sync; keep it that way when adding entries (from BRIEF.md §folder
rules the mirrored copy must stay in sync).

## 2. Section map — present vs missing

Main.tex active outline (real line numbers):

| Element | Line | State |
|---|---|---|
| Author/affiliation/correspondence block | 42–56 | present, 4 authors (see §3 flag) |
| Title (short + long) | 64 | STALE — see §3 |
| Abstract + keywords | 72–82 | present, early register (§3) |
| Introduction | 84–95 | present (611 words), gaps (§3) |
| Methods and Materials | 98 | present (1,500 words) |
| — Overview | 100 | present |
| — Actuator Force Model | 104 | present (ε, ε620, F*, Fmax, F — eqs 1–6) |
| — Isometric Knee Torque | 160 | present |
| —— Robot Model | 162 | present |
| —— Knee Torque Data Collection and Equipment | 197 | present |
| —— Hybrid Torque Calculation Method | 207 | present (thin; ends at L221) |
| Results | 228 | present (2,395 words) |
| — Flexor on Pin Joint | 231 | present (χ0–χ2 model derivation, eqs 7–22) |
| — Extensor on Pin Joint | 380 | present (χ3 wrap-loss, eqs 23–25) |
| — Extensor on Biomimetic Knee | 420 | present |
| — Flexor on Biomimetic Knee | 443 | present (human-torque comparison) |
| Discussion | 456 | present but SHORT — 547 words, 4 paragraphs + a dangling scratch fragment (L469–474) |
| Resource Identification Initiative | 476 | present (MATLAB/OpenSim/Arduino RRIDs) |
| Conflict of Interest / Author Contributions / Funding / Acknowledgments / Supplemental Data / Data Availability | 486–507 | all present |
| Bibliography | 513–515 | Frontiers-Harvard + BolenFrontiers22 |
| Figures (trailing Frontiers stub) | 521 | empty stub (figures are embedded inline in the body) |
| Tables (trailing) | 536 | holds `tab:OptCoeff` (L538–548) |
| **Conclusion** | — | **MISSING ENTIRELY** |

Missing vs the editorial contract (BRIEF.md §Editorial contract, Actuators-paper
structure): **no Conclusion section**, and the Discussion does not carry the LONG-form
additional results the contract requires — there is **no two-bracket flexor/extensor
identification (Xi0–Xi3 system-identification) content, no route-redesign
(mesh-optimization attachment redesign meeting human isometric torque targets) content,
and no fictive-strain (ε620) framing of the correction terms** anywhere in Main.tex.
The only ε620 content is the definitional eq. (2) at L114–120 — a natural hook for the
fictive-strain rewrite. Methods likewise has no two-bracket method, no route-redesign
method, and no human-benchmark method (the human-values paragraph survives only as the
commented block at L225).

Also absent: any 20 mm extensor biomimetic-knee results (only flexor 20 mm at L445);
the results narrative stops at the flexor biomimetic knee.

## 3. Stale elements

1. **Title block (L64):** `\title[Torque Values About Robotic Knee]{Isometric Torque
   Values About Robotic Knee Using Braided Pneumatic Actuators}` — old title. Working
   title per BRIEF.md: "Design Factors for Isometric Robotic Knee Torque Using Braided
   Pneumatic Actuators" (short: "Design Factors for Isometric Robotic Knee Torque").
2. **Abstract register (L72–82, 231 words):** early-register draft. Enumeration style
   ("has resulted in 1) ... and 2) ..."), present-tense procedure narration ("test
   setups are created"), and a future-tense promise ending ("This will lead to faster
   design processes..."). Keywords (L80, 8 of them) include the odd pair
   "Function Fit". Nothing in the abstract hints at the design-factor/system-ID scope
   of the target framing.
3. **Introduction gaps (L84–95):** intro ends at the one-bracket χ0–χ2/χ3 story and the
   attachment-optimization pointer. No setup for (a) the two-bracket identification and
   cross-configuration consistency requirement, (b) the route-redesign results meeting
   human isometric targets, (c) the ε620 fictive-strain expression of correction
   factors, (d) the "design factors" framing of the new title.
4. **Dangling Discussion fragment (L469–474):** an unnumbered, unpunctuated scratch
   derivation ("In equation $m\ddot\varepsilon^* + c(\dot\varepsilon^*)\dot\varepsilon^*
   + k(\varepsilon^*)\varepsilon^* = F$ given $k(\varepsilon^*,P^*)= ...$ then if we set
   mass m and external force F to zero we get ...") — reads as leftover working math,
   not manuscript prose. Must be replaced by real Discussion content.
5. **Author-count flag (from BRIEF.md, verified):** Main.tex:44 lists 4 authors
   (Bolen, Burgess, Morrow, Hunt). The dissertation publication note
   (`Dissertation/ProofFinal/chapters/30-torque.tex:5`) names **7** contributors —
   Bolen, Burgess, Morrow, **Elzein, Pang, Scharzenberger**, and Hunt — for
   bolen_isometric_2025. Missing from Main.tex: Elzein, Pang, Scharzenberger.
   Ben's call; do NOT change (BRIEF.md instructs flag-only).
6. **Author Contributions defect (L490):** "Investigation, B.B., C.M., and C.M." —
   C.M. listed twice; L.B. absent from Investigation.
7. **In-text numbers disagreeing with the supplementary GoF table** (never harmonize
   silently — flag per contract). Active-text vs `tab:TorqueFit`
   (frontiers_SupplementaryMaterial.tex:163–266):
   - Main L445 (flexor biomimetic 10 mm, l0 = 41.5 cm): text RMSE 0.3007 / FVU 0.0102 /
     max 0.4647 vs table 0.680 / 0.052 / 1.055.
   - Main L445 (flexor biomimetic 20 mm @ 620 kPa): text 1.8024 / 0.0499 / 3.5678 vs
     table 2.868 / 0.126 / 3.550.
   - Main L426 (extensor biomimetic, "52.0 cm"): text 0.8920 / 0.3396 / 1.4188 vs
     table (51.8 cm row) 1.461 / 1.041 / 2.755.
   - Commented-out older numbers also disagree with the table: Main L235 (pre-opt
     4.8196/0.7082/9.8652 and 2.9424/0.1545/6.6968) and L376 (1.616/"0.0.0958"/3.2527
     and 2.1214/0.0878/5.1690) vs table rows 4.871/0.723/9.865, 3.024/0.163/6.753,
     1.608/0.095/3.242, 1.915/0.072/4.695.
   - Resting-length inconsistency for the extensor biomimetic test: Main L426 says
     "l0 = 52.0 cm", its own figure caption L432 and the supp table say 51.8 cm
     (file name `Ext10mm_52cm.eps` says 52).
8. **Commented-out salvage blocks in Main.tex** (50 full-line comments; active text
   unaffected): old RMSE blocks L235/L376, Arduino/Windows note L203, the
   human-benchmark paragraph L225 (cites delp_interactive_1990, hoy_musculoskeletal_1990,
   yamaguchi_planar_1989, sherman_what_2013, millard_flexing_2013, thelen_adjustment_2003,
   seth_opensim_2011 — all keys EXIST in the bib, currently inert), extensor intro L382.

## 4. Figures — \includegraphics vs disk

`Figures/` (14 files): Ext10mm_52cm.eps, ExtPin_group.eps, F-P_contraction.eps,
FlexText4_pinned_FishScale_leak.png, FlxGrp.eps, FlxPin_angleMisalign.png,
FlxPin_group.eps, Flx_brktBend.png, Gait2392_rendered.PNG, KneeExt_10mm_passthru.png,
KneeICR.eps, bktFrame.pdf, extensor_misalign.png, testJigs.pdf.

Main.tex active targets (script check, all EXISTS):

| Line | Target | Label |
|---|---|---|
| 169 | Figures/testJigs.pdf | fig:testJigs |
| 180 | Figures/KneeICR.eps | fig:ICR |
| 239 | Figures/FlxPin_group.eps | fig:FlxPin_group |
| 367 | Figures/bktFrame.pdf | fig:bktFrame |
| 388 | Figures/ExtPin_group.eps | fig:ExtPin_group |
| 430 | Figures/Ext10mm_52cm.eps | fig:Ext10mm_52cm |
| 437 | Figures/KneeExt_10mm_passthru.png | fig:Ext_passthru |
| 449 | Figures/FlxGrp.eps | fig:FlxGrp |

Supplementary active targets (both EXISTS): Figures/extensor_misalign.png
(supp L321, fig:extensor_misalign), Figures/KneeExt_10mm_passthru.png (supp L329,
fig:passThru — the SAME image as main fig:Ext_passthru, presented twice).

**No missing targets; no broken figure references.** Files on disk not referenced by
any active includegraphics (candidate material, currently orphaned):
F-P_contraction.eps, FlexText4_pinned_FishScale_leak.png, FlxPin_angleMisalign.png,
Flx_brktBend.png, Gait2392_rendered.PNG. Supplementary figures referenced only in
comments point at files that DO NOT exist: Figures/LoadCellCalibration.png,
Figures/RobotSim_frames.pdf (supp L313, L67 — commented).

Caption defects to fix when touching figures: panel-D mislabeled "(A)" in the
bktFrame caption (Main L369 — reads "(A)...(B)...(C)... and (A) flexor configuration
on the biomimetic knee with the 20 mm BPA"); L445 "it is \ref{fig:bktFrame}\textbf{D}"
missing "Fig." before \ref. None of the existing result figures states the QC-required
"negative values denote knee flexion" caption note (checked all captions L171, L182,
L241, L369, L390, L432, L439, L451 — none mention it; BRIEF.md figure_qc requires it).

## 5. Citations — \cite keys vs bib

Script check result: **Main.tex has 37 active cite instances, 21 unique keys —
ALL 21 EXIST in BolenFrontiers22.bib. ZERO missing keys.** Supplementary has 0 active
cites (but keeps `\bibliography{BolenFrontiers22}` at supp L345 → empty References
section there). Bib has 62 entries, no duplicate keys.

Active keys (usage counts): asano_design_2017 (1), asano_musculoskeletal_2019 (2),
aschenbeck_design_2006 (1), bolen_2026 (4), bolen_determination_2019 (4),
delp_interactive_1990 (2), delp_opensim_2007 (2), goldsmith_investigating_2021 (1),
hunt_development_2017 (2), hunt_modeling_2017 (1), ijspeert_amphibious_2020 (1),
morrow_optimization_2020 (4), schilling_adaptive_2022 (1), seth_opensim_2018 (2),
shin_constructive_2018 (1), steele_biomimetic_2018 (2), steele_biomimetic_2018-1 (1),
steele_development_2017 (1), wu_for_1995 (1), wu_isb_2002 (1), young_analyzing_2019 (2).

Keys cited only inside Main.tex comments (inert, all present in bib):
hoy_musculoskeletal_1990, millard_flexing_2013, seth_opensim_2011, sherman_what_2013,
thelen_adjustment_2003, yamaguchi_planar_1989.

35 bib entries are currently cited nowhere (main + supp incl. comments) — headroom for
the Discussion expansion without adding new bibliography: aksoz_design_2019,
al-rawi_joint_1985, andrews_general_1985, arnold_model_2010, brand_model_1982,
chen_realizing_2019, chou_measurement_1996, derrick_isb_2020, festo_2026,
he_mechanism_2020, hitzmann_anthropomorphic_2018, hunt_modeling_2017-1,
john_contributions_2012, kurumaya_musculoskeletal_2016, lang_robert_musclesim_2005,
leibach_development_2020, li_influence_2020, liang_comparative_2020, liu_robotic_2018,
madden_artificial_2004, mahapatra_multi-legged_2020, martens_modeling_2017,
menegaldo_moment_2004, mulholland_activities_2001, romero_four-bar_2020,
sarosi_comparative_2017, sarosi_investigation_2013, sarosi_new_2012,
scharzenberger_design_2019, suderman_moving_2012, suzumori_trends_2018,
tondu_static_2022, veale_towards_2016, wu_isb_2005, zajac_determining_1989.

Note: `bolen_2026` is an `@unpublished` entry (bib L1085) — the Actuators 15(5):230
paper per BRIEF.md; check its fields say that before submission.

## 6. Cross-document structure (xr)

- Main.tex:37 `\externaldocument{frontiers_SupplementaryMaterial}`; `latexmkrc`
  (whole file, John Collins 2023-03-29 xr recipe) auto-compiles the supplementary
  into `output-subdoc/` and copies its .aux next to its .tex. BRIEF.md's note applies:
  if xr complains about a missing aux, compile the supplementary once first.
- **`tab:TorqueFit` — cited twice in Main.tex active text (L417, L426: "goodness of
  fit values ... found in table \ref{tab:TorqueFit}") — is defined ONLY in the
  supplementary** (supp L165). Ref-vs-label check: no undefined references anywhere
  (script output "Main.tex ACTIVE refs with NO label anywhere: NONE"), i.e. it resolves
  through xr. But a load-bearing GoF table living in the supplement while both results
  subsections lean on it is a structural choice to revisit (Frontiers main-text budget:
  currently 8 figures + 1 table = 9 of the permitted 15 combined — comment at Main L523;
  there is room to move it into the main text).
- Supplementary content: "Supplementary Methods and Materials" (supp L38) is ENTIRELY
  commented out (moment-arm derivation, transformation matrices, bracket-frame
  construction — supp L40–150); active supp content = tab:TorqueFit, `table:muscle_attachment`
  (supp L274–300), and the two figures. Supp class `frontiers_suppmat` with the same
  preamble duplication (double `\DeclareUSUnit\inch{in}` at supp L17–18 as well).

## 7. Class/preamble constraints (what may and may not be added)

Class: `FrontiersinHarvard` v2.6 (cls L2–3), utf8 option (Main L14). The class
already loads (cls L32–40): inputenc, crop, graphicx, amsmath, array, color, amssymb,
flushend, stfloats, amsthm, chngpage, times, datetime, epstopdf, parskip [=v1],
xcolor [usenames,dvipsnames], caption [labelsep=quad,indention=10pt],
subcaption [labelfont=bf,list=true]. It defines the \helvetica* font commands (cls
L34–37) and prints `./logo1.eps` on the title page (cls L486) — logo1.eps is present.
Body font is Times (`times`, cls L32); page layout is class-controlled (crop + class
geometry) — do not add `geometry`.

Main.tex preamble (L16–37) adds: url, hyperref, lineno, microtype, xr, amssymb, amsmath,
bm, siunitx, graphicx, caption (plain re-load — safe, no option clash), tabularx,
booktabs, multirow, array, setspace [onehalfspacing]; plus
`\let\DeclareUSUnit\DeclareSIUnit` / `\let\US\SI` (L21–22, siunitx-v2 idiom),
a duplicated `\DeclareUSUnit\inch{in}` (L23–24 — harmless), `\providecommand`
\abs/\norm (L25–26), `\renewcommand{\arraystretch}{2}` (L34), `\linenumbers` (L35).

- **MAY add:** `\usepackage{mathabx}` — currently NOT loaded (verified by preamble read
  + grep "mathabx" → no hits). BRIEF.md requires the engineering `\diameter` symbol and
  says to add mathabx and flag it: **FLAGGED — mathabx is missing and must be added**
  before using `\diameter`. Today the manuscript uses `$\phi$` for diameter in 8 places
  (Main L165 [comment], 185, 369, 382 [comment], 422, 432, 445, 451) plus the supp
  table header `\textbf{\phi}` (supp L171) — all need conversion.
- **MAY add:** booktabs-style tables are already supported (booktabs L30, used in
  tab:OptCoeff L542–547).
- **MUST NOT / beware:** do not load `caption`, `graphicx`, `amsmath`, `amssymb`,
  `array`, `xcolor` WITH NEW OPTIONS (class already loaded them; option clash error);
  do not add `geometry`; keep natbib-style `\citet/\citep` (class + Frontiers-Harvard.bst
  pairing); keep `logo1.eps` in the compile folder; keep the Frontiers trailing
  "Figures"/"Tables" stub sections and ≤15 combined figures+tables (comment L523);
  tables must be real LaTeX, not images (comment L550).
- **siunitx version risk (compile-pending):** the `\DeclareUSUnit/\US` shims are v2
  idioms; if the MiKTeX install ships siunitx v3 they are deprecated-but-present
  aliases — untested here because TeX is not installed (§9). Verify on first compile.
- Figures embedded inline (review style) with `\begin{center}` wrappers; Frontiers
  wants EPS-class figures like the existing set when replacing files (from BRIEF.md
  §Figure specification).

## 8. TODO markers

Script scan (TODO/FIXME/XXX/todo/fillme/TBD/PLACEHOLDER/???, case-insensitive):
**0 hits in Main.tex, 0 hits in frontiers_SupplementaryMaterial.tex.** The unfished
work expresses itself as the §2/§3 gaps, not as markers.

## 9. Toolchain / compile status

`ls D:\MiKTeX\miktex\bin\x64\{pdflatex,bibtex,latexmk}.exe` → "No such file or
directory" for all three; `where pdflatex / latexmk / bibtex` → not on PATH either.
**No TeX on this machine yet — every compile is PENDING, per BRIEF.md ("treat missing
binaries as compile pending, never as a manuscript failure").** No build was attempted
in this session. The PDF in the parent folder (2026-10-17 21:39) is Ben's Overleaf
export, not a local build.

## 10. Typo / unit-convention ledger (line-cited, for the edit pass)

Style-contract violations (BRIEF.md §Editorial contract):
- Bare units outside siunitx: L174 `$6 \, \text{ft}$`; L199 `$3/16\; \text{inch}$`;
  L201 `$0-120\; \text{psi}$` and `$3/8 \; \text{inch}$`; L445 `$l_{0}=41.5 cm$` (×2),
  `$0.4647$ $N\cdot m$`, `$3.5678$ $N\cdot m$`; L462 `$\pm \; 20\,\text{mm}$`.
- `$\theta_\text{k}=0 \deg$` (L174) — `\deg` is a math operator, renders "0deg";
  should be `\ang{0}` for consistency with L171/L182.
- `width=1\textwidth` (L239, 367, 388, 430, 449) — works, but sloppy; use `1.0`.

Typos/grammar (line-cited): L88 "robot robot"; L90 "actuator force model force";
L128 "it's total"; L195 "given din \citep"; L373 "the the model"; L384 "offet";
L415 "gamultobj" (correct "gamultiobj" at L373); L445+L451 "acheivable"; L462
"the  hybrid" (double space), "BPA This is why" (garbled splice), "a change the force
vector", and the unit error `$l_{0}=\qty{42}{\mm}$` (the extensor BPA is 42 cm —
cf. L384's correct `\qty{42}{\cm}`); L369 caption panel "(A)" duplicated where "(D)"
belongs; L445 "for the $\phi$\qty{20}{\mm} it is \ref{...}" missing "Fig.".

## 11. Summary judgment

The manuscript is a complete, self-consistent OLD draft: every citation resolves, every
figure target exists, no TODO markers, no undefined refs. What it lacks is exactly what
the bridge session must add — the new title/abstract/intro register, a Conclusion,
a LONG Discussion carrying the Xi0–Xi3 two-bracket identification and route-redesign
additional results in ε620 fictive-strain form, the `\diameter`/mathabx and siunitx
unit cleanup, the caption QC pass ("negative values denote knee flexion"), and
reconciliation (flag, not silent) of the in-text vs supplementary GoF numbers.
