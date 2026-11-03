# v21 mathematical review

Reviewed 2 November 2026 from a mathematics committee perspective. This is an independent content review, not an impersonation of a committee member. No manuscript, result, plot, or deck-authoring file was changed.

## Verdict

No blocking numerical or equation-transcription error was found in the requested force and optimization-backup content. Two precision refinements are recommended before delivery. Full-slide PowerPoint visual review is still pending; this report does not certify the assembled deck.

## Corrections recommended

1. `backup_content.mjs:6`: replace the note that the maximum-shortfall term “prevents a good average from hiding one failed angle” with “The large weight prioritizes the worst-angle shortfall; acceptance is checked against the pointwise torque requirement.” A finite weighted sum is not a lexicographic optimizer or a feasibility theorem. The 100:1 weighting alone cannot prove that no trade-off is possible, particularly with additional penalties. The displayed title, “Route redesign prioritizes the worst torque shortfall,” is already appropriately scoped. The submitted source has similar strong wording; do not edit that protected source in this deck pass.
2. `build_identification_assets.py:29` (`shortfall.png`): the numerical floor in `tau = max(1, max |M^H|)` implicitly means 1 N·m. State this unit either in the equation or in its immediate explanation. All other dimensional terms in the inspected equations are consistent.

Neither refinement requires refitting, changing recorded results, or editing the dissertation.

## Checked and passed

| Item | Check and source evidence |
|---|---|
| Maximum-force arctangent | `F620 = b1 atan[b2(l0-b3)]` and both coefficient sets match frozen submitted `20-muscles.tex:95–100`: 303.5 N, 19.03 m^-1, 0.0075 m; 922.4 N, 15.37 m^-1, 0.013 m. |
| Pressure-dependent maximum force | `Pmax` appears in both prefactor and arctangent argument. The 10 mm coefficients 0.4895 N/kPa and 0.03068/(kPa·m) match the submitted source. The plotted pressure surface is correctly limited to the 10 mm dataset. The 7.5 mm offset is preserved rather than silently replaced by the separate 20 mm/620 kPa offset. Higher-pressure extrapolation is not presented as validated. |
| Normalized force | `max(0, c0[exp(-c1 epsilon*)-1] + P* exp[-c2(epsilon*)^2])` is equivalent to the published positive-force piecewise form. All six coefficients match `20-muscles.tex`, Table `tab:Fstar`. |
| Normalization | `F/F620`, `P/(620 kPa)`, and `epsilon/epsilon620` retain the original reference conditions. The proposed Pmax notation is not incorrectly reused to normalize pressure or force a second time. |
| Model-comparison tables | Independently parsed the frozen submitted table and compared all 12 rows / 36 numerical values with `model_comparison_table.json`: zero mismatches. All new-model values are lower than their corresponding comparison values, supporting the two table headings. These are reported comparison cases, not newly held-out-actuator experiments. |
| Identification objective | The mean of per-test, zero-correction-normalized RMSE/FVU/maximum-error components matches frozen `30-torque.tex:314–338`. Notes correctly recognize that FVU is proportional to squared RMSE within a fixed test and does not independently isolate waveform shape. |
| Identification bounds | The physical bounds and mappings agree with frozen `94-AppendixC.tex:240–258`. Flexor chi0: 0–20 mm; chi1: 30,000–1,000,000 N/m; chi2: 5,000–20,000 N/m. Extensor chi0: -20–5 mm; chi3: 0–1; stiffnesses locked. |
| Adopted parameters | The displayed values are consistent with the adopted front members and submitted rounded table: flexor/extensor chi0 +3.937/-6.378 mm, shared chi1 39,978.7 N/m and chi2 14,733.8 N/m, chi3 0.15829. Effective system stiffness is not presented as a literal material modulus. |
| Fold bookkeeping | Flexor row 77: held-out file indices 2 and 4; training 1 and 5; test 3 absent from the authoritative noT3 identification. Extensor row 32: held-out 1,2,8; training 5,6,7; excluded 3,4,9. Notes distinguish nominal file labels from measured lengths. |
| Candidate plots | Plotted columns match the original front compilation layout: flexor training 7–9, validation 10–12; extensor training 9–11, validation 12–14. Plots use only candidates from the adopted fold. The source hashes in `pareto_provenance.json` still match both authoritative MAT files. The plots are labeled as 2D views of retained 3-objective candidates, not as newly computed 2D Pareto fronts. |
| Validation scope | Notes explicitly state that validation entered candidate filtering/ranking, so no blind-validation claim is made. This distinction is important for defending generalization. |
| Saved solver settings | The displayed flexor population 150 / maximum generations 600 agree with the saved scalar inventory. Extensor saved `opts`: population 50, maximum generations 150, crossover 0.8, scattered crossover, adaptive-feasible mutation, tolerance 0.001, Pareto fraction 0.35, maximum stall generations 100. Missing flexor options are not inferred from today's defaults. |
| Route objective and penalties | Coefficients and selected penalties match frozen `30-torque.tex:376–446`: worst-angle 1e5, mean-square 1e3, rest 1e6, strain 1e7/1e4, extensor contraction 1e6, flexor off-axis/wrap 1e3, geometry 1e-2, length 1e-3. The table is a selected summary, not the complete objective. |
| Route search budgets | Flexor 1,000/15,000; direct extensor 3,000/15,000; cam 1,000/1,500 match frozen `30-torque.tex:453`. Notes correctly distinguish maximum evaluation budgets from actual evaluations, and this search from gamultiobj identification. |
| Five-percent allowance | Correctly described as a design requirement, not an experimental error interval. Required torque direction and the selected monoarticular-muscle targets remain explicit. |
| Chi4 | Correctly identified as a future compliance hypothesis, with no fitted or validated value asserted. This agrees with the submitted Discussion. |

## Asset inspection

Visually inspected the standalone general arctangent, Pmax arctangent, normalized-force equation, identification-objective equation, and extensor candidate-projection PNGs. Their symbols and mathematical grouping are legible and correct at source resolution. Actual slide scaling, equation legibility at projection size, and table wrapping require the parent’s assembled-slide render.

## Scope and source notes

- Frozen dissertation source: `C3NS_Revision_20261031_v18/submitted_source/chapters/20-muscles.tex`, `30-torque.tex`, `94-AppendixC.tex`, and the supplied Discussion context.
- Force source extracts: `assets/force/SOURCE_NOTES.md`, `force_equations.json`, and `model_comparison_table.json`.
- Identification evidence: `build/identification_inventory.json`, `saved_solver_options.json`, `pareto_provenance.json`, and the read-only front-assembly source used to verify column semantics.
- The historical BEN_TODO and campaign style/ruling files named in the dissertation skill are absent at their old paths. Current `DISS_STATUS.md` and the newest digest were read; no stale campaign was restored.
- No optimizer, controller, experiment, manuscript build, or publication upload was run for this review.

## Assembled PowerPoint render review

Individually inspected final slides 27, 28, and 130–138 at 1600 × 900 pixels. Image hashes and exact reviewed paths are recorded in `committee_jiang_v21_visual.json`. These eleven slides all pass visual checks: no clipped values, cut-off equations, overlapping labels, or unreadable table wrapping were found. The two earlier corrections are implemented in slide 130's equation and its notes.

| Final slide | Individual verdict |
|---|---|
| 27 | Pass. The 10 mm comparison table is spacious, values and units are legible, the two new-model rows are clearly emphasized, and the source is present. |
| 28 | Pass. The 20 mm table preserves all values, including 1321.1 N, without clipping. Column alignment and diameter/initial-length notation are clear. |
| 130 | Pass. All three equations are legible. The torque-scale floor visibly reads 1 N·m. The 5% allowance and required torque direction are explicit. Corrected notes distinguish finite penalties from final feasibility checks. |
| 131 | Visual pass. Table exponents, paired regularizer values, and the cam explanation are legible. Minor scientific wording refinement remains: prefer “Route search penalizes violations of actuator and packaging limits” to “Feasibility penalties enforce the actuator and packaging limits.” Finite penalties alone do not enforce hard feasibility. |
| 132 | Pass. Surrogate/global search and local pattern refinement are clearly separated, and numbers are labeled maximum function evaluations rather than observed run lengths. |
| 133 | Pass. Baseline normalization and the three objective meanings are readable. The table does not misleadingly separate shape from level error. |
| 134 | Pass. Bounds, units, signed offsets, and shared stiffness values are clear and unclipped. The fixed extensor stiffness pair is explicit. |
| 135 | Pass. gamultiobj identification is clearly distinguished from route redesign. Saved settings and missing-option status are legible; no absent flexor values are invented. |
| 136 | Pass. Both axes, candidate/adopted markers, fitting/validation labels, and the same-fold/three-objective qualification are readable. Notes retain the non-blind-validation qualification. |
| 137 | Pass. The extensor projection and adopted star are clear. The visible caption states which parameters are fixed and refitted. No new-front or unconditional-generalization claim is made. |
| 138 | Pass. Current fitted results are separated from the chi4 hypothesis. The absence of a fitted or validated chi4 is explicit. |

There is no visual blocker in this reviewed set. The remaining slide 131 wording refinement can be made without changing data or layout. This coverage is limited to the eleven named slides; the parent editor and other reviewers cover the remaining deck.
