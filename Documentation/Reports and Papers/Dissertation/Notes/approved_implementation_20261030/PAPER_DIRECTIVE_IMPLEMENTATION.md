# Ben's paper corrections: source implementation

30 October 2026. Implemented only Ben's U01, U02, U06, and U08–U16 instructions in `Knee_Torque_Test/ProofFinal/Main.tex` and `frontiers_SupplementaryMaterial.tex`. U10 was already satisfied and its citations remain unchanged. No unmarked T/OP finding was implemented independently. No source bibliography, figure, model input, result file, metric, or optimization was edited/run by this reviewer.

The two exact pre-edit files are preserved under `evidence/paper_source_baseline`. Reviewable diffs and checks are `evidence/Main.tex.diff`, `evidence/frontiers_SupplementaryMaterial.tex.diff`, `PAPER_IMPLEMENTATION_CHANGES_20261030.json`, and `PAPER_IMPLEMENTATION_VERIFICATION_20261030.json` (the JSON files are also under `evidence`).

| User instruction | Source disposition / current lines |
|---|---|
| U01 — remove supplementary Table S1 | Removed its table and both references. Main:409–413 retains the pre-correction comparison, all common-domain sample counts, the hybrid/full-test distinction, and the authoritative main fit table. The removed table remains in the exact baseline. |
| U02 — remove unused42-cm slippage photograph | Removed the supplementary figure/callout and its inventory entry. Main:541 retains the independently relevant path-sensitivity discussion and apparatus reference; the photo-dependent slippage/±20-mm sentence is removed with the excluded-test illustration. The original is preserved in baseline. |
| U06 — scientific segue and 27-site contribution | Main:101 now connects corrected joint prediction, both physically tested knee routes, and the existing 27-site modeled comparison. The 37.28→65.19% values are copied from the unchanged result table and scoped to inherited routing/sizing/transmission assumptions. Main:108 introduces RoM for its later use. |
| U08/U09/U14 — replace and relocate pulley | Main:397 retains a Methods pointer. Main:621–647 places the single new biomimetic-pulley include in §4.4 after physical capacity evidence, with the force/travel equation in text, total input2F20, origin reactionF_R, and short “Example of…” caption. There is no mount-reaction calculation. Root supplies the new asset. |
| U10 — Morrow cost citation | Already satisfied. Existing Morrow-citing lines are unchanged. The present knee-route objectives remain explicitly this study's formulation. |
| U11 — hybrid removal from current Figure11 | Main:405/494–500 removes the biomimetic-extensor hybrid references, caption phrase, and premature figure callout; pinned-knee hybrid evidence remains. Basename `Ext10mm_52cm_current_publication.pdf` is retained for root's corrected graphic. |
| U12 — shorten §4.3.1 | Main:576–596 uses two short outcome paragraphs and a native two-row capacity table. The section decreases from899 to245 regex-counted tokens, including table/LaTeX tokens under the same method. Existing objectives and unique direct-route/cam geometry values move to supplement:177–192; none is recalculated. |
| U13 — §4.3.2 results first | Main:600 opens with all-nine/all-fifteen measured capacity outcomes, then the figure, ranges, and margins. Main:602–606 retains pressure-specific prediction comparisons, hung loads, pressure ranges, adjoint procedure, dimensions, and the discrete-static scope. |
| U15 — Figure13 wording/style | Caption at Main:615 says “cam-assisted extensor route.” Both original panel basenames remain for root's matched-font/line-width/scatter-order/opacity/edge export. This reviewer did not alter either graphic. |
| U16 — exact transmission selection account | Main:647 describes the saved maximum-meeting-cell winner and points to supplement:245–275. The supplement states the executed search objective, a term/weight table, length bounds, pattern-search setup, escalation rule, and integer lexicographic winner ranking, distinguishing the search's positive-strain secondary term from the subsequently reported positive-insertion-force coverage. |

## Existing-code provenance for U16

The current canonical replay loads the saved winner files identified by `tmp/paper_bridge/capability_table_20261022/BiPulley_survey2_perMuscle_20261023.csv` (`regrade_capability_active_20261028.m:6–15`). The result assembler reads that roster and the frozen replay (`assemble_capability_active_20261028.py:12–30`). The dated survey driver in that same folder, `Opt_run_BiPulley_survey2_20261023.m`, defines:

- Executed cost at432–477: `−(1000*meetFrac + activeFrac) + 1e−3*deadFrac`; active search condition is at458 and is not silently replaced by the later reported force-active metric.
- Bounds at428–429; pattern-search seeds/settings at480–517.
- Winner comparison at599–612: greatest integer meeting-cell count, then greatest search-active count, then smallest actuator-count-plus-travel-multiplier; a remaining tie preserves the earlier rung.
- The actual ladder's force/travel diagnosis is in `tmp/paper_bridge/opt_improve_20261022/escalation_advisor_20261022.m:120–159,200–228`.

These are recovered existing operations, not an invented objective or retrospective new assessment method. No code or saved result was executed to recompute a scientific quantity.

## Conservation and integration checks

The main pinned/biomimetic fit tables, adopted coefficients, final capability table, all metric equations, common knee-route objective, and every retained supplementary coordinate/size table compare exactly to baseline. The configuration-specific objective block moved to the supplement is exact. Common-domain sample counts and all physical-test values are retained. Unique direct-route release angles, endpoint/length values, and unadjusted cam comparisons are preserved in the supplement. The authorized Table S1 values and excluded-photo displacement statement are the only removed numerical evidence; both remain archived in baseline.

Static source checks find no duplicate labels, no missing reference target across the two source files, no dangling removed-table/photo reference, and exactly one new pulley include after the physical-capacity section. The root reviewer owns figure exports and the fresh integrated build/visual check; this note does not claim a compile or final figure verification has yet occurred. No upload or commit was performed.
