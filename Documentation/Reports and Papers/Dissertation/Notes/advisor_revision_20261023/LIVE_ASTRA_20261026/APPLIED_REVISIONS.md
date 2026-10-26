# Applied Overleaf revisions — 26 October 2026

The live dissertation was revised through bounded selections in Reviewing mode. No chapter was replaced wholesale. No existing tracked change was accepted or rejected, and no advisor comment was resolved or reply posted. The peer-reviewed isometric-force chapter was not edited.

## Scientific argument

The revised argument identifies the need to distinguish neural-circuit behavior, sensory contribution, actuator mechanics, and numerical implementation in the human-scale and BPA models. The tools provide specified, variable model components and controlled tests. The representative results establish selected implementation agreement, corrective reflex behavior, and contact/Ib-feedback dependence in a supported W2L configuration. These conclusions do not claim unrestricted controller transfer, physiological validation of every model, or harness-free walking.

## Applied passages

| Patch | File and current local line | Purpose |
|---|---|---|
| A00 | `after/60-sns_LIVE.tex`:80 | Ben's requested heading; already applied live, verify rather than repeat. |
| A01 | `after/60-sns_LIVE.tex`:14 | Replace vague missing-link paragraph with the concrete experimental need. |
| A02 | `after/60-sns_LIVE.tex`:29 | Address advisor knowledge-gap request with concrete questions, contributions, and demonstrated answers; retain approved list order. |
| A03 | `after/60-sns_LIVE.tex`:43 | Credit existing tools and replace unsupported blanket statements about missing software capabilities. |
| A04 | `after/60-sns_LIVE.tex`:47 | Replace incomplete-work table framing. |
| A05 | `after/60-sns_LIVE.tex`:49 | Use approved tools-and-demonstrations table; Gait2392 second. |
| A06 | `after/60-sns_LIVE.tex`:70 | Explain the actual conversion problems and repairs, including Ben's reply about the initial robotbody attempt. |
| A07 | `after/60-sns_LIVE.tex`:78 | Define the BPA interface explicitly and connect it to the preceding conversion account. |
| A08 | `after/60-sns_LIVE.tex`:84 | Include the author-confirmed coupled SolidWorks assembly; distinguish representative component examples. |
| A09 | `after/60-sns_LIVE.tex`:108 | Explain what the reflex result establishes, rather than appending a generic outstanding-work warning. |
| A10 | `after/60-sns_LIVE.tex`:132 | Clarify bridge tests and add the reviewed fixed-setting feedback proof of concept. |
| A11 | `after/60-sns_LIVE.tex`:144 | Lead the discussion with completed contributions. |
| A12 | `after/60-sns_LIVE.tex`:146 | Connect each demonstrated result to the knowledge gap. |
| A13 | `after/60-sns_LIVE.tex`:152 | End on usable contributions, with specific scope rather than an incomplete-project narrative. |
| A14 | `after/60-sns_LIVE.tex`:16 | Restore correct Connor--Ben model and visible A/B labels. |
| A15 | `after/10-introduction_LIVE.tex`:12 | Repair grammar and duplicated text; state evidence and neuromechanical knowledge gaps without claiming no prior human model exists. |
| A16 | `after/10-introduction_LIVE.tex`:30 | Replace vague framework objective with concrete tools and demonstrations. |
| A17 | `after/70-overall_LIVE.tex`:16 | Make the whole-work discussion reflect the completed SNS contributions and distinct evidence. |
| A18 | `after/70-overall_LIVE.tex`:65 | Remove unsupported claim that arbitrary full-body model transfer is already cheap, automatic, and verified. |
| A19 | `after/80-conclusion_LIVE.tex`:10 | State the completed tools and what their demonstrations establish. |
| A20 | `after/80-conclusion_LIVE.tex`:12 | Keep the scope without presenting the dissertation as unfinished integration. |
| A21 | `after/80-conclusion_LIVE.tex`:14 | Conclude with delivered scientific and methodological contributions. |
| A22 | `after/15-background_LIVE.tex`:208 | Replace unsupported field-wide novelty claim with the scoped methodological gap; retain cited prior work. |
| A23 | `after/95-AppendixD_LIVE.tex`:56 | Substantiate the Rybak table entry with the verified resetting/non-resetting mechanism; identify Python and neural-only scope. |

Additional insertion: Appendix D section `sec:app_conversion_repairs`, detailing the human-reference conversion repairs and the separate robotbody conversion attempt. Uploaded `connor_ben_shared_rg_original.png` and `ablation_pair_metrics_20261026.pdf` to `figs/60-sns/`. Figure 6.1 now has visible A/B labels and the original shared-RG diagram.

## Advisor-comment and author-response matrix

| Concern / recorded response | Action | Status |
|---|---|---|
| “Describe the knowledge gap and how the tools I created fill them.” (Ben relaying advisor, this session) | Rewrote SNS gap/contributions/discussion and aligned introduction, background, overall discussion, and conclusion. | Applied; no field-wide absence claim. |
| “huh?” / “This paragraph confuses me.” | Replaced the vague five-system prose and “connections” descriptions with named models, blocks, integrations, and experiments. | Applied; live threads retained. |
| “need a better word” | Used neural pathways, actuator integration, signal exchange, and model implementation as appropriate. | Applied; no blanket substitution with “framework.” |
| “What BPA interface? Non seq” / Ben explained the initial robot import | Explained length/pressure-to-force calculation, tendon route, and sensory return; separated the failed robotbody conversion from human-reference repair. | Applied in chapter and Appendix D. |
| Advisor requested specific conversion problems / Ben described small muscles and robotbody routing | Added small-muscle deactivation, moving-path-point constraints, driver limits, patella tracking, and contact repairs with scope. | Applied. |
| “Huh? Just delete all this?” / Ben replied “deleted” | Preserved the tracked deletion of the suspended-stepping qualifier. | Not restored; comment remains open. |
| “for each leg?” / Ben described hip and knee/ankle PF groups and shared original RG | Correct original image and caption explain the bilateral PF groups and distinguish later per-side RG development. | Applied; reply preserved. |
| Ben: Researchers; Gait2392 second; tools and selected proof-of-concept tests | Opening starts with Researchers; Gait2392 second in list and table; no incomplete-work column. | Applied. |

## Evidence and number changes

- The 2.00/7.85-degree loaded-arm comparison is retained and interpreted as a functional model test.
- The new contact/Ib-ablation passage uses 33/28 versus 1/0 cycles and 1.16 versus 0.09 m from the saved paired metrics. The scored window is 25.499 s within each 30 s run; Ia/II remain active. It is the standard muscle-actuated W2L model, with identical harness configuration but differing measured support contribution.
- Table values, including 1.000 s rhythm-core agreement, correlation at least 0.99996, the 12-actuator BPA example and 372 kPa peak pressure, come from the previously reviewed packet. Their scopes are retained.
- The Rybak Python implementation is supported by `Code/MuJoCo_SNS/spinal/replication/verify_rybak2006.log`; the added explanation describes the distinct PF and RG perturbation responses behind its 13 checks.
- No new optimization or simulation was run. The held CPG demonstration figure and uncorrected clipped demo figures were not uploaded. The new Deng experiment and further W2L perturbation packet remain available for a separate bounded insertion; this pass did not import them wholesale.

## Verification and remaining limits

- Full local build: pdflatex, BibTeX, and two pdflatex passes exited zero; 266 pages; no errors, undefined references, or missing figures.
- Final Overleaf build: Errors 0, Warnings 0; downloaded PDF verified as 268 pages and contains the new gap text and Rybak subsection. `overleaf_compile_status.txt` preserves the UI status.
- Whole-file readback of all six revised source files matches the prepared content after normalizing line endings. The source/PDF are preserved here alongside the before snapshot and diff.
- The protected force chapter, existing bibliography, and main document are byte-identical between the local baseline and revised build. No commit or push was made.
- Figure 6.1, the new table, gap prose, feedback subsection, and conversion appendix pages were visually inspected. Box-overflow messages remain elsewhere in the complete document; this is not a claim that every page or advisor thread has been fully polished.
- The downloaded Overleaf PDF is the live compiled artifact. Local pagination differs between the build environments; do not replace the live project with the local tree.
- Advisor comments and replies remain in Overleaf. The matrix reports content addressed, not automatic thread resolution.
