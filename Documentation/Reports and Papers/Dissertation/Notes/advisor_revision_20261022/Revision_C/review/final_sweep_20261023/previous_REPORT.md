# Revision C

An edited, independently compilable local dissertation copy. Started from the committed round-two chapter proposals at 39e3f707 and unchanged supporting ProofFinal sources. Original proposal chapters, ProofFinal sources, bibliography bytes, and protected chapter-3 figures passed the integrity check. This pass made no Overleaf writes or Git commits/pushes. A merge to 39de1ccc was observed during the work; the dissertation, companion torque-paper, workbook, and flexor-endpoint sources used here did not differ across those commits.

## Main editorial decisions

- SNS now follows the AnimatLab construction and Ben's interconnected left/right RG extension → unresolved coordination and verification questions → integration requirements → tools → representative verification checks. The main chapter has five figures and one capability table. Larger reference sweeps, pathway removals, spiking comparisons, and supported-walker evaluations are in Appendix D. GUI inventories and debugging chronology were removed from the compiled argument; original text and figure assets remain recoverable.
- SADb leads with the resource and distinguishes the paper-only application from the five-table working base. Coverage details, linked-pathway tables, and the citation neighborhood are in Appendix E. Actual application screenshots replace retired plots; the pivot is framed through LaTeX to remove its empty viewport and enlarge the matrix excerpt.
- Torque retains final knee-design results in the main chapter. Earlier multi-site batch diagnostics, whole-body route-builder history, and transmission surfaces are in Appendix C. The newer companion-paper EPS exports were used for the cam route and two transmission figures: the similarly named old PDFs were stale. The false “two pulley positions” caption is corrected to actual joint-angle, pressure, and transmission variables. The two surfaces occupy separate pages. The second flexor path is separately solved, as confirmed in the endpoint and prediction code; the old shared-route statement is corrected.
- Introduction, selected background explanations, Discussion, and Conclusion were aligned with those decisions. The erroneous cross-product moment-arm definition was removed. Captions were shortened without adding unmeasured claims.

## Latest instruction: completed-study wording

Ben reported that the extensor data arrived last night and flexor collection is in progress, and directed writing as if the intended result were established. The abstract outcome sentence is therefore retained, and the torque summary and Conclusion use completed-result wording. This is an explicit drafting assumption. No missing flexor values, error scores, or confidence intervals were invented. Before release, the final extensor/flexor comparison must substantiate the stated human-reference outcome. The spreadsheet-derived cam table and torque range were removed after Ben clarified that Excel rows 16 and 17 are sanity checks; final torque values and presentation plots come from the test-specific analysis scripts. See `review/07_Bilateral_and_torque_source.md`. The active bipulley optimization on easteregg2 is a new campaign, not a completed result represented by the older sweep figures.

## Length

Approximate narrative source-word counts, excluding captions and tables, using the campaign helper with escaped-percent and end-of-file fixes. LaTeX labels/arguments can contribute tokens, so these are comparative counts rather than publication word counts.

| Chapter | Round two | Revision C | Reduction |
|---|---:|---:|---:|
| 10-introduction.tex | 2,227 | 1,447 | 35% |
| 15-background.tex | 5,387 | 5,000 | 7% |
| 30-torque.tex | 9,958 | 8,960 | 10% |
| 50-database.tex | 2,085 | 743 | 64% |
| 60-sns.tex | 6,865 | 1,244 | 82% |
| 70-overall.tex | 3,382 | 3,307 | 2% |
| 80-conclusion.tex | 519 | 441 | 15% |

SNS plus new Appendix D has 2,410 narrative tokens versus 6,865 in the old SNS chapter. SADb plus new Appendix E has 1,019 versus 2,085. Across changed files, including added appendices, narrative/caption/table tokens decrease from approximately 39,243 to 31,804. Movement to appendices was accompanied by selection and shortening.

## Per-directive check

| Directive | Source evidence |
|---|---|
| Separate local copy; no Overleaf writes | [main.tex:194](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/main.tex:194>) |
| SNS begins with AnimatLab attempt | [60-sns.tex:6](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/chapters/60-sns.tex:6>) |
| SNS has tools-and-verification structure | [60-sns.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/chapters/60-sns.tex:25>) |
| Detailed SNS results moved into appendix | [95-AppendixD.tex:60](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/chapters/95-AppendixD.tex:60>) |
| SADb shortened; coverage/pathways in appendix | [96-AppendixE.tex:6](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/chapters/96-AppendixE.tex:6>) |
| Pulley axes/configuration caption corrected | [94-AppendixC.tex:289](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/chapters/94-AppendixC.tex:289>) |
| Updated cam geometry shown | [30-torque.tex:542](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/chapters/30-torque.tex:542>) |
| New interconnected bilateral RG construction credited | [60-sns.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/chapters/60-sns.tex:12>) |
| Final torque results sourced from analysis scripts | [30-torque.tex:578](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/chapters/30-torque.tex:578>) |
| Completed-study wording per latest instruction | [30-torque.tex:578](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/chapters/30-torque.tex:578>) |
| Conclusion uses that drafting assumption | [80-conclusion.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/chapters/80-conclusion.tex:8>) |
| Unsupported-walking claims scoped | [60-sns.tex:115](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/chapters/60-sns.tex:115>) |

The protected isometric-force chapter is byte-identical. All 19 entries in the protection manifest (active/retired figures and chapter source) match their recorded hashes; the copied active figures also match. The bibliography is byte-identical. See `review/source_integrity.json`.

## Verification

- Full MiKTeX build: pdflatex → BibTeX → three pdflatex passes, all successful in a fresh temporary directory. No undefined references/citations, duplicate labels, changed-label rerun warning, or oversized-float warning in the final log.
- Compiled output: `Revision_C.pdf`, 254 physical PDF pages. SADb occupies printed pages 112–116; SNS 117–127. Supporting appendices remain part of the complete document.
- Rendered and inspected the changed SADb/SNS sections, new appendices, and mechanical figures. Page renders and overview sheets are in `review/rendered`; the manifest records the exact pages.
- Existing wide equations/tables in inherited Appendices A–C still produce overfull-box warnings. These pre-existing layout issues were not declared repaired; the new moved tables were fitted to the text width. This is an editorial revision with a successful build, not a claim that all inherited typesetting is publication-ready.
- Final torque results must use the test-specific analysis-script outputs. Excel rows 16 and 17 are sanity checks, not presentation data. The spreadsheet-derived table was removed. `review/torque_sheet_evidence.json` is retained only as a superseded inspection record and is explicitly marked accordingly; no reduction script was run during this editorial correction.

## Review and integration

Open `Revision_C.pdf` to read, or `main.tex` to compile this separate copy. `chapters` contains the edited TeX; `figs` contains existing assets and the three newly converted existing torque EPS exports. Preserve that folder layout.

`review/changes.diff` compares the revised chapters against the saved round-two baseline and includes the new appendix registrations. `review/figure_dispositions.json` lists retained, moved, and omitted figure labels. `review/numeric_change_review.json` contains every removed/replaced source line containing numerals, including relocated lines: it is an inspection ledger, not a claim that each number disappeared. Detailed section notes explain scientific corrections, explicit omissions, and moves. Do not apply the entire copy over an advisor-edited Overleaf project blindly; integrate the selected chapter/section changes and required assets against the current live source.

The new data campaign, the expanded SADb app export, unsupported walking, and hardware/network deployment were not completed or evaluated by this editorial pass.
