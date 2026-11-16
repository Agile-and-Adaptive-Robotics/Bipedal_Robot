# Chapter 6 and Appendix D revision

The canonical working dissertation is `Final/main.tex`, with the rewritten Chapter 6 and Appendix D in `Final/chapters/60-sns.tex` and `95-AppendixD.tex`. Active figures are in the sibling `Figures/60-sns/` and `Figures/95-AppendixD/` directories. At Ben's explicit direction, the generated review ZIP was removed; the manuscript is delivered as ordinary source files and figure folders. The other sixteen chapter/front-matter TeX files are synchronized with the latest downloaded source, with their previous versions preserved in Archive. Bibliographies, classes, and support files were not replaced.

The Chapter 7 import retains the verified distinction between supplementary hand drive and retained shared descending drive as one disclosed exception to the downloaded text. All other Chapter 7 content follows the latest supplied source. The synchronization copies individual chapter files while preserving the revised Chapter 6 and Appendix D.

## Synchronization of the other chapters

Source: `Bolen_Dissertation.zip`, downloaded November 15 at 21:58:20, 36,977,754 bytes. `FINAL_SOURCE_SYNC.json` records the input hash and every before/ZIP/final chapter hash. The sixteen imported files are preserved beforehand under `Archive/final_source_sync_20261116/chapters/`. The rewritten `60-sns.tex` and `95-AppendixD.tex` were excluded and their hashes verified unchanged.

Only three imported files differ textually from the previous local Final versions: the author's latest Introduction, the torque publication note (`lumbosacral actuator-placement` becomes `actuator placement`), and Chapter 7. The published muscle chapter and the other twelve files have line-ending differences only. The author's Introduction is imported verbatim. The ZIP's removal of the higher-level command-interface and embedded-scaffold paragraphs from Chapter 7 is retained; the previous versions remain in the backup. `OTHER_CHAPTERS_SYNC.diff` shows the textual changes without line-ending noise.

Final/main.tex is already textually identical to the downloaded main file. The synchronized chapters require no additional figure or bibliography imports. Static checks against the actual active `thesis.bib` and all Final sources report zero missing citations, references, duplicate labels, or graphic files.

## Canonical figure organization

Thirty-five superseded files were moved, preserving their hashes: 22 into `Figures/60-sns/Deprecated/20261116/` and 13 into `Figures/95-AppendixD/Deprecated/20261116/`. Current figure previews and alt text remain with their active assets. The missing active `LI_DENG_FALL_FRAME.png` was recovered byte-for-byte from the saved authoritative baseline. `FIGURE_ORGANIZATION.json` records every move and the recovery hash. All active graphic references in every Final chapter now resolve.

Historical ProofFinal sources were not edited. Several retain old figure paths; those historical builds need the preserved Deprecated assets. Final is the working source of record for this revision.

## Changes and directive checks

| Directive | Implemented result | Evidence |
|---|---|---|
| Coherent dissertation chapter using the Aim 3 and workflow/calibration slides | Introduction, Methods, Results, Discussion, and conclusion; functional design, Hunt calibration, and physical curriculum have distinct purposes | `Final/chapters/60-sns.tex`; `CHAPTER_REVIEW.md` gives section/line anchors |
| Independent judgement of structure, prose, grammar, captions, and evidence | Separate reviewer read both complete revised sources and current evidence; identified two caption clarifications, both applied and read back | `SUPERVISOR_REVIEW.md` |
| Separate Appendix D revision | Dedicated appendix author reorganized supporting checks, added relevant script map/current results, consolidated duplicate/irrelevant displays, and retained source evidence | `APPENDIX_REVIEW.md`, `APPENDIX_CHANGE_MANIFEST.json` |
| Current relevant results | BiW2Lvar+ air/noise/pulse evidence, newer W2L108 supported-ground test, CAD reference-tension trial, superseding PF timing audit, and explicitly scoped timing comparisons | Claim-to-source tables in the chapter/appendix reviews |
| Current SNS-Simscape graphics | Three new vector library/knee/cup figures, verified six-path topologies, actual saved traces, current shared muscle glyphs, panel-aware captions, and native 14 cm/10 pt sizing | `GRAPHICS_REVIEW.md`, `GRAPHICS_QC.json`, `figscripts/` |
| Remove irrelevant appendix sections, scripts, and figures | Removed active catalogues/duplicate plots, reduced a historical control table, and replaced old CAD chronology; original scientific source files and recordings remain preserved | `APPENDIX_REVIEW.md` lists all removals and retained anchors |
| Include relevant scripts and figures | Appendix implementation/source map; retained original shared-RG circuit; readable current stimulus and PF vector panels | Appendix D and its change manifest |
| Follow-up on both remote machines after Ben commits | Two machine-specific prompts with exact evidence identities, source locations, verification gates, and remaining graphics checks | `EB475WS4_FOLLOWUP.txt`, `EASTEREGG2_FOLLOWUP.txt` |

## Scientific scope and numerical changes

The chapter and appendix review ledgers disclose all replaced, added, and removed numerical statements. No source MAT/NPZ/SLX recording, model parameter, force equation, optimizer result, or simulation trajectory was changed, and no new scientific simulation was run.

The most consequential clarifications are:

- The current BiW2Lvar+ candidate passed two suspended air stages, not a complete five-stage physical curriculum. Stage conditions differ in commands and sensory return, so the period change is not attributed solely to feedback.
- The syn6 output fit did not verify its intended PF timing: zero of twelve channels passed the superseding audit. The s3k audit likewise does not establish an intended timing design. Failed bounded repairs remain failed.
- The newer W2L108 and CAD ground examples are separate configurations. Support windows and definitions remain explicit. CAD uses ideal reference tensions and an engineering allocator, not 54 dynamically calibrated BPAs. No harness-free walking is claimed.
- The knee transient uses a simplified actuator-force element. The cup uses the BPA20 block; its quoted sag values occur at 10 s with an actual 0.625 kg ramp endpoint, rather than a settled 0.5 kg condition.
- Complete-process, neural-call, warm-kernel, and full diagnostic timings retain their distinct denominators. The latest kernel gain is not substituted for whole-simulation acceleration.

## Graphics and implementation conservation

The icons represent the same twelve block types shown on defense slide 69. That slide contains native editable PowerPoint drawings; the new chapter figure derives from the actual local library mask definitions. The shared muscle emitter had three display-only coordinate errors, corrected in `Code/Matlab/SNS_Simscape/sns_icon_muscle.m`. MATLAB R2025b verified the emitted geometry against the shared vertices. Existing primary `.slx` files retain their stored expanded graphics and were not overwritten; the supplied refresh script makes a separate library copy for normal-session inspection.

MATLAB read and exported the exact saved numeric data and mask primitives, but graphics initialization and the Simulink localhost helper failed within this restricted execution environment. Final publication figures use a standard vector backend on those exported values. Both native MATLAB and fallback renderers are retained. This is not reported as a successful MATLAB figure export or a newly simulated result.

The main air-motion figure rearranges exact vector strips from the two-page portable recording; no curves are digitized, smoothed, or recomputed. Its positive-flexion display is declared in the caption while the model retains flexion-negative coordinates. Full original plots remain in the handoff. Appendix stimulus plots are split into readable landscape strips, preserving all seven rows.

## Verification and delivery

`FINAL_STATIC_CHECK.json` records the actual Final-tree scan: no missing source references, duplicate labels, or missing citation keys. `FIGURE_ORGANIZATION.json` records zero missing Final graphic files after reorganization. These static checks do not establish a clean LaTeX build. The earlier `STATIC_PACKAGE_RESULT.json` is retained only as historical evidence for the superseded package workflow.

Compilation and rendered-page verification are separate checks. The earlier staged build did not finish: it stopped during scratch font-cache preparation with Windows API error 5. Its 646-byte PDF is invalid. The subsequent bounded build of actual Final/main.tex also failed before reading the manuscript because MiKTeX attempted automatic format generation and its child encountered registry-access error 5. `FINAL_BUILD_RESULT.json` records that blocked status. All owned jobs were stopped and the temporary Final/pdflatex.fmt was hash-verified and removed. No new dissertation PDF or compiled-page verification is claimed. The normal-session no-ZIP driver is `build_final.py`.

Full Windows source/figure paths and both machine-specific prompts are provided directly in chat. The prompt text is also retained beside this report. No Overleaf upload, live comment disposition, commit, or push was performed. Ben's pre-existing working-tree changes were preserved. Scratch staging, rendering, and build probes remain in ignored `tmp/ch6_20261116/`.

## Files to commit from this revision

- `Final/chapters/60-sns.tex` and `95-AppendixD.tex`.
- New `*20261116*` assets, the recovered active fall image, and the recorded moves into Deprecated under `Figures/60-sns/` and `Figures/95-AppendixD/`.
- `Notes/ch6_revision_20261116/`, including the evidence/review ledgers, generation scripts, and follow-up prompts.
- The display-only `Code/Matlab/SNS_Simscape/sns_icon_muscle.m` correction.

After synchronization, compile Final/main.tex directly using its sibling Figures folders. Do not use the superseded packaging helpers to replace the canonical Final sources. The runtime note documents the laptop's isolated TeX setup; other machines should use their normal machine-specific compiler configuration.
