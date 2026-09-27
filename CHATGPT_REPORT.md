# ChatGPT report for ZCode

**Last updated:** 2026-09-26

**Current priority:** Ben's dissertation review

**Detailed historical record:** [`CHATGPT_REPORT_ARCHIVE.md`](CHATGPT_REPORT_ARCHIVE.md)

This is the current operational index. It is organized by durable item number and
status rather than by an ever-growing chronological transcript. Newest session updates
appear first. Older detail is retained in the archive and must not be deleted when an
item is summarized or superseded.

## Read this first: current state

- The authoritative source for the current review is Ben's downloaded Overleaf archive:
  `Documentation/Reports and Papers/Dissertation/Bolen_Dissertation (1).zip`.
- The ZIP itself is unchanged. Its isolated working copy is:
  `Documentation/Reports and Papers/Dissertation/Overleaf_review_20260926_yellow/`.
- The live Overleaf project was not modified during this review pass.
- The review PDF is:
  `output/pdf/Bolen_Dissertation_proposed_edits_yellow.pdf`.
- Yellow material in that PDF is proposed, not approved or merged into the live
  dissertation.
- The 160-page review PDF compiles with zero LaTeX errors, undefined citations,
  undefined references, or rerun warnings. The edited pages and inserted figures were
  visually inspected.
- Ben is reviewing the PDF now. Do not produce a clean merge until he identifies the
  yellow changes he approves.

## Status vocabulary

- **READY FOR REVIEW:** implemented in the isolated yellow review copy and awaiting Ben.
- **APPROVED:** Ben explicitly approved the item; safe to place in a clean integration copy.
- **PROPOSED:** drafted but not yet included in the yellow review PDF.
- **OPEN:** requires additional work, evidence, or a decision.
- **REVISIONS REQUESTED:** Ben reviewed the item and specified changes.
- **REJECTED:** Ben explicitly declined the proposal; retain only for provenance.
- **NOT TOUCHED:** explicitly outside the session scope; no current verification implied.
- **SUPERSEDED:** preserved for provenance but replaced by a newer item or decision.

## Active dissertation items

### Item D1.1: Current-Overleaf yellow review build — READY FOR REVIEW

Ben requested one PDF based on the current Overleaf source, with proposed changes shown
in yellow. The downloaded `(1).zip` was extracted to a separate review tree, edited, fully
compiled, and visually inspected. The live project and source ZIP remain unchanged.

**Primary artifact:**
`output/pdf/Bolen_Dissertation_proposed_edits_yellow.pdf`

**Decision needed:** Ben should approve, reject, or revise each yellow item before a clean
Overleaf-ready source copy is produced.

### Item D1.2: Abstract replacement — REJECTED

The current abstract was replaced in the review copy with a four-paragraph proposal that
summarizes the BPA force characterization, knee-torque correction and validation, human
torque comparison, Sensory Afferent Database, OpenSim-to-MuJoCo workflow, BPA interface,
and SNS-Toolbox controller tools. It does not invent numerical results for simulations
whose final configurations and metrics have not been frozen. Ben rejected this replacement
in review comment 9: the PSU abstract must be one paragraph, double-spaced, and no longer
than one page. Retain the original abstract as the base for Ben's revision. Review comment
10 also requests a missing comma in the original wording.

**File:** `chapters/02-abstract.tex` in the isolated review tree.

### Item D1.3: Dissertation organization section — READY FOR REVIEW

The one-sentence chapter summaries in `sec:organization` were expanded into full
paragraphs. The prose points to actual chapter and section labels for the isometric BPA
force study, knee-torque study, neuromechanical database, modeling tools, Results,
Discussion, appendices, and Future Work. Future Work is described as concrete uses of
completed research rather than unfinished dissertation requirements.

**File:** `chapters/10-introduction.tex` in the isolated review tree.

### Item D1.4: Steele knee mechanism figure — REVISIONS REQUESTED

The review copy uses only the two Steele mechanism panels selected by Ben. Panel
descriptions and the Steele citation are in the caption, not embedded beside replacement
panel headings. The figure is stored under the Background chapter-specific figure folder.

**Files:**

- `chapters/15-background.tex`
- `figs/Background/steeleknee.pdf`
- `figs/Background/steeleknee.png`
- `thesis.bib` (`steele_experimental_2023`)

### Item D1.5: Xi frame and two-bracket figures — REVISIONS REQUESTED

The two-rotation frame construction and tall two-bracket force-path projection figures
are included in Methods. The transformation and projection equations remain in the text.
The two-bracket figure uses frames `{br,1}` and `{br,2}`, equal-and-opposite force vectors,
`delta_tendon` along the path, `eta` for three-dimensional bracket displacement, and a
color-vision-deficiency-safe palette. Editable SVG and PowerPoint versions remain beside
the Methods figure assets.

**Files:**

- `chapters/20-methods.tex`
- `chapters/92-AppendixA.tex`
- `figs/Methods/xiFrameTransform.pdf`
- `figs/Methods/xiProjection.pdf`
- `figs/Methods/xiProjection.svg`
- `figs/Methods/xiProjection.pptx`

### Item D1.6: Completed-work neuromechanical framing — REVISIONS REQUESTED

Methods, Discussion, Future Work, and Conclusion now distinguish completed tools from
future experiments. The database, model-conversion workflow, BPA coupling, controller
interfaces, and two-layer spinal-network implementation are presented as reusable
research infrastructure. Future Work gives actionable experiments that unnamed
researchers can perform with those contributions.

**Files:**

- `chapters/20-methods.tex`
- `chapters/40-discussion.tex`
- `chapters/50-futurework.tex`
- `chapters/60-conclusion.tex`

### Item D1.7: Chapter-specific figure folders and paths — READY FOR REVIEW

Every chapter has a dedicated folder under `figs/`, and every chapter source begins with
a local `\graphicspath` entry. Existing `Aim1`, `Aim2`, and `Preliminary` paths remain as
fallbacks so the current document compiles without relocating all legacy figures in this
review pass.

**Folders:** `Introduction`, `Background`, `Methods`, `Results`, `Discussion`,
`FutureWork`, `Conclusion`, `AppendixA`, `AppendixB`, and `AppendixC`.

### Item D1.8: Weight-bearing simulation results — OPEN

Full-body closed-loop and weight-bearing simulations are running on multiple machines,
but the yellow review copy deliberately does not invent or freeze their outcomes. A
yellow Results note requests the final configuration, run duration, gait cycles or time
to failure, joint ranges of motion, duty factor, trunk and pelvis orientation, contact or
ground-reaction behavior, and cross-machine reproduction. Once Ben selects the runs of
record, those measurements should strengthen the Abstract, Results, and Discussion.

**File:** `chapters/30-results.tex` in the isolated review tree.

## Review round 1: PDF comments 1–44

These comments refer to
`output/pdf/Bolen_Dissertation_proposed_edits_yellow.pdf`. Do not overwrite or repaginate
that PDF while Ben continues reviewing it. Implement the accepted corrections in a new
review round after Ben finishes commenting.

### Item D2.1: Force-characterization figures and tables — REVISIONS REQUESTED

**PDF comments:** 1, 31–35, and 42.

- Resize the affected figures so their captions remain on the same page.
- Preserve the accessible ranked-result palette: indigo is the primary result, followed
  by progressively lighter colors for less important series.
- Increase the clearance between y-axis titles and tick labels.
- Figure 4.3 should separate the maximum-force result from maximum contraction. Add the
  10 mm maximum-contraction-versus-resting-length result and Moe's 20 mm data. The 20 mm
  result may be a new dissertation result, with an explicit single-production-batch
  limitation and a batch-number placeholder for Ben.
- Tables 4.1 and 4.2 need 10 pt text and should share one PSU-compliant sideways page.
- Table 5.1 is shifted into the left margin; use 10 pt text and move it to a sideways page
  if it cannot fit correctly at that size.

**Dependency:** ZCode should regenerate the figures from the current data and apply the
project-wide figure standards. Do not merely scale raster images.

### Item D2.2: Current torque-identification and biomimetic-knee results — OPEN

**PDF comments:** 2–7, 24, 36, and 37.

The pinned-knee figure, coefficient table, extensor paragraph, extensor figure, biomimetic
knee paragraphs, goodness-of-fit table, and complete biomimetic torque figure all contain
stale results. The Methods statement that the fit used the 48.5 cm flexor BPA has also
changed. ZCode must rerun or retrieve the current results of record, regenerate every
affected plot, update table values, and rewrite the associated Results prose together so
the numbers cannot drift between text, tables, and captions.

**Do not patch isolated numbers from the annotated PDF.** Confirm the current scripts,
input files, held-out tests, selected solution, and generated artifacts first.

### Item D2.3: Appendix provenance and accuracy audit — OPEN

**PDF comment:** 11.

Audit every appendix for current numerical results, filenames, script callouts, and model
provenance. In particular, confirm that the recorded workflow calls
`minimizeFlxPin10mm`, not `minimizeFlxPin10_2bracket`; the experiment used two brackets,
but the filename must match the actual script of record. Cross-check every appendix table
against the current Results tables after Item D2.2 is resolved.

### Item D2.4: Abstract constraints and punctuation — IMPLEMENTED LOCALLY; OVERLEAF SYNC PENDING

**PDF comments:** 9 and 10.

Reject the proposed multi-paragraph abstract. Restore Ben's abstract as the base, retain a
single paragraph, keep it to one double-spaced PSU page, and insert the requested comma.
Ben will revise the substantive abstract language.

**Implementation, 2026-09-26:** restored Ben's single-paragraph abstract from the current
Overleaf download, including the parenthetical comma pair around “from the trunk down.”
Removed the unapproved four-paragraph rewrite. Compile and Overleaf synchronization are
pending the end of this implementation pass.

### Item D2.5: Biomimetic-project rationale — IMPLEMENTED IN OVERLEAF

**PDF comment:** 12.

Revise the project rationale to express the iterative program explicitly: use knowledge
from biomechanics and neuroscience to build robots that perform complex tasks in
unstructured environments with improved efficiency, balance, and related capabilities;
use biomimetic robots to test, validate, and improve understanding of biological systems;
iterate between those directions; and translate the resulting knowledge to prosthetic and
orthotic devices.

**Implementation, 2026-09-26:** rewrote the opening project paragraph in
`chapters/10-introduction.tex` to state the iterative biomechanics/neuroscience-to-robot,
robot-to-biological-hypothesis, and translation-to-assistive-devices cycle.

### Item D2.6: Background figure program — PARTIALLY IMPLEMENTED IN OVERLEAF

**PDF comments:** 13–17, 19, and 20.

Add or replace Background figures as follows:

- Section 2.1: Kengoro and Kenshiro, or a more current biomimetic biped from Asano et al.
- Section 2.2: representative torque motor, BPA, hydraulic, and dielectric-elastomer
  actuator images.
- Section 2.3: BPAs from Ben's prior reports, Hunt's work, or Festo.
- Section 2.4: Ben's isometric-force-study graphical abstract.
- Section 2.5: a two-panel comparison of Steele's 3D-printed robot skeleton and Ben's
  `gait2392_robotbody.osim`, showing visible force actuators, the 27 right-side BPAs, and
  the left-side OpenSim Gait2392 musculature.
- Section 2.7: a multi-panel synthesis of important neural pathways from Rybak 2006, Deng
  2019, Shevtsova/Rybak 2016–2026, Shinohara, and Ijspeert.
- Section 2.9: use an enhanced still from Ben's supplied zombie-cat video rather than a
  sanitized schematic; document the source and citation/permission basis.

**Dependency:** source, citation, permission, image-quality, and alt-text review are
required before insertion. Prefer original or author-provided assets over web thumbnails.

**Asset check, 2026-09-26:** Ben's graphical abstract is available as
`Figures/Aim1/00_GraphicalAbstract.png`. No local Kengoro, Kenshiro, Asano, or Steele
skeleton source image was found by filename search. Existing `gait2392` screenshots are
only in deprecated folders, so the OpenSim comparison needs a new render from the current
model rather than reuse of a legacy screenshot.

**Implementation checkpoint, 2026-09-26:** added Ben's existing BPA graphical abstract
to the force-model background and added the existing literature-synthesis spinal-circuit
figure after the locomotor-control background, with dissertation captions and alt-text
comments. Robot, actuator-technology, OpenSim/Steele comparison, and zombie-cat panels
remain open pending source/permission or new-render work.

### Item D2.7: Caption construction, test-jig figure, and Xi figure order — OVERLEAF SYNC IN PROGRESS

**PDF comments:** 8, 18, and 21–23.

- The Xi and Steele captions are LaTeX captions, not text embedded in the figures. The
  current yellow `\colorbox`/`\parbox` wrapper forced “Figure 3.x:” onto a separate line.
  Remove that wrapper or highlight the complete caption in a way that preserves the normal
  inline label and first sentence.
- Replace Figure 3.2 with the correct three-panel `TestJigs` figure from
  `Documentation/Reports and Papers/Knee_Torque_Test`; check for duplicate filenames
  before selecting the source.
- Keep the two-rotation frame explanation, which Ben marked “ok for this I think.” Place
  the current Figure 3.6 immediately after the new two-panel frame figure.
- Reconsider the two-bracket projection equations and figure placement because they may
  precede the point at which the two-bracket model is properly introduced.

**Asset check, 2026-09-26:** the correct three-panel source is
`Knee_Torque_Test/Figures/Figure components/testJigs2/testJigs.pdf`; the similarly named
`testJigs1.pdf` files are single-panel force-test-stand figures and are not substitutes.

**Implementation, 2026-09-26:** replaced the dissertation asset with the confirmed
three-panel `testJigs2/testJigs.pdf` source; placed the bracket-frame orientation figure
immediately after the two-rotation construction figure; and moved the two-bracket series
projection equations and figure to follow the explicit introduction of the two-bracket
identification. The clean `ProofFinal` source has ordinary LaTeX captions, so the yellow
review wrapper that split the caption label is not carried forward.

### Item D2.8: Chapter and Results structure — OPEN

**PDF comments:** 25, 30, 38, 39, and 41.

- Rewrite the clunky AnimatLab section introduction.
- Move and rewrite the pasted Sensory Afferent Database, inverted-pendulum, and dynamic-BPA
  paragraph so it opens the neuromechanical-study Results rather than reading as an
  Introduction/Methods transplant.
- Move “Follow-Up Identification and Route Redesign” before the neuromechanical modeling
  section.
- Consider placing the Morrow/Bolen optimization-study results before BPA force
  characterization.
- Reserve the end of the route-redesign Results for Ben's final biomimetic torque data and
  plots. Follow it with the inverted-pendulum test results. Discussion should interpret
  those results alongside the updated `MonoPam_pulley` work.

**Implementation checkpoint, 2026-09-26:** moved “Follow-Up Identification and Route
Redesign” ahead of the balance and neuromechanical-recording sections. The current clean
Introduction already places the Morrow/Bolen optimization foundation before force
characterization, and the AnimatLab introduction has been rewritten. Final torque plots,
the precise neuromechanical-results opening, and `MonoPam_pulley` synthesis remain tied
to the current runs-of-record.

### Item D2.9: Neuromechanical figures and completed results — OPEN

**PDF comments:** 26–29 and 44.

- Add schematics for Ben's multiple ongoing neuromechanical models.
- Replace the poor knee-reflex demonstration image with a visually improved version and
  include the beer-cup arm reflex model where it supports the narrative.
- Replace the obsolete Simulink library and reflex figures with ZCode's newer schematics
  and library implementation.
- The MuJoCo–SNS Toolbox, SNS Simscape, MuJoCo–Simulink bridge, and related workflows are
  completed in substantial form. Present verified implementations and results in Results
  or Discussion instead of describing them only as future benchmark opportunities.

**Dependency:** ZCode must identify the current figures and runs of record and distinguish
validated results from active tuning.

**Asset check, 2026-09-26:** current source material exists as
`Code/Matlab/SNS_Simscape/figures/KneeReflex_circuit.*`,
`Code/Matlab/SNS_Simscape/results/pictures/sns_beer_results.png`, the corresponding
beer-cup Simulink export, and the dissertation's `CPG_airstepping_figs/circuit_literature.*`
and `sns_diagram_panels.*`. These are inputs for a publication-quality composite, not an
automatic final selection; in particular, the raw Simulink wiring export needs redesign.

### Item D2.10: Use “model,” not “plant” — IMPLEMENTED LOCALLY; OVERLEAF SYNC PENDING

**PDF comments:** 40 and 43.

Replace “plant” with “model” throughout the proposed dissertation prose unless the text is
quoting a source that specifically requires control-theory terminology. Apply the same
preference to captions, tables, and figure alt text.

**Implementation, 2026-09-26:** replaced every standalone authored occurrence of “plant”
or “plants” in the included chapter files with “model,” “biomechanical model,” or
“musculoskeletal model,” including captions and Future Work. A word-boundary search of
`chapters/*.tex` now returns no remaining occurrences.

## Items not touched in this session

- **Live Overleaf:** the iterative project rationale, BPA graphical overview, literature-synthesis
  locomotor-circuit figure, and reviewed three-panel knee-test figure have been applied. The Xi
  figure-order edits are in progress. Uploaded project-root assets are
  `circuit_literature.pdf`, `inverted_pendulum_photo.jpg`, `steeleknee.pdf`,
  `testJigs_reviewed.pdf`, `xiFrameTransform.pdf`, and `xiProjection.pdf`.
- **Canonical repository copies:** `ProofFinal/`, `upload/`, and `ZCode_drafts/` were not
  replaced with the yellow review sources.
- **Simulation and optimization code:** no MATLAB, MuJoCo, SNS-Toolbox, OpenSim,
  AnimatLab, or Simulink code or parameters changed.
- **Scientific data:** no experimental data, optimization outputs, or running simulation
  results changed.
- **Git publication:** no commit, push, rebase, history rewrite, or branch operation was
  performed for this review copy.
- **Clean dissertation integration:** no yellow proposal has been treated as approved.

## Verification record

### Item V1.1: LaTeX build — PASSED

- MiKTeX `pdflatex` and `bibtex` were run directly because this installation's `latexmk`
  wrapper requires an unavailable Perl interpreter.
- Final page count: 160.
- Final log: zero LaTeX or package errors, undefined citations, undefined references,
  rerun warnings, or fatal errors.
- Existing legacy overfull and underfull box messages remain outside the narrow yellow
  review scope; the build is not represented as a complete ETD typography audit.

### Item V1.2: Visual review — PASSED

Rendered pages were inspected for the Abstract, dissertation organization, Steele knee
figure, Xi frame figure, Xi projection figure, simulation-status block, Results insertion
note, Discussion synthesis, Future Work, and Conclusion. A two-line Conclusion orphan was
removed by tightening the proposed prose without changing its substance.

### Item V1.3: Source-scope checks — PASSED

- Every chapter has a chapter-local `\graphicspath` declaration.
- Malformed placeholder references such as `\ref(ch}` and `\ref{app}` are absent from the
  proposed material.
- Bracket displacement uses `eta`, avoiding conflict with Ben's use of epsilon elsewhere.
- The source ZIP and delivered review PDF were kept separate.

## Files created or changed by the current review session

- `AGENTS.md`: added Ben's standing Oxford-comma preference.
- `CHATGPT_REPORT.md`: reorganized current operational index.
- `CHATGPT_REPORT_ARCHIVE.md`: byte-for-byte preservation of the former 1,244-line report.
- `Documentation/Reports and Papers/Dissertation/Overleaf_review_20260926_yellow/`:
  isolated edited source and build tree.
- `output/pdf/Bolen_Dissertation_proposed_edits_yellow.pdf`: delivered review PDF.

## Session updates, newest first

### 2026-09-27 05:33 PDT: Annotated text corrections complete in Overleaf

The remaining terminology corrections are now saved in Conclusion and Appendix A, completing the
live `plant`-to-`model`/`system` pass. Results now also has the corrected 620-kPa caption spacing,
the `chi_0, chi_1, chi_2` sequence, the `(B)` panel label, and “achievable” spelling. All requested
source edits from this finishing pass are live. The next step is a manual recompile, log review,
targeted PDF-page inspection, and prohibited-name verification.

### 2026-09-27 05:31 PDT: Terminology pass in progress

The live Discussion wording now uses “biomechanical model.” Future Work now uses
“biomechanical system,” “mechanical model,” “MuJoCo model,” and “feedforward model” in the
four annotated locations. These edits are saved and Overleaf is recompiling automatically.
Conclusion and Appendix A are next, followed by the four Results corrections and final QA.

### 2026-09-27 05:30 PDT: Final live-edit and verification pass started

The live Results restructuring and Xi-method insertions remain saved in Overleaf. The current pass
is applying the remaining narrow terminology corrections in Discussion, Future Work, Conclusion,
and Appendix A, followed by the four small Results caption/notation corrections. A clean Overleaf
recompile and targeted PDF inspection will follow immediately. Plot/data items that depend on the
other coding assistant remain tracked as external dependencies rather than being represented as
completed edits. No assistant or tool names are being added to the dissertation source or comments.

### 2026-09-27 05:27 PDT: Live Overleaf implementation checkpoint

Overleaf now contains the revised iterative biology--robot--biology rationale, the existing BPA
graphical abstract placed in Background, the literature-synthesis locomotor-circuit figure, and
the reviewed three-panel test-jig figure. The new Xi frame-construction text, transformation
matrix, and frame-orientation figures were inserted next to their defining text. The old duplicate
bracket-frame float has now been removed, and the two-bracket deflection/projection equations and
figure were inserted immediately after the two-bracket model description. Results now opens with a
separate Completed Foundation Studies section, including the completed placement optimization,
Sensory Afferent Database, balance-platform, and dynamic-BPA studies. A direct source-file upload
was refused by the browser, so the Results synchronization continued through focused source edits.
The Results tail is now reordered to place the completed identification and route redesign before
the balance-platform and preliminary-simulation records; the route status was updated from work in
progress to the adopted shared-route result, and the balance demonstrations now have their own
section. Remaining live work is the other `plant` replacements, a clean recompile, and live-PDF
inspection. No assistant/tool names were added to Overleaf source or comments.

### 2026-09-26: Review round 1 comments captured

Ben supplied 44 annotated PDF comments. They are indexed under Items D2.1–D2.10 with
their original comment numbers, dependencies, and requested outcomes. The delivered PDF
was not overwritten or repaginated while review continued. The project-wide figure rules
now record Ben's preferred indigo-to-lighter importance hierarchy and required y-axis
title clearance.

### 2026-09-26: Report reorganized

The former chronological report was archived intact. This working report now uses durable
item numbers, explicit statuses, a current-state summary, verification evidence, an
explicit not-touched section, and newest-first updates.

### 2026-09-26: Yellow dissertation review delivered

The current Overleaf ZIP was used as the base for an isolated 160-page review PDF. All
visible proposals are yellow. The PDF compiled cleanly, and the changed pages were visually
inspected. Ben began reviewing it.

### 2026-09-26: Current Overleaf snapshot adopted

Ben downloaded `Bolen_Dissertation (1).zip`. It superseded stale local dissertation copies
as the base for this review. The original ZIP remains unchanged.

## Historical record

The detailed September 8 through September 26 chronological history, including earlier
Overleaf work, toolchain repair, neuromechanical figures, MuJoCo/SNS verification, Xi work,
figure-tree checkpoints, and prior dissertation consistency passes, is preserved verbatim
in [`CHATGPT_REPORT_ARCHIVE.md`](CHATGPT_REPORT_ARCHIVE.md).

When an archived item becomes active again, summarize its current state here under a new
or existing durable item number and link to the archived heading. Do not copy the entire
historical narrative back into this operational report.
### 2026-09-27 05:36 PDT: Live compile-path repair in progress

- The clean Overleaf rebuild identified missing-graphic errors rather than prose or equation failures.
- Confirmed that the four newly uploaded reviewed figures were placed under `chapters/`; their source paths are being corrected accordingly.
- Confirmed the existing Results artwork remains under `figs/Aim1/`, `figs/Aim2/`, and `figs/Preliminary/`, while the discussion comparison figure remains under `figs/Aim1/`. Restoring those search paths now, then recompiling twice to settle cross-references.

### 2026-09-27 05:40 PDT: Overleaf figure paths repaired

- Restored the Results chapter search paths to `figs/Aim1/`, `figs/Aim2/`, and `figs/Preliminary/` and the Discussion comparison path to `figs/Aim1/`.
- Added the `chapters/` search path to Background and Methods so the newly uploaded reviewed figures resolve.
- Reinstated the Steele biomimetic-knee figure with its source attribution and `fig:steele_knee` label; the old commented placeholder remains inert and will be cleaned after the build check.
- Next: explicit two-pass Overleaf compile, log review, and final project-wide terminology/prohibited-name audit.

### 2026-09-27 05:47 PDT: Recovered two truncated live sections

- Detected that the live Introduction had been reduced to only the revised project paragraph. Restored the complete reviewed Introduction, including Motivation, the labeled project section, Problem Statement, Research Objectives, and Dissertation Organization; the requested iterative biology-to-robot-to-biology rationale remains intact.
- Detected that the Methods balance-platform section was absent. Restored the full hardware/control-method description, the `sec:balance_testbeds` label, and the uploaded inverted-pendulum photograph using an explicit Overleaf path.
- All uploaded reviewed figures now resolve in the current build. Remaining work is formatting cleanup, a final compile/log check, and project-wide name/terminology audit.

# 05:53 PDT — Recovery complete; final Overleaf audit underway

- Restored the complete reviewed Introduction and the missing balance-platform / inverted-pendulum methods section after detecting that both live Overleaf files had been truncated during earlier edits.
- The dissertation is back to its expected 167-page length.
- Reduced the reviewed test-jig figure to `0.73\textwidth` to resolve its oversize-float warning, and repaired two insertion boundaries with explicit paragraph breaks.
- Final work now in progress: one remaining Background insertion-boundary cleanup, two clean compilation passes, project-wide prohibited-name and terminology searches, and spot checks of the rendered PDF.

