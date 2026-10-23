# Defense revision handoff

23 October 2026. Target: approximately 35 minutes within Ben's stated 30–45-minute range, followed by audience and committee questions. Advisor favors assertion/evidence. No additional committee format requirements were supplied.

## Assessment and story

The 19 October source has 54 slides, including 11 appendix slides. Its useful plots are diluted by text panels, a broad simulation catalogue, and claims that precede the latest physical reductions. The revised sequence is: project gaps and attribution; characterized actuator force; transmission-model diagnosis and held-out tests; optimized routes and physical torque capacity; sensory evidence; selected neural-interface verification; conclusions and next experiments. A completed implementation appears only when it supports a stated gap or a verification boundary.

## Completed sections

- Rebuilt an editable presentation copy using existing figure components, Calibri, the original brown/teal palette, and 13.333 × 7.5-inch dimensions. Original deck and builder are preserved. `reference/Original_Defense.pptx` is an additional exact source copy.
- 31 main slides, including title and questions, plus 6 hidden backups. Sentence assertions lead the substantive slides. Technical explanations, attribution, scope qualifications, and delivery cues are in notes.
- Force: retained existing measured plots; scoped the comparison to reported cases; added an editable chart from the three reported 112 mm-case force values. Its embedded workbook is a literal snapshot of those reported values.
- Torque: separated characterization, correction identification (`gamultiobj`), route redesign (`surrogateopt` and `patternsearch`), and physical tests. Added final extensor measurements and the latest flexor reduction. Retained per-test holdout counts and exact validation means.
- Flexor: regenerated three plots from a local copy of the 23 October MATLAB generator. The numeric calculation, source arrays, reaction-point transform, correction values, and design-result loading are unchanged. Changes are plot size, font size, legend position, and output destination. No experimental data were digitized. Numerical rerun: 15/15 margins nonnegative; minimum 0.415038 N·m at −105°, maximum 17.004549 N·m at −14°. Displayed rounding follows the paper: +0.42 to +17.0 N·m.
- SADb: reduced to two main slides, emphasizing source-linked records, use, incomplete classification, and the distinction between literature evidence and controller parameters.
- SNS: distinguished initial one-RG suspended stepping from the later interconnected bilateral extension; added an editable conceptual bilateral diagram; retained selected reflex, matched rhythm-core, and bridge checks. Explicitly exposed the full-network tonic-drive discrepancy. Removed unsupported walking and Hunt/Young comparison framing.
- Broader pulley, supported-motion, and spiking results are hidden backups with scope qualifications. Equations and an earlier biomimetic-extensor configuration are backups too.

## Figure provenance and decisions

Exact paths, hashes, and crops: `build/asset_manifest.json`. Slide assertions, notes, and asset keys: `slide_plan.json`. A readable slide-by-slide source index is in `SLIDE_REVIEW.md`.

| Figure group | Source | Decision |
|---|---|---|
| Force jig, length curves, normalized curves | Existing `deck_builder/deck_assets/figs` components of the published Actuators study | Reuse and enlarge; comparison image selects the two upper 10 mm cases |
| Knee apparatus and earlier single-flexor tests | Existing deck figure components, checked against Revision_C Ch. 4 | Reuse; move interpretation to notes |
| Hybrid and corrected fit | Current Revision_C `FlxPin_preopt-eps-converted-to.pdf` and `FlxPin_optimized-eps-converted-to.pdf` | Replace older deck plots with current adopted-fold validation panels C/D; preserve data and add native shared ordinate label |
| Cam route and assembly | Revision_C `figs/30-torque/extensor_cam_route.pdf` and `cam_extensor_hardware.jpg` | Crop existing panel A, pair with original photograph |
| Predicted cam envelope | Revision_C `Vas_Pam_20mm_cam_winner_1_r2.png` | Add before measured tests; label as prediction |
| Final extensor reduction | Revision_C `extensor_20mm_physical.png`, manifest source `Code/Matlab/Mesh_Optimization/Results/Ext20mm_compare_20261022.png` | Add; nine holds, variable pressure |
| Latest flexor plots | `tmp/paper_bridge/figscripts/gen_Flx20mm_3panel_20261023.m`; local copied calculation `build/flexor_defense.m` | Regenerate presentation layouts only; final values in `assets/flexor_final_values.csv` |
| SADb citation and pivot views | Revision_C `figs/50-database` application captures | Reuse; crop pivot to selected columns and remove blank screenshot area |
| Initial AnimatLab body and traces | Revision_C `AnimatLab_reference/OurWalker.png` and `jointAngles.PNG`; Ch. 6 establishes chronology | Replace GUI capture of later two-CPG project and regenerated phase-1 traces with historical references |
| Later bilateral architecture | Editable conceptual diagram based on Revision_C Ch. 6 | Replace illegible full-connectome image in main talk |
| Library, reflex, rhythm traces, pipeline | Existing deck components and Revision_C `figs/60-sns` | Select representative checks; proposed hardware remains dashed in source pipeline |

## Source hierarchy and unresolved factual questions

1. **Flexor test pressure is material.** Recorded pressures are 226.5–638 kPa, with five holds above 620 kPa. All 15 exceed the selected human-muscle reference, but these are discrete holds at different pressures, not a measured continuous 620 kPa envelope. The pressure range is visible on the main slide. Small positive margins are point estimates; this pass did not establish repeatability or a confidence interval for the minimum margin.
2. **Latest flexor result supersedes the earlier Revision_C figure.** The 23 October MATLAB generator uses measured reaction-point geometry and confirms the paper's 15/15 claim. The older 22 October physical plot is retained as a reference asset but is not used for the final flexor claim.
3. **Extensor status conflicts across sources.** Revision_C has a completed nine-point MATLAB reduction, whereas a paragraph in the torque paper still describes six points and a pending reduction. The main talk uses Revision_C's final reduced figure. Confirm any subsequent rig/geometry changes before the defense.
4. **Nominal versus actual lengths:** the current figure caption explains that nominal 48, 41, 46, and 40 cm titles correspond to actual 48.5, 41.0, 45.7, and 40.6 cm lengths. The copy now uses the current validation panels (nominal 46 and 40 cm tendon), rather than the older deck's training/validation pair. Numeric validation means come from the matching adopted-fold table. No experimental labels were silently changed.
5. **Earlier 10 mm biomimetic extensor:** Revision_C retains an earlier 0.8920 N·m fit description and a later adopted shared-pair comparison of 1.17 N·m. The main talk omits that numerical comparison. Backup notes explain the different scopes instead of silently selecting a number.
6. **Target scope:** the final flexor target is maximum-isometric biceps femoris short-head torque from OpenSim. It is neither total human knee torque nor an experimental human-subject result. Keep the selected extensor reference distinct from it.
7. **SNS validation boundary:** matched rhythm-core results do not resolve the full-network discrepancy. Signal exchange, mechanism import, reflex behavior, controller equivalence, and locomotion are separate claims. No harness-free walking or completed embedded deployment is asserted. No fresh simulation benchmark was run for this presentation pass.
8. **SADb snapshot:** latest Chapter 5 describes 1,393 records (943 papers + 450 linked non-paper records). Old deck curation/coverage statistics were removed from the main talk. Captures are static snapshots and do not prove live behavior or curated accuracy.
9. **Historical AnimatLab capture:** original left-joint traces do not state angular units explicitly. The notes preserve that uncertainty and make no numeric amplitude claim. The historical outcome and initial one-RG organization are supported by the current chapter text. The older GUI capture names a later two-CPG project and was removed from the initial-model slide.

## Pacing and remaining work

Suggested rehearsal budget: opening/context 3 minutes; force 6; torque and physical validation 15; SADb 2; SNS 7; conclusion/transition 2. Total approximately 35 minutes. To reach 30, shorten apparatus descriptions and use one flexor prediction/mismatch explanation. To use 45, add discussion of the identification fold, static/dynamic boundary, and committee questions, rather than the full simulation catalogue.

Remaining author decisions: rehearse timing aloud; confirm the earlier extensor-fit scope and historical trace units if needed; confirm collaborator roles and preferred spoken credits; confirm whether later measurements supersede the current extensor figure. Repeat-pressure tests, uncertainty quantification, neural-controller reconciliation, pulley hardware tests, and walking are research work beyond this deck edit.

## Verification

The deck was exported by the presentation-authoring tool, structurally validated, opened by PowerPoint, and rendered to PNG/PDF for visual inspection. Visual repairs included external MATLAB legends, larger plot layouts, a corrected route crop, and corrected schematic arrowheads. The final report records the final inspection status below. All original presentation/dissertation files were treated as references. No dissertation edits, commits, pushes, or Overleaf actions were performed.

### Final inspection status

Deliverable: `output/Bolen_Defense_Assertion_Evidence_20261023_v6.pptx`. Companion PDF contains the 31-slide main talk; the six backups remain hidden in the PPTX and are available during questions. Notes exist on all 37 slides. All slides were rendered in PowerPoint and inspected individually, with changed views inspected again after repair. The package, slide geometry, Calibri font policy, first-party re-import, and embedded chart-data checks passed (`build/validation_v6.json`). The final PowerPoint text-bound check identified one title whose glyph bounds extend about 5.6 pt past its nominal text box while remaining inside the slide; visual inspection found no clipped text. No unresolved visual clipping or unintended overlap was observed. Some historical screenshots retain small interface text; the talk relies on their visual organization, not on reading that interface text aloud.

Original-deck SHA-256 comparison passed. The original deck, its builder/components, and dissertation sources remain preserved. Full source-slide dispositions are in `SLIDE_REVIEW.md`. Editable replacement content, embedded notes, provenance manifest, copied MATLAB generator, and numeric reduction CSV are retained for the next handoff.
