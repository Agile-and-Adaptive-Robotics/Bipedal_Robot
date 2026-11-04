# v21 neural sources and main-talk sequence

Prepared on easteregg2, 2026-11-02. This is a deck source/implementation plan, not new scientific work. Original documents, controllers, figures, and results remain unchanged.

## Confirmed assets

Staged 37 byte-identical new-campaign assets, 7.36 MB, plus 13 historical tool figures under assets/neural. ASSET_MANIFEST.json records source paths, timestamps, sizes, and SHA-256 for these 50 copies. A separate 13-file EB475WS4 delivery is recorded in eb475ws4_motor_loops/SOURCE_MANIFEST.json.

- ground9_rc_s4/PRESENTATION_REPLAY.gif and PRESENTATION_STILL.png: three actual recorded cycles at 1× playback, 7.663–13.000 s, selected ground9 RC-enabled S4.
- ground9_rc_s4/PRESENTATION_NEURAL.svg/.pdf/.png: actual bilateral RG, four-cell RF/HAM PF, and hip/knee/ankle MN traces, with legends outside the axes. Preserve colors and amplitudes. This portrait plate needs selected rows or a faithful landscape rearrangement to be readable in the main talk.
- ground9_rc_s4/ALIGNED_MECHANICS.pdf/.png: matched joint, foot-loading, and apparatus-force traces. The displayed knee ordinate is explicitly positive flexion; do not silently relabel as raw OpenSim knee coordinates.
- rfham_selected/RFHAM_FOURCELL_SELECTED.svg/.pdf/.png: all 12 incident edges of the actual four right-leg PF cells. RF integrates hip flexor and knee extensor PF input; HAM integrates hip extensor and knee flexor PF input. Each integration cell excites its relay, which inhibits the opposing integration cell. Original SVG retains vector text/shapes; full-precision fields, exact IDs, source ledger, and edge audit are adjacent.
- task_feedback/HIGH_VS_LOW_TASK_FEEDBACK.svg/.pdf/.png: saved excessive-versus-reduced task-feedback comparison, useful for a technical question.
- contact_feedback/STAGE5_CONTACT_STRENGTH.svg/.pdf/.png: saved full-versus-quarter contact-PF comparison, useful for the main experimental inference or a question.
- benchmark/BACKEND_AND_VALIDATOR_COMPARISON.md: paired runtime numbers and exact trajectory equivalence.
- sources/RESULTS_AND_FILES.md and COMMITTEE_QA.md: exact controller/configuration scope and committee explanations.

The user supplied the authoritative new EB475WS4 folder: Code/MuJoCo_SNS/spinal/functional_subnetwork_review/drawings/full_motor_loops. README.md and VERIFICATION.json identify the historical human w2lvar standing export, whose actual local ledger SHA-256 matches 49b946c627222e030aa7bc741d7de20faeb804646d791b51548e6d3d704d4e2f. Its zero-effective PF-to-MN routes are gray/dashed. This is not the selected ground9 controller. Two readable vector crops of its RF motor/force and Ia/II sensory maps are included as question slides, with scope visible in their source captions and full off-viewport context in notes. The full original SVG/PDFs remain alongside the copies. Main 65 retains the actual selected RF/HAM subnetwork, whose source identity matches the new result.

## Scientific scope to carry into notes and visible labels

The new result is the **human Gait2392 model with graded biw2lvar-Hunt-PT, four-cell RF/HAM PF modules, and a separate autogenic Renshaw repair**. It is neither the 108-population W2L robot nor historical Syn6. It has 1,033 actual cells and 92 muscle outputs; four extra logging names are aliases.

Ground9 RC S4 records five right/four left qualifying unloaded, whole-foot-mesh-clear swings during the 5–13 s assessment. Mean net superior assistance is 48.848356% BW, comprising 40.051043% prescribed support and 8.797312% signed vertical restraint. Horizontal restraint RMS is 115.560892 N (15.03% BW). These are guided steps in place. Do not label this unassisted walking or a minimum required support threshold.

The same RC-enabled saved controller passes nominal S1–4. S5 at contact-PF scale 0.25 passes; full-strength S5 and the next lower-support test fail clearance. S4 and quarter-strength S5 repeats have identical saved traces and compiled models. Eight additional S4 ±1% trials perturb ten selected scalars with fixed RC settings and pass. This is finite selected-parameter evidence, not a stability theorem or an all-controller result.

The paired benchmark is 199.624 s native versus 29.284 s accelerated, **6.82× complete-process speed-up**. Neural calls alone are 179.079 versus 9.780 s, **18.31×**. Fully instrumented ground trials remain slower than simulated real time. Preserve the distinction from the old suspended-W2L benchmark.

Hunt A–D are functional calibration steps, distinct from physical validation stages 1–5. This campaign reused and checked intrinsic/timing fits, then refined four shared local feedback targets and three paired motor settings. Cite Hunt's meeting presentation, not a journal paper: A. J. Hunt, “Real-time control of dynamic walking in a dog-like robot,” INRC2 Meeting, 15 February 2019.

The RF/HAM engineering integration module is not an exact reproduction of the public rat circuit. Sartorius and gracilis remain separate. Local FSA targets are engineering calibration choices, not measured human synaptic values. These facts belong in notes unless they change the visible comparison.

## Proposed replacement of v20 logical slides 57–68

Keep v20 slide 56's historical suspended/unsupported movies and its failure-to-transfer transition. The next slide establishes the knowledge gap; scientific results follow the tool sequence. Numbers below identify v20 replacement positions, not final v21 pagination.

| v20 position | Proposed title or purpose | Authentic evidence and transition |
| --- | --- | --- |
| 57 | Detailed neuromechanical models require a faster, repeatable testing workflow | One clean model/controller/sensory-loop schematic using historical controller and actual Gait2392 sources. The missing automated workflow is explicit without implying that AnimatLab lacked a graphical circuit editor. No results plot here. |
| 58 | The model library links simplified walkers with detailed human mechanics | Actual W2L and Gait2392 models, with their controller/body boundary clearly identified. Covers the first two rows of Table 6.1 without claiming identical controllers. |
| 59 | Neural blocks connect sensory input to actuator commands | Actual SNS-Simscape library and one readable connected reflex/actuator subsystem. Put the loaded-arm response and detailed agreement plot in questions. |
| 60 | CAD assemblies bring robot mechanics into the workflow | Actual SOLIDWORKS/Simscape assembly and BPA block. Distinguish the existing CAD component example from the laptop's new 3-DoF-foot work. |
| 61 | The bridge exchanges commands and sensory measurements | Actual MuJoCo–Simulink interface diagram or screenshot, with two clear signal directions. No unsupported new bridge benchmark. |
| 62 | Staged experiments separate timing, motor output, and feedback | Hunt's actual four-part method mapped to the experiment drivers. Keep air/standing/ground environments visibly separate from functional fitting steps. |
| 63 | Literature circuits support targeted perturbation experiments | Actual reduced Deng circuit and a clearly positioned stimulus; one short result may demonstrate what the tool enables. Move multi-model and spiking taxonomy to questions. |
| 64 | Compiled execution reduced complete-process time by 6.82× | Native editable timing table with native/accelerated seconds and separate complete-process and neural-only rows. Benchmark configuration is explicit in notes; it is not the ground9 trial. |
| 65 | Four-cell RF and hamstring modules combine hip and knee input | Actual selected RF/HAM vector schematic. The new EB475WS4 drawing describes a different saved controller and is therefore retained separately in questions. |
| 66 | The selected controller maintains alternating neural activity | Readable RG/PF/MN selections from PRESENTATION_NEURAL, preserving shared scales, original legends, and recorded-time labels. Do not shrink the entire portrait plate onto one slide. |
| 67 | The human model takes repeated steps with measured support | Large actual GIF plus legible matched foot-load evidence. Visible: guided steps in place; 5 right/4 left swings; 48.8% BW vertical assistance; 115.6 N horizontal guidance. Exact assessment interval in notes. |
| 68 | Quarter-strength contact feedback restored clearance in this test | Actual S5 full-versus-quarter comparison. This turns a limitation into a specific sensory-control experiment and leads naturally to future work. If main-slide budget is tight, put this in questions and end the main results on 67. |

This maintains Table 6.1 order through models, library, CAD, bridge, experiment drivers, literature circuits, and acceleration. The new selected-controller mechanism and results appear afterward. It replaces twelve positions with twelve slides rather than appending a second results section.

## Questions and existing results

Retain the authentic original SADb offline movie and SADb questions at the end. Detailed original W2L108 results, matched NumPy/Simulink tests, loaded-arm reflex verification, complete neural/mechanics plates, all-controller distinctions, finite repeat/parameter checks, RF/HAM synapse fields, FSA equations, and support decomposition belong in questions. Keep historical submitted results and the post-submission campaign separately identified in notes and source captions.

New question evidence should answer: what makes this biomimetic; which parameters were fitted; what the support apparatus contributes; whether the same controller passed the environments; what acceleration includes; why full contact feedback failed; and why reduced-kernel, spiking, Syn6, and CAD results cannot inherit the selected human result.

Do not replace the final conclusion with a catalog of repaired bugs. Tie the actuator, joint-design, and neural-tool contributions back to constructing and testing biomimetic robots.

## Delivered module and visual handoff

neural_content.mjs is complete and passes node --check. It replaces exactly twelve existing main positions (57–68), updates 71, 72, and 74 to follow the new evidence, and adds six question slides: same-controller stages, assistance decomposition, historical NumPy–Simulink agreement, the loaded-arm reflex, RF motor/force maps, and RF Ia/II maps. No controller, model, original source figure, or dissertation file was edited, and no simulation was run.

The neural RG/PF/MN and contact-strength figures are faithful landscape arrangements of original vector axes and data paths. Original chart colors, dashed lines, ordinates, and time ticks are preserved. DERIVATIVE_MANIFEST.json records source hashes and the selected axes or exact crop viewport. PNG previews were rendered using Sharp because the installed MuPDF SVG conversion dropped dash patterns. The final PowerPoint render must retain the source dashes as well.

Root is rendering the packaged deck. The independent Reddy-perspective reviewer owns the final visual pass. Priorities: confirm main 65's complete scientific schematic is sufficiently large; check the native crop on main 63 retains or adds an explicit stimulus/baseline legend; confirm the new EB cropped vector diagrams render on white and remain clearly labeled as the standing export. These are visual handoff checks, not verified final-render failures. I inspected all derived scientific-panel PNGs individually before handoff; the final deck itself has not yet been rendered by this subtask.

## Bounded render fixes delivered after review

Root requested four replacement assets, now available as SVG, PDF, and PNG in assets/neural/slide_panels. Each PNG was inspected individually. The first asset has aspect ratio 3.125:1; the other three are exactly 16:9. RENDER_FIX_PROVENANCE.json records exact source/output hashes, clocks, scope, and unchanged amplitudes. No model/controller construction or integration occurred.

- NEURAL_RIGHT_HAM_HIERARCHY uses actual aligned RG, RF/HAM PF, and their downstream BFLH motor-neuron voltage. BFLH retains its original orange line and the original bilateral hip-MN panel ordinate range. These are matched post-update records; no inferred RF voltage is added.
- CONTACT_PF_AND_FOOT_LOADING_WIDE uses the original comparison's FIGURE_DATA.npz arrays. The title is “Ankle pattern-formation activity,” Plantarflexor and Dorsiflexor are spelled out, and both legends are outside axes. Shared row scales, dashed lines, and the distinct native PF/force clocks are preserved.
- FOOT_LOADS_AND_ASSISTANCE_WIDE uses the same aligned arrays and 768.756726 N denominator as the published source panel. The legend is outside the graph and distinguishes foot force, vertical restraint, and prescribed lift.
- DENG_RG_STIMULUS_COMPLETE restores all time ticks and the Baseline/Stimulus legend using the original saved baseline.npz and RG_plus50_2cycles.npz. The existing cropped source PNG had already lost its x ticks. New plotting preserves measured/simulated saved voltages, source-relative time, the +50 nA pulse, and its 3.876 s interval.

Channel-availability correction: the compact ground9 records lack RF MN, but the separate ground9_rc_s4.mechanics.npz stores all 1,033 actual cells, including RF, BFLH, and semimembranosus MNs, at 50 ms pre-update clocks. This was reported to root immediately after finding the full sidecar. Those clocks must not be silently mixed with the 2 ms aligned post-update trace. saved_mn_channel_inventory.json records both scopes. The requested BFLH-only replacement consistently uses the aligned record.

## Projection-sized replacements after the second PowerPoint pass

The independent review found small text when the figures were placed in the actual slide frames. Root requested larger source fonts and matched aspect ratios. New SVG/PDF/PNG files are staged, individually inspected, and recorded in PROJECTION_RENDER_PROVENANCE.json:

- DENG_RG_STIMULUS_PROJECTION: 6.4 × 3.8 in, 16 pt axes and key, target frame 525 × 268 pt.
- CONTACT_PF_AND_FOOT_LOADING_PROJECTION: 10 × 4 in, 18 pt axes and key, target frame 870 × 333 pt. Full plantarflexor/dorsiflexor labels remain outside axes. The two plots per condition retain their original scales and clocks.
- FOOT_LOADS_AND_ASSISTANCE_PROJECTION: 5.2 × 4.2 in, 18 pt axes and 17 pt key, target frame 380 × 295 pt. “Vertical forces” scopes the shorter Right foot, Left foot, Restraint, and Lift keys. Time still means the actual recorded force-interval end.
- W2L_PATHWAY_REMOVALS_PROJECTION: 10.5 × 4.5 in, 18 pt axes and external Left/Right key. Recovered the original generator and exact JSON in Code/MuJoCo_SNS/spinal/reports_20261027_codex. The original seven labels/order, Left blue bars, Right striped red bars, and contact-event definition are retained. Left counts are 33, 24, 17, 35, 25, 19, and 20; Right counts are 26, 22, 15, 34, 2, 20, and 31. These are historical supported-W2L events, not current ground9 steps.

No shared deck module was edited. Root owns insertion and final PowerPoint verification.
