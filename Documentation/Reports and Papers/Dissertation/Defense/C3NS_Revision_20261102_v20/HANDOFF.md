# Defense v20 handoff

Date: 2 November 2026. Machine: easteregg2. Status: **v20 is the current delivery; final PPTX/PDF, package checks, static rendered review, and bounded media checks are complete**. The parent editor promoted `Defense/current_deck.json` after media verification. Spoken timing remains unrehearsed.

## Scope and intended result

The final RoseGreen/C3NS defense presentation contains **127 slides: 74 main slides and 53 hidden question slides**. Main slides are positions 1–74; the final package verifies hidden positions 75–127. Pacing is **37 minutes at the 30-second-per-slide baseline, a 39-minute delivery target, and a hard 45-minute limit**. These are planning values; the talk has not received a spoken timing rehearsal.

The main narrative now follows:

1. Terrain and robotics motivation; BPA rationale; force/moment-arm mechanics; sensory feedback; anatomical correspondence.
2. Limited contraction constrains joint range, and historical BPA routes do not reproduce human-model torque.
3. The reciprocal biology–robot research loop, potential prosthetic application, and three main-talk gaps.
4. Published BPA force characterization: model deficiencies, experiment, maximum force, pressure/contraction normalization, and comparisons.
5. Isometric knee torque: experimental mismatch, offset/compliance/bending methods, bracket frames, inner displacement solve, correction-factor identification, corrected predictions, physical shortfall, route redesign, actual route GIFs, and physical knee results.
6. Models and tools in the order of dissertation Table 6.1, with separate model identities and test conditions.
7. Mechanical synthesis, dynamic and supported-control next steps, and the three-contribution conclusion.

**SADb is excluded from the main talk**, including the roadmap and rewritten conclusion notes. Its seven question slides remain at the end, positions 121–127. The native offline walkthrough on source slide 123/final 127 is preserved; it must not be replaced by its imported poster image.

## Governing files and source identity

- `revise.mjs` is the entry-point authoring source. It imports v19, applies the initial edits, then calls `pass2_content.mjs` and `pass3_content.mjs` in that order.
- `pass2_content.mjs` supplies the new opening, mechanics sequence, identification slide, enlarged native GIFs, and scientific scope labels. It calls `sadb_questions_fix.mjs` for the bounded question-slide changes.
- `pass3_content.mjs` applies the latest first-use wording, figure definitions, source-aware notes, legend/layout repairs, and rewritten conclusion.
- `make_equations.py` supplies the displayed equation images. In particular, the displacement root uses the reference factor `F620(l0)` times the dimensionless BPA response, avoiding a second pressure normalization.
- `package_current.py` is the current packaging path. **Do not use the older `package_deck.py` as the authority for the current order.** The current packager retains original native chart/workbook, movie, master, and layout parts and asserts their preservation.
- `build/slide_order.json` governs order. `build/package_report.json` records the source-to-final mapping produced by packaging. `build/changed_slides.json` records changed authoring source IDs, not final positions.
- Scientific source identity is the submitted dissertation and the frozen `../C3NS_Revision_20261031_v18/submitted_source` / `submitted_figs` copies. The submitted PDF is `C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_Dissertation.pdf`. Do not substitute stale historical AGENTS.md numbers or overwrite the manuscript source.
- Advisor evidence: `../AJH_Review_20261030/comments.json` and the reviewed advisor PDF. The opening argument and evidence map are in `build/OPENING_NARRATIVE_SUPERVISOR.md`.

The latest pass-3 source has been rebuilt into the final output. All 127 images from the final-output rerender exactly match the inspected pass-3 images (`build/final_render_comparison.json`: no differences). The final PPTX is **90,314,963 bytes**, SHA-256 **`69ea30c4fe7ab533f1c8a74b258cf431ad3025c4c4b8a946215916393aa9d272`**. The review PDF contains **127 pages**.

## ZCode source and path guide

```text
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\
├── Defense\
│   ├── current_deck.json                  [promoted to verified v20]
│   ├── SESSION1_NEURAL_CONTROLLER_HANDOFF.md
│   ├── AJH_Review_20261030\comments.json
│   ├── C3NS_Revision_20261031_v18\
│   │   ├── submitted_source\chapters\    [frozen submitted scientific source; retain]
│   │   ├── submitted_figs\               [original submitted figure assets; retain]
│   │   └── transplant_v18.py             [packaging helper read by v20; retain]
│   ├── C3NS_Revision_20261102_v19\
│   │   ├── assets\                       [extracted figures and approved cutout; retain]
│   │   └── output\Bolen_PhD_Defense_RoseGreen_20261102_v19.pptx
│   │                                      [v20 import/native-preservation baseline]
│   └── C3NS_Revision_20261102_v20\
│       ├── HANDOFF.md
│       ├── DIRECTIVE_COMPLIANCE.md
│       ├── run_build.mjs                  [calls revise.mjs]
│       ├── revise.mjs                    [imports v19, then pass 2 and pass 3]
│       ├── pass2_content.mjs
│       ├── pass3_content.mjs
│       ├── sadb_questions_fix.mjs         [source 117/118/120/121 only]
│       ├── make_equations.py
│       ├── package_current.py            [current packager; not package_deck.py]
│       ├── finalize.mjs
│       ├── render_review.ps1
│       ├── verify_delivery.py
│       ├── check_video_playback.ps1
│       ├── assets\                       [v20 diagrams/equations and sourced images]
│       │   ├── slide17_comparison_*       [data, provenance, and regeneration source]
│       │   ├── xiFrameTransform.png
│       │   ├── xiProjection.png
│       │   └── force_balance.png
│       ├── build\
│       │   ├── supervisor_sources\       [older-talk evidence used by the opening]
│       │   ├── render\slide_001.png … slide_127.png
│       │   ├── slide_order.json / package_report.json / changed_slides.json
│       │   ├── validation_delivery.json / delivery_checks.json
│       │   ├── final_render_comparison.json / pass3_all_hashes.json
│       │   └── committee_*_review.md / SUPERVISOR_ACCEPTANCE_CHECKLIST.md
│       └── output\
│           ├── Bolen_PhD_Defense_RoseGreen_20261102_v20.pptx
│           └── Bolen_PhD_Defense_RoseGreen_20261102_v20.pdf
├── Archive\session1_submission_20261031\source\
├── Figures\                             [original dissertation figure sources]
└── ProofFinal\                           [reference only for this deck pass]

D:\GitHub\Bipedal_Robot\Knee_Flexor_20mm.gif
D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Knee_Extensor_20mm.gif
D:\GitHub\Bipedal_Robot\Pictures\PictureAnalysis\20220307_222908.jpg
C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_Dissertation.pdf
```

The standalone v18 final PPTX has been deleted; its original source/figure assets and packaging helper are retained. **Retain the v18 source directory.** The current v20 build reads the retained v19 PPTX, v19 assets, v18 submitted figures, the v18 helper, v20 assets, and the external repository GIF/photo paths shown above. It does not require the deleted standalone v18 PPTX. Rebuilding v19 from its old source is a separate task and may still have that historical dependency; rebuild v20 directly from the retained v19 baseline.

The v20 folder is therefore not a standalone source bundle. Preserve the named sibling/source dependencies. Cleanup is complete: the parent removed the v18 final PPTX, v20 `build/artifact_draft.pptx`, v20 `build/candidate.pptx`, and the stale `build/render/slide_128.png`. `build/cleanup.json` records the four deletions, totaling 307,614,121 bytes. **Final PPTX revisions 19 and 20 are retained; original v18 source assets remain present.**

## Critical source-to-final mapping

The final order omits source slides **14 and 50** to remove redundant bridges. Six new authoring IDs, 124–129, are inserted into the main sequence. The mapping below agrees with the final package report.

| Source ID | Final position | Content |
| --- | ---: | --- |
| 1–13 | 1–13 | New opening, section divider, and BPA-length definition |
| 15–27 | 14–26 | Force characterization and comparisons |
| 28–37 | 27–36 | Torque section, setup, discrepancies, and correction terms |
| 124 | 37 | Two-rotation bracket frame |
| 125 | 38 | Inner force/displacement solve |
| 129 | 39 | Fit/validation stages and shared stiffness pair |
| 38–42 | 40–44 | Corrected pinned results, fit table, and biomimetic-knee comparison |
| 126 | 45 | Single-BPA physical shortfall |
| 43 | 46 | Route design variables, torque requirement, and constraints |
| 127 / 45 | 47 / 48 | Flexor GIF / physical results |
| 128 / 46 / 47 | 49 / 50 / 51 | Direct extensor GIF / cam explanation / physical results |
| 44 / 48 | 52 / 53 | Separate 27-site model comparison / mechanical synthesis |
| 49, 51–64 | 54–68 | Neural section and Table 6.1 model/tool flow |
| 65–70 | 69–74 | Synthesis, future tests, and conclusion |
| 77 / 78 | 81 / 82 | Corrected duty-side labels / own-peak force comparison |
| 101 / 103 | 105 / 107 | Neural control labels / separate-model distinction |
| 117–123 | 121–127 | SADb questions; source 123 preserves native offline video |

## Scientific distinctions to retain

- χ1 and χ2 are effective system-stiffness parameters. They are not uniquely identified bracket material properties. The extensor retains the flexor stiffness pair and fits its own χ0 and χ3.
- The inner scalar displacement solve, factor identification, and geometric route optimization are three different operations. The selected validation tests also informed Pareto acceptance; do not describe them as blind external validation.
- Free contraction still requires individual measurement. The published scatter fit does not establish a reliable prediction from initial length.
- Force-model sample counts and the 80%/20% split stay in notes, per Ben's direction. The selected-length comparisons do not establish a separately held-out-actuator study.
- The physical knee claims apply at the measured isometric angles. The route GIFs illustrate geometry; they are not physical motion tests or proof of byte-identical measured-test geometry.
- The 27-site result is a separate model comparison at 620 kPa using actuator sizing and pulley ratios on inherited routes. Its percentages are equal-site means of angular-domain coverage, not fractions of equally weighted sample points or physical whole-leg results.
- W2L108's 9/8 swings and 40.5% body-weight assistance belong to the newer supported demonstration. The recovered-W2L 0.498/0.223 feedback comparison is a separately specified experiment. Its controller, drive, and support are fixed across conditions.
- The +50 nA, 3.876 s stimulus plot is a reduced fourteen-neuron Deng-derived neural-response test. It is not the response of the complete imported rat model or the larger W2L108 model.
- No harness-free walking result is claimed. Dynamic hardware actuation and reduced-support control remain future tests.

## Review history and latest repairs

Four independent committee-perspective reviews are preserved in `build/committee_hunt_review.md`, `committee_jiang_review.md`, `committee_reddy_review.md`, and `committee_turcic_review.md`. They are analytical perspectives, not statements by the actual committee. They refer to different candidate snapshots; use their source IDs and dated/iteration sections rather than treating every historical issue as still open.

`build/SUPERVISOR_ACCEPTANCE_CHECKLIST.md` records the second rendered inspection. That inspection covered final positions 1–12, 38, 39, 42, 44, 45, 46, 47, 49, 52, 60, 74, and 82. It confirmed the causal opening and identified specific defects. The latest pass-3 source addresses the pinned-photo/legend collision, the own-peak-normalization/legend collision, the ambiguous biomimetic prediction legend, equation definitions, fit-objective notes, route-direction wording, first-use BPA expansion, reciprocal-loop wording, and conclusion evidence.

The completed second-pass static coverage was: **1–12 supervisor, 13–53 Turcic-perspective reviewer, 54–105 Reddy-perspective reviewer, and 106–127 parent editor**. After pass 3, the parent re-opened all **15 changed images: 3, 9, 23, 38, 39, 42, 44, 46, 47, 50, 57, 74, 82, 105, and 127**. The supervisor independently rechecked **38, 39, 42, 44, 46, 74, and 82**; the Reddy-perspective reviewer independently rechecked **57, 74, 82, and 105**. The specified photo/legend collisions, normalization-label collision, and model-scope defects are resolved in those images. The final-output rerender matches all 127 inspected pass-3 images exactly.

This is combined review coverage, not a claim that every reviewer opened every slide. Static review does not establish animation playback or spoken timing. `build/SUPERVISOR_FAILURE_LOG.md` preserves the withdrawal of the former v19 blanket design judgment.

## Verification and delivery state

| Check | Final recorded state |
| --- | --- |
| Order/count and hidden flags | Verified: 127 total; main 1–74; hidden 75–127 |
| Latest source authoring and packaging | Complete; final PPTX hash and size recorded above |
| Main-talk SADb exclusion and rejected wording | Final main text and linked notes scan clear; no missing slide relationships |
| Package integrity and first-party import | Pass, 127 slides; no structural findings |
| Native charts and embedded workbooks | All **7 charts and 7 workbooks** unchanged from the baseline, on owner slides 15, 25, 30 (two), 52, 101, and 102 |
| Template/title preservation | Original title, masters, layouts, themes, and template-image bytes unchanged |
| Full final PowerPoint render | Complete; all 127 images exactly match inspected pass-3 images |
| Latest repaired slides visually inspected | Complete under the combined coverage and explicit rechecks above |
| GIF crop mask | Narrow check recorded: 199 flexor frames, zero colored geometry pixels inside the small axis-fragment mask; this does not certify all crop bounds or playback |
| Native video/media preservation | Four original MP4 assets retained; three active clips remain attached to slides 56 (two) and 127 (one); one was already inactive. Seven GIF parts retain multi-frame content |
| Route-GIF movement | `build/native_gif_playback.json` verifies movement on final 46, 47, and 49 in a native PowerPoint video export, sampled at six times across a source GIF cycle. The parent also checked two exported frames per slide for clarity. This is bounded to those three slides |
| Native MP4 playback | `build/native_video_playback.json` verifies both final-56 clips advance and restart their loops. Final 127's offline SADb walkthrough advances with its original click-to-play, no-loop behavior |
| Submitted PDF conservation | Verified unchanged, SHA-256 `447c2cc5a9dd219fa0857e298727fb68d4797af8f406251865eb22f4b126361f` |
| PDF includes all 127 slides | Verified: 127 pages |
| Actual spoken duration | Unrehearsed; 37-minute baseline / 39-minute target / hard 45-minute limit |
| Current-delivery pointer | Promoted by the parent to v20, with the final hash and evidence paths |
| Cleanup | Complete; final v19 and v20 PPTX files retained, v18 source assets retained, and obsolete final/draft/candidate/stale-render files removed |
| Commit/push or Overleaf upload | None authorized or performed by this handoff task |

The validation helpers now use the 127-slide order and mapped chart-owner positions. `validation_delivery.json` reports zero blocking findings and three nonblocking layout warnings on slides 36, 41, and 60; the PowerPoint-rendered versions are covered by the static inspection above. The arithmetic helper skipped tables without unambiguous totals; it is not a recomputation of the experimental values. The completed native-video report separately verifies playback on actual owner slides 56 and 127.

Current deliverables are `output/Bolen_PhD_Defense_RoseGreen_20261102_v20.pptx` and its matching `.pdf`. `build/validation_delivery.json`, `build/delivery_checks.json`, and `build/final_render_comparison.json` provide package/static evidence. Separate `build/native_gif_playback.json` and `build/native_video_playback.json` provide the bounded movement/playback results above. All requested delivery checks are complete; spoken timing remains unmeasured. A future new physical photograph would be an optional user-directed revision, not an outstanding delivery requirement.

See `DIRECTIVE_COMPLIANCE.md` for the per-directive implementation and remaining checks. No manuscript, model, database, or result file was edited by the handoff/documentation task.
