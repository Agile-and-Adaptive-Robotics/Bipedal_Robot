# ZCode pickup: defense deck and W2L recording models

## Start here

This handoff is for continuing the 31 October 2026 defense-deck and W2L work on **easteregg2**. Run `hostname` first, then read the repository `AGENTS.md`. These absolute paths are valid on easteregg2. Do not interpret a missing `D:` path on another machine as a stale file.

Do not clean, reset, move, or regenerate the whole working tree. Several active results are ignored by Git, and the repository also contains unrelated dirty work from parallel dissertation and BPA sessions. Ben commits and pushes. Agents do not.

Recommended pickup order:

1. Read this file and the records named below.
2. Treat defense revision 17 as the current slide baseline. Do not rebuild it merely to prove reproducibility.
3. Open the requested AnimatLab project only after confirming Ben no longer has an unsaved AnimatLab window. The prior GUI-open check was deliberately skipped because the existing window contained unsaved edits.
4. If Ben records a corrected walker or passive-guide run, make a new `v18` deck from `v17`. Do not overwrite `v17`.
5. Render and inspect every slide in any new revision, verify native video playback in PowerPoint, export the full PDF including hidden question slides, and retain only the last two validated revisions. While building `v18`, keep `v16` until `v18` passes; then retain `v17` and `v18`.

## Current defense deck

Authoritative pointer:

```text
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\current_deck.json
```

Current deliverables, both intentionally Git-ignored:

```text
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\C3NS_Revision_20261031_v17\
├── output\
│   ├── Bolen_PhD_Defense_RoseGreen_20261031_v17.pptx
│   ├── Bolen_PhD_Defense_RoseGreen_20261031_v17.pdf
│   ├── AJH_Review_Response_v17.txt
│   └── sadb_app.html
├── build_deck.py
├── add_recordings.py
├── finalize.mjs
├── render_review.ps1
├── check_video_playback.ps1
├── assets\                         copied figures, posters, and embedded media
└── build\
    ├── main_story.json
    ├── slide_order.json
    ├── build_summary.json
    ├── validation_delivery5.json
    ├── preservation.json
    ├── visual_review.json
    ├── comment_compliance.json
    ├── native_chart_owners.json
    ├── native_recordings.json
    ├── native_video_playback.json
    ├── native_text_geometry.json
    ├── geometry_flags.json          expected: []
    ├── render\slide_01.png ... slide_123.png
    └── contact\sheet_001.jpg ... sheet_118.jpg
```

Final `v17` facts:

- PPTX SHA-256: `9d9375e2f83cec31a40ef782edc530c4c4b1173aa3765bc29cb6d00c5ba61e09`
- PPTX size: 94,467,329 bytes, below GitHub's 100 MB per-file limit but intentionally ignored.
- 123 slides: 70 main-talk slides and 53 hidden question slides.
- Planned main-talk time: 2,031 seconds, or 33 minutes 51 seconds. Hard cap: 45 minutes.
- All 123 slides were rendered and visually inspected. The PDF has 123 pages.
- Rose-green C3NS template, original honeycomb title treatment, smaller original logo/footer, and original plot colors are preserved.
- SADb is last, on hidden question slides 117–123.
- Slide 52 contains two native MP4s with autoplay, mute, and loop verified in PowerPoint.
- Slide 74 contains Ben's newer bilateral-RG Root-fixed GIF.
- The media already in `v17` predate the corrected body/joint geometry. Do not label them as corrected-geometry evidence.
- The new flexor test photo remains pending.
- The separate dissertation revision of the Hunt 2019 citation remains deferred. The correct source is A. J. Hunt, *Real-time control of dynamic walking in a dog-like robot*, 2nd Intel Neuromorphic Research Community meeting, Intel Labs, 15 February 2019.

Previous validated revision, retained under the two-revision rule:

```text
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\C3NS_Revision_20261030_v16\output\Bolen_PhD_Defense_RoseGreen_20261030_v16.pptx
SHA-256 cfa20a90ee6216bb36acb9737586d91d6e93ff3a8a7a62e1c28c693fd1da59cd
```

Advisor review inputs and extracted comments:

```text
C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_PhD_Defense_C3NS_20261030_v13_review_AJH.pdf
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\AJH_Review_20261030\
├── comments.json
└── REVIEW_STATUS.txt
```

Local slide references:

```text
C:\Users\Ben Bolen\Downloads\19_02_15 - INRC.pptx
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Movement and Control of Biomimetic Humanoid Robots.pptx
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\
D:\GitHub\Bipedal_Robot\Pictures\PictureAnalysis\20220308_032827.jpg
D:\GitHub\Bipedal_Robot\Pictures\PictureAnalysis\20230419_225929.jpg
D:\GitHub\Bipedal_Robot\Pictures\PictureAnalysis\20220308_031822.jpg
D:\GitHub\Bipedal_Robot\Pictures\PictureAnalysis\20220307_222908.jpg
D:\GitHub\Bipedal_Robot\Code\Matlab\Colors.m
```

For a new revision, copy the `v17` source tree to a new `v18` folder, make `v17` the baseline, change all hard-coded revision paths and receipt names, and write a new final file. `finalize.mjs` will not safely reuse an existing final/receipt as a new validation event. Preserve the package resources and native charts. Do not flatten the deck through PDF or screenshots.

## Local recordings

These files are outside the repository and will not survive a clone:

```text
C:\Users\Ben Bolen\Videos\Screen Recordings\
├── frozen-root treadmill.mp4          embedded on slide 52
├── free_drop.mp4                      embedded on slide 52
├── w2l_bilateralRGfreeze.gif          used on slide 74
├── frozen-root treadmill.gif          source alternative, not the final embedded movie
├── free_drop.gif                      source alternative
└── Screen Recording 2026-10-31 005831.mp4   inspect before using; not identified as final media

C:\Users\Ben Bolen\Downloads\2021-02-17-20-58-24.mp4   older walker recording
```

The copies actually embedded into `v17` also exist under `C3NS_Revision_20261031_v17\assets\`, so the final PPTX is self-contained. The external originals are still useful for recutting or rebuilding.

## AnimatLab geometry models

These are new copies. The original models and Ben's recordings were not overwritten. The `.aproj`, `.asim`, `.aform`, traces, and most build artifacts are ignored locally by each folder's `.gitignore`. Preserve the folders before cleaning or recloning.

Frozen-Root corrected geometry, original shared-RG controller:

```text
D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_Geometry_20261031\
├── W2L_Geometry.aproj                 open this in AnimatLab
├── W2L_Geometry.asim
├── W2L_Geometry_Standalone.asim       native runner input
├── README.txt
├── ANATOMICAL_LIMITS.md
├── GEOMETRY_REPORT.json
├── NATIVE_TEST.json
├── PRESERVATION.json
└── tools\
    ├── build_geometry.py
    ├── inspect_geometry.py
    └── test_geometry.py
```

Frozen-Root corrected geometry, newer type-II/contact/separate-LH-RH-RG controller:

```text
D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_Geometry_Newer_20261031\
├── W2L_Geometry.aproj                 open this in AnimatLab
├── W2L_Geometry.asim
├── W2L_Geometry_Standalone.asim
├── README.txt
├── ANATOMICAL_LIMITS.md
├── GEOMETRY_REPORT.json
├── NATIVE_TEST.json
└── PRESERVATION.json
```

Geometry result of record:

- Left hip 10° extended with 10° knee flexion. Right hip 20° flexed with 55° knee flexion. Ankles begin at 0° anatomically.
- The old tibial knee anchor was 11.514 mm anterior to the distal femur. Correct placement uses `p_child = p_parent_anchor - R_child q_child_anchor`.
- Both knees begin with coincident markers. Maximum native dynamic marker separation is 0.549 mm for the original and 0.585 mm for the newer model.
- Solid-foot clearance above the moving belt is 10 mm.
- Native 10-second tests pass. The corrected right knee extends from 55° to about 2.2°, with a small soft-limit overshoot below the nominal 3° anatomical limit.
- Root stays frozen in these two geometry projects. These are recording models, not unsupported-walking claims.

## AnimatLab passive-support models

Simple hip-wing/rail copies requested for recording:

```text
D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_20261031\
├── W2L_VirtualWalker.aproj            original shared-RG controller
├── W2L_VirtualWalker_Standalone.asim
├── README.txt
├── TEST_RESULTS.json
├── VERIFICATION.json
├── BUILD.json
├── PRESERVATION.json
├── verification_comparison.png
├── SUBNETWORK_EXPERIMENT_RECORD.md
├── SUBNETWORK_REBUILD.json
├── SUBNETWORK_VERIFICATION.json
├── subnetwork_review\                 SVG review of saved GUI pages
└── tools\                             all builders, runners, and verifiers

D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\
├── W2L_VirtualWalker_Newer.aproj      newer 96-neuron controller
├── W2L_VirtualWalker_Newer_Standalone.asim
├── README.txt
├── TEST_RESULTS.json
├── VERIFICATION.json
├── BUILD.json
├── PRESERVATION.json
├── SUBNETWORK_REBUILD.json
├── SUBNETWORK_VERIFICATION.json
└── subnetwork_review\
```

The simple copies have an unfrozen pelvis, passive lateral hip wings, and invisible horizontal rails. Native 10.1-second tests completed at normal playback speed. They remain mechanically supported and drift backward with the treadmill. Do not call this unsupported walking. The editable project pages were subsequently reorganized into the historical uniarticular hierarchy. The runtime `.asim` files remained byte-identical during that layout-only pass. Read `SUBNETWORK_EXPERIMENT_RECORD.md` before rebuilding, because the old physical builder recreates flat pages.

More realistic finite-clearance native guide companions:

```text
D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Guided_20261031\
├── W2L_VirtualWalker_Guided.aproj     canonical 10 mm, source clock
├── W2L_VirtualWalker_Guided_Standalone.asim
├── EXPERIMENT_RECORD.md               read first
├── BUILD_GUIDE.json
├── NATIVE_GUIDE_SMOKE.json
└── GUIDE_VERIFICATION.json

D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Guided30_20261031\
└── W2L_VirtualWalker_Guided30.aproj   30 mm clearance

D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_GuidedT2_20261031\
└── W2L_VirtualWalker_GuidedT2.aproj   10 mm, doubled neural time constants
```

The guide carriage uses physical contacts, a finite ±1.9 m fore-aft slider, and no root weld or continuous balance controller. The robot Root remains free. These native runs are finite, supported trials with no sampled knee hyperextension, not a completed walking or curriculum pass. If Ben wants a screen recording with bounded wobble, start with `W2L_VirtualWalker_Guided.aproj`; if he wants the simplest visual hip-wing support, start with `W2L_VirtualWalker.aproj`.

Newer-controller engineering candidates and diagnostics:

```text
D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_SIContact_20261031\
├── W2L_VirtualWalker_Newer_SIContact.aproj
└── EXPERIMENT_RECORD.md               20 nA/contact engineering candidate; falls

D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_GUI1nA_20261031\
├── W2L_VirtualWalker_Newer_GUI1nA.aproj
└── EXPERIMENT_RECORD.md               1 nA/contact interpretation; falls

D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_NewerV3Guided_20261031\
├── W2L_VirtualWalker_NewerV3Guided.aproj
├── V3_DIAGNOSTIC_BUILD.json
├── MUJOCO_REFERENCE_V3_DIAGNOSTIC.json
└── NATIVE_GUIDE_SMOKE.json
```

Do not silently substitute these candidates for Ben's supplied newer model. The editable newer project has 96 neurons and 204 synapses, while the supplied standalone had 96 neurons and 208 synapses with different contact paths, types, tonic drive, and polynomial coefficients. The two contact-gain candidates are explicit engineering interpretations, not identified optima. The newer branch has not passed the full supported-stepping gate.

## MuJoCo port and continuation work

The corresponding W2L MuJoCo package is local and partly ignored:

```text
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\
├── README.md                          commands and output layout
├── EXPERIMENT_RECORD.md               authoritative claims and limitations
├── ORIGINAL_TIMING2_SAGITTAL.json
├── ORIGINAL_TIMING2_PHYSICAL_GUIDES.json
├── run_virtual_walker.py
├── native_kernel.py
├── import_virtual_walker.py
├── contact_guides.py
├── assess_trial.py
├── draw_trial.py
├── render_trial.py
└── runs\                              raw XML, NPZ, GIF, PNG, and JSON; Git-ignored
    ├── GUIDE_STUDY_V1.json
    ├── GUIDE_EXTREMES_V1.json
    ├── GUIDE_DEADBAND_VERIFICATION.json
    ├── GUIDE_LONG_COMPARISON_V1.json
    ├── NEWER_GUIDE_PAIRS_V1.json
    ├── BENCHMARK_FINAL.json
    └── FRESH_GUIDE_CLI_PARITY.json
```

Read `EXPERIMENT_RECORD.md` in full before continuing. Current bounded result:

- The original 80-neuron/150-edge source controller has a verified direct port and supported treadmill-stepping checkpoints.
- Physical guides reduce wobble but still provide about 85% body-weight vertical assistance through rails. This is not independent standing or harness-free walking.
- The corrected five-stage curriculum and Hunt four-part training sequence are incomplete.
- The newer GUI-authoritative branch remains latched or fails the full supported-stepping gate under the tested candidates.
- The record's current next steps are a basin test of the retained physical-guide checkpoint, explicit small reductions in rail support, standalone standing/Hunt gates, and preservation of every failed case. Keep the knee nonhyperextension gate.
- A source-2/0-tonic plus V3-output-zero control was still being distinguished from tied-0/0 diagnostics. Check the newest JSON in `runs\newer_rg_isolation_v1\` before rerunning or claiming causation.

Useful ignored media from the MuJoCo package:

```text
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\timing2_fresh_ground\SUPPORTED_STEPPING_REPLAY.gif
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\timing2_fresh_ground\SUPPORTED_STEPPING_STILL.png
```

Those files show the continuous sagittal-guidance checkpoint, not the later physical-guide condition. Preserve that distinction if they enter the deck.

## Skills and runtimes

Shared ZCode skill source:

```text
C:\Users\Ben Bolen\Documents\GitHub\ZCode_Skills\animatlab\SKILL.md
C:\Users\Ben Bolen\Documents\GitHub\ZCode_Skills\mujoco\SKILL.md
C:\Users\Ben Bolen\Documents\GitHub\ZCode_Skills\neuroscience\SKILL.md
C:\Users\Ben Bolen\Documents\GitHub\ZCode_Skills\python\SKILL.md
```

The active junction for the AnimatLab skill is `C:\Users\Ben Bolen\.agents\skills\animatlab\SKILL.md`. Use the host's built-in `Presentations` skill for any deck revision.

Executables on easteregg2:

```text
D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatLab2.exe
D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe
D:\Anaconda\envs\myo\python.exe
C:\Users\Ben Bolen\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\bin\node.exe
C:\Users\Ben Bolen\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe
C:\Users\Ben Bolen\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe
D:\Anaconda\envs\myo\Library\bin\ffmpeg.exe
```

PowerPoint COM is required for the native deck render, PDF export, and video playback check.

## Git and preservation warnings

- The defense `output\`, `build\`, and copied `assets\` trees contain ignored local artifacts.
- Each W2L model folder has a local `.gitignore` that ignores `.aproj`, `.asim`, `.aform`, traces, and other generated files while allowing selected records such as `README.txt`.
- The raw MuJoCo `runs\` tree is ignored.
- The working tree contains unrelated modified and untracked dissertation, SNS, and BPA files. Do not use a repository-wide clean, reset, checkout, or deletion.
- `D:\GitHub\Bipedal_Robot_Ignored\` contains older moved build trees and scratch. It is archival, not the active source for this pickup.
- No commit or push was performed for this work.

## Definition of a safe completed pickup

For a deck-only continuation, completion means a new revision beside `v17`, all slides rendered and individually inspected, geometry and package checks passed, PDF page count matched, media played and looped in PowerPoint, current pointer updated, and only the last two validated final PPTX revisions retained.

For an AnimatLab continuation, completion means the chosen `.aproj` opens without changing Ben's source projects, reset/play works at normal speed, root-Y and joint/MN telemetry remain available, the exact controller/support condition is named, and any recording or plot avoids unsupported-walking claims.

For a MuJoCo continuation, completion means fresh-process verification, exact configuration and hashes saved, failed cases retained, no knee hyperextension, support forces stated explicitly, and no claim beyond the gates documented in `EXPERIMENT_RECORD.md`.
