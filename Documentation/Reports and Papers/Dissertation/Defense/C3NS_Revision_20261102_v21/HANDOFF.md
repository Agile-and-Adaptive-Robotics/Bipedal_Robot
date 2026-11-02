# Defense v21 handoff

Machine: easteregg2. The authoritative completion status and artifact hashes are in `../current_deck.json`. Do not promote a partial build by changing that pointer manually.

The deck has 77 main slides and 68 hidden questions. The 30-second baseline is 38.5 minutes; allow about 40 minutes, with a hard 45-minute limit. A spoken rehearsal is still needed. Suggested checkpoints: finish the force section at about 13–14 minutes (slide 28), the mechanics/AnimatLab transition by about 29 minutes (slide 59), and the neural demonstrations by about 35 minutes (slide 71). The conclusion is slide 77. SADb is only on 139–145.

## Relevant files

```text
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\
├── current_deck.json
├── SESSION1_NEURAL_CONTROLLER_HANDOFF.md
├── AJH_Review_20261030\comments.json
├── C3NS_Revision_20261031_v18\
│   ├── submitted_source\chapters\     frozen submitted scientific source
│   ├── submitted_figs\                submitted figures
│   └── transplant_v18.py              packaging dependency; retain
├── C3NS_Revision_20261102_v19\assets\ source assets; retain
├── C3NS_Revision_20261102_v20\
│   ├── assets\                       retained baseline assets
│   └── output\Bolen_PhD_Defense_RoseGreen_20261102_v20.pptx
└── C3NS_Revision_20261102_v21\
    ├── output\
    │   ├── Bolen_PhD_Defense_RoseGreen_20261102_v21.pptx
    │   └── Bolen_PhD_Defense_RoseGreen_20261102_v21.pdf
    ├── DIRECTIVE_COMPLIANCE.md
    ├── HANDOFF.md
    ├── revise.mjs                     imports v20 and applies the modules
    ├── force_content.mjs
    ├── mechanics_content.mjs
    ├── neural_content.mjs
    ├── backup_content.mjs
    ├── visual_fixes.mjs                final review corrections
    ├── package_current.py             preserves original native PPT parts
    ├── finalize.mjs                   first-party export validation
    ├── render_review.ps1              native render + PDF
    ├── verify_delivery.py
    ├── check_video_playback.ps1
    ├── export_gif_test.ps1
    ├── assets\
    │   ├── force\                    source/style conservation and plots
    │   ├── mechanics\                native figures, data, and generators
    │   ├── identification\           saved Pareto candidates and equations
    │   └── neural\
    │       ├── ASSET_MANIFEST.json
    │       ├── eb475ws4_motor_loops\  historical standing-export drawings
    │       ├── rfham_selected\       actual ground9 selected circuit
    │       ├── ground9_rc_s4\        preserved animation/records
    │       └── slide_panels\         presentation derivatives + provenance
    └── build\
        ├── render\slide_001.png … slide_145.png
        ├── slide_order.json           baseline/source IDs in final order
        ├── slide_map.csv              original and final positions
        ├── structure.json
        ├── MECHANICS_SOURCES.md
        ├── validation_delivery.json
        ├── delivery_checks.json
        ├── visual_review_coverage.json
        ├── committee_*_v21_review.md
        ├── native_video_playback.json
        ├── native_gif_playback.json
        └── final_delivery_receipt.json
```

The output, build, and assets folders are gitignored. Do not commit a review ZIP or the media-heavy PowerPoint. Only the latest two final PowerPoint revisions are retained. Older source/helper/asset trees and archived v12/v13 review PDFs remain as source evidence. Do not delete those trees as though they were duplicate presentations.

## Sources outside this build

- Submitted committee copy: `C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_Dissertation.pdf`, SHA-256 `447c2cc5a9dd219fa0857e298727fb68d4797af8f406251865eb22f4b126361f`.
- Advisor review: `C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_PhD_Defense_C3NS_20261030_v13_review_AJH.pdf`.
- Advisor’s original lecture: `C:\Users\Ben Bolen\Downloads\19_02_15 - INRC.pptx`. Cite Hunt’s February 15, 2019 INRC2 address, *Real-time control of dynamic walking in a dog-like robot*.
- New neural campaign: `D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\RESULTS_AND_FILES.md`, `COMMITTEE_QA.md`, and `figures\`.
- EB475WS4 drawings: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\functional_subnetwork_review\drawings\full_motor_loops`.
- Scientific generators, original MAT files, and measurements are mapped in `build/MECHANICS_SOURCES.md`, `assets/force/SOURCE_NOTES.md`, and the neural provenance JSON files.

## Scientific boundaries

The new human-model demonstration is the selected graded Hunt–PT ground9 RC controller. It is guided stepping in place: five right/four left qualifying clear, unloaded swings over 5–13 seconds, 48.848% BW net vertical assistance, and 115.561 N horizontal restraint RMS. It is not unsupported walking. The displayed GIF is saved-state replay, not new integration. The S5 quarter-strength contact comparison is separately identified. The native RF/HAM schematic and the BFLH motor-neuron trace use the selected controller’s actual edges and recorded channels.

The EB475WS4 full-motor-loop delivery instead depicts a historical w2lvar standing export. Gray dashed zero-conductance descending routes remain faithful to that source. The slides do not attribute that drawing to the new controller.

The 6.82× result is complete-process acceleration in a paired source benchmark; 18.31× covers neural updates. It is not a claim that the fully instrumented ground experiment runs faster than real time. The recorded numerical fields agree in that paired comparison.

No force or torque fits, route searches, controllers, physical models, paper, or dissertation sources were changed. Regenerated plots preserve saved values. Their rendering and provenance reports are retained. Older backup screenshots remain source screenshots; detailed interface labels are better examined in the accompanying source or offline demonstration.

## Rebuild

Run from this v21 directory, using the bundled Node and Python executables recorded in the scripts. Order: `run_build.mjs`, `package_current.py`, `finalize.mjs`, `render_review.ps1 -MakePdf`, and `verify_delivery.py`. PowerPoint rendering and final Windows file-link creation may require the normal unrestricted local process. Never kill an existing PowerPoint session. The scripts close only presentations they open.

Then inspect every changed render, check the review closure, and run the native movie and GIF checks. `finish_delivery.py` updates the pointer only after their checks pass. The finalizer writes an immutable receipt: archive the preceding final and `build/validation_delivery.json` before reusing their paths for a changed candidate. Intermediate PPTX copies and playback-test videos are regenerable and may be removed after delivery verification. Do not rerun `prepare_pipeline.py` without reconciling the later manual improvements to the generated scripts.

No commit or push was made. Ben controls commits through GitHub Desktop.
