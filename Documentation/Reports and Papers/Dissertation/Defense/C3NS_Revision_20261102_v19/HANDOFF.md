# Defense v19: delivery and continuation

Final deck: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\C3NS_Revision_20261102_v19\output\Bolen_PhD_Defense_RoseGreen_20261102_v19.pptx`

Final PDF: same folder, `Bolen_PhD_Defense_RoseGreen_20261102_v19.pdf`.

The presentation has 70 main slides and 53 hidden question slides. A 30-second average gives a 35-minute main talk, with 10 minutes available before Ben's 45-minute ceiling. This is a pacing plan, not a timed rehearsal. SADb detail remains at the end of the question slides, 117–123. Slides 12–48 emphasize the actuator and physical torque work; slides 49–64 follow the implemented models/tools in Table 6.1. Do not reorganize SADb into the neural section.

## What changed

The opening now establishes the research goals, mechanical benchmark, and four knowledge/tool gaps before presenting numerical results. Section transitions reuse the same robot and its actuator, knee, and neural-model insets. The mechanical sequence proceeds from existing force-model errors, through characterization, to knee-torque mismatch, correction, route redesign, and physical validation. The neural sequence begins with the actual controller-transfer failure, then presents models and scoped tool verification.

Removed the rejected “transmission model” and “transmission path” phrasing. No visible slide text contains “transmission.” Corrected normalization to epsilon620, with diameter-dependent coefficients. Corrected l0 and l620 to the inner boundaries of the active BPA section. Chi0 is a signed correction to current predicted BPA length; the positive flexor example points inward. The extensor's negative offset is explained in notes. Chi1/chi2 are effective system stiffnesses; the photograph does not invent coordinate-axis assignments. The bending-loss arc and equation remain.

Ben accepted the background-removed knee photograph shown in this chat. `assets\knee_cutout_approved.png` is used upright. Original photos remain in their original locations and in extracted source assets. The C3NS honeycomb, rose-green theme, small original branding, and original plot colors remain. The obsolete oversized title-page logo is covered without changing the honeycomb.

## Scientific references and boundaries

- Submitted PDF, unchanged: `C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_Dissertation.pdf`.
- Submitted PDF SHA256: `447c2cc5a9dd219fa0857e298727fb68d4797af8f406251865eb22f4b126361f`.
- Neither exact rejected phrase occurs in that 311-page PDF. Technical uses of “transmission” occur in the pulley/tackle discussion; `build\dissertation_wording_audit.json` records the actual contexts and PDF-page numbers. This is a bounded wording check, not a new full dissertation audit.
- Advisor comments: `C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_PhD_Defense_C3NS_20261030_v13_review_AJH.pdf`; extracted comments are in `Defense\AJH_Review_20261030\comments.json`.
- Latest handoff rechecked: `Defense\SESSION1_NEURAL_CONTROLLER_HANDOFF.md`, including its November 2 update. GitHub/local commit checked: `37b558f1`, “audits and defense.” The local update after that commit is also included in the source review.
- Neural audit sources: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\GAIT_DISSERTATION_RECIPE_REVIEW.md`, `DISSERTATION_TRIAL_PROVENANCE.md`, `DISSERTATION_FIRST_PASS_FIX.md`, `DISSERTATION_MODEL_MAP.md`, `ROUTE_BUILDER_EXPORT_CHECK.md`, and `DEFENSE_DESIGN_DECISIONS.md`.
- Historical full-s3k had active central feedback. The later 410-neuron standing export is a different condition. Six current-code static reconstructions are not byte-identical historical replays or new whole-body results. Optional repairs have not all been composed and dynamically validated. No full Hunt training or support-free walking is claimed.
- W2L80 suspended runtime measurements are explicitly distinguished from the newer W2L108 supported-ground demonstration. The original AnimatLab recordings predate the geometry correction. New design audits do not retroactively change the conditions of those recordings.
- Hunt attribution: A. J. Hunt, *Real-time control of dynamic walking in a dog-like robot*, INRC2 Meeting, Intel Labs, February 15, 2019. Ben's separate dissertation citation revision remains deferred.

No dissertation, Overleaf, simulation, controller, or experimental-result changes were made during this deck revision. No training was run. No commits or pushes were performed by the deck agents.

## Relevant file tree

```text
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\
├── Defense\
│   ├── current_deck.json
│   ├── SESSION1_NEURAL_CONTROLLER_HANDOFF.md
│   ├── AJH_Review_20261030\comments.json
│   ├── C3NS_Revision_20261031_v18\
│   │   ├── output\Bolen_PhD_Defense_RoseGreen_20261031_v18.pptx
│   │   └── submitted_figs\                     [source figures used by v19]
│   └── C3NS_Revision_20261102_v19\
│       ├── HANDOFF.md
│       ├── output\                            [final PPTX and PDF; Git-ignored]
│       ├── assets\                            [extracted figures and approved cutout]
│       ├── build\
│       │   ├── render\slide_001.png … slide_123.png
│       │   ├── supervisor_review.md
│       │   ├── mechanics_review.md
│       │   ├── neural_review.md
│       │   ├── validation_delivery.json
│       │   ├── delivery_checks.json
│       │   ├── native_video_playback.json
│       │   ├── slides_test.log
│       │   ├── final_slide_text.json
│       │   └── dissertation_wording_audit.json
│       ├── rebuild.mjs                         [main editable-slide authoring]
│       ├── mechanics_annotations.mjs
│       ├── supervisor_fixes.mjs
│       ├── neural_fixes.mjs
│       ├── package_deck.py                     [preserves native charts/media]
│       ├── finalize.mjs
│       ├── render_review.ps1
│       ├── check_video_playback.ps1
│       └── verify_delivery.py
├── Archive\session1_submission_20261031\
│   ├── Overleaf_source.zip
│   └── source\60-sns.tex, 95-AppendixD.tex
├── Figures\                                  [current source figures]
└── ProofFinal\                                [reference only in this revision]

D:\GitHub\Bipedal_Robot\Pictures\PictureAnalysis\
    20220308_031822.jpg                         [pinned flexor]
    20220307_222908.jpg                         [pinned extensor; 180° view]
    20230419_225929.jpg                         [biomimetic knee; CCW90° view]

D:\Users\Ben Bolen\NeuroCampaign_20261031\
├── review_20261101\                           [separate later research media]
└── fsa_design\                               [audit evidence and configuration records]

C:\Users\Ben Bolen\Documents\Github\ZCode_Skills\
└── assertion-evidence\SKILL.md                [updated shared design/verification guidance]

D:\GitHub\Bipedal_Robot_Ignored\               [older moved scratch/build trees]
```

## Verification and reviewer coverage

All 123 slides were rendered at 1600 × 900 in PowerPoint and individually reviewed. The supervisor inspected 1–27 and 71–88. The mechanics/controls/mathematics reviewer inspected 28–48 and 89–106, plus the corrected BPA-dimension slide. The neuroscience reviewer inspected 49–70 and 107–123. Corrected slides were re-rendered and individually rechecked. Their expertise lenses do not represent statements by the actual committee members.

Final package/import/geometry checks passed; the padded-canvas test reported no overflow. All seven native chart parts and seven embedded workbooks are byte-identical to v18. Embedded MP4 hashes are preserved. The two native AnimatLab videos on slide 52 were tested in PowerPoint: autoplay, advancement, repeat looping, and muted playback passed. GIF frame counts are retained in `delivery_checks.json`.

Final PPTX SHA256: `3495a9edd586b2c463b39d430dfac227b055b0de7eeb5dcfdc7750fd183cb0a9`.
Final PPTX: 94,941,742 bytes. PDF: 123 pages. Keep generated presentations, PDFs, renders, and intermediate packages ignored. Do not create another large review ZIP. Final PowerPoint revisions retained are v18 and v19; original reference talks and the advisor's PDF are not disposable revisions.

For another edit, use the current v19 source and the retained v18 package needed by its preservation step. Rebuild and re-render every changed slide. The disposable artifact draft/candidate and temporary test files are removed after delivery and regenerate from these sources. Ben commits source/skill changes himself.
