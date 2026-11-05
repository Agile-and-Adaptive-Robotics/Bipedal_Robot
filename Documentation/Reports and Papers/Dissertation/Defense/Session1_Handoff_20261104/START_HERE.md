# Start here: Session 1 handoff for the defense-deck builder

Prepared on **easteregg2**, 4 November 2026. This package puts the previously external handoff documents and selected presentation assets inside the repository. It is ready for Ben to commit and push; creating this package did not commit or upload anything.

**Later calculation update:** [new Gait2392 air GIF, controller results, postural tests, and full experiment paths](experiments_20261104/EXPERIMENT_UPDATE.md). Ben subsequently authorized training. The new graded Hunt candidate passed both suspended-air stages; all four families passed intrinsic/timing blocks. This update supersedes the preparation-only status below. Raw calculations remain on easteregg2; selected media and reports travel in this package.

## Open this on EB475WS4

Repository: `D:\Github\Bipedal_Robot`.

Full handoff path:

```text
D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\START_HERE.md
```

Ben should commit/push this folder on easteregg2, then pull the same branch on EB475WS4: **KneeTestSetup_BenBo_stw**. The current repository branch is not master. All links in this index point to repository files and work after that pull. Original-source paths under `D:\Users\Ben Bolen\...` inside the copied reports refer to archives on easteregg2; do not assume they exist on EB475WS4. Use the media and figure links below for building the deck.

## Read in this order

Builder task: revise the neuromechanical-tools slides around the two worked workflows, then use selected neural/mechanical plots and supported-motion GIFs as evidence that the components operate. Explain the conversion repairs, graphical neural construction, actuator interfaces, accelerated testing, and staged calibration as a connected research toolchain. Preserve the deck's established style and source citations. Use the support and controller identities below; do not present the untested new postural pathway as a completed result. No new simulations or manuscript edits are needed to use this handoff.

For committee preparation, use the math/design reports to address Turcic's control-design and timing questions, Reddy's circuit/afferent interpretation questions, Jiang's objective/constraint and robustness questions, and Hunt's staged calibration and circuit-function questions. Distinguish a six-synergy basis from the conductance-based circuit implementing it, and distinguish the seven-step spiking FSA design procedure from the physical curriculum.

1. [Defense tools brief](docs/DEFENSE_TOOLS_BRIEF.md): contributions, implemented tools, and evidence.
2. [Gait repair, training, and deployment walkthrough](docs/GAIT_REPAIR_TRAINING_AND_DEPLOYMENT_20261103.md): scaled OpenSim conversion/repair, FSA, Hunt/curriculum, speed-ups, robustness, and hardware interfaces.
3. [GUI workflow handoff](docs/GUI_WORKFLOW_HANDOFF.md): the human-model and CAD/SNS-Simscape routes.
4. [Implementation clarifications](docs/CLARIFICATIONS_AND_IMPLEMENTATION_20261104.md): exact muscle families, parameter-test scope, drive naming, manuscript omissions, and the Rybak comparison.
5. [Script/function evidence inventory](docs/TOOLCHAIN_EVIDENCE_INVENTORY.md): concrete implementations and checks.
6. [All-GIF index](docs/GIFS_AND_TOOLS_INDEX.md) and [fresh walker result paths](records/OTHER_MODELS_RESULTS_AND_FILES.md): identities, support, and provenance. Their external links are archival; the portable copies are below.

## Current presentation media in this package

| Model/condition | Measured vertical assistance | GIF | Still |
|---|---:|---|---|
| Gait2392, graded Hunt-PT + autogenic RC | 48.85% BW, 5–13 s window | [GIF](media/Gait2392_Hunt_PT_RC_48p85pctBW.gif) | [PNG](media/Gait2392_Hunt_PT_RC_48p85pctBW.png) |
| Corrected CAD model, reference-tension actuators | 77.50% BW, 16–28 s window | [GIF](media/CAD_reference_tension_77p5pctBW.gif) | [PNG](media/CAD_reference_tension_77p5pctBW.png) |
| Graded W2L108 | 39.24% BW, whole 10 s trial | [GIF](media/W2L108_graded_39p24pctBW.gif) | [PNG](media/W2L108_graded_39p24pctBW.png) |
| Mixed spiking/graded W2L | 42.91% BW, whole 10 s trial | [GIF](media/W2L704_hybrid_42p91pctBW.gif) | [PNG](media/W2L704_hybrid_42p91pctBW.png) |

These are saved-state replays of the documented supported ground conditions, not runs of the newly prepared four-command or future postural-premotor revisions. Guides and attitude assistance remain. Gait and CAD step approximately in place. Support fractions describe their stated time windows, not matching conditions across every model. CAD's current local result is 5/8 perturbation cases; the other three stop at the predictive knee margin.

The CAD GIF uses 54 ideal reference tensions with SNS timing and an engineering allocator. Do not label it as a calibrated dynamic BPA/pressure-actuator experiment. The graded and spiking W2L GIFs have 12 uniarticular routes and differing retained mechanics/poses; they are not the separate biarticular rat or a matched comparison proving a benefit of spiking.

The older hopping GIF is already tracked elsewhere in the repo: [historical hopping GIF](<../../../../../Code/MuJoCo_SNS/slide_media/session1_20261030/Robot_supported_hopping_65pct.gif>). Its [original three-GIF handoff](<../../../../../Code/MuJoCo_SNS/slide_media/session1_20261030/SIMULATION_GIF_EDITOR_HANDOFF.md>) gives that condition's parameters and assistance. Keep its condition separate from the newer CAD clip above.

## Figures ready for slides

- [Human neural traces](figures/GAIT_ALIGNED_NEURAL.png), [PDF](figures/GAIT_ALIGNED_NEURAL.pdf).
- [Human mechanics/contact traces](figures/GAIT_ALIGNED_MECHANICS.png), [PDF](figures/GAIT_ALIGNED_MECHANICS.pdf).
- [CAD neural and mechanical traces](figures/CAD_NEURAL_AND_MECHANICS.png), [editable SVG](figures/CAD_NEURAL_AND_MECHANICS.svg), [PDF](figures/CAD_NEURAL_AND_MECHANICS.pdf).
- [Research toolchain workflow](figures/RESEARCH_TOOLCHAIN_WORKFLOW.png), [editable SVG](figures/RESEARCH_TOOLCHAIN_WORKFLOW.svg).
- [FSA/Hunt/curriculum workflow](figures/FSA_HUNT_CURRICULUM_WORKFLOW.png), [editable SVG](figures/FSA_HUNT_CURRICULUM_WORKFLOW.svg).
- [Paired Gait speed comparison](figures/PAIRED_GAIT_SPEEDUP.png), [editable SVG](figures/PAIRED_GAIT_SPEEDUP.svg).

The workflow diagrams include proposed orchestration as described in their source reports. The runtime bridge exchanges signals with a MuJoCo mechanism; it is not a universal Simscape-mechanism exporter. Keep kernel-only speed ratios distinct from complete simulator timing and hardware latency. [Acceleration record](records/NEURAL_ACCELERATION_20261104.md) gives the more recent incremental-kernel comparison as well.

## Dissertation source and committee preparation

- [Final revision source entry](<../../Final/main.tex>), [revision overview](<../../Final/README.md>), [build/review report](<../../Final/BUILD_REVIEW.md>).
- [Chapter 6](<../../Final/chapters/60-sns.tex>): the first rewrite follows the human-model and graphical/CAD workflows.
- [Chapter 7](<../../Final/chapters/70-overall.tex>) and [Appendix D](<../../Final/chapters/95-AppendixD.tex>): now explicitly describe the embedded-interface scaffold and unfinished SNS/hardware binding.
- [Controller math and evidence](<../../../../../Code/MuJoCo_SNS/spinal/fsa_design_20261101/DEFENSE_CONTROLLER_MATH_20261102.md>), [design decisions](<../../../../../Code/MuJoCo_SNS/spinal/fsa_design_20261101/DEFENSE_DESIGN_DECISIONS.md>). They describe the versions identified in those reports, not an already validated new postural controller.

Final is a separate review draft. ProofFinal and live Overleaf were not overwritten by Session 1. Source/support files are in the repo; the large compiled 310-page review PDF remains on easteregg2 at `D:\Users\Ben Bolen\ControllerPreparation_20261104\final_tex_build_v1\main.pdf`. Build locally from Final/main.tex if a PDF is needed on EB475WS4.

## Corrections and current research status: read before making claims

1. The historical human hand-drive scale reaching zero did **not** eliminate all descending drive. `wean_tune_20261016.py` retained shared DRIVE=2.0773198435110425 and scaled only supplementary left/right channels. Final Chapter 7 and Appendix D correct that interpretation.
2. The latest four-command Hunt/s3k/w2lvar/syn6 candidates passed static construction/wiring checks. They have not undergone new neural replay, retraining, or body trials. [Preparation record](records/CONTROLLER_PREPARATION_20261104.md).
3. Ben's latest instruction is to **test one postural-premotor candidate before deploying it into all four families**. The older report's pending postural replacement remains pending; do not present it as an implemented or trained new controller.
4. Contact policy clarification: flag **any explicitly inhibitory contact route**, even to one channel. Toe-specific feedback to RG is excluded from the proposed baseline; toe-to-ankle-PF feedback is a separately testable hypothesis. Distinguish explicit sensory wiring from normal reciprocal inhibition inside the CPG. The stronger audit/test revision is next work; the copied static report predates this clarification.
5. Training remains held while Ben's MATLAB parallel calculations run. This handoff copied documents/media only; it ran no training or simulation and did not modify those MATLAB jobs.

## Provenance and what to commit

Commit **this entire Session1_Handoff_20261104 folder** so its documents, media, and relative links travel together. There are no raw simulation NPZ files, optimizer databases, model mesh trees, or trial logs in the package. [Package manifest](PACKAGE_MANIFEST.json) records exact source paths, sizes, and SHA256 values; files were copied without re-encoding. Source files elsewhere in the repo remain linked rather than duplicated.

The copied reports retain their original dates and archival machine paths. This index's status/corrections take precedence over older statements such as “current,” “next,” or pending postural work inside those snapshots.
