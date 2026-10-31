# SNS narrative and figure repair — Session 1

**Historical proposal.** The bounded implementation on October 31 is recorded in APPLIED_20261031.md. Statements below that the active manuscript is unchanged describe the earlier proposal stage, not the current working source.

Machine: easteregg2. Review proposal, 2026-10-30. No active chapter, bibliography, figure, or Overleaf content changed by this pass.

## Inputs and limits

Ben confirmed the ZIP at `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Bolen_Dissertation.zip`. Its Chapter 6 and Appendix D differ from the current working sources in ProofFinal. The ZIP has been inspected without extracting it over those sources. See INPUT_MANIFEST.json for hashes and snapshots.

Requested additional inputs:

- `D:\GitHub\Bipedal_Robot\output\pdf\Dissertation_Audit_PENDING_ONLY_v5_20261030.pdf`
- `D:\GitHub\Bipedal_Robot\output\pdf\Dissertation_Postdoc_Comment_Review_20261030.pdf`

Neither PDF was present at that path on easteregg2 during this pass. Their contents have NOT been reviewed. Do not substitute audit v4 or claim this proposal reconciles the missing documents. Ben has been asked for their accessible location/machine.

Read current Chapter 6, Appendix D, relevant CHATGPT_REPORT history, dissertation campaign requirements, and saved experimental records. Advisor comments include M105, “Results? Lessons?”, and M122, “This all still feels a bit awkward.” These support restoring the missing results-to-requirements transition. They do not authorize deleting unrelated material.

Live Overleaf was inspected read-only: 0 errors, 1 warning, 19 informational typesetting notices. The warning is the lone `\begin{table}[h]` in Appendix D, currently line 412 locally. Proposed correction: `[htbp]`. No numerical or scientific change is needed to address that warning.

## Proposed chapter order

Retain the section numbers and use Table 6.1 to order the methods/results narrative.

| Section | Proposed content | Evidence and purpose |
|---|---|---|
| 6.1 | Original AnimatLab model and controller | Mechanical model, shared RG, PF subsystems, MNs, and sensory connections. Clearly distinguish the original course model from later bilateral W2L revisions. |
| 6.2 | **AnimatLab Results and Lessons** | Short suspended neural/joint recording; brief separate unsupported ground-test result; explain why neural oscillation is insufficient for supported locomotion. |
| 6.3 | Knowledge gaps and requirements | Transfer controllers between models; efficient headless testing; GUI circuit editing; high-fidelity human mechanics; BPA robot mechanics; measured verification. Introduce Table 6.1 as the chapter roadmap. |
| 6.4, row 1 | Li–Deng/W2L implementation | Compact engine agreement and feedback examples. Relocate the short MATLAB implementation account here or to the graphical-tools subsection. |
| 6.4, row 2 | Gait2392-based models | Human and robot-route variants, conversion checks, muscle-force limits, existing IK comparison, and limits of supported motion. |
| 6.4, row 3 | SNS–Simulink/Simscape | Block editing, neuronal variants, verified reflex examples, and the advantage of connecting circuits to mechanics. |
| 6.4, row 4 | CAD-derived robot | CAD import and neural/BPA coupling. Identify the actual actuator law used in each demonstration. |
| 6.4, row 5 | MuJoCo–Simulink bridge | Activation/sensor exchange and tested closed-loop cases. Distinguish signal exchange from automatic translation of arbitrary block diagrams. |
| 6.4, row 6 | Curriculum and experiment drivers | Hunt's four operations and Ben's five physical tasks, harness/contact measurements, fixed-configuration checks, and scoped speed measurements. |
| 6.4, row 7 | Literature circuits and neuronal options | Demonstrate circuit adaptation, stimulation at different layers, and spiking/nonspiking interfaces; state which translations have actually run. |
| 6.5 | Discussion and transition | Explain what the tools enable, then finish with the CAD/BPA supported-standing example and its connection to the robot-design chapter. |

The separate draft file supplies replacement prose for 6.2 and 6.3 and proposed transition text. It is not a replacement for the whole chapter. Retain useful MATLAB implementation results and move them into the relevant tool demonstration.

## Figure decisions

| Figure | Proposed action | Data and safeguards |
|---|---|---|
| Original model/controller | Keep a readable mechanical rendering beside a hierarchical neural schematic. | Cite the actual original configuration; shared RG and later bilateral RG must not be conflated. |
| New AnimatLab result | Add a compact RG/PF/MN + hip/knee/ankle example, about three cycles. | Preview generated from the same saved graded recording, 1–2.5 s. Original AnimatLab joint signs retained; radians converted to degrees after checking `JointRotation` chart types. Representative left-leg channels, not a new bilateral or ground-walking result. |
| 6.3 | Replace four isolated HC panels with two leg panels, E/F paired, NumPy and Simulink overlaid. | Preview uses existing traces at 2–5 s. No time shift or amplitude adjustment. Voltages are relative to the network's resting reference. |
| 6.4 | Replace the dense 30 s display with a few common cycles of neural and mechanical response. | Retain full-record outcome metrics in a small results table. Match RG, PF, MN, and joint traces to one run and one time axis. If a channel was not recorded, identify that gap rather than infer it. Feedback-on/off comparisons must use the same physical interval or disclose their different states. |
| Existing Gait2392 IK overlay | Retain and improve readability if needed. | It already exists in the main chapter. Preserve absolute-angle error information in Appendix D; cycle normalization must not conceal offsets or unilateral motion. |
| D1 | Replace “motifs” with **neural subsystems**. | Caption prefixes: `\textbf{(A)} Rhythm generation`, `\textbf{(B)} Pattern formation`, and `\textbf{(C)} Sensorimotor pathways`. Preserve references when changing descriptive wording. |
| D6 | Shorten the neural-only graded/spiking display. | Preview uses 2–4.5 s and spike-event marks derived from positive cumulative-count increments. Original saved graded signal is E–F; it is labeled as a readout, not fabricated individual-cell voltages. No unrelated kinematic data appended. |
| D7 | Short mechanical comparison, with an independently identified benchmark overlay where valid. | Show hip, knee, and ankle in aligned panels. State differences in motion, tuning, support, and recording. An untuned neuronal replacement is not an equal-performance test. |
| D8 | Recommend retiring the old startup comparison from the active argument; retain source/results. | This five-second startup test differs from the newer settled-pulse experiment. It should not be called the same matched operating-point experiment. Ben's disposition is pending. No file has been moved or deleted. |
| Stimulation example | Use one clear RG/PF/MN intervention with pulse bar, baseline, and short joint response. | Existing layer-stimulus records contain same-run neural and joint arrays; select a response that supports the stated reset/non-reset distinction. |
| Reflex demonstrations | Keep a concise knee and/or loaded-arm demonstration, with complete feedback connections. | Preserve verified results and corrected six-feedback-path diagrams. Explain the experimental question before naming software. |
| Final robot example | Finish with neural/BPA supported standing if the matching recording is available. | Existing MuJoCo robot GIFs use ideal reference tensions plus an engineering allocator. They cannot be relabeled as calibrated BPA-force-law standing. The author-confirmed Simscape/BPA contribution may be described, with its recording/status identified honestly. |

### Review previews already generated

All in the adjacent `previews` folder, with sources and hashes in PREVIEW_PROVENANCE.json:

1. `animatlab_short_neural_mechanical_preview.png`
2. `fig6_3_compact_engine_overlay.png`
3. `figD6_short_neural_preview.png`

These are review images, not active chapter assets. No new simulation was run to obtain them. Every included panel was visually inspected. The compact engine figure's legend was moved outside the plots to avoid covering the curves.

## Evidence for the missing AnimatLab narrative

The saved graded suspended example is documented in `Code\MuJoCo_SNS\spinal\reports_20260924\animatlab_w2l_airwalk.md`. Its median onset period is 0.4432 s; the FFT estimate is 2.222 Hz. These are different estimators, so use approximately 2.2 Hz in prose rather than implying exact reciprocity. The E/F traces are graded activity, not action potentials.

The separate unsupported ground-test series is documented in `Code\MuJoCo_SNS\spinal\reports_20261020\w2l_animatlab\W2L_ANIMATLAB_TUNED_DROP_20261020.md`. In the stock 30 s example, pelvis height falls from approximately 0.99 m to 0.10 m by 1.5 s while the neural rhythm continues. This supports the distinction between an oscillating circuit and a controller that supports the body. It is not a matched platform-release continuation of the suspended recording, and the text must not imply that it is.

## Literature and citation corrections

- The supplied memory-function and animal-kinematics comparison examples match Deng et al. (2019), https://doi.org/10.3390/biomimetics4010021 . Use their presentation principles: short shared time axis, stimulus timing, neural layers, and mechanical output. Do not claim these figures are our results.
- Incorporate Deng et al. (2023), *Biarticular Muscles Improve the Stability of a Neuromechanical Model of the Rat Hindlimb*, https://doi.org/10.1007/978-3-031-39504-8_2 . It motivates treating biarticular mechanics explicitly. It does not by itself validate the present RF/hamstring four-cell implementation or establish that this implementation exactly reproduces that paper.
- Cite Hunt's four operations to the user-supplied 2019 INRC2 presentation: A. J. Hunt, *Real-time control of dynamic walking in a dog-like robot*, February 15, 2019. The current attribution to the 2017 article needs correction. Inspect the presentation for the exact four operations before finalizing their wording and mapping.
- Do not claim Python lacks neural GUI tools. Nengo provides a GUI (https://www.nengo.ai/getting-started/). State the narrower practical requirement: the chosen SNS workflow needed graphical circuit editing connected to the mechanical modeling tools used here.

## Production schematics: drawing contract

Create one consistent diagram style and separate configuration-specific exports for s3k, w2lvar, and syn6, plus a clearly labeled experimental biw2lvar-Hunt-PT diagram. Do not imply the experimental controller is the validated production baseline.

Each diagram needs: configuration identifier; neuronal model; RG/PF/MN hierarchy; commissural connections; Ia/II/Ib and separate heel/toe feedback; mechanical inputs/outputs; active actuator groups; and an expanded representative subsystem. Use separate symbols for excitation and inhibition and show any sign-changing/gated feedback explicitly. Distinguish grouped populations from individual cells. Populate edges and strengths from the saved run's actual built network and parameters. A zero-gain pathway must not look like an active synapse. Preserve separate sartorius and gracilis groups and identify whether the RF/hamstring four-cell PF subsystem is present in that specific run.

Keep engineering balance/force allocation visually distinct from biological reflex circuitry. No unexplained “Hunt bank” box. If auxiliary fitting features are present, name their signal-processing role and cite their implementation separately from Hunt's method. Do not draw literature-derived connections that the selected run did not actually contain.

## Speed, tendon, curriculum, and deployment claims

The saved backend benchmark at `Code\MuJoCo_SNS\spinal\hunt_training\runs\backend_benchmarks.json` includes an 8 s case with wall times 122.1196 s and 13.2528 s: approximately 9.21× end-to-end acceleration with identical saved numerical outputs. Its isolated neural calculation improves by approximately 26.6×. Keep those scopes separate.

Do not publish the suggested 240× full-simulation or 8× compliant-tendon figures without matched benchmark records. The compliant-tendon muscle-bank benchmark is an isolated component; dividing one component ratio by an unrelated speedup does not establish total simulation speed. State whether rendering, logging, compilation/warm-up, and physics are included. “Faster than real time” is allowed only for the measured computation that achieves it. Hardware export/readiness does not demonstrate real-time operation of the complete detailed model.

Separate the four training operations from the five physical tasks. Stage 4 requires ground contact, stance support, unloading, and swing clearance while the harness assistance is measured. Suspended oscillation does not satisfy it. Report which exact unchanged configuration passes which stages; do not splice stage successes from incompatible configurations. Preserve the no-knee-hyperextension requirement for new runs.

The latest experimental Gait2392 record has supported stepping at approximately 0.637 body-weight measured assistance, while stage 3 remains unresolved for that force-model configuration. Earlier successful standing under a different configuration cannot be used to close that gap. No harness-free walking claim is supported. This is background for selecting a demonstration, not authorization to restart training during the dissertation repair.

### Clarification requested by Ben: ideal reference tension and hybrid allocation

The robot run's engineering controller calculates joint torques from desired posture/motion and distributes these across the routed actuators using bounded least squares and their moment arms. It then supplies target motor-neuron currents. Crucially, `motor_target_current` cancels the instantaneous background synaptic current at those MNs to make them track the engineering targets. Thus, merely showing active PF and MN cells does not demonstrate that their ordinary synaptic coordination generated the forces. RG timing and central sensory interactions remain active, but this is a mechanics/routing diagnostic with SNS pacing and imposed muscle-force control.

The ideal actuator has approximately `tension = activation × force capacity`, with a 30 ms first-order activation filter and force caps. It is not instantaneous or unlimited, but it omits the BPA pressure/contraction force law and pneumatic dynamics. The 65% GIF also retains substantial vertical support and lateral/angular guidance. Appropriate short label: **SNS-paced robot motion with engineering muscle-force allocation and harness support.** Do not label it autonomous SNS control or a BPA-actuated standing result. The source and archived runner are identified in SIMULATION_GIF_EDITOR_HANDOFF.md.

## Completion and remaining decisions

### Ben's clarification: neural coordination with an approximate BPA model

An output conversion between a motor neuron and an actuator is expected and does not invalidate neural coordination. AnimatLab's polygon mapping, SNS–Simscape conversion blocks, and muscle excitation/activation interfaces serve this role. Millard et al. (2013), https://pmc.ncbi.nlm.nih.gov/articles/PMC3705831/, distinguish excitation, activation dynamics, and state-dependent contraction dynamics. A static conversion alone is not the same as those dynamic processes; neural output need not encode a literal force target.

Ben's intended architecture is neural coordination with a statically characterized BPA model that will be improved dynamically. Sensory feedback should let the SNS compensate for actuator/model errors. A suitable signal chain is SNS output → command conversion (possibly inverse BPA force map) → pressure dynamics → BPA force → robot mechanics → sensory feedback to SNS. The conversion must preserve the SNS's authority over timing and recruitment. The issue in the historical robot diagnostic is the independent torque allocator and MN-current cancellation, not the presence of a conversion block.

Robustness to dynamic-model error must be evaluated with errors actually represented: alternative pressure response, force-map parameters, loading/unloading behavior, and external disturbances, followed by held-out tests or hardware measurements. A controller cannot be shown to compensate for an effect absent from its evaluation model. Do not imply that better neural feedback requires a perfect actuator model, or that static BPA characterization is unsuitable as the initial forward model. No controller code was changed during this conceptual clarification.

Completed: source snapshots and provenance; active narrative/figure inspection; relevant advisor-comment history; saved-data previews; live error/warning check; exact future controller prompt saved in `Code\MuJoCo_SNS\spinal\hunt_training\NEXT_SESSION_PROMPT.md`.

Pending: accessible v5 audit and colleague-review PDFs; Ben's D8 disposition; selection/status of the ending BPA-standing recording; configuration-specific production schematics; the full short-window replacement for 6.4; final source reconciliation and manuscript edits. After edits, perform a fresh full build, visual checks, and a per-directive compliance table with source-line evidence. Do not claim a completed dissertation repair before these steps.
