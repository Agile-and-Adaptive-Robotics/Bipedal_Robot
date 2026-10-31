# Session 1 handoff — 31 October 2026

## Current priority and completed results

Ben's latest instruction is **newer W2L ground contact and training**, testing
whether contact input is routed to too many layers. Rat-to-Gait2392/CAD transfer
is paused. The selected canonical W2L now has **108 neurons/216 chemical edges**:
12 new homonymous II-to-MN relays cover flexors/extensors at both hips, knees,
and ankles, alongside the four retained phase-modulation II channels.

Actual108 fast-MuJoCo baseline and matched RG/PF/MN stimulus graphs are complete:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\CANONICAL108_STIMULUS_GRAPHS.md`.
The 12 s afferented fixed-root air baseline is stable (period approximately
1.072 s). All three -10 nA/100 ms trials have exactly matching pre-stimulus
trajectories. PF stimulation changes downstream output and joint motion while
largely preserving RG timing; measurements are in the figure index and raw data.
All 12 new homonymous synaptic currents were measured nonzero.

The actual108 ground baseline is finite with positive knee flexion, but produces
3 left/1 right qualified swings and 33.3% alternation, with 89.25% BW mean measured
assistance. This is distinct from older96 passing results below. Ground contact
class ablations and Rustuna calibration are the active next work. The contact
audit resolves 16 direct pathways, including a heel-to-ankle-MN maximum
conductance of 0.5 uS versus 0.01 uS for the other contact pathways.

Rustuna driver/protocol:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\RUSTUNA_CAPABILITY_RECORD.md`.
Rustuna runs in an isolated Python3.11 process, communicating with a persistent
myo Python3.10 simulation worker; the existing Optuna databases remain untouched.
No controller-training result is implied by the optimizer capability check.

Chapter 6 evidence/placement audit:
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\CH6_TRANSFER_EVIDENCE_GAPS.md`.
Interpret its proposed rat-transfer next step as paused by the newer instruction.

## Ben's decoded synapse names

Ben supplied this convention on 31 October: `{destination}_{Source}{Configuration}{type}{exc/inh}`.
Example: `RG_HF2i` denotes inhibitory type-II hip-flexor feedback to the rhythm
generator: RG = destination, H = hip source, F = flexor configuration, 2 = type II,
and i = inhibitory. Use this intended-path label alongside the actual saved
endpoints, reversal potential, maximum conductance, threshold, and saturation.
Do not silently rename original synapse types or infer physical attachment from
the name alone.

**Correction after Ben's model-selection clarification:** W2L means the newer
bilateral-RG controller with contact and type II feedback. Its canonical editable
project is `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer.aproj`.
The older 80-neuron results below are historical comparisons, not the selected
W2L paper result. The six generated GUI companion folders were **deleted** from
the active model directory after source snapshots and provenance were rehomed.
Their pre-cleanup contents are preserved in the single recovery ZIP
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\archives\w2l_generated_copies_before_cleanup_20261031.zip`.
Old generated-project paths below are historical and resolve inside that ZIP,
not active folders. Use the current experiment record for settings and cleanup map.

Importable public biarticular rat project:
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\Biarticular_net.aproj`.

Current W2L reference-drawing comparison:
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\subnetwork_review\LH_reference_comparison.png`.
Named type fields:
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\subnetwork_review\NAMED_SYNAPSE_TYPES.csv`.
Current work explicitly copies reference-local named type values as well as
node spacing, saved control points, and arrowheads; numerical changes have a
separate ledger and require fresh tests. Earlier parameter results are not
silently assigned to the updated canonical export.

The selected newer model has four actual II adapters: bilateral hip- and
knee-flexor-labelled receptor routes. It has bilateral RGs and heel/toe contacts.
The positive-V3 older-parameter benchmark passes both air checks and guided ground
(8/8 qualified swings; 15/15 alternating transitions), with 73.53% BW mean
assistance. The 35.2% BW result belongs specifically to the V3-zero diagnostic.
Keep those conditions distinct. Both have positive knee flexion throughout.

Rat import baseline: native AnimatLab completed 5.01 s. Left-leg RG/PF/MN and
hip/knee/ankle results are at
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\Native_baseline_neural_joint.png`.
Right RG latches in this public-source baseline. Mechanical MuJoCo import is
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\Deng_Biarticular_Rat.xml`.
It compiles with 16 actuators; coupled MuJoCo controller equivalence is a separate
next step. This is the public 2020 biarticular source, not a verified final-2023 release.

---
The earlier handoff below is preserved as history; the correction above takes
precedence for model selection and retired project locations.

**Machine: easteregg2.** Repository: `D:\GitHub\Bipedal_Robot`.
Python: `D:\Anaconda\envs\myo\python.exe`; MuJoCo 2.3.7.
Work is saved on disk. Usage ending does not remove these files. Raw run folders
are gitignored and local to this machine; include them when backing up experiments.
No commits, pushes, or Overleaf upload were made by Session 1.

## Start here

| Work | Full-path documentation |
|---|---|
| Virtual-walker experiments and current results | `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\EXPERIMENT_RECORD.md` |
| Virtual-walker run, plot, and GIF commands | `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\README.md` |
| Native AnimatLab guides and clock verification | `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Guided_20261031\EXPERIMENT_RECORD.md` |
| Dissertation changes and checks | `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\APPLIED_20261031.md` |
| BPA neural-coordination tests | `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\EXPERIMENT_RECORD.md` |
| Extended model-by-model instructions | `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\how_to_run_code_extended.md` |
| Cross-session summary for ZCode | `D:\GitHub\Bipedal_Robot\CHATGPT_REPORT.md` |

## Virtual walkers: current results

These are direct source-specific ports, using the supplied neural equations,
muscle paths, and Hill muscles. They are separate from Gait2392 and the CAD/BPA robot.
The original circuit has **80 neurons/150 chemical connections**; the newer GUI
branch has **96/204**. The older newer-model standalone has 208 edges and differs
from its GUI project; do not substitute it silently.

### Original: slower supported stepping

Matched 20 s trials, physical guides with 30 mm clearance, rails lowered 20 mm,
and a phase-triggered 5 N·m yaw perturbation:

| Neural clock scale | RG period | Qualified swings L/R | Alternating transitions |
|---|---:|---:|---:|
| 2 | 0.882 s | 21/22 | 100% |
| 4 | 1.760 s | 11/11 | 100% |
| 6 | 2.640 s | 8/11 | 77.8% |

Double support was approximately 27–28%. Mean **pelvis-only rotational energy**
decreased from 0.0941 to 0.0799 to 0.0722 J. Slower timing therefore reduced that
energy measure without guaranteeing better coordination. Clock 4 produced about
146/141 mm belt-relative strides, with measured net assistance 85.4% BW.
Knees stayed flexed. The guide is a passive finite carriage/collision apparatus;
the continuous sagittal harness is off in these conditions.

Comparison data:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\GUIDE_LONG_COMPARISON_V1.json`.
Each case saves its model XML, parameters, trajectory, and assessments.

### Newer: successful explicitly altered diagnostic

The supplied newer circuit latches. Symmetric zero RG-extensor tonic currents
**and** zero output gain on the two mirrored V3 → contralateral RG-extensor-IN
connections restore rhythmic coordination. Both changes condition this result;
V3 ablation alone does not restore all four RG cells under source tonic 2/0.
V3 neurons and c1 pathways remain present. Source files and default gains remain
preserved, with default-path trajectory equivalence checked.

Same saved diagnostic settings pass deafferented air (7/8 full knee cycles,
phase 0.488), afferented air (6/7, phase 0.496), and the bounded guided-ground
check (6/8 unloaded-clear swings, 11/13 alternating transitions = 84.6%).
Ground period is approximately 1.094 s; knees remain flexed. Measured **net**
vertical assistance is approximately 35.2% BW. Belt-relative progress is modest:
0.117 m over 10 s. Preserve the distinction between supported treadmill stepping
and forward travel through the room.

Support breakdown confirms this is almost entirely upward rail support, with
zero cap contacts across the saved trial. Full signed components are in
SUPPORT_BREAKDOWN.json in the diagnostic ground-result folder below.

Diagnostic configuration:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\NEWER_GUI1NA_TIED0_V3ZERO_DIAGNOSTIC.json`.
Ground result folder:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer_GUI1nA\newer_v3_diag_v1_ground_guides`.
This is a diagnostic variant, not a final optimized V3 architecture or a completed
Hunt/five-stage curriculum.

## Open these native AnimatLab models

Original guide, 10 mm:
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Guided_20261031\W2L_VirtualWalker_Guided.aproj`

Original guide, 30 mm:
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Guided30_20261031\W2L_VirtualWalker_Guided30.aproj`

Original guide, 10 mm, clock 2:
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_GuidedT2_20261031\W2L_VirtualWalker_GuidedT2.aproj`

Newer diagnostic, 30 mm, clock 2:
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_NewerV3Guided_20261031\W2L_VirtualWalker_NewerV3Guided.aproj`

Matching `.asim` and `_Standalone.asim` files are beside each project. All four
complete native 10 s diagnostics with finite states and no recorded knee
hyperextension. The newer native RG period is about 1.074 s, with approximately
half-cycle bilateral separation. Native foot-force/clearance qualification was
not inferred from MuJoCo. Native diagnostic verification:
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_NewerV3Guided_20261031\V3_NATIVE_VERIFICATION.json`.

Full newer native record and commands:
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_NewerV3Guided_20261031\EXPERIMENT_RECORD.md`.

Historical uniarticular subnetwork rebuild documentation:
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_20261031\SUBNETWORK_EXPERIMENT_RECORD.md`.
Page previews are in that model folder's `subnetwork_review` directory.
The guided projects are separate companions; existing unsaved GUI windows were
not touched. GUI opening itself was not verified by automation.

## GIF, still, and short results plot

Full paths, parameters, program versions, and caption for the original clock-4
example are in:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\GUIDED_WALKER_MEDIA.md`.

Animation:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\guides_long_v1_t4_c0.03_h0_p5\SUPPORTED_STEPPING_REPLAY.gif`

Still:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\guides_long_v1_t4_c0.03_h0_p5\SUPPORTED_STEPPING_STILL.png`

RG/PF/MN, joint angles, and ground-force plot:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\guides_long_v1_t4_c0.03_h0_p5\SHORT_NEURAL_MECHANICAL.png`

The newer diagnostic folder named above also contains a GIF, still, and short
neural/mechanical plot. Keep their controller and support captions distinct.

## Reproduce the MuJoCo examples

Use unique tags to preserve previous raw results. Commands run on easteregg2.

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py' --config 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\RUN_GUIDE_30MM_CLOCK2.json' --neural-time-scale 4 --duration 20 --tag my_original_clock4

& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py' --source 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_GUI1nA_20261031\W2L_VirtualWalker_Newer_GUI1nA_Standalone.asim' --config 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\NEWER_GUI1NA_TIED0_V3ZERO_DIAGNOSTIC.json' --guide-clearance .03 --duration 10 --tag my_newer_v3_diagnostic
```

Warm matched benchmark: 4 simulated seconds took **17.208 s scalar reference
versus 0.970 s Numba**, a 17.75× improvement, with physics and observers included.
Import, compilation, saving, and rendering are excluded. Contact-rich 20 s runs
took about 9.5–9.6 s. These are scoped measurements, not a 240× whole-model claim.

## Dissertation working files and review build

Working chapter:
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\60-sns.tex`

Working appendix:
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\95-AppendixD.tex`

Current Session 1 PDF:
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Archive\build_session1_sns_20261031_v5\main.pdf`

301 pages; zero LaTeX errors, warnings, or undefined citations/references.
Revised figure pages were rendered and checked. This is a local review build;
coordinate with the other editor before replacing delivery files or live sources.

Applied changes: brief AnimatLab results before MATLAB; Table 6.1 workflow;
short RG/PF/MN and kinematics displays; matched short spiking ground comparison;
archived superseded knee-startup D8; correct Hunt lecture citation; Deng 2023;
scoped speed measurements; positive simulated-BPA-control narrative and future
dynamic/hardware research. Chapter 20 and torque-paper sources were not edited.

Final check record:
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\FINAL_CHECKS_20261031.json`

Audit/comment reconciliation:
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\EB475WS4_RECONCILIATION_20261031.md`

Original bytes, requested PDFs retrieved from EB475WS4, and retired material:
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Archive\session1_sns_20261031`

## BPA tests, literature, and physics followups

BPA record above covers bounded one-joint neural coordination with pressure-state
delays, force-map errors, and loading changes. Feedback reduced tracking error in
the four tested conditions. The full CAD interface and force-sign checks are in:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\CAD_INTEGRATION_AUDIT.md`.

Public biarticular rat source and differences from Deng 2023:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\LITERATURE_MODEL_SOURCE_AUDIT.md`.
The located source is a 2020 precursor; final-2023 equivalence is not established.

Compliant-tendon/damping work remains separately documented under:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity`.
Contact initializer results:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\CONTACT_INITIALIZER_AND_ALIGNMENT.md`.
Do not substitute constant-activation holds for SNS curriculum results.

## Resume from here

1. Read the virtual-walker experiment record, then repeat the retained candidate
   in a fresh process and assess small perturbations before expanding training.
2. Reduce measured support in explicit small rail-height steps while retaining
   the same controller, knee guards, full-foot clearance, and alternation checks.
   Keep clock 2 and clock 4 results distinct.
3. Treat V3-zero as a causal diagnostic; investigate positive symmetric V3 gains
   separately rather than silently changing the supplied network.
4. Complete standing/Hunt gates before labeling the full curriculum complete.
   A biarticular derivative follows the retained uniarticular baseline.
5. Coordinate dissertation/deck integration with the other editor using the
   applied-change report and media manifest, preserving exact support captions.

Ben's saved Gait2392 continuation prompt remains at:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\NEXT_SESSION_PROMPT.md`.
That campaign is separate from these virtual walkers. Preserve the instruction
to prevent knee hyperextension, use speed-ups, and keep its training record separate.

## Stop status

Session 1's numerical and native AnimatLab jobs have finished; no background
simulation from these tasks remains running. The referenced documentation and
model/media paths were checked on disk before this handoff. Agent completion
records remain with their experiments. Resume explicitly from the steps above.
