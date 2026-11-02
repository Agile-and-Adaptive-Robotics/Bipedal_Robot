# Session 1 handoff — 31 October 2026

## Editor coordination: submitting manuscript stays unchanged

Ben supplied the latest 304-page PDF at
`C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_Dissertation.pdf`
and explicitly instructed Session 1 to answer the editor's coordination request,
leaving working manuscript sources unchanged. Reply:
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\EDITOR_SNS_COORDINATION_REPLY_20261031.md`.
It resolves the 1.938 s full-cycle versus approximately 0.969 s half-cycle question,
documents the old D4 ankle's gross limit violations, proposes an existing-data
canonical108 bilateral PF preview with visible hip limits, and gives the exact
source/run hashes and v5-versus-captured-remote diffs. Preview and verification
JSON are in the adjacent `editor_D4_verification` folder. No simulation was rerun
for this coordination analysis. The editor owns integration and compilation.
Native-style current subnetwork vectors and parameter tables are staged in
`current_W2L_results\native_schematic`; the main plate now preserves at least
10 pt lettering at nine-inch placement. Original Figure 6.1B still requires the
editor's own-page sideways placement. Ben prioritizes Gait2392 and CAD/BPA robot
integration over extended W2L or weak production-controller comparisons.

## Latest: lower hip effort after the Rustuna fit

Machine: **easteregg2**. Preferred lower-effort replay/configuration:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\GENTLE_HIP_0P6_FLEX_0P5.json`.
This post-fit adjustment uses bilateral hip PF-to-MN factors 0.6 and an additional
flexor-only factor 0.5: hip extensor input gains are 60%, and flexor gains 30%,
of the original trained setting. Other recruitment, contact, posture, and physics
settings are retained. In the 10 s trial, qualified swings remain 9/8 with 100%
alternation; mean net assistance is 40.50% BW. Measured hip tension RMS decreases
from about 403 to 277 N, and positive hip muscle work from 683 to 413 J.
Knee minima are 2.453/2.256 degrees. A fresh-process repeat is numerically identical.
Hip flexion-limit dwell remains 34.5%/39.9%; this is a lower-effort result, with
remaining angle plateaus and no uniform acceleration improvement.

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py' --config 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\GENTLE_HIP_0P6_FLEX_0P5.json' --duration 10 --tag my_fresh_lower_hip_w2l
```

Replay folder:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_hip0p6_flex0p5_fresh_v1`.
It contains `SUPPORTED_STEPPING_REPLAY.gif` (3.00–7.24 s),
`SUPPORTED_STEPPING_STILL.png` (4.84 s), and
`HIP_ANGLE_TORQUE_COMPARISON.png/.pdf/.svg`.
Settings, derived mechanics, and trial comparisons are recorded separately in
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\EXPERIMENT_RECORD.md`
and `CANONICAL108_GROUND_GRAPHS.md` in that directory.
The original 39.24% BW Rustuna winner and staged paper figures below retain their
original data and captions; the 40.50% result is an explicitly separate refinement.

## Current priority and completed results

**Earlier ground calibration:** the corrected physical initial pose and four
explicit ankle-PF target swaps produce 9/8 qualified swings with 100% alternation
and 42.24% BW mean net assistance over 10 s. A fresh-process repeat is byte-exact.
The untuned crossed reference also passes 20 s (18/17 swings); the finite-belt
audit records adequate loaded-contact margins.

Rustuna completed 24/24 accepted trials in 113.89 s. Its selected shared PF-to-MN
factor is 1.305863499088313, and the bilateral absolute toe-to-ankle-MN conductance
is 0.031839851437587975 uS. The selected 10 s run retains 9/8 swings and 100%
alternation while reducing mean net assistance to 39.24% BW. Fresh-process and
actual CLI replay checks are byte-identical. Training and rendering jobs have
finished. This supersedes the earlier within-day failed-ground status below.

**Original Rustuna configuration used in the staged paper package:**

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py' --config 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\RUSTUNA108_GROUND_WINNER.json' --duration 10 --tag my_fresh_trained_w2l
```

The JSON includes the actual 17-element physical initial pose, phase assignment,
both fitted gains, source filename, and source-hash guard. The `.aproj` remains
the editable reference; **its default home/routing is not the complete tested
runtime configuration**. Use the JSON/command above for this working MuJoCo run.
Effective executed graph, including the target swaps and fitted conductances:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\PRODUCTION108_EFFECTIVE_CONNECTOME.json`.

**Paper/editor pickup:**
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results\build_review\REVIEW.pdf`.
Adjacent INSERT.tex, FIGURES.tex, RESULTS_TABLE.tex, figures, numerical results,
and provenance form the review package. This is staged insertion material;
the active dissertation and Overleaf were not changed by this result-packaging pass.

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


## Session 1 — final live submission revision, 2026-10-31

Ben authorized direct live Overleaf edits. Final PDF: 310 pages; Errors 0 / Warnings 0. Figure 6.1 pairs original walker and fall; original controller is standalone sideways Figure 6.2. Added W2L108 contact/type-II and rat subnetwork plates, supported treadmill still, human geometry, headless/GUI explanation, scoped speed table, and actuator fidelity paragraph. Detailed engine and IK plots moved to Appendix D. Current complete source ZIP and PDF preserved in Dissertation/Archive/session1_submission_20261031. Full change record: D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\submission_live_20261031\SESSION1_FINAL_HANDOFF.md. No new training in this final publication pass.


## Session 1 — urgent chapter opening and final knee orientation correction

Live Ch. 6 now has 6.1 compact AnimatLab motivation, 6.2 compact gap/framing, and 6.3 Gait2392. Bold Verification labels removed; Background circuit figures cross-referenced; clean vector SNS library diagram replaces distorted icons/tiny wiring; curriculum rows spaced. Neural-only 100.17x/56.55x and Python-to-Numba 584.33x timings added with explicit scopes. Final knee panel C: flexor LEFT, extensor RIGHT. Live compile Errors 0 / Warnings 0. Use live download for final PDF; last browser download attempt timed out. Full handoff: D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\structure_revision_20261031\SESSION1_STRUCTURE_HANDOFF.md.

## Hourly continuation checkpoint — 2026-10-31

Four workstreams remain active or have been resumed with bounded next tasks. Verified external syn6 biological standing checkpoint is fresh-process bit-exact and passes 8/8 corners of the three fitted parameters at ±1%; posture remains laterally/rotationally constrained and is not the earlier Hunt controller. Source: D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\STANDING_CHECKPOINT.json and standing_basin\BASIN_SUMMARY.json. Sensory-spiking W2L and Gait stepping checkpoints pass with graded central RG/PF/MN retained; central spiking feasibility is still failing its neural gate. Full-Gait compliant saved-excitation replay is validated and runs 2.05x real time, but live SNS integration remains next. Editable W2L II subnetwork parity is verified; rat adapter audit found effective ampere-scale injection and the agent is deriving isolated corrected mappings. No manuscript changes. Raw outputs external. Process CIM inspection unavailable; Get-Process fallback used. No duplicate campaign launched by root.

## Hourly continuation — source-consistent live mechanics and actuator uncertainty

Verified SOURCE_DYNAMIC_CHECKPOINT.json: exact same constrained syn6 controller passes native, source-rigid, and source-compliant standing; compliant fresh repeat exact. Knees protected, retained rig. Hunt clearance fails both source-rigid and compliant while native passes; next calibration must be force-law-consistent, not blame compliance alone. Verified BPA_HELDOUT_FRESH_CHECKS.json: live-feedback virtual-BPA standing passes held-out20% forward-map loss and16ms pressure delay with exact fresh traces; corresponding nominal sensory replay fails characterized equilibrium domain. Geometry/dynamics remain declared assumptions, local inverse remains active, and composite stressor still fails. Resumed physics motor calibration and rat source-circuit work; spiking and BPA agents continuing. Central spiking W2L ground success reported by agent, compact record requested before public promotion. No manuscripts changed; no duplicate root simulations launched. Get-Process inspection returned no visible matching processes; agent status/checkpoint inventory used for coordination.

## Hourly continuation — calibrated compliant stepping and central-spiking W2L

Disk-verified HUNT_CALIBRATION_CHECKPOINT: same mirrored1.5/1.5/1.5 functional recruitment scales pass both air stages and source-compliant supported ground;5/3 sustained cleared swings,68.8%BW assistance, exact fresh repeat,8/8±1% gain corners across air/ground. Stage3 balance remains unverified for this calibrated controller; next work assigned, no completed-curriculum claim. Disk-verified CENTRAL_RESET_REST_CHECKPOINT: central-spiking W2L10/10 clear swings,46.02%BW net support,192 actual spiking cells with96 graded remaining,16/16 perturbations per record; lower rails increase support, not weaning. Prepared virtual-BPA standing passes20s with rotational scale.8; reset-pressure controls identify startup reserve dependence. Physics, prepared hybrid/fair family comparisons, rat distal-PF causal diagnosis and combined spiking work resumed. All sources/manuscripts remain isolated from experiments; external durable records remain authoritative.

## Hourly continuation — combined spiking and remaining biomechanical limits

Disk-verified COMBINED_SPIKING_CHECKPOINT:704 actual RG/PF+sensory spiking cells,80 remaining graded source cells; nominal air/ground pass and complete-state/fresh event parity retained. Latest agent load probes show simply increasing ankle gain increases support, and rail release fails posterior pitch; mechanically valid CoM-home initialization is being tested separately. Physics agent reports a true50%BW source-compliant stepping clearance passer with exact repeat/1ms check, but severe ankle DF outside saved IK envelopes; not promoted as normal gait, biomechanical constraints requested before further fits. Prepared biological/BPA sensory-spiking fresh comparisons saved and release-state audited; s3k standing passes, w2lvar still fails broad three-gain search and now needs current/group diagnosis. Rat converter diagnostic restores PF rhythms but right motor output wiring remains deficient; source-supported experimental ipsilateral reconstruction assigned. All four workstreams continue, originals/manuscripts unchanged. No duplicate root simulations launched.

## Hourly continuation — all-path-preserving w2lvar standing mapping

Verified NORMALIZED_STANDING_CHECKPOINT: w2lvar standing variant normalizes each IIIN→MN group by actual antagonist fan-in, keeping all those connections nonzero; PF/aggregate-Ib walking drive remains disabled in this standing mode. Fresh repeat exact,20s and8/8 small-parameter checks pass; retained lateral/rotational guidance,4.17mm visual penetration,97.70%BW feet. Source-compliant gait still fails the predeclared subject-specific IK envelope in8/8 Ib trials, despite clearance; next assignment is source-consistent motor/afferent calibration, not further support reduction. W2L CoM-valid pooled spiking checkpoint is stable but rail-free balance fails; next phase-dependent neural body-balance pathway explicitly experimental. Rat ankle defects are real supplied-source settings, so the next comparison must separate enabled-state/home/force-length corrections from neural fits. Four bounded workstreams resumed; no manuscript edits or duplicate root physics runs.

## Hourly continuation — longer-hold failures and bounded next diagnostics

Verified normalized w2lvar virtual-BPA inverse-feasibility audit:8s frozen prepared transfer repeats exactly, but both20s runs later fail declared pect_r pressure/sleeve domain. Nonzero requested force can be infeasible; current inverse flag3 zero sentinel is consumed as a pressure command and vents. Assigned isolated explicit feasibility-boundary policies with unchanged forward guards, not geometry relaxation. Hunt ankle-current reconstruction found direct ankle drive active via circuit_options.json, correcting the initial missing-drive suspicion; local Ia ablation has tiny effect and amplitude/transient-only candidates still fail anatomical gates. Next loaded hip/knee feasibility and Hunt C motor coordination, not more blind ankle fits. Rat source-derived biomechanics comparisons all fail and need constraint/current audit. W2L pitch/capture relays not promoted; landing/timing analysis continues. Raw external artifacts retained, no manuscript edits or new root simulations. No user action required at this checkpoint.

## Hourly continuation — contact-feasible teacher and predictive reserve work

Reviewed BALANCE_CHECKPOINT, source-compliant causal audit and inverse-policy comparison; no duplicate root simulations. W2L0.2m/s belt condition retains frozen704/80 neural controller and improves measured supported load sharing, with fresh exact/8-seed robustness; this is an environmental condition change, not neural training/rail-free gait. Original Gait fixed-root IK target is not contact/no-slip feasible, so contact-consistent supported-step teacher construction assigned before further Hunt C fitting. Static nearest-feasible pressure policy does not resolve delayed actuator domain failure; explicit anticipatory neural reserve experiment assigned with strict forward guards and original failures preserved. Rat true-marker constraints supersede misleading SDK accumulated angles; efficient full-input PF diagnosis continues. No manuscript/Overleaf edits; experiment outputs remain external.

## Hourly continuation — no new promoted controller

Reviewed failed reserve→MN nominalfit and W2L outputgaincheckpoint. The nominal lower-support W2L knee2.3 candidate fails one heldoutgroundseed while matchedbaselinepasses; robust0.2baseline stayspreferred. Pressure reserve relay computations/parity are verified but allsixnominal20sfitsfail; separate explicit valve-reference relay assigned. Rat six source2019-supported missingRH PF→Ib gates now authorized as isolated132/270 reconstruction, not a gain tweak. Contact-consistent Gaitteacher passes positional/no-slip checks but source inverse rejects some dynamic/contact states; actual MTP capacity makes initial50/50heel/toe support infeasible, so capacity-constrainedheel-dominant teacher and separate contactgeometryaudit continue before any neural fit. No source anatomy, forwardguard, manuscript, or Overleaf changes. Root launched no duplicate numericaljobs; activeworkdelegated with existingcheckpoints.

## Hourly continuation — nominal valve-reserve success, no robustness promotion

Disk-verified pressure_reserve_valve/NOMINAL_CHECKPOINT: explicit neural pressure-reference relay permits20s nominal/pulse supportedstanding with exactfreshrepeat, but withheld delay/map/combined tests fail; globalcheckpoint unchanged. New preregistered uncertaintytraining assigned, nopostheldoutreselection. Mechanics teacher92-source single-supportholds at60% support have actualunloading/clearance and exact/timesteprepeat, but prunedproductionvariant fails; no neuralcurriculumclaim. Smooth native-contacttransition teacher assigned before HuntC. W2L terminationhypothesis was contradictedbymeasured timing, so nofurther gainsearch; symmetricoffline stanceprimitive target now assigned. Rat132/270 teachingbilateral timing improves but native/mechanicalgatesremainnegative; higher-rateclampedinputverification and boundedexpected-afferentmethodcontinue. No manuscriptchanges or rootduplicates; externalstorageaudit initiated.

## Hourly continuation — expected-input rat timing verified, mechanical teacher not complete

Verified strictHunt270 native checkpoint: mirrored7paired sensory groups, unpairedsourcepaths frozen, all132native neuralcadence/withinleg/bilateral gates pass with~4.86e−5mV replayRMS. This is prescribedengineeringaff input, not endogenousmechanicalsuccess; genuine HuntC waveform stage assigned with2019 homologues and explicitbiarticularlimits. Verified Gaitleft lift-return source92teacherpasses no-slip/contact/sourceROM withfresh/substep evidence, butrightreturnfails so bilateralcycle remainsunaccepted; nextroot-aware worldfootteacher, no neuralfit. BPAuncertainty8trialstudyhas0alltrainingpassers, validationunused; internalpressure-capacityweight variant assigned ratherthanraisingirrelevantoutputgain. W2L preservedrobustbaseline stillpreferred aftermechanicalfeasibilitydiagnostics. No manuscripts or rootduplicatejobs; compactnewoutputs policyretained.

## Continuation checkpoint — 2026-11-01 07:03 UTC

Machine: easteregg2. No duplicate simulations or manuscript edits.

- BPA observer comparison finished: pressure observer nominal 4/4, moderate uncertainty 3/4, strong uncertainty 0/4; measured-tension observer 4/4, 0/4, 0/4. No training finalist, and independent validation remains unused. Full neural and estimator state restoration and default physical parity pass. Next bounded experiment explicitly increases shared neural valve correction authority from 5 kPa to a registered 6–10 kPa range, retaining physical pressure limits and effort gates. Records: external gait_bpa/pressure_reserve_capacity/EXPERIMENT_RECORD.md and PROGRESS_SUMMARY.json.
- Rat Hunt C six-function waveform fit improved drive RMS 0.2245 to 0.1594, but VA/TA gates fail. Native neural parity is excellent while actual left ankle chatters near 483 Hz and exceeds its stop by 6.87 degrees. Neither neural parity nor smooth RG plots establish mechanical success. Next task isolates mechanical chatter before further fitting. Records: external subnetworks/Rat132_bounded_diagnostic/diagnostics/RH_Ib_gate_diagnosis/huntC_waveform/HUNTC_WAVEFORM_RECORD.md.
- W2L actual source-muscle/contact landing witnesses survive short bilateral mechanical tests at about 9–10% BW upward support. These are isolated 0.5-second stance witnesses, not reduced-support neural walking. Neural stance prototype retained alternating swings but worsened support to 47.47% BW, so the previous 36.48% baseline remains preferred. Next task examines an observable pre-contact hip-length pathway.
- Gait2392 teacher audit found that earlier two-second endpoint holds drift and are not equilibria. Full native follower/contact/source-muscle equilibrium is now the prerequisite for new continuous teachers. Protective knee-stop compensation was separately identified and preserved as a failed diagnostic; physical stops and all acceptance gates remain enabled.

All four agents have bounded next tasks. Exploration uses compact records; raw outputs remain in owned external campaign folders. No new closed-loop checkpoint promoted.

## Continuation checkpoint — 2026-11-01 08:03 UTC

Verified a full-native Gait2392 double-support endpoint with all 105 residuals retained, constant source controls, and five-second stationarity. Fresh process saved fields are bit-identical; 1 ms refinement also passes. Maximum joint drift is 0.01595 degrees and final static residual is 0.02727 N. This is original-source92 mechanics with 60% prescribed support, not an SNS stepping result. Left single support has also passed its initial five-second test with the opposite foot 13.319 mm clear and zero force; independent repeats are pending. Right endpoint remains required before a continuous bilateral teacher.

Rat chatter is now isolated to mechanical integration: original 2019 mechanics are quiet, reconstructed passive muscles can chatter with every command off, and identical active neural outputs yield 32.56 degrees left-ankle RMS disagreement between step sizes. Installed AnimatLab recognizes PhysicSubsteps, not the supplied PhysicsSubsteps tag; corrected substeps alone do not cure the problem. Further PF fitting is gated. Agent next audits the existing MuJoCo LinearHill import and a source-consistent stable tension integrator, with explicit model provenance.

BPA valve-authority study passed nominal and moderate fixtures in all four trials but none passed strong uncertainty. At the last failure, the neural target is 17.39 N while the assumed sleeve geometry requires at least 67.94 N. No validation or controls were run without a training passer. Next work audits posture, moment arms, and state-dependent neural recruitment rather than another fixed pressure-cap sweep. W2L interlock variants retain alternating motion but offer no meaningful support improvement; a load-gated additional spiking PF primitive is being tested. Prior preferred checkpoints remain unchanged.

All records and raw results remain in the external owned campaign folders; no manuscript edits or duplicate root simulations.

## Continuation checkpoint — 2026-11-01 09:07 UTC

All three original-source92 Gait2392 stationary endpoints now pass five-second, fresh-process, and 1 ms checks with unchanged physical parameters. Right endpoint required native acceleration-based solver conditioning; its final residual is 0.000198 N. Continuous lift/return sequences clear both swing feet and recover double support, but touchdown slip still fails (37.82/17.89 mm/s versus 5 mm/s). Next correction preserves calibrated loaded follower FK in Cartesian targeting. Hunt C remains gated.

The reconstructed rat now has an accepted analytic LinearHill mechanical alternative in MuJoCo: native force-equation replay RMS about 0.000014 N, 250-to-125 microsecond joint convergence RMS 0.02423 degrees, stops within 0.1645 degrees, and exact fresh repeats. Coupled mechanics is 3.1–3.2 times real time; 331 times applies only to compiled versus scalar force calculations. Muscle coefficients and source inertia are preserved; stop response is explicitly calibrated to a native torque fixture. Implicit velocity prototype remains rejected. Next task is actual receptor/strict-B neural coupling, not claiming it already works. Full-path instructions are in external subnetworks/Rat132_bounded_diagnostic/diagnostics/RH_Ib_gate_diagnosis/mujoco_muscle_fidelity/HOW_TO_RUN.md.

BPA measured-deficit-to-MN architecture has a frozen 5 nA paired candidate passing all three training fixtures and fresh 30-second nominal repeats. Independent validation is running, with no reselection; this is provisional guided standing. W2L load-conditioned PF mechanisms operate as intended but do not materially improve baseline support/clearance. Whole-cycle mechanically feasible landing references are the next task, not another strength sweep.

All agents have bounded continuation tasks. No manuscript edits, no new root simulation, and no production checkpoint replacement.

## Continuation checkpoint — 2026-11-01 10:06 UTC

Verified BPA neural compensation milestone: the frozen 5 nA paired pectineus MN feedback passes all three known guided-standing training conditions, two exact fresh nominal 30-second repeats, and an exact fresh strong 20-second repeat. Its matched added-feedback-off control fails the strong case at 17.294 seconds. Independent validation is 3/4, with the weaker 0.78 force-map failure preserved at 18.072 seconds and no reselection. Strong-case upward support averages 1.137% BW; lateral/rotational guidance remains. Geometry is assumed virtual BPA geometry. Full record and reproduction: external gait_bpa/source_domain_recruitment/EXPERIMENT_RECORD.md and HOW_TO_REPEAT.md. Next bounded task is a spiking replacement of only the two added relays, not a claim of an entirely spiking controller.

Rat strict-B now couples in MuJoCo through all 44 actual adapters and 32 muscle/receptor states. Mechanical air gates pass and neural expected-input replay is exact, but right/distal rhythm latches under actual afferents. No supported drop was attempted. Next task separates Ib/II and central PF competition and fits only real homologous paired connections, preserving unpaired paths and the original strict-B checkpoint.

Gait2392 predictive actuator audit reproduces native activation to 1.11e−16 and forces to about 2.7e−10 N on both saved sides. A causal 10 ms forecast predicts torque substantially better than the old instantaneous static command model. Live prediction remains experimental; warmup command blending was identified before any promotion. W2L full-cycle source dynamics reconstruction passes numerical checks; periodic closure/landing references remain under feasibility analysis.

All outputs remain external and compact where possible; independent evidence archives are preserved. No manuscript edits or global checkpoint changes.

## Continuation checkpoint — 2026-11-01 11:09 UTC

Spiking BPA error-relay mirror passes the nominal eight-second physical run with exact fresh repetition. The strong run strictly fails at 19.374 seconds, with the sleeve-domain guard unchanged; the graded successful result remains preferred. Only the two added relays are spiking. CAD mapping confirms neither pectineus route exists among the CAD54 routes, so direct reuse is not claimed. Next CAD task starts with explicit paired-route pressure mapping and feasibility, retaining mixed-actuation provenance if only the verified vas_med pair is usable.

Gait2392 loaded-FK geometry/reference repairs preserve legacy/default-off replay exactly, but no continuous teacher passes. Full q/velocity/activation/fiber/filter replay and directional dynamics demonstrate why accurate force-only prediction is insufficient: command-to-motion Jacobians differ by 153–264 percent. The next bounded inverse method propagates the complete augmented state; no horizon/strength sweep or Hunt C fitting is authorized before mechanical acceptance.

Rat actual-afferent replay search restores timing on the teaching tape but fails amplitude and real bilateral-input gates. Deafferented startup requires a declared symmetry-breaking pulse; sustained autonomous period is 0.5172 seconds versus the 0.292-second target. Next task explicitly starts at Hunt A with shared parameter fitting, rather than treating expected-input timing as complete training.

W2L mirrored motor replay has a physically valid whole cycle with 5/5 alternating cleared swings, but full-state periodic closure remains false and support is worse than the live baseline. Automatic approval review rejected the premature phase-only Rustuna search; no optimizer ran and no retry or workaround is authorized. Unaffected saved-state closure diagnosis may continue. The review reason and blocked search are retained for Ben's return.

No manuscript edits, no accepted checkpoint replacement, and no root duplicate simulations. All four agents have bounded continuation tasks; at most three heavy jobs.

## Continuation checkpoint — 2026-11-01 12:06 UTC

Verified CAD mixed pressure interface: two modeled vas_med BPAs carry about 221/223 N mean tension, peak 239/253 N, and peak pressure 499 kPa in a two-second test. Two fresh runs reproduce the saved physical/neural/pressure arrays exactly. Other 52 routes retain ideal reference tension; the existing 262-cell hybrid controller still imposes allocator targets through MN current cancellation. Harness support is 91.7% BW, and both knees stay flexed. The zero-coactivation pressure case produced zero tension and is preserved separately. Next experiment isolates direct SNS recruitment of the pressure pair without its online MN cancellation; no claim of all-BPA or autonomous balance.

Rat Hunt A now passes all declared timing, amplitude, phase, and stability gates with exact fresh parity. Its autonomous period is 0.6464 seconds, consistent with the preregistered Hunt target near twice the 0.292-second final cadence. Stage-A deafferentation explicitly masks outgoing afferent edges, because zero input alone leaves prototype resting currents active. Restored B baseline still fails. A new staged 12+12 trial B study is authorized; A remains frozen.

Gait2392 full augmented-state derivative audit identifies coupled branch changes at equilibrium and inadequate local ten-millisecond correction authority; no live controller or Hunt C was promoted. Next bounded nonlinear solve uses full state propagation and shared functional corrections with physically separate residual scales. No horizon sweep. W2L fixed-wave full-state closure still fails; rejected phase optimization remains disabled. Its agent is now preparing source-derived editable functional subnetworks and aligned baseline figures, an unaffected priority.

No manuscript edits, no original source/production physics changes, no duplicate root simulations. Raw results stay in external owned folders; separate records remain current.

## Continuation checkpoint — 2026-11-01 13:07 UTC

Gait2392 small load-transfer primitive now passes full nonlinear source-state dynamics: 1% BW redistribution over the source-response-derived 302 ms horizon, using eleven shared muscle-family corrections. Fresh final states are exact, 1 ms refinement passes without refitting, all-window slip stays below 0.44 mm/s, and knees remain flexed. This is a double-support primitive, not swing or a trained neural teacher. Next bounded continuation is 3% then 5% redistribution in both directions with continuous states and unchanged gates.

W2L editable functional package is verified in Neuromechanical_Models/Mujoco_SNS_models/W2L_704_80_baseline/README.md. Exact census is 108 logical channels, 216 chemical pathways, and 704 spiking plus 80 graded membranes; RG and PF populations are spiking. Nested editor JSON preserves all edge fields; vectors and a four-cycle live-baseline figure are rendered and inspected. No user browser state changed. The next figure task uses the actual saved Gait2392 family configurations rather than substituting generic build diagrams.

CAD PF-only unclamped recruitment is infeasible: additive RG drive through IBEXC dominates the supposed stance gate. More than 99.8% of its positive input is RG feedforward in this record. Selected-pair calibration is now explicitly extended to existing PF plus IBEXC output weights, tied bilaterally, with no MN background cancellation on the pressure pair and all original models preserved.

Rat staged B fails its dual-signal RG gate, so conditional PF trials and body/drop tests were not run. Next task reviews whether indefinite failed-load replay is an appropriate criterion versus finite perturbation recovery, using primary literature before registering any changed experiment. Old failures and gates remain intact.

No manuscript edits or production checkpoint promotion. Four bounded tasks assigned; no root duplicate simulations.

## Continuation checkpoint — 2026-11-01 14:12 UTC

CAD modeled BPA pair now uses literal SNS motor-neuron dynamics with zero applied clamp currents and no imposed 5% pair floor. Six existing PF/IBEXC conductances are normalized identically across sides. Matched two-second reference and pressure trials pass with exact fresh repeats, pressure peak 422 kPa, tension peaks 250/262 N, and protected knees. Other 52 routes remain hybrid-controlled; upward harness load is 91.8% BW. Next test extends this frozen pair to eight seconds and explicitly checks stationarity before any harness reduction. Record: external cad_bpa/pair_neural_recruitment/normalized_unclamped/EXPERIMENT_RECORD.md.

Gait2392 continuous 1%→3% load transfer passes in both directions, with exact fresh replay and 1 ms verification. Five-percent transfer fails near the diagnostic 0.2 mm lateral heel-body hold, while actual slip remains small. A separately preregistered rolling-contact task will distinguish physiological no-slip from that artificial body-origin hold; original failures will not be relabeled.

Editable Gait2392 family figures are ready at Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families/README.md. They match actual saved standing candidate neural steps exactly, show active versus zero-conductance paths, and preserve full synapse fields. All three selected candidates are graded; no spiking claim is inferred from these diagrams. Next task is actual short spiking Ia/Ib reflex demonstrations with measured downstream force and joint responses.

Rat teaching audit found asymmetric clipping and insufficient excitation in three fitted channels. Further searches on the failed-body-derived tape stop. Next tapes must propagate each side's own muscle/receptor states from validated kinematics and explicitly declared commands. Source coefficients and accepted Hunt A stay preserved.

No manuscript edits or production changes; bounded tasks remain active with no duplicate root jobs.

## Continuation checkpoint — 2026-11-01 15:15 UTC

CAD pressure pair passes the frozen eight-second comparison and steady-tail gates; ideal and pressure cases both repeat exactly. Literal pair SNS currents remain unclamped, pressure peak 442.22 kPa, tail tensions 201.98/202.26 N, and knees remain flexed. Harness support is still 91.87% BW. Startup shifts COM posterior by approximately 5–6 cm at settlement, shown explicitly in the context figure. Next single experiment modestly reduces measured vertical assistance while retaining all neural/pressure settings and lateral/rotational guides.

Gait2392 rolling-stance 1→3→5% BW transfers pass both directions, fresh replay, and 1 ms gate checks; actual slip stays below 1.17 mm/s and penetration below 1.17 mm. Right MTP acceleration is timestep-sensitive around 0.806 seconds, so no claim of pointwise dynamical convergence. Next bounded continuation is 7.5% then 10%, preserving continuous states and physical constraints. No swing or Hunt C promotion yet.

Rat source-consistent teaching motion now maps at 0.013 mm RMS, but earlier experimental receptor conversions flatten three channels. Original 2019 telemetry verifies real modulated Ia input and identifies receptor/adapter-domain mismatch, not missing biology or a unit-prefix guess. Next explicit calibration variant preserves all previous source files and accepted Hunt A.

Reflex-demo provenance audit caught cross-family mutable parameter leakage in diagram-export metadata (syn6 heel/toe defaults). Candidate matrices and recorded neural parity were correct, but the first reflex batch used altered contact gains. Those v1 results are disqualified and retained. Independent-process extraction now reproduces the accepted biological physical prefix exactly; corrected v2 pulse cases are being quantified. Source models and manuscripts remain unchanged.

External records remain authoritative; no duplicate root simulations or excessive new logs.

## Continuation checkpoint — 2026-11-01 16:12 UTC

Verified corrected Gait2392 sensory-spiking reflex results: two 16-cell populations emit 64 actual events per 100 ms pulse. Ia increases soleus MN by 1.843 mV and tension by 908.97 N, with a 2.346-degree ankle effect. Ib initially inhibits by 0.1565 mV and 56.14 N, with later rebound. Outgoing-path-off controls preserve spikes but give exactly zero downstream response. Both fresh repeats are exact; corrected baseline matches the archived biological physical prefix. These are short functional tests, not an all-spiking controller or new standing winner. Figures and record: external spiking/Gait2392_reflex_demo_v2.

Gait rolling-stance load transfer reaches 10% BW in both directions with preserved earlier prefixes, exact fresh replay, and all 1 ms acceptance gates. Donor feet still carry 8.1%/11.8% BW, so no swing. Toe acceleration remains timestep-sensitive. Next bounded continuation is 12.5% then 15%.

CAD frozen pressure pair passes one support reduction: actual upward apparatus force falls from 91.87% to 87.34% BW, with exact fresh repetition, knees flexed, and peak pressure 479.81 kPa. Ideal-tension comparison fails its COM-drift gate and is retained. Posterior COM excursion remains substantial; next task audits actual support geometry/AP harness forces before reducing assistance further.

Rat receptor variant qualifies all eight RG input waveforms across seven conditions, but full teaching remains blocked by eight other mismatches and target force-length domain defects. Next audit uses actual biarticular gastrocnemius function and source provenance; no gain fitting is permitted on flat or physically inconsistent teaching channels.

No manuscripts or accepted production controllers changed; separate records and external raw data remain authoritative.

## Continuation checkpoint — 2026-11-01 17:17 UTC

Frozen-family spiking reflex comparisons are complete for s3k, normalized w2lvar, and syn6. Positive Ia/Ib responses repeat exactly; outgoing-source-off controls retain spikes but eliminate downstream effects. Family-specific pathways are documented rather than treated as matched topology-only experiments. Source-derived figures/table: external spiking/Gait_family_reflex_v1/EXPERIMENT_RECORD.md. Next spiking task targets actual RG/PF population replacement with spike-derived chemical outputs, first offline.

Gait2392 reaches 12.5% dynamic redistribution in both directions, with fresh and 1 ms gate checks, but cannot hold the left milestone: even unchanged commands continue lowering the pelvis. Fifteen-percent attempts remain failed. No unloading/swing or Hunt C claim. Next task returns to full-native stationary load-bearing solutions at the achieved split instead of extending an unstable transient.

CAD posture audit finds centimeter-scale initial four-bar loop separation, 5.5–6.5 mm tail gaps, approximately 29–30 degree toe-out, and non-grounded initial left foot. The trunk is nearly upright; posterior translation and harness moment balance explain the apparent lean. A feet-fixed 10 mm anterior proposal fails actual static support/route-force feasibility and remains rejected. Existing pressure-interface successes remain scoped to their compliant, guided model. Next task repairs initialization/stance constraints in a distinct mechanical reference without changing CAD geometry or neural weights.

Rat author-method curve/rest reference passes all seven fixed-kinematic domains, and the actual current graph lacks the published GA integrator. Next explicit 134/274 variant adds the documented knee/ankle GA integration pathway and requalifies all sensory inputs before training. Prior source132, accepted Hunt A, and all failed maps remain preserved.

No manuscript edits, production replacements, or duplicate root jobs; records retain all material limitations.

## Continuation checkpoint — 2026-11-01 18:18 UTC

Both Gait2392 12.5% load-split references now pass full 105-coordinate equilibrium, five-second stationarity, exact fresh repeats, and 1 ms checks. Direct activation interpolation between equilibria fails, confirming that static capacity alone does not define a dynamic motor primitive. Next task is a bounded augmented-dynamics/equilibrium-manifold transition at a smaller split with terminal stationarity; source physics and all prior gates stay fixed.

Central-spiking Gait prototype has 256 central LIF cells with 778 graded cells; none of 24 candidates passed both standing and rhythmic replay, so no body trial ran. Next isolated mechanism test lets NaP inactivation observe only the spike-derived population envelope, preserving actual spike chemical outputs and the fixed numerical clock. This tests the identified reset-voltage distortion rather than repeating gains blindly.

CAD stance repair finds exact grounded, zero-heading geometry without CAD edits, but source soft constraints require small loaded deflection. Exact zero-gap statics and the frozen virtual BPA domain both fail. Next distinct mechanical-only ideal-force reference allows at most 0.5 mm physical compliance while retaining the native constraint law; no existing sleeve geometry is resized or pressure result reclassified.

Rat editable134/274 GA integration variant passes Hunt A/native neural parity, but only32/44 physical input qualifications. Two LH receptor rest lengths exceed their full route range and remain tensionless; prototype Ib/II labels actually bind Ia. Source-binding provenance and a declared intrafusal-only geometry normalization are the next prerequisite, not another gain search. No full B or body/drop result is claimed.

Four bounded tasks continue, with no duplicate root simulations, manuscript edits, or production replacements.

## Continuation checkpoint — 2026-11-01 19:15 UTC

No new controller promotion. Gait manifold study adds stationary 1% references with exact fresh and 1 ms checks, but its best dynamic candidate still fails settling. The audit shows the global six-mode torque reduction discards most distal row authority; this is not evidence against the full 46-family basis. Next work audits full/task-normalized derivatives and contact-sensitive objective blocks before any additional fit.

Central-spiking mean-field inactivation restores oscillations and h-phase, but all 24 refinements fail combined standing/rhythmic fidelity; zero body trials. Next structural test derives a conductance-aware spike-only RG readout instead of further gain sweeps.

CAD single 0.5 mm-compliance initialization still fails actual equilibrium and loop bounds; five-second tests were withheld. Next audit measures native equality compliance, row rank, and force accounting on a small known-wrench fixture. Any tighter rigid-link numerical approximation will be a separate derived variant, not an overwrite or pressure-capacity repair.

Rat receptor-only geometry/calibration variant qualifies 35 informative inputs across seven cases; four low-drive baselines and five unused signals remain explicitly classified. All original ordinary muscle forces, bindings, and central graph remain unchanged. Unfitted closed neural replay latches only LH hip PF. Next causal target is its two unmatched prototype-Ia returns, with source provenance checked before any topology policy change.

Separate records and snapshots are current. No manuscript edits, duplicate root jobs, or continuation past failed physical gates.

## Continuation checkpoint — 2026-11-01 20:13 UTC

No new stepping result. Gait full46/1ms audit passes physical guards but left force/velocity derivatives remain irregular. Verified zero-effect shadow localizes predictor/corrector FV-branch/tangent disagreement near MTP spikes without proving causality. Next explicitly separate integration variant seeks within-step branch/tangent consistency under unchanged source laws; original A403 backend remains protected.

CAD fixture demonstrates that old soft equality settings can permit millimeter loaded gaps. A single derived numerical variant gives tens-of-micrometers static gaps; the stronger load peaks at115µm, missing its100µm design target but satisfying the0.5mm condition. Native EOM and timestep checks pass. Next full-model mechanical static audit uses that frozen variant, ideal54 forces, and repaired initial geometry; no pressure/hardware claim or new stiffness selection.

Rat two-prototype-return causal diagnosis is decisive: both contribute6.74nA to LH hip PF; either alone removed fails, both removed restore all basic RG/PF cycling. A separately declared bilateral-role policy masks those two existing LH edges (none exist RH), retaining localIa/allII. Source files remain untouched. Full unfittedA/B gates run before any conditional12+12 study or actualair run.

Central-spiking inverse readout improves static fidelity but fails dynamic envelopes. Actual slow burst cadence is2.326s; rapid peaks are readout splitting, not fast biological bursts. Next mechanism uses causal conductance-matched spike history rather than another gain search. No physical central-spiking trial is authorized before offline gates.

All four agents have bounded next tasks; no running jobs were duplicated, manuscripts edited, or failed candidates promoted.

## Continuation checkpoint — 2026-11-01 21:21 UTC

Reviewed central-spiking burst scoring against source code and saved traces. The old local-maximum counter splits within-burst ripple; a common 100 ms envelope/hysteresis detector gives source2.132 s versus candidate2.130/2.149 s, while the original unfiltered waveform error0.14977 already passes. Root approves a separately versioned burst metric, preserving all original-negative records. Missing/doubled/irregular bursts, isolated events, duty, and within-leg phase checks are required before short physical standing. This is an after-the-fact scorer correction, not retroactive preregistration or a walking result.

Whole-step consistent source integration closes forces/energy but leaves the MTP peaks and settling failure unchanged at1/0.5 ms. This refutes the proposed predictor/corrector mismatch as a sufficient cause. Further work returns to bounded actual nonlinear trajectory scoring, retaining distal task authority and rejecting all failed terminal conditions; no more tangent-cause tuning.

CAD tighter-constraint full-model static reference still fails native acceleration and force balance despite acceptable loop gaps. The old5%-MIF force floor materially alters feasibility and is not5%physiological activation. Next audit separates force-only mechanical fixtures from source-activation coactivation semantics before one new bounded equilibrium attempt. No neural or pressure result is inferred from LP feasibility.

Rat role-policy B has stable cycles and reduced amplitude error but still fails waveform. Same-config A remains exact. Next limited review checks teaching phase/actual sensor lag before a bounded refinement; no independent per-neuron warping or relabeling of old failures.

No new physical controller promotion, manuscript edit, duplicate job, or permission needed from Ben.

## Ben returned for review — 2026-11-01

Hourly automation continue-neuromechanical-campaigns is PAUSED. All numerical campaigns stopped after current bounded work was saved. Review media were rendered from saved poses, not new simulations: D:\Users\Ben Bolen\NeuroCampaign_20261031\Review_20261101\GAIT2392_REVIEW.md and cad\README.md.

New central-spiking Gait2392 eight-second standing passes the unchanged standing criterion; all82 fresh saved arrays are identical. Direct measured support is upward0.5937%BW, downward0.0440%BW, net0.5497%BW, feet99.3228%BW. This supersedes the earlier0.719%BW scorer estimate. Selected256central-spike/778graded configuration and source hashes are in Gait_dynamic_readout_v4/TIMING_CANDIDATE.json. Lateral/AP/rotational apparatus remains.

Review stepping GIF uses the corrected-initial-pose graded Hunt-PF checkpoint at63.7263%BW assistance,5/4 sustained mesh-cleared swings. Earlier/lower-support runs are not silently substituted. Same configuration has not passed its preceding standing stage. Source-compliant50%BW stepping remains a separate anatomical-ankle-envelope failure.

CAD pressure-pair GIF is the complete support085 record:2modeledBPAs/52hybrid routes,87.337%BW upward support. Visual review confirms distal crossing/overlap; the historical settling gate did not test uncrossed stance. Manifest now flags this alongside toe-out and loop gaps. It must not be presented as corrected autonomous robot balance.

The old8765 localhost server serves an empty directory. A new loopback-only server at8766 serves the permanent Code\MuJoCo_SNS\spinal folder; editor returnsHTTP200 and is open as an IAB deliverable. Existing tabs/localStorage were not changed. Port8766 has a distinct origin, so older browser-saved tabs are not migrated automatically.

Source/literature review: virtual_walker_20261031/LITERATURE_SYNERGY_3D_REVIEW.md. Submitted311-page Defense PDF lacks full NMF/conductance/objective derivation; reviewed Chapters6 andD describe procedures, not enough numerical controller-design detail. No manuscript or Zotero edits were made. Exact Micklem2026/Schwarzenberger identity requested from Ben; no unverified substitute attribution.

## Supplied-paper and scaling review — 2026-11-01

Ben supplied Micklem2026, Scharzenberger/Hunt encoding/CPG papers, Scharzenberger/Hays ROA, Szczecinski2020, and Sutton2021 PDFs. Earlier unresolved-name status is superseded by direct PDF review; no guessed substitution is needed. User expanded Table2audit to all neural controller families, with a separate syn6/DiRusso comparison.

Completed read-only records:
- Code\MuJoCo_SNS\spinal\virtual_walker_20261031\SZCZECINSKI2020_SEVEN_STEP_CONTROLLER_AUDIT.md: actual seven-step/full-family matrix, exact code evidence, GLIF versus graded applicability, 1000x capacitance-unit adapter.
- Code\MuJoCo_SNS\spinal\virtual_walker_20261031\SUTTON_HUMAN_SCALING_REVIEW.md: size/time dimensionless analysis, fastest running caches/harmonics, and proposed scaling process.
- D:\Users\Ben Bolen\NeuroCampaign_20261031\Review_20261101\CONTROL_DESIGN_LITERATURE_REVIEW.md: supplied Scharzenberger/Micklem/ROA equations and page-level evidence.

Fastest usable cached stride is subject11_Run_50002 cachedl0.559s (1.78891strides/s); pairedr0.564s. All43caches load;32runningtrials/64cycles. Original raw cycle completeness/anatomical side labeling remain unverified. The100point phase grid is not raw sampling rate. Fast ankle99%AC modulation reaches harmonic4/~7.156Hz. Prior'Arnold' labels are not established dataset attribution; catalog points towardHamner/Delp.

No whole current spiking controller has evidence of full2020Table2completion; gradedmodels require corresponding2017functional audit, not literal adaptive-threshold variables. Syn6 directlyuses2017Eq18, while six synergy factors are NMF-derived, not DiRusso's five prescribedcosine primitives. No models, fitted parameters, physics, manuscripts, or Zotero records changed. Hourlycampaign remainsPAUSED.
