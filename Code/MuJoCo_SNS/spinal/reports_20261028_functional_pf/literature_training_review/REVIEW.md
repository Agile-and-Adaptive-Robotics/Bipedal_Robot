# Deng 2023 and Hunt's INRC training method: Session 1 correction

Machine: **easteregg2**. Repository: `D:\GitHub\Bipedal_Robot`.
This is a controller-method review and bounded experiment, not a dissertation
revision or a new successful walking result.

## Sources actually inspected

- Alexander Hunt, *Real-time control of dynamic walking in a dog-like robot*,
  INRC, 15 February 2019. User-supplied source:
  `C:\Users\Ben Bolen\Downloads\19_02_15 - INRC.pptx`.
  Extracted all 57 slides and their speaker notes in presentation order.
  Training guidance is on slides 39-46; inverse BPA force solving is on
  slides 48-49. The embedded circuit and training plots on slides 44-46
  were inspected separately. PowerPoint could not open the deck for native
  slide export; full slide composites and embedded movies were not reviewed.
  Extraction does not modify the presentation.
- Deng, Hunt, Chiel, and Quinn (2023),
  [*Biarticular Muscles Improve the Stability of a Neuromechanical Model of
  the Rat Hindlimb*](https://doi.org/10.1007/978-3-031-39504-8_2), Living
  Machines 2023, pp. 20-37. Read the user-supplied author manuscript at
  `D:\GitHub\_zcode_backups\lit_20261026\kaiyu2023.docx`, identified by
  `Code\Matlab\SNS_Simscape\dev\DENG2023_SOURCES_20261025.md`.
  Extracted paragraphs, both parameter tables, and all seven figures;
  inspected the original vector circuit figure after rasterization.
  Paragraph references below use the local extraction, not published page
  numbers. This is the supplied manuscript, not a verified final-typeset PDF.
- Rechecked [Deng et al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6477610/),
  Sections 2.2, 3.1, 3.2, and 4.2. Its shared knee/ankle PF precedent remains
  valid for w2lvar. Its fitting prioritizes MN activity, then kinematics,
  emphasizing peak timing and magnitude. It also acknowledges transient
  feedback-related phase shifts after PF perturbation in the intact model.

Raw source extractions and artwork are local, gitignored inspection files.
Source hashes and reproducible experiment records accompany this review.

## The missing training distinction

Hunt's slide 44 specifies four functional training steps:

1. Choose a free-running CPG with approximately **twice the target stepping
   period without drive**, equivalently half the target stepping frequency.
2. Replay expected afferent feedback and fit the input connections so that
   transitions occur at the target walking speed.
3. Fit CPG/PF-to-MN output to the desired activation peaks.
4. Fit Ia/Ib feedback to improve the remaining activation waveform.

Slide 43 places most nominal MN activity in the feedforward component, with
sensory modulation providing adjustment. Its 70%+ figure is presentation
design guidance. It is not a minimum coactivation level, a universal human
physiology value, or an automatically additive fraction of SNS output.
Conductance-based inputs shunt one another; the contribution measure must
be defined before turning that guidance into a numerical objective.

Slide 46 explicitly warns that training the whole network together can
produce a local solution that depends strongly on the expected feedback.
Slides 40-42 supply the preceding chain: desired motion -> muscle tension
-> desired MN activity. Slide 48 revisits the inverse force problem for
robot actuators. Therefore IK alone is not the MN-training target.

The completed 500-trial campaign changed physical environments and searched
many gains together. It did **not** implement these four functional fits.
It did not verify the undriven intrinsic period, independently fit expected
afferent transitions, or fit the new PF-to-MN activation peaks before closing
the loop. Its short inverse-force estimate and NMF analysis remain exploratory
engineering targets, not measured EMG or identified biological synergies.
The requested five physical environments remain useful; they are a different
dimension of the training and validation procedure.

## What Deng 2023 actually adds

The rat model has five monoarticular muscle groups and three biarticular
groups: BFP, RF, and gastrocnemii. The model is sagittal. The manuscript says
that the authors first manually tuned weight-bearing forward movement, then
suspended the pelvis for the reported perturbation comparisons. The air tests
do not establish unsupported 3D Gait2392 or robot walking.

Figure 1 separates **joint PF oscillators** from **biarticular integration**:

| Element | Published-model description | Difference from the previous functional prototype |
|---|---|---|
| Hip, knee, and ankle joint PF | Separate extensor/flexor oscillator populations with inhibitory relays; PF oscillators have persistent sodium dynamics | w2lvar's original knee/ankle drive is shared; the added ankle and other PF cells are leaky, not independent NaP oscillators |
| BFP intermediate population | Combines hip extension and knee flexion PF signals | HAM channels previously received direct generic RG stance/swing drive |
| RF intermediate population | Combines hip flexion and knee extension PF signals | RECT channels previously received direct generic RG stance/swing drive |
| GA intermediate population | Combines knee flexion and ankle extension PF signals | GAS channels previously received direct generic RG stance/swing drive |

Figure 1E has two intermediate populations for BFP/RF together, one for each
muscle, not an independent four-cell oscillator per biarticular muscle.
Table 2 gives hip/knee BFP input conductances 0.967/0.365 microS, hip/knee RF
0.236/0.872 microS, and knee/ankle GA 1/0.3 microS. Their thresholds also differ.
These numbers cannot be copied directly into the current normalized-voltage
human network. Nor is rat BFP anatomically identical to all human hamstrings.

RG-to-joint-PF conductance decreases from hip to knee to ankle
(0.1, 0.05, and 0.02 microS), allowing different phase offsets and disturbance
responses while retaining a common approximately 0.59 s period. That is a
phase-locking/flexibility tradeoff, not a prescription for human cadence.
The paper's voltage-dependent NaP time constant also differs from the current
fixed-tau implementation. A future paper-faithful NaP-PF branch requires its
own solver verification; the existing fixed-tau fast backend must not silently
be used to change those equations.

### Results that change our interpretation

- **PF perturbations can affect RG timing in the complete loop.** A knee-PF
  current injection changes limb motion and sensory feedback, advancing the
  RG and other joints. Figure 5 reports approximately 5.2% persistent phase
  advance after the initial disturbed cycle. Thus RG/PF separation does not
  imply that perturbing PF can never reset a closed-loop neuromechanical model.
- **Mechanical coupling and active recovery are distinct.** The authors compare
  BFP/RF active, active-force-off with passive tension retained, and both forces
  disabled. Active biarticular muscles reduce the disturbance and restore timing
  faster than the disabled condition. Passive tension can limit disturbance
  while failing to restore the original step timing.
- **The comparisons have limitations.** The nominal RoM differs across conditions,
  the conditions are not independently optimized for ground walking, and the
  muscle model differs from MuJoCo's. Discussion explicitly limits generalization
  to realistic ground locomotion. The existing PF-to-MN edge deletion is a
  neural-pathway ablation, not this mechanical active/passive experiment: reflex
  and postural inputs can still activate the same MNs.

Evidence: manuscript P016-P017, P029-P046, P048-P068, Figures 1, 3, 5, and 7,
and Tables 1-2. The manuscript occasionally says BFA where the stated experimental
conditions say BFP, and the torque sign differs between body text and caption.
Use the explicitly defined BFP/RF conditions and torque **magnitude** 0.1 Nm;
do not silently treat those inconsistencies as independently verified values.

## Changes implemented

1. Added an opt-in `deng2023_joint_pf.py` experimental variant. Each leg has
   three single premotor integration populations for HAM/BFP-like, RECT, and
   GAS functions, driven by the corresponding joint PF signals above. Input
   ratios are normalized within a configurable human-network gain budget.
   Other biarticular families retain the previous experimental routing.
2. Explicitly retain the adaptation boundary: joint PFs remain leaky, sodium
   kinetics remain the current human implementation, source
   threshold offsets are not reproduced, and MN gains are inherited rather
   than retrained. The three changed families retain base-W2L reflex paths,
   including its IaIN phase gating, rather than the prototype's extra direct
   Ia/II inputs to biarticular PF. Changing both integration and those pathways
   makes this an architecture probe, not a single-edge causal comparison.
   This is a partial human adaptation, not a Deng replication.
3. Added pelvis-frame ankle/knee ordering to future evaluation acceptance.
   The 20 mm body-origin separation threshold is an engineering screen, not a
   full foot-edge clearance/self-collision test or a value from Deng's paper.
4. Future curriculum runs stop after a failed fresh replay. Starting at a later
   stage requires a passed predecessor; ground promotion also requires the
   frontal gate. `--diagnostic-continue` explicitly permits diagnostic exploration
   without certifying progression. Final selection prefers a feasible trial
   when one exists, rather than choosing a higher-scoring failed trial.
5. Updated the shared neuroscience skill with the four training steps, the
   biarticular integration distinction, and the qualified phase-reset rule.
   Production family builders, their selectors, and the completed 500-trial
   campaign's stored evidence are unchanged.

## Anatomical-axis correction

A neutral-pose regression test exposed an error in the previous frontal audit: it rotated a world-right vector as if it were pelvis-local, making the measurement nearly vertical. The corrected method derives right from the hip anchors and verifies local +z. Neutral ankle/knee widths are 0.18284/0.167 m, and a crossed-pose fixture is rejected. Earlier crossing percentages and widths are withdrawn. Saved finalist trajectories and all 16 final perturbations still fail the corrected gate; earlier corrective-probe rankings need separate re-audit. `corrected_frontal_reaudit.json` records the corrected values without rerunning physics. The table below uses the corrected axis.

## Bounded verification and probes

The new variant has **962 neurons and 400 inputs** before the finite initialization
port, compared with the earlier combined prototype's 980/400. Its six intermediate
populations each receive exactly two appropriate joint-PF projections. Those muscle
MNs receive their intermediate drive rather than the old hip/knee pattern drive.
The accelerated and preallocated reference calculations agree exactly for every
neural voltage and NaP h-gate over 1,000 steps. Building the new variant between
two original-variant builds leaves the original traces exactly unchanged.

Two untuned, 16 s air probes inherit the old air finalist's other parameters:

| Joint integration gain | Knee cycles R/L | Period (s) | Mean error from half-cycle separation | Hip RoM R/L (degrees) | Minimum ankle separation (m) | Accepted |
|---|---|---:|---:|---|---:|---|
| 0.3 | 5/5 | 1.752 | 0.0123 cycles | 171.4/158.3 | -0.537 | No |
| 1.0 | 5/5 | 1.755 | 0.0507 cycles | 146.2/146.5 | -0.314 | No |

Both remain finite but retain crossing and excessive excursions. They do not
improve the old finalist or validate the new architecture. Inheriting output
gains is not equivalent to Hunt's dedicated output-pattern training. No later
ground stage was promoted from these failed air trials. The stance-fraction
parameters for HAM, RECT, and GAS are explicitly recorded as inherited but inert
in this variant, rather than presented as active optimization variables.

Verification: `joint_integration_verification.json` and
`promotion_gate_verification.json`. All three tested promotion failures
(missing predecessor, failed predecessor, and missing frontal gate) stop before
any simulation trial begins. Probe configurations, metrics, and trajectories
are under `joint_integration_probes/`. A repeated gain-1.0 condition checks the
complete trajectory against the slower evaluator. All 13 saved arrays agree
exactly, including full joint states, activations, commands, neural voltages,
and contact. This is validation, not a third tuning condition. See
`joint_integration_trajectory_parity.json`.

## Revised experiment sequence

The machine-readable training specification is in `TRAINING_PROTOCOL.json`.
First establish model-consistent motion, force, afferent, and MN targets;
then fit intrinsic timing, sensory transitions, PF/MN output, and reflex
residuals as separate blocks. Freezing preceding blocks is the implementation
choice used to preserve the four-step separation. Failed closed-loop gates
return to the responsible block and require earlier gates to be rechecked.

Keep spatial recruitment, temporal activation primitives, and joint-function
classification separate. Secondary actions and hip-rotation sign changes must
use the measured pose-dependent moment arms. An NMF factorization alone does
not identify a neural circuit. Human actuators and the CAD robot's BPAs require
different force-to-MN mappings; do not transfer a human activation fit unchanged.

After those fits, validate the five requested environments, including stable
standing, uncrossed legs, physical alternating cycles, contact, and assistance.
Only then compare phase-controlled RG-current, PF-current, and mechanical
perturbations under active, passive-only, and disabled biarticular conditions.
Each mechanical condition needs its own unperturbed reference. Scale disturbance
torque to the human model explicitly instead of assuming rat 0.1 Nm is comparable.
Record phase error, amplitude disturbance, recovery cycles, contact changes,
and support forces separately. Fresh-process all-scalar +/-1% basin testing and
held-out afferent/perturbation tests remain required before promotion.

The four-step functional fitting is **not complete**. These changes correct the
architecture interpretation, trial gates, and next training procedure; they do
not turn the earlier 500 trials into a graduated controller. No dissertation
chapter, biomimetic result, defense figure, or live Overleaf source was changed.
