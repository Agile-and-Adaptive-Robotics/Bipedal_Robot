# Session 1: force provenance, symmetry, and supported stepping

Machine: **easteregg2**. Repository: `D:\GitHub\Bipedal_Robot`.
Interpreter: `D:\Anaconda\envs\myo\python.exe`. MuJoCo 2.3.7,
implicitfast, 2 ms physics/neural steps, accelerated sparse/pairwise evaluator.
This is the **human Gait2392** campaign, not the CAD/BPA robot.

## Strength provenance and physics

- Source OpenSim bilateral force specifications agree. MyoConverter Step 3
  independently fitted four force-law parameters for each muscle. Its force
  multiplier changes the stored Fmax. Soleus, for example, became 3395.5 N
  right and 3186.36 N left from the same 3549 N source specification.
- Full-activation fitting can conceal opposite active/passive errors. Neutral
  BFSH has similar full-activation totals (736.19/732.60 N), but at 5% activation
  the original surrogates produce 241.60/157.34 N. Simply copying right-side
  parameters worsened held-out source agreement in four of six initial families.
- The existing runner also repairs eight converter active scales that arrived
  as 1 N: ercspn 2500 N, intobl/extobl 900 N, and ext_hal 162 N, bilaterally.
  Their suspicious converted length ranges remain a separate alignment issue.
- `physics_fidelity/shared_force_law_functional_v3.json` fits common parameters
  for 18 families/36 muscles against matched OpenSim static forces. All 92
  muscles have NOT been recalibrated. RF held-out RMSE falls from 259.69 to
  22.78 N after fitting the actually repaired route. This is a static surrogate,
  not an elastic-tendon model or a full-ROM/dynamic validation.
- Coordinate adapters matter: matched pose/derivative comparisons use hip
  flexion and adduction sign -1, knee/ankle +1. The apparently reversed signs
  disappear after this explicit adapter; remaining moment-arm errors are
  quantified in `../physics_fidelity/MOMENT_FORCE_AUDIT.md`.
- The physics subagent's compliant-tendon prototype agrees with OpenSim static
  equilibrium to 2.16e-6 N. A warmed 92-muscle bank takes 2.47 s per 16 simulated
  seconds versus 0.0776 s for its rigid comparator. These are bank-only timings.
  Dynamic force-velocity and added damping remain experimental. The compliant
  bank is NOT integrated into the accepted SNS walker.
- An unused legacy nonlinear-contact option misindexes solimp and has equal
  impedance endpoints. It does not provide its claimed depth-dependent damping;
  this does not establish that the current default contact law is defective.

## Neural symmetry

The old left/right sensory-timing fit was a mistake: hip-extensor length gains
were independently fitted at about 13.0181 and .0026777. New B fits pool both
legs. PF/BI/MN/reflex learning pools mirrored rows and columns, fitting one
conductance per canonical pair. New Hunt training defaults to tied fitting;
asymmetric timing inputs are rejected. Explicit HUNT_BILATERAL_TIED=0 is reserved
for historical reproduction. Saved strict candidates reject mismatched weights.
A regression distinguishes the pooled optimum 3.5987 from post-fit averaging
(3.0), and confirms one solver call for the pair. Balance recruitment can also
be pooled sagittally and mirrored laterally with bilateral_recruitment=true.
Global midline directional channels still need an anatomical-axis-specific
learning audit before claiming unrestricted 3-D control.

The both-cadence tied fit in
`runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact` passed stages 1/2
and fresh-process repeats exactly. Subsequent force-law and recruitment changes
require new validation; do not substitute that older early-stage checkpoint.
The added strong E/F commissural pair caused double knee pulses in the newer
air trial. Removing motor reflex strength or timing feedback did not fix it;
removing that added pair did. The tested alternative keeps crossed hindfoot
landing feedback to the contralateral RG-F, with shared timing_scale=.01.

## Acceptance corrections and current ground evidence

Separate hindfoot/forefoot forces are partitioned from actual contact positions;
normal forces are counted once. Landing events require prior whole-foot unloading
and produce a decaying pulse rather than a sustained RG reset.

Force unloading is insufficient: some nominal 30-40% support runs passed the
old motion/force gates while dragging both feet. New stage-4/5 acceptance also
requires at least three sustained swings per foot in 5-13 s: GRF <.025 BW,
the entire visual foot mesh >1 mm above the plane for >=80 ms, and peak >=10 mm.
These are explicit engineering margins, not a published biological standard.
Historical results are not overwritten; independent reassessments are saved as
`*.ground_clearance_gate.json` and `GROUND_CLEARANCE_REASSESSMENT.json`.

Ben's no-hyperextension directive is enforced in the evaluator on BOTH knees:
0-degree upper range, protective stop engaging 2 degrees before neutral,
solref .004/1, solimp .999/.9999/.001/.5/2. No position or velocity clipping.
Every saved physics state and final state is checked; ANY positive knee angle
rejects the run. Earlier passers lacking this check are historical only.

`runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/zero_hyperextension_v1`:

| Trial | Measured vertical assistance | Sustained swings R/L | Maximum extension R/L | All current ground gates |
|---|---:|---:|---:|---|
| ground_original | 76.43% BW | 3/4 | -1.324/-1.423 deg | pass |
| ground_no_comm | 61.74% BW | 5/4 | -1.315/-1.408 deg | pass |
| ground_half_support | 43.61% BW | 5/4 | -1.324/-1.398 deg | pass |

The no-commissural-pair controller also passes both air stages with protected
knees in `zero_hyperextension_air_v1/air_pair/SUMMARY.json`. Its standing check
still fails. These remain supported stepping-in-place trials: horizontal and
rotational restraints, prescribed support, and protective knee stops remain.
There is no harness-free walking claim and no completed five-stage curriculum.

## Contact-feasible starting pose: validation in progress

The physics follow-up found the old initial feet penetrate the plane by roughly
35-58 mm, causing initial reactions around 6.6-7 BW. Raising height alone leaves
heels about 20 mm above the floor, creating a tiptoe pose. Alignment must be
performed on repaired model axes before constructing actual harness references.
The tested helper is `physics_fidelity/ground_alignment.py`,
align_via_repaired_probe. A common native ankle adjustment of about +5.62 deg
and vertical shift about +38-39 mm brings the regional foot surfaces near the
plane. The exact shift depends on requested initial gap and pose.

The subagent obtained stable eight-second CONSTANT-ACTIVATION holds with this
pose and its new contact-aware initializer; these are NOT SNS standing passes.
SNS tests are in `contact_feasible_start_v1`. Do not choose a final defense
animation or claim completed curriculum until the corrected initial pose,
zero hyperextension, and geometric clearance are checked together.

## Patella/inertia clarification

Verified in the actual repaired model: vas_med-P4 is a moving surrogate body
under the femur, driven by knee-angle equality polynomials. RF routes through
it to P3 under the tibia. The route transfers muscle tension, not tibia/foot
inertia into the via point. Those bodies retain their own inertia. Follower
armature adds reflected inertia A(ds/dq)^2; measured right-knee contributions
are .03256 kg m^2 neutral and .02570 at swing, 6.19%/5.35% of reduced inertia.
The old 36 kg m^2 statement was wrong. These percentages describe follower
ARMATURE, not every possible surrogate-body mass contribution.

## Reproduction and next work

Use `physical_eval.py` with a copied JSON in a NEW ignored run directory.
Retain absolute weight/force-profile paths. Never overwrite original evidence.
The entry point is
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\physical_eval.py`.
Track follow-up results in TRAINING_RECORD.md and this session checkpoint,
not CHATGPT_REPORT.md. Freeze v3 and accepted weights; physics prototypes are
separate. No manuscript/defense-source edits, commits, or pushes were made.


### Corrected-pose checkpoint (supersedes the unaligned lower-support claim)

Fresh repeat: aligned_support_ladder_v1/repeat_0p5, same saved controller as
protected-knee air passers. All current ground gates pass, including actual
foot-mesh clearance and no positive knee angle. Measured assistance63.73%BW;
sustained swings5right/4left; maximum knee angles-1.201/-1.460degrees.
Recorded qfull/ctrl/neuro are bit-identical to the first corrected-pose run.
The initial feet are aligned near the floor before building rig references.
This also rebases the existing ankle spring rest under the runner convention.

The earlier43.61%assistance run used the bad deeply penetrating initial pose;
it is NOT the current checkpoint. After correcting the pose, lower prescribed
support.3/.1 produces measured.556/.598BW assistance and ZERO qualifying
swings per foot. Those runs are rejected; the guide takes more load as the
model collapses. No successful further support wean is claimed.

Corrected-pose standing with the original coupled initializer improves to
12.69degrees tilt and93.6mm sway, load.477BW, but still fails the60mm/.5BW
limits. The contact-aware initializer succeeds in constant-activation holds
but fails in the actual SNS interface (about80mm sway/.17BWfootload), and its
ground test also fails phase criteria. It is an experimental opt-in variant,
not the selected checkpoint. Next investigate delivered MN activation versus
initializer targets, gait-trained reflex activity in standing, and stance
support capacity. Do not sweep insensitive gains or count static holds as SNS
curriculum passes. No stage5 advancement yet.

Full details: SESSION1_PHYSICS_GROUND_CHECKPOINT_20261030.md. Canonical ground
checkpoint: tonic_0p02_exact/aligned_support_ladder_v1/CHECKPOINT.json.
