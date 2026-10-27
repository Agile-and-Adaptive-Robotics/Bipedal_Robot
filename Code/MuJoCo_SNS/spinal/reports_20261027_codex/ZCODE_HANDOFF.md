# Prompts for continuing the simulation recovery

Use these prompts in ZCode with this checkout. Read `REPORT.md` first. The
verified W2L result and actual CAD model are available; Di Russo replication and
CAD robot neural locomotion remain open. Do not reclassify either as successful
on the basis of a long duration, a favorable scalar score, or a supported
mechanical test.

## Di Russo original baseline and gait2392 adaptation

Continue task 2 in `D:\GitHub\Bipedal_Robot`. Check hostname first, read AGENTS
and the Python, OpenSim, MuJoCo, and neuroscience skills, and preserve the dirty
working tree. Do not commit or push.

Use the actual original sources:

- `https://github.com/DiRussoAndrea/Spinal_controller`
- `https://gitlab.com/simgait/SCONE`

The first is staged under `tmp/codex_dirusso_original/`; the second is a sparse
checkout under `tmp/codex_simgait_scone/`, revision
`58ea68762791e133b0b18f2f0db83e3eb2a18fe3`. The neumovelab comparison framework
is not the author's source. Do not replace the published phase oscillators with
a Matsuoka controller or the Thelen muscles with an Ogihara model on that
attribution.

First establish a reproducible replay of the author's optimized `spinal.par`
on `gait9dof18musc_Thelen_20170320.osim` using the custom SegmentalController
runtime. Record the source revision, executable, dependency versions, complete
parameter file, initial state, objective terms, and unassisted trajectory. The
model initializes in this machine's OpenSim 4.4.1, but that does not establish
controller compatibility. The fork's supported build route is Ubuntu 18.04;
ordinary SCONE installations need not contain SegmentalController. If the native
runtime cannot be built or run, report that exact prerequisite and preserve the
source; do not spend another optimization campaign compensating for a guessed
controller.

Audit all 256 optimized parameters against original source and scenario:
one frequency, 16 state offsets, two phase offsets, nine balance parameters,
45 pattern-formation weights, 174 reflex connections, and nine length offsets.
Verify delayed neuron state histories and update ordering, muscle-specific
delays, fiber-length and velocity normalization, afferent scaling, saturation,
Renshaw-to-antagonist-IaIN wiring, balance gates, and heel-reset semantics.
Use recorded original neural and actuator traces as regression evidence before
retargeting the controller to gait2392. Keep the original baseline and the
gait2392 adaptation clearly distinguished.

For the MuJoCo adaptation, preserve 2.3.7, implicitfast, the approved 2 ms
physics step, initialization, contact repairs, and support-rig definitions.
Use compiled neural evaluation and persistent process workers only after
trajectory parity; keep one BLAS thread per worker. Do not use the old FastTrial
loop or the old sparse backend as an equivalent physics runner.

The previous 960-trial survival/speed search produced a somersault and is
rejected. `optimize_dirusso.py` deliberately requires `--debug-surrogate` for
that experiment. Replace its objective with the original objective terms and
explicit pose, contact, and progression acceptance checks after baseline
reproduction. A candidate must pass fresh-process replay, ±1% basin checks,
finite-state and numerical-warning checks, per-term diagnostics, and visual
inspection. Record external support forces; zero nominal support is not
harness-free if PD forces remain active.

Deliver the replay command, parameter file, numerical comparisons, plots,
failure cases, and an honest statement of which published results reproduce.

## Actual robot with lumbar hinge

Continue task 5 from `Code/MuJoCo_SNS/robot_torso/README.md`. Use
`ah001_with_torso.xml`, not the intermediate rigid mesh-group URDF. The source
assembly is `Solid_Models/00_BipedalRobot_Redesign/10_AH_001.SLDASM`; it was
opened read-only and was not saved. Both four-bar knees and both spherical hips
are represented. Preserve the two loop constraints and their rank-four
Jacobian. Do not replace the knees with arbitrary single hinges or weld them
through generic CAD constraints.

The user approved a 250 mm rigid spine, a sagittal lumbar hinge at the adjustable
upper rear center of the pelvis, total mass 75 kg, and payload center 100 mm
anterior and 75 mm above the hinge. Payload mass is the remainder after all
modeled lower-body, spine, and added hardware masses. Current estimate:
4.307078 kg lower body, 0.12 kg spine, and 70.572922 kg payload. Read the material
and hardware assumptions before treating this as a measured physical mass.

The user selected adjustable routes from the existing robot-body OpenSim model.
Use `ah001_with_opensim_routes.xml` and editable `opensim_routes.json`. All 92
routes are transferred; 20 lose joint leverage without foot articulation.
Bilateral add_mag2, bilateral grac, and left tfl reverse their neutral sagittal
action and need attachment review. Moving and conditional points are frozen
at neutral; original definitions remain in `source_route_inventory.json`.
Implement their intended dynamic behavior before claiming dynamic equivalence.
The 65-pose constrained motion audit must continue to pass after adjustments.

`trunk_routes.json` retains the earlier six-route template separately. Its
static demand is about 1394.6 N per erector; the OpenSim-derived routes require
about 1354.8 N each to share 69.23 N·m. These are geometric force requirements.
Verify BPA force/pressure/stroke data and neural mapping before adding active
forces. All 92 current interfaces deliberately have zero gain.
Preserve the distinction between retargeted OpenSim routes and measured robot
routes. The source assembly has no separate feet or ankle joints; do not invent
feet without labeling their geometry and mass.

Then couple the specified SNS controller, verify action signs dynamically,
inspect moment arms over the motion range, check force and stroke limits, and
complete supported stance before stepping and support weaning. Apply the
trajectory-parity speed discipline above. The existing two-second test uses
explicit gravity compensation and is a mechanical diagnostic only.

## Existing human and robot families

The user confirmed task 3 means the existing `s3k`, `w2lvar`, and `syn6` human
and robot families. Begin this campaign after tasks 1, 2, and 5 meet their
completion criteria, as requested. Read `spinal/DESIGN.md` before modifying
runner, builder, or parameters. Preserve the default census and the documented
bit-exact anchor evaluation; do not silently substitute older baselines.

Run fresh-process baselines for every family and model pairing, audit saved
knob loading and conditional topology, compare per-term metrics and trajectories,
and perform dynamic sign, sensitivity, basin, and support-weaning checks before
claiming improvement. Reuse existing studies safely; never overwrite a study
with a different search-space distribution or schema. Save every tuned knob
and implement its loader branch. Report unsupported stepping and unsupported
walking separately from harness-supported motion.
