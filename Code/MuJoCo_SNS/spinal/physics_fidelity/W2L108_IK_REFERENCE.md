# Newer W2L108 supported stance reference

Machine: easteregg2. This audit concerns the newer canonical 108-neuron W2L
model, not Gait2392 or AH001. Source is the exact
`virtual_walker_20261031/runs/W2L_VirtualWalker_Newer/canonical108_ground_baseline_v1/model.xml`
and its saved trajectory. No source geometry, joint limits, masses, routes,
passive springs, runner, or canonical neural connections were changed.

## Pose to trial

Durable output directory:
`physics_fidelity/runs/canonical_12II/ik_reference/`.
Load `flat_existing_rails.npz`, key `qpos`, as the complete 17-element initial
position vector. Use the exact corresponding source XML; unrelated models may
have different joint ordering. Set initial velocities to zero.

Root position is `[-3.454, 0.9408845688307746, 0]` m, quaternion `[1,0,0,0]`.
This lowers the root 42.000003 mm from the XML assembly pose and puts its
existing wings on the existing rail tops. It does not lower the rails.
No guide/rail relocation was needed. Root pitch is zero.

| Coordinate | Left native q (degrees) | Right native q (degrees) |
| --- | ---: | ---: |
| Hip | 24.8929964 | -5.1057105 |
| Knee | 3.9431478 | -41.0353262 |
| Ankle | -0.6857989 | -0.6655573 |
| Toe | -1.5185417 | -1.5185203 |

After the source assembly offsets, thigh pitch is 14.893/14.894 degrees and
signed knee flexion is 13.943/13.965 degrees. Positive native ankle q is
dorsiflexion, with zero anatomical offset. Both ankles remain within the
unchanged approximately -12 to +3 degree limits. Minimum source knee flexion
allowed by the existing ranges is approximately 3 degrees; no hyperextension
was introduced.

The main feet are approximately horizontal: pitch +0.264 degrees and sole
height span 0.559/0.557 mm. Main-foot minimum belt gaps are -0.021/-0.020 mm;
toe minimum gaps are -0.530/-0.529 mm. These small penetrations are explicit
soft-contact initialization, not altered geometry. All selected sole-face
corners are within 0.538 mm of the belt top. Wing/rail gaps are -3.1 nm,
numerical solver tolerance. The body is not suspended above the rails.

Exactly coplanar main-foot and toe faces cannot be obtained with these fixed
attachments. An exact main-foot-horizontal/toe-edge-touch solution instead
rotates toes about +25 degrees, lifting their fronts about 18 mm and creating
large passive spring torques. That diagnostic is saved as
`exact_mainsole_existing_rails.npz` but is **not the recommended reference**.

## Inverse recruitment prior

`STATIC_RECRUITMENT.json` uses the actual pulling actuator moment matrix,
native passive springs, six root coordinates, eight leg/toe coordinates,
nonnegative regional foot/toe normal reactions, and actual overlapping
wing/rail support footprints. It includes unavoidable positive parallel Hill
tension as a lower force bound. Reactions are optimized unknowns, not measured
simulation forces, and are never injected into the simulation.

The projection allocates 0.20974 BW to feet and 0.79026 BW to rails. Root
equilibrium residual is below 1.3e-9 in its respective N/Nm coordinates;
hip/knee/ankle residuals are below 6.2e-7 Nm. The two passive toes retain
0.329 Nm residual apiece. Therefore this is not an exact static equilibrium
of every coordinate. Toe-region reactions were included; allowing negative
normal forces to conceal that residual would be unphysical.

Original Hill steady-state recruitment uses

`F = clip((Kpe*(L-L0) + TL*(Bsig/(1+exp(C*(A-V)))+D))/(1+Kpe/Kse), 0, Fmax)`

with `TL=max(0,1-((L-L0)/width)^2)` and measured route lengths. The JSON records
the inverted stimulus voltage, original MN binding, adapter polynomial, and
candidate MN voltage. These are a static inverse prior, **not a fitted MN
trajectory or a validated network gain set**. At a passive lower bound,
voltage inversion is ill-conditioned; do not treat an arbitrarily very
negative value as a precise biological target. At zero clipped tension, the
inversion describes a voltage bound rather than a unique solution.

The native binding itself matters: `R Hip MN ext` drives `hip_R_flx`, whereas
`R Hip MN flx` drives `hip_R_ext`; ankle MN ext drives the physical
dorsiflexor and ankle MN flx drives the plantarflexor. Preserve measured
functional signs when using this prior.

## Actual forward check

The 2 s test applied the solved twelve constant muscle tensions through their
original tendon actuators in actual MuJoCo, with all native contacts and
passive springs enabled and a stationary belt. It was not an SNS run or a
dynamic Hill feedback run. Actual reactions settle near 0.204 BW feet and
0.796 BW rails. At 1.982 s, root Y was 0.940883876 m and pitch -0.000338
degrees. The knees stayed flexed at every one of the 1,000 physics states:
left 13.957–14.791 degrees, right 13.978–14.816 degrees. No clipping was used.

The toes relaxed from -1.519 to +0.686 degrees and lifted about 6.1 mm. Main
feet retained actual contact, with final minimum gaps about -0.086/-0.081 mm.
Their final pitch was +3.36/+3.34 degrees, with approximately 7 mm sole height
span. Thus this probe maintained supported heel/main-foot contact, **not the
initial flat-foot orientation or persistent toe contact**. Optimized static
normal reactions need not match the reactions selected by compliant contact.
The raw hold and settled
pose are `FIXED_FORCE_HOLD.json/.npz` and
`fixed_force_settled_qpos.npz`. The complete 2 s fixed-force physical loop
took 0.219 s, approximately 9.1 times real time; this excludes neural
computation and does not benchmark the full controller.

Next: initialize the authored SNS from this actual pose, measure heel and toe
contact separately, and compare its MN outputs with the inverse prior before
gain search. Passing this short assisted hold is not evidence of stepping,
reduced assistance, or harness-free walking.

Generator: `physics_fidelity/w2l_ik_reference.py`. The script uses compiled
MuJoCo kinematics/dynamics, bounded SciPy least squares, preallocated Jacobian
buffers, and direct original Hill-law inversion. Full pose, raw measurements,
toe torque scan, and static projection are retained in the ignored run
directory.

## Independent anatomical action check

`w2l_fk_action_audit.py` performed central finite differences of 1e-5 radians
at the saved near-flat pose. Hip flexion was defined independently by the
knee anchor moving anterior to its hip anchor, knee flexion by the relative
thigh/shank directions, ankle dorsiflexion by foot pitch relative to shank
orientation, and toe extension by toe pitch relative to the main foot.
No additional IK or hold was run for this check.

Positive native q increases the named anatomical action on **both sides**;
all eight diagonal derivatives are approximately +1. Assembly offsets are:

| Joint | Left home offset | Right home offset | Positive native q |
| --- | ---: | ---: | --- |
| Hip | -10 degrees | +20 degrees | Hip flexion; knee moves anterior |
| Knee | +10 degrees | +55 degrees | Knee flexion |
| Ankle | 0 degrees | 0 degrees | Dorsiflexion |
| Toe | +2 degrees | +2 degrees | Toe extension |

For example, left native hip 24.893 degrees is anatomical hip flexion
14.893 degrees; right native hip -5.106 degrees is anatomical flexion
14.894 degrees. The offsets describe the source assembly rather than
different angle signs for the two sides.

Signed pulling moment arms at this pose and actual adapter bindings are:

| Side(s) | MN | Actual muscle | Anatomical action | Moment arm |
| --- | --- | --- | --- | ---: |
| Left | Hip MN flx | hip_L_flx | Hip flexion | +0.09649 m |
| Right | Hip MN flx | hip_R_ext | Hip flexion | +0.09649 m |
| Left | Hip MN ext | hip_L_ext | Hip extension | -0.07941 m |
| Right | Hip MN ext | hip_R_flx | Hip extension | -0.07939 m |
| Both | Knee MN flx | knee_*_flx | Knee flexion | +0.05393/+0.05396 m |
| Both | Knee MN ext | knee_*_ext | Knee extension | -0.05929 m |
| Both | Ank MN flx | ankle_*_ext | Plantarflexion | -0.03978/-0.03977 m |
| Both | Ank MN ext | ankle_*_flx | Dorsiflexion | +0.07616/+0.07617 m |

Thus the right hip **muscle names** are reversed relative to their action,
but the original hip **MN labels and bindings are functionally correct** on
both sides. The ankle MN labels require the explicit dorsiflexor/
plantarflexor interpretation above. All twelve MN-to-muscle adapter
polynomials are identity mappings in volts. There are no active toe MN
bindings in this model. Full joint-position derivatives and bindings are
saved in `FK_ACTION_AUDIT.json` in the same durable run directory.
