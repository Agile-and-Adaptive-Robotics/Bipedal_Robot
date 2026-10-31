# Contact-aware initialization and repaired-model ground alignment

easteregg2, 30 October 2026. New isolated functions only. The existing
`standing_solver.py` and all v1/v2/v3 force profiles were preserved.
Evidence is under local ignored `runs/contact_standing_v1/`, `v3/`, and `v4/`.

## Corrected diagnosis

The current runner's solve copy does **not** disable ground or springs:
`model_solve=apply_harness(...,kxy=1500)` defaults to ground enabled.
The reconstructed model has six predefined ground pairs and disableflags=0.
The old corrected initializer instead freezes zero-activation constraint
forces and fits seventeen internal coordinates. That is a different limitation.

The supplied failed standing recording reproduces semimembranosus 0.8536,
biceps-long-head 0.7977, and gluteus-maximus-3 0.5326 initializer activations.
Its initial feet penetrate the actual ground plane by approximately:

- Right hindfoot 34.9 mm, forefoot 56.2 mm.
- Left hindfoot 36.4 mm, forefoot 57.6 mm.

Actual initial MuJoCo foot reactions are about 6.6–7.0 BW, not quiet standing.
An optimization allowed to invent approximately 1 BW of static reactions can
have a small algebraic residual while failing this actual contact law badly.
The existing −20 mm height offset is therefore not a validated ground-placement
choice. This observation concerns the initial pose; a force-supported model
can subsequently rise/settle, so it does not prove constant late penetration.

## New root/contact-aware prototype

`contact_standing_solver.py` provides:

```python
solve_contact_standing_activations(
    model, data, Fmax=None, preload=False,
    bilateral=True, support_fraction=None,
    contact_gap_tolerance_m=0.003, prior_weight=0.001,
    return_diagnostics=False,
)
```

Returns 92 activations, or `(activations, diagnostics)`. It builds the full
equality reduction including the six root coordinates (23 independent scalar
coordinates total), with 46 bilateral activation parameters and up to four
nonnegative hindfoot/forefoot normal-force variables. Existing joint springs,
native muscle passive bias, actual pose-dependent active capacities, gravity,
and already-applied external forces enter equilibrium. Optional support_fraction
replaces the scratch pelvis upward load; it does not apply a new force to the
caller or actual simulation. Caller model/data remain unchanged.

Optimized ground reactions are **never injected** into dynamics. The routine
also evaluates actual MuJoCo forward contact forces at the same pose. Gap
eligibility within 3 mm is only a candidate-contact convention, not proof that
the actual contact law generates the optimized load. Unknown tangential forces
are omitted. Root force rows and torque rows use declared BW/length scaling;
the prior is weaker than in the old initializer. Comparisons therefore test a
whole formulation variant, not the isolated effect of adding root rows.

On the original deeply penetrating pose the new linear residual is small,
but actual forward force residual remains large; two-second open-loop holds
still lean substantially. Contact optimization alone does not fix the setup.

## Tested initial-pose/reference convention

Use `ground_alignment.py`, explicitly opt-in:

```python
from physics_fidelity.ground_alignment import align_via_repaired_probe

metadata = align_via_repaired_probe(
    raw_model, raw_data, original_apply_harness,
    flatten_ankles=True, target_gap_m=-0.002,
)
```

Call after selecting native hip/knee/etc coordinates and before runner captures
the pose and constructs its solve/live support rigs. Pass the captured ORIGINAL
`runner.apply_harness`, not a recursive wrapper.

The helper constructs a repaired probe, seeds the current native-coordinate
pose, finds one common ankle angle that equalizes mean hindfoot/forefoot
minimum heights, shifts the pelvis along the ground normal, and copies the
aligned native coordinate values back. Final actual rig references are then
constructed from this aligned pose. No runtime qpos/qvel clipping occurs.

**Do not align on the raw conversion geometry.** Runner repairs flip hip axes;
the same numbers produce different geometry before and after repair. A first
raw-model attempt produced a misleading approximately 20° ankle/73 mm shift.
Those `runs/contact_standing_v2/` results are not a valid alignment convention.

For the supplied recording, repaired-model alignment gives:

| Choice | Common native ankle | Pelvis raise along ground normal |
| --- | --- | --- |
| Height only, minimum gap −1 mm | Unchanged 0° | 56.634 mm |
| Flatten feet, minimum gap −1 mm | +5.62137° | 38.90395 mm |
| Flatten feet, minimum gap −2 mm | +5.62137° | 37.90395 mm |

Height-only adjustment leaves heels 20–22 mm above ground: it preserves the
tiptoe posture. The −2 mm flat result has approximately right regional gaps
−0.558/−0.563 mm and left −2.000/−1.995 mm. The residual 1.44 mm R/L height
difference comes from the borrowed pelvis roll/list; the helper does not
silently level the pelvis or force equal leg loads. Recompute alignment for
each pose/model; do not hard-code these shifts as universal dimensions.
"Flatten" here means sagittal hindfoot/forefoot minimum-height alignment,
not a claim that the entire nonplanar mesh sole lies in the ground plane.

Existing `runner.patch_xml` also sets ankle ligament-surrogate springref from
the initial pose. Rebuilding the rig after changing ankle angle consequently
rebases that rest angle too, approximately a 0.98 N·m spring-torque difference
at stiffness 10 N·m/rad. This is existing pipeline behavior, but must be disclosed
as part of the pose variant. To isolate a pure pose/vertical-reference change,
the old ankle spring reference must be restored explicitly after construction.
No such additional comparison is claimed here.

## Actual eight-second holding checks

`contact_standing_audit.py --aligned-long` tests constant activation/control
without an SNS. It reconstructs the supplied v3 model, joint restraints, and
support settings, and tests both the old initializer and the new prototype.
The knee upper stop is 0°, with 2° pre-neutral margin, solref [0.004,1], and
solimp [0.999,0.9999,0.001,0.5,2]. Every integrated knee state is checked; no
trajectory clipping is performed.

For flat alignment at −2 mm and the root/contact initializer:

| Holding condition | Spring-supported standing | Constant 70% upward support |
| --- | --- | --- |
| Actual initial foot normal force | 0.786 BW | 0.329 BW |
| Mean later foot force | 1.011 BW | 0.299 BW |
| Maximum absolute pelvis tilt | 2.68° | 2.18° |
| Maximum COM displacement from initial COM | 12.1 mm | 6.15 mm |
| Runner-style foot-referenced maximum sway | 31.8 mm | 31.9 mm |
| Maximum right/left knee extension | −1.846° / −1.885° | −1.950° / −1.967° |

Both complete eight seconds with no positive knee angle. The stops, harness,
joint restraints, and constant activation remain assistance. The torso's
absolute lean is a separate metric from the runner's pelvis-tilt gate.
Holding force means use 1–8 s; geometric/extension/tilt extrema cover all saved
holding steps. Parent ground-curriculum support uses its separate 5–13 s window.
These are **holding diagnostics, not SNS curriculum passes**. The old internal
initializer also holds much better after geometric alignment, especially at
70% support; the contact-aware prototype is not universally necessary/superior.

## Actual SNS transfer: important negative result

Parent-session results were checked separately under
`hunt_training/runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/`:

- `contact_feasible_start_v1/standing.result.json`, aligned −1 mm with the
  previous coupled initializer, still fails: pelvis tilt 12.69° and sway
  93.6 mm; independently measured scored foot load sums to about 0.423 BW.
- `contact_aware_initializer_v1/standing_gap_0p001.result.json` and the −2 mm
  counterpart also fail: sway about 80 mm and independently measured scored
  foot load about 0.13 BW. The attractive fixed-activation holding result does
  not transfer through the present POST/MN/reflex interface.
- `contact_aware_initializer_v1/ground.result.json` passes contact timing but
  fails movement phase: knee phase error 0.15266 versus the 0.15 bound.
- `contact_feasible_start_v1/ground_guide.result.json`, using the previous
  coupled initializer, passes the evaluator's current full ground checks.
  Its independently measured support in the 5–13 s scoring window is
  0.63726 BW, and contact phase error is 0.06293. It remains a guided,
  assisted ground demonstration; parent owns fresh repeats and support weaning.

Therefore repaired-model pose alignment is justified by geometry, startup-force,
holding, and actual SNS evidence. The new contact-aware initializer remains
experimental and must not replace the verified guided-ground winner by default.
Further work must inspect how its desired activations are expressed by the
POST/MN/reflex network, rather than treating an algebraic or fixed-control
equilibrium as a neural-controller success. No causal walking improvement is
attributed to the new initializer alone.
