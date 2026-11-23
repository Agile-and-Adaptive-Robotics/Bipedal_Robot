# Pressure-driven robotbody experiments

The model uses the 27 bilateral routes of the mirrored Gait2392 robotbody.
Its spinal network has one motor output per route. Each output sets a valve
target between 0 and 620 kPa. The baseline from Table C13 uses 102 sleeves and
the selected design from Table C14 uses 138. Several sleeves can share a valve
and pull on the same output cable. The route geometry, survey sleeve dimensions,
actuator counts, tackle ratios, and cable and bracket compliance come from the
same [sizing transfer](../robotbody_survey.py) as the fixed-pelvis experiments.

At 50.47% body-weight lift in stage 5, the selected design sustained
alternating loading and clearance in four independent 16 s trials. Every run
produced five clear swings on the right and six on the left. Nominal and
perturbed starts recorded five right and six left loading onsets. Halving the
timestep or reducing the contact time constant from 20 to 10 ms gave five
onsets on each side. In the nominal trial the feet carried 29.9%/19.6% of
body weight on average, the knees stayed flexed about 100–118°, and the
center of mass advanced 1.1 cm. Knee excursions of 14°/18° fell below the
20° criterion. [Reduced-lift evidence](selected_reduced_lift_evidence.json).

In that trial right vastus tension reached 48 N during right stance. The
chamber reached 131 kPa, whereas the same route lengths provided a median
modeled capacity of 1.94 kN at 620 kPa. Right medial gastrocnemius remained
slack. [Pressure, force, and capacity verification](selected_reduced_lift_force_check.json).

The separate [selected ground trial](selected_ground.json) reproduced seven loading
onsets per foot with 79% body-weight lift in four independent 16 s simulations.
The feet also cleared while the other foot carried load. The nominal clearance
peaks were 6.2 cm on the right and 3.2 cm on the left. The knees stayed deeply
flexed, at 92–119°, and flexed twice per loading cycle. The center of mass
advanced 1.9 cm under horizontal guidance. These results demonstrate supported
load transfer and clearance in a crouch. The knee timing failed the stepping
criterion. [Per-foot and per-route evidence](selected_ground_evidence.json).

In the nominal 5–13 s scoring window, when the right foot carried more than
5% of body weight, the right vastus chamber reached 91 kPa and its tension
reached 36 N. Evaluating the documented design at 620 kPa at those same route
lengths gives a median capacity of 2.08 kN. A weak stance response in this
trial therefore reflects, in part, how little the circuit recruited that
available extensor. The medial gastrocnemius was slack over the same interval.

| Frozen condition | Right/left loading onsets | Right/left clear swings | Load timing and duty | Knee timing |
| --- | ---: | ---: | --- | --- |
| Nominal, 0.5 ms physics | 7/7 | 7/7 | Passed | Failed |
| Knee starting angles perturbed ±1.25° | 7/7 | 7/5 | Left foot stayed loaded for 87% of the window | Failed |
| 0.25 ms physics | 7/7 | 7/7 | Passed | Failed |
| 10 ms contact time constant | 7/7 | 7/7 | Passed | Failed |

## Installations and control

The [kinematic atlas](locomotion_atlas.json) gives the route lengths observed
in 4,360 poses of the converted robotbody. The default installation reserves
10 mm of extra output travel beyond the longest sampled length on either side
of each mirrored pair. It retains each source sleeve, tendon, and fitting.
The required route offsets correspond to as much as 0.61 m of additional
fixed proximal path, or 0.17 m (baseline) and 0.74 m (selected) of extra cable.
They would require a packaging study on a physical robot. Six baseline and ten
selected routes are slack at the neutral pose in this installation. The
additional suspension trials also evaluate specified neutral trunk contractions.
Each saved installation identifies the mounting datum used on each route.

The 10 and 20 mm coefficients derive from force measurements, and the 40 mm
coefficients and force scale derive from manufacturer data. Pressure fills
with a 40 ms and vents with a 60 ms time constant. The force surface is
quasistatic and uses the same coefficients on both branches. Each bundle
pulls through a compliant cable. At maximum sleeve shortening the remaining
cable stretch sets the transmitted force. A route that would require extension
beyond the characterized domain at its present pressure stops the trial.

The controller updates every 2 ms. Its rhythm generator drives pattern-forming
neurons and the named motor pools. The first 180 ms of drive also supplies a
finite asymmetric start pulse. The suspended deafferented test removes
peripheral feedback and interleg connections. The afferented condition enables
length, velocity, and tension feedback. Ground trials add foot load and regional
contact. The last condition can recruit additional contact inputs to the
pattern-forming layer. A static equilibrium fit supplies constant initial
posture currents. Certain refinement trials calibrate those currents with the
pose held for one second. The model then fills its chambers for 0.6 s before
release. During the physical trial the neural outputs, chamber states, joint
positions, and contact forces evolve together.

## Fixtures and acceptance

These five conditions are the dissertation's physical tests. Its separate
four-stage neural training study uses human reference muscle recruitment. The
pneumatic trials tune the controller against joint motion, foot loading,
clearance, and measured support in the robot model.

The torso suspension in stages 1 and 2 applies nominal full body-weight lift
and records its six-component wrench. There is no floor contact in suspension.
On the ground, the pelvis has a specified vertical lift, forward/lateral spring
stiffnesses of 1,500/2,000 N/m, and an angular stiffness of 400 N·m/rad.
Separate passive hip guides constrain adduction and axial rotation, with
450/400 N·m/rad stiffness. Ground trials also guide the three lumbar rotations
at 300 N·m/rad. Foot contacts, passive end stops, and gravity remain active.
The force and moment histories quantify the pelvic support supplied in each
trial. A separate no-rig condition sets pelvic lift, applied pelvic guide
forces and moments, and the added hip and lumbar guides to zero. The passive
joint-limit model and its small ligament surrogates remain in place.

| Stage | Physical experiment | Acceptance |
| --- | --- | --- |
| 1 | Legs move with torso suspension and peripheral feedback disabled | Bilateral knee rhythm and timing, hip and knee excursion, frontal ordering, six active trunk routes, and lumbar excursion at most 15°. |
| 2 | Repeat in suspension with sensory return | The same mechanical criteria, now with bilateral connections and peripheral input. |
| 3 | Standing under 35% vertical lift, with a 35 N forward pelvic disturbance for 150 ms | Complete 8 s. Pelvic tilt below 15°, horizontal center-of-mass excursion below 60 mm, a mean combined foot load above half body weight, and each foot carrying 30–70% of that load. |
| 4 | Ground trials with pelvic support | Complete 16 s. Pass the movement and frontal criteria, alternating foot-load criteria, and at least three genuine clear swings on each side while the opposite foot carries load. |
| 5 | Ground trials with lower lift and additional contact feedback | Repeat stage 4 under the specified reduction in vertical lift. |

The motion assessment uses the [original hip/knee and timing criteria](../spinal/functional_gait2392/curriculum_metrics.py).
The [force check](../spinal/hunt_training/ground_contact_check.py) uses 5–13 s
of data, an unloading threshold of 2.5% and a loading threshold of 5% of body
weight, and minimum unloaded/loaded intervals of 80/100 ms. The position check
uses the lowest point of each talus, heel, or toe mesh, at 20 ms intervals.
One clear swing requires at least 80 ms above 1 mm, a peak of at least 10 mm,
and simultaneous contralateral load above 5% of body weight. Later stages
still run as diagnostic trials when an earlier prerequisite fails.

## Repeating the search

To replay any of the ten design and stage combinations individually, pass
`--config Code/MuJoCo_SNS/robotbody_curriculum_reference/study.json` together
with `--design`, `--stage`, the atlas, `--dt .0005`, and an output prefix.
Use `--rank 0` or `--rank 1` for the two frozen finalists in that record.
Use `--fit local` to replay the final targeted standing or supported-stepping
configuration. Use `--fit coarse --dt .001` to replay a first-round winner.

Use fresh output folders. The first search uses deterministic Sobol starts and
a bounded local refinement. The second tries narrower parameter changes,
transfers gains between designs, and repeats the two finalists per design in
fresh processes. The held-out condition starts the right and left knees
1.25° on opposite sides of their nominal initial angles. Separate checks halve
the physics timestep or reduce the contact solver time constant from 20 to 10 ms.
The [exploratory ground seeds](ground_search_seeds.json) supply starting gains
only. Every reported physical outcome uses the shortening-limited pressure law.

The ground check then adjusts hip, knee, and ankle recruitment around one
well-timed and one kinematically promising candidate per design. It compares
changed lift, external guidance, and actuator design with the exact 54 named
motor-current presets from the original trial. It includes a second held-out
knee asymmetry, 1.75° in the opposite direction. A nominal full gate pass
also triggers eight additional trials with neural gains varied by up to 1%.

The same search tests a connection from the existing ipsilateral load
interneuron to the vastus motor pool because the model developed little vastus
force in loaded crouch despite available actuator capacity. Its conductance is
recorded as `load_to_vastus_uS`. Optional postural feedback drives the spinal
ankle, hip-abductor, and trunk balance inputs from center-of-mass and torso
position and velocity. An additional standing comparison fits this feedback
at two initial knee poses (±0.5° asymmetry), then retests the same gain at
the nominal pose and the two held-out initial poses. These interventions act
through motor neurons and the same pressure and tension equations.

```powershell
$python = 'D:/Anaconda/envs/myo/python.exe'
$models = 'Code/MuJoCo_SNS'
$atlas = "$models/robotbody_curriculum_reference/locomotion_atlas.json"
& $python "$models/tune_robotbody_curriculum.py" --output D:/temp/robotbody_coarse --atlas $atlas --workers 4 --trials 24 --dt .001
& $python "$models/refine_robotbody_curriculum.py" --initial D:/temp/robotbody_coarse --output D:/temp/robotbody_refined --atlas $atlas --workers 4 --trials 22 --dt .0005
& $python "$models/check_robotbody_ground.py" --refined D:/temp/robotbody_refined --output D:/temp/robotbody_ground --atlas $atlas --workers 4 --trials 18
```

## Results across the curriculum

The searches ran 520 training evaluations and 152 independent validation or transfer evaluations. The table gives finite full-duration training runs and nominal passes of every criterion for a stage. Each physical condition was attempted on both designs.

| Stage | Baseline trials | Completed | Nominal passes | Selected trials | Completed | Nominal passes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1. Suspended, deafferented | 46 | 7 | 0 | 46 | 7 | 0 |
| 2. Suspended, afferented | 46 | 17 | 0 | 46 | 26 | 0 |
| 3. Disturbed standing | 58 | 38 | 4 | 58 | 32 | 11 |
| 4. Pelvic-supported stepping | 64 | 9 | 0 | 64 | 51 | 1 |
| 5. Lower lift | 46 | 1 | 0 | 46 | 17 | 0 |

The final local selected-design ground configuration passed its nominal 16 s trial at 78.5% lift. It also passed 4 of 8 one-percent gain perturbations. The opposite initial knee offsets, halved timestep, and stiffer contact all retained six right and seven left loading onsets but failed at least one knee-motion criterion. The knees remained deeply flexed. Neither design passed all five physical stages across the independent conditions.

Baseline fell at 1.70 s with zero vertical lift and the other guides retained. Removing the added guides as well exceeded the calibrated sleeve-extension domain at 0.056 s (intobl_r, intobl_l). Selected fell at 2.66 s with zero vertical lift and the other guides retained. Removing the added guides as well exceeded the calibrated sleeve-extension domain at 0.346 s (ercspn_l). All of these support comparisons held the 54 named postural motor currents and the circuit gains fixed. [Complete results and runnable finalist settings](study.json).
