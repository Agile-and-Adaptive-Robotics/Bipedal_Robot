# s3k, w2lvar, and syn6: derivation, wiring, and improvement attempts

Session 1, 28 October 2026, **easteregg2**. Repository:
`D:\GitHub\Bipedal_Robot`. Python: `D:\Anaconda\envs\myo\python.exe`.

## Findings

These are literature-informed engineering controllers, not three reproductions
of a single published circuit. They share an RG/PF hierarchy and graded SNS
neurons, but differ in their PF representation, sensory pathways, and crossed
connections. The current CAD robot's successful defense demonstration uses
an engineering MN target-current clamp. That clamp cancels incoming PF/reflex
MN current, so its motion is not evidence that the native PF muscle patterns
work on the robot. This audit separately tests the conventional Gait2392 runner,
where those neural pathways remain active.

The recorded s3k anchor reproduces exactly, including its **frozen left leg**.
Its numerical stability must not be reported as bilateral walking. No wiring
change from the first bounded human-model screen improves its family's
baseline kinematic score. Existing winners and production builders are retained.

## What was actually built

The following are fresh default-family inventories. The 54-actuator column
uses the existing CAD preset's contact gains. Human inventories use the 92-name
basis/prune union; the files save the exact name list. Configured trial counts
differ because conditional pathways are enabled.

| Family | PF excitatory cells per leg | Default human neurons / inputs / builder edges | CAD neurons / inputs / builder edges |
|---|---|---:|---:|
| s3k | E1, E2, F1, F2; shared E/F inhibition | 410 / 376 / 1186 | 262 / 230 / 672 |
| w2lvar | HIP-E/F and shared KNEE/ANKLE-E/F; inhibition per pair | 888 / 382 / 7866 | 544 / 230 / 3292 |
| syn6 | Six NMF channels S1–S6; shared family inhibition | 794 / 382 / 3033 | 488 / 230 / 1481 |

w2lvar's four ankle watch indices alias knee cells; they are not four extra
neurons. Its CAD census is 544, not the inherited runner's 548 index count.
Builder-edge counts are not always unique active conductances: default
w2lvar has 7864 nonzero matrix entries, and CAD s3k has 650. The inventory
records duplicate endpoints and zero edges rather than silently treating all
builder entries as independent synapses.

### s3k

The RG has persistent-sodium half-centers, with RG-E -> InE -> RG-F and the
mirror inhibitory chain. Weak direct RG-to-RG excitation is optional. Four
graded PF cells have staggered time constants; these are phase-window cells,
not a separate antagonistic pair at every joint. E cells excite PF_IN_E,
which inhibits F cells, and vice versa. MN drive is a functional-group table,
including reduced secondary-group weights for some biarticular muscles.

Ia excites the homonymous MN; reciprocal inhibition uses IaIN when enabled.
II and Ib routing depend on `full_rules`; the reduced and full-rules builds
must not be conflated. Enabled central pathways project extensor Ib to RG-E,
InE, and E PF cells, and flexor Ia/II to RG-F, InF, and F PF cells.
The reference anchor explicitly uses `full_rules=0` and has 602 configured
neurons, rather than the 410-neuron all-default census. The configured
kinematic anchor is -161.56754173676563.

The default PF table supplies no phase drive to twelve human adductors, or
to the CAD robot's bilateral gracilis routes. These MNs may still receive
posture/reflex inputs. Primary-group-only antagonist pairing is an engineering
classification, not a complete map of each actuator's pose-dependent 3D action.
In particular, the map suppresses secondary sartorius, TFL, and gracilis actions
for historical tuning reasons and has no separate internal/external rotation
groups. CAD route mechanics must be checked using moment arms, not those labels.

### w2lvar

Its strongest direct precedent is the hip pair plus shared knee/ankle pair
in [Deng et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6477610/).
That paper used six muscles per rat hindlimb, sagittal mechanics, and body
support. The 54/92-muscle extension, 3D CAD body, trunk, and V-class crossed
circuitry are additional modeling choices. Its PF mapping uses a fitted
joint-layer table where available, with anatomical and biarticular fallbacks.
There is no independent PF pair exclusively for each biarticular muscle.

Both pairs have laminated reciprocal inhibition. Per-muscle Ia/II/Ib central
feedback reaches the appropriate F/E PF cells, RG, and RG inhibitory interneuron,
with normalization by each side's functional-muscle count. Ia, II, Ib, antagonist
interneurons, and Renshaw circuits are more extensive than the rat PF template.
Renshaw cells exist for every actuator in this variant, even when their output
gain is zero. The six trunk actuators have **no PF input by design**; the rat
template supplied no trunk PF circuit. They rely on posture/balance/reflex drive.

The V2a/V0V/INI, V0D, and V3-E pathways draw on Shevtsova-style rules, but use
calibrated gains and a lumbar reduction. `w2l_comm` is a **positive/zero topology
switch**, not a proportional gain: changing 1 to 0.5 changes nothing. The
conductance-reduction probe therefore changes the actual nine crossed-pathway
gain keys. This prevents another nominal search over an inactive parameter.

The suspected fitted-zero fallback issue was checked against the current
weight file: **no explicit fitted zeros are currently reintroduced**. It is a
potential future file-semantics concern, not a demonstrated present defect.

### syn6

The six channels come from NMF of this project's back-solved activations from
subject01_walk1 IK and measured GRF. They are not six anatomically identified
neuron populations or a multi-subject EMG validation. Channel families are
F/E/E/F/F/E on both sides. S5 is mixed, receives both RG drives, and is excluded
from the default reciprocal PF inhibition to avoid a measured latch. Therefore
the six channels are **not three conventional antagonistic half-center pairs**.

The basis covers every current CAD actuator. It omits only the six deliberately
pruned human short-rotator actuators. The initial concern about missing CAD
coverage was disproved. Synergy gains use the analytical conductance mapping
from [Szczecinski et al., 2017](https://www.frontiersin.org/journals/neurorobotics/articles/10.3389/fnbot.2017.00037/full).
The mapping does not establish exact linear NMF reconstruction when several
conductance-based inputs simultaneously drive an MN; synaptic shunting and
activation clipping remain part of the dynamics.

Ia/II/Ib reflex interneurons and antagonist inhibition are unconditional, but
per-muscle sensory projections into RG/PF are **absent at defaults**.
`ib_rge` enables a group load path, while the newer `ib_e_central` enables
per-muscle extensor force input to RG-E/InE and relevant E synergy channels.
The older curriculum does not load that newer key; the experimental loader
explicitly applies it and records the effective dictionary. `iain_pf_gate`
and `pf_mixed_recip` similarly enable previously absent PF gating/lamination.
Default syn6 has no Renshaw cells. The optional shared brainstem drive is off
because it previously caused tonic latching.

### Connection does not imply mechanical authority

The separate `robot_torque_coverage.json` projects the actual dynamic tendon
Jacobians onto the knee-loop tangent at five recorded poses. **Twenty of the
54 CAD actuators have no resolved free-joint torque at all five poses**:
the ten monoarticular distal routes per side, including soleus, tibialis,
peroneals, and digital flexors/extensors. Their distal joints are welded.
Gastrocnemius retains knee leverage. Thus a fully populated PF-to-MN map is
not a fully actuated 54-muscle locomotor model. Internal tension and Ib
feedback can still exist in a welded route, so force feedback there is not
automatically evidence of useful stance-load support. This is another reason
to adapt the sensory/group map to the robot's actual free joints.

## Literature similarity and deviations

The Deng 2023 biarticular follow-up and Hunt's 2019 INRC training slides were
reviewed after this baseline audit. See
`reports_20261028_functional_pf/literature_training_review/REVIEW.md` for the
joint-PF integration architecture, active/passive perturbation distinction,
and four-step functional training method. Physical-environment gain searches
are not equivalent to that method. The inventories and baseline results in
this audit remain the original recorded measurements.

| Source | Supported similarity | Material deviations in these families |
|---|---|---|
| [Rybak et al., 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1890439/) | Separate RG timing and PF distribution, inhibitory interneuron pathways, resetting versus non-resetting perturbations as a test | Graded single-cell surrogates, different PF counts, reduced MN dynamics, fitted engineering gains, and additional postural control. Wiring similarity alone does not reproduce the perturbation results. |
| [Deng et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6477610/) and [Deng et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9776051/) | RG/PF hierarchy, hip and knee/ankle synergy precedent, sensory interaction with mechanics | Human/CAD morphology and many-muscle expansion, additional crossed pathways, new support rig and actuator model, and normalized voltage/conductance calibration. |
| [Deng et al., 2023](https://doi.org/10.1007/978-3-031-39504-8_2) | Joint PF signals feed intermediate BFP/RF/GA populations; active/passive mechanical perturbation comparisons; PF perturbations can shift RG through the closed loop | None of the three baseline families reproduces that full architecture or experiment. The new opt-in human joint-integration variant implements only part of it and remains untuned. |
| [Shinohara et al., 2025 preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC12642437/) | Flexor length/velocity and extensor force feedback can reach RG, inhibitory interneurons, PF, and MNs; C1/V3 coordination | That cat model has seven muscles and two PF populations per side, voltage-dependent sodium kinetics, subthreshold NaP in PF/MNs, nonlinear velocity/length encoding, and a flexor-driven regime. These families retain fixed-tau RG NaP, lack PF/MN NaP, use different sensory encodings, and mix several reflex motifs. |
| [Shevtsova et al., 2026](https://elifesciences.org/articles/107480) | Named V2a/V0V/V0D/V3 and descending-drive motifs | Lumbar-only reduction of a quadruped network, folded interneurons in syn6, recalibrated gains, altered drive, no injury-recovery replication, and fixed tau rather than the paper's voltage-dependent bell-shaped tau. |

The fixed-tau approximation is useful for the existing calibrated family, but
it is not a general literature requirement. The paper-specific Shevtsova SNS
implementation already uses a distinct bell-tau backend; applying a fixed-tau
fast evaluator there would change the equations. The new preallocated backend
rejects gated classes other than `SNS_NumpyFixedTau` rather than silently doing so.

The old assertion that the original Deng circuit had no sensory input to RG/PF
is explicitly **retracted** in `WIRING_RULINGS_20261003.md`, line 395. Earlier
paragraphs in that same file and the September feedback audit retain the wrong
claim. This report follows the corrected ruling, not the stale paragraph.

## Fast calculations and acceptance boundary

All improvement probes retain the original MuJoCo model, physics step, contact,
support, and initialization for the model being tested. The existing exact-fast
backend computes only live conductances with Numba and preallocates its dense
workspaces, retaining NumPy's original synaptic reduction order. This matters
near bifurcations. The older sparse/Numba trial implementation's speed claims
must not be transferred to gated networks or to a different harness loop.

The additional `preallocated_exact_eval.py` removes unnecessary dense one-hot
input/output products and reuses voltage/reduction buffers. Ion-channel
equations and synaptic sum order remain unchanged. Some channel arithmetic
still allocates arrays; this is not a zero-allocation whole-simulation claim.
It is an optional evaluator, not a change to production SNS builders.

The audit checks every default family on both actuator sets, including voltage
and h-gate agreement. Full human-model baseline replays additionally compare
the actual recorded state, activation, and neural arrays. Timings are warmed
measurements on this shared machine; only numerical parity is a strict gate.
The small s3k networks do not consistently benefit from the additional
preallocation, so the earlier exact-fast backend remains appropriate there.

## Experiments and results

The bounded screen uses three frozen baseline configurations and ten targeted
changes, not a new large optimizer search. The s3k probes add opposite RG-IN
inhibition to PF (ratio inspired by Shinohara, calibrated engineering scale),
halve extensor central force input, or halve flexor central inputs. The w2lvar
probes reduce central feedback or actual crossed-pathway conductances. The
syn6 probes enable central force feedback and/or PF reflex gating/mixed-channel
lamination. Each runs in a fresh process and saves the effective gains. There
is no positive locomotor result to promote from that screen.

Separately, CAD trials replace cancellation-based MN tracking with a bounded
postural equilibrium bias. Complete immediate release fails lumbar posture in
all three families. A 5–15 s ramp to 10% natural synaptic-current authority is
also tested; it preserves a substantial engineering clamp and is not a pure
neural controller. The surviving w2lvar partial-release run is tested with a
delayed PF-to-MN deletion after 15 s, so histories can be compared before the
intervention. The geometry and provisional attachment assumptions from the
previous task-5 report remain in force.

<!-- GENERATED_RESULTS -->

### Human-model screen

| Family | Configuration | Kinematic score | R/L cycles | Fell |
|---|---|---:|---:|---|
| s3k | afferent_half | -262.355 | 5/0 | False |
| s3k | baseline | -161.568 | 9/0 | False |
| s3k | ib_half | -165.451 | 9/0 | False |
| s3k | rg_pf_gate | -176.770 | 8/0 | False |
| w2lvar | baseline | -174.848 | 5/7 | False |
| w2lvar | central_half | -194.998 | 3/4 | False |
| w2lvar | central_quarter | -258.544 | 2/2 | False |
| w2lvar | comm_half | -187.132 | 6/8 | False |
| syn6 | baseline | -200.380 | 4/0 | False |
| syn6 | central05 | -216.966 | 4/0 | False |
| syn6 | central15 | -301.042 | 2/3 | False |
| syn6 | central_gate | -222.038 | 8/0 | False |
| syn6 | pf_gate | -226.509 | 9/0 | False |

Higher score is better. These short supported tests do not establish walking. Every targeted change scored below its own frozen baseline. No gait winner was replaced.

### Full-trajectory numerical checks

- s3k versus reference: all saved arrays bit-exact.
- s3k versus preallocated: all saved arrays bit-exact.
- w2lvar versus reference: all saved arrays bit-exact.
- w2lvar versus preallocated: all saved arrays bit-exact.
- syn6 versus reference: all saved arrays bit-exact.
- syn6 versus preallocated: all saved arrays bit-exact.

### Preallocated neural evaluator

| Family/model | Reference / preallocated speed ratio | Previous fast / preallocated ratio | State parity |
|---|---:|---:|---|
| s3k/human92 | 20.73 | 0.93 | True |
| s3k/robot54 | 2.70 | 0.95 | True |
| syn6/human92 | 11.39 | 1.05 | True |
| syn6/robot54 | 9.49 | 1.19 | True |
| w2lvar/human92 | 14.53 | 1.20 | True |
| w2lvar/robot54 | 7.42 | 1.41 | True |

Ratios above one indicate faster preallocated stepping. These are warmed neural-step timings on a shared machine, not whole-trial speedups or controlled throughput benchmarks. For s3k, retain the earlier exact-fast evaluator if it is faster.

### CAD neural-authority probes

| Trial | Simulated duration (s) | Outcome |
|---|---:|---|
| s3k_authority0p1_ramp | 22.242 | pelvis height failure |
| s3k_authority1p0 | 0.126 | lumbar posture failure |
| syn6_authority0p1_ramp | 13.830 | lumbar posture failure |
| syn6_authority1p0 | 0.066 | lumbar posture failure |
| w2lvar_authority0p0_ramp | 25.000 | Completed requested duration |
| w2lvar_authority0p0_ramp_pf_deleted | 25.000 | Completed requested duration |
| w2lvar_authority0p1_ramp | 25.000 | Completed requested duration |
| w2lvar_authority0p1_ramp_pf_deleted | 25.000 | Completed requested duration |
| w2lvar_authority0p1_ramp_preallocated | 25.000 | Completed requested duration |
| w2lvar_authority1p0 | 0.186 | lumbar posture failure |

Delayed PF deletion: identical history before 15 s = **True**; maximum subsequent control difference = **0.95006**. This is a causal effect of PF transmission inside a largely engineering-controlled system, not an accepted gait result.

Within the first 100 ms, the exposed-current case changes control by **0.0661464**.
The fully clamped negative control changes control by **6.95949e-11** in that same first 100 ms (maximum over the full post-deletion interval: **2.5945e-09**).

Fresh preallocated replay of the 25 s partial-neural CAD run: **all saved arrays bit-exact**.

## Recommended next changes

1. Use the frozen baseline, compiled edge inventory, and effective-gain loader
   as the starting point for each family. Keep the default and configured
   neuron censuses distinct. Reject nominal gain changes that do not reach
   actual matrix entries.
   Require bilateral physical motion, finite states, measured assistance,
   and contact/clearance checks before ranking locomotor scores. The frozen-left
   s3k anchor demonstrates why the scalar score alone is insufficient.
2. For the CAD robot, develop neural authority with the air/afferent-air/stance
   curriculum and a separate, explicitly bounded trunk-posture contribution.
   Check a delayed PF deletion before claiming that a PF layer controls muscles.
3. Audit the 27-route robot's pose-dependent torque vectors before fitting
   functional groups. Hip rotation and biarticular actions should follow actual
   routes, rather than inheriting the human-name map or the old splay-suppression
   edits. The foot-mount datum is still provisional.
4. Refit the synergy output for CAD geometry if syn6 is retained. Its present
   human IK-derived basis is fully connected but has not been validated for
   these routes, welded feet, or BPA force/pressure limits.
5. Preserve fixed-tau calibrated families as engineering models. Use each
   paper's actual channel dynamics and its own verified evaluator for a claim
   of equation-level replication. No fixed-tau fast wrapper should replace
   a bell-tau paper model.

No production winner, optimization database, dissertation chapter, or defense
deck was changed by this audit. New source files, edge inventories, trials,
verification records, and this report are separate reviewable artifacts.
