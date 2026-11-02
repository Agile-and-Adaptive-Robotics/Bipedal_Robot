# Optional graded Gait design profile

Machine: easteregg2. Implementation: `graded_profile.py`; builders: `build_network.py`, `build_network_w2lvar.py`, and `build_network_syn6.py`. All hooks default to `graded_profile=None`, with a lazy import only when supplied. Existing numerical defaults, source NMF files, and historical configurations are preserved.

The three-family results are electrical design checks of canonical builders, not new standing or walking results. The authoritative `CAMPAIGN_GATE.json` remains closed for Hunt training, curriculum execution, and whole-body optimization. `assert_campaign_allowed` delegates those operations to the hashed all-family `design_gate.require_calculation_release` check. Local pathway results cannot release that gate.

## Applied units and timing

The implemented SNS equation uses seconds, mV, nA, microSiemens, and numerical capacitance in microfarads. Thus `C=0.03`, `G=1` means 30 ms and physical capacitance 30 nF. The physical conversion is `C_nF=1000*C_backend`; the profile changes actual compiled `c_m` and `time_factor_membrane` fields.

The engineering context is the shortest usable cached stride estimate, 0.559 s, and its ankle 99%-AC-energy fourth harmonic, 7.155635 Hz. It is not a universal human maximum or a measured neural bandwidth. A declared 10-degree **per-encoder** budget gives a continuous-RC bound of 3.921848 ms. At the unchanged 2 ms update, that capacitance produces 10.064927 degrees of clocked phase lag, so the applied value is instead 3.896229 ms.

For the explicit update, let `a=1-dt/tau`, `theta=2*pi*f*dt`, and use

`H=(1-a)*exp(-j*theta)/(1-a*exp(-j*theta))`.

The output timestamp is the end of the update; the input is sampled at its beginning and held during that update. Solving this discrete transfer function for the declared phase budget determines the actual capacitance. The corresponding discrete pole time constant is reported separately. Actual SNS-Toolbox probes measure 2.508619 degrees at the stride fundamental and 10 degrees at the fourth harmonic. An additional stale-input age of 0–2 ms adds 0–1.288014 degrees or 0–5.152057 degrees, respectively. This extra age is not hidden inside the per-encoder budget.

Only the 276 isolated Ia/II/Ib input encoders receive the new time constant. They must be passive, zero-rest, zero-bias neurons with no incoming synaptic conductance. NaP parameters, RG/lamination/PF/MN capacitances, PF shape and phase assignments, and actuator dynamics remain unchanged. There is no whole-loop phase-budget or intrinsic-period claim.

## Function-derived signed terminal pools

The [graded functional-subnetwork approach](https://www.frontiersin.org/journals/neurorobotics/articles/10.3389/fnbot.2017.00037/full) distinguishes desired signal transformations from raw conductance values. This profile derives four existing sensory-to-MN terminal pools: Ia homonymous excitation, II homonymous excitation, Ib autogenic inhibition, and aggregate Ia reciprocal inhibition. It creates no neurons, edges, or missing homologues.

The chosen engineering calibration is a +/-0.5 mV change, or +/-0.1 of the 5 mV activation span, about an operating voltage of 2.5 mV. These are declared functional test targets, not physiological findings or Hunt-fitted recruitment. A 1 nA pulse is a calibration amplitude, not a physical maximum. For relay paths, the actual isolated afferent-to-relay response determines the presynaptic voltage used in the terminal-pool calculation.

For each existing source, `a_i=clip((V_i-Elo_i)/(Ehi_i-Elo_i),0,1)`. With existing relative pool weights `gbar_i`, desired output `U*=U0+delta`, and operating current `I0=G*U0`, solve

`scale = G*delta / sum(gbar_i*a_i*(E_i-U*))`.

Each actual pool conductance becomes `scale*gbar_i`. This accounts for sign, reversal, source voltage range, postsynaptic operating point, and conductance-dependent shunting. It is not division by fan-in. Unreachable signs, absent contracts, asymmetric derived homologues, and unknown fields are rejected before committing the profile.

The isolated calibration endpoints pass to roundoff for all 368 terminal pools per family. Additional amplitude tests include 0, 0.25, 0.5, 1, 2, 5, and 20 times the calibration pulse at four postsynaptic operating points. The response is a rational/piecewise-saturating map, not a global linear gain. Simultaneous signed-pool checks show nonadditivity; even two positive 0.5 mV targets combine to approximately 0.9167 mV rather than exactly 1 mV. Operating-point errors and activation-range crossings are preserved in the records.

The admissible **engineering** encoder-current range is derived from actual outgoing voltage windows and port conductance (0–5 nA for these builders). Physical signal maxima remain distinct and unresolved. No current clamp is introduced. Actual isolated 20 nA probes show membrane voltage above 5 mV while outgoing synaptic activation saturates at 1, followed by a measured release delay. This is bounded synaptic representation, not evidence that physical inputs cannot exceed the declared engineering interval.

## Actual consumption and historical mirror projection

`verify_profile(net)` hashes all relevant final electrical fields, including ionic parameters, capacitances, synapses, windows, and input mappings. It also checks accelerated edge caches. A later gain or mirror pass invalidates this ledger; it must be followed by explicit re-derivation and verification. The consumed profile passes 128-step dense-versus-sparse neural trajectory identity. Unknown top-level or nested configuration keys are rejected.

The existing `gait_standing_family.tie_weights` changes no pairs in canonical s3k/w2lvar. For syn6 it changes 189 pairs and adds 67 nonzero compiled edges. That historical behavior is audited, not approved as a new biological symmetry policy. Its changes invalidate the old profile hash; reapplication verifies electrical consumption only.

In the unchanged K6 source, same-index right/left W cosine similarities for S5/S6 are about 0.587/0.193, and their per-leg phase-profile correlations are about 0.438/0.152. S5 peaks at 26% versus 58.5% of its respective leg cycle. S3 membrane tau is 72/48 ms, while S6 is 48/72 ms. Averaging same-index PF outputs does not establish functional homology. The component ledger includes original H peak times, per-leg phase, stance fraction, tau assignment, and added source targets. No permutation, NMF refit, or canonical-side mirrored candidate is applied. Strict-symmetry release remains blocked pending an explicit verified policy.

## Use after final network construction

```python
from fsa_design_20261101.graded_profile import (
    default_profile, apply_profile, verify_profile, adapt_inputs)

profile = default_profile("syn6")
net = build_network.build(actuator_names, graded_profile=profile)
# Apply any explicitly chosen final legacy gain/tie transformations here.
apply_profile(net, profile, family="syn6")
# Install the existing accelerated wrapper after applying, then verify again.
verify_profile(net)
currents_nA, measured_range = adapt_inputs(
    net, supplied_currents, unit="nA", enforce_domain="engineering")
```

Calling `apply_profile` on a baseline built with `graded_profile=None` after its final transformations is also supported. Prefer this order for externally assembled networks. `enforce_domain=True` requires independently supplied physical bounds; `"engineering"` checks the window-derived interval; `False` measures/converts without clipping. Each accepted field has an explicit consumer. The source builder graph retains baseline numerical labels, so diagrams of the applied variant must use its final compiled ledger rather than inherited weight labels.

## Dimensional coverage and remaining scope

Gait2392 has six floating-base and 17 internal anatomical coordinates. Six primary sagittal hip/knee/ankle reference coordinates do not establish control of the remaining hip rotations, subtalar/MTP, or lumbar coordinates. The observed native 105-coordinate conversion adds 82 follower coordinates, not 82 anatomical joints. Soft apparatus guides provide reactions without deleting anatomical degrees of freedom. Zero vertical support does not mean free balance.

`gait_coverage_ledger` distinguishes actual population groups, distributed muscle routes, passive contributions, and apparatus restraint. The separate read-only `syn6_torque_span.py` diagnostic found rank17 for the full92/represented86 muscle activation maps and rank12 for the bilateral K6 activation basis at four postures. It excludes extra reflex/posture inputs and is not the compiled SNS Jacobian, a cone-feasibility result, or a stability proof.

Unchanged PF/NaP/stance-load and other reflex branches still need their appropriate functional and interaction checks. Canonical electrical checks do not transfer to saved standing winners or release body calculations. Full physical domains, whole-loop performance, syn6 component homology, and dimensional control coverage remain separate review items.

External evidence on easteregg2: `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\graded_gait\`. Each family has an unchanged pre-hook baseline, a consumable profile JSON, and a final electrical-design report. No body simulation, Hunt training, curriculum, NMF modification, or manuscript edit occurs in this work.
