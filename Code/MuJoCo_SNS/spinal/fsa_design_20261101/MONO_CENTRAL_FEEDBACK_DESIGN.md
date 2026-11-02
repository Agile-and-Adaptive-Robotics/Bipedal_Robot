# Optional sagittal monoarticular central feedback

Machine: **easteregg2**. This is a **new unmasked canonical-template proposal**, not reconstruction or repair of the dissertation's historical trials. In particular, the historical s3k stage-5 recipe has central afferent pathways; its wiring must not be inferred from a later standing candidate with zero PF gain. Saved standing comparisons remain separate. No existing model, winner, dissertation result, source K6 basis, or manuscript was edited.

`mono_central_feedback.py` implements a strict opt-in compiled overlay for s3k and syn6. Passing `None` changes nothing. The full template is built before saved standing gain suppression; its original 128-step default anchor is verified. Syn6 uses the explicit right-canonical policy before adding the new routes, without averaging unrelated left/right components.

## Explicit routing scope

| Reviewed sagittal role | Muscle bases, on each side | New positive targets |
| --- | --- | --- |
| Flexor | iliacus, psoas, bifemsh, tib_ant | Ia and II → ipsilateral RG-F, InF, and supported F-family PF cells |
| Extensor | glut_max1/2/3, vas_med/int/lat, soleus | Ib and II → ipsilateral RG-E, InE, and supported E-family PF cells |

PF targets require an actual positive compiled PF→that-muscle-MN connection in the correct phase family. Within each muscle, these actual conductances determine relative PF feedback weights. No absent-path or all-PF fallback is invented. The matrix records eligible and rejected PF support, including cross-phase support. Syn6 uses its actual family/mixed metadata; mixed channels are excluded. After canonical projection, raw left-side NMF labels are not used as evidence of functional homology.

Twenty-two muscles receive these new routes. Sixty-four other active routes remain unreviewed for this proposal, and six force-pruned routes remain reserved. Biarticular, sartorius/gracilis/TFL, frontal/rotational, trunk, and distal multi-joint decisions require their own role-specific design. Existing arcs on those routes remain untouched. Tibialis anterior is represented only by its reviewed sagittal dorsiflexor role, not by an independently controlled inversion pathway.

Source inheritance is explicit in each `RULE_MATRIX.json`, with hashes. The local `replication/shinohara_rules.json` eq10 target set supports flexor velocity/length input to RG-F/IN-F/PF-F; eq11 supports extensor force input to RG-E/IN-E/PF-E. Our actual bounded SNS Ia/II encodings differ from those paper equations, and no paper gain is copied. Extensor II is the local same-group rule and requested engineering extension; it is not attributed to the paper's force equation. No universal extensor-Ia sign, new contralateral connection, state-dependent reflex reversal, or paper-exact connectome is claimed.

## Current function

Each postsynaptic cell and afferent type gets an independent **+0.5 nA** maximum-current contract at reference membrane voltage **2.5 mV**. Excitatory reversal is 8 mV, and the new source activation window is 0–5 mV. Actual canonical input maxima are Ia 2.88, II 2.92, and Ib 1.125 nA, conditional on gamma0, vestibular-proprioception0, and muscle-afferent scale1. Actual input weights and encoder leak convert these to source voltage bounds.

For relative weights `w_i` and maximum activation fractions `a_i`, the new conductance is

`G_i = 0.5 nA × w_i / Σ_j[w_j a_j (8 − 2.5) mV]`.

The complete fan-in participates in the denominator. One right half determines its verified left mirror. Different types are not pooled: Ia0.5 plus II0.5 can supply 1 nA to an F target. Current away from the reference voltage follows the actual conductance driving force; this is not a 0.5 mV equilibrium target, an activation floor, or a trained physiological gain. Direct central Ib excitation and retained local autogenic inhibition coexist; this addition alone does not implement a phase-dependent sign reversal.

Existing nonzero selected connections cause rejection rather than addition or silent replacement. Duplicate installation, unused profile fields, corrupted sparse caches, and out-of-domain source currents are rejected. The overlay preserves neurons, initial/dynamic states, other synaptic entries, local reflexes, ionic parameters, and clocks. Its explicit input validator must be used by a future caller; no production body runner is silently modified.

## Electrical verification

| New full-template proposal | New edges | Independent typed target pools | Final active edges | Largest native current error |
| --- | ---: | ---: | ---: | ---: |
| s3k, 410 neurons |168|32|1,354|1.44e−14 nA|
| syn6, 794 neurons |148|36|3,188|1.43e−14 nA|

At 0/25/50/100% of the declared source envelope, isolated actual sparse on/off membrane increments reproduce the analytical pool currents. Default-None preservation, the original default anchor, source state preservation, nonselected entries, and the negative guards pass. Independent fresh processes reproduce source hashes, graph-operator fingerprints, counts, and current errors exactly. These isolated tests clamp presynaptic voltages and set the reference postsynaptic voltage; there is no physical-model clamp or body run.

All authoritative review paths are in:

`D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design_20261101\mono_central_feedback\full_template\REVIEW_CHECKPOINT.json`

Each family's preferred `review_final` folder contains `PROFILE.json`, `RULE_MATRIX.json`, `RESULTS.json`, and **`FINAL_COMPILED_GRAPH.json`**, plus a fresh-process proof. The graph exports actual final neuron parameters, input-current mappings, every nonzero g/E/window, motor readouts, and unresolved coverage. Its edge-operator fingerprint is exactly the accelerator's. Use this graph for drawings; the original builder graph and older exported JSON do not contain these added connections.

## Application boundary

Construct the intended unmasked template and its declared symmetry policy first. Apply other static design profiles in a reviewed composition, then install this overlay and verify its final caches. Do not apply an existing-edge central-budget normalizer again to these new, already budgeted connections. Later static edits invalidate the recorded fingerprint and require an explicitly verified composed design.

```python
import mono_central_feedback as M
M.install(net, M.default_profile("syn6"), P=params)
M.verify(net, P=params)
u = M.validate_input_currents(net, u, P=params)
voltage = net.step(u)
```

The completed work is the new monoarticular wiring/current contract and its electrical proofs. Whole-neural phase behavior, muscle-waveform quality, biarticular/3D coverage, and body performance are not established. No Hunt training, optimization, curriculum, tendon/damping change, or global training release occurred.
