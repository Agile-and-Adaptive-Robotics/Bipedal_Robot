# Repaired Syn6: explicit composed candidate

Machine **easteregg2**. This implements Ben's selected heel policy in a separately named candidate. Historical builders, source weights, saved results, K6 basis, and `syn6_graph_ledger.json` remain unchanged. No body, training, tendon, or damping change is included.

## Exact runnable profile and drawing source

Profile:
`D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\syn6_single_in\composed\PROFILE.json`

Final executed graph:
`D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\syn6_single_in\composed\FINAL_COMPILED_GRAPH.json`

Composition proof:
`D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\syn6_single_in\composed\COMPOSED_VERIFICATION.json`

Use **this final graph** for new drawings. It exports every actual node's membrane/ionic fields, all 3,052 active chemical operators with source/target indices and names, Gmax/reversal/threshold/saturation, 382 current-input bindings, and MN/actuator identities. Removed branches are recorded separately. The source builder's original connections and old drawing files do not describe the repaired final graph.

## Construction and function

The existing composed-family builder now accepts this named profile through `build_named_profile`. It reconstructs the actual full dissertation Syn6 recipe with current source, applies the explicit canonical-right policy, removes the two ipsilateral heel→InF conductances, then applies graded terminal/encoder design, budgets on eligible existing central routes, and actual phase/load coincidence. Acceleration follows these graph changes. All final cache/current operators are verified together. After electrical checks, the builder restores the original neuronal/gate states and returns the actual composed network for subsequent use.

The remaining heel policy is **HEEL_s→InE_s excitation**, with the actual source 0.5 µS/+8 mV/0–5 mV window. InE inhibits RG-F through the preserved source relay. Heel→InF has zero conductance and emits exactly zero direct current, including above synaptic saturation. Heel→PF_IN_E remains unchanged. No contralateral reset, inhibitory heel→InF replacement, or new event encoder is added. Existing direct/crossed reset combinations are rejected rather than stacked. Neurons may still interact indirectly through retained recurrence; zeroing one branch does not claim that every possible heel influence on InF disappears.

Excitatory/inhibitory labels describe reversal relative to rest. Actual current remains `g*s*(Erev−Vpost)`; an +8 mV reversal does not force positive current above +8 mV. Engineering current budgets are evaluated at their declared 2.5 mV operating point.

This candidate has **794 cells / 382 inputs / 3,052 active pairs**. It is distinct from the isolated single-IN proof's 3,035 pairs because the composed candidate also uses the declared canonical projection and coincidence design. It does **not** silently include the separate 22-muscle monoarticular-central overlay. Its source Syn6 central Ia/II omissions are not relabeled as fixed by a zero-group inventory check.

## Rebuild without body integration

```powershell
& 'D:\Anaconda\envs\myo\python.exe' `
  'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\gait_design_composition.py' `
  --family syn6 `
  --profile 'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\syn6_single_in\composed\PROFILE.json' `
  --out 'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\syn6_single_in\NEW_EXPLICIT_REBUILD_FOLDER'
```

Use a new output folder; completed graphs are protected from overwrite. The same module's Python API is `net, evidence = build_named_profile(profile_path)`. `net.make_inputs()`, `net.input_index(name)`, and `net.step(inputs)` retain the actual source controller interface. This is a neural builder, not an automatic body/curriculum launcher. Source hashes are checked before construction; changed sources require explicit review/rebinding, not silent acceptance. The previous `run(family,out)` and default builders retain their old behavior.

## Verification scope

The composed check passes 368 signed terminal-pool contracts (largest equilibrium error 4.44e−16 mV), 318 actual PF→MN operator slots, declared muscle-input bounds and out-of-domain rejection, phase/load coincidence truth/release tests, 120 actual input-domain network steps, and 128-step exact dense/sparse parity. Final whole-graph and cache operators match. The repaired branch remains zero after all transformations, and all initial dynamic states are restored before returning the network.

The isolated repair's fresh `RESULTS.json` is byte-identical under `D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\syn6_single_in\{proof,fresh}`. The separate composed fresh rebuild is under `D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\syn6_single_in\composed_fresh`. Both `FINAL_COMPILED_GRAPH.json` and `COMPOSED_VERIFICATION.json` are byte-identical across the fresh processes; `D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\syn6_single_in\COMPOSED_FRESH_PROOF.json` records their hashes.

These tests establish the selected graph and its scoped electrical consumption. Full-network phase-reset behavior, body ranges, Hunt fitting, and curriculum performance remain separate. No old standing or walking result is inherited.
