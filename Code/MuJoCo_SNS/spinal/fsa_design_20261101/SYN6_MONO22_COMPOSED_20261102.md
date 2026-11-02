# Syn6 composed central-feedback candidate

Machine: **easteregg2**. Named candidate: `Syn6_singleInE_canonical_FSA_mono22_extII`. This is a separate electrical design, not a replacement for the historical Syn6 controller, the earlier no-overlay composition, or Hunt1033. No body, fitting, curriculum, physics, tendon, damping, or manuscript changes are included.

## Runnable definition and final graph

- [PROFILE.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/syn6_single_in_mono22/composed/PROFILE.json>)
- [FINAL_COMPILED_GRAPH.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/syn6_single_in_mono22/composed/FINAL_COMPILED_GRAPH.json>)
- [COMPOSED_VERIFICATION.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/syn6_single_in_mono22/composed/COMPOSED_VERIFICATION.json>)
- [FRESH_PROOF.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/syn6_single_in_mono22/FRESH_PROOF.json>)

Builder: [syn6_mono_composed.py](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/syn6_mono_composed.py).

```python
import syn6_mono_composed as S
n, evidence = S.build_named_profile(profile_path)
u = n.make_inputs()
u[n.input_index('DRIVE')] = 2.5
voltage = n.step(u)
```

The existing `make_inputs`, `input_index`, and graded `step` interfaces remain. The new step wrapper validates input domains and compiled consumption. It adds validation overhead; its wall time must not be reported as sparse-kernel-only performance. A declared, fully prevalidated electrical stimulus schedule may instead use `n.compiled.forward(u)` with full `S.verify(n)` checks around each case; state/parameter ablations must be explicit and restored, not hidden by bypassing validation.

For a fresh reconstruction, with no body integration:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' `
  'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\syn6_mono_composed.py' `
  --profile 'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\syn6_single_in_mono22\composed\PROFILE.json' `
  --out 'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\syn6_single_in_mono22\NEW_EXPLICIT_REBUILD_FOLDER'
```

Completed outputs are protected from overwrite. Source/profile hashes are validated before building. The old no-overlay profile and all its referenced modules remain unchanged.

## Composition and exact scope

The builder reconstructs the reviewed full-recipe, canonical-right, single-heel-IN Syn6 composition, including its graded encoder/local-pool design, existing central-current budgets, and IBEXC coincidence. It then installs the **unchanged reviewed** 22-muscle mono-role overlay. No existing-edge normalizer is reapplied to the new, already budgeted routes.

All original numerical operators remain exactly unchanged. The new pairs are disjoint from both previously designed local terminal pools and existing central-budget selections. Because the graded verifier fingerprints the entire matrix, the new candidate explicitly records an approved disjoint fingerprint extension only after proving that every old entry and every membrane/ionic/input parameter is preserved. It does not weaken a verifier or rebind an old profile silently. The source canonical-right symmetry, single-IN heel policy, local functions, and IBEXC checks all pass on the final graph.

There are **794 neurons, 382 input ports, and 3,200 active chemical pairs**: 3,052 retained base pairs plus 148 new pairs. No neurons are added.

| New source type | Distinct source nodes | RG targets, edges | IN targets, edges | PF targets, edges | Total new edges |
| --- | ---: | ---: | ---: | ---: | ---: |
| Ia | 8 | 8 | 8 | 12 | 28 |
| II | 22 | 22 | 22 | 30 | 74 |
| Ib | 14 | 14 | 14 | 18 | 46 |

On each side, iliacus, psoas, bifemsh, and tib_ant provide Ia/II input to RG-F, InF, and their actually supported F-family PF cells. Glut_max1/2/3, vas_med/int/lat, and soleus provide Ib/II input to RG-E, InE, and actually supported E-family PF cells. PF allocation uses positive **final compiled** PF→MN support and canonical-right phase metadata; mixed channels are excluded rather than assigned an invented phase. All 22 reviewed muscles have eligible PF support in this composition.

Extensor II routing is the previously reviewed **engineering extension**, not an attribution to the literature's extensor-force equation. Sixty-four other active routes and six force-pruned reserved routes remain outside this new overlay. Biarticular, sartorius/gracilis/TFL, trunk, secondary-axis, and independent 3D central control are not silently completed. The existing Hunt four-cell RF/hamstring controller is untouched; no four-cell mechanism is newly claimed for Syn6.

## Current and input contracts

Each postsynaptic cell/type pool has a declared +0.5 nA full-envelope current at 2.5 mV postsynaptic voltage. The conductance law, including the actual source activation fraction and +8 mV reversal, determines weights. Relative PF contributions are preserved and the complete typed fan-in is normalized. This is a current contract, not a 0.5 mV equilibrium target, a physiological gain estimate, or a firing/phase guarantee. Current can reverse above the +8 mV reversal.

All 276 muscle-afferent inputs retain explicit clipped-caller domains: Ia 0–2.88 nA, II 0–2.92 nA, and Ib 0–1.125 nA, conditional on gamma coactivation zero, vestibular-proprioceptive addition zero, and muscle-afferent scale one. The 44 distinct new central sources are a subset of those guarded ports. Bounds derive from the caller's formulas and clamps, not measured biological maxima. The other 106 ports retain their existing source interface; this overlay adds no new reader or physiological range claim for them. The domain test holds those ports at zero, and the separate persistence fixture varies DRIVE over 0/2.5/0 nA.

**The 46 new direct central Ib connections bypass the separate IBEXC phase/load gate.** They excite RG-E, InE, and eligible E-family PF cells directly. Existing local autogenic Ib inhibition and IBEXC-gated excitatory MN output also remain. This composition does not make all load feedback phase-gated and does not implement a universal feedback sign reversal. The caller's muscle-force signal is distinct from heel/contact force. Those distinctions are explicit in the profile, rule matrix, final graph, and proof.

## Completed verification

- Every old operator, local/ionic/membrane/input parameter, and initial state is preserved; source/target pairs are unique, with zero collision against existing active pairs.
- All 368 local signed terminal-pool contracts pass; largest equilibrium error is 4.44e−16 mV.
- All 36 new typed central pools pass actual compiled on/off current checks at 0/25/50/100% input envelope; largest current error is 1.43e−14 nA.
- Single-IN/canonical-right/coincidence checks and final dense/sparse caches pass together. The final graph exports actual executed numerical fields, not stale source-drawing labels.
- Actual source-network stepping passes 120 domain/transient checks, followed by 128 exact dense/sparse steps and explicit command-persistence checks. Duplicate installation, corrupted cache, unknown profile fields, and out-of-domain input reject.
- Initial neuronal/gate states are restored before returning the runnable network.
- A separate fresh process produces byte-identical final graph and verification files: graph SHA256 `697678b3281d7fb3b2a9322c2bc77b0671e7d05beb932fcb4364b8aa5ddaddea`; proof SHA256 `c562ecebd872e7eccd580df815eafa0e404612454cd242445b2632d58ff2a029`.

These are completed electrical and routing proofs. The separate full-network phase probe, motor-waveform behavior, body-domain measurements, and walking performance are not inferred from them. The candidate remains explicitly excluded from automatic body/training release.
