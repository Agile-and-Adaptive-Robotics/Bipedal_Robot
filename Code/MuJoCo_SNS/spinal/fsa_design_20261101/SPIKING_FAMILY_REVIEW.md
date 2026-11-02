# Spiking design review

Machine: **easteregg2**. Design and electrical verification only. No body simulation, optimizer, Hunt training, curriculum, or training release occurred.

## New canonical central-spiking candidate

The accepted TimingZeroSFA v4 standing controller remains preserved. Its original left/right S3/S6 membrane-time assignments differ, even though its chemical weights are mirrored. The new candidate applies the explicit right-canonical policy before accelerator construction, followed by the graded functional profile and the frozen v4 spiking core. It does not average left/right weights, refit K6, copy dynamic phase states, or add the optional sensory-spiking port.

Actual verification covers 382 input bindings, 2,906 final nonzero chemical edge operators, and every spike-derived central readout. Maximum edge-operator error is `2.78e-17`. All 368 modified graded signed-pool contracts match their declared operating-point targets within `4.45e-16 mV`; extracted transient and off-point/fan-in behavior is recorded. These terminal maps are conductance-dependent and saturating, not globally linear.

With the same frozen central parameters, the 3 s recorded standing-input fixture passes: chemical-fraction RMSE 0.02582 and selected-MN RMSE 0.003531 mV. The 10 s deafferented fixture passes: chemical RMSE 0.16455; RG-E/F periods 2.130/2.149 s versus the graded reference's 2.132 s. Both sides have identical spike counts. Seven saved arrays, including names, repeat exactly in a fresh process. Disabled, shadow-poison, and complete-state controls pass. This is an electrical standing/slow-rhythm contract, not proof of alternating walking, a 0.559 s rhythm, or the new candidate's physical standing.

The caller-domain ledger derives all 382 current bounds from the actual runner formulas and the saved, fresh-verified bilateral static activation calibration. It does not infer bounds from observed trace maxima. Conditional muscle maxima are Ia 2.88, II 2.92, and Ib 1.125 nA at afferent scale 1, gamma 0, and vestibular proprioception 0. The 92 POST commands use `6 × static activation × 0.9619588149071685`, with that posture scale explicitly fixed. Standing DRIVE/POSTURE are 0/1; reflex-stand BAL commands and the selected zero-gain contact commands are zero. The separate deafferented fixture has DRIVE 2.5 nA and other input ports zero, plus its declared intrinsic startup kick.

All 382 endpoints pass unchanged, every out-of-domain upper probe is rejected, and changed caller parameters are rejected. A validation-only wrapper is numerically identical for 64 actual fixture steps and restores its complete state exactly. The graded profile receives only the 276 muscle-afferent domains; the full 382-port ledger stays separate. Adding this metadata and validation changes no neural numerical field. The generic graded ledger's body-domain flag remains unpromoted: the separate caller-mathematics proof is not a new mechanical-domain experiment.

Evidence, with source/profile hashes and eight scope-qualified checks awaiting root review:

- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design_20261101\spiking_family_review\CANONICAL_CENTRAL_FAMILY_EVIDENCE.json`
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design_20261101\spiking_family_review\canonical_central\BOUND_CANONICAL_CENTRAL_PROFILE.json`
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design_20261101\spiking_family_review\canonical_central\CALLER_DOMAIN_382.json`

Broader posture scales, drive schedules, cadence, gamma, vestibular, or contact gains require a new declared domain and appropriate functional qualification. All 92 MN pools remain, with 86 active biological routes and six force-pruned routes; six sagittal reference targets do not establish independent control of all 17 internal coordinates. The old v4 physical result belongs to its original configuration only.

## Optional sensory-spiking components

The typed N128 profile replaces six bilateral soleus Ia/II/Ib input pools. It has 768 actual spiking cells; the other 270 afferent pools remain graded. It is neither an all-sensory conversion nor a required addition to v4.

All selected passive outgoing receivers pass the registered 5% modulation-relative RMS and 10-degree extra-phase limits: 14 s3k, 18 w2lvar, and 12 syn6 edges. The two additional ionic bindings are **Ib_soleus_r → RG_E_r** and **Ib_soleus_l → RG_E_l**. Each produces 4,095 actual events in its isolated check; peak input current is 0.01712 nA. Actual membrane-increment binding error is `4.23e-16 mV`, and same-step ionic gate updates are exact against matched outgoing-off states. This establishes current binding; it does not establish longitudinal NaP phase equivalence under arbitrary modulation.

`D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design_20261101\spiking_family_review\SENSORY_FAMILY_EVIDENCE.json` records all checks, profiles, source hashes, component scope, and remaining base-family/composition requirements. Earlier failed designs remain preserved. The implementation rejects silently stacking onto an opaque central wrapper.

## Legacy selector

The 422-cell `AARL_NET=spiking` selector remains a reachable historical implementation: 310 spiking and 112 graded cells. Its calibration, SFA, and additive capped synapses are not the accepted v4 mechanism or a literal completed Table 2 design. The dissertation explicitly presents its untuned graded/spiking comparison and timing mismatch. Retaining that comparison does not require adopting it for new training or pretending the new v4 candidate repairs it.

`D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design_20261101\spiking_family_review\LEGACY_422_USE_CLASSIFICATION.json` links the source audit and manuscript locations. All five audited legacy source hashes still match. Any future active use requires a separate redesign decision; the audit does not grant training eligibility.

The root-owned global registry is unchanged. Both new family dossiers retain `reviewed_complete=false`, `training_allowed=false`, and `curriculum_allowed=false`.
