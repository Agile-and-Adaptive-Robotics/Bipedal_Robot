# Wiring review: quick disposition

Full decisions: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\DISSERTATION_WIRING_DECISIONS_20261102.md`.

- **Historical s3k:** central Ia/II/Ib were active in the full dissertation recipe; later standing exports cannot illustrate that graph. Pruned central Ib is not removal of local Ib.
- **Heel sign:** InE inhibits RG-F; InF inhibits RG-E. W2L's contralateral InF path therefore inhibits contralateral RG-E. Syn6's two positive ipsilateral IN inputs influence both half-centers inhibitory. Actual current still changes sign with `Erev−Vpost`.
- **Ben's selected Syn6 repair:** retain ipsilateral heel→InE excitation, remove heel→InF entirely. This is a new default-off compiled repair, not either earlier A/B policy or an edit to historical results. Two conductances change; no new edge or encoder is added. Direct/crossed reset stacking is rejected.
- **Overwrite:** inspect final SNS matrices when duplicate source/target definitions exist. Full s3k's final heel→InE uses `ib_e_central`, not the originally declared heel gain.
- **Contact timing:** generic tonic contact plus optional event additions differs from the selected Hunt caller's finite hindfoot landing pulse. The integrated Hunt profile retains direct cross-heel reset and adds no heel-B stack.
- **Syn6:** six NMF/half-center channels differ from Di Russo's five analytic primitives and phase oscillators. K8 inverse reconstruction is not a new forward/controller/body pass.
- **Composed:** existing-family graded/current/coincidence candidates and Hunt1033 local/timing/coincidence buffers have composition evidence. Historical recipe terminal checks and the new 22-muscle central overlay are different scopes.
- **Non-Gait:** W2L108/704 remain configuration-specific; original rat120, new122, and excluded diagnostic134 are distinct; CAD262 pressure-pair remains hybrid; Simulink and reduced mechanisms retain their actual implementation scopes.
- **Still required:** actual body-domain currents, few-cycle hierarchy/force measurements, PT calibration, and individual pilot acceptance. An isolated electrical pass does not complete these tests.

Historical models, parameters, and results are preserved. The explicitly selected repaired candidate changes two compiled conductances; no physics or manuscript change is included.
