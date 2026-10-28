# AH001 supported-control figure and results

Session 1, easteregg2, 28 October 2026. Repository: `D:\GitHub\Bipedal_Robot`.

## Defense slides

The animation is slide 29 in `Defense\AE_Revision_20261023\output\Bolen_Defense_Assertion_Evidence_20261028_AH001_v8.pptx`. The still and editable comparison table are hidden backup slide 39. The source v6 deck, its 37 original slides, and their notes are preserved. The GIF is embedded byte for byte. Native PowerPoint rendered both added slides for inspection.

## Dissertation placement

The current source has the CAD integration paragraph in `ProofFinal\chapters\60-sns.tex`, under Repeatability and Mechanical Support, and the detailed CAD-Derived Robot Integration subsection in `95-AppendixD.tex` (`sec:app_cad_supported_stance`). The SNS discussion ruling, `Notes\revision_campaign_20261015\SNS_DISCUSSION_REDESIGN_20261017.md`, calls for a tools demonstration with explicit limits.

Place the still next to that scoped CAD integration paragraph in Chapter SNS. Put the matched results table and assumptions in the existing Appendix D subsection. In `70-overall.tex`, use a brief cross-reference if needed; a duplicate figure or table adds little. Appendix B describes a different W2L connectome and is not the place for this CAD model's census or results.

The older 20 s stance and 94.68% support record should remain separately identified. The new 40 s trial is not a matched support-reduction comparison with that older result. No master dissertation source was edited in this pass.

## Proposed caption

Anterior (A) and oblique (B) views of the CAD-derived robot model at the same recorded simulation state. The right and left legs are blue and gold, respectively, and selected hip-abducting routes are green. The 75 kg model has 54 reference-tension actuators and a hybrid neural and postural controller. Vertical assistance supplies 80% of body weight in the displayed interval; lateral and angular pelvis assistance also remain active. Route attachments and the welded foot mount are provisional, and passive knee-stop surrogates are assumed. This figure demonstrates supported model integration, rather than validated BPA actuation or independent walking.

## Matched comparison

Recorded interval: 15–40 s. Both runs use the same geometry, 2 ms physics step, initial pose, external assistance, reference tensions, and other control settings. The control bundle adds separation feedback (160 N/m, capped at 20 N) and raises the floor on currently abducting routes from 5% to 6%.

| Quantity | Extra abductor control off | On |
|---|---:|---:|
| Reversed ankle ordering | 67% | 0% |
| Minimum lateral ankle separation | −38.3 mm | 142.4 mm |
| Median lateral ankle separation | −6.2 mm | 159.0 mm |
| Mean right TFL activation | 22.3% | 81.8% |
| Mean left TFL activation | 99.6% | 99.8% |
| Mean vertical assistance | 80% BW | 80% BW |
| Maximum sampled knee-loop position error | 4.485 mm | 0.387 mm |
| Maximum sampled opposite-leg mesh penetration | 2.700 mm | 0 mm |

Ankle ordering and sampled mesh contacts are separate checks. Toe-out remains approximately 25–28 degrees in the displayed case. Knee motion after the support ramp is small; regular stepping is not established. The exact-fast and reference backends reproduce all saved arrays bit for bit in fresh 40 s runs. This establishes numerical replay within the stated model, not robustness or hardware validity.

## Attachment interpretation

Every route site belongs to a named rigid body, including moving and conditional sites. A site need not lie on a CAD surface to transmit tendon force. Such an offset is mechanically equivalent to an ideal rigid bracket or guide at the prescribed body-local location. There is no automatic nearest-body connector, bracket mass, bracket compliance, collision-aware wrapping, or physical bracket design in this diagnostic. Consequently, visible floating routes are an engineering packaging assumption that must be audited before hardware claims.

The foot is welded using a top-of-assembly datum, not a verified mating interface. The initial mount-to-sole height is approximately 103.7 mm. Ben's suspected 60 mm discrepancy is a hypothesis, not a measured fabrication error. A separate proximal-foot-translation experiment is recorded in `Code\MuJoCo_SNS\robot_torso\foot_mount_sensitivity_spec.json`; it does not shorten or redesign the tibia.

## Evidence

Model and traces: `Code\MuJoCo_SNS\robot_torso\ah001_abduction_diagnostic_v2.xml`, `runs\abd_defense_exact.*`, `runs\abd_defense_control_off.*`, their abduction assessments and geometry audits, and `abduction_defense_verification.json`. Asset hashes and recording times are in `asset_provenance.json`. These paths are relative to the repository root above.
