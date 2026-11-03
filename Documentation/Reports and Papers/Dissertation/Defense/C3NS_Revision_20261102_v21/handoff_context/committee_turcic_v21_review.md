# v21 mechanical and controls review

Independent review from a mechanical/controls committee perspective, completed 2 November 2026. This is not a statement from or impersonation of David Turcic. I did not author the reviewed slides, figures, or source data. No deck, manuscript, model, or optimizer was edited. Host verified as easteregg2.

## Coverage and verdict

Individually viewed the actual PowerPoint raster renders for **all main slides 29–59 and all question slides 78–99: 53 images**. The parent is replacing the already identified slide 47 crop and adding the slide 35 hybrid definition; those corrections require the parent's second-render check. Slide 55 already showed the shortened title when inspected.

No additional blocking numerical or mechanics contradiction was found. The mechanics sequence now establishes muscle force and geometry, shows the torque/length discrepancy, introduces corrections and fitting, demonstrates their transfer, and uses the remaining hardware/design limitation to motivate route redesign. The flexor and cam results distinguish model capacity from measurements. The transition into neural work shows a body/controller combination and its failed released condition before the new tools.

Three small wording corrections were sent to the parent:

1. **Main 37:** replace “Bracket stiffness is expressed in the local attachment frame” with “Effective stiffness is expressed in the local attachment frame.” The fitted parameters lump bracket, fixture, and test-stand effects; the following slide correctly calls them system stiffness. This makes the headings consistent with the authoritative interpretation.
2. **Main 50 footer:** replace “measured pressure used for hybrid prediction” with “model prediction at measured pressure.” The pinned-test hybrid on 35 uses measured BPA length, pressure, and moment arm. The redesigned flexor's yellow prediction substitutes measured pressure, so reusing an undefined broad “hybrid” label can confuse the comparison.
3. **Question 96:** “Complete pathway design / Still under review” is vague after the newly mapped selected circuit. “Physiological equivalence / Not established” would state the remaining limit directly, without implying the selected implemented pathway is missing.

Optional polish: main 42 and 44 retain tiny embedded “Unoptimized/Optimized” legends in addition to the much clearer native legends. Removing the duplicate legend through a native figure export would eliminate inconsistent terminology. This is less consequential than the slide 47 ambiguity immediately before route redesign. Question 86 and 91 duplicate the editor screenshot, but this redundancy is confined to hidden questions. The internal labels in the editor/Simulink screenshots are small; they remain usable as supporting overviews, not as an instruction to read every label from the back of the auditorium.

## Scientific checks

- The OpenSim illustration identifies a force vector, position vector, knee axis, and cross-product torque. The accompanying provenance computes the instantaneous axis for the CustomJoint and verifies the signed geometric moment arm against OpenSim; it does not assume a fixed hinge.
- The human and robot BFLH plots share torque and color scales. Their historical-model scope is cited. I did not independently rerun the source equations.
- The paired pinned-example plots use the same test and distinguish measured BPA length from rigid-geometry length. The correction is motivated by a mechanical discrepancy rather than merely by an aggregate error number.
- The root equation on 40 is dimensionally consistent: dimensional reference force times normalized force balances equivalent stiffness times displacement. The accompanying text states the normalized force is dimensionless.
- The split diagram preserves two flexor fitting/two validation tests and three extensor fitting/three validation tests, with the same stiffness pair retained. The source documentation correctly limits the interpretation of validation used in candidate selection; no blind-validation claim is visible in the inspected slides.
- The pinned-extensor panel's “44 cm” label is a rounding of the 43.6 cm measured length, not a new test. Checked against `Plot_KneeExtPin10mm_43cm.m` (`rest{1}=0.436`) and the dissertation's `30-torque.tex` training and figure caption. The file-label “43 cm” is distinct from this rounded figure label.
- The route criterion explicitly includes the required torque direction and a 5% magnitude allowance. It is presented as a design criterion, not as an experimental confidence interval.
- Slide 53 explicitly compares total BPA force at 620 kPa on the left and pressure for the same human torque target on the right. It does not incorrectly call the left curve the force required for the same torque.
- The physical cam comparison identifies vastus medialis, rather than total human extensor torque. Blue measurements remain visible above the yellow model points.
- The ideal reverse-pulley relation is correct for two parallel BPAs and two cable parts: tendon travel doubles, and tendon force is half the combined BPA force. The slide labels the relation as ideal.
- Whole-leg coverage is explicitly predicted, averaged across 27 sites, and distinct from a physical walking result.
- The released W2L failure and fixed-root condition are explicitly separated. Native video playback was not tested by this review; the parent owns that check.
- Questions 83–85 separate supported-cycle discrepancies, contact timing, and matched-input conversion checks. The own-peak normalization on 85 is visible, so shape agreement is not represented as absolute force agreement.

## Per-slide visual coverage

| Final slide | Viewed result |
|---:|---|
| 29 | Pass: mechanics divider connects BPA force to the knee-torque question; no irrelevant neural panel. |
| 30 | Pass: human paths and robot routes are visually paired and labeled. |
| 31 | Pass: actual model geometry, F, r_k, axis, and M_z relation are visible. |
| 32 | Pass: BFSH and VM targets are separately identified with signed torque axes. |
| 33 | Pass: common scales make the historical surface comparison meaningful. |
| 34 | Pass: apparatus labels identify load cell, swing arm, winch, and joint. |
| 35 | Pass pending known hybrid-definition insertion; one test is shown clearly. |
| 36 | Pass: matched measured/predicted active lengths support the mechanism. |
| 37 | Pass with effective-stiffness wording correction recommended. |
| 38 | Pass: offset, system stiffness, and short bending-loss segment are readable. |
| 39 | Pass: two orientation operations are shown using source geometry. |
| 40 | Pass: coupled force/displacement equation and physical meaning are visible. |
| 41 | Pass: fitting/validation split and shared stiffness are explicit. |
| 42 | Pass; optional removal of tiny duplicate embedded legend. |
| 43 | Pass: full GoF table with units and fitting/validation distinction. |
| 44 | Pass; optional removal of tiny duplicate embedded legend. |
| 45 | Pass: migrating rotation center and knee angles are shown. |
| 46 | Pass: physical extensor/flexor comparisons and photos identify the test. |
| 47 | Known crop and legacy-legend issue; parent has clean replacement asset. |
| 48 | Pass: design variables, required direction, target, and constraints are shown. |
| 49 | Pass as static poster; animation playback belongs to parent test. |
| 50 | Pass: measured blue points visible; footer terminology fix recommended. |
| 51 | Pass as static poster; hardware-overload claim is supported by source/narration, not the GIF alone. |
| 52 | Pass: geometry illustrates moment-arm/cam/clearance changes. |
| 53 | Pass: force and pressure scopes are distinguished on their axes. |
| 54 | Pass: VM reference, model capacity, and physical points distinguished. |
| 55 | Pass: ideal travel/force tradeoff is consistent with the mechanism. |
| 56 | Pass: mean predicted coverage values and averaging scope are visible. |
| 57 | Pass: neural section transition links mechanics and control. |
| 58 | Pass: separate controller and body origins are explicit. |
| 59 | Pass as static posters: fixed root versus released failure clearly labeled. |
| 78 | Pass: biological-engineering motivation without animal-study imagery. |
| 79 | Pass: two diameter datasets visibly support the length-force dependence. |
| 80 | Pass: conceptual feedback transition; not a detailed implementation schematic. |
| 81 | Pass as static poster: newer controller remains explicitly root-fixed. |
| 82 | Pass: force-law replacement is separated from mechanical coordinates. |
| 83 | Pass: supported human-model hip errors are visible. |
| 84 | Pass: contact timing and duty differences remain visible. |
| 85 | Pass: matched-input comparison explicitly states own-peak normalization. |
| 86 | Pass: editor overview; small UI labels should not be read individually. |
| 87 | Pass: current-balance implementation visible as a supporting schematic. |
| 88 | Pass: summed conductance/reversal-potential current equation is consistent. |
| 89 | Pass: implemented BPA block and empirical normalized-force law connected. |
| 90 | Pass: length/velocity versus force sensing distinguished. |
| 91 | Pass: duplicate editor overview, confined to questions. |
| 92 | Pass: six feedback paths and excitation/inhibition legend are visible. |
| 93 | Pass: stimulus, response, recovery, units, and sign convention visible. |
| 94 | Pass: mechanical state, controller, and actuator interface are separated. |
| 95 | Pass: separate CAD model identity and two views are explicit. |
| 96 | Pass with more precise remaining-limit wording recommended. |
| 97 | Pass: forward command and return measurement directions are explicit. |
| 98 | Pass: declared bridge-check outcomes are readable; scope remains short interface tests. |
| 99 | Pass: historical full-s3k and later standing configurations are explicitly distinct. |

## Limits

This review covers the above actual rendered images and selected source checks. It is not a new experimental validation, re-optimization, full manuscript audit, animation test, or rehearsed timing guarantee. The parent must re-render and inspect any corrected slides before delivery.
