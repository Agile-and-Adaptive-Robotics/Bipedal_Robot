# Independent Chapter 6 and Appendix D review

Reviewed the complete revised `Final/chapters/60-sns.tex` and `95-AppendixD.tex`, the latest downloaded ZIP baseline, the rewritten Introduction, the relevant November 4–5 experiment records, and the new standalone figure previews. This review is independent of the chapter, appendix, and graphics authors. Manuscript sources were not edited by the reviewer.

## Required correction identified and resolved

The reviewed Chapter 6 workflow caption at line 37 initially described proposed orchestration and the MuJoCo–Simulink bridge, neither of which is explicitly drawn in the new figure. The visible diagram shows two mechanical preparation routes, the neural–actuator–mechanical feedback loop, sensory encoding, and calibration/testing/inspection. The chapter author applied a replacement describing those visible elements. The required caption correction is resolved. The knee caption also now distinguishes the zoomed pulse time axis from the longer response window.

## Structure and prose

- The Introduction now states the controller-transfer problem before introducing the motivating AnimatLab experiment. It connects that example to dissertation objective 4 and identifies the author's contributions separately from the software foundations (`60-sns.tex:6–26`).
- Methods establish model identity and interfaces before calibration and acceptance (`30–143`). Results follow the corresponding sequence: mechanical/interface checks, circuit/reflex behavior, suspended coupled tests, interventions, separate ground demonstrations, and execution time (`150–251`). This resolves the baseline's interleaving of tool inventory, fitting procedures, and unrelated controller results.
- Functional design, Hunt calibration, and the physical curriculum have distinct purposes (`108–133`). Hunt's method is described once. The curriculum table states physical tasks rather than implying that a single controller has completed all stages.
- The Discussion interprets what the experiments establish, explains the distinction between implementation agreement and physiological validation, and identifies the immediate integration requirements (`253–277`). It now answers the research objective rather than repeating the software inventory.
- Terminology and grammar are consistent with an engineering dissertation. Models, controllers, signal interfaces, and actuator formulations remain distinct. No substantial grammatical repair was identified in the revised chapter.
- Appendix D follows construction checks, interface definitions, implementation tests, staged experiments, supported evaluations, spiking comparisons, and reproducible sources. Its historical material is identified separately from the revised controller. The duplicate inverse-kinematics figure is included once with both required labels.

## Evidence and inference checks

| Claim or comparison | Review finding |
|---|---|
| BiW2Lvar+ air stages | Same learned weights and expression scale are reported for the two actual air tests. The changed sensory and command conditions are explicit, so the period change is not attributed solely to feedback (`60-sns.tex:182`). Standing and ground stages are not inherited from older controllers. |
| Air support | The 99.860% body-weight statistic is tied to the stage-2 suspension apparatus and the post-5 s window. It is not treated as a minimum ground-walking requirement (`191`). |
| Current noise and fixed perturbations | Ongoing RG-current noise, repeated noisy seeds, fixed parameter changes, and mechanical pushes remain separate tests (`141`, `193`). The six noisy air trials do not confer ground-balance robustness. |
| Layer stimulation | Silent-phase transmission failure is distinguished from the active-phase response. Later RG changes after PF/MN stimulation are interpreted within the coupled loop rather than as a direct synaptic gain (`195`; Appendix D stimulus section). |
| PF timing failure | The revised text incorporates the superseding November 5 audit: zero of twelve syn6 channels passed, variance components are not neuron counts, and s3k lacked declared numerical PF windows (`208–210`). The bounded repair remains unsuccessful. |
| W2L108 ground result | Assistance windows, retained rails/guides, treadmill-relative progression, and separate controller identity are explicit (`216–218`). Appendix D also records the large impact-force transient and the unsuccessful shuffle criterion. |
| CAD result | The current condition uses 54 ideal reference-tension actuators with SNS timing and engineering allocation, not 54 dynamically calibrated BPAs. The 77.5% support statistic is attached to 16–28 s; the limited 5/8 local perturbation result is retained (`220–229`). |
| Knee reflex | The reduced simplified force law is distinguished from the characterized Festo model. The response is a numerical circuit test with broadened ramps, not a physiological tendon tap (`162`). Appendix D retains the settled baseline, matched control, step sensitivity, and recovery thresholds. |
| Loaded cup | The corrected 2.00°/7.85° values are explicitly the 10 s endpoints of a ramp reaching 0.625 kg. No equilibrium sag or stiffness ratio is claimed (`171`). Independently checked against `figscripts/graphics_metrics.json`. |
| Force comparison | Matched pose/activation calculations are distinguished from dynamics; negative rectus-femoris and biceps-femoris-short-head R² values remain visible (`152`; Appendix D force section). |
| Speed | The 6.82× complete-process and 18.31× neural-call ratios share one declared benchmark. The later 1.613× warm-kernel comparison and 1.0156× full diagnostic comparison use the already accelerated sparse baseline and are not conflated (`233–251`). |

## Figures, captions, and source map

Independently inspected the new library, knee, cup, workflow, air-stage, and CAD previews. No clipped labels, obscured paths, or illegible legends were identified in those standalone images.

- The library figure shows twelve elements present in the current native library masks. It no longer contains a separate graded-neuron internal-flow panel, and its revised caption correctly omits that old panel.
- Knee and cup schematics each show all six actual connections: two homonymous Ia excitatory, two reciprocal Ia inhibitory, and two homonymous Ib inhibitory. The current builder files independently confirm those signs and targets. The old cup artwork's incorrect cross-Ib connection is not reused.
- Knee panel mapping is A circuit, B stimulus, C angle recovery, and D MN response. Cup mapping is A circuit, B arm geometry, C load ramp, D elevation, and E pressure command. Revised captions match these mappings.
- The air-stage display retains positive flexion magnitude while the model coordinate is flexion-negative. Its caption explicitly identifies that display transformation. The four selected rows are aligned over 7–11 s.
- The CAD image is a static sagittal/oblique view of the selected recorded state, not a trajectory or a force comparison. Its caption does not imply otherwise.
- Appendix D's script map now includes conversion/repair, CAD route and derivative checks, bridge verification and loop demonstrations, neural builders/tests, air/noise/stimulus drivers, PF timing checks, drawing checks, and the compact backend. Independently checked the added conversion, CAD, bridge, and controller-profile paths: they exist in the laptop checkout.

## Review limits

This is a completed source/evidence/caption review with standalone-image inspection. No unresolved material source-quality finding remains from this review. It is not a new simulation, a full physiological validation, or independent recovery of raw recordings that remain on the other machines. The coordinator reports that the assembled-package checks found no unknown references, duplicate labels/citation keys, or missing figures and preserved the unrelated original Introduction, published muscle chapter, torque chapter, and bibliography bytes. Those package checks are distinct from this independent prose review. Final compiled-page layout and figure placement remain pending because MiKTeX startup currently fails on registry permissions. Remote-machine prompts should request a clean full build, source-bound replays, and any newer records before replacing these explicitly scoped results.
