# Independent final air-evidence review

Reviewed the edited canonical Chapter 6 and Appendix D air/noise/stimulus/W2L108 passages against AIR_EVIDENCE.md, AIR_RESULTS.json, AIR_CYCLE_RESULTS.json, raw arrays, saved configurations, and executable source. Review was performed after the primary edits, before the coordinator's final figure substitutions and full build. No canonical source was edited by this reviewer.

## Material findings communicated for correction

1. **W2L108 native coordinates:** the reviewed Appendix D line 289 added “corresponding to native maximum coordinates of −2.087 and −1.946.” This conversion is not the W2L108 native coordinate convention and must be removed. `virtual_walker_20261031/run_virtual_walker.py:140` stores `angles=initial+degrees(qpos[qadr])`, with initial anatomical knees left10°/right55°. The native positions are pose-relative hinge offsets. They are not the converted Gait2392 model's absolute flexion-negative anatomical coordinates. Retain minimum anatomical flexion **right2.087°/left1.946°**, with explicit side labels. Communicated immediately to coordinator.
2. **Air cycle-scoring window:** the reviewed air table reports run duration16s but does not state that counts, period, and phase come from **5≤t_state<14s**. Add this to its caption or adjoining paragraph. Maximum knee angles use the full retained physical trajectory; support uses legacy samples from5s onward. These distinct masks should not be inferred from one another. Communicated to coordinator.
3. **Stimulus figures:** the inherited figure strips fall below10pt at their scaled size. They show the **7.65s active-phase series**, not the7.1s silent controls. New vector panels were regenerated directly from raw aligned arrays and visually inspected, with10.2pt type at14cm width. They are ready for the coordinator's portrait-figure substitution; active files and canonical references remain the coordinator's responsibility.

## Results confirmed

- Stage1 F/E commands0 and sensory return disabled versus stage2 F/E0.15 on both sides and live sensory return. Common100ms initialization, expression0.1, and disabled postural/balance/function channels remain correctly distinguished. The revised text no longer attributes the period change solely to feedback.
- Recomputed principal cycles4/4 and7/7, periods2.466/1.300s, phase errors0.0336667/0.0932768 cycles, and period CV agree exactly with the archived assessment. The detector's rejection controls pass under the declared available base environment.
- Adopted stage2 support99.860% retains legacy total weight769.737726N and mean768.661982N. The revised denominator disclosure correctly distinguishes the0.1kg fixed ground body and aligned pelvis-subtree denominator768.756726N. Figure force units in newtons prevent silent normalization changes. The additional99.9877% is explicitly a diagnostic ratio using the same legacy numerator, not a replacement result.
- Stimulus targets are right RG_F, right hip-flexor PF, and right iliacus MN. Earlier onset7.1s and later onset7.65s,100ms duration,0/5mV outgoing ramp, and3.4s onset-relative maximum-difference windows are correct. Rounded1.14/0.16/1.79/6.04mV results retain their original values and context. The very small silent-PF MN response is correctly scoped.
- Six independent noisy cases, one zero control, and one repeated seed are distinguished correctly. Raw repeated noise, neural, and qpos arrays agree exactly.
- W2L108 explicit side labels are correct: right2.087°/left1.946°. Upward29.519% uses independently nonnegative rail/guide/harness forces over1≤t<10s; signed net35.069% uses0≤t<10s. The411.2254N model-weight denominator is correct. Rails/slack guides and moving treadmill remain explicit.
- Separate accepted supported stepping and failed short-shuffle refinement remain distinct; no continuous BiW2Lvar+ curriculum or unsupported-walking claim follows.

## New stimulus figures

Generator: `AIR_stimulus_figures.py`; provenance: `AIR_stimulus_provenance.json`.

- `AIR_stimulus_neural.pdf`:14×15cm, four rows for current/rightRG-F/right hip-FPF/right iliacusMN.
- `AIR_stimulus_mechanics.pdf`:14×12.5cm, three rows for right hip/native flexion-negative knee/ankle.
- Both use7.5–11.05s display, all1776 saved samples per condition, actual7.65–7.75s pulse shading, distinct line styles, baseline plus separateRG/PF/MN trials, and minimum10.2pt text. Native clocks are retained; no interpolation, smoothing, phase shift, amplitude rescaling, or decimation is applied. The only joint conversion is radians to degrees.

The PNG previews were inspected for label overlap, legend clipping, sign convention, and curve identity; no material defect found. Final manuscript scaling and rendered-page inspection remain the coordinator's build review.

## Remaining source limitation

The exact launch-time air `body_probe.py` hash in SOURCE.json was not recovered among later snapshots. The captured physical evaluator hash is matched, raw measurements are independently verified, and the later stable command implementation is consistent with the saved manifests. This source-binding limitation remains documented; do not claim complete launch-source recovery. The initial myo SciPy import error was resolved for cycle calculation using the available base environment and is not an unresolved cycle-verification gap.

## Closure after canonical readback and figure installation

All three material findings above are resolved in the current canonical sources:

- `95-AppendixD.tex:197` explicitly gives the principal-cycle, period, and phase window as **5≤t<14 s**, within each 16 s run. The distinct support statistic and clock/denominator explanation remain at line 211.
- `60-sns.tex:216` and `95-AppendixD.tex:271` explicitly bind W2L108 minimum anatomical knee flexion to **1≤t<10 s**, right **2.087°**, and left **1.946°**. The appendix now explains initial-pose offsets and distinguishes this convention from Gait2392 coordinates. The incorrect native-coordinate conversion is removed.
- `95-AppendixD.tex:219–231` includes the two new stimulus PDFs at **14 cm**, in ordinary portrait figures. Captions explicitly bind them to the **7.65–7.75 s active-phase** pulses, retain the common **7.50–11.05 s** display, identify the right-side targets, and separate the silent-phase textual control. Knee coordinates remain native flexion-negative.
- The voltage reference is explicit: `60-sns.tex:187` labels the new air-recording voltages rest-relative, line 195 names the rest-relative 0 mV stimulus transmission threshold, and `95-AppendixD.tex:235` defines the PF timing equation's V as rest-relative membrane voltage.

Active PDF bytes match the inspected Notes originals exactly:

| Asset | SHA256 of both active and Notes copies |
|---|---|
| `AIR_stimulus_neural.pdf` | `859087a4d0469922700c469c518d14e1a1700a10a06e69ccbc9fbcd18de306a1` |
| `AIR_stimulus_mechanics.pdf` | `acc06620bb9685a09a68a31dba33621f1b52a82a169e344cfd3b6920b6f712f7` |

Read-only PDF text extraction confirms a **10.1999998 pt** minimum in both active PDFs. Page dimensions are **396.8504×425.1968 pt** (neural) and **396.8504×354.3307 pt** (mechanics), consistent with their intended 14×15 cm and 14×12.5 cm sizes. No generator was rerun during closure. The matching active hashes bind this check to the previously inspected previews.

This closure verifies source/evidence consistency and figure installation. The coordinator's full build and rendered-manuscript-page inspection are separate checks and are not claimed by this review.

### Additional command-provenance qualification

The final provenance follow-up checked whether the exact-hash physical evaluator and saved air-run configurations alone prove the consumed F/E command contrast. They do not: the body configurations omit the numeric command values, the exact evaluator delegates neural construction to body_probe.py, and the recorded launch-time adapter hash remains unrecovered. The fitted C manifest records F/E=0.15; the retained later adapter specifies stage1=0 and stage2=those fitted values. All three Git-history versions were checked read-only and also fail to match the original launch hash (AIR_SOURCE_CHECK.json). No numeric F/E command samples were separately logged in the aligned capture.

Accordingly, **0 versus 0.15 is the declared stage contrast reconstructed from the retained adapter and calibration configuration, not independently recovered consumed-command telemetry**. The coordinator was advised to narrow certainty accordingly in the manuscript and preserve this limitation in the evidence report. This qualification does not change the independently recomputed joint motion, support forces, pulse samples, or noise-repeat results. It reinforces the decision not to interpret the two air-stage periods as an isolated causal test of feedback.

### Final command-provenance prose closure

Read back the final air-stage paragraphs in both canonical files. Chapter 6 now attributes the command/sensory definitions to the **retained adapter**, explicitly states that complete applied-command histories and the original launch adapter were not recovered, and rejects a controlled feedback-only interpretation. Appendix D independently distinguishes the fitted manifest's 0.15 command value, the retained adapter's zero/0.15 stage definitions, and the unrecovered consumed samples. The common initialization seed and remaining command flags are likewise explicitly scoped to that adapter, rather than presented as recovered launch telemetry.

The **prose finding is resolved**. Recovery of the exact launch adapter and full consumed-command histories remains an **evidence gap**, plainly disclosed in both manuscript passages. No additional source search, plot generation, or numerical reanalysis was performed for this final readback.
