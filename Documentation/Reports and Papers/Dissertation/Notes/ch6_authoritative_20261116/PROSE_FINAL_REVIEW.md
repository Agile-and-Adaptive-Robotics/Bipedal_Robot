# Independent final prose and evidence review

Both edited canonical sources were read again in full on easteregg2 on 16 November 2026, after the primary scientific revisions and before the coordinator's final figure substitutions. This is a manuscript-source and evidence review, not a claim of successful compilation or rendered-page verification. This reviewer did not edit either canonical source.

## Verdict

The chapter now has a coherent scientific argument: controller transfer exposes the need for inspectable mechanical/neural interfaces; the methods define construction, calibration, and test conditions; the results demonstrate bounded capabilities and meaningful failures; the discussion distinguishes implementation agreement, functional response, engineering acceptance, and physiological validation. No broad structural rewrite is needed.

The contribution statement credits the existing software and published neural formulations while identifying the author's preparation, repairs, editable blocks, actuator integration, and experiment procedures (`60-sns.tex:26, 43-77, 110-112, 257-263`). Separate controller/body/actuator/support configurations remain identifiable. The successful air, W2L108, CAD, and offline PF cases are not presented as one completed curriculum (`60-sns.tex:45-61, 132, 210-229, 267`). Failure and limitation statements remain scientifically useful rather than being replaced by administration-only status language.

## Findings communicated to the coordinator

| ID | Finding at reviewed source lines | Needed disposition |
|---|---|---|
| P1 | `60-sns.tex:216` and `95-AppendixD.tex:277` bind W2L108 knee minima correctly to right 2.087/left 1.946 degrees, but omit their post-startup scoring window. AIR_EVIDENCE binds them to `1<=t<10 s`; these are not whole-trial minima. | Add that window to the minimum-flexion sentence in both locations. This clarifies interpretation without changing values. |
| P2 | Air/PF passages use a 0/5 mV transmission ramp and negative graded PF voltages without explicitly naming the rest-relative reference (`60-sns.tex:187,195`; `95-AppendixD.tex:223,241`). The reduced Simulink demonstrations use absolute voltages near -45 mV. | Explicitly label the SNS-Toolbox air/PF voltages as rest-relative, at least at the first relevant definition/caption. This prevents interpreting zero as an absolute physiological spike threshold. No numerical change. |
| P3 | The compact script-map row still says “whole-process timing” (`95-AppendixD.tex:465`), while the revised performance text correctly identifies evaluator diagnostic timing. | Replace with “scoped evaluator timing” or “timing scopes.” Do not reintroduce a full-process claim for the compact saved wall-time field. |
| P4 | The text says Hunt supplies the calibration method, rather than controller authorship (`95-AppendixD.tex:197`), but timing uses “Hunt controller” and the unexplained retained label “Hunt-PT+RC.” | Prefer the archived `hunt` implementation label and identify the later postural/Renshaw variant once. This is a naming clarification, not a new model or result. |

These findings were sent directly to the coordinator for resolution. Line numbers describe the fully read version and can move with final substitutions.

### Resolution readback

The coordinator's corrections were subsequently read back in the canonical files:

- P1 is resolved at `60-sns.tex:216` and `95-AppendixD.tex:271`: knee minima now explicitly use `1<=t<10 s`.
- P2 is resolved at `60-sns.tex:187,195` and `95-AppendixD.tex:235`: the air caption, pulse threshold, and PF transmission definition identify rest-relative voltages.
- P3 is resolved at `95-AppendixD.tex:459`: the compact script-map row now says numerical agreement and scoped evaluator timing.
- P4's attribution concern is resolved at `60-sns.tex:233,237,251` and `95-AppendixD.tex:425`: the older implementation is selected by the archived `hunt` name, and the later implementation is expressly labelled Hunt-PT+RC. Retaining that recorded configuration label is acceptable without inventing an unsupported acronym expansion.

No material prose/evidence finding remains open from this review. Figure substitutions and the separate build/render limitations retain their own verification status.

## Numerical and interpretation checks

The timing, cup, and knee values were checked against this reviewer's independently recomputed raw-array audit, not only against the prior laptop report. The other new result passages were compared with the independent AIR, PF, and CAD evidence reports and their result ledgers; this reviewer did not rerun those reviewers' scientific calculations.

| Passage | Reviewed outcome |
|---|---|
| Main cup (`60-sns.tex:171-176`) | Correct 2.00/7.85 degree comparison at 10 s/0.625 kg; increasing-ramp endpoint remains distinct from equilibrium. |
| Spiking cup (`95-AppendixD.tex:401`) | Correct matched endpoints 1.78/2.00/7.85 degrees at 10 s/0.625 kg. Original 2.02 remains explicitly a different maximum at 9.688 s. Raw maximum window is `2<t<=10 s`; endpoint and maximum are not conflated. |
| Knee (`60-sns.tex:161-167`; `95-AppendixD.tex:119-125`) | Correct 2 nA/30 ms pulse at 30 s; 1.09-degree extension, 111 ms latency, 6.355 s recovery. Six ramp changes and simplified force law remain explicit. Fine-step discrepancy now correctly rounds to 5.20e-5 degrees. |
| First performance pair (`60-sns.tex:233-251`; `95-AppendixD.tex:431`) | Correct 1,025 cells; native SNS versus sparse, one 16 s stage-1 protocol with 0.6 s held-pose preparation/15.4 s physics. Full process 199.624/29.284 s is kept distinct from neural calls 179.079/9.780 s. 45 arrays/model identity are correctly scoped. |
| Compact performance (`60-sns.tex:237-251`; `95-AppendixD.tex:433`) | Correct 1,033-cell later condition; median of three 3,000-update blocks, 1,200 comparison updates; warm 1.613 ratio is separate from single evaluator pair 92.944/91.516 s and 62-array/model identity. Exclusions around the evaluator timer are stated correctly. |
| Historical W2L performance (`95-AppendixD.tex:326`) | 1.425 ratio is bound to a bilateral-RG+Ib-port core, one 10,000-update block per backend after 5,000 parity steps. Separate 30 s body agreement is not used as its timing denominator. |
| Air stages (`60-sns.tex:182-195`; `95-AppendixD.tex:197-223`) | F/E commands 0 versus 0.15, sensory differences, shared onset seed, and disabled legacy drive are separated from unchanged learned weights. Principal-cycle window is in the table caption. 99.860% is retained with original total-model denominator, and 99.9877% is explicitly diagnostic, not a substituted adopted result. |
| PF follow-up (`60-sns.tex:208-210`; `95-AppendixD.tex:241-273`) | Original 0/12 remains; later 48-channel-cycle checks have an explicitly different denominator/criteria. Offline cadence success, live transmitted timing, and rejected body stages remain separate. Obsolete blanket “held” status is removed. |
| W2L108 (`60-sns.tex:216-218`; `95-AppendixD.tex:277-279`) | Side labels, anatomical versus native-coordinate distinction, 29.519% upward versus 35.069% signed net windows, rails/guides, and rejected shuffle refinement are correct, subject to P1's missing knee-window clause. |
| CAD (`60-sns.tex:220-229`; `95-AppendixD.tex:337-353`) | Ideal reference tension/engineering allocation remains distinct from dynamically calibrated BPAs. Final-rung 77.5% assistance and 7R/6L intervals are separate from whole-trial 2 ms knee extrema, 10 ms plots, and 5/8 local cases. Missing historical compiled-binary identity is disclosed. |
| Historical drive (`95-AppendixD.tex:319`) | Zero supplementary hand drive remains distinct from retained shared input 2.0773 nA. No harness-free walking is claimed. |

No unexpected change to an independently audited timing/cup/knee number, sign conversion, side assignment, or averaging window was found beyond the findings above. Numerical additions in the new air, PF, and CAD passages agree with their audit ledgers.

## Prose, captions, and conservation

Methods and results correspond: the model-preparation checks lead to force/interface results; neural blocks lead to matched calculations and reflex responses; declared intermediate timing leads to preserved PF failure and subsequent partial repairs; acceptance definitions lead to separately identified air and supported-ground cases. The chapter conclusion stays within these demonstrated capabilities.

The appendix supplies the configuration and scoring detail required for the main chapter without asserting physiological calibration from qualitative resemblance. The detailed tables are appropriate to the appendix. Repeated support and actuator distinctions are retained where they prevent conflation between cases; they do not require broad deletion for concision.

The current main air caption uses native flexion-negative knees, and the CAD trace caption distinguishes its negative knee guard from positive anatomical flexion magnitudes in prose. Original AnimatLab positive knee conventions remain explicitly identified. These conventions should not be globally converted under a single sign rule.

Final compact knee/cup and stimulus figure substitutions were still being coordinated during the read. Their actual panel letters, legends, source samples, final embedded fonts, and caption fit require the separate final figure review. This source review does not certify the later rendered pages or resolve the documented LaTeX runtime blockage.
