# Appendix D revision and evidence record

This pass edits `Final/chapters/95-AppendixD.tex` from the most recent downloaded source ZIP, not from the older Final appendix. Baseline SHA256: `85b1feb584807723c293d0f1b08cb29f3e061457fdc655e37f97342076b2f2a3`. The original ZIP, result files, model sources, and archived graphics are unchanged.

## Structure and editorial decisions

The appendix now follows mechanical/interface checks, neural implementation checks, staged calibration and actual-body tests, controlled interventions, spiking comparisons, and reproduction sources. Chapter 6 owns the main argument; Appendix D supplies protocols and limitations.

- Retained the specific conversion repairs, twelve-muscle force table, settled knee-transient protocol, matched W2L rhythm-core comparison, reduced Rybak mechanism, and decisive reduced-Deng layer controls.
- Retained the historical Gait2392 absolute-angle comparison once. Both `fig:gait2392_ik_overlay` and `fig:app_gait2392_ik_overlay` now identify the same float. The historical controller is distinguished from BiW2Lvar+.
- Retained the independently replayed W2L condition and its matched pathway removals. Its 96-neuron circuit, later W2L108 body, and revised Gait2392 controller are not treated as the same configuration.
- Added the latest portable BiW2Lvar+ air-stage, current-noise, and phase-specific-stimulus evidence. Explicitly preserved the distinction between offline activation fits and actual-body acceptance.
- Added the failed syn6/s3k PF-timing audit as a verification result. Passing output fitting no longer implies verified PF timing.
- Updated CAD integration to the selected 77.5% measured-support condition, with ideal reference tensions and engineering force allocation stated explicitly.
- Added a purpose-based source-script map covering conversion, route verification, graphical block construction, bridge execution tests, training, interventions, and timing audits. Long listings are avoided.
- Reduced the hardware material to its interface boundary and unfinished status; no hardware throughput is claimed.

## Removed or consolidated active material

User authorization covers removing irrelevant sections, scripts, and figures from the appendix. No underlying scientific result or script was deleted.

1. Removed the three low-resolution legacy subsystem screenshots (`RGnetwork.PNG`, `PFnetwork.PNG`, and `sensoryMotor_Network.PNG`). Retained the requested original full shared-RG source architecture as one appendix figure.
2. Removed the composite converted-robot/reduced-knee overview, which identified two different bodies but added no verification beyond the retained construction/interface explanation and Chapter 6 model descriptions.
3. Removed the older bilateral hip-PF inhibition figure in favor of the revised Gait2392 matched-current experiment. The reduced-Deng neural-only controls remain explicitly separate.
4. Consolidated the two identical IK figure inclusions into one. No angle-error table values changed.
5. Removed the four-family reference-score catalogue and multirow pathway-score catalogue. Their source results remain preserved, while the retained text explains transfer limitations and the danger of the `-25` no-cycle fallback. The matched W2L intervention remains the detailed pathway example.
6. Replaced the lengthy old CAD stance/foot-mount chronology and seven-row abductor table with the current supported-clearance experiment. The old matched 80%-support comparison's 67%-to-zero reversed-ankle result remains identified as a separate causal test.
7. Reduced the stronger-current Deng table from nine rows to four decisive rows: RG and PF two-cycle controls, MN two-cycle inhibition, and the one-second PF case separating motor shift from RG reset. Five shorter/weak control rows are preserved in the source ZIP and source results. No retained value changed.

Removed labels are enumerated in `APPENDIX_CHANGE_MANIFEST.json`. Root coordination confirmed the rewritten Chapter 6 will not retain those references. All other chapters are outside this agent's edit scope.

## New numerical results and source bindings

All paths in this table are under `Defense/Session1_Handoff_20261104/` unless stated otherwise.

| Added or updated result | Authoritative portable record | Scope applied in the prose |
|---|---|---|
| Air stages: 4/4 cycles, 2.466 s, 0.0337 phase error; 7/7, 1.300 s, 0.0933 | `experiments_20261104/GAIT_AIR_VALIDATION.json` and `EXPERIMENT_UPDATE.md` | Same learned weights/expression 0.1, 16 s per run; native-muscle human model; only physical stages 1 and 2 |
| 99.860% BW suspension and zero foot load | Same stage-2 record | Measured suspension accounting from 5 s onward; not minimum ground-support requirement |
| Six noisy trials, 0.01/0.025 nA, three seeds each, 20 ms correlation; exact repeated seed | `experiments_20261104/BiW2LvarX_RG_NOISE_VERIFICATION.json` | Frozen weights and suspended-air gates; separate from parameter perturbations |
| 2 nA/100 ms stimuli; paired differences 1.14/0.16/1.79/6.04 mV | `experiments_20261104/BiW2LvarX_STIMULUS_PAIRED_METRICS.json` and `EXPERIMENT_UPDATE.md` | Maximum response-window differences, not direct synaptic gains; silent and active phases distinguished |
| syn6 0/12 timing passes; variance rank 3 versus intended 6; S5 duties 1.000 versus about 0.325/0.330; s3k rank 2 | `experiments_20261104/circuit_trace_20261105/PF_DESIGN_AUDIT_20261105.md` and family timing JSONs | Offline PF design failure; engineering thresholds; rank is variance description, not neuron count |
| Gain-only mirrored-target repair loss 2.51171 to 2.15643, still 0/12 passes | Same audit | Only the failed acceptance is retained in prose; detailed loss retained here and in source |
| W2L108 8L/7R clear swings; 29.519% upward assistance after 1 s; 35.069% whole-trial signed assistance; peak force 11.727 BW | `experiments_20261104/LiDeng_W2L108_RESULT.md` | Separate native-muscle body, rails/guides retained; belt-relative displacement is not overground progression |
| Slower candidate 7/7 swings, 10.105% assistance, 0.351/0.417 m excursion, 1.929% sustained double support versus 8% target | `experiments_20261104/SHORT_SHUFFLE_RECORD_20261104.md` | Rejected short-shuffle condition despite reduced assistance |
| CAD 28 s; 77.5% measured BW support in 16–28 s; 7R/6L clear intervals; 5.10° knee margin; 5/8 perturbation passes | `records/CAD_CONTACT_HANDOVER_EXPERIMENT_RECORD.md` | Guided stepping approximately in place; reference tensions; not a calibrated 54-BPA or harness-free result |
| CAD contact criteria, 1/5 N confidence thresholds, 10 ms filter, 16.13/12.91 mm transient penetrations | Same record | Controller task and soft-contact limitations retained; no physical damping change claimed |
| Complete-process 199.624/29.284 s = 6.82; neural 179.079/9.780 s = 18.31 | `docs/TOOLCHAIN_EVIDENCE_INVENTORY.md`, acceleration section | Matched historical complete process and its neural portion; not transplanted to other benchmarks |
| Compact warm-kernel increment 1.613; whole diagnostic pair 92.944/91.516 s = 1.016; all 62 arrays exact | `records/NEURAL_ACCELERATION_20261104.md` | Increment versus an already accelerated backend; one instrumented body pair |
| Retained shared drive 2.0773 nA | Existing Final correction and `Code/MuJoCo_SNS/spinal/wean_tune_20261016.py` | Zero additional hand-drive scale is not zero total descending drive |

No bibliography additions are needed. Existing cited keys were checked against the downloaded ZIP's bibliography.

## Figure provenance and bounded inspection

The original full twelve-channel PF plot remains in the handoff. `crop_pf_timing_panels.py` composes representative right-side S1, S2, and S5 vector panels without changing their time coordinates, transmission coordinates, or phase. It replaces the original small legend with a common 10-point legend and adds a 10-point x label and ticks. The derived one-page figure was rendered and inspected. It supplies a readable representative display; the 0/12 verdict comes from the full machine-readable audit.

Three other dated appendix-root copies preserve exact source bytes: matched W2L rhythm-core traces, original shared-RG architecture, and archived AnimatLab aligned recording. Each has an alt-text companion. There are eleven active graphic inclusions versus fourteen in the ZIP baseline, with two inclusions forming one continued stimulus figure.

The full BiW2Lvar+ seven-row stimulus figure has 10-point fonts at its original 597.6-point width. `split_biw2l_stimulus.py` retains all seven rows in two landscape vector strips, removing only the original heading and inter-row whitespace. The first strip shows actual current and RG/PF/MN responses; the continued strip shows hip/knee/ankle responses. Common 10-point time ticks and labels are added to each. Both strips were rendered and inspected; their time and data coordinates are unchanged. Landscape placement keeps the source fonts near their original physical size instead of reducing the seven-row figure into a small portrait page. Raw aligned recordings remain identified on easteregg2 for subsequent replotting if needed.

## Checks and reviewer response

`APPENDIX_STATIC_CHECK.json` reports no unmatched environments, duplicate labels, missing graphics, or missing/new bibliography keys. This is a structural check, not a substitute for the root agent's full dissertation build and rendered-page review.

The independent supervisor found no material air/PF/CAD overstatement and requested five improvements. Applied: clearer flexion-negative knee terminology, wrapped air-table header, conversion/CAD/bridge-test entrypoints in the source map, retained explicit neural-only interpretation for the reduced circuit with joint-limit violations, and removal of the incidental repair-loss values from prose. Figures and page layout remain subject to the final compiled review.
