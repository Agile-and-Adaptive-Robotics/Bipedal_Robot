# Final independent CAD and figure review

The revised canonical `Final/chapters/60-sns.tex` and `95-AppendixD.tex` were read back after integration. This is a bounded source, saved-array, vector-font, and image review. It does not assert that a full dissertation build or compiled-page inspection succeeded.

## Evidence and prose closure

| Check | Final source evidence | Verdict |
|---|---|---|
| CAD remains separate from the human-controller curriculum | Chapter 6, lines 220–229 and 267 | Pass: reference tensions, engineering allocation, retained guidance, and limited local outcome remain explicit |
| Support denominator and scoring window | Appendix D, line 334 | Pass: 570.20625 N / 735.75 N, 75 kg, and `16 <= t < 28 s` identified |
| Right/left knee identity and sampling | Chapter 6, line 220; Appendix D, line 334 | Pass: 11.20° right / 5.10° left, whole-trial online 2 ms bound distinct from the saved 10 ms trace |
| Acceptance criteria | Appendix D, line 334 | Pass: 5 mm mesh height, 7.36 N ipsilateral unloading, 15 N contralateral loading, and 80 ms duration |
| Exact saved repeat comparisons | Appendix D, line 345 | Pass: all 26 arrays; current source/geometry checks kept distinct from unavailable historical compiled binary |
| Local perturbation scope | Chapter 6, line 229; Appendix D, line 345 | Pass: five of eight; predictive-knee failures; no global-stability or all-connection claim |
| Soft contacts and actuation fidelity | Appendix D, lines 332–347 | Pass: retained penetrations and ideal reference tensions are explicit |
| New CAD trace | Appendix D, lines 338–340 | Pass: native signs, units, simultaneous load/clearance, threshold mapping, and all 1,200 final-window samples |
| Knee figure mapping | Chapter 6, lines 166–168 | Pass: A circuit, B pulse, C angle, D MN; simplified actuation retained |
| Cup figure mapping | Chapter 6, lines 175–177 | Pass: A circuit, B geometry, C load, D elevation, E normalized pressure command; exact plotted elevation is the negative of saved positive sag |
| No unsupported-walking claim | Chapter 6, line 267; Appendix D, line 313 | Pass |

Line numbers identify the readback at this review and can move during final layout corrections. Source hashes and all CAD calculations are in `CAD_AUDIT_RESULTS.json`; figure input/output bindings are in `CAD_FIGURE_PROVENANCE.json` and `CAD_SNS_FIGURE_PROVENANCE.json`.

No material CAD evidence or interpretation error remains in the reviewed revision. The minimum-knee claim remains explicitly dependent on the retained online accumulator because the entire 2 ms trajectory was not saved. The missing historical compiled CAD binary remains disclosed.

## Vector font and scaling check

`CAD_FIGURE_QC.py` reads the actual canonical graphics options and embedded PDF text. The class delegates to installed `book.cls`/`bk12.clo`, giving a source-derived text width of 468 PDF points and text height of 546.4508 PDF points. The font sizes below include the specified manuscript scaling, including limiting heights and rotated figures.

| Active vector figure | Minimum scaled font size (PDF pt) | Figure height (cm) |
|---|---:|---:|
| Workflow | 10.50 | 14.00 |
| SNS block library | 10.00 | 15.20 |
| Compact knee reflex | 10.00 | 15.00 |
| Compact cup comparison | 10.00 | 15.00 |
| Aligned air-stage trace | 10.00 | 15.00 |
| W2L feedback comparison | 11.61 | 14.27 |
| Neural stimulus responses | 10.20 | 15.00 |
| Mechanical stimulus responses | 10.20 | 12.50 |
| PF timing | 10.00 | 11.70 |
| Historical IK overlay | 11.00 | 15.11 |
| New CAD trace | 10.00 | 15.00 |
| Graded/spiking coupled responses | 10.71 | 15.42 |

All active vector figures pass the 10 pt floor. The script records unrounded values in `CAD_FIGURE_QC.json`. The source-derived scaling is an independent check, not a claim that the complete manuscript compiled.

The inherited force-sequence, appendix knee, W2L multiperturbation, and full shared-RG bitmap displays that failed readability were removed from active inclusions by the coordinating editor, while their useful results and original files remain preserved. The retained W2L cross-engine and archived AnimatLab response rasters were visually inspected: their plot labels are readable at the specified size. Embedded-font metadata cannot independently certify the text size in a raster. The motivating model/fall images and CAD still are illustrations rather than substitutes for quantitative traces.

## Compact knee and cup figures

The prior 18.7 and 19.1 cm heights left inadequate caption space in the approximately 19.28 cm text area. The new figures use 14 × 15 cm canvases with unchanged 10 pt minimum text. Panel captions remain A–D and A–E. Redundant descriptive subplot titles were replaced by panel letters, diagram labels were separated, and geometry was repositioned to remove overlaps. All six circuit pathways and their excitation/inhibition signs remain. No numerical sample, waveform amplitude, or timing was changed.

The figures were regenerated from `TIMING_trace_arrays.npz`, whose reproducible MAT reader and original-source hashes are retained in the timing audit. Knee angle is converted from the saved internal positive-flexion coordinate to the dissertation's flexion-negative coordinate. Cup elevation is the negative of saved positive sag. The pulse-relative knee clock subtracts 30 s; it does not shift the response relative to its pulse. Cup series retain each saved native time array without interpolation. No smoothing, added resampling, or decimation was applied. The current shared BPA glyph expressions are retained from `sns_icon_geom.m` and hashed.

The final compact PNGs were visually inspected after correcting title/axis overlaps. Diagram labels, legends, physical units, panel letters, negative knee coordinates, and cup endpoint geometry are readable. The selected figure files and companion alt text are now active in the canonical chapter folders. No further generator run is needed for this review.

## Remaining publication check

Actual full-page review is contingent on a successful canonical build. The new CAD caption is relatively long for the 4.28 cm remaining beneath its figure; the graded/spiking figure retains approximately 3.86 cm for its caption. These must be checked in the actual compiled pages. If a float is oversized, shorten repeated caption explanation while retaining it in adjacent prose, rather than shrinking the figure below the required font size. The coordinating build record owns the compilation verdict and exact failure or success evidence.

No model parameter, actuator law, route, scientific result array, controller configuration, or original generator was changed by the final review. No simulation, optimization, retraining, commit, push, or upload was performed.
