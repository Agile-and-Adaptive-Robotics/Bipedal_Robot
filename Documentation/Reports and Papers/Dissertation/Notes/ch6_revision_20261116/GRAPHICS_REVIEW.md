# SNS-Simscape graphics revision

The Chapter 6 library and two reflex figures are regenerated from the current library drawing definitions and existing recordings. Original graphic files and all source/result SLX/MAT files are preserved. No model parameter or experimental result was changed, and no simulation was rerun.

## Delivered figures

| Figure in `Figures/60-sns/` | Evidence and scope |
|---|---|
| `sns_library_20261116.pdf` | Twelve actual `SNS_Library.slx` masks corresponding to the block selection on defense slide 69. Muscle glyphs use the corrected shared emitter. This is a library inventory, not a validation result. |
| `sns_knee_reflex_20261116.pdf` | Existing `results/knee_transient_reflex_20261028/replay_summary.mat`; all six actual sensory-to-MN synapses, input pulse, paired knee trajectories, and MN voltages. The graded-ramp configuration gives 1.09079283 deg extension, 111 ms peak latency, and 6.355 s recovery. Its actuator is the simple `BPAForce` element, not the Festo pressure-length model. |
| `sns_cup_reflex_20261116.pdf` | Existing `results/sns_beer_results.mat`; six actual synapses, side-view BPA glyphs, added-mass ramp, paired elevation trajectories, and normalized pressure commands. The force element is `BPA_20mm`. Endpoint sag is 1.99653505 deg with feedback and 7.84641287 deg without feedback. |

Each PDF has a PNG preview and an adjacent `.alt.txt`. Numerical traces are displayed without smoothing, resampling, or filtering. The knee convention is the recorded sign conversion `theta_paper = -theta_internal`. Cup elevation is `-sag`, so a lower curve represents greater droop.

**Loading clarification:** the saved cup trial is 10 s, with an unclamped 0.0625 kg/s ramp. It reaches 0.5 kg at 8 s and **0.625 kg at the 10 s endpoint** associated with the quoted sag values. The builder's “0.5 kg over 8 s” comment does not describe the complete recording. The revised graphic displays the full 10 s. A startup transient remains visible; this is not a settled load-step recovery experiment.

## Topology and provenance checks

`figscripts/extract_sns_graphics_sources.py` reads the current SLX XML without modifying it. Both knee and cup contain exactly six sensory-to-MN synapses: each Ia source excites its own MN and inhibits the antagonist MN, and each Ib source inhibits its own MN. That is two excitatory and four inhibitory connections. The cup graphic does not reuse the inaccurate older v5 diagram that sent biceps Ib to both MNs and omitted triceps Ib.

The simplified pathway panels show sensory-neuron channels and their six synapses. Fixed drive, transducer internals, activation filters, and the mechanical return are omitted from these panels; the underlying models retain them. Excitation uses an open terminal triangle, with its base toward the target, matching the library convention. Inhibition uses a filled terminal circle and dashed line. Crossing lines without a terminal or junction dot are not joined.

`figscripts/sns_graphics_sources.json` preserves masks, connectivity, and SHA256 source hashes. `figscripts/graphics_metrics.json` records endpoint values. `GRAPHICS_QC.json` records final PDF hashes, dimensions, embedded font sizes, absence of raster-image XObjects, and source conservation.

## Display-only helper correction

`Code/Matlab/SNS_Simscape/sns_icon_muscle.m` had three display errors (corrected at lines 20, 27, and 40):

- The outer body used `g.bodyY + edgeSign(...)` although `edgeSign` already returns the expanded y coordinates. It therefore doubled the body height. It now uses `edgeSign(...)` alone.
- Braid rows store `[x1 y1 x2 y2]`; the emitter passed this ordering into `plot([x1 x2],[y1 y2])`. It now reorders to `[1 3 2 4]`.
- The biological-muscle fiber emitter had the same coordinate-order issue and now uses the same correction.

MATLAB R2025b verified all 14 braid segments and the outer-body coordinates against `sns_icon_geom.m` to the emitter's numeric precision (maximum tolerance 6e-5 icon units). No force, neural, timing, unit, or mask-value calculations changed. Existing SLX masks contain expanded commands, so changing the helper does not update those saved masks automatically. `figscripts/refresh_sns_library_icons.m` produces a separate corrected library copy under `tmp/sns_library_icons_20261116/`, leaving the original library untouched. Its runtime execution remains to be checked outside this restricted session.

## Execution and review limits

MATLAB **25.2.0.2998904 (R2025b)** successfully read the MAT recordings, exported numerical arrays, evaluated the actual mask drawing commands into vector primitives, and passed the geometry check. Its graphics initialization then failed in this execution sandbox (`hgrc` screen initialization followed by a graphics-driver access violation). Simulink also failed to start its localhost helper. These failures occurred before any source-model change; they are not evidence of failed model dynamics.

The delivered figures therefore use the standard Matplotlib Agg vector backend on the exact MATLAB-exported arrays and glyph primitives. This fallback is explicit rather than described as a successful MATLAB graphics export. The full MATLAB renderer is retained as `figscripts/sns_graphics_20261116.m` for execution in a normal MATLAB session. The numerical bridge is `figscripts/export_sns_numeric_data.m`; its generated data are temporary, not added to the repository.

All three final PNGs were inspected. Library labels and descriptions fit; six-edge diagrams have no hidden sensory sources; knee panels have separate input, motion, and MN axes; the cup panel shows the complete load history. Quantitative traces use Tol blue/red with solid/dashed differentiation. Figure text is at least 10 pt at native 14 cm width; final dissertation sizing must preserve that floor. Final compiled-page inspection remains the coordinating editor's responsibility.

An independent reviewer inspected the three final graphics and found the pathways correct, labels readable, and panels free of clipping or material overlap. MATLAB Code Analyzer reported no syntax issues: the isolated library refresher has zero messages, and the numeric icon-capture script has five advisory messages about its temporary global capture state. An obsolete suppression in the native renderer was removed after inspection.

## Reproduction

1. Run `figscripts/extract_sns_graphics_sources.py` from the repository using Python.
2. In MATLAB, add the `figscripts` directory and call `sns_graphics_20261116` for native MATLAB rendering.
3. If the session has the same restricted graphics problem, call `export_sns_numeric_data` in MATLAB, then run `figscripts/render_sns_graphics.py` with Python. This is the path used for the delivered PDFs.
4. To preview corrected saved library masks separately, call `refresh_sns_library_icons`; inspect the returned copy before replacing any working library.
