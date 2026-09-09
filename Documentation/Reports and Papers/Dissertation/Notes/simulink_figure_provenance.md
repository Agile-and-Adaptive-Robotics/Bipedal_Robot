# Preliminary Simulink figure provenance

Generated 2026-09-09 with MATLAB R2025b (25.2.0.2998904) using
`export_preliminary_simulink_figures.m` in MATLAB `-batch` mode. Inputs are the
current saved `Code/Matlab/SNS_Simscape/SNS_Library.slx` and
`Code/Matlab/SNS_Simscape/KneeReflexDemo.slx`. No simulation was run and neither
source model was saved. SHA-256 comparison in `simulink_model_hash_verification.json`
confirms both files remained byte-identical to their pre-export versions.

Outputs in `ProofFinal/figs/Preliminary/` are native Simulink diagram exports,
each in vector PDF and 300-dpi PNG:

- `simulink_sns_library_grid`: recommended compact library figure. The same seven
  saved masked blocks were repositioned into three columns in memory solely for
  export; their definitions and parameters were unchanged.
- `simulink_sns_library`: original saved single-column library layout.
- `simulink_knee_reflex`: original saved top-level knee-reflex topology. Its saved
  layout has overlapping labels and blocks; this is an implementation overview,
  not a polished circuit schematic. Several primitive icons render blank before
  model initialization, while block labels and connectivity remain visible.
- `simulink_knee_reflex_arranged`: recommended full topology export (8167 by
  5954 pixels and vector PDF). The same saved model was opened through the
  Simulink API and automatically arranged in memory. This removes the saved
  block overlaps and renders primitive icons correctly. An explicit comparison
  of every input port's source block and source port before/after arrangement
  passed. No block definitions, mask names, parameters, or connections changed.
  Use a dedicated landscape page or a full-resolution supplemental view because
  the complete flat model contains many blocks.
- `simulink_neuron_detail`: interior of the saved non-spiking neuron block,
  opened through `open_system(path, 'force')` before export so the primitive
  icons render correctly. The final visual check confirmed Vrest/Thr constants,
  the 1/s integrator, gains, and saturation icon. Gain expressions too long for
  the icon use Simulink's standard generic gain symbol. No initialization,
  simulation, parameter edits, or model saves were required.

All PNGs were visually inspected. The compact library has readable labels and no
overlaps; the arranged topology removes the original block overlaps. These are application-generated model diagrams, not fabricated GUI
screenshots. They establish implemented blocks/topology only. BPA force and
reduced-order plant parameters remain placeholders; no claim is made about
quantitatively valid knee regulation or a completed CAD-coupled Simscape
Multibody plant.

The successful final export log records the connectivity check and explicit
discard of all in-memory display changes. Both saved source hashes were checked
again after the final export and remained unchanged.
