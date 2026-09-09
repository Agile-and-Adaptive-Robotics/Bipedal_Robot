# Neuromechanical figure integration - 2026-09-09

The preliminary workflows section is now after the historical AnimatLab walker
section, at the end of Methods. The actuator-force section remains the first
technical Methods section after the overview.

The figure sequence includes the AnimatLab body/hierarchy, controller reference
schematics, real working-project GUI captures, native MuJoCo model views,
Simulink library/neuron implementation, and the complete Simulink reflex topology.
Results uses the new bilateral hip/knee/ankle plot. The earlier neural/contact
diagnostic PDF is retained as a supporting asset and is no longer the included
Results figure.

## Provenance and limits

- Nine user-supplied images were copied unchanged from the explicitly supplied
  temporary paths into `ProofFinal/figs/Preliminary/AnimatLab_reference/`.
  The body, hierarchy, RG, PF, and sensorimotor images are used as architecture
  illustrations. Their precise historical dates were not established. The
  supplied joint and CPG plots are not relabeled as the new phase-1 run.
  `footMechano.png` and `fullNetwork.PNG` are retained as reference material.
- `animatlab_body_gui.png` and `animatlab_network_gui.png` are genuine window
  captures of `Biped_2xCPG_wSubs.aproj` opened in AnimatLab Pro. No simulation
  was run or project saved. They show the project GUI, not a captured phase-1
  time point. The phase-1 numerical export comes from the separately instrumented
  `.asim` named in `preliminary_animatlab_verification.json`.
- MuJoCo model figures and their exact view settings are reproducible with
  `export_mujoco_models.py`; render provenance is alongside the figure assets.
  The converted robot is at keyframe 0. The separate knee is statically posed
  at -35 degrees for visibility. Neither is a new dynamic result.
  An actual MuJoCo viewer was also launched at keyframe 0, but the installed
  desktop capture tool hung. No viewer screenshot was saved. The viewer exited
  cleanly without stepping; the included MuJoCo figures are native renderer
  outputs, not GUI screenshots.
- Simulink uses native exports with display-only layout changes. See
  `simulink_figure_provenance.md` and the source-hash verification. No model
  definition, scientific parameters, or port connectivity was changed.
- `build_preliminary_animatlab_figure.py` checks complete finite aligned
  recordings and excludes nine padded rows beyond 10 seconds. It generates
  both diagnostic and joint-motion figures from the saved data. Joint values
  are the recorded `JointRotationDeg` channels without sign reversal or
  smoothing. Contact-sensory exports are voltages, not force measurements.

The five original claim limits remain explicit: no converted-robot SNS control,
no validated Xi-corrected MuJoCo torque, no completed CAD-coupled Simscape
Multibody plant, no stable AnimatLab walking or proven contact gating, and no
quantitatively established Simulink knee regulation from the saved run.

## Review and synchronization

`build_neuromechanical_review.py` creates a review packet with captions and
native application exports. Simulink's high-resolution PNG exports are included
because its PDF printouts retain excess paper margins and unembedded-font
dependencies; original vector exports remain available beside them.
This is not a LaTeX compilation or an Overleaf layout
proof. Local edits are in `ProofFinal/chapters/20-methods.tex` and
`30-results.tex`; referenced new assets are below `figs/Preliminary/`.
The Overleaf project, dissertation ZIP, and existing full dissertation PDF
have not been updated by this figure-integration work. No git commit or push.

Before overwriting the online chapters, compare them with the current Overleaf
source to preserve edits made elsewhere. Upload the referenced assets with
their relative paths and recompile there to verify final float pagination.
