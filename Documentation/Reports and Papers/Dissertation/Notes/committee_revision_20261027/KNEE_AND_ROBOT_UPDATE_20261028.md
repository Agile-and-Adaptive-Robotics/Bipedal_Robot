# Knee transient and Session 1 robot results

Completed in the established ProofFinal and Figures folders, and uploaded to the existing dissertation Overleaf project. No commits or pushes. No old ZIP extracted over working edits.

## Knee demonstration

Replaced the startup trace in Figure 6.4(A) with the verified transient sensory-current response. Panels B and C retain the model illustration and six implemented feedback connections. The illustrative knee pose is identified as such. The old composite is preserved in Figures/60-sns/Deprecated and in this pass's Archive.

The new result is an explicitly modified engineering configuration. Independent comparison of the loaded original and final saved block parameters confirms exactly six saturation-ramp changes: both MN widths 1 to 10 mV; four Ia synaptic widths 0.5 to 10 mV. The original configuration maintained oscillations, so its old startup trace cannot substantiate recovery from a settled baseline. No torque-model identification term was changed.

Independent calculation from the saved complete control/pulse/fine arrays reproduces:

- Positive extension excursion 1.0907928305 degrees, peaking 0.111 seconds after the 2 nA, 30 ms extensor-Ia current pulse.
- Recovery at 6.355 seconds under simultaneous 0.01-degree displacement and 1e-4 rad/s velocity thresholds.
- Maximum angle difference 5.2027343e-5 degrees when the maximum ODE45 step is halved.
- Exact control/pulse equality before stimulus onset.

The supplied fresh-process replay has zero difference in its complete resampled pulse arrays. This pass independently checked those saved arrays and parameter changes; it did not rerun the MATLAB simulation. Original integrity records, failure traces, successful comparison, fresh replay, solver checks, and raw MAT files are preserved together.

Chapter 6 reports the response and points to the new Appendix D subsection, Transient Knee Sensory-Current Response. Appendix D specifies the pulse location, settling, ramp changes, simple actuator law, inactive baseline flexor, recovery thresholds, and solver checks. Figure D4 shows stimulus, angle, velocity, MN voltages, activations, and forces. The current perturbation is not a tendon impact, physiological calibration, or proof of suppression of a tonically active antagonist. The old 44.4/43.0-degree startup comparison is superseded in the main narrative, with its previous source preserved in Archive.

## Session 1 robot integration

Read the current remote CHATGPT_REPORT.md Session 1 block through its final foot-mount/w2lvar continuation, FOOT_MOUNT_W2LVAR_RESULTS.md, media placement/provenance records, recorded assessments, and backend verification. The captured report is preserved separately; it was not copied over the local report.

Figure 6.9 adds the accepted anterior/oblique CAD-model still. Its SHA256 matches the supplied provenance. The existing Appendix D CAD subsection now adds the matched 40-second separation-control comparison (Table D8), while keeping the earlier 20-second/94.68%-support stance result separately identified.

The new comparison uses 80% vertical assistance over 15–40 seconds in both conditions. Added muscle-mediated separation feedback and a 6% abducting-route activation floor change reversed ankle ordering from 67% to zero and minimum separation from -38.3 to 142.4 mm. These values agree with the saved assessments. The supplied full-trajectory reference/accelerated verification is bit-exact; this pass reviewed that evidence rather than rerunning the robot simulation.

The important control limitation is explicit: compensated MN target currents cancel incoming PF/reflex synaptic currents. Effective control combines RG timing, engineered motion targets and force allocation, and the support rig. This is not full PF-driven muscle-pattern control. The 262-neuron baseline and separate 544-neuron/230-input W2L-organization timing probe are distinguished; four named ankle aliases are not neurons.

The provisional welded foot datum and ideal offset route brackets are documented. The separate 60 mm proximal-foot hypothesis changes median spacing to 57.3 mm and introduces 1.316 mm sampled penetration under the baseline hybrid controller. Its shifted w2lvar timing comparison passes unloading/lifting counts but still fails alternating-apex and progression tests. No attachment hypothesis is represented as a verified fabrication correction. No new gait-success claim is added.

## Verification and delivery

- Edited sources: 60-sns.tex and 95-AppendixD.tex only. Before/after copies and diffs are in Archive/knee_and_robot_update_20261028.
- Five-pass fresh delivery build: all return codes zero, no errors or undefined references, 319 pages. The existing 0.669-point protected Appendix C overflow remains.
- Rendered and inspected final Figure 6.4, Figure 6.9, Figure D4, and Table D8. A cropped knee label was caught and repaired before the delivery build. The full transient plot was arranged in two columns to improve label size.
- Overleaf final compile: 320 pages, Errors 0 / Warnings 0. Both edited sources were read back from the live editor and match the local sources after line-ending normalization. Active figures were uploaded into figs/60-sns and figs/95-AppendixD.
- Root Bolen_Dissertation.pdf/.zip updated; 89 active graphics packaged. ZIP integrity and source/support byte comparisons pass.
- Published Chapter 20, torque Chapter 30, Appendix C, and torque-paper Main.tex match their pre-pass hashes. Biomimetic results and torque curves remain untouched.
- Active files: Figures/60-sns/knee_transient_current_layout.pdf and ah001_supported_control.png; Figures/95-AppendixD/knee_transient_response.pdf. Generators, provenance, and alternative text are in the corresponding Components folders; complete incoming results remain in Archive.

Directive evidence: Chapter 6 Coupled Block-Model Tests and Figure 6.4; Chapter 6 Repeatability and Mechanical Support and Figure 6.9; Appendix D Transient Knee Sensory-Current Response (source lines 58–72); Appendix D CAD-Derived Robot Integration, matched comparison and table (source lines 335–360). Independent array/parameter checks, source readbacks, protected hashes, visual proofs, and delivery verification are in the pass Archive.

Final accessibility check: the reused knee illustration carried hidden, clipped text from the discarded startup plot. Its visible illustration was embedded as a high-resolution crop, while the new plot and circuit remain vector graphics. Text extraction from the active figure and final delivery page confirms the old startup/NumPy labels are absent. Fresh five-pass build after this cleanup passes with unchanged page count and results.
