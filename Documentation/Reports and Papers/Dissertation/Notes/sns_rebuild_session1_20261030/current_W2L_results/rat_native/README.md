# Public rat model: native subnetwork figure

Machine: **easteregg2**. Figure ready at:

`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results\rat_native\DENG_RAT_NATIVE_SUBNETWORKS.pdf`

SVG and PNG with the same stem are alongside; `CAPTION.tex` supplies a proposed caption. The **final figure is reflowed on a 9 × 5.5 inch canvas with 10 pt minimum type**, fitting the ETD sideways text block plus caption. The body panel is rasterized from the static imported mechanics; the actual native neural drawings remain vector. The final PNG was visually inspected.

Source: `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\Biarticular_net.aproj`, imported from `https://github.com/vipamol/Rat-SNS-design-and-tuning/tree/main/Biarticular_net`. The public source is a **bilateral hindlimb rat precursor**, with 120 neurons and 16 physical muscles. It is not identified as the final 2023 published configuration. The source pelvis is frozen; this figure is a model/circuit description, not a new ground-walking result.

The figure reproduces actual native pages `LH`, `LH_RG`, `LH_Hip`, and `Biar`. The last page contains the source BFP and RF coupling interneurons and their real hip/knee PF inputs; it is not replaced by a newly invented four-cell circuit. Native node centers and internal route control points are preserved. Neural and most off-page glyphs are enlarged for publication type; original names, shortened labels, source IDs, and hashes are in `PROVENANCE.json`. Existing source cross-side connections and asymmetries are not repaired or hidden by a newly synthesized circuit. Detailed right-side pages remain in the source project.

The static body panel uses `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\Deng_Biarticular_Rat.xml`, the already compiled MuJoCo mechanical import. Only a publication camera, limb colors, visible site groups, and line width change in memory. There is **no integration step**, controller fitting, source-file write, or physics change. It depicts the imported zero-coordinate source pose, with the author's fixed pelvis preserved. Neural data for this precursor remain the native recordings documented separately in the model's README and `EXPERIMENT_RECORD.md`.

`C:\Users\Ben Bolen\Pictures\Saved Pictures\Bipedal_Walker.png` was inspected: it is the CAD humanoid with muscle routes in anterior/oblique views, **not this rat model**. That identity must remain distinct in the dissertation.

Generator: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results\make_rat_native_figure.py`; Python: `D:\Anaconda\envs\myo\python.exe`. No manuscript, Overleaf project, native model, or original result was modified.
