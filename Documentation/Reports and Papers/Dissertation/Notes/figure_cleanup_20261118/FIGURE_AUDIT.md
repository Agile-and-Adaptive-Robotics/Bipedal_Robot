# Dissertation figure-folder audit, November 18, 2026

The working manuscript uses 30 graphic files for 27 numbered figures in Chapters 6 and 7 and Appendices A–E. Figure 6.1 and Figure 7.2 each use two source files. Appendix D places the neural and mechanical panels of Figure D3 on consecutive pages. The chapter-folder roots contain the exact graphic inputs and their available accessibility descriptions.

The [figure transfer ZIP](Revised_chapter_figures.zip) contains the 30 graphics and 25 descriptions under Overleaf’s existing `figs/` paths. Its source is the local revised TeX. Compared with the existing project export, 26 graphics match byte for byte and four PDFs have newer artwork under the same names.

| Figure | PDFs to replace |
|---|---|
| 7.2 | `70-overall/bpa_construction.pdf` and `70-overall/bpa_exploded_cutaway.pdf` |
| C4 and C5 | `94-AppendixC/extensor_cam_AF_bones.pdf` and `94-AppendixC/extensor_cam_GH_bones.pdf` |

Appendix A uses equations and a parameter table. Appendix E points to the citation map in Figure 5.1 and the selected SADb record in Figure 5.3. Its former `sadb_explorer.png` was the one unused standalone figure at a chapter-folder root.

## Printed figures and source lines

Page numbers are printed manuscript pages. Lines refer to the chapter's TeX in `Final/chapters`. The [compiled figure inventory](COMPILED_FIGURES.json) also records PDF page numbers, source callouts, and the transfer archive’s SHA-256.

| Figure | Page | Subject | Active graphic | TeX line |
|---|---:|---|---|---:|
| 6.1 | 134 | Original controller transfer and fall on release | `OurWalker.png` and `LI_DENG_FALL_FRAME.png` | 19, 23 |
| 6.2 | 136 | Implemented modeling and testing workflow | `simulation_workflow_authoritative.pdf` | 41 |
| 6.3 | 141 | Editable SNS-Simscape blocks | `sns_library_20261116.pdf` | 102 |
| 6.4 | 147 | Transient knee reflex | `sns_knee_reflex_authoritative.pdf` | 168 |
| 6.5 | 148 | Loaded-cup feedback | `sns_cup_reflex_authoritative.pdf` | 177 |
| 6.6 | 150 | BiW2Lvar+ in the afferented air stage | `biw2lvar_air_aligned_authoritative.pdf` | 188 |
| 6.7 | 152 | W2L sensory-pathway removal | `w2l_feedback_short_20261031.pdf` | 205 |
| 6.8 | 155 | CAD robot under supported hybrid control | `cad_supported_20261116.png` | 227 |
| 7.1 | 161 | BPA joint design workflow | `design_workflow.pdf` | 12 |
| 7.2 | 164 | BPA end transition and proposed loading experiments | `bpa_construction.pdf` and `bpa_exploded_cutaway.pdf` | 33, 34 |
| 7.3 | 170 | Proposed additions to the control stack | `sensorium_hierarchy.pdf` | 77 |
| 7.4 | 171 | Simulation-to-hardware workflow | `pipeline.pdf` | 88 |
| B1 | 215 | AnimatLab spiking-conversion variants | `animatlab_spiking_conversion.pdf` | 194 |
| B2 | 218 | W2L population connectivity | `w2l_connectome_current.pdf` | 213 |
| B3 | 219 | W2L synaptic conductances | `w2l_conductances_current.pdf` | 222 |
| C1 | 235 | Extensor route on the bone meshes | `bone_context_extensor20mm_20261030.png` | 325 |
| C2 | 236 | Direct-route contact and release | `direct_route_contact_detail_20261030.pdf` | 333 |
| C3 | 242 | Torque at five saved pressure levels | `pulley_torque_isobars_20261030.pdf` | 407 |
| C4 | 243 | Cam-extensor poses A–F | `extensor_cam_AF_bones.pdf` | 417 |
| C5 | 244 | Cam-extensor poses G and H and release-angle key | `extensor_cam_GH_bones.pdf` | 424 |
| D1 | 256 | Matched-input W2L rhythm-core comparison | `w2l_crossengine_20261116.png` | 108 |
| D2 | 263 | Archived AnimatLab neural and mechanical response | `animatlab_short_response_20261116.png` | 187 |
| D3 | 267 | Neural response to active-phase interventions | `AIR_stimulus_neural.pdf` | 232 |
| D3 continued | 268 | Hip, knee, and ankle response in the same intervention trials | `AIR_stimulus_mechanics.pdf` | 240 |
| D4 | 270 | Measured syn6 PF transmission versus intended primitives | `PF_timing_authoritative.pdf` | 253 |
| D5 | 275 | Bilateral joint angles and stance against Gait2392 inverse kinematics | `gait2392_ik_overlay.pdf` | 314 |
| D6 | 280 | Supported CAD model traces | `CAD_supported_traces_20261116.pdf` | 355 |
| D7 | 284 | Graded and spiking responses in matched ground trials | `spiking_ground_neuromechanics_short_20261031.pdf` | 400 |

## Archived files and verification

Fifteen files moved out of chapter roots. Fourteen were alternate PNG or SVG exports of figures included as PDFs, and one was the Appendix E screenshot. Thirty-five additional files belonged to earlier layouts, previews, and their companions in `Components`. All 50 moved to their owning chapter’s `Deprecated/20261118/`. Thirty-seven earlier files found in the existing incoming project archive also have byte-matched copies there. Three accessibility descriptions for the current Appendix B graphics were recovered from the earlier archive and copied beside their PDFs.

The [file manifest](FIGURE_FILES.json) records original paths, current paths, reasons, and SHA-256 digests. The [retired artwork and companions](Archived_figure_files.zip) preserve all 87 newly archived files under their repository-relative dissertation paths. This ZIP carries them through Git because the repository ignores new files inside `Deprecated`. Source arrays, original exports, and the November 17 manuscript delivery remain available.

MiKTeX compiled `Final/main.tex` with BibTeX in a fresh directory. All 313 pages match the previous delivered proof in extracted text. All 28 figure placements were inspected and rendered pixel for pixel against that proof at 144 dpi. The recorder lists each of the 30 expected graphic files directly under the correct chapter. The final logs contain zero LaTeX warnings, missing graphic files, undefined references or citations, overfull boxes, and BibTeX warnings.
