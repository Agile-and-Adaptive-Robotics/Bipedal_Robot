# Figure 6.1, controller history, and conversion details

Ben corrected the historical identity of Figure 6.1(B). The prior draft incorrectly used the later controller hierarchy beside a description of the original Connor–Ben model. The exact supplied attachment is now the B asset, preserved byte-for-byte. Labels (A) and (B) are added in LaTeX. No circuit connections or raster content were redrawn.

The original model has a shared four-neuron RG with bilateral PF groups and MN-level afferent feedback. Years later Ben split the RG into one per side and added RG/PF sensory feedback and foot-contact feedback. The text and caption distinguish those stages explicitly. The claim that Ben constructed the right-side subsystems only during that later extension was removed.

The main chapter already contained specific MyoConverter problems. Before this correction, Appendix A also discussed force comparisons and rectus-femoris routing, and Appendix C discussed robot geometry; there was no single consolidated repair account in Appendix D. New §D.1 supplies it: hip-direction corrections, retained knee convention, six disabled entries within 92, restoration of eight active-force scales, coupled-path limits versus 14 anatomical driver hinges, rectus-femoris patella routing, and the 19-to-six contact-pair restriction. These numerical additions are traceable to the implementation in conversion_source_evidence.json; no new simulation measurements were generated.

The main text was also corrected to say that anatomical joint limits were re-enabled separately from coupled path-point limits. Passive torques were not presented as the sole replacement for those limits. The driver-joint and follower-joint distinction prevents an ambiguous account of “removing joint limits.”

## Per-directive evidence

- **Editable source**, `chapters/60-sns.tex:1`: The complete editable chapter and supporting TeX project remain in Revision_E.
- **Original versus later controller**, `chapters/60-sns.tex:12`: Original shared RG and Ben’s later per-side RG/feedback extension are explicitly separated.
- **Figure 6.1 panels**, `chapters/60-sns.tex:23`: A and B labels added in TeX; B replaced with Ben’s exact attachment; old image retained but unused in Figure 6.1.
- **Section heading**, `chapters/60-sns.tex:63`: Heading changed verbatim to Block Programming Neural Pathways.
- **MyoConverter summary**, `chapters/60-sns.tex:43`: Retained readable main-text problems and added appendix cross-reference.
- **MyoConverter details**, `chapters/95-AppendixD.tex:6`: Consolidated source-supported repairs in Appendix D.1; existing Appendix A force comparisons and Appendix C robot geometry details remain.

## Build and preservation

Fresh build `build_20261026_051358`: 270 pages, zero errors, zero unresolved references/citations, and zero missing figures. The 26 document-wide overflow warnings remain. No oversized-float warning was produced. Rendered Figure 6.1 (printed p.121) and new Appendix D.1 (printed p.229 onward) were checked. Output: `Revision_E_corrected.pdf`. The earlier Revision_E_local.pdf was open and locked, so it was retained and the correction saved under this new filename.

Exact prior sources are preserved in before_60-sns.tex and before_95-AppendixD.tex. changes.diff and conservation.json record the changes. Old figure assets and all original labels remain. Diagram density and the supplied raster resolution still limit the smallest text; the attachment itself has not been altered. Overleaf, repository model code, and remote machine files were not edited. No commit or push was made.
