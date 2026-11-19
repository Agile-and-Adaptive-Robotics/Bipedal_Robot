# Workflow arrows and the illustrated BPA cutaway

The previous manuscript placed the Chapter 7 hierarchy at 6.24 inches wide, where its arrowheads measured 3.84 points in both directions. The revised filled heads print 9.60 points long and 7.68 points wide. All six workflow figures now use the same design for signal direction. The archived exploded BPA illustration is active again, and the axial profile and the three loading experiments retain their own full-size panels.

## Printed workflow figures

The proof supplied for the review is `tmp/ch6_advisor_20261119/proof_r9/main.pdf`. The table uses its PDF page placement matrices and its figure numbering. Arrow measurements are **length along the path × full width across it**, in printed PostScript points (72/in). The original open chevrons were measured between their rear vertices and the tip, and the original circuit-symbol arrows were measured on a 720 dpi raster. The Appendix D source outlined its triangle with a 1.5-point stroke. The table gives the underlying vertices. [Print measurements and file hashes](FIGURE_PRINT_QA.json)

| Earlier figure / PDF page | Graphic | Earlier head (pt) | Revised head (pt) | Smallest printed type (pt) | New tips |
| --- | --- | --- | --- | ---: | ---: |
| 6.2 / 157 | Two model routes | 5.2 × 6.1 | 10.04 × 8.03 | 10.46 | 4 |
| 6.6 / 169 | Controller retest and optimization | 5.3 × 6.2 horizontal and 6.3 × 7.3 vertical. The circuit-symbol heads were about 2.5–2.7 × 4.8–5.0. | 10.04 × 8.03 | 11.05 | 8 |
| 7.1 / 184 | Design workflow | 4.0 × 4.0 | 10.00 × 8.00 | 10.70 | 7 |
| 7.3 / 193 | Sensorium hierarchy | 3.84 × 3.84 | 9.60 × 7.68 | 10.37 | 8 |
| 7.4 / 194 | Simulation-to-hardware pipeline | 4.0 × 4.0 | 10.00 × 8.00 | 10.50 | 10 |
| D8 / 310 | Research toolchain | 4.91 × 4.91 | 10.24 × 8.19 | 10.24 | 5 |

The original Chapter 7 PowerPoint generator specified a 1 pt shaft and the smallest triangular end available in Office. Revised workflow shafts are 1.54–1.64 pt at the six inclusion widths. White space around the larger tips keeps the direction readable at each receiver. In the hierarchy, both proposed descending signals enter the CPG independently. The outer body-orientation signal now enters the vestibular/ocular box from its outside edge. In the hardware diagram the joint-controller return and the proposed transfer from SNS-Toolbox enter the Orin through separate ports. The mechanical conversion, simulation feedback, spinal loop, and proposed hardware paths retain their original endpoints. Dashed outlines and paths preserve the proposed extensions. [Editable vector generator](../../Figures/70-overall/Components/workflow_review_20261119/revise_workflow_figures.py)

All 42 vector tips were counted in the revised PDFs. Extracted words and their multiplicities match each pre-edit diagram after normalizing case and punctuation. I examined every revised diagram at manuscript scale, including the four arrows inside Figure 6.6’s two neural examples and each return path in the two Chapter 7 hierarchies. The fonts and network drawings in Chapter 6 and Appendix D remain those of the original graphics. All three redrawn Chapter 7 diagrams embed Arial. Re-running the generator reproduces the six active PDF files. [Print-size proof](../../../../../tmp/toolchain_sadb_20261119/figure_readability/FIGURE_PRINT_PROOF.pdf)

## Restored cutaway

The 468 × 246 pt file under `Figures/70-overall/Deprecated/20261117` is the original labeled half-section. The retained fitted-end illustration contains the fabric covering, outer and inner rubber, two oppositely wound aramid fiber families, and the translucent longitudinal plane. Its editable PowerPoint file places labels on a prompted raster illustration. The drawing describes a possible model geometry. Specimen measurements would establish the wall thicknesses, angles, and bonded or sliding interfaces. The cut plane identifies the illustration view. A symmetry boundary for an actual reinforcement pattern would require its own justification. [Illustration](../../Figures/70-overall/bpa_exploded_cutaway.pdf) · [Original labeling source](../../Figures/70-overall/Components/bpa_exploded_cutaway/bpa_exploded_cutaway.pptx) · [Source image prompt](../../Figures/70-overall/Components/bpa_exploded_cutaway/bpa_exploded_cutaway.prompt.txt)

| Active file | Disposition |
| --- | --- |
| `bpa_construction.pdf` | The current axial diameter profile remains panel A at 14 cm. |
| `bpa_loading_conditions.pdf` | The current panel B with free, loaded, and fixed-length tests was copied byte-for-byte under this descriptive name, with its alt text. |
| `bpa_exploded_cutaway.pdf` | Restored the original PDF byte-for-byte, SHA-256 `f7075e81d7f3c4fdb62959a7743718d669b3afe8d7eb7c029f41e62fc16b7da0`. Alt text describes what the source illustration supports. |
| `design_workflow.pdf`, `sensorium_hierarchy.pdf`, and `pipeline.pdf` | Chapter 7 vectors were redrawn with the same labels and signal destinations. Each now has active alt text. |
| `simulation_workflow_authoritative.pdf`, `controller_retest_and_optimization.pdf`, and `research_toolchain_workflow.pdf` | Replaced the directional glyphs in the original artwork. The separate masks follow the narrow corridors between the existing nodes. |

**Placement supplied to the manuscript editor.** In the A/B experimental float, change only the second included filename from `bpa_exploded_cutaway.pdf` to `bpa_loading_conditions.pdf`. Both graphics use `width=14cm`, and both retain their original letters. Introduce the material cutaway where Festo’s rubber wall and aramid fabric are discussed, and place the cutaway in its own float immediately after the experimental figure. Its archival panel heading ends 17.9 bp from the top, and its first callout begins 26.5 bp from the top. The following inclusion shows the entire illustration and all callouts at their native size. All callouts print at 10.8 pt or larger. [Rendering at full manuscript width](../../../../../tmp/toolchain_sadb_20261119/figure_readability/bpa_restored_print.png)

```tex
\begin{figure}[htbp]
\centering
\includegraphics[trim=0 0 0 20bp,clip,width=\textwidth]{bpa_exploded_cutaway.pdf}
\caption[Exploded half-section of a proposed BPA wall model]{Illustrative exploded
longitudinal half-section of the proposed BPA wall model. The fixed fitting
constrains the end. The separated surfaces identify the rubber, the two oppositely
wound aramid fiber families, and the outer covering. Sectioning a specimen would
establish the layer dimensions, fiber angles, and interfaces used in a finite-element model.}
\label{fig:bpa_exploded_wall}
\end{figure}
```

## Preservation and verification

Each affected chapter has a `Deprecated/20261119_toolchain` directory containing its former active PDFs, available alt companions, earlier editable source, and a SHA-256 manifest. The original cutaway raster, label-editable PowerPoint, EPS, prompt, and descriptive text are also in the active Chapter 7 Components folder. The portable ZIP preserves 58 entries, including the archived PDFs and sources, restored illustration sources, revised active assets, and the current generator. Every member was decompressed and compared with its on-disk source. [Archival manifest](FIGURE_ARCHIVE.json) · [Portable graphics and sources](WORKFLOW_AND_CUTAWAY_SOURCES.zip)

## Figure 7.5 in “From Design Tools to Neuromechanical Experiments”

Ben requested a Chapter 7 adaptation of the original defense workflow. The added figure joins the controller and its compatible mechanical model at the muscle commands and sensory inputs. It then shows how to run a common baseline and repeat it after changing either one sensory pathway or one actuator attachment or route. Each branch names the settings held fixed. The two traces are aligned at stimulus onset or a comparable foot-contact event to identify where the first difference appears in neural activity, muscle force, joint motion, or foot loading. [Chapter 7 PDF](../../Figures/70-overall/neuromechanical_experiment_workflow.pdf) · [Alt text](../../Figures/70-overall/neuromechanical_experiment_workflow.alt.txt)

The portrait vector is 396 × 405 pt (5.50 × 5.625 in). At the intended 5.5-inch width, body labels print at 11 pt, headings print at 12 pt, and each of six directional tips measures 9 × 8 pt. At a 5-inch width, type remains 10.0–10.9 pt. Visual checks at both widths showed clear branches and at least 8.93 pt of native horizontal clearance between the longest label and its frame. The PDF embeds Arial and contains zero raster images. Its SVG retains 24 separately editable text objects. A second run produced byte-identical SVG and PDF files. The Components folder preserves the original defense SVG, the adapted SVG, the script, and their hashes. [Editable SVG](../../Figures/70-overall/Components/neuromechanical_experiment_workflow/neuromechanical_experiment_workflow.svg) · [Generator](../../Figures/70-overall/Components/neuromechanical_experiment_workflow/make_figure.py) · [Source and print verification](../../Figures/70-overall/Components/neuromechanical_experiment_workflow/SOURCE.json) · [5-inch print check](../../../../../tmp/toolchain_sadb_20261119/figure_readability/7_5_print_reduced_5in.png)

## Placement in the combined dissertation

Figure 7.2 has one List of Figures entry at printed page 168. Its axial profile and loading tests appear together there. The exploded wall illustration continues on page 169, with its own cross-reference and a caption describing the proposed layer arrangement and the measurements needed to establish it. Figure 7.3 is on page 175. Its solid arrows retain the implemented spinal paths and its dashed arrows identify the proposed additions. Section 7.3 and its entire explanatory paragraph fit on page 177. Figure 7.5 and its caption share page 178 at the original 5.5-inch width with 11–12 pt type. Chapter 8 opens on page 179. The LoF and section links agree with these pages in the 323-page eighth proof. [Final page spread](../../../../../tmp/toolchain_sadb_20261119/figure_readability/final_ch7/proof8_section_7_3.png) · [Combined proof](../../../../../tmp/toolchain_sadb_20261119/proof_8/main.pdf)

The controller-retesting PDF also received the figure-font embedding pass. Its arrow geometry, text, and all original page and form streams compare exactly. [Font repair and examination](FIGURE_FONT_EMBEDDING.md)

In the delivered 324-page proof, the four-point side allowance places panels 7.2A/B on p. 167 and 7.2C on p. 168. Figure 7.3 appears on p. 174, with arrowheads measuring 9.44 × 7.55 printed points. Across the six revised workflows the smallest printed label is 10.06 points. The complete comparison paragraph is on p. 176, and Figure 7.5 and its caption occupy p. 177. That figure retains its fixed 5.5-inch width. Chapter 8 starts on p. 178. [Final figure dimensions and folios](FINAL_VERIFICATION.json)
