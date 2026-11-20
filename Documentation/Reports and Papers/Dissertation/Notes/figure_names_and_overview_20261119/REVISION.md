# Live dissertation revision

The existing Bolen-Dissertation project contains the revised chapters, front matter, bibliography, and figures. Its final full PDF has 325 pages. The native compilation reported zero errors, zero warnings, and sixteen informational underfull horizontal boxes. All six of Dr. Hunt’s November 19 Chapter 6 comments remained in the live review overview after the changes. The editor was left on Chapter 6 in Reviewing mode. [Overleaf project](https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb) · [Historical native PDF comparison](OVERLEAF_PROOF.json)

## Chapter 6

The overview on printed p. 134 states the need for systematic controller tuning, explains why the lab robot required its own mechanical model, and ends with what researchers can change and measure in a repeated experiment. Section 6.1 on the same page opens with Ben Bolen and Connor Morrow’s course project, identifies Deng’s circuit and Li’s biped, introduces AnimatLab, and reports what happened when support was removed. The earlier review of the opening was reopened after the advisor identified the repeated introduction. The review table distinguishes his original deductions from the assessment of the revised passages. [Overview and §6.1](../../Final/chapters/60-sns.tex#L4) · [Live page](chapter_opening.jpg) · [Historical review matrix](../toolchain_sadb_20261119/REVIEW_TABLE.html)

## Chapter 3 and the figures

| Printed page | Revision and verification |
| --- | --- |
| 42–43 | Chapter 3 opens with the publication citation and the authors’ published contribution statement. The overview and scientific text begin on the next page. The 10 mm reference points to Figure 3.2(a). [Publication and panel checks](../toolchain_sadb_20261119/CHAPTER3_PUBLICATION_VERIFICATION.md) |
| 105 | The six panels of Figure 4.9 share angle and torque limits and show every measured and hybrid point. Redrawing the figure from the saved MATLAB tests restored the four legend swatches. All 226 markers and 1,083 finite curve vertices were verified in the vector PDF. [Sources and scientific checks](../../Figures/30-torque/Components/extensor_figures/SOURCE.md) |
| 114 | In Figure 4.12, panels A and B show the route at −120° and +10°. Panel C identifies the photograph of the assembled two-BPA extensor in the isometric knee test stand. Both redrawn figures use embedded 10–11 pt Arial. [Live figure](cam_route_overleaf.jpg) |
| 117–118 | Figures 4.13 and 4.14 share plot dimensions, angle limits, typography, and line weights. The flexor caption retains the result across the 15 measured points and the installed lengths. The graph supplies the key. [Measurements and plot checks](VALIDATION_FIGURES.md) |

Panel callouts in figure captions are bold. The caption review retained information needed to interpret trial conditions, sample sizes, and uncertainties. Graphical source files are in Components and superseded artwork is in Deprecated. The combined live figure tree contains 186 files with undated basenames. The cleanup removed 17 standalone Chapter 6 PDFs across easteregg2 and EB475WS4 and retired the extraction steps. [Caption review](../toolchain_sadb_20261119/manuscript_captions.py) · [Artwork inventory](ARTWORK_CHECK.md) · [Archived filenames](ARCHIVED_ART_NAMES.json) · [Chapter PDF cleanup](CHAPTER_PDF_CLEANUP.md)

## Native source and print verification

The project’s own source download contains 212 files. All files match the independently rebuilt transfer tree after normalizing text line endings. The 166 byte-identical files include all binary artwork. The other 46 files differ only in text line endings. The downloaded source contains the 15 revised original documents and the final Figure 4.9 PDF and alt text. [File comparison](OVERLEAF_EXPORT.json)

The initial live build used the project’s TeX Live 2022 setting. Four duplicate citation rectangles reached the paper edge on printed pp. 23, 29, 30, and 36. The project now uses TeX Live 2026 with pdfLaTeX. The full native proof passed Poppler’s visible-ink check at 144 and 300 dpi. All 760 citation borders are visible, all internal link destinations resolve, and all 182 font resources are embedded. At 300 dpi, the smallest left, right, top, and bottom clearances are 73.92, 72.48, 74.16, and 76.80 PDF points. An inch is 72 PDF points. [Print measurements and citation borders](OVERLEAF_PRINT.json)

The updated live and local full PDFs have identical extracted text and folios on all 325 pages. At 72 dpi, the Chapter 3 opening, the four revised Chapter 4 figure pages, and the Chapter 6 opening are pixel-identical. Across the full book, 28 pages have small rendering differences, at most 124 of 484,704 pixels on any page. [Native proof comparison](OVERLEAF_PROOF.json) · [Independent local and transfer rebuilds](../toolchain_sadb_20261119/FINAL_VERIFICATION.json)
