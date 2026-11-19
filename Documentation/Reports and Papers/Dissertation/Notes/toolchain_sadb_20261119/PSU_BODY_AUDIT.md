# PSU body-format check before repair, November 19, 2026

This historical audit covers `proof_5`. The release proof has a separate all-page margin and hyperlink check, and the figure repair has its own font-embedding record. [Final print margins and links](FINAL_PRINT_MARGINS.json) · [Figure font embedding](FIGURE_FONT_EMBEDDING.md)

The 324-page proof had portrait US Letter pages, 12-point body text, and double-spaced prose. Nine figures displayed text in an unembedded font, three others declared a font they never displayed, and some viewers drew citation boxes into the margins. All twelve figure resources already existed in the incoming PDF and the earlier local proof.

**Standard and scope.** Portland State requires 8.5 × 11 inch portrait pages, embedded fonts, and one-inch margins on all four sides, including folios. Ordinary text and headings are 12 pt. Tables, captions, and notes may use 10 pt. Figures and appendices must be legible. Prose is double spaced, with exceptions for such elements as captions, tables, references, and appendices. Geometry and font resources were checked on every page. Type and leading were measured in representative passages. The separate front-matter review covers the preliminary-page text and order. [PSU Graduate School formatting requirements](https://www.pdx.edu/gradschool/etd-formatting-requirements)

## Margins and type

- All 324 MediaBoxes and CropBoxes measure 612 × 792 PDF pt with zero page rotation. Poppler at 144 dpi with link annotations hidden measured minimum visible margins of 71.0 pt left, 70.0 pt right, 99.5 pt top, and 112.5 pt bottom. MuPDF independently measured 71.5, 70.0, 99.5, and 112.5 pt. The TeX text block is 6.5 inches wide. Hanging punctuation, hyphens, and occasional letter edges extend up to about 2 pt past that block. With TeX text suppressed in memory, artwork and rules stayed at or inside both 72 pt side margins. Strict one-inch ink clearance would also require adjusting the optical protrusion or adding a small allowance to the text margins. [Page geometry and type settings](../../Final/psutex.cls#L55)
- With annotations shown, Poppler draws a green sliver **at the left paper edge** and a box around the full running-header line on printed pp. **13–14** (PDF pages 34–35). A Suzumori–Faudzi citation crossing a page break generated two link rectangles on each page. One rectangle starts at x = −0.996 pt. In total, 1,660 link annotations have a one-point border. The right margin with borders visible reaches 68.5 pt on printed pp. 11, 19, 35, 72, 97, 133, 156, 172, and 270. The edge marks and header boxes also occur in `proof_r9`. The incoming PDF contains bordered links, but has no paper-edge or full-header link rectangles. The release repair preserves the visible one-point borders and excludes the page edge and running header from split links. [Release preamble](../../Final/main.tex#L135)
- The embedded URW Palladio body font uses 11.9552 PDF pt, equivalent to 12 TeX pt. Successive prose baselines were 23.91 PDF pt on printed pp. 1, 42, 68, 121, 133, 164, and 179. Chapter and section headings use 12 TeX pt. The Fig. 6.4 caption on p. 145 uses 12 TeX pt and 14.45 PDF pt leading. Table 6.1 on p. 136 uses 10.95 TeX pt body type. Tables 3.1–3.2 and 4.5–4.7 use 10 TeX pt for ordinary entries, with smaller mathematical indices. [Class options](../../Final/main.tex#L1) · [Heading definitions](../../Final/psutex.cls#L125)
- Content extraction includes material clipped from the source artwork. For example, Fig. 2.9 on p. 29 retains a source article paragraph, Fig. 2.11 on p. 32 retains a journal heading and adjoining source text, and Figs. 4.12, C4, and C5 retain out-of-crop axes. These produced no visible ink in the reported margins. Page-render measurements include the ink the reader actually sees.

## Font resources in placed figures

Pages below use printed folios. Each row names the active PDF within its chapter figure folder. On all three PDF versions, Helvetica and Helvetica-Bold were the only unembedded families. Each resource was traced through its placed form. To determine visibility, its text was temporarily made invisible in memory while every other drawing instruction was preserved, and the page was rendered before and after. Twelve of the 18 font objects affected visible pixels.

| Figure | Proof 5 page | Earlier proof / incoming page | Source figure | Unembedded use |
|---|---:|---:|---|---|
| 1.1 | 10 | 10 / 10 | `roadmap.pdf` | Helvetica declared, no visible text |
| 2.4 | 22 | 22 / 22 | `steeleknee.pdf` | Helvetica-Bold panel letters |
| 4.7 | 99 | 99 / 99 | `FlxPin_preopt_publication.pdf` | Helvetica-Bold torque-axis label |
| 4.8 | 101 | 101 / 101 | `FlxPin_optimized_publication.pdf` | Helvetica-Bold torque-axis label and a small axis annotation |
| 4.9 | 104 | 104 / 104 | `ExtPin_complete_publication.pdf` | Two Helvetica-Bold resources, neither selected |
| 4.12 | 113 | 113 / 113 | `cam_current_legend.pdf` | Helvetica legend labels |
| 4.13 | 115 | 115 / 115 | `Ext20mm_physical_20261026_publication.pdf` | Helvetica-Bold torque-axis label |
| 5.1 | 122 | 123 / 123 | `sadb_bubble.pdf` | Helvetica and Helvetica-Bold callout text and letters |
| 6.4 | 145 | 145 / 145 | `sns_knee_reflex_authoritative.pdf` | Helvetica declared, no visible text |
| 6.6 | 149 | 149 / 149 | `controller_retest_and_optimization.pdf` | Helvetica and Helvetica-Bold panel headings and circuit labels |
| C4 | 249 | 246 / 247 | `extensor_cam_AF_bones.pdf` | Helvetica coordinate-axis labels, plus an unused Helvetica in the callout overlay |
| C5 | 250 | 247 / 248 | `extensor_cam_GH_bones.pdf` | Helvetica coordinate axes, key, and release-angle text, plus an unused Helvetica in the callout overlay |

The affected sources span Chapters 1, 2, 4, 5, and 6, and Appendix C. The release repair embeds all eighteen Helvetica resources. The dissertation text fonts were already embedded. [Figure repair and verification](FIGURE_FONT_EMBEDDING.md)

## Historical figure diagnostics

Poppler 26.07.0 omits part of two Chapter 4 legends in both `proof_5` and the final `proof_11`. Each original drawing gives the PDF line-dash operator an array but omits its required phase. A test supplied a zero phase in memory and rendered the full pages at 144 dpi. All changed pixels were confined to the following legend entries. MuPDF 1.28.2 displayed the complete legends before and after that test.

- **Figure 4.9, printed p. 104.** The `Unoptimized` swatch uses `[6 6] d`. Poppler shows `Measured` and omits the `Unoptimized` and `Optimized` swatches and labels. The affected form in the final PDF is object 3101. [Pinned-knee extensor source figure](../../Figures/30-torque/ExtPin_complete_publication.pdf)
- **Figure 4.12, printed p. 113.** The `Bend radius` swatch uses `[4 3] d`. Poppler shows the key through `Tibia upper` and omits the six remaining swatches and labels for `Bend radius`, `Optimized route`, `Cam contacts`, `Origin p1`, `Insertion p9`, and `Reference bones`. The affected nested form in the final PDF is object 3311. [Cam-geometry legend source figure](../../Figures/30-torque/cam_current_legend.pdf)

**Compared files.** SHA-256 identifiers below distinguish the precise builds even if the delivery PDF is replaced.

| Copy | Pages | SHA-256 |
|---|---:|---|
| `tmp/toolchain_sadb_20261119/proof_5/main.pdf` | 324 | `5D271D4B860EA95F9774D27BFE1417263CE1E40530EFE0519A996782DEFFFAD3` |
| `tmp/ch6_advisor_20261119/proof_r9/main.pdf` | 319 | `CD550ABCE2FA0984E2415D8E58C59D039A673D0E92F46A6BAF07971FB43A13E3` |
| Incoming `Dissertation/Bolen_Dissertation.pdf`, 11:55 a.m. | 320 | `F8255114C6F84428D0722106470A6F9A28848A2E9B52EBE8A91B538BDF842B17` |
