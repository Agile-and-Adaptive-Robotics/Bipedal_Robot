# Embedded fonts in the dissertation figures

**Current 325-page proof.** Nine earlier embedding repairs are directly included under their original names. Three earlier artwork paths have been superseded. The current audit also checks the four final Chapter 4 figure files, including the renamed extensor validation figure. Its font-repaired bytes match the preceding version. All 68 PDF figure inputs and 182 compiled fonts passed. [Active paths and hashes](ACTIVE_FIGURE_FONTS.json) · [Current audit](ACTIVE_FIGURE_FONT_AUDIT.json)

## Earlier twelve-figure embedding audit

Twelve included PDFs contained 18 unembedded Helvetica resources. Nine used regular Helvetica with WinAnsi encoding, seven used bold Helvetica with WinAnsi encoding, and two used bold Helvetica with Adobe StandardEncoding. The repair embeds the regular and bold Nimbus Sans font programs that MuPDF uses to display Helvetica. Each edited PDF begins with every byte of its original. An incremental update replaces the affected font dictionaries and appends their font programs, descriptors, and character widths. The existing graphics, figures’ page boxes, word positions, and previously embedded font programs compare exactly. [Twelve source and output hashes](FONT_EMBEDDING_FILES.json) · [Object and rendering audit](FIGURE_FONTS.json)

| Figure-folder root | PDFs | Font dictionaries |
| --- | ---: | ---: |
| 10-introduction | 1 | 1 |
| 15-background | 1 | 1 |
| 30-torque | 5 | 7 |
| 50-database | 1 | 2 |
| 60-sns | 2 | 3 |
| 94-AppendixC | 2 | 4 |

Across all twelve pairs, the 319 original image, page-content, form, and other streams have identical compressed bytes. All 36 MuPDF before-and-after render comparisons had zero differing pixels at 72, 144, and 216 dpi. Poppler independently confirmed an embedded program for every font. Its 216-dpi renders at each figure’s manuscript scale were identical outside the previously unembedded text. Inside those text spans, Poppler had previously used its system font substitute. The protected SADb overview preserves the screenshot, the citation-map layout, and the seven A–G callouts. [Audit details, including font names and text positions](FIGURE_FONTS.json)

The 323-page combined proof contains 194 distinct font resources, each embedded. All 70 PDF figure inputs recorded by LaTeX also contain their fonts. The SHA-256 of the combined PDF is `2a2c88916486555ef4d3e8d45bfe220390a60f2a15e71bdcceb60aef169dd8ce`. The twelve repaired PDFs in the Overleaf transfer ZIP each match their active figure’s hash. A separate extraction of the forensic archive reproduced all twelve two-renderer comparisons. [Combined proof](../../../../../tmp/toolchain_sadb_20261119/proof_8/main.pdf) · [Recorded checks](FIGURE_FONTS.json)

## Reproduction

The final production step is independent of the original drawing application. After regenerating a figure, retain its pre-embedding PDF and pass it to the font tool. The twelve original inputs are stored in each chapter’s `Deprecated/20261119_font_embedding` directory. The archive preserves those inputs, the corresponding active outputs, the scripts, and the protected overview’s editable annotation source and cropped screenshot. The controller workflow and knee-response plotting source are also supplied. [Postprocessing source](embed_figure_fonts.py) · [Verification source](audit_figure_fonts.py) · [Portable original and repaired figures](FIGURE_FONT_SOURCES.zip)

From this Notes directory, the following PowerShell example recreates the exact delivered protected figure in a separate temporary file. PyMuPDF 1.28.2 and fontTools 4.25.0 produced the recorded hashes. The verification script also uses Pillow, NumPy, Poppler’s `pdftoppm`, and `pdffonts`.

```powershell
$before = '..\..\Figures\50-database\Deprecated\20261119_font_embedding\sadb_bubble.pdf'
$check = Join-Path ([System.IO.Path]::GetTempPath()) 'sadb_bubble_font_check.pdf'
& 'D:\Anaconda\python.exe' -B .\embed_figure_fonts.py $before $check
Get-FileHash -Algorithm SHA256 -LiteralPath $check
```

The expected SHA-256 is `6ada03600eb535f377a86bc8eada4f91a266b938ccf36a5e7cabdec5cf0d3d25`. The archived object and rendering report preserves the original 12-figure check. The current audit repeats the byte-identical embedding checks for the nine retained paths and verifies the four current Chapter 4 inputs and combined dissertation. [Historical manifest](FONT_EMBEDDING_FILES.json) · [Current audit](ACTIVE_FIGURE_FONT_AUDIT.json)
