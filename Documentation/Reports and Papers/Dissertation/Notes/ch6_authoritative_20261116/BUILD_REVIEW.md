# Canonical build and final page review

The accepted build is **r7**, compiled from the complete canonical `Final/main.tex`. The existing artifact was retained after the connection returned: all 23 files in `Final` still match the build's frozen SHA256 values, and the PDF hash is unchanged. The resumed review required no manuscript edit, figure edit, or rebuild.

## Accepted artifact

- Entry: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Final\main.tex`
- PDF: `D:\GitHub\Bipedal_Robot\tmp\ch6_authoritative_20261116\build_final_r7\main.pdf`
- SHA256: `03314775553c05483786f42950e991512e18e8958169d6cc6c29a920a97ce104`
- Extent: **318 PDF pages**, 18,375,404 bytes; successfully opened and parsed with PyMuPDF.
- Full command, environment, source-hash, timing, and exit-code record: `BUILD_RESULT_r7.json`.
- Canonical result pointer and final review closure: `BUILD_RESULT.json`.

The PDF remains an ordinary local file in the repository's ignored `tmp` area. No PDF was installed in `Final`, substituted for manuscript sources, or copied over the unrelated root dissertation PDF.

## Compilation and source conservation

The build used `D:\MiKTeX\miktex\bin\x64\pdflatex.exe`, BibTeX with absolute `BIBINPUTS` and `BSTINPUTS`, and four complete pdflatex passes. All five commands exited **0**. Their measured command/verification phase was **34.984 s**, excluding preparation of the support cache. The individual command timings were 8.969 s, 0.109 s, 9.187 s, 8.438 s, and 8.078 s. The complete fresh build directory is outside `Final`.

The installed runtime initially encountered registry/format access errors. The successful run used scratch-local MiKTeX startup/configuration/data/install roots and byte-identical copies of readable installed support files and the existing format. Scratch format timestamps were refreshed before each pass to avoid MiKTeX's post-BibTeX format-freshness failure; the format bytes remained identical. Installed files, the original format, registry settings, and manuscript files were not changed by this recovery. `BUILD_RUNTIME_DIAGNOSIS.md`, the driver, per-attempt records, and runtime manifests retain the details. Earlier interrupted or failed attempts are not accepted builds.

Every file under `Final` was hashed before and after r7; all 23 hashes were conserved. On resumption, all 23 current hashes were checked again against that baseline and matched. In particular:

| File | SHA256 |
|---|---|
| `main.tex` | `b8d2ce675f7fce05f6ad392f1a9039c44783e76af612794930e5fc6ad8610b6c` |
| `chapters\60-sns.tex` | `eaf8a81eea28276b980438f1edb93aced8fd7e1fc3612b0ab1c37f7b2a1032ed` |
| `chapters\95-AppendixD.tex` | `e3a19ee49d35764a220d706121348ee518c6fcf8187d57901c97ea502b4e0ea0` |

The final-pass gates found **no build errors, undefined references or citations, duplicate TeX labels, missing inputs or graphics, oversized floats, or request for another reference pass**. `BUILD_STATIC_CHECK.json` separately checks 19 TeX sources, 331 labels, and 158 cited keys. Static checks are supporting checks; the actual compiled PDF and successful bibliography/pass sequence establish compilation.

## Final rendered-page coverage

Full-page renders are 1,224 by 1,584 pixels at 144 dpi, corresponding to the complete publication-page frame. Review included tables, figure text, panel mappings, legends, native signs, captions, page transitions, and margins. Unchanged body pixels were compared exactly, without resampling, against previously inspected full r6 pages; this carryover is identified explicitly below and is not described as a fresh direct inspection of every final page.

| Final PDF pages | Printed pages | Coverage and record |
|---|---|---|
| 153-183 | 133-163 | Complete Chapter 6. Final body regions on 153, 156, 160-162, 168-172, and 174 match the fully inspected same-number r6 pages byte for byte at pixel level (`BUILD_PAGE_COMPARISON_r6_r7.json`). The coordinating reviewer directly inspected final 154-155, 157-159, 163-166, 173, and 175-178. The root reviewer closed final 167. `COMPLETION_PROSE_REVIEW.md` directly closes final 179-183. |
| 273-294 | 253-274 | Complete first Appendix D segment. `AIR_COMPILED_REVIEW.md` identifies exact unchanged-body carryover from r6 and full final inspection of all changed/new pages. It also directly inspects the final transition 289-296. |
| 295-314 | 275-294 | Complete remaining Appendix D segment. Every final full page was directly inspected individually; final closure is persisted in `CAD_COMPILED_REVIEW.md`. |

Thus all **31 Chapter 6 pages** and **42 Appendix D pages** have final visual coverage. The complete document was compiled; this report does not claim a new full-page visual review of unrelated chapters.

Material layout findings from the earlier rendering were corrected and closed in r7: the OurWalker viewport now retains the readable 3D view without its tiny internal joint-label grid; the CAD trace caption fits; earlier-section floats no longer interrupt the main or appendix PF discussion; and the Appendix D source-map handoff path renders correctly. Final Figure D3 pages 291-292 precede the PF subsection on 293, and Figure D4 fits on 294. The main knee, cup, aligned-air, CAD, PF, and spiking figures and their captions remain legible at the final placement. Table 6.4 and the multi-page appendix source map fit.

The last six assigned Chapter 6 pages, 166, 173, and 175-178, were inspected individually after reconnection. No clipping, unreadable labels, caption collision, or material float defect was found. Page 173 retains white space before the flushed W2L figure. A single paragraph-opening line below Figure 6.8 on page 177 continues on 178; it is a minor typographic artifact, with complete readable text. No claim of eliminating every widow or club is made.

## Figure-font check and remaining notices

`BUILD_ACTUAL_FIGURE_FONT_CHECK_r7.json` measures the actual transformed vector-figure text spans in the accepted PDF, cross-checked against source figure records. The minimum detected figure text is **10.000599861 PDF bp**, with **zero spans below the 9.98 bp tolerance threshold**. This is a compiled-placement check, rather than a source-size assumption. Raster text has no recoverable PDF font size and was assessed visually; a numeric raster font floor is not claimed.

No material visual finding remains in the requested chapters. The successful build is not warning-free. The final log retains an unrelated Appendix C overfull box of **0.9562 pt** at its source lines 15-16, underfull-box notices, older included-PDF version/page-group warnings, and a duplicate PDF destination `page.1` from existing front-matter pagination. This PDF destination warning is distinct from duplicate TeX labels, of which none were detected. Those unrelated sources were preserved. The earlier first-pass Appendix C 18.94823 pt box resolved after reference convergence; it is not a final-pass defect.

This build/page review does not recover missing historical launch inputs or compiled binaries, establish physiological validation, or remove the scientific evidence limitations disclosed in the evidence reports and revision ledger. No simulations, training, optimization, model points, transforms, stiffnesses, bounds, flags, numerical samples, or actuator equations were changed during this verification.
