# Dissertation version comparison

**EB475WS4 · November 22, 2026, Pacific time**

**The live export and finalized working PDF give the reader the same text, figures, tables, and pagination.** The November 20 reconstruction differs on the abstract, acknowledgments, one contents entry, and the first part of Chapter 6. All three have 353 PDF pages. The reconstructed file establishes an earlier editorial state. The identity of the PDF uploaded to ProQuest remains unverified.

## Identification

Paths below are relative to `C:\Users\Ben Bolen\Dissertation_Handoff`. Sizes and SHA-256 hashes were computed directly from the files on EB475WS4.

| Origin | File | Bytes | SHA-256 |
|---|---|---:|---|
| Live Overleaf project, PDF saved by Ben | `Live_Overleaf\Bolen_Dissertation.pdf` | 20,175,010 | `aa51a1d7c060b8f4286797fde36c68d2c0e04fe8703c6a047b876d8516994ce6` |
| Matching live project source, 170 ZIP members | `Live_Overleaf\Bolen_Dissertation.zip` | 35,509,883 | `b7fdeb12d6fa07dfc6c4ee5b84cd98d5ced090874d5d85ae871788a1d693337e` |
| Finalized working reference in Ben’s handoff | `Bolen_Dissertation_Review.pdf` | 21,035,014 | `deaba99867b93382f8f9f4c4c11ecf5b584cbd0b8c7b5f4c3083a6197f9735e2` |
| Provisional rebuild from the recovered November 20 pre-submission sources | `ETD_Candidate\ETD_Candidate_Rebuilt.pdf` | 21,034,004 | `2a16b8ca9b6edef387c7c427f46f5913dfa2b02c9e70f094ba4d6dfb78dd6d7a` |

All four hashes agree with their respective handoff records. The repository-root deliveries **within the dissertation folder**, `Documentation/Reports and Papers/Dissertation/Bolen_Dissertation.pdf` and its ZIP, also match the live pair byte for byte in this checkout.

**Working source provenance.** Git identifies `477ae53231288271cd1fe1516c8decc3cf98bb91` as the November 22 commit “finalized dissertation and bibliography.” Its parent is `0212c64cc36a229402491bce822740db3fda9f85`, the base recorded in the overlay manifest. All 33 individual overlay files match their manifest hashes. For alignment with the Git sources, the export’s `chapters/` and support files map to `Final/`, and its `figs/` maps to sibling `Figures/`. Of the 170 ZIP files, 160 match the finalized commit exactly, eight differ only in line endings, and the finalized abstract differs only in whitespace after its last sentence. The remaining file is an extra nested `chapters/dissertation_upload.zip`, outside the active TeX inputs. Every figure file used from the live export therefore has the same source content as its working counterpart. Ben’s finalized abstract and acknowledgments are present in both.

**Submission provenance.** The candidate’s provenance describes reversing a verified three-file source diff against an archive. It records source-archive SHA-256 `48dbf8a254e08091bc09050b1c2694b6eaab4f1464046c4beaaae6e780a03138`, diff SHA-256 `82b51c59f81a645d11594f8893f6c2b704b2bbc2c2c073f221553488a7b4fe46`, and an earlier Overleaf PDF hash of `8be76780d626f0a1b91d280e8f425eda458aacc9338969c079159d7884bae7bb`. The last value is documentary evidence supplied in both the manifest and candidate provenance. The handoff contains the rebuild. I found no independently identified ProQuest copy or separately identified `_ETD_Submission` export in the handoff, local dissertation tree, or searched user download/document locations. The older release ZIP predates the author edits. Its existence does not date or invalidate the live sources.

## Page evidence

I compared extracted text on every corresponding page, including the text inside vector figures, and rendered all three PDFs in RGB at 90 dpi with one MuPDF version. I compared the rendered pixels and inspected sample changed regions on each of the 34 live/working exception pages. I also compared the multiset of decoded embedded bitmap hashes, widths, and heights per page, and the decoded PDF form streams for the two local builds. These complementary checks distinguish a changed figure, a figure moving on a page, and a renderer drawing the same glyph slightly differently.

| Pair | Same extracted text | Identical RGB page raster | Interpretation |
|---|---:|---:|---|
| Live / working | 353/353 | 319/353 | The 34 exception pages differ by 6–461 pixels each in a 765 × 990 image, 2,832 pixels across the document. Changed marks are edges of the same printed letters, punctuation, numerals, or rules. Embedded bitmap contents agree on every page. The live build uses PDF 1.5 and pdfTeX, the working build PDF 1.7 and MiKTeX pdfTeX. The producers, timestamps, and some internal object encodings differ. |
| Working / ETD rebuild | 336/353 | 336/353 | Exactly 17 pages have changed text and pixels. All per-page embedded bitmap identities and decoded form-stream hashes agree. The principal figure difference is placement relative to the revised prose. |
| Live / ETD rebuild | 336/353 | 302/353 | The same 17 edited pages, plus the 34 small compiler-rendering differences seen between live and working. The two sets are disjoint. |

Printed labels follow the paper, where i is physical PDF page 3 and Arabic 1 is physical PDF page 21. **The complete set of edited printed pages is i, iii, vii, 136–145, 147, 149, 152, and 153.** No change to a plotted experimental trace, panel legend, result number, bibliography entry, or appendix page emerged from the version comparison. This measures agreement across editions. The technical review evaluates claims shared by them.

## Editorial differences

The source lines below refer to the matching live ZIP, relative to `chapters/`. The ETD handoff provides the rebuilt PDF and provenance, including hashes of its three reconstructed source members. It does not supply those three source files individually. Candidate wording and page positions below were checked in that PDF. Historical Git copies that resemble a passage do not establish its ETD source line.

| Printed page(s) | ETD reconstruction → live and finalized working reference | Source location in live ZIP |
|---|---|---|
| i, abstract | “Graphical neural blocks and simulation interfaces support” becomes “I developed graphical neural blocks and simulation interfaces to support.” The final sentence changes from “potential prosthetic and orthotic applications” to the explicit use of simulation to refine mechanical designs and neural controllers before physical robot tests. | `02-abstract.tex` 7 |
| iii, acknowledgments | Funding stays intact. Ben addresses the committee directly, specifies collaborators’ roles, recalls Alex Steele sharing prototypes, places Stu McNeal’s conceptual contribution next to the Chapter 6 work, and links Jasmine Bradley’s figure work to Chapters 3 and 4. The live version adds “It is hard to imagine reaching this point without you two.” Both versions already credit Stu. | `06-ack.tex` 7, 9, 11, 13 |
| 136–137, Chapter 6 opening and §6.1 | The earlier tools inventory gives way to the interaction between robot geometry, controller performance, and feedback. Both give slow simulation and manual tuning as obstacles. The fall when the walker is released and its conference attribution now immediately follow the tuning description. The circuit and joint details follow those observations. | `60-sns.tex` 4, 8, 20, 24, 26 |
| 137–138, §6.2 | The live text describes the closed loop explicitly, from controller commands to muscle force and joint motion and back through sensory signals. The harness explanation identifies pelvis and foot forces. W2L’s contact-feedback question moves to the start of its paragraph. The OpenSim human reference, rerouted counterpart, CAD geometry, and adapted MathWorks interface are introduced by their respective contributions. | `60-sns.tex` 30, 32, 34, 36, 38 |
| 139, Table 6.1 | The CAD row changed from “A torque controller distributes the load across 54 actuator routes” to “An engineering controller converts torque requests to forces for 54 ideal actuators.” This makes the ideal-force assumption explicit in the table, where the adjacent robot-counterpart row permits pressure-driven BPAs. | `60-sns.tex` 58–59 |
| 140–142, §6.2.1 | The live account specifies corrected hip axes, the patellar route, contact conditions, matching joint states, and holding the converted paths fixed while comparing muscle equations. The active-force errors of 27.8/39.4 N and 3.6/1.9 N, and the compliant-tendon result within 0.19 N RMS at 2 ms, are the same. The reference to the clearance experiment moves to the clearance sentence. A new transition distinguishes matching stance time from matching joint motion. | `60-sns.tex` 68, 70, 79, 81, 83, 85, 87 |
| 142–144, §6.2.2 | The robot counterpart description relates muscle length to the rerouted paths at the same joint angles. Both distinguish suspended stepping with biological model actuators from the pressure-driven configuration. The pressure paragraph now points to Tables C13 **and** C14 for site capabilities and Appendix D.1.5 for the bilateral transfer. The earlier paragraph only identified Table C13 for sleeve sizes and counts. | `60-sns.tex` 91, 93, 102 |
| vii and 144–153, later Chapter 6 | Prose flows around the same tools table and the same figures. The §6.4.3 heading moves from 149 to 152; the printed contents entry follows it. §§6.4.1–6.4.2 text moves across the fixed table/figure pages, and the §6.4.3 discussion and Figure 6.8 shift within pp. 152–153. The §6.5 printed heading still begins on 154. | `60-sns.tex` 104–118, 139, 154–156, 174–184 and the Figure 6.8 block |

**What moved in the artwork.** Figure 6.1 on 137, Figure 6.4 on 143, and Figure 6.8 on 153 retain their figure numbers and captions. Their caption top positions change by approximately 76 points upward, 24 points downward, and 32 points downward, respectively, from ETD to working/live. Figures 6.2 and 6.3 remain at the same caption height on 140 and 141. The knee and cup reflex figure pages 150 and 151 have exactly the same local-build pixels. There is no figure replacement associated with moving §6.4.3 past these figure pages.

## Electronic navigation and contents layout

One reproducible defect crosses versions. The live and working PDF’s section links can land before the advertised heading. Among 167 printed contents entries on iv–x, 34 links open one to three printed pages early. The two longest entries required visual matching. The printed page for §6.4.3 was updated in the live file, but its link still opens 149. In the ETD reconstruction that particular link and heading both lead to 149. Other early links occur in all three editions.

| Linked TOC row | Where the live link lands | Where the heading prints |
|---|---:|---:|
| §3.2, contents v | 46 | 48 |
| §6.2.1, contents vii | 138 | 140 |
| §6.4.3, contents vii | 149 | 152 |
| §6.5, contents vii | 153 | 154 |
| §B.5, contents ix | 231 | 234 |

The PDF outline also sends §6.4.3 to 149. The first repair to test is to tie the hyperlink anchor to the typeset heading across float/page breaks. The active source uses `main.tex` 6, the custom section handling in `psutex.cls` 401–421, and the representative heading at `60-sns.tex` 174. Rebuild a temporary copy and repeat the destination-versus-printed-page comparison before deciding whether to change the production class.

Separately, the printed contents on x visibly runs the titles of D.5.1 and D.5.3 into the page numbers (“Harness317” and “Parameters319”). Short optional contents titles at `95-AppendixD.tex` 438 and 473 would restore the gap and preserve the full body headings. The page numbers and destination pages for these two rows are correct.

These observations do not establish a need to replace a submitted file. A file retrieved from the submission system can be identified against the hashes above, then compared using the same page checks if its hash is different.
