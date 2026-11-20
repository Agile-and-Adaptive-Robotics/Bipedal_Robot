# Chapter 6 standalone PDF cleanup

Seventeen generated Chapter 6 PDFs were removed across the two lab machines. Twelve were whole-chapter extracts, four were seven-page reading packets, and one was a proof of a proposed chapter opening. Every validated path was rechecked after removal. Full dissertations and scientific figure files were retained.

## easteregg2

The table records repository-relative paths, PDF page counts, printed folios, and the SHA-256 of the deleted bytes. Its eleven whole-chapter copies start with the Chapter 6 opening and run to the end of the chapter. The remaining two files contain selections from that chapter.

| Repository path | PDF pages | Printed folios | SHA-256 before removal |
|---|---:|---|---|
| `Documentation/Reports and Papers/Dissertation/output/pdf/Bolen_Dissertation_Chapter_6.pdf` | 30 | 133–162 | `a08339af8663116678d4dbca69f2a78a291fb3e295b79a47353995b15f312f76` |
| `Documentation/Reports and Papers/Dissertation/output/pdf/Bolen_Chapter_6_Toolchain_20261119.pdf` | 30 | 133–162 | `60e4ba425ac49495063f7e70258c6323d0179dd2a2d38ed9b89c3c51c4be18cb` |
| `Documentation/Reports and Papers/Dissertation/Notes/ch6_sequence_20261118/Chapter_6_reading_pages.pdf` | 7 | 133, 135, 136, 139, 143, 147, 160 | `b9befaf0b080314dfd9ee9e85efffa47b1d617b4399c1b7f5d4251c5bb052b39` |
| `Documentation/Reports and Papers/Dissertation/Notes/ch6_first_use_20261118/Chapter_6_terms_reading_pages.pdf` | 7 | 133, 134, 135, 139, 151, 153, 155 | `c8b55d2b4c8a7b0b207c7fdd274b513a923177165999cb10a2b50bb44fac48d5` |
| `tmp/pdfs/chapter6_standalone_20261118/chapter6.pdf` | 31 | 133–163 | `fd35000ed2d32aa16bb63d41e6250c15d0e67effa1113a2e66f4117b77ab4ff6` |
| `tmp/ch6_advisor_20261119/proof_r9/chapter6.pdf` | 30 | 133–162 | `a08339af8663116678d4dbca69f2a78a291fb3e295b79a47353995b15f312f76` |
| `tmp/ch6_advisor_20261119/proof_r8/chapter6.pdf` | 30 | 133–162 | `14d17bdf989cdfd0ce80abd840117ed22492a4e3f3abeae851b7c4ea4119f41b` |
| `tmp/ch6_advisor_20261119/proof_r7/chapter6.pdf` | 30 | 133–162 | `53af19ba5afb17063b1cf78b4d67f677e3c0dbaceccf6b026ef1e15a6ed71af6` |
| `tmp/ch6_advisor_20261119/proof_r4/chapter6.pdf` | 29 | 133–161 | `6c09767d2383fa6ac32e92b036a0f5936172092d7eb1d6916901dfdcc98a1051` |
| `tmp/ch6_advisor_20261119/proof_r3/chapter6.pdf` | 29 | 133–161 | `eb10aa44b2b32d8a489f24ac514aae8c41eaf43966cd4d87e339859aece3e25c` |
| `tmp/toolchain_sadb_20261119/chapter6_link_trial.pdf` | 31 | 133–163 | `657e272c12a0fe0fa042158dbb3a751d5a0b1fd17c581a777cc3fcfa4c790889` |
| `tmp/toolchain_sadb_20261119/proof_11/chapter6_checked.pdf` | 30 | 133–162 | `60e4ba425ac49495063f7e70258c6323d0179dd2a2d38ed9b89c3c51c4be18cb` |
| `tmp/toolchain_sadb_20261119/proof_10/chapter6_checked.pdf` | 30 | 133–162 | `5ffd86529c73e3a13404bc57f42abf43e29c3cd66ece848c11688cc2f80f3079` |

## EB475WS4

The EB475WS4 repository is rooted at `D:/Github/Bipedal_Robot`. The November 18 public file was identical to the 31-page easteregg2 scratch copy. The two reading packets were identical across machines. The October two-page PDF typeset a proposed Chapter 6 opening and its Table 6.1. Its editable TeX source is still present.

| Repository path | PDF pages | Printed folios or content | SHA-256 before removal |
|---|---:|---|---|
| `Documentation/Reports and Papers/Dissertation/output/pdf/Bolen_Dissertation_Chapter_6.pdf` | 31 | 133–163 | `fd35000ed2d32aa16bb63d41e6250c15d0e67effa1113a2e66f4117b77ab4ff6` |
| `Documentation/Reports and Papers/Dissertation/Notes/ch6_sequence_20261118/Chapter_6_reading_pages.pdf` | 7 | 133, 135, 136, 139, 143, 147, 160 | `b9befaf0b080314dfd9ee9e85efffa47b1d617b4399c1b7f5d4251c5bb052b39` |
| `Documentation/Reports and Papers/Dissertation/Notes/ch6_first_use_20261118/Chapter_6_terms_reading_pages.pdf` | 7 | 133, 134, 135, 139, 151, 153, 155 | `c8b55d2b4c8a7b0b207c7fdd274b513a923177165999cb10a2b50bb44fac48d5` |
| `Documentation/Reports and Papers/Dissertation/Archive/Notes/advisor_revision_20261023/Revision_E/review/targeted_reconciliation/tools_inventory/tex_check/P01_check.pdf` | 2 | Proposed opening and Table 6.1 | `dbbb3ffb26e7bcc9a086547a9efd0ad1f6143d2a21c4f2cde1a5108c91d8085b` |

## Search and verification

- **easteregg2.** Screened 2,787 loose PDFs under the repository, the ignored repository, `D:/temp`, `D:/Users/Ben Bolen`, `D:/Download`, `D:/Documentation`, and Ben’s local Downloads, Desktop, OneDrive, Documents, and temporary directory. Eight additional PDFs in the academic-writing/project/worktree locations were unrelated. The reader checked page counts, opening and ending text, section headings, bookmarks, and printed folios wherever an extract was plausible. The archive check read 1,758 PDF members in 101 ZIPs. It found no Chapter 6-only PDF. Two files with ZIP extensions under `Code/Matlab/SNS_Simscape/dev` are invalid ZIPs. The unrelated temporary installer 7z was excluded.
- **EB475WS4 and redirected drives.** Verified the remote hostname and searched the repository, `D:/Bipedal_Robot_Local_Artifacts`, `D:/temp`, Ben’s research directories and Downloads, the AARL and Codex documents folders, and the available F: and G: media. A fingerprint and member-name pass checked 6,645 PDF entries (1,015 distinct CRC/length combinations) in 392 readable ZIPs. It found no match to a removed file. The two non-ZIP Simscape filenames and an unrelated FPV download account for the three invalid ZIPs. The redirected F: held defense presentations and two unrelated 2018 PDFs. G: contained the SIMS directory. E: was inaccessible and was not a mounted share on easteregg2. This inventory covers the accessible working copies and archives. Git history was preserved.
- **Misleading filenames.** Read the actual contents of the compact `tmp/ch6_reader_20261117/ce_proof/main.pdf` and `table_proof/main.pdf` (Appendices C/E and E) and the 33–34-page placement-test `main.pdf` documents (Chapter 7 with preliminary material). The seven-page `neuromechanical_figure_review.pdf` presents separate scientific figures. All were retained. The November 16 legacy scripts copy a complete dissertation into filenames containing “Chapter6.” Their code compiles and verifies the full project. [November 16 builder](../ch6_revision_20261116/build_final.py) · [November 16 package check](../ch6_revision_20261116/verify_and_package.py)
- **Removal.** Each operation resolved the exact path, confirmed its repository location and SHA-256, and deleted a single named PDF. A second path check confirmed that all 17 files were absent. No archive required rewriting. Both additional academic-writing ZIPs contained no PDFs. The rewritten first-use verifier had the same SHA-256 through EB475WS4’s local path and the redirected `\\tsclient\D` path.

## Generators and report links

- The current verifier publishes the complete dissertation and still compares local and transfer builds, checks fonts, validates figure placement and sources, and records chapter folios from the full PDF. The previous Nov 19 packager retains its complete-dissertation, source-ZIP, chapter-bookmark, and raster checks. Their standalone writers, output paths, and file-hash slots were removed. [Current verifier](../toolchain_sadb_20261119/verify_delivery.py) · [Earlier Nov 19 packager](../ch6_advisor_20261119/package_proof.py)
- Both Nov 18 verifiers now read the reviewed pages in the full dissertation. In separate temporary output directories the historical verifiers produced exactly one PDF each, the 315-page and 317-page complete dissertations. The chapter and appendix checks ran. The four applicable verifier/report templates on EB475WS4 were updated against their exact prior text, as were their published links and file lists. [Sequence proof verification](../ch6_sequence_20261118/verify_proof.py) · [First-use proof verification](../ch6_first_use_20261118/verify_first_use.py)
- Older reading-page links in the review templates, HTML, November 18–19 reports, and file lists were removed or redirected to the appropriate full PDF. The November 18 links open physical page 154 and the older Nov 19 link opens physical page 153. Historical measurements and source packages remain. A read-only smoke check of the current verifier on a complete 324-page build found 194 embedded font resources and zero undefined references, oversized floats, or overfull boxes.
