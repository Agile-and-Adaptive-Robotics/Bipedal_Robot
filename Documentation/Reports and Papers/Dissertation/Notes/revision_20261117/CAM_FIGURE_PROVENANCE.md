# Appendix C4 and C5 cam-route callouts

**Proof.** In the 313-page `build_r5/main.pdf`, Figure C4 occupies printed p. 243 (PDF p. 263) and Figure C5 occupies printed p. 244 (PDF p. 264). The PDF SHA-256 is `10126c24720774c651db780ab91a6b88d6a18e4d5262bb3723b9e73b54255409`. At full page size, all 51 active points have distinct labels. C4 panel A shows p3, p4, and p5 individually at the end of the femoral segment and p6, p7, and p8 across the tibial cam. C5 panels G and H distinguish the adjacent p7 and p8 contacts.

## Geometry and material values

The original C4 and C5 PDFs contain the same registered 415 × 663 pt MATLAB plot as an embedded vector form. This is the artwork used below. The registered MATLAB generator reads `ctxCamBest` and `predBest.routeInfo` from `Code/Matlab/Mesh_Optimization/Results/Vas_Pam_20mm_Result_cam_20261014_0416.mat` (SHA-256 `caaf5e4d6a0ce1abc1179a3c3943389c2c00c6468c0a17424e3012f74219623f`). The generator applies `T_ICR_t1`, then `T_Pam`, to carry p6–p9 into the femur frame. The registration file describes placement of the reference bones. Source locations, relative to the dissertation directory, are below.

| Source | SHA-256 |
| --- | --- |
| `Figures/94-AppendixC/Components/registration_20261028/diss_cam_c_registered.m` | `9239e989b2a883fb6f615fd3bf330429958fcbda4c9d9b9cebe1c6999a7b3877` |
| `Figures/94-AppendixC/Components/registration_20261028/appendix_c_bone_background.m` | `e13655e07ac51e3ad5df275e365f0b98d0048e7c10407c3813117edda4ab8172` |
| `Figures/94-AppendixC/Components/registration_20261028/registration.txt` | `ce6020b2d0553f2d278fea8fe3b5ae8bf6ad2b9e453b146fb5e1021dbce47872` |
| `Figures/94-AppendixC/Components/registration_20261028/source_constants.json` | `68b4ea7cbf6821bce55cd1bd6d73b87c455aac87435fabeaa816d72acc95793b` |

The A–F figure occupies its original 490 × 480 pt page. G–H occupies its original 490 × 360 pt page. The revised PDFs carry the complete original plotting stream after removal of the prior point-label text blocks. The resulting stream in each PDF hashes to `d2fe6d031c28ef04d4032070a38a8fd96b47ccd4c53eb4cc9797cfcf6fc1ec20`. This exact binary comparison covers every route, bone outline, marker, clearance envelope, bend-radius line, and coordinate axis. Both original clipping windows and placement matrices also compare exactly. The older full-page vector in `Deprecated/Components/cam_route` has a different bone drawing. The archived active figures supply the complete registered artwork for this edit.

Each new leader points to an original route vertex and stops at the marker perimeter. MATLAB independently exported the model coordinates and active row identities for the following poses. The final column is the maximum Euclidean distance between any transformed model coordinate and the corresponding original, plotted vertex, after fitting each panel's original equal x/y scale and translation. The original and revised vector vertices compare exactly in all panels.

| Panel | Displayed knee angle | Active p rows | Maximum model-to-artwork difference (pt) |
| --- | --- | --- | ---: |
| A | −120° | 1, 2, 3, 4, 5, 6, 7, 8, 9 | 0.00000475 |
| B | −102° | 1, 2, 3, 4, 6, 7, 8, 9 | 0.00000998 |
| C | −83° | 1, 2, 3, 6, 7, 8, 9 | 0.00001301 |
| D | −65° | 1, 2, 3, 6, 7, 8, 9 | 0.00001672 |
| E | −45° | 1, 2, 6, 7, 8, 9 | 0.00001619 |
| F | −27° | 1, 2, 6, 7, 8, 9 | 0.00002260 |
| G | −8° | 1, 7, 8, 9 | 0.00001380 |
| H | +10° | 1, 7, 8, 9 | 0.00001583 |

The canonical MAT gives `routeInfo.wrap.radius(2,1) = 0.0969033864406247 m`. Its `eliminatedAngleD` values round to p5 at −110.8°, p4 at −89.8°, p3 at −54.3°, p2 at −21.5°, and p6 at −18.9°. All five release values and the 96.9 mm radius are visible in the final C5 key.

## Published artwork and originals

The active chapter directory is `Figures/94-AppendixC/`. The two `.alt.txt` companions follow the PDF basenames. `Figures/94-AppendixC/Deprecated/20261117_cam_callouts/` holds both superseded PDF files and `manifest.json` with their original paths, byte lengths, and SHA-256 hashes. `CAM_FIGURE_ORIGINALS.zip` in these dated notes provides a second copy under the original `Figures/94-AppendixC/` paths and contains the same manifest. The archive SHA-256 is `b660710213ad5a86919b74dd05a4b2b87376b8d336ecd5766e275e78f3ad9fdf`. The repository's `.gitignore` line 111 covers every chapter's `Deprecated` directory. Git picks up the Notes ZIP normally. Tracking the three unpacked archive files separately requires an explicit forced add.

| Active file | Archived SHA-256 | Revised SHA-256 |
| --- | --- | --- |
| `extensor_cam_AF_bones.pdf` | `108bdda44425e7343be2e20958d43d4009ca16f7193c5e7fcb01f69efc9b1a0d` | `18aa0b081a4f75270d9bf6f07d9ccbea510fd429b94f588f219b6433bd768f30` |
| `extensor_cam_GH_bones.pdf` | `76da6e7b41bff1c940a715b2e2b4e259b04c512c83749cadef01378bfae4e275` | `20379eb28766cfdc27644f263a3bbd50f5ef36d029eb8e1dd4d73167858ed5fd` |

## Reproduction and final visual check

From the repository root, call MATLAB `addpath('Documentation/Reports and Papers/Dissertation/Notes/revision_20261117'); export_cam_label_geometry` to regenerate `CAM_LABEL_GEOMETRY.json`. Run `make_cam_label_figures.py OUTPUT_DIRECTORY` with the document Python, which provides `pypdf`, `pdfplumber`, and ReportLab. The script reads the archived PDFs and MATLAB reference values, and writes both revised figures and `CAM_LABEL_VERIFICATION.json` to the specified directory. That verification file gives all 51 original on-page coordinate pairs, the new label positions, the placement matrices, the byte hashes, and the model residuals. The accepted copy is stored alongside the generator. A second run reproduced both active PDFs and the verification JSON byte for byte.

The final dissertation embeds each figure at a scale of 0.95514. The 11.3 pt Arial Bold point labels measure 10.79308 pt on the PDF page. Both axis titles and the C5 key use 11 pt artwork text (10.50654 pt placed). Thirteen former 10 pt C5 key and note text runs were enlarged. One original C5 legend swatch encoded its dash pattern as `[4 3] d`, omitting the required phase. On the old artwork Poppler stopped after the fifth key entry. The corrected `[4 3] 0 d` prints the intended swatch and all following entries.

Both final pages were rendered with Poppler at 144 dpi. Visual inspection and extracted PDF text confirmed all 11 C5 keys, including the dashed bend-radius swatch, and both lines of release angles. The original marker centers, axes, and material values can be checked against the records above. The proof images are in `tmp/ch6_reader_20261117/cam_trials/v2/placed-263.png` and `placed-264.png` at the repository root.
