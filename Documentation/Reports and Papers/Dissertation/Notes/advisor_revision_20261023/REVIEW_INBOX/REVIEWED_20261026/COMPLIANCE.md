# Directive compliance and verification

These line numbers refer to the delivered patch files, not stale Overleaf positions. Live advisor-thread mapping still requires checking the application anchors.

- Researchers and Gait2392 second: `patches/P01_replacement.tex`, line(s) 1.
- Actual referenced tools table: `patches/P01_replacement.tex`, line(s) 11.
- Table Gait2392 second: `patches/P01_companion_table.tex`, line(s) 10.
- Correct Ib versus commissural claim: `patches/EG_controlled_tests.tex`, line(s) 2.
- Mass count/height/tilt limits: `patches/EG_controlled_tests.tex`, line(s) 4.
- Partial feedback and scored-window scope: `patches/EB_feedback_ablation.tex`, line(s) 1.
- Original engine and original model scope: `patches/DENG_original_engine.tex`, line(s) 3.
- Actual panel letters and preserved labels: `patches/EG_extensor_merged_figure.tex`, line(s) 3, 5.
- Captionless figure counter: `patches/EG_physical_figures_captionless_OPTION.tex`, line(s) 3, 8.
- Two-part tackle: `patches/EG_pulley_discussion.tex`, line(s) 5.

Proof build: pdflatex, BibTeX, and two further pdflatex passes all exited zero. Final log: no errors, undefined citations/references, or overfull boxes. Six references to dissertation material outside this fragment packet are explicitly mocked as X, as recorded in `proof_result.json`. This is a fragment proof, not the assembled dissertation.

Proof length: 14 pages. The new Deng plot was visually inspected before and after layout repair. Figure/model simulations were not rerun.

Overleaf, active chapter source, and existing bibliography were not edited. No commit, push, or outgoing ZCode message was made.
