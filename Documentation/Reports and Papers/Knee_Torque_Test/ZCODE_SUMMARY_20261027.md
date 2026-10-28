# ZCode campaign summary — 2026-10-27/28

Status of the torque paper, the dissertation torque chapter, and the supporting
optimization studies, as of the close of the 2026-10-27 session. Machine:
easteregg2. Nothing committed; everything below is in the working tree.

## Deliverables of record

| Item | Path | State |
|---|---|---|
| Dissertation build (final, merged numbers) | `Documentation/Reports and Papers/Dissertation/ProofFinal/build_r11_20261027/main.pdf` | 320 pp, 0 errors, 0 undefined refs/citations |
| Dissertation torque chapter | `Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex` | Restructured per directive; placeholders filled with study numbers |
| Dissertation appendix (new C.7-C.9) | `Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/94-AppendixC.tex` | +133 lines: contact/release angles, pulley surfaces (correct captions), cam trims, restored fit table |
| Torque paper source | `Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex` | 38 pp clean; 10 figs + 4 tables = exactly 15 display items; objectives formalized |
| Paper build of record | `Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.pdf` | In-place double compile, 0 errors, 0 undefined |
| Supplement | `Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/frontiers_SupplementaryMaterial.tex` | Cleaned; S1 = RouteVariants, S2 = alignment photo |
| Return package (diss-side patches from the 10-26 audit round) | `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/REVIEW_INBOX/easteregg2/20261026_1700_torque-physical-audit/` | REVIEW.md + 11 patches + figures + 4 audit CSVs |

## Study results of record

**THE GOAL TABLE** (Ben's ruled two columns; CSV:
`tmp/paper_bridge/capability_table_20261022/goal_table_torque_active_20261028.csv`;
the 21-of-27 variant excludes the six updated rows carrying the incoherent
zero-active signature Ben excluded by ruling in the October capability campaign:
ercspn, semimem, med/lat gas, glut_max, add_mag):

| Routing | % of RoM torque meets or exceeds human | % of RoM with active force |
|---|---|---|
| Masters-era (40 mm sites as 4x 20 mm) | 37.3% (38.5% on 21) | 87.8% (87.1% on 21) |
| Updated optimization | 65.2% (70.2% on 21) | 74.3% (95.5% on 21) |

The updated routing trades nothing on torque and, on the trustworthy rows,
also wins active-force coverage; the all-27 active average is dragged by the
six suspect rows, not by real slack.

**Vastus medialis minimum-BPA** (driver: `tmp/paper_bridge/figscripts/study1_vasmed_canon_20261027.m`;
summary CSV: `tmp/paper_bridge/capability_table_20261022/VasMed_minBPA_study1_summary_20261027.csv`):
no configuration (nBPA <= 3 x tackle <= 3, both stiffness modes, both arm methods)
meets all 5 graded angles. Best: 3 BPAs x 2 tackle parts, stiffness-aware,
cable moment arm — 4/5 angles, margins +56.3 / +40.1 / +35.4 / +1.0 N*m;
fails only at +10 deg where every configuration delivers ~0 torque
(routing-family geometry collapse at full extension).

**Masters-era vs updated, capability v3** (script:
`tmp/paper_bridge/figscripts/study2_capability_v3_20261027.m`; CSV:
`tmp/paper_bridge/capability_table_20261022/capability_v3_three_metric_20261027.csv`):
over the 27 actuators (masters 40 mm sites substituted 4x 20 mm on the original side),
updated meets 65.19% vs original 37.28% of the graded domain (area share; updated
wins 23/27). Average-of-axis-means equals area share to 1e-14 on all 23 full grids;
they differ (up to 10.7 points) only at the 4 hip-envelope sites. NOTE: the third
column (mean margin, the average capped sign-aware torque ratio) is an
AGENT-ADDED DIAGNOSTIC from the 10-27 metric discussion, not one of Ben's ruled
metrics; Ben's ruled columns are the torque share and the RoM/active share, and
the ruled four-row structure (masters/updated x stiffness-aware/unaware) was NOT
built for this comparison — the updated side is the as-solved (stiffness-aware)
designs, and the aware/unaware split exists only in the vastus study. The margin
column reads as depth-of-shortfall: updated -3.63 vs original -1.93 (updated wins
coverage 23/27 but its misses are taut-but-short rather than slack-dead).

**4x20 substitution run** (original side, the eight 40 mm sites):
`tmp/paper_bridge/capability_table_20261022/BiPulley_survey3_orig_4x20_20261027.csv`
+ mats in `Code/Matlab/Mesh_Optimization/Results/BiPulley_surveys/`.

## Open flags for Ben

1. Appendix C.8 dropped the stale "40.2 N*m at 804 N" clause (record of record
   = 27.23 N*m peak); replacement wording is Ben's call.
2. MonoPam port residual: 0.99% at 25 collapsed poses (Funit quirk). The
   stiffness-aware study mode evaluates through MonoPam_pulley itself
   (gate rel-err exactly 0); the 1-D port residual is flagged, not adopted.
3. The staged `Notes/advisor_revision_20261022/Revision_C` package still carries
   the OLD pA/pB swap; the live chapter is correct — do not apply that package blindly.
4. `BiPulley_survey2_perMuscle_20261023.csv`'s matFile column predates the
   Results/BiPulley_surveys consolidation; resolve by basename.
5. Figures folder triage (10-26): used files at root, components in
   `Figures/Figure components/`, unused in `Figures/deprecated/` — all moves, nothing deleted.

## Build recipes (verified this session)

- Paper: in-place double compile inside ProofFinal with
  `D:/MiKTeX/miktex/bin/x64/pdflatex.exe --disable-installer -interaction=nonstopmode -halt-on-error`.
  NEVER use -output-directory (stale in-place Main.aux breaks forward refs).
- Dissertation: fresh dated build dir + pdflatex -output-directory, then COPY the
  prebuilt `main.bbl` from the previous build dir into it (bibtex is not in the
  loop; the aux carries no \bibdata), then pdflatex x2.
- MATLAB: `-batch "cd(...); <functionName>"` with file named to match its
  primary function; interleaved nested functions defeat script/function
  classification — use primary + local-functions-after-end layout.
