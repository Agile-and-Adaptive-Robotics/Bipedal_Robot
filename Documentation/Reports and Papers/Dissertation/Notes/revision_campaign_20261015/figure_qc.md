# Figure QC — verdicts only (2026-10-15)

Verdicts-only legibility QC for the torque writer's new figures, gating dissertation
inclusion. Checked per figure: axis labels, font size at print scale, legend contrast,
and sign convention vs the dissertation knee convention (**negative = flexion**).
No figure was edited, regenerated, or restyled by this pass.

## Verdicts

- **Vas_Pam_20mm_cam_winner_1.png — FAIL** — legend contrast: the orange solid
  (`vas_med_r`) and yellow dashed (`Required (+5%)`) curves are close in hue and nearly
  overlap across the mid-range (dash pattern is the only separator), and the yellow
  dashes read weakly on white; fonts (~8 pt effective), labeled axes, and the −113..+8°
  x-range are otherwise acceptable.
- **Vas_Pam_20mm_cam_winner_2.png — FAIL** — legend contrast: the `zero` (gray dashed)
  and `required +5%` (gray dotted) reference lines share the same gray and differ only
  by dash pattern (indistinguishable in grayscale print), and the dotted +5% line
  visually merges with the x-axis baseline; fonts and labeled axes otherwise acceptable.
- **Vas_Pam_20mm_cam_winner_3.png — FAIL** — legend contrast: `0 (no stretch)` and
  `KMAX` are rendered with identical gray dotted styling (distinguishable only by
  legend text, and the "0" line is not clearly visible in the plot); fonts and labeled
  axes otherwise acceptable.
- **Vas_Pam_20mm_cam_winner_4.png — PASS** — single curve with labeled, unit-bearing
  axes (`Knee angle, deg` / `Moment arm, mm`), no legend conflicts, and ~10 pt MATLAB
  tick fonts giving ≈8 pt effective at ~\textwidth print; nothing needs regeneration.
- **route_check_run8_20261014.png — FAIL** — print-scale font size: the 1817×900 px
  (~16.5 in native) 8-panel layout collapses to ~0.38× at a ~6.2 in print width,
  leaving ~10 pt native tick text at ~2.5–4 pt printed — sub-legible at any placement
  short of a full landscape page, and even then marginal.
- **route_check_run8b_20261014.png — FAIL** — print-scale font size plus legend
  contrast: at ~0.37× print scale the ~10 pt native ticks fall to ~2.5–4 pt, the
  p3–p8 route cluster becomes an unreadable blob, several panel pairs rely on
  orange-vs-yellow-only separation (p4/p6), and θk values live only in tiny titles
  with the flexion-negative convention stated nowhere.
- **Knee_Extensor_20mm.gif (first frame) — FAIL** — the extracted first frame is
  155×807 px with the panels' aspect destroyed and no text, axis labels, or legend
  legible at that resolution (<100 DPI even at 1.5 in wide); unusable as a print still.

## Cross-figure note (sign convention)

All four cam-winner panels and both route-check figures plot knee angle over the
−113..+10° range, which is **numerically consistent** with the dissertation convention
(negative = flexion, knee RoM [−120, +10]), but **no figure states the convention** in
any label, legend, or annotation. This is a caption-level fix (state "negative values
denote knee flexion" in each caption), not a regeneration defect — it does not change
the verdicts above, which rest on legend contrast and print-scale font size.

## Method and limits

- Dimensions/DPI read programmatically with PIL (`tmp\_qc_gif_frame.py`): cam winners
  1118–1124×747 px @150 DPI; route figures 1817/1835×900 px @110 DPI; GIF frame 0 of
  199 extracted to `tmp\_qc_gif_first_frame.png`.
- Print-scale assumption grounded in ProofFinal usage: `findstr /s /r /i
  "includegraphics" *.tex` shows figures placed at 0.65–1.0×\textwidth (≈6–6.5 in);
  verdicts computed at 6.2 in print width. Cam winners at 6.2 in = 0.83× scale; route
  figures = 0.37–0.38×.
- Content (labels, legends, titles, sign statements) judged with an AI vision model
  over each image; tick-glyph pixel heights cross-checked with a pixel-projection
  script (`tmp\_qc_text_px.py`, `tmp\_qc_text_px2.py`). Vision estimates for the cam
  set converged on 14–18 px tick digits (≈10 pt MATLAB default); the script's
  bottom-strip medians confirmed ~10–11 px digits for the route figures at 110 DPI.
- **Not checked:** the actual compiled PDF/LaTeX rendering of these files (no
  `\includegraphics` for them exists yet — these are candidates), color appearance
  under a specific print provider's grayscale conversion (judged from RGB only), and
  any figure content beyond the first GIF frame.
