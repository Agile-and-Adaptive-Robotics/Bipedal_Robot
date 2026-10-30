# Torque figure polish —2026-10-30

**Final follow-up:** The physical-flexor export now has distinct dissertation/paper font profiles, yielding placed legends of 10.142/10.200 pt at the existing inclusion scales. See `PHYSICAL_FONT_PROFILE_FOLLOWUP.md` and `PHYSICAL_FONT_PROFILE_VERIFICATION.json`. These supersede the initial physical-font and identical-PDF statements below; the Figure 4.11 findings remain current.

Only plot generators and the corresponding four active PDF assets were changed. No TeX, experimental data, model/evaluator, coefficient, fit, optimization, or simulation changed. Figure 4.14 is the physical **flexor**, `fig:flexor_physical_final`; the initially named Ext20 asset is Figure 4.13 and was left untouched.

| Requested change | Result |
|---|---|
| Figure 4.11 A/B | One row/two panels retained. A/B sit outside the graph, left of the respective rotated y-axis label, at exactly the title's vertical anchor. Actual native PDF rasters were inspected. |
| Figure 4.11 legend | Southwest in B. Final order: Unoptimized; Optimized; Human model (max); Measured; Measured at reduced pressure. The two measured entries are contiguous and last. |
| Figure 4.14 fonts/legend | All native text ≤12 pt. MATLAB legend FontSize=10; the renderer emits 9.75 pt native legend glyphs. Legend is padded 8 points inside the axis, clearing ticks and the predictions. |
| Figure 4.14 draw/style | Pressure-specific gold predictions draw before measured blue points. Every scatter has MarkerFaceAlpha 0.75 and MarkerEdgeColor none. The current paper's 1.2 pt vector curve weight is retained for consistency with its companion extensor panel. |

## Generator identity and safe regeneration

- `Testing_Data/2022_02_Festo/Plot_KneeFlxBio_all.m` is the proper Figure 4.11 generator. Its default path now loads retained full-precision plotted series and exports the current A/B layout. The former three-panel 2×2 plotting block is replaced. An explicit `BPA_FLX_BIO_REBUILD=1` opt-in retains the older numerical path; it was **not executed**. Both paths use the same shared layout/legend helper, including the current A/B x-limits.
- `Testing_Data/2022_02_Festo/Plot_KneeFlxPhysical_final.m` is the plot-only physical-flexor entry point.
- `Testing_Data/2022_02_Festo/Polish_KneeTorque_Figures.m` supplies shared formatting and synchronized dissertation/paper export. Normal calls update active assets; `BPA_FIGURE_STAGE_ONLY=1` stages review candidates.

MATLAB version 25.1.0.2973910 (R2025a Update 1) was verified. All three .m files passed checkcode with zero messages. Drawnow completes the actual native legend geometry before vector export; the PDF itself, rather than only a MATLAB preview, was rendered and inspected. The PDF-skill edit-operation marker was invoked once before the first export.

## Data conservation

The retained dissertation and paper FIGs have identical scientific X/Y/marker-size datasets. Those FIGs were also checked against the **current active PDFs before editing**: visible curve vertices and marker centers agree to less than 0.00005 pt after the graphic-axis coordinate transform. This includes all 4 reduced-pressure points in Figure 4.11 and all 15 gold/15 blue points in Figure 4.14. The Figure 4.11 B measurement on the axis boundary is retained in the full series; its clipping at that boundary is inherited. All retained XData, YData, and SizeData are asserted bit-identical before/after formatting.

See CURRENT_PDF_VS_SAVED_SERIES.json, SAVED_FIGURE_INSPECTION.json, the per-figure verification JSON, and FINAL_TORQUE_ASSET_MANIFEST.json. Original PDFs/generator are backed up under Dissertation/Archive/figure_polish_20261030/before with their repository-relative paths. The proper generator diff is Plot_KneeFlxBio_all.m.diff. New editable FIG outputs are retained in each manuscript's Figures/Components/figure_polish_20261030 directory.

## Native fonts versus document placement

The physical-flexor PDF is 470×279 pt. Its largest native text is 12 pt and nominal 10 pt legend exports as 9.75 pt because of renderer sizing. With the pre-existing dissertation 0.97textwidth inclusion, the legend is approximately 9.42 pt; the paper's scale 0.85 gives 8.29 pt. These measurements were reported to the lead and build agents for actual placed-page assessment; no TeX scale change is hidden here.

The Figure 4.11 PDF is 699×281 pt; 15 pt native legend text places at approximately 10.04 pt at 468 pt width. Both manuscripts' corresponding assets are synchronized byte-for-byte. Final hashes and paths are in the manifest. The lead/build agents own final compilation, placement QA, captions, and uploads.
