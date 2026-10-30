# Physical flexor: final placed-font correction

The initial shared PDF reduced the legend below 10 pt after manuscript placement. This follow-up supersedes the physical-flexor font settings and identical-PDF statement in the initial release record. Figure 4.11 and the extensor asset were unchanged.

| Setting | Dissertation | Paper |
|---|---:|---:|
| Native legend size | 10.5 pt | 12 pt |
| Final placed legend size | 10.142 pt | 10.200 pt |
| MATLAB tick font | 11 pt | 12 pt |
| Maximum native text size | 12 pt | 12 pt |
| Maximum placed text size | 11.590 pt | 10.200 pt |
| Inclusion | Existing 0.97 textwidth | Existing scale 0.85 |

`Testing_Data/2022_02_Festo/Polish_KneeTorque_Figures.m` now exports distinct physical profiles from the same retained full-precision series. `Plot_KneeFlxPhysical_final.m` remains the normal plot-only entry point. The per-target font compensation is documented directly in the shared generator. Labels are reasserted after changing tick sizes, preventing MATLAB's automatic label scaling from exceeding 12 pt natively.

The legend starts 14 pt inside the axis, providing a visible gap from the major y ticks. Both final native PDFs were rendered and visually inspected: no legend/tick, curve, or marker overlap. Curve linewidth remains 1.2 pt natively, matching the unchanged companion paper panel; pressure predictions render before measured points. All scatter faces remain at alpha 0.75 with no marker edges.

Each profile asserts exactly unchanged XData, YData, and SizeData. No data file, fit, evaluator, optimization, or simulation was run or changed. The shared MATLAB generator passes checkcode with zero messages. Pre-follow-up helper and physical PDFs are preserved under `Dissertation/Archive/figure_polish_20261030/profile_before`.

Final hashes, exact native/placed font measurements, active paths, and per-profile evidence are in `PHYSICAL_FONT_PROFILE_VERIFICATION.json` and the updated `FINAL_TORQUE_ASSET_MANIFEST.json`. Actual native-PDF rasters and editable FIGs are under `physical_profiles/dissertation` and `physical_profiles/paper`; each manuscript also retains its corresponding FIG in `Components/figure_polish_20261030`.

The top-level physical PDF and `release_torque_assets.py` belong to the initial shared export. For the final styles, use the proper MATLAB entry point above, or the two candidates and `release_physical_font_profiles.py` recorded here. The final two PDF styles differ; their scientific arrays are identical.
