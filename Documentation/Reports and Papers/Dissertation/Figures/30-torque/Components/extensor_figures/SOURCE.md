# Pinned-knee extensor torque and cam-assisted extensor route

Both figures are drawn for a 468 pt (6.5 in) manuscript width. Figure 4.9 is 351 pt tall and Figure 4.12 is 343 pt tall. Both use embedded Arial regular and bold at 10–11 pt. The vector SVGs retain separately editable type. Poppler rendered both PDFs at print scale, including the four torque-series swatches and all eleven route-key entries. [Recorded checks and hashes](verification.json)

## Figure 4.9

The six test configurations use row 32 of `Testing_Data/2022_02_Festo/minimizeExt10mmX3_results.mat`. The published six-panel generator calls `minimizeExtX3` with the four adopted terms and plots the hybrid torque, measured torque, unoptimized model torque, and corrected model torque. The three training tests are 5, 6, and 7. The three validation tests are 1, 2, and 8. The paper and dissertation use the same six tests. The older paper figure `Ext10mm_pin_all.fig` supplies a different layout. The raw observations and the torque curves for the present comparison come from the six-test evaluator. [Original six-test generator](../../../../../../../tmp/paper_bridge/figscripts/gen_ExtPin6Panel_20261026.m) · [MATLAB evaluator](../../../../../../../Testing_Data/2022_02_Festo/minimizeExtX3.m)

The panel headings give the measured resting lengths. Gold triangles and blue circles distinguish hybrid and measured points, including in grayscale. Dashed and solid lines distinguish model predictions. All six panels span −135° to +50° and −2 to 25 N·m. A fine horizontal line passes through the zero torque tick. Panel C shows measured and hybrid values at −125°, −123°, and −120.5°, and measured torque of −0.057265164 N·m at +2°. The earlier figure spanned −120° to +45° and 0 to 25 N·m. Its visible marker centers agreed with the MATLAB records within 0.006 native PDF pt. Its model-curve vertices above 0.2 N·m agreed with the interpolated MATLAB curves within 0.002 N·m. [Earlier six-panel figure](../../Deprecated/pinned_extensor_and_cam/ExtPin_complete_publication.pdf) · [Archived first redraw](../../Deprecated/pinned_extensor_and_cam/initial_axis_limits/pinned_extensor_torque.pdf) · [Current editable graphic](pinned_extensor_torque.svg)

An independent MATLAB replay recovered the same source arrays and predictions in all six tests. In the current PDF, the centers of all 113 measured circles and 113 hybrid triangles agree with the source coordinates within 0.00018 native pt. The 1,083 finite model samples across both curves agree within 0.000023 pt. The closest marker center is 5.65 pt inside its frame, with 3.17 pt of clear space between its outer stroke and the frame. Both zero ticks and the zero torque reference were checked separately on each panel. [MATLAB replay](check_extensor_replay.m) · [Per-panel measurements](verification.json)

## Figure 4.12

The cam data come from `Vas_Pam_20mm_Result_cam_20261014_0416.mat`, the optimization result also used by Appendix C. Panels A and B use the first and last records at the reported angles of −120° and +10°. The plotted route, active-point identities, circle centers and radii, wall and condyle boundaries, bend-radius segments, and bone outlines follow the Appendix C frame transforms and landmark registration. The leader lines separate labels in the full-flexion contact cluster. Panel C is the original photograph of the assembled two-BPA extensor. The full camera frame and orientation are preserved. The eleven-entry key includes the reference bones and the origin and insertion symbols. [Appendix C drawing code](../../../94-AppendixC/Components/registration_20261028/diss_cam_c_registered.m) · [Bone registration](../../../94-AppendixC/Components/registration_20261028/appendix_c_bone_background.m) · [Original hardware photograph](cam_hardware_photo.jpg) · [Editable graphic](cam_extensor_route.svg)

## Reproduction

From this Components directory, run the following commands on easteregg2. The first repeats the adopted extensor calculation in a temporary MAT file and compares the arrays with the stored plotting data. The second draws both figures. The third checks the original inputs and final vector paths, reads the embedded font programs, and renders both figures with Poppler.

```powershell
matlab -batch "check_extensor_replay"
& 'D:\Anaconda\envs\myo\python.exe' -B .\draw_extensor_figures.py
& 'D:\Anaconda\python.exe' -B .\check_figure_sources.py
```

Run `export_extensor_data` in MATLAB to rebuild the plotting data, or pass it an output filename for a separate copy. The plotting script accepts `--figure torque` or `--figure cam` to redraw one figure. Each PDF is stored at the chapter root with its matching alt-text file. Previous plots and the unmodified photograph are in the chapter archive. [MATLAB data preparation](export_extensor_data.m) · [Drawing code](draw_extensor_figures.py) · [Scientific and rendering checks](check_figure_sources.py) · [Previous assets and hashes](../../Deprecated/pinned_extensor_and_cam/manifest.json)
