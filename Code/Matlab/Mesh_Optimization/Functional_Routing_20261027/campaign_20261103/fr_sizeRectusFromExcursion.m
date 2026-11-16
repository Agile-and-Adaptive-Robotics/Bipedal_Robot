function fr_sizeRectusFromExcursion(outputDir)
% Turn human travel into minimum rest length and a parallel count screen.
C = fr_setup();
grid = readtable(fullfile(outputDir,'rectus_pose_grid.csv'));
summary = jsondecode(fileread(fullfile(outputDir,'rectus_excursion.json')));
T=fr_excursionSizing(C,grid.human_path_m,grid.bay_m,summary.human_fmax_N);
T.parallel_at_1x_human_Fmax=T.parallel_count_screen;
T.parallel_at_1_25x_human_Fmax=ceil(1.25*T.human_Fmax_with_scale_N./T.minimum_cable_force_per_BPA_N);
T.parallel_at_1_25x_human_Fmax(~isfinite(T.parallel_at_1_25x_human_Fmax))=NaN;
T(:,{'parallel_count_screen','initial_10mm_half_Fmax_N','human_Fmax_with_scale_N'})=[];
writetable(T,fullfile(outputDir,'rectus_sizing.csv'));
disp(T);
end
