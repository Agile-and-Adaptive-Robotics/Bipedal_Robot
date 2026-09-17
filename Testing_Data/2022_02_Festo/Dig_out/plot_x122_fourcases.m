% plot_x122_fourcases.m — training/validation figures for all four runs
% (figures are left OPEN on screen; nothing is saved from the GUI)
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
tags = {'pick1_h18','pick107_h18','pick1_all_h3479','pick107_all_h3479'};
for mm = 1:4
    fprintf('--- plotting %s ---\n', tags{mm});
    load(['D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo/minimizeExt10mmX3_results_20260916_' tags{mm} '.mat']);
    run('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo/Dig_out/plotblock_ext_x122.m');
end
fprintf('plots done\n');
