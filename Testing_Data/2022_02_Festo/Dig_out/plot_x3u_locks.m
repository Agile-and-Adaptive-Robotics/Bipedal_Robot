% plot_x3u_locks.m — run the driver's plotting section on both finished x3u mats
% (figures are left OPEN on screen; nothing is saved from the GUI)
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
mats = {'minimizeFlxPin10mmX3_2brkt_results_20260915_L107_DIAG_x3u.mat', ...
        'minimizeFlxPin10mmX3_2brkt_results_20260915_L1trans_DIAG_x3u.mat'};
for mm = 1:2
    fprintf('--- plotting %s ---\n', mats{mm});
    load(['D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo/' mats{mm}]);
    run('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo/Dig_out/plotblock_x3u.m');
end
fprintf('plots done\n');
