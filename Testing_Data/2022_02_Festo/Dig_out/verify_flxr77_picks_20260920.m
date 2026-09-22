% verify_flxr77_picks_20260920.m
% Pre-launch check for the 2026-09-20 Opt_run / Opt_run_Ext campaign:
% print the flexor row-77 pick and the flxr77 row-32 extensor pick that
% the context builders are about to be pointed at, plus the current
% builder picks for comparison. Read-only.
thisDir = fileparts(mfilename('fullpath'));
festoDir = fileparts(thisDir);
cd(festoDir);

fprintf('--- flexor front: minimizeFlxPin10_results_20260908_2brkt_2trans_noT3 ---\n')
S = load('minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat', ...
    'filtered_results', 'xCols');
fprintf('front size = %d rows x %d cols\n', size(S.filtered_results,1), size(S.filtered_results,2))
fprintf('xCols = [%s]\n', num2str(S.xCols))
g77 = S.filtered_results(77, S.xCols);
fprintf('row 77 (NEW pick): Xi0 = %.6g m, Xi1 = %.6g N/m, Xi2 = %.6g N/m\n', ...
    g77(1), g77(2), g77(3))
g107 = S.filtered_results(107, S.xCols);
fprintf('row 107 (current builder pick): Xi0 = %.6g m, Xi1 = %.6g N/m, Xi2 = %.6g N/m\n', ...
    g107(1), g107(2), g107(3))
g1 = S.filtered_results(1, S.xCols);
fprintf('row 1:  Xi0 = %.6g m, Xi1 = %.6g N/m, Xi2 = %.6g N/m\n', g1(1), g1(2), g1(3))

fprintf('\n--- extensor front: minimizeExt10mmX3_results_20260920_flxr77 ---\n')
E = load('minimizeExt10mmX3_results_20260920_flxr77.mat', ...
    'filtered_results', 'xCols');
fprintf('front size = %d rows x %d cols\n', size(E.filtered_results,1), size(E.filtered_results,2))
fprintf('xCols = [%s]\n', num2str(E.xCols))
g32 = E.filtered_results(32, E.xCols);
fprintf('row 32 (NEW pick): Xi0 = %.6g m, Xi1 = %.6g N/m, Xi2 = %.6g N/m, Xi3 = %.6g\n', ...
    g32(1), g32(2), g32(3), g32(4))

fprintf('\n--- lock consistency: flxr77 lock pair vs flexor row 77 ---\n')
fprintf('flexor row 77 pair  = %.6g / %.6g\n', g77(2), g77(3))
fprintf('extensor row 32 pair = %.6g / %.6g\n', g32(2), g32(3))
if abs(g77(2) - g32(2)) < 1 && abs(g77(3) - g32(3)) < 1
    fprintf('LOCK MATCH: extensor flxr77 front carries the flexor row-77 Xi1/Xi2 pair.\n')
else
    fprintf('LOCK MISMATCH -- check which flexor pair this front was locked to.\n')
end

fprintf('\n--- seed files referenced by Opt_run.m ---\n')
resDir = fullfile(fileparts(festoDir), 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
seedFiles = { ...
    'Bifemsh_20mm_Result_20260910_1234.mat', ...
    'Bifemsh_20mm_Result_20260917_1927.mat', ...
    'Bifemsh_20mm_Result_20260917_2004.mat', ...
    'Bifemsh_20mm_Result_20260918_1052.mat', ...
    'Bifemsh_20mm_Result_20260918_1138.mat'};
for k = 1:numel(seedFiles)
    f = fullfile(resDir, seedFiles{k});
    if isfile(f)
        Sd = load(f, 'xBest');
        fprintf('OK  %s  (xBest %dx%d)\n', seedFiles{k}, size(Sd.xBest,1), size(Sd.xBest,2))
    else
        fprintf('MISSING  %s\n', seedFiles{k})
    end
end
fprintf('\nverify done.\n')
