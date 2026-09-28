% reeval_gof_20260926.m - Xi/GoF audit: re-evaluate pick-77 (flexor) and
% flxr77 pick-32 (extensor) per-test GoF [RMSE, FVU, MaxResidual] with the
% CURRENT evaluators of record. Read-only: no mats written, no RNG.
% Run: matlab -batch run("...reeval_gof_20260926.m")
thisDir = fileparts(mfilename('fullpath'));
root = thisDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot'), break, end
    if strcmp(parent, root), error('repo root not found'), end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization')); % must win shadowing
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));

fprintf('MATLAB %s | %s\n', version, string(datetime('now')));

%% Flexor: rows 77 (live pick) and 107 (dissertation row) via minimizeFlxPin
S = load('minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat', ...
    'filtered_results', 'xCols');
g77  = S.filtered_results(77,  S.xCols);
g107 = S.filtered_results(107, S.xCols);
fprintf('\n--- flexor row 77  Xi = [%.6g %.6g %.6g] (minimizeFlxPin, mixed convention) ---\n', g77)
f77 = minimizeFlxPin(g77(1), g77(2), g77(3), 1:5);
disp('test  RMSE      FVU        MaxRes');
for i = 1:5, fprintf('%4d  %.6g  %.6g  %.6g\n', i, f77(i,:)); end
fprintf('\n--- flexor row 107 Xi = [%.6g %.6g %.6g] ---\n', g107)
f107 = minimizeFlxPin(g107(1), g107(2), g107(3), 1:5);
for i = 1:5, fprintf('%4d  %.6g  %.6g  %.6g\n', i, f107(i,:)); end

%% Extensor: flxr77 row 32 via minimizeExtX3 (default 2trans, K=[X1,X2,X2])
E = load('minimizeExt10mmX3_results_20260920_flxr77.mat', ...
    'filtered_results', 'xCols');
g32 = E.filtered_results(32, E.xCols);
fprintf('\n--- extensor flxr77 row 32 Xi = [%.6g %.6g %.6g %.6g] (minimizeExtX3, 2trans default) ---\n', g32)
f32 = minimizeExtX3(g32(1), g32(2), g32(3), g32(4), 1:9);
disp('test  RMSE      FVU        MaxRes');
for i = 1:9, fprintf('%4d  %.6g  %.6g  %.6g\n', i, f32(i,:)); end

fprintf('\nREEVAL DONE (no mats written)\n');
