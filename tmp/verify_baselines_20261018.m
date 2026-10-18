% verify_baselines_20261018.m - zero-correction (improved-BPA-char-alone)
% baselines for torque tables 4.2/4.3 via the evaluators' documented rigid
% escape paths, plus numeric front-row dumps. Read-only; no mats written.
root = 'D:\Github\Bipedal_Robot';
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));

fprintf('--- flexor zero-correction (rigid, Xi0=0, X1=X2=Inf) ---\n');
z = minimizeFlxPin(0, Inf, Inf, 1:5);
disp('test  RMSE      FVU        MaxRes');
for i = 1:5, fprintf('%4d  %.6g  %.6g  %.6g\n', i, z(i,:)); end

fprintf('--- extensor zero-correction (rigid, Xi0=Xi3=0, X1=X2=Inf) ---\n');
ze = minimizeExtX3(0, Inf, Inf, 0, 1:9);
disp('test  RMSE      FVU        MaxRes');
for i = 1:9, fprintf('%4d  %.6g  %.6g  %.6g\n', i, ze(i,:)); end

S = load('minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat', 'filtered_results');
fprintf('\nFLX filtered_results row 77  = %s\n', mat2str(S.filtered_results(77,:), 6));
fprintf('FLX filtered_results row 107 = %s\n', mat2str(S.filtered_results(107,:), 6));

E = load('minimizeExt10mmX3_results_20260920_flxr77.mat');
fprintf('\nEXT top-level fields: %s\n', strjoin(fieldnames(E)', ', '));
fprintf('EXT filtered_results row 32  = %s\n', mat2str(E.filtered_results(32,:), 6));
if isfield(E, 'baselineScores')
    disp('EXT baselineScores:'); disp(E.baselineScores);
end
lf = {'LOCKSRC','xi1lock','xi2lock','flxRow'};
for k = 1:numel(lf)
    if isfield(E, lf{k}), fprintf('%s = %s\n', lf{k}, mat2str(E.(lf{k}), 10)); end
end
fprintf('\nBASELINE VERIFY DONE\n');
