%dump_ext_pool_20261026.m -- decisive pool-size evidence for the extensor
% identification of record: numBPA / scores_cv / results_cv / baselineScores
% dimensions and the saved harness labels, plus the flexor mat's NUMHOLD.
root = 'D:/GitHub/Bipedal_Robot';
cd(fullfile(root, 'Testing_Data/2022_02_Festo'));
warning('off', 'all');

E = load('minimizeExt10mmX3_results.mat', 'numBPA', 'scores_cv', ...
    'results_cv', 'baselineScores', 'baselineScores1', 'baselineScores2', ...
    'subT', 'titles', 'validLabels', 'labels', 'allBPA', 'fold', ...
    'flxRow', 'tag', 'pick', 'filtered_results', 'xCols', 'trainIdx', ...
    'holdoutIdx', 'trainTests', 'holdTests', 'trainCols', 'valCols');
if isfield(E, 'numBPA'), fprintf('numBPA = %s\n', mat2str(E.numBPA)); end
if isfield(E, 'allBPA'), fprintf('allBPA = %s\n', mat2str(E.allBPA)); end
if isfield(E, 'labels'), disp('labels:'); disp(E.labels); end
if isfield(E, 'validLabels'), disp('validLabels:'); disp(E.validLabels); end
for f = {'scores_cv','results_cv','baselineScores','baselineScores1','baselineScores2'}
    if isfield(E, f{1})
        v = E.(f{1});
        fprintf('%s: size %s\n', f{1}, mat2str(size(v)));
        if isnumeric(v) && max(size(v)) <= 12
            disp(v);
        end
    end
end
if isfield(E, 'subT'), disp('subT:'); disp(E.subT); end
if isfield(E, 'titles'), disp('titles:'); disp(E.titles); end
if isfield(E, 'flxRow'), fprintf('flxRow = %s\n', mat2str(E.flxRow)); end
if isfield(E, 'tag'), fprintf('tag = %s\n', char(string(E.tag))); end
if isfield(E, 'pick'), fprintf('pick = %s\n', mat2str(E.pick)); end
if isfield(E, 'trainCols'), fprintf('trainCols = %s\n', mat2str(E.trainCols)); end
if isfield(E, 'valCols'), fprintf('valCols = %s\n', mat2str(E.valCols)); end
if isfield(E, 'trainIdx'), fprintf('trainIdx = %s\n', mat2str(E.trainIdx)); end
if isfield(E, 'holdoutIdx'), fprintf('holdoutIdx = %s\n', mat2str(E.holdoutIdx)); end

S = load('minimizeFlxPin10_results.mat', 'NUMHOLD', 'TAG', 'PICK', 'numBPA', 'labels');
if isfield(S, 'NUMHOLD'), fprintf('FLX NUMHOLD = %s\n', mat2str(S.NUMHOLD)); end
if isfield(S, 'TAG'), fprintf('FLX TAG = %s\n', char(string(S.TAG))); end
if isfield(S, 'PICK'), fprintf('FLX PICK = %s\n', mat2str(S.PICK)); end
if isfield(S, 'numBPA'), fprintf('FLX numBPA = %s\n', mat2str(S.numBPA)); end
if isfield(S, 'labels'), disp('FLX labels:'); disp(S.labels); end

% evaluate the adopted extensor terms over all nine tests and compare with
% the published tab:ExtPinGof corrected column, to confirm score provenance
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
g = E.filtered_results(32, E.xCols);
fprintf('row32 terms: Xi0=%.6f Xi1=%.1f Xi2=%.1f Xi3=%.6f\n', g);
[af, ~] = minimizeExtX3(g(1), g(2), g(3), g(4), 1:9);
fprintf('corrected scores at row-32 terms, all nine tests (RMSE/FVU/Max):\n');
for i = 1:9
    fprintf('  test %d: %8.4f %8.4f %8.4f\n', i, af(i,:));
end
disp('DONE_dump_pool');
