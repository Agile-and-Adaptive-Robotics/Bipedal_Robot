% dump_holdouts2_20261018.m - extensor holdout fields + front column structure.
% No evaluator calls, no parpool, read-only.
root = 'D:\Github\Bipedal_Robot';
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));

E = load('minimizeExt10mmX3_results_20260920_flxr77.mat', ...
    'exthold','holdTests','trainTests','holdoutIdx','validLabels', ...
    'train','val','pick','labels','TAG','numBPA','NUMHOLD','holdCols','trainCols','fold');
fprintf('EXT flxr77: pick=%s TAG=%s exthold=%s numBPA=%s NUMHOLD=%s\n', ...
    mat2str(E.pick), mat2str(E.TAG), mat2str(E.exthold), mat2str(E.numBPA), ...
    mat2str(E.NUMHOLD));
fprintf('EXT labels: %s\n', strjoin(string(E.labels), ' | '));
fprintf('EXT holdTests=%s trainTests=%s holdoutIdx=%s fold=%s\n', ...
    mat2str(E.holdTests), mat2str(E.trainTests), mat2str(E.holdoutIdx), mat2str(E.fold));
fprintf('EXT train=%s val=%s holdCols=%s trainCols=%s\n', ...
    mat2str(E.train), mat2str(E.val), mat2str(E.holdCols), mat2str(E.trainCols));
fprintf('EXT validLabels: %s\n', strjoin(string(E.validLabels), ' | '));

F = load('minimizeExt10mmX3_results_20260920_flxr77.mat', 'filtered_results');
cols = unique(F.filtered_results(:, 2:4), 'rows');
fprintf('EXT unique (col2,col3,col4) across all %d front rows:\n', size(F.filtered_results, 1));
disp(cols);

S = load('minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat', 'filtered_results');
colsF = unique(S.filtered_results(:, 2:3), 'rows');
fprintf('FLX unique (col2,col3) across all %d front rows:\n', size(S.filtered_results, 1));
disp(colsF);
fprintf('\nHOLDOUT DUMP 2 DONE\n');
