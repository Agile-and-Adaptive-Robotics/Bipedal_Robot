% dump_holdouts3_20261018.m - guarded extensor holdout dump. Read-only.
root = 'D:\Github\Bipedal_Robot';
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));

E = load('minimizeExt10mmX3_results_20260920_flxr77.mat');
gf = @(nm) local_get(E, nm);
for nm = {'exthold','holdTests','trainTests','holdoutIdx','validLabels', ...
           'train','val','pick','labels','numBPA','NUMHOLD','holdCols', ...
           'trainCols','fold','validLabels'}
  if isfield(E, nm{1})
    v = E.(nm{1});
    if isnumeric(v) || islogical(v) || isstring(v) || ischar(v)
      fprintf('EXT %s = %s\n', nm{1}, mat2str(v));
    else
      fprintf('EXT %s = %s\n', nm{1}, strjoin(string(v), ' | '));
    end
  else
    fprintf('EXT %s : (absent)\n', nm{1});
  end
end

F = load('minimizeExt10mmX3_results_20260920_flxr77.mat', 'filtered_results');
cols = unique(F.filtered_results(:, 2:4), 'rows');
fprintf('EXT unique (col2,col3,col4) across all %d front rows:\n', size(F.filtered_results, 1));
disp(cols);

S = load('minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat', 'filtered_results');
colsF = unique(S.filtered_results(:, 2:3), 'rows');
fprintf('FLX unique (col2,col3) across all %d front rows:\n', size(S.filtered_results, 1));
disp(colsF);
fprintf('HOLDOUT DUMP 3 DONE\n');

function s = local_get(~, ~)
s = [];
end
