%dump_ext_folds2_20261026.m -- extensor mat fold metadata: translate row-32
% holdCols values (1,2,8) into resting lengths via the mat's own
% enumerations (holdTests, trainTests, grpTests, validLabels, holdoutIdx,
% trainIdx, grpNames, titles/subtitles) + flexor mat bpa resting lengths.
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
addpath(fullfile(root, 'Code/Matlab/Functions'));
cd(festodir);
warning('off', 'all');

E = load('minimizeExt10mmX3_results.mat');
for f = {'holdTests','trainTests','grpTests','validLabels','grpNames', ...
         'holdoutIdx','trainIdx','train','exthold','numHold','titles','subtitles'}
    if isfield(E, f{1})
        v = E.(f{1});
        fprintf('--- %s (%s) ---\n', f{1}, mat2str(size(v)));
        if isnumeric(v) || islogical(v)
            disp(v);
        elseif iscell(v)
            disp(v(:));
        else
            disp(v);
        end
    end
end
% bpa structs: resting lengths of the extensor campaign in index order
if isfield(E, 'bpa')
    try
        for i = 1:numel(E.bpa)
            fprintf('bpa(%d): rest=%.4f ten=%.4f\n', i, E.bpa(i).rest, E.bpa(i).ten);
        end
    catch me
        fprintf('bpa field read failed: %s\n', me.message);
    end
end
if isfield(E, 'allBPA')
    try
        for i = 1:numel(E.allBPA)
            fprintf('allBPA(%d): rest=%.4f ten=%.4f\n', i, E.allBPA(i).rest, E.allBPA(i).ten);
        end
    catch me
        fprintf('allBPA field read failed: %s\n', me.message);
    end
end

S = load('minimizeFlxPin10_results.mat');
if isfield(S, 'bpa')
    for i = 1:numel(S.bpa)
        fprintf('FLX bpa(%d): rest=%.4f ten=%.4f\n', i, S.bpa(i).rest, S.bpa(i).ten);
    end
end
if isfield(S, 'subtitles'), disp('FLX subtitles:'); disp(S.subtitles); end
if isfield(S, 'titles'), disp('FLX titles:'); disp(S.titles); end
if isfield(S, 'labels'), disp('FLX labels:'); disp(S.labels); end
disp('DONE_dump2');
