%dump_fold_memberships_20261026.m -- record-derived training/validation
% memberships for the adopted picks, BEFORE rewriting the identification
% prose (Ben: the assumed mappings are wrong; recompute from the mats).
% READ-ONLY: loads Testing_Data/2022_02_Festo result mats and prints the
% fold columns of the adopted rows. Flexor: minimizeFlxPin10_results.mat
% row 77 (leading columns of filtered_results = held-out validation test
% indices of that member's fold). Extensor: minimizeExt10mmX3_results.mat
% row 32, holdCols=[2 3 4]. Family/test enumerations printed alongside.
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');
setenv('FLX_T5_YMM', '');

%% flexor mat
S = load('minimizeFlxPin10_results.mat');
fprintf('=== minimizeFlxPin10_results.mat ===\n');
fprintf('fields: %s\n', strjoin(fieldnames(S).', ', '));
if isfield(S, 'filtered_results')
    FR = S.filtered_results;
    fprintf('filtered_results: %d x %d\n', size(FR,1), size(FR,2));
    fprintf('row 77 (all cols): '); fprintf(' %g', FR(77,:)); fprintf('\n');
    fprintf('rows 76-78 leading 8 cols for context:\n');
    for r = 76:78
        fprintf('  r%02d: ', r); fprintf(' %g', FR(r,1:min(8,end))); fprintf('\n');
    end
end
if isfield(S, 'xCols')
    fprintf('xCols: %s\n', mat2str(S.xCols));
    fprintf('row77 xCols values: %s\n', mat2str(S.filtered_results(77, S.xCols), 10));
end
for f = {'testNames','test_names','tests','family','fold','folds','holdCols','LOCKSRC'}
    if isfield(S, f{1}), fprintf('%s: ', f{1}); disp(S.(f{1})); end
end
% the flexor family enumeration, from the evaluator: 1=48.5 2=45.7 3=47.9
% 4=40.6t 5=41.0 (gen_FlxPin_allmodels comment). Print the ke/kf campaign
% resting lengths if present.
if isfield(S, 'kf')
    try
        for i = 1:numel(S.kf)
            fprintf('kf(%d): rest=%.4f ten=%.4f\n', i, S.kf(i).RestingL, S.kf(i).Ten);
        end
    catch
        fprintf('kf present but not a struct with RestingL\n');
    end
end

%% extensor mat
E = load('minimizeExt10mmX3_results.mat');
fprintf('\n=== minimizeExt10mmX3_results.mat ===\n');
fprintf('fields: %s\n', strjoin(fieldnames(E).', ', '));
holdCols = [2 3 4];
if isfield(E, 'filtered_results')
    FR2 = E.filtered_results;
    fprintf('filtered_results: %d x %d\n', size(FR2,1), size(FR2,2));
    fprintf('row 32 (all cols): '); fprintf(' %g', FR2(32,:)); fprintf('\n');
    fprintf('row 32 holdCols: %s\n', mat2str(FR2(32, holdCols)));
    fprintf('rows 31-33 leading 10 cols for context:\n');
    for r = 31:33
        fprintf('  r%02d: ', r); fprintf(' %g', FR2(r,1:min(10,end))); fprintf('\n');
    end
end
if isfield(E, 'xCols')
    fprintf('xCols: %s\n', mat2str(E.xCols));
    fprintf('row32 xCols values: %s\n', mat2str(E.filtered_results(32, E.xCols), 10));
end
for f = {'testNames','test_names','tests','family','fold','folds','holdCols','LOCKSRC','usedTests'}
    if isfield(E, f{1}), fprintf('%s: ', f{1}); disp(E.(f{1})); end
end
if isfield(E, 'ke')
    try
        for i = 1:numel(E.ke)
            fprintf('ke(%d): rest=%.4f ten=%.4f\n', i, E.ke(i).RestingL, E.ke(i).Ten);
        end
    catch
        fprintf('ke present but not a struct with RestingL\n');
    end
end
disp('DONE_dump_folds');
