% verify_cleanmats_20261018.m - do the clean-named mats carry the adopted picks?
% Read-only. Flexor: minimizeFlxPin10_results.mat row 77 vs adopted
% [0.00393737, 39978.7, 14733.8]. Extensor: minimizeExt10mmX3_results.mat
% row 32 vs [-0.00637802, 39978.7, 14733.8, 0.158289].
root = 'D:\Github\Bipedal_Robot';
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));

S = load('minimizeFlxPin10_results.mat');
fprintf('CLEAN FLX fields: %s\n', strjoin(fieldnames(S)', ', '));
if isfield(S, 'filtered_results')
    fprintf('CLEAN FLX filtered_results size: %s\n', mat2str(size(S.filtered_results)));
    if isfield(S, 'xCols'), fprintf('CLEAN FLX xCols = %s\n', mat2str(S.xCols)); end
    fprintf('CLEAN FLX row 77  = %s\n', mat2str(S.filtered_results(77, :), 8));
    fprintf('CLEAN FLX row 107 = %s\n', mat2str(S.filtered_results(107, :), 8));
    fprintf('CLEAN FLX unique (col2,col3) = \n'); disp(unique(S.filtered_results(:, 2:3), 'rows'));
end
if isfield(S, 'TAG'), fprintf('CLEAN FLX TAG = %s\n', mat2str(S.TAG)); end
if isfield(S, 'TRANSMODE'), fprintf('CLEAN FLX TRANSMODE = %s\n', mat2str(S.TRANSMODE)); end

E = load('minimizeExt10mmX3_results.mat');
fprintf('\nCLEAN EXT fields: %s\n', strjoin(fieldnames(E)', ', '));
if isfield(E, 'filtered_results')
    fprintf('CLEAN EXT filtered_results size: %s\n', mat2str(size(E.filtered_results)));
    if isfield(E, 'xCols'), fprintf('CLEAN EXT xCols = %s\n', mat2str(E.xCols)); end
    fprintf('CLEAN EXT row 32 = %s\n', mat2str(E.filtered_results(32, :), 8));
    fprintf('CLEAN EXT unique (col2,col3,col4) = \n'); disp(unique(E.filtered_results(:, 2:4), 'rows'));
end
for f = {'LOCKSRC', 'xi1lock', 'xi2lock', 'flxRow', 'baselineScores', 'labels', 'holdTests', 'trainTests'}
    if isfield(E, f{1})
        v = E.(f{1});
        if isnumeric(v), fprintf('CLEAN EXT %s = %s\n', f{1}, mat2str(v, 8));
        else, fprintf('CLEAN EXT %s = %s\n', f{1}, mat2str(v)); end
    end
end
fprintf('\nCLEAN MAT VERIFY DONE\n');
