% reeval_bio_20260926.m - bio-extensor (52cm) GoF at pick-32 values + old-value
% reproduction check. Read-only: evaluator calls only, no mats saved.
thisDir = fileparts(mfilename('fullpath'));
root = thisDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot'), break, end
    if strcmp(parent, root), error('repo root not found'), end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));
fprintf('MATLAB %s | %s\n', version, string(datetime('now')));

% exact old adopted extensor values (appendix of record until 09-20):
% front 20260910_noT3 pick 1
E0 = load('minimizeExt10mmX3_results_20260910_noT3.mat', 'filtered_results', 'xCols');
ep = E0.filtered_results(1, E0.xCols);
fprintf('\nold adopted extensor pick 1: Xi0=%.6g Xi1=%.6g Xi2=%.6g Xi3=%.6g\n', ep);

base_bio = minimizeExt(0, Inf, Inf, 0, 1);
old_bio  = minimizeExt(ep(1), ep(2), ep(3), ep(4), 1);
% current of record: flxr77 row 32, pure substitution on the bio 52cm test
new_bio  = minimizeExt(-6.37802211e-03, 39978.72850713, 14733.75228587, 1.58288685e-01, 1);
% pure substitution with CURRENT flexor pair but OLD Xi0/Xi3, for attribution
mix_bio  = minimizeExt(ep(1), 39978.72850713, 14733.75228587, ep(4), 1);

disp('                RMSE       FVU      MaxRes');
fprintf('baseline     %9.4f %9.4f %9.4f\n', base_bio);
fprintf('old adopted  %9.4f %9.4f %9.4f  (expect ~2.20/2.35/4.73)\n', old_bio);
fprintf('pick32 pure  %9.4f %9.4f %9.4f\n', new_bio);
fprintf('old Xi0/Xi3 + new pair %9.4f %9.4f %9.4f\n', mix_bio);
fprintf('BIO DONE (no mats written)\n');
