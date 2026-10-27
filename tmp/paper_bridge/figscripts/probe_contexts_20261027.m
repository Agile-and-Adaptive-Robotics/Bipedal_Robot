% probe_contexts_20261027.m -- staged probe (easteregg2, 2026-10-27)
% Read-only inspection of the data of record for the two studies:
%  (1) MovePoints.mat contents (does the pipeline carry a t1-ICR table?)
%  (2) adopted extensor chi values (minimizeExt10mmX3_results_20260920_flxr77
%      pick 32, xCols order) -- buildKneeExtContext20mm.m's own load
%  (3) cam design of record mat variables (Vas_Pam_20mm_Result_cam_20261014_0416)
%  (4) helper availability on path (festo4, maxBPAforce, RpToTrans)

R = 'D:\GitHub\Bipedal_Robot';

fprintf('--- (1) MovePoints.mat (tmp/old_code_copy_20261022) ---\n');
w = whos('-file', fullfile(R, 'tmp', 'old_code_copy_20261022', 'MovePoints.mat'));
for k = 1:numel(w)
    fprintf('  %-16s %-10s %s\n', w(k).name, w(k).class, mat2str(w(k).size));
end

fprintf('--- (2) adopted extensor chi (flxr77 pick 32) ---\n');
S = load(fullfile(R, 'Testing_Data', '2022_02_Festo', ...
    'minimizeExt10mmX3_results_20260920_flxr77.mat'));
fprintf('  variables: %s\n', strjoin(fieldnames(S)', ', '));
xc = S.xCols(:).';
fprintf('  xCols: [%s]\n', num2str(xc));
g = S.filtered_results(32, xc);
fprintf('  pick 32 g = [%.6g %.6g %.6g %.6g]\n', g);
fprintf('  (Xi0=%.6g Xi1=%.6g Xi2=%.6g Xi3=%.6g)\n', g(1), g(2), g(3), g(4));

fprintf('--- (3) cam result of record ---\n');
C = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'));
fprintf('  variables: %s\n', strjoin(fieldnames(C)', ', '));
if isfield(C, 'x')
    fprintf('  x  = [%s]\n', num2str(C.x(:).', 6));
end
if isfield(C, 'xBest')
    fprintf('  xBest = [%s]\n', num2str(C.xBest(:).', 6));
end
if isfield(C, 'finalX')
    fprintf('  finalX = [%s]\n', num2str(C.finalX(:).', 6));
end
if isfield(C, 'meta'); disp(C.meta); end
if isfield(C, 'result'); disp(C.result); end

fprintf('--- (4) helpers on path ---\n');
addpath(genpath(fullfile(R, 'Code', 'Matlab')));
for fn = {'festo4', 'maxBPAforce', 'RpToTrans', 'RowVecTrans'}
    fprintf('  %-12s -> %s\n', fn{1}, which(fn{1}));
end

fprintf('--- (5) surfaces CSV act 11 grid (vas_int convention) ---\n');
T = readtable(fullfile(R, 'tmp', 'paper_bridge', 'capability_table_20261022', ...
    'human_surfaces_5x5_grading_20261023.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
ix = T.actuator_id == 11;
fprintf('  pri dof %s / sec dof %s; p nodes: %s\n', T.primary_dof{find(ix,1)}, ...
    T.secondary_dof{find(ix,1)}, mat2str(unique(T.p_node_deg(ix))));
fprintf('  s nodes: %s\n', mat2str(unique(T.s_node_deg(ix))));

fprintf('PROBE DONE\n');
