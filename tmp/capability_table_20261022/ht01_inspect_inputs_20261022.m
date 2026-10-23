% Hand-tuned BPA evaluation, step 1: inspect inputs on eb475ws4 (READ-ONLY).
% Staged 2026-10-22 in tmp/capability_table_20261022. Prints inventories
% and values of the old-code .mats needed for the Connor-vs-Ben value
% comparison, and checks the Curve Fitting Toolbox (fit) that Connor's
% Bi_Hip knee/patella functions need.

root = 'D:\Github\Bipedal_Robot';
oldc = fullfile(root, 'tmp', 'old_code_copy_20261022');
repo = fullfile(root, 'Code', 'Matlab', 'Previous Optimization Code');

fprintf('=== fit availability (Curve Fitting Toolbox) ===\n');
try
    fobj = fit([0;1],[0;1],'cubicspline');
    fprintf('fit cubicspline: OK\n');
catch e
    fprintf('fit cubicspline FAIL: %s\n', e.message);
end
clear fobj

mats = { ...
    fullfile(oldc,'AttachPoints_Robot.mat'), ...
    fullfile(oldc,'MovePoints.mat'), ...
    fullfile(oldc,'PAM_Prop.mat'), ...
    fullfile(oldc,'AttachPoints.mat'), ...
    fullfile(oldc,'MuscProp.mat'), ...
    fullfile(repo,'Human_Data','MovePoints.mat'), ...
    fullfile(repo,'Human_Data','AttachPoints.mat'), ...
    fullfile(repo,'Human_Data','MuscProp.mat'), ...
    fullfile(repo,'Human_Data','PAM_Prop.mat')};

for k = 1:numel(mats)
    fprintf('\n=== %s ===\n', mats{k});
    try
        w = whos('-file', mats{k});
        for v = 1:numel(w)
            fprintf('  %-18s %-10s %-14s', w(v).name, w(v).class, mat2str(w(v).size));
            if strcmp(w(v).class, 'double') && prod(w(v).size) <= 12
                S = load(mats{k}, w(v).name);
                fprintf('  = %s', mat2str(S.(w(v).name), 6));
            end
            fprintf('\n');
        end
    catch e
        fprintf('  FAIL: %s\n', e.message);
    end
end

% Full values of the robot attachment points (small numeric vars)
fprintf('\n=== AttachPoints_Robot.mat values ===\n');
try
    APR = load(fullfile(oldc,'AttachPoints_Robot.mat'));
    fn = fieldnames(APR);
    for v = 1:numel(fn)
        val = APR.(fn{v});
        if isnumeric(val)
            fprintf('%-6s = %s\n', fn{v}, mat2str(val, 8));
        end
    end
catch e
    fprintf('FAIL: %s\n', e.message);
end

fprintf('\n=== PAM_Prop.mat values (old-code copy) ===\n');
try
    PP = load(fullfile(oldc,'PAM_Prop.mat'));
    fn = fieldnames(PP);
    for v = 1:numel(fn)
        val = PP.(fn{v});
        if isnumeric(val)
            fprintf('%-12s = %s\n', fn{v}, mat2str(val, 8));
        else
            fprintf('%-12s : %s %s\n', fn{v}, class(val), mat2str(size(val)));
        end
    end
catch e
    fprintf('FAIL: %s\n', e.message);
end

fprintf('\n=== MovePoints.mat (old-code copy) numeric values ===\n');
try
    MP = load(fullfile(oldc,'MovePoints.mat'));
    fn = fieldnames(MP);
    for v = 1:numel(fn)
        val = MP.(fn{v});
        if isnumeric(val) && prod(size(val)) <= 12
            fprintf('%-14s = %s\n', fn{v}, mat2str(val, 8));
        else
            fprintf('%-14s : %s %s\n', fn{v}, class(val), mat2str(size(val)));
        end
    end
catch e
    fprintf('FAIL (needs cfit toolbox to deserialize?): %s\n', e.message);
end

fprintf('\nINSPECT DONE\n');
