% Dig_OptRun_frommat_ext_20260921.m
% Ben's workflow, verified FAST: load Results\Vas_Pam_20mm_Result.mat
% (now the 1519 content) into the BASE workspace, execute Opt_run_Ext.m
% from the '%% Run with adjusted seed' section to end of file. Per Ben's
% instruction the refine patternsearch is capped for THIS run only:
% optsP is pre-set in the base workspace (MaxFunctionEvaluations 50,
% serial) so the section's ~exist guard no-ops and the capped options
% are used. Asserts no new dated result mat is minted.

thisDir = fileparts(mfilename('fullpath'));
root = thisDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root)
        error('Bipedal_Robot repo root not found from %s', thisDir)
    end
    root = parent;
end
meshDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization');
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(meshDir);
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));

resDir = fullfile(meshDir, 'Results');
matFile = fullfile(resDir, 'Vas_Pam_20mm_Result.mat');
before = {dir(fullfile(resDir, 'Vas_Pam_20mm_Result_2*.mat')).name};

fprintf('=== loading %s ===\n', matFile)
load(matFile)   % plain load -> base workspace, like Ben's workflow
fprintf('loaded: fBest = %.6g, XiUsed = [%.6g %.6g %.6g %.6g]\n', ...
    fBest, XiUsed)

% Verification-only patternsearch cap (Ben, 2026-09-21): section guard
% sees optsP already defined and keeps these options.
optsP = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'UseParallel', false, ...
    'MaxFunctionEvaluations', 50, ...
    'MeshTolerance', 1e-4, ...
    'StepTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-6);

% Extract the section from the driver and stage it as a runnable script.
L = splitlines(string(fileread(fullfile(meshDir, 'Opt_run_Ext.m'))));
i0 = find(L == "%% Run with adjusted seed", 1);
assert(~isempty(i0), 'section marker not found in Opt_run_Ext.m')
tmpF = fullfile(tempdir, 'frommat_ext_section_20260921.m');
fid = fopen(tmpF, 'w');
fwrite(fid, char(join(L(i0:end), newline)));
fclose(fid);
addpath(tempdir);

fprintf('=== running Opt_run_Ext adjusted-seed section -> EOF from mat ===\n')
frommat_ext_section_20260921

after = {dir(fullfile(resDir, 'Vas_Pam_20mm_Result_2*.mat')).name};
assert(isequal(before, after), ...
    'a new dated result mat was minted by the from-mat run')
fprintf('FROMMAT EXT PASS (50-eval cap): section ran clean; no new dated mat minted.\n')
