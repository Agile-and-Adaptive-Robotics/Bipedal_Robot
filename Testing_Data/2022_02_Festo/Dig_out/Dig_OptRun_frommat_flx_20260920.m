% Dig_OptRun_frommat_flx_20260920.m
% Reproduce Ben's display-from-mat workflow for the flexor driver:
% load Results\Bifemsh_20mm_Result.mat into the BASE workspace, then
% execute Opt_run.m from the '%% Run with adjusted seed' section to the
% end of the file (local functions included so the extracted text is a
% valid standalone script). Verifies the from-mat guards and that no new
% dated result mat gets minted. Read-only for Results\.

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
matFile = fullfile(resDir, 'Bifemsh_20mm_Result.mat');
before = {dir(fullfile(resDir, 'Bifemsh_20mm_Result_2*.mat')).name};

fprintf('=== loading %s ===\n', matFile)
load(matFile)   % plain load -> base workspace, like Ben's workflow
fprintf('loaded: fBest = %.6g, XiUsed = [%.6g %.6g %.6g %.6g]\n', ...
    fBest, XiUsed)

% Extract the section from the driver and stage it as a runnable script.
L = splitlines(string(fileread(fullfile(meshDir, 'Opt_run.m'))));
i0 = find(L == "%% Run with adjusted seed", 1);
assert(~isempty(i0), 'section marker not found in Opt_run.m')
tmpF = fullfile(tempdir, 'frommat_flx_section_20260920.m');
fid = fopen(tmpF, 'w');
fwrite(fid, char(join(L(i0:end), newline)));
fclose(fid);
addpath(tempdir);

fprintf('=== running Opt_run adjusted-seed section -> EOF from mat ===\n')
frommat_flx_section_20260920

after = {dir(fullfile(resDir, 'Bifemsh_20mm_Result_2*.mat')).name};
assert(isequal(before, after), ...
    'a new dated result mat was minted by the from-mat run')
fprintf('FROMMAT FLX PASS: section ran clean; no new dated mat minted.\n')
