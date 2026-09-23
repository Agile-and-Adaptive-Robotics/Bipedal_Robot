% run_Opt_flxr77_20260920.m
% Ben, 2026-09-20: flexor Opt_run (2brkt 2trans noT3 row-77 pick) then
% extensor Opt_run_Ext (flxr77 front row-32 pick), both at the 5% torque
% margin target with 7000 surrogateopt evals. Sequential in one session,
% single log. Launched headless via matlab -batch from this file.
% Both drivers begin with `clear`, which wipes this runner's variables --
% each stage re-derives what it needs from mfilename. The diary and the
% parallel pool survive a clear.

logFile = 'Opt_flxr77_20260920.log';
diary(logFile); diary on;
fprintf('[%s] ===== run_Opt_flxr77_20260920 started =====\n', string(datetime('now')));

%% Preflight: paths + ctx builds so both picks are verified before the
% multi-hour compute. The drivers rebuild ctx themselves.
% mfilename('fullpath') INCLUDES this script's filename -- strip it first,
% then walk up to the repo root.
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
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
% Mesh_Optimization must win any shadowing contest against data subfolders.
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));

ctxF = buildKneeFlexorContext20mm();
fprintf('PREFLIGHT FLX ctx Xi = [%.6g %.6g %.6g %.6g] margin %.4g\n', ...
    ctxF.Xi0, ctxF.Xi1, ctxF.Xi2, ctxF.Xi3, ctxF.requiredTorqueMargin)
ctxE = buildKneeExtContext20mm();
fprintf('PREFLIGHT EXT ctx Xi = [%.6g %.6g %.6g %.6g] margin %.4g\n', ...
    ctxE.Xi0, ctxE.Xi1, ctxE.Xi2, ctxE.Xi3, ctxE.requiredTorqueMargin)
clear ctxF ctxE

%% Stage 1: flexor (Opt_run.m, row-77 ctx)
% Call the driver BY NAME, never run(<fullpath>): run() cds to the
% script's folder while it executes, which breaks the drivers' relative
% reads (OpenSim_*_Results.txt, the front mats). By name, the script is
% found on the path and executed in place with cwd = 2022_02_Festo, which
% is how it is normally launched.
try
    root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    addpath(genpath(fullfile(root, 'Code', 'Matlab')));
    addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
    cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));
    fprintf('[%s] --- Flexor Opt_run (pick 77) starting, cwd = %s ---\n', ...
        string(datetime('now')), string(pwd))
    Opt_run
    fprintf('[%s] --- Flexor Opt_run DONE ---\n', string(datetime('now')))
catch MEf
    fprintf(2, '[%s] FLEXOR RUN ERROR:\n%s\n', string(datetime('now')), getReport(MEf))
end

%% Stage 2: extensor (Opt_run_Ext.m, flxr77 row-32 ctx)
try
    root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    addpath(genpath(fullfile(root, 'Code', 'Matlab')));
    addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
    cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));
    fprintf('[%s] --- Extensor Opt_run_Ext (flxr77 pick 32) starting, cwd = %s ---\n', ...
        string(datetime('now')), string(pwd))
    Opt_run_Ext
    fprintf('[%s] --- Extensor Opt_run_Ext DONE ---\n', string(datetime('now')))
catch MEe
    fprintf(2, '[%s] EXTENSOR RUN ERROR:\n%s\n', string(datetime('now')), getReport(MEe))
end

fprintf('[%s] ===== run_Opt_flxr77_20260920 finished =====\n', string(datetime('now')))
diary off
