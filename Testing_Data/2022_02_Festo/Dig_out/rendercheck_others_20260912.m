% rendercheck_others_20260912.m
% Verify the plot-section edits in minimizeFlxPin10mmX3 / minimizeFlx10mm /
% minimizeExt10mm: checkcode, live runs of the two light drivers, and X3
% plot-section renders from the K2allX1 workspace in both legend branches.
% NOTE: the drivers call `clear` at top, which wipes script variables -- every
% block re-sets the path constants after its run().
DIGO = 'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo\Dig_out';
FEST = 'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo';
repo = 'D:\GitHub\Bipedal_Robot';
cd(FEST);
addpath(FEST); addpath(DIGO);
addpath(genpath(fullfile(repo, 'Code', 'Matlab'))); %Robot_Data etc. for biomimetic evaluators

%% 1) checkcode all touched files
files = {'minimizeFlxPin.m', 'minimizeFlxPin10mm.m', 'minimizeFlxPin10mmX3.m', ...
         'minimizeFlx10mm.m', 'minimizeExt10mm.m'};
for iF = 1:numel(files)
    r = checkcode(files{iF});
    fprintf('checkcode %-24s : %d findings\n', files{iF}, numel(r));
    for k = 1:numel(r), fprintf('   L%d: %s\n', r(k).line, r(k).message); end %#ok<AGROW>
end

%% 2) live run minimizeFlx10mm (fast, no solver)
try
    run(fullfile(FEST, 'minimizeFlx10mm.m'));
    DIGO = 'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo\Dig_out'; %driver's clear wiped these
    figs = struct('fig', {figT, figL, figMA, figS}, 'name', {'Torque','MuscleLength','MomentArm','RelStrain'});
    for k = 1:numel(figs)
        exportgraphics(figs(k).fig, fullfile(DIGO, sprintf('plotcheck_Flx10mm_%s_20260912.png', figs(k).name)), 'Resolution', 100);
    end
    fprintf('Flx10mm live run + export OK\n');
catch ME
    fprintf(2, 'Flx10mm FAILED: %s\n', getReport(ME, 'basic'));
end
close all force;

%% 3) live run minimizeExt10mm (fast, no solver)
FEST = 'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo'; %re-set: Flx10mm's clear wiped it
try
    run(fullfile(FEST, 'minimizeExt10mm.m'));
    DIGO = 'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo\Dig_out'; %driver's clear wiped these
    figs = struct('fig', {figL, figMA, figT, figE}, 'name', {'MuscleLength','MomentArm','Torque','RelStrain'});
    for k = 1:numel(figs)
        exportgraphics(figs(k).fig, fullfile(DIGO, sprintf('plotcheck_Ext10mm_%s_20260912.png', figs(k).name)), 'Resolution', 100);
    end
    fprintf('Ext10mm live run + export OK\n');
catch ME
    fprintf(2, 'Ext10mm FAILED: %s\n', getReport(ME, 'basic'));
end
close all force;

%% 4) X3 plot sections from the K2allX1 workspace -- odd branch (5 tests)
FEST = 'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo';
DIGO = [FEST '\Dig_out'];
matFile = fullfile(FEST, 'minimizeFlxPin10_results_20260912_2brkt_2trans_K2allX1_2folds.mat');
S = load(matFile);
fn = fieldnames(S);
for k = 1:numel(fn), eval([fn{k} ' = S.(fn{k});']); end %#ok<SAGROW>
try
    run(fullfile(DIGO, 'plotsectionsBody_X3_20260912.m'));
    figs = struct('fig', {figTpre, figTpost, figL, figMA, figS}, 'name', {'TorquePre','TorquePost','MuscleLength','MomentArm','RelStrain'});
    for k = 1:numel(figs)
        exportgraphics(figs(k).fig, fullfile(DIGO, sprintf('plotcheck_X3_%s_20260912.png', figs(k).name)), 'Resolution', 100);
    end
    fprintf('X3 odd-branch render OK\n');
catch ME
    fprintf(2, 'X3 odd-branch FAILED: %s\n', getReport(ME, 'basic'));
end
close all force;

%% 5) X3 plot sections -- even branch (4 tests: legend must land in tile (1,2))
FEST = 'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo';
DIGO = [FEST '\Dig_out'];
matFile = fullfile(FEST, 'minimizeFlxPin10_results_20260912_2brkt_2trans_K2allX1_2folds.mat');
S = load(matFile);
fn = fieldnames(S);
for k = 1:numel(fn), eval([fn{k} ' = S.(fn{k});']); end %#ok<SAGROW>
try
    run(fullfile(DIGO, 'plotsectionsBody_X3even_20260912.m'));
    exportgraphics(figTpre, fullfile(DIGO, 'plotcheck_X3even_TorquePre_20260912.png'), 'Resolution', 100);
    fprintf('X3 even-branch render OK\n');
catch ME
    fprintf(2, 'X3 even-branch FAILED: %s\n', getReport(ME, 'basic'));
end
close all force;
disp('RENDERCHECK OTHERS DONE');
