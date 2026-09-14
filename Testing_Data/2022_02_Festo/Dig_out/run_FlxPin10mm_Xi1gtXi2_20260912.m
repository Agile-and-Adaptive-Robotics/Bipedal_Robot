%% run_FlxPin10mm_Xi1gtXi2_20260912.m
% Launch wrapper for minimizeFlxPin10mm (legacy driver, as committed) with the
% two-bracket minimizeFlxPin evaluator and the +5.3 deg encoder offset on
% test 4 (40cm-tendon). Result mat captured in the family format:
%   minimizeFlxPin10_results_20260912_2brkt_2trans_offT4_Xi1gtXi2.mat
% Log: Dig_out\FlxPin10mm_2brk_Xi1gtXi2_20260912.log

here = fileparts(mfilename('fullpath'));              % ...\Dig_out
root = fileparts(fileparts(fileparts(here)));         % ...\Bipedal_Robot
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));
addpath(pwd);   % so pool workers find the data mats via path regardless of worker cwd

logfile = fullfile(here, 'FlxPin10mm_2brk_Xi1gtXi2_20260912.log');
diary(logfile); diary on;
fprintf('=== run_FlxPin10mm_Xi1gtXi2 | %s ===\n', string(datetime('now')));

if isempty(gcp('nocreate'))
    parpool(min(10, feature('numcores')));
end

try
    minimizeFlxPin10mm;
    fprintf('=== driver finished normally | %s ===\n', string(datetime('now')));
    drvErr = '';   %assigned AFTER the driver: its opening `clear` wipes anything set earlier
catch ME
    drvErr = getReport(ME);
    fprintf(2, 'ERROR in driver (saving partial results): %s\n', drvErr);
end
if ~exist('drvErr', 'var'), drvErr = ''; end   %belt and braces

% Results capture -- built AFTER the driver call because its `clear` wipes
% workspace variables. Campaign date = run-start date per family convention.
NOTES = {
 'evaluator: minimizeFlxPin.m (2-bracket port from minimizeFlxPin2brk, 2026-09-11)'
 'frames: two-rotation (Z then Y) both brackets; K=[X1,X2,X1], K2=[X1,X1,X2] (2trans, fixed)'
 'Pbri=[-48.11,-107.81,13.8]mm; Pbr2=[-52.61,0,75.06]mm'
 'encoder correction: +5.3 deg on EXPERIMENTAL angles of test 4 (40cm-tendon) ONLY'
 '(the 20260908 2brk campaign had it on test 3 -- wrong test, Ben 2026-09-11)'
 'driver: minimizeFlxPin10mm.m as committed (allBPA=1:5, numHold=2, POP=150, MAXGEN=600)'
 'bounds: Xi0 [0,2] cm, Xi1 [5e3,5e7], Xi2 [5e3,5e7] (legacy wide bounds, NOT re-centered)'
 'NONLINEAR CONSTRAINT (Ben, 2026-09-12): Xi1 > Xi2 via the drivers own nonlcon2 in the'
 '  gamultiobj call (x(3) < x(2) in log10 space) -- the unconstrained offT4 front landed'
 '  at Xi1 ~= Xi2 (pooled pick=1) and Xi2 > Xi1 throughout the fold-1 cluster'
 ['driver error: ' drvErr]
};
saveFile = 'minimizeFlxPin10_results_20260912_2brkt_2trans_offT4_Xi1gtXi2.mat';
vars = {'results_cv','all_candidates','results_sort','results_sort_actual', ...
        'filtered_results','xCols','a0','baselineScores','f','k1','k2','k3', ...
        'allBPA','numHold','labels','NOTES'};
%Plain loop, NOT cellfun/anonymous: who()/exist() inside an anonymous function
%see the anonymous workspace, not this script's, and report everything missing.
saveVars = {};
for iv = 1:numel(vars)
    if exist(vars{iv}, 'var') == 1
        saveVars{end+1} = vars{iv};   %#ok<AGROW>
    end
end
if isempty(saveVars)
    fprintf(2, 'Curated list empty -- dumping ENTIRE workspace as emergency capture.\n');
    save('flxpin_Xi1gtXi2_workspace_dump_20260912.mat');
else
    save(saveFile, saveVars{:});
    fprintf('Saved %s (%d vars: %s)\n', saveFile, numel(saveVars), strjoin(saveVars, ','));
end
diary off;
