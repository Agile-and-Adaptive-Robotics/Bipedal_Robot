%% run_FlxPin10mm_K2allX1_20260912.m
% Launch wrapper for minimizeFlxPin10mm with:
%   - K2 = [X1,X1,X1] (isotropic origin bracket, Ben 2026-09-12)
%   - exactly TWO folds: train {1,2,5}/validate {3,4}; train {3,4,5}/validate {1,2}
%   - Xi1 > Xi2 constraint still active (nonlcon2 in the gamultiobj call)
% SAVE BEHAVIOR (Ben, 2026-09-12): the ENTIRE driver workspace is saved to the
% results mat (save with no variable list) so that after `load <mat>` any
% section of minimizeFlxPin10mm can be Run Section'd / Run-and-Advanced.
% Result mat: minimizeFlxPin10_results_20260912_2brkt_2trans_K2allX1_2folds.mat
% Log: Dig_out\FlxPin10mm_2brk_K2allX1_20260912.log

here = fileparts(mfilename('fullpath'));              % ...\Dig_out
root = fileparts(fileparts(fileparts(here)));         % ...\Bipedal_Robot
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));
addpath(pwd);   % so pool workers find the data mats via path regardless of worker cwd

logfile = fullfile(here, 'FlxPin10mm_2brk_K2allX1_20260912.log');
diary(logfile); diary on;
fprintf('=== run_FlxPin10mm_K2allX1 | %s ===\n', string(datetime('now')));

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

% Provenance notes (included in the workspace dump)
NOTES = {
 'evaluator: minimizeFlxPin.m (2-bracket, 2trans frames, +5.3 deg on test 4 experimental angles)'
 'K2 = [X1,X1,X1] origin bracket ISOTROPIC (Ben, 2026-09-12): tibia bracket ~15mm y-defl. at'
 '  low force; femur bracket ~25mm only at much higher force -> origin bracket stiff, isotropic'
 'K = [X1,X2,X1] insertion bracket (unchanged)'
 'driver folds (Ben, 2026-09-12): exactly two -- train {1,2,5}/validate {3,4}; train {3,4,5}/validate {1,2}'
 'Xi1 > Xi2 nonlinear constraint ACTIVE (nonlcon2 in the gamultiobj call)'
 'bounds: Xi0 [0,2] cm, Xi1 [5e3,5e7], Xi2 [5e3,5e7] (legacy wide bounds)'
 ['driver error: ' drvErr]
};

saveFile = 'minimizeFlxPin10_results_20260912_2brkt_2trans_K2allX1_2folds.mat';
save(saveFile);   %ENTIRE workspace (Ben, 2026-09-12: load-then-Run-Section must work)
fprintf('Saved %s (entire workspace, %d vars)\n', saveFile, numel(whos));
diary off;
