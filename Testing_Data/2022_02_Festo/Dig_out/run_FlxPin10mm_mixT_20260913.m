% run_FlxPin10mm_mixT_20260913.m
% Two sequential campaigns on the MIXED-convention evaluator (Ben, 2026-09-13):
%   tibia/insertion bracket ONE-TRANSFORM (pitch-only), K = [X1,X2,X1];
%   hip/origin bracket TWO-ROTATION, K2 = [X1,X1,X2];
%   folds: holdout {1,5} (train 2,3,4) then holdout {3,4} (train 1,2,5) -- none left out;
%   encoder +5.3 deg on kf(3) (Ben's own 2026-09-12 move).
% Run 1: FLX_T5_YMM = 0 (no test-5 pB shift)
% Run 2: FLX_T5_YMM = 5 (test-5 insertion point +5 mm along bracket-frame y,
%        equivalent to thetabrB + ~4.56 deg -- plastic-deformation allowance)
% Full-workspace saves (Ben directive 2026-09-12): load mat -> run any section.
DIGO = 'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo\Dig_out';
FEST = fileparts(DIGO);
REPO = 'D:\GitHub\Bipedal_Robot';
cd(FEST); addpath(FEST); addpath(DIGO);
addpath(genpath(fullfile(REPO, 'Code', 'Matlab')));   %MonoPamDataExplicit + Colors etc.

diary(fullfile(DIGO, 'run_FlxPin10mm_mixT_20260913.log'));
fprintf('=== mixed-convention campaign, start %s ===\n', datestr(now));
if isempty(gcp('nocreate')), parpool(min(10, feature('numcores'))); end

%% ---------------- Run 1: baseline mixed convention ----------------
setenv('FLX_T5_YMM', '0');
try
    minimizeFlxPin10mm;
    drvErr = '';
catch ME
    drvErr = getReport(ME, 'basic');
    fprintf(2, 'RUN1 DRIVER ERROR:\n%s\n', drvErr);
end
if ~exist('drvErr', 'var'), drvErr = ''; end
FEST = 'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo';   %re-set: driver's clear wiped it
NOTES = {
 'evaluator: minimizeFlxPin.m, MIXED convention (Ben, 2026-09-13)'
 'tibia bracket ONE-TRANSFORM (pitch-only) Tkbr=RpToTrans(RkbrZ,Pbri), K=[X1,X2,X1]; pbrBnew 1-rotation line'
 'hip bracket TWO-ROTATION Thbr=RpToTrans(RhbrZ*Ryh'',Pbr2), K2=[X1,X1,X2] (supersedes 2026-09-12 isotropic test)'
 'folds (driver list): holdout {1,5} train {2,3,4}; holdout {3,4} train {1,2,5} -- leave none out'
 'encoder +5.3 deg on kf(3) EXPERIMENTAL angles only (Ben, 2026-09-12 final attribution)'
 'FLX_T5_YMM = 0: NO test-5 pB shift in this mat'
 ['driver error: ' drvErr]
};
saveFile = fullfile(FEST, 'minimizeFlxPin10_results_20260913_2brkt_mixT_2folds.mat');
save(saveFile);
fprintf('RUN 1 saved: %s\n', saveFile);
close all force; clearvars -except DIGO FEST REPO

%% ---------------- Run 2: + test-5 bracket-y shift ----------------
FEST = 'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo';
setenv('FLX_T5_YMM', '5');
try
    minimizeFlxPin10mm;
    drvErr = '';
catch ME
    drvErr = getReport(ME, 'basic');
    fprintf(2, 'RUN2 DRIVER ERROR:\n%s\n', drvErr);
end
if ~exist('drvErr', 'var'), drvErr = ''; end
FEST = 'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo';   %re-set: driver's clear wiped it
NOTES = {
 'evaluator: minimizeFlxPin.m, MIXED convention (Ben, 2026-09-13)'
 'tibia bracket ONE-TRANSFORM (pitch-only) Tkbr=RpToTrans(RkbrZ,Pbri), K=[X1,X2,X1]; pbrBnew 1-rotation line'
 'hip bracket TWO-ROTATION Thbr=RpToTrans(RhbrZ*Ryh'',Pbr2), K2=[X1,X1,X2]'
 'folds (driver list): holdout {1,5} train {2,3,4}; holdout {3,4} train {1,2,5} -- leave none out'
 'encoder +5.3 deg on kf(3) EXPERIMENTAL angles only (Ben, 2026-09-12 final attribution)'
 'FLX_T5_YMM = 5: test-5 (42cm) insertion point +5 mm along INSERTION-BRACKET-frame y (bending axis)'
 '  = permanent ~5mm bend of bracket arm; equivalent thetabrB rotation +4.56 deg; knee-frame shift [-4.997,-0.159,0] mm'
 '  (pure knee-frame +5mm y is unreachable by rotation: point B sits nearly straight above the bracket, theta=91.8 deg)'
 ['driver error: ' drvErr]
};
saveFile = fullfile(FEST, 'minimizeFlxPin10_results_20260913_2brkt_mixT_T5y5mm_2folds.mat');
save(saveFile);
fprintf('RUN 2 saved: %s\n', saveFile);

setenv('FLX_T5_YMM', '0');   %leave the environment clean
fprintf('=== campaign complete, %s ===\n', datestr(now));
diary off;
