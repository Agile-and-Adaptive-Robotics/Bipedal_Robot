% run_ExtPinX3_1translock_20260916.m
% Extensor CV with Xi1/Xi2 locked to the 1trans flexor pair (EXTX3_PASS=2 env);
% GA searches Xi0 and Xi3. Full-workspace save (Ben 2026-09-12).
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');

pool = gcp('nocreate');
if isempty(pool), parpool(8); end

run('minimizeExt10mmX3.m');   % driver clears this workspace at its start

RESULTFILE = 'minimizeExt10mmX3_results_20260916_1translock.mat';
EXTPASS = getenv('EXTX3_PASS');
EXTXI1LOCK = getenv('EXTX3_XI1');
EXTXI2LOCK = getenv('EXTX3_XI2');
LOCKSRC = 'minimizeFlxPin10_results_20260908_2brkt_1trans_noT3.mat filtered_results row 1 (pick 1)';
RUNTAG = '1translock';
save(RESULTFILE);   % full workspace
fprintf('\nSAVED FULL WORKSPACE: %s\n', RESULTFILE);
