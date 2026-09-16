% run_ExtPinX3_pk107lock_20260914.m
% Runner: pinned extensor CV with Xi1/Xi2 locked to flexor 2trans_noT3 pick 107
% (row 107: Xi1 = 56237.84753, Xi2 = 18542.25769) via EXTX3_PASS=2 env lock.
% Saves FULL workspace when done (Ben directive 2026-09-12).

cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath(genpath('D:/GitHub/Bipedal_Robot/Code/Matlab'));
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization');   % win shadowing contests

pool = gcp('nocreate');
if isempty(pool), parpool(8); end   % EB475WS4 cap ~6-8 workers

% Driver clears the workspace at its start and runs to completion in this
% workspace; RESULTFILE is (re)assigned AFTER run() so it survives the clear.
run('minimizeExt10mmX3.m');

RESULTFILE  = 'minimizeExt10mmX3_results_20260914_pk107lock.mat';
EXTPASS     = getenv('EXTX3_PASS');
EXTXI1LOCK  = getenv('EXTX3_XI1');
EXTXI2LOCK  = getenv('EXTX3_XI2');
LOCKSRC     = 'minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat filtered_results row 107 (pick 107)';
RUNTAG      = 'pk107lock';
save('minimizeExt10mmX3_results_20260914_pk107lock.mat');   % full workspace
fprintf('\nSAVED FULL WORKSPACE: %s\n', RESULTFILE);
