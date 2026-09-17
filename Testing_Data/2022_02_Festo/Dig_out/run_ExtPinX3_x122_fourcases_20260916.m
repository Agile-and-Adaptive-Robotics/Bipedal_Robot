% run_ExtPinX3_x122_fourcases_20260916.m
% Four single-fold extensor CVs with K = [X1,X2,X2] (one bracket):
%   pick1_h18         lock = 20260910 front pick 1 pair (4.35e4/1.70e4), holdout {1,8}
%   pick107_h18       lock = pick 107 pair (5.62e4/1.85e4),        holdout {1,8}
%   pick1_all_h3479   lock = pick 1 pair,  allBPA = all 9 tests,   holdout {3,4,7,9}
%   pick107_all_h3479 lock = pick 107 pair, allBPA = all 9 tests,  holdout {3,4,7,9}
% The GA solves Xi0 and Xi3 (EXTX3_PASS=2). Full-workspace saves.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');

m1 = load('minimizeExt10mmX3_results_20260910_noT3.mat', 'sol_actual');
m107 = load('minimizeExt10mmX3_results_20260914_pk107lock.mat', 'sol_actual');
cases = { ...
    'pick1_h18',         m1.sol_actual(2),   m1.sol_actual(3),   '1,8',     ''; ...
    'pick107_h18',       m107.sol_actual(2), m107.sol_actual(3), '1,8',     ''; ...
    'pick1_all_h3479',   m1.sol_actual(2),   m1.sol_actual(3),   '3,4,7,9', '1'; ...
    'pick107_all_h3479', m107.sol_actual(2), m107.sol_actual(3), '3,4,7,9', '1'};
fprintf('pick1 pair: (%.10g, %.10g) | pick107 pair: (%.10g, %.10g)\n', ...
    m1.sol_actual(2), m1.sol_actual(3), m107.sol_actual(2), m107.sol_actual(3));

pool = gcp('nocreate');
if isempty(pool), parpool(8); end

for c = 1:size(cases, 1)
    runOne(cases{c,1}, cases{c,2}, cases{c,3}, cases{c,4}, cases{c,5});
end
fprintf('\nALL FOUR RUNS DONE\n');

function runOne(tagEnv, lk1, lk2, hold, allt)
    setenv('EXTX3_PASS', '2');
    setenv('EXTX3_XI1', num2str(lk1, '%.17g'));
    setenv('EXTX3_XI2', num2str(lk2, '%.17g'));
    setenv('EXTX3_HOLD', hold);
    setenv('EXTX3_ALLTESTS', allt);
    setenv('EXTX3_RUNTAG', tagEnv);
    RUN_ERROR = '';
    try
        run('minimizeExt10mmX3.m');   % driver clears THIS function's workspace only
    catch ME
        RUN_ERROR = ME.message;
        fprintf('\nRUN %s raised: %s\n', getenv('EXTX3_RUNTAG'), RUN_ERROR);
    end
    RESULTFILE = sprintf('minimizeExt10mmX3_results_20260916_%s.mat', getenv('EXTX3_RUNTAG'));
    LOCKSRC = 'pick1 = 20260910 front pick-1 pair; pick107 = 2trans_noT3 pick-107 pair';
    KMODE = 'K = [X1, X2, X2], single bracket';
    save(RESULTFILE);   % full workspace
    fprintf('\nSAVED FULL WORKSPACE: %s\n', RESULTFILE);
end
