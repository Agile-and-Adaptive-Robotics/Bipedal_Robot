% run_FlxFrX3_2brkt_x3u_20260915.m
% Two DIAG lock CVs with the UNITLESS Xi3 (wrap-loss, [0,1]) evaluator.
% Xi1/Xi2 locked per mat (FX3B_LOCK env); GA searches Xi0 and Xi3.
% Full-workspace saves with _x3u tag. ENFORCE legs cancelled (Ben 2026-09-15).
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
setenv('FLXPX3_TANGENCY', 'DIAG');

m1 = load('minimizeFlxPin10_results_20260908_2brkt_1trans_noT3.mat', 'k2', 'k3');
m2 = load('minimizeExt10mmX3_results_20260914_pk107lock.mat', 'sol_actual');
locks = { ...
    'L1trans', m1.k2, m1.k3; ...
    'L107',    m2.sol_actual(2), m2.sol_actual(3)};
fprintf('Lock pairs: L1trans = (%.10g, %.10g) | L107 = (%.10g, %.10g)\n', ...
    locks{1,2}, locks{1,3}, locks{2,2}, locks{2,3});

pool = gcp('nocreate');
if isempty(pool), parpool(8); end

order = [2 1];   % L107 first
for L = order
    runOne(locks{L,1}, locks{L,2}, locks{L,3});
end
fprintf('\nALL RUNS DONE\n');

function runOne(lockTagEnv, lk1, lk2)
    setenv('FX3B_LOCK1', num2str(lk1, '%.17g'));
    setenv('FX3B_LOCK2', num2str(lk2, '%.17g'));
    setenv('FX3B_LOCKTAG', lockTagEnv);
    RUN_ERROR = '';
    RUN_STACK = '';
    try
        run('minimizeFlxPin10mmX3_2brkt.m');   % driver clears THIS function's workspace only
    catch ME
        RUN_ERROR = ME.message;
        RUN_STACK = strjoin(arrayfun(@(s) sprintf('%s L%d', s.name, s.line), ME.stack, ...
            'UniformOutput', false), ' | ');
        fprintf('\nRUN %s raised: %s\n    at: %s\n', getenv('FX3B_LOCKTAG'), RUN_ERROR, RUN_STACK);
    end
    lockTag = getenv('FX3B_LOCKTAG');
    RESULTFILE = sprintf('minimizeFlxPin10mmX3_2brkt_results_20260915_%s_DIAG_x3u.mat', lockTag);
    SRCMAT1 = 'minimizeFlxPin10_results_20260908_2brkt_1trans_noT3.mat';
    SRCMAT2 = 'minimizeExt10mmX3_results_20260914_pk107lock.mat';
    LOCKPAIR = [sscanf(getenv('FX3B_LOCK1'), '%g'), sscanf(getenv('FX3B_LOCK2'), '%g')];
    TANGMODE = 'DIAG';
    X3SEMANTICS = 'unitless wrap-loss [0,1]';
    save(RESULTFILE);   % full workspace
    fprintf('\nSAVED FULL WORKSPACE: %s\n', RESULTFILE);
end
