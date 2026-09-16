% run_FlxFrX3_2brkt_locks_20260915.m
% Four CV runs of minimizeFlxPin10mmX3_2brkt (two-bracket evaluator):
%   Xi1/Xi2 LOCKED to each mat's pair (GA searches only Xi0 and Xi3):
%     L1trans = pick pair of minimizeFlxPin10_results_20260908_2brkt_1trans_noT3.mat
%     L107    = lock pair carried in minimizeExt10mmX3_results_20260914_pk107lock.mat
%   x tangency mode: DIAG vs ENFORCE (default).
% Each run saves its FULL workspace (Ben 2026-09-12). Driver crashes (e.g. the
% pick section on an empty ENFORCE front) are caught and the partial state saved.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');

m1 = load('minimizeFlxPin10_results_20260908_2brkt_1trans_noT3.mat', 'k2', 'k3');
m2 = load('minimizeExt10mmX3_results_20260914_pk107lock.mat', 'sol_actual');
locks = { ...
    'L1trans', m1.k2, m1.k3; ...
    'L107',    m2.sol_actual(2), m2.sol_actual(3)};
fprintf('Lock pairs: L1trans = (%.10g, %.10g) | L107 = (%.10g, %.10g)\n', ...
    locks{1,2}, locks{1,3}, locks{2,2}, locks{2,3});

pool = gcp('nocreate');
if isempty(pool), parpool(8); end

order = [2 1];                 % L107 first, then L1trans
modes = {'DIAG', 'ENF'};
for mi = 1:2                   % DIAG runs first, ENFORCE after
    for L = order
        runOne(locks{L,1}, locks{L,2}, locks{L,3}, modes{mi});
    end
end
fprintf('\nALL RUNS DONE\n');

function runOne(lockTagIn, lk1, lk2, modeTagIn)
    setenv('FX3B_LOCK1', num2str(lk1, '%.17g'));
    setenv('FX3B_LOCK2', num2str(lk2, '%.17g'));
    setenv('FX3B_LOCKTAG', lockTagIn);
    setenv('FX3B_MODETAG', modeTagIn);
    if strcmp(modeTagIn, 'DIAG')
        setenv('FLXPX3_TANGENCY', 'DIAG');
    else
        setenv('FLXPX3_TANGENCY', '');
    end
    RUN_ERROR = '';
    RUN_STACK = '';
    try
        run('minimizeFlxPin10mmX3_2brkt.m');   % driver clears THIS function's workspace only
    catch ME
        % driver's clear wiped the arguments -- identity comes from env, not them
        RUN_ERROR = ME.message;
        RUN_STACK = strjoin(arrayfun(@(s) sprintf('%s L%d', s.name, s.line), ME.stack, ...
            'UniformOutput', false), ' | ');
        fprintf('\nRUN %s_%s raised: %s\n    at: %s\n', ...
            getenv('FX3B_LOCKTAG'), getenv('FX3B_MODETAG'), RUN_ERROR, RUN_STACK);
    end
    % driver final (or partial) state is in this function's workspace now;
    % re-derive identity from env because the driver cleared the arguments
    lockTag = getenv('FX3B_LOCKTAG');
    modeTag = getenv('FX3B_MODETAG');
    RESULTFILE = sprintf('minimizeFlxPin10mmX3_2brkt_results_20260915_%s_%s.mat', lockTag, modeTag);
    SRCMAT1 = 'minimizeFlxPin10_results_20260908_2brkt_1trans_noT3.mat';
    SRCMAT2 = 'minimizeExt10mmX3_results_20260914_pk107lock.mat';
    LOCKPAIR = [sscanf(getenv('FX3B_LOCK1'), '%g'), sscanf(getenv('FX3B_LOCK2'), '%g')];
    TANGMODE = modeTag;
    save(RESULTFILE);   % full workspace
    fprintf('\nSAVED FULL WORKSPACE: %s\n', RESULTFILE);
end
