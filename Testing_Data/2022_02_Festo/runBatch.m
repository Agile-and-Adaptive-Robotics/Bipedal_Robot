function runBatch()
%RUNBATCH Overnight sequence (call from matlab -batch):
%   Step 0: one-fold A/B (abSolver.m) — gamultiobj vs surrogateopt. Winner
%           sets FLX2BRK_SOLVER for the night.
%   Step 1: full 2brk CV  (minimizeFlxPin10mm_2brk, FLX2BRK_MODE='full')
%   Step 2: cross-prediction on the new 2brk solution (crossPredictFlx)
%   Step 3: extensor high-Xi1 sweep (sweepExtX3, trimmed pool)
% Each step is try/catch-isolated so a late failure keeps earlier results.
% Run from bash like:
%   matlab -batch "cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo'); runBatch"

base = 'D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo';
fnPath = 'D:/GitHub/Bipedal_Robot/Code/Matlab/Functions';
mrPath = 'D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics';
t00 = tic;
batchLog = struct('step', {}, 'note', {}, 'minutes', {});

%% Step 0: solver A/B (skipped if FLX2BRK_SOLVER already set in the environment)
if isempty(getenv('FLX2BRK_SOLVER'))
try
    t0 = tic;
    fprintf('\n########## STEP 0: solver A/B (gamultiobj vs surrogateopt) ##########\n');
    [winner, abStats] = abSolver();
    setenv('FLX2BRK_SOLVER', winner);
    batchLog(end+1) = struct('step', 0, 'note', sprintf('A/B winner: %s (ga %.0fs/d=%.3f | surr %.0fs/d=%.3f)', ...
        winner, abStats.gaTime, abStats.gaDist, abStats.surrTime, abStats.surrDist), ...
        'minutes', toc(t0)/60); %#ok<SAGROW>
    save(fullfile(base, 'batchAB_solver.mat'), 'abStats', 'winner');
catch ME
    fprintf('A/B failed (%s) — defaulting to gamultiobj\n', ME.message);
    setenv('FLX2BRK_SOLVER', 'gamultiobj');
    batchLog(end+1) = struct('step', 0, 'note', sprintf('A/B failed: %s', ME.message), 'minutes', toc(t0)/60); %#ok<SAGROW>
end
else
    fprintf('STEP 0 skipped: FLX2BRK_SOLVER="%s" already set (A/B result: gamultiobj 333s/d=0.309 vs surrogateopt 231s/d=0.350 on 2026-09-07)\n', getenv('FLX2BRK_SOLVER'));
end
%% Step 1: full 2brk CV (script runs in the BASE workspace; it clears it)
try
    t0 = tic;
    fprintf('\n########## STEP 1: full 2brk CV ##########\n');
    evalin('base', sprintf( ...
        'cd(''%s''); addpath(''%s''); addpath(''%s''); setenv(''FLX2BRK_MODE'',''full''); minimizeFlxPin10mm_2brk', ...
        base, fnPath, mrPath));
    batchLog(end+1) = struct('step', 1, 'note', 'full 2brk CV done', 'minutes', toc(t0)/60); %#ok<SAGROW>
catch ME
    batchLog(end+1) = struct('step', 1, 'note', sprintf('FAILED: %s', ME.message), 'minutes', toc(t0)/60); %#ok<SAGROW>
    fprintf('STEP 1 FAILED: %s\n', ME.message);
end

%% Step 2: cross-prediction on the fresh 2brk solution
try
    t0 = tic;
    fprintf('\n########## STEP 2: cross-prediction ##########\n');
    d = dir(fullfile(base, 'minimizeFlxPin10_2brk_results_*.mat'));
    if isempty(d)
        error('no minimizeFlxPin10_2brk_results_*.mat found (step 1 output missing?)');
    end
    [~, ix] = max([d.datenum]);
    crossPredictFlx(fullfile(d(ix).folder, d(ix).name), 1);
    batchLog(end+1) = struct('step', 2, 'note', sprintf('crossPredictFlx on %s', d(ix).name), 'minutes', toc(t0)/60); %#ok<SAGROW>
catch ME
    batchLog(end+1) = struct('step', 2, 'note', sprintf('FAILED: %s', ME.message), 'minutes', toc(t0)/60); %#ok<SAGROW>
    fprintf('STEP 2 FAILED: %s\n', ME.message);
end

%% Step 3: extensor sweep
try
    t0 = tic;
    fprintf('\n########## STEP 3: extensor high-Xi1 sweep ##########\n');
    sweepExtX3();   %defaults: pool {1,2,5,6,7,8}, numHold 3 and 2, MAXHOURS 6
    batchLog(end+1) = struct('step', 3, 'note', 'sweepExtX3 done', 'minutes', toc(t0)/60); %#ok<SAGROW>
catch ME
    batchLog(end+1) = struct('step', 3, 'note', sprintf('FAILED: %s', ME.message), 'minutes', toc(t0)/60); %#ok<SAGROW>
    fprintf('STEP 3 FAILED: %s\n', ME.message);
end

%% Wrap up
fprintf('\n########## BATCH SUMMARY (%.1f h total) ##########\n', toc(t00)/3600);
for i = 1:numel(batchLog)
    fprintf('  Step %d (%6.1f min): %s\n', batchLog(i).step, batchLog(i).minutes, batchLog(i).note);
end
save(fullfile(base, sprintf('batch_summary_%s.mat', char(string(datetime('now'),'yyyyMMdd')))), 'batchLog');
end

