function Collect_batch()
%COLLECT_BATCH Overnight sequence (call from matlab -batch):
%   Step 1: full 2brk CV  (minimizeFlxPin10mm_2brk, FLX2BRK_MODE='full')
%   Step 2: cross-prediction on the new 2brk solution (Dig_crossPredict)
%   Step 3: extensor high-Xi1 sweep (Collect_ExtPinX3_sweep, trimmed pool)
% Solver comes from FLX2BRK_SOLVER ('gamultiobj' default | 'surrogateopt').
% A/B verdict 2026-09-07: gamultiobj kept (faster per quality at full scale,
% keeps the Pareto front; see Dig_out logs). Each step is try/catch-isolated
% so a late failure keeps earlier results.
% Run from bash like:
%   matlab -batch "cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo'); Collect_batch"

base = 'D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo';
fnPath = 'D:/GitHub/Bipedal_Robot/Code/Matlab/Functions';
mrPath = 'D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics';
t00 = tic;
batchLog = struct('step', {}, 'note', {}, 'minutes', {});

if isempty(getenv('FLX2BRK_SOLVER'))
    setenv('FLX2BRK_SOLVER', 'gamultiobj');
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
    Dig_crossPredict(fullfile(d(ix).folder, d(ix).name), 1);
    batchLog(end+1) = struct('step', 2, 'note', sprintf('Dig_crossPredict on %s', d(ix).name), 'minutes', toc(t0)/60); %#ok<SAGROW>
catch ME
    batchLog(end+1) = struct('step', 2, 'note', sprintf('FAILED: %s', ME.message), 'minutes', toc(t0)/60); %#ok<SAGROW>
    fprintf('STEP 2 FAILED: %s\n', ME.message);
end

%% Step 3: extensor sweep
try
    t0 = tic;
    fprintf('\n########## STEP 3: extensor high-Xi1 sweep ##########\n');
    Collect_ExtPinX3_sweep();   %defaults: pool {1,2,5,6,7,8}, numHold 3 and 2, MAXHOURS 6
    batchLog(end+1) = struct('step', 3, 'note', 'Collect_ExtPinX3_sweep done', 'minutes', toc(t0)/60); %#ok<SAGROW>
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

