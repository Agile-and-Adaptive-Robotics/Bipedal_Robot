function nightBatch_20260908()
%NIGHTBATCH_20260908 STAGED 2026-09-08 by laptop zcode session: NEVER EXECUTED.
% Owner going forward: easteregg2 zcode chat (see HANDOFF_laptop_20260908.md).
% Rename/erase/adapt freely. Self-locating paths: safe on either machine.
%
% Queued sequence (Ben's asks, 2026-09-08):
%   1 picksScan on the 20260907 2trans SMOKE results (different picks; ~15 min)
%   2 full CV with the 1TRANS evaluator = rerun of the run killed early 09-07
%   3 crossPredictFlx + picksScan on the fresh 1trans full result
%   4 full CV with the 2TRANS evaluator (method of record, production-grade)
%   5 crossPredictFlx + picksScan on the fresh 2trans full result
%
% Solver: surrogateopt for the full runs (A/B on 09-07: gamultiobj 333 s /
% dist 0.309 vs surrogateopt 231 s / dist 0.350; the abSolver.m rule picks
% surrogateopt). Full gamultiobj CV at POP 150 x MAXGEN 600 x 10 folds is a
% DAY-length run on the laptop's 6 workers - do not launch casually; easteregg2
% (10 workers) is the place for it if Ben wants a gamultiobj production run.
% Expected wall time for this batch: roughly 1-2 h per full CV + ~15 min per
% picks scan (6 workers), so plan ~4-5 h total.

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
cd(here);
addpath(fullfile(root, 'Code', 'Matlab', 'Functions'));
addpath(fullfile(root, 'Code', 'Matlab', 'Functions', 'ModernRobotics'));
addpath(fullfile(root, 'Code', 'Matlab', 'Robot_Data'));
t00 = tic;
batchLog = struct('step', {}, 'note', {}, 'minutes', {});
smokeMat = fullfile(here, 'minimizeFlxPin10_results_20260907_2brkt_2trans_smoke.mat');
stamp = char(string(datetime('now'), 'yyyyMMdd'));

%% Step 1: picks on the 2trans smoke results
t0 = tic;
try
    picksScan_2brkt(smokeMat, 12);
    batchLog(end+1) = struct('step', 1, 'note', 'picksScan on 2trans smoke', ...
        'minutes', toc(t0)/60); %#ok<SAGROW>
catch ME
    batchLog(end+1) = struct('step', 1, 'note', sprintf('FAILED: %s', ME.message), ...
        'minutes', toc(t0)/60); %#ok<SAGROW>
    fprintf('STEP 1 FAILED: %s\n', ME.message);
end

%% Step 2: full 1trans CV (rerun of the killed run; script runs in base workspace)
t0 = tic;
try
    evalin('base', sprintf(['setenv(''FLX2BRK_MODE'',''full''); ', ...
        'setenv(''FLX2BRK_SOLVER'',''surrogateopt''); cd(''%s''); ', ...
        'minimizeFlxPin10mm_2brk_1trans'], here));
    batchLog(end+1) = struct('step', 2, 'note', 'full 1trans CV done (surrogateopt)', ...
        'minutes', toc(t0)/60); %#ok<SAGROW>
catch ME
    batchLog(end+1) = struct('step', 2, 'note', sprintf('FAILED: %s', ME.message), ...
        'minutes', toc(t0)/60); %#ok<SAGROW>
    fprintf('STEP 2 FAILED: %s\n', ME.message);
end

%% Step 3: cross-prediction + picks on the fresh 1trans full result
t0 = tic;
try
    f1 = newestResult(here, '*_2brkt_1trans.mat');
    crossPredictFlx(f1, 1);
    renameTodaysCrossPredict(here, stamp, '_1trans_full');
    picksScan_2brkt(f1, 8);
    batchLog(end+1) = struct('step', 3, 'note', sprintf('crossPredict+picks on %s', ...
        char(string(f1))), 'minutes', toc(t0)/60); %#ok<SAGROW>
catch ME
    batchLog(end+1) = struct('step', 3, 'note', sprintf('FAILED: %s', ME.message), ...
        'minutes', toc(t0)/60); %#ok<SAGROW>
    fprintf('STEP 3 FAILED: %s\n', ME.message);
end

%% Step 4: full 2trans CV (method of record)
t0 = tic;
try
    evalin('base', sprintf(['setenv(''FLX2BRK_MODE'',''full''); ', ...
        'setenv(''FLX2BRK_SOLVER'',''surrogateopt''); cd(''%s''); ', ...
        'minimizeFlxPin10mm_2brk'], here));
    batchLog(end+1) = struct('step', 4, 'note', 'full 2trans CV done (surrogateopt)', ...
        'minutes', toc(t0)/60); %#ok<SAGROW>
catch ME
    batchLog(end+1) = struct('step', 4, 'note', sprintf('FAILED: %s', ME.message), ...
        'minutes', toc(t0)/60); %#ok<SAGROW>
    fprintf('STEP 4 FAILED: %s\n', ME.message);
end

%% Step 5: cross-prediction + picks on the fresh 2trans full result
t0 = tic;
try
    f2 = newestResult(here, '*_2brkt_2trans.mat');
    crossPredictFlx(f2, 1);
    renameTodaysCrossPredict(here, stamp, '_2trans_full');
    picksScan_2brkt(f2, 8);
    batchLog(end+1) = struct('step', 5, 'note', sprintf('crossPredict+picks on %s', ...
        char(string(f2))), 'minutes', toc(t0)/60); %#ok<SAGROW>
catch ME
    batchLog(end+1) = struct('step', 5, 'note', sprintf('FAILED: %s', ME.message), ...
        'minutes', toc(t0)/60); %#ok<SAGROW>
    fprintf('STEP 5 FAILED: %s\n', ME.message);
end

%% Wrap up
fprintf('\n########## NIGHTBATCH SUMMARY (%.1f h total) ##########\n', toc(t00)/3600);
for i = 1:numel(batchLog)
    fprintf('  Step %d (%6.1f min): %s\n', batchLog(i).step, batchLog(i).minutes, batchLog(i).note);
end
save(fullfile(here, ['nightBatch_summary_', stamp, '.mat']), 'batchLog');
end

function f = newestResult(here, pattern)
%Newest non-smoke match (smoke files end in _smoke.mat, so these patterns skip them)
    d = dir(fullfile(here, ['minimizeFlxPin10_results_', pattern]));
    assert(~isempty(d), 'no %s result found (earlier step failed?)', pattern);
    [~, ix] = max([d.datenum]);
    f = fullfile(d(ix).folder, d(ix).name);
    fprintf('using result file: %s\n', f);
end

function renameTodaysCrossPredict(here, stamp, suffix)
%crossPredictFlx always saves crossPredict_<stamp>.mat - rename before the next
%call overwrites it (same-day collision gotcha).
    src = fullfile(here, ['crossPredict_', stamp, '.mat']);
    if isfile(src)
        dst = fullfile(here, ['crossPredict_', stamp, suffix, '.mat']);
        movefile(src, dst, 'f');
        fprintf('renamed to %s\n', dst);
    end
end
