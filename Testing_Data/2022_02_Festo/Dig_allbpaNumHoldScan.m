function Dig_allbpaNumHoldScan()
%DIG_ALLBPANUMHOLDSCAN Ben's questions (2026-09-08), answered from data.
% Laptop-session script (2026-09-08), renamed into the Collect/Dig scheme same
% day. Still UNTESTED end-to-end (never executed); smoke before trusting.
% Superseded sections were REMOVED 2026-09-08 (Ben: easteregg2 owns this
% folder's compute):
%   - old (0) provenance check -> CLOSED: the canonical 09-07 mat is the
%     restored-2trans chain run, and 1trans==2trans is proven for these
%     y-symmetric K arrays, so "which evaluator" is moot.
%   - old (B) angle-lag scan -> RESOLVED: the shifted test is the flexor 47cm
%     (kf(3)); minimizeFlxPin2brk adds +5.3 deg in its kf-build.
%
%   (A) Per-test held-out patterns from the stored full result (best candidate
%       per fold, leave-2-out, 2brk evaluator on all 5 tests).
%   (C) New smoke-scale gamultiobj folds with the CURRENT 2brk evaluator:
%       E1 ALLBPA=[1..5] NUMHOLD=1 (5 folds) - numHold contrast vs (A)
%       E3 ALLBPA=[2,3,4,5] NUMHOLD=2 (6 folds) - does dropping test 1 hurt?
%       E4 holdout [2,3,4] train [1,5] (Ben's split) + reverse [1,5]/[2,3,4]
%       E5 holdout [1,2,3,4] train [5] - train=1 extreme
% Runtime: (A) ~1 min; (C) 14 ga folds at 25x30, ~15-30 min on 6 workers.
% Self-locating paths. Results -> Dig_out\Dig_allbpaNumHoldScan_20260908.mat

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));   % 2022_02_Festo -> Testing_Data -> repo
cd(here);
addpath(fullfile(root, 'Code', 'Matlab', 'Functions'));
addpath(fullfile(root, 'Code', 'Matlab', 'Functions', 'ModernRobotics'));
addpath(fullfile(root, 'Code', 'Matlab', 'Robot_Data'));   %REQUIRED: 2brk kf-build reads MonoPamDataExplicit objects
set(0, 'DefaultFigureVisible', 'off');
warning('off', 'fortz:LengthBalanceMismatch');
if isempty(gcp('nocreate'))
    try
        parpool(6);
    catch
        fprintf('NOTE: no parallel pool (running serial).\n');
    end
end

labels = ["48cm", "46cm", "47cm", "40cm-tendon", "41cm"];   %tests 1..5
lb = [0, log10(5e3), log10(5e3)];
ub = [2, log10(5e7), log10(5e7)];
opts = optimoptions('gamultiobj', 'UseParallel', true, 'Display', 'off', ...
    'PopulationSize', 25, 'MaxGenerations', 30, ...
    'MutationFcn', {@mutationadaptfeasible}, 'CrossoverFraction', 0.8, ...
    'CrossoverFcn', {@crossoverscattered}, 'FunctionTolerance', 4e-3);
fullMat = 'minimizeFlxPin10_results_20260907_2brkt_2trans.mat';  %canonical 09-07 full CV
prov = '2trans-verified';   %restored-2trans chain run; 1trans==2trans proven for these arrays
Sm = load(fullMat);
a0_5 = minimizeFlxPin2brk(0, Inf, Inf, [], true);   %baseline on all 5 tests

%% ---------- (A) per-test held-out patterns from the stored full result ----------
fprintf('\n===== (A) STORED FULL RESULT: leave-2-out, best cand per fold (2brk evaluator) =====\n');
fprintf('%-16s %7s %7s %7s %7s %7s\n', 'holdout', labels);
nF = numel(Sm.results_cv);
normA = zeros(nF, 5);
holdA = false(nF, 5);
distA = zeros(nF, 1);
for k = 1:nF
    rc = Sm.results_cv{k};
    [~, b] = min(rc.distance_all);
    xr = rc.optParams_all(b, :);
    F = minimizeFlxPin2brk(xr(1)/100, 10^xr(2), 10^xr(3), [], true);
    normA(k, :) = (F(:, 1) ./ a0_5(:, 1))';
    holdA(k, :) = ismember(1:5, unique(rc.foldIdx(b, :)));
    distA(k) = rc.distance_all(b);
    fprintf('%-16s', strjoin(labels(holdA(k, :)), '+'));
    fprintf(' %7.2f', normA(k, :));
    fprintf('\n');
end
fprintf('Per-test held-out mean (norm RMSE vs baseline):\n');
for j = 1:5
    if any(holdA(:, j))
        fprintf('  test %-12s : %.2f  (n=%d folds)\n', labels(j), ...
            mean(normA(holdA(:, j), j)), sum(holdA(:, j)));
    end
end

%% ---------- (C) new folds ----------
fprintf('\n===== (C) NEW SMOKE FOLDS (2brk evaluator, ga 25x30) =====\n');

fprintf('\n--- E1: ALLBPA=[1..5], NUMHOLD=1 (5 folds) ---\n');
E1 = runFolds(1:5, num2cell(1:5), a0_5, lb, ub, opts);

fprintf('\n--- E3: ALLBPA=[2,3,4,5], NUMHOLD=2 (6 folds; does dropping test 1 hurt?) ---\n');
E3 = runFolds(2:5, num2cell(nchoosek(2:5, 2)), a0_5, lb, ub, opts);

fprintf('\n--- E4: holdout [2,3,4] train [1,5] (Ben''s split) ---\n');
E4a = runFolds(1:5, {[2 3 4]}, a0_5, lb, ub, opts);
fprintf('\n--- E4rev: holdout [1,5] train [2,3,4] ---\n');
E4b = runFolds(1:5, {[1 5]}, a0_5, lb, ub, opts);

fprintf('\n--- E5: holdout [1,2,3,4] train [5] (train=1 extreme) ---\n');
E5 = runFolds(1:5, {[1 2 3 4]}, a0_5, lb, ub, opts);

fprintf('\n===== PER-TEST HELD-OUT NORM RMSE SUMMARY (best candidate per fold) =====\n');
fprintf('%-34s %8s %8s %8s %8s %8s\n', 'experiment \ test', labels);
sumE = {E1, E3, E4a, E4b, E5};
namesE = {'E1 leave1 [1-5]', 'E3 leave2 [2-5]', 'E4 Ben split', 'E4rev', 'E5 train=1'};
for e = 1:5
    foldsE = sumE{e};
    holdMask = false(numel(foldsE), 5);
    for k = 1:numel(foldsE), holdMask(k, foldsE{k}.holdout) = true; end
    for k = 1:numel(foldsE)
        rc = foldsE{k};
        [~, b] = min(rc.distance_all);
        pt = perTestNorm(rc.optParams_all(b, :), rc.holdout, a0_5);
        fprintf('%-34s', sprintf('%s h=[%s]', namesE{e}, num2str(rc.holdout)));
        row = nan(1, 5); row(rc.holdout) = pt;
        fprintf(' %8.3f', row);
        fprintf('\n');
    end
    fprintf('%-34s', sprintf('  ^ %s held-out MEAN', namesE{e}));
    for j = 1:5
        if any(holdMask(:, j))
            m = collectMean(foldsE, holdMask(:, j), j, a0_5);
            fprintf(' %8.3f', m);
        else
            fprintf(' %8s', '-');
        end
    end
    fprintf('\n\n');
end
fprintf(['Reference, legacy leave-2 (0817, mining yesterday): 48cm 0.29 | 46cm 0.15 | ', ...
    '47cm 0.49 | 40cm-tendon 0.21 | 41cm 0.17\n']);

outDir = fullfile(here, 'Dig_out');
if ~exist(outDir, 'dir'), mkdir(outDir); end
save(fullfile(outDir, 'Dig_allbpaNumHoldScan_20260908.mat'), 'prov', 'normA', 'holdA', ...
    'distA', 'E1', 'E3', 'E4a', 'E4b', 'E5', 'a0_5', 'labels');
fprintf('\nSaved %s\n', fullfile(outDir, 'Dig_allbpaNumHoldScan_20260908.mat'));
end

%% ================= helpers =================
function out = runFolds(allbpa, foldList, a0, lb, ub, opts)
    nf = numel(foldList);
    out = cell(nf, 1);
    for k = 1:nf
        holdoutIdx = foldList{k}(:)';
        trainIdx = setdiff(allbpa, holdoutIdx);
        fprintf('  fold: holdout [%s], train [%s] ... ', num2str(holdoutIdx), num2str(trainIdx));
        [x, fvals] = gamultiobj(@(X) min1(X, trainIdx, a0), 3, [], [], [], [], lb, ub, opts);
        valF = zeros(size(x, 1), 3);
        parfor i = 1:size(x, 1)
            valF(i, :) = min1(x(i, :), holdoutIdx, a0);
        end
        d = vecnorm(fvals - valF, 2, 2);
        out{k} = struct('holdout', holdoutIdx, 'train', trainIdx, 'optParams_all', x, ...
            'trainScores_all', fvals, 'validation_all', valF, 'distance_all', d);
        [~, b] = min(d);
        fprintf('done. best: Xi0=%.4f Xi1=%.3e Xi2=%.3e | val RMSE-ratio %.3f FVU-ratio %.3f\n', ...
            x(b, 1)/100, 10^x(b, 2), 10^x(b, 3), valF(b, 1), valF(b, 2));
    end
end

function ff = min1(x, trainIdx, kompare)
    Xi0 = x(1)/100; Xi1 = 10^x(2); Xi2 = 10^x(3);
    try
        f_all = minimizeFlxPin2brk(Xi0, Xi1, Xi2, trainIdx, true);
        ff = mean(f_all(trainIdx, :) ./ kompare(trainIdx, :), 1, 'omitnan');
    catch
        ff = [Inf, Inf, Inf];
    end
end

function m = collectMean(foldsE, mask, j, a0)
    vals = [];
    for k = find(mask(:))'
        rc = foldsE{k};
        [~, b] = min(rc.distance_all);
        f_j = minimizeFlxPin2brk(rc.optParams_all(b, 1)/100, 10^rc.optParams_all(b, 2), ...
            10^rc.optParams_all(b, 3), j, true);
        vals(end+1) = f_j(j, 1) / a0(j, 1); %#ok<AGROW>
    end
    m = mean(vals);
end

function r = perTestNorm(xrow, holdoutIdx, a0)
    Xi0 = xrow(1)/100; Xi1 = 10^xrow(2); Xi2 = 10^xrow(3);
    r = zeros(numel(holdoutIdx), 1);
    for t = 1:numel(holdoutIdx)
        j = holdoutIdx(t);
        f_j = minimizeFlxPin2brk(Xi0, Xi1, Xi2, j, true);
        r(t) = f_j(j, 1) / a0(j, 1);
    end
end
