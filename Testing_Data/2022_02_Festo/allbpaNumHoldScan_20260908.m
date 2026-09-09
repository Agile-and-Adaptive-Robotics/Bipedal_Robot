function allbpaNumHoldScan_20260908()
%ALLBPANUMHOLDSCAN_20260908 Ben's questions (2026-09-08), answered from data:
% STAGED 2026-09-08 by laptop zcode session (DESKTOP-5Q16KE9): UNTESTED, never
% executed. Rename/erase/adapt freely - see HANDOFF_laptop_20260908.md.
%
%   (0) PROVENANCE: the file "..._2brkt_2trans_smoke.mat" is actually the FULL
%       09-07 production CV (ALLBPA 1-5, NUMHOLD=2, POP 150, MAXGEN 600; no
%       TRANSMODE field inside). Which evaluator produced it - restored 2trans
%       or superseded 1trans? Re-evaluate one stored candidate with both and
%       compare to the stored training scores. ~1 min.
%   (A) Per-test held-out patterns from that stored full result (best candidate
%       per fold, leave-2-out, 2brk evaluator on all 5 tests).
%   (B) Angle-lag scan per test, several solutions/models: is one test's
%       measured curve angle-shifted (~+5 deg encoder concern)? Shift applied
%       to the MEASURED angles Aexp before interpolating the model torque.
%   (C) New smoke-scale gamultiobj folds with the CURRENT 2brk evaluator:
%       E1 ALLBPA=[1..5] NUMHOLD=1 (5 folds) - numHold contrast vs (A)
%       E3 ALLBPA=[2,3,4,5] NUMHOLD=2 (6 folds) - does dropping test 1 hurt?
%       E4 holdout [2,3,4] train [1,5] (Ben's split) + reverse [1,5]/[2,3,4]
%       E5 holdout [1,2,3,4] train [5] - train=1 extreme
% Runtime: (0)+(A)+(B) ~2 min; (C) 14 ga folds at 25x30, ~15-30 min on 6 workers.
% Self-locating paths. Results -> allbpaNumHoldScan_20260908.mat

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));   % 2022_02_Festo -> Testing_Data -> repo
cd(here);
addpath(fullfile(root, 'Code', 'Matlab', 'Functions'));
addpath(fullfile(root, 'Code', 'Matlab', 'Functions', 'ModernRobotics'));
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
fullMat = 'minimizeFlxPin10_results_20260907_2brkt_2trans_smoke.mat';  %MISNAMED: full run
Sm = load(fullMat);
a0_5 = minimizeFlxPin2brk(0, Inf, Inf, [], true);   %baseline on all 5 tests

%% ---------- (0) provenance ----------
fprintf('\n===== (0) PROVENANCE: which evaluator produced the stored full result? =====\n');
rc1 = Sm.results_cv{1};
[~, b0] = min(rc1.distance_all);
xr0 = rc1.optParams_all(b0, :);
tr1 = setdiff(1:5, unique(rc1.foldIdx(b0, :)));
fNow = evalRatio(xr0, tr1, @minimizeFlxPin2brk);
fOld = evalRatio(xr0, tr1, @minimizeFlxPin2brk_1trans);
fprintf('stored fold-1 best train scores : [%.4f %.4f %.4f]\n', rc1.trainScores_all(b0, :));
fprintf('recomputed with 2trans (current): [%.4f %.4f %.4f]\n', fNow);
fprintf('recomputed with 1trans (supersd): [%.4f %.4f %.4f]\n', fOld);
if max(abs(fNow - rc1.trainScores_all(b0, :))) < 1e-6
    fprintf('-> stored result MATCHES the current 2trans evaluator.\n');
    prov = '2trans';
elseif max(abs(fOld - rc1.trainScores_all(b0, :))) < 1e-6
    fprintf('-> stored result MATCHES the superseded 1trans evaluator (treat as 1trans data!).\n');
    prov = '1trans';
else
    fprintf('-> matches NEITHER exactly (different evaluator revision?) - inspect manually.\n');
    prov = 'unknown';
end

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

%% ---------- (B) angle-lag scan ----------
fprintf('\n===== (B) ANGLE-LAG SCAN: shift applied to measured angles =====\n');
fprintf('delta* = shift (deg) that best aligns measured torque with model; ratio = RMSE(delta*)/RMSE(0)\n');
sols = struct('name', {}, 'x', {});
bpas = {};

[fL0, bpaL0] = minimizeFlxPin(0, Inf, Inf);                          %#ok<NASGU>
sols(end+1) = struct('name', 'legacy baseline', 'x', [0 Inf Inf]); %#ok<SAGROW>
bpas{end+1} = bpaL0; %#ok<SAGROW>
[fL1, bpaL1] = minimizeFlxPin(0.01, 2e4, 0.8e4);                     %#ok<NASGU>
sols(end+1) = struct('name', 'legacy moderate', 'x', [0.01 2e4 0.8e4]); %#ok<SAGROW>
bpas{end+1} = bpaL1; %#ok<SAGROW>

bestLegacy = [];
if isfile('minimizeFlxPin10_results_20260817_all.mat')
    S17 = load('minimizeFlxPin10_results_20260817_all.mat', 'results_cv');
    dAll = []; for k = 1:numel(S17.results_cv)
        dAll = [dAll; S17.results_cv{k}.distance_all]; %#ok<AGROW>
    end
    [~, ib] = min(dAll);
    cum = 0; xr = [];
    for k = 1:numel(S17.results_cv)
        nk = numel(S17.results_cv{k}.distance_all);
        if ib <= cum + nk
            xr = S17.results_cv{k}.optParams_all(ib - cum, :);
            break;
        end
        cum = cum + nk;
    end
    if ~isempty(xr)
        [~, bpaLb] = minimizeFlxPin(xr(1)/100, 10^xr(2), 10^xr(3));
        sols(end+1) = struct('name', 'legacy best0817', 'x', xr(1:3)); %#ok<SAGROW>
        bpas{end+1} = bpaLb; %#ok<SAGROW>
        bestLegacy = xr;
        fprintf('legacy best0817 candidate: Xi0=%.4f m, Xi1=%.3e, Xi2=%.3e\n', xr(1)/100, 10^xr(2), 10^xr(3));
    end
end

[fB0, bpaB0] = minimizeFlxPin2brk(0, Inf, Inf, [], true);            %#ok<NASGU>
sols(end+1) = struct('name', '2brk baseline', 'x', [0 Inf Inf]); %#ok<SAGROW>
bpas{end+1} = bpaB0; %#ok<SAGROW>
kB = [];
kB = [Sm.k1, Sm.k2, Sm.k3];   %stored pick-1 of the full result
[~, bpaB1] = minimizeFlxPin2brk(kB(1), kB(2), kB(3), [], true);
sols(end+1) = struct('name', '2brk fullPick1', 'x', kB); %#ok<SAGROW>
bpas{end+1} = bpaB1; %#ok<SAGROW>
fprintf('2brk fullPick1: Xi0=%.4f m, Xi1=%.3e, Xi2=%.3e\n', kB(1), kB(2), kB(3));

testIdx = 1:5;
lagDelta = NaN(numel(sols), numel(testIdx));
lagRatio = NaN(numel(sols), numel(testIdx));
for s = 1:numel(sols)
    for j = testIdx
        [lagDelta(s, j), lagRatio(s, j)] = angleLagScan(bpas{s}(j));
    end
end
fprintf('\n%-16s', 'solution \ test');
fprintf(' %9s', labels);
fprintf('\n delta* (deg):');
for s = 1:numel(sols)
    fprintf('\n%-16s', sols(s).name);
    fprintf(' %9.2f', lagDelta(s, :));
end
fprintf('\n RMSE ratio   :');
for s = 1:numel(sols)
    fprintf('\n%-16s', sols(s).name);
    fprintf(' %9.3f', lagRatio(s, :));
end
fprintf('\n(If one test shows a consistent |delta*| ~ 5 deg across ALL solutions, suspect its encoder.)\n');

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

save('allbpaNumHoldScan_20260908.mat', 'prov', 'normA', 'holdA', 'distA', ...
    'sols', 'lagDelta', 'lagRatio', 'E1', 'E3', 'E4a', 'E4b', 'E5', ...
    'a0_5', 'bestLegacy', 'kB', 'labels');
fprintf('\nSaved allbpaNumHoldScan_20260908.mat\n');
end

%% ================= helpers =================
function r = evalRatio(xr, trainIdx, ev)
%Baseline-normalized training metrics for a stored candidate, given an evaluator.
    f = ev(xr(1)/100, 10^xr(2), 10^xr(3), [], true);
    a0 = ev(0, Inf, Inf, [], true);
    r = mean(f(trainIdx, :) ./ a0(trainIdx, :), 1, 'omitnan');
end

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

function [dStar, ratio] = angleLagScan(bpaJ)
%Shift MEASURED angles by delta so the model curve matches; find best delta.
    Ak = bpaJ.Ak; Mp = bpaJ.M_p(:, 3);
    valid = ~isnan(Mp) & ~isnan(Ak);
    [Aks, is] = sort(Ak(valid));
    Vs = Mp(valid);
    F = griddedInterpolant(Aks, Vs(is), 'linear', 'nearest');
    Aexp = bpaJ.Aexp(:); Mexp = bpaJ.Mexp(:);
    deltas = -10:0.25:10;
    rmse = zeros(size(deltas));
    for di = 1:numel(deltas)
        MpAt = F(Aexp + deltas(di));
        v = ~isnan(MpAt);
        rmse(di) = sqrt(mean((Mexp(v) - MpAt(v)).^2));
    end
    [rmBest, ib] = min(rmse);
    i0 = find(abs(deltas) < 1e-9);
    rm0 = rmse(i0);
    dStar = deltas(ib);
    ratio = rmBest / rm0;
end
