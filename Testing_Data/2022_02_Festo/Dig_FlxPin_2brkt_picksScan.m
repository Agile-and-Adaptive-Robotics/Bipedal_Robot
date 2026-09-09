function summary = Dig_FlxPin_2brkt_picksScan(srcMat, nPicks, modeOverride)
%DIG_FLXPIN_2BRKT_PICKSSCAN Score the WHOLE filtered Pareto front of a 2brk
%results .mat, then cross-predict a shortlist of picks in detail.
% Laptop-session script (2026-09-08), renamed into the Collect/Dig scheme same
% day after Ben's adjudication: the laptop's separate-file 1trans evaluator was
% deleted -- the transMode flag on minimizeFlxPin2brk is the vehicle of record.
% Still UNTESTED end-to-end (the laptop session never executed it); smoke on a
% small nPicks first.
%
% FAST PASS (all filtered candidates, ~2-3 min for ~265: one 5-test evaluator
%   call each): flxPin = mean RMSE ratio vs baseline (<1 beats rigid baseline).
%   Prints the plateau structure and writes the full table to CSV + .mat.
%   The front is nearly flat in validation distance, so "pick" is effectively a
%   choice of Xi1 - this pass shows whether flexor fit quality discriminates
%   at all along the front.
% DEEP PASS (picks 1..nPicks in val-distance order, ~1-2 min each):
%   bioFlx  : biomimetic flexor (minimizeFlx) mean RMSE
%   extPin  : pinned extensor pool {1,2,5,6,7,8} mean RMSE, pure substitution
%             (published Xi0/Xi3 + this pick's Xi1/Xi2)
%   rfXi0/rfXi3/extPinRF : coarse refit of (Xi0,Xi3) with Xi1/Xi2 LOCKED to the
%             pick (fminsearch 100 evals, scan-grade; use Dig_crossPredict for
%             an exact refit once a final pick is chosen)
%   extBio  : biomimetic extensor (52cm) mean RMSE, pure substitution
%
% Usage:
%   Dig_FlxPin_2brkt_picksScan('minimizeFlxPin10_results_20260907_2brkt_2trans.mat', 12)
%   Dig_FlxPin_2brkt_picksScan()                            % newest 2brkt file
%   Dig_FlxPin_2brkt_picksScan(srcMat, nPicks, '1trans')    % force evaluator variant
% Evaluator variant defaults to the TRANSMODE stored in the results file (if
% present; older files without it default to '2trans'). The canonical 09-07
% mat has no TRANSMODE field (pre-416f08f driver) -- provenance VERIFIED on
% easteregg2: it is the restored-2trans chain run (corrected Pbr2), so the
% default is correct for it.

here = fileparts(mfilename('fullpath'));
cd(here);
root = fileparts(fileparts(here));
addpath(fullfile(root, 'Code', 'Matlab', 'Functions'));
addpath(fullfile(root, 'Code', 'Matlab', 'Functions', 'ModernRobotics'));
addpath(fullfile(root, 'Code', 'Matlab', 'Robot_Data'));
set(0, 'DefaultFigureVisible', 'off');
outDir = fullfile(here, 'Dig_out');
if ~exist(outDir, 'dir'), mkdir(outDir); end

if nargin < 3, modeOverride = ''; end
if nargin < 2 || isempty(nPicks), nPicks = 8; end
if nargin < 1 || isempty(srcMat)
    d = dir(fullfile(here, 'minimizeFlxPin10_results_*_2brkt_*.mat'));
    d = d(~contains({d.name}, '_smoke'));
    assert(~isempty(d), 'no non-smoke 2brkt results found - pass srcMat explicitly');
    [~, ix] = max([d.datenum]);
    srcMat = fullfile(d(ix).folder, d(ix).name);
end
S = load(srcMat, 'filtered_results', 'xCols', 'results_sort_actual', 'TRANSMODE');
fprintf('=== Dig_FlxPin_2brkt_picksScan | %s ===\n', srcMat);

mode = '2trans';                       %default for pre-TRANSMODE files (verified for the 09-07 mat)
if isfield(S, 'TRANSMODE') && ~isempty(S.TRANSMODE), mode = S.TRANSMODE; end
if ~isempty(modeOverride), mode = modeOverride; end
switch lower(mode)
    case '2trans'
        evalFlx = @(a, b, c) minimizeFlxPin2brk(a, b, c, [], true);
    case '1trans'
        %Ben 2026-09-08: the transMode arg is the 1-trans vehicle (the laptop's
        %separate-file evaluator was deleted); proof: 1trans==2trans for these
        %y-symmetric K arrays, so this branch exists for completeness only.
        evalFlx = @(a, b, c) minimizeFlxPin2brk(a, b, c, [], true, '1trans');
    otherwise
        error('unknown transmode "%s"', mode);
end
fprintf('evaluator variant: %s | candidates: %d filtered, %d total\n', ...
    mode, size(S.filtered_results, 1), size(S.results_sort_actual, 1));

pool = S.filtered_results;
if isempty(pool)
    fprintf('filtered_results empty - scanning results_sort_actual instead (NO baseline-pass guarantee)\n');
    pool = S.results_sort_actual;
end
nAll = size(pool, 1);
nPicks = min(nPicks, nAll);

% Baselines (fresh, on ALL 5 tests, so ratios align regardless of the run's ALLBPA)
a0full = evalFlx(0, Inf, Inf);
baseFlx = minimizeFlx(0, Inf, Inf);
a0e = minimizeExtX3(0, Inf, Inf, 0);
be0 = minimizeExt(0, Inf, Inf, 0, 1);
POOL = [1 2 5 6 7 8];   % tests 3/4/9 excluded per Ben

%% ---------- FAST PASS: whole front ----------
fprintf('\n--- FAST PASS: flxPin RMSE ratio for all %d filtered candidates ---\n', nAll);
fastXi = pool(:, S.xCols);             % [Xi0 m, Xi1, Xi2] per candidate
fastDist = pool(:, end);
fastR = zeros(nAll, 1);
for p = 1:nAll
    fAll = evalFlx(fastXi(p, 1), fastXi(p, 2), fastXi(p, 3));
    fastR(p) = mean(fAll(:, 1) ./ a0full(:, 1));
end
[rSort, rOrd] = sort(fastR);
rBest = rSort(1);
plateau = fastR <= rBest * 1.01;       %within 1% of the front's best flexor fit
fprintf('best flxPin ratio on front: %.4f (candidate %d) | worst: %.4f\n', ...
    rBest, rOrd(1), rSort(end));
fprintf('candidates within 1%% of best: %d | Xi1 range there: %.3e .. %.3e\n', ...
    sum(plateau), min(fastXi(plateau, 2)), max(fastXi(plateau, 2)));
fprintf('Xi1 range over whole front: %.3e .. %.3e | Xi2: %.3e .. %.3e | Xi0: %.4f .. %.4f\n', ...
    min(fastXi(:, 2)), max(fastXi(:, 2)), min(fastXi(:, 3)), max(fastXi(:, 3)), ...
    min(fastXi(:, 1)), max(fastXi(:, 1)));
fprintf('\n top 25 candidates BY FLEXOR FIT (re-ranked; "pick" col = val-distance rank):\n');
fprintf('pick  flxPin     Xi0(m)        Xi1        Xi2     dist\n');
for q = 1:min(25, nAll)
    p = rOrd(q);
    fprintf('%4d  %6.4f  %8.4f  %9.3e  %8.3e  %6.3f\n', ...
        p, fastR(p), fastXi(p, 1), fastXi(p, 2), fastXi(p, 3), fastDist(p));
end
[~, srcName] = fileparts(srcMat);
csvFile = fullfile(outDir, sprintf('Dig_FlxPin_2brkt_picksScan_%s_fast.csv', srcName));
outTbl = [(1:nAll)', (1:nAll)', fastXi, fastR, fastDist];
fid = fopen(csvFile, 'w');
fprintf(fid, 'pick,valDistRank,Xi0_m,Xi1,Xi2,flxPinRMSEratio,valDist\n');
fprintf(fid, '%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g\n', outTbl');
fclose(fid);
fprintf('full front table written to %s\n', csvFile);

%% ---------- DEEP PASS: shortlist cross-prediction ----------
E = load(fullfile(here, 'minimizeExtPin10_results_20260819_2transforms_Z2.mat'), ...
    'filtered_results', 'xCols');
ep = E.filtered_results(1, E.xCols);
fprintf('\n--- DEEP PASS: picks 1..%d (val-distance order) | published extensor ref: ', nPicks);
fprintf('Xi0=%.4f Xi1=%.3e Xi2=%.3e Xi3=%.4f ---\n', ep(1), ep(2), ep(3), ep(4));

rows = zeros(nPicks, 10);
fprintf('\npick     Xi0(m)        Xi1        Xi2   flxPin  bioFlx  extPin   rfXi0    rfXi3  extPinRF  extBio\n');
for p = 1:nPicks
    g = pool(p, S.xCols);
    bF = minimizeFlx(g(1), g(2), g(3));
    pureE = minimizeExtX3(ep(1), g(2), g(3), ep(4));
    rPure = mean(pureE(POOL, 1));
    [rf0, rf3, refitE] = refitCoarse(ep, g, a0e, POOL);
    bioE = minimizeExt(ep(1), g(2), g(3), ep(4), 1);
    rows(p, :) = [g, fastR(p), mean(bF(:, 1)), rPure, rf0, rf3, mean(refitE(POOL, 1)), mean(bioE(:, 1))];
    fprintf('%3d  %8.4f  %9.3e  %8.3e  %6.4f  %6.2f  %6.2f  %7.4f  %7.4f  %7.2f  %6.2f\n', ...
        p, rows(p, 1), rows(p, 2), rows(p, 3), rows(p, 4), rows(p, 5), rows(p, 6), ...
        rows(p, 7), rows(p, 8), rows(p, 9), rows(p, 10));
end
fprintf(['\nflxPin <1 beats baseline; bioFlx/extPin/extBio are raw mean RMSE (N*m) - ', ...
    'baselines: bioFlx %.2f | extPin pool %.2f | extBio %.2f\n'], ...
    mean(baseFlx(:, 1)), mean(a0e(POOL, 1)), mean(be0(:, 1)));

outFile = fullfile(outDir, sprintf('Dig_FlxPin_2brkt_picksScan_%s.mat', srcName));
save(outFile, 'rows', 'fastXi', 'fastR', 'fastDist', 'csvFile', ...
    'srcMat', 'mode', 'nPicks', 'POOL', 'ep', 'baseFlx', 'a0e', 'be0', 'a0full');
fprintf('Saved %s\n', outFile);
summary = struct('rows', rows, 'fastR', fastR, 'fastXi', fastXi, 'src', srcMat, 'mode', mode);
end

function [Xi0out, Xi3out, fitAll] = refitCoarse(ep, g, a0e, POOL) %#ok<INUSD>
    fun = @(z) refitObj(z, g, a0e, POOL);
    z0 = [ep(1), max(min(ep(4), 1), 0)];
    opts = optimset('Display', 'off', 'MaxFunEvals', 100, 'MaxIter', 60);
    [z, ~] = fminsearch(fun, z0, opts);
    Xi0out = z(1); Xi3out = z(2);
    fitAll = minimizeExtX3(Xi0out, g(2), g(3), Xi3out);
end

function f = refitObj(z, g, a0e, POOL)
    Xi0 = z(1); Xi3 = z(2);
    if Xi0 < -0.02 || Xi0 > 0 || Xi3 < 0 || Xi3 > 1
        f = 1e6; return
    end
    try
        a = minimizeExtX3(Xi0, g(2), g(3), Xi3);
        f = mean(a(POOL, :) ./ a0e(POOL, :), 'all');
    catch
        f = 1e6;
    end
end
