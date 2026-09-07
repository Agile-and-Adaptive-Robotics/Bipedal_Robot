function sweepResults = sweepExtX3(configs, XI1MIN, XI3MIN, MAXHOURS)
%SWEEPEXTX3 Hunt for pinned-extensor solutions with higher Xi1 and non-negligible Xi3.
% Runs the minimizeExt10mmX3 cross-validation pipeline over a small explicit
% config list, with:
%   - test pool trimmed to {1,2,5,6,7,8} (40cm, 40cm-tendon, 43cm, 43cm-tendon,
%     46cm, 47cm); tests 3 (42cm), 4 (42cm-tendon), 9 (48cm) are EXCLUDED
%   - Xi1/Xi2 bounds widened by a per-config factor x12factor of the flexor
%     solution (1 = locked to flexor values, as in minimizeExt10mmX3.m)
%   - budget guard: times the first fold, then truncates the fold list so the
%     config fits inside MAXHOURS
%   - solver from FLX2BRK_SOLVER (gamultiobj default | surrogateopt)
%   - post-filter: Xi1 >= XI1MIN (default: flexor Xi1) and Xi3 >= XI3MIN (0.05)
%   - confirmation of top survivors on all 9 pinned tests + biomimetic 52cm
%
%   sweepResults = sweepExtX3();              % default 2 configs
%   sweepResults = sweepExtX3(configs)        % struct array: allBPA,numHold,x12factor,pop,maxgen,surrvals
%   sweepResults = sweepExtX3(configs, XI1MIN, XI3MIN, MAXHOURS)

if nargin < 1 || isempty(configs)
    configs = struct( ...
        'allBPA', {[1 2 5 6 7 8], [1 2 5 6 7 8]}, ...
        'numHold', {3, 2}, ...
        'x12factor', {[1], [1]}, ...
        'pop', {50, 50}, 'maxgen', {150, 150}, 'surrvals', {4000, 4000});
end
if nargin < 2 || isempty(XI1MIN)
    S = load('minimizeFlxPin10_results_20260730_2transforms_Z2.mat', 'filtered_results', 'xCols');
    g = S.filtered_results(1, S.xCols);
    XI1MIN = g(2);
end
if nargin < 3 || isempty(XI3MIN), XI3MIN = 0.05; end
if nargin < 4 || isempty(MAXHOURS), MAXHOURS = 6; end

SOLVER = getenv('FLX2BRK_SOLVER');
if isempty(SOLVER), SOLVER = 'gamultiobj'; end
W = [0.5, 0.3, 0.2];   %scalarization weights for surrogateopt

labels9 = ["40cm","40cm-tendon","42cm","42cm-tendon","43cm","43cm-tendon","46cm","47cm","48cm"];
POOL = [1 2 5 6 7 8];

%Flexor reference solution (bounds anchor) and baseline (all 9 tests)
S = load('minimizeFlxPin10_results_20260730_2transforms_Z2.mat', 'filtered_results', 'xCols');
g = S.filtered_results(1, S.xCols);
a0 = minimizeExtX3(0, Inf, Inf, 0);

fprintf('=== sweepExtX3 | SOLVER=%s | XI1MIN=%.3e | XI3MIN=%.3f | MAXHOURS=%.1f ===\n', ...
    SOLVER, XI1MIN, XI3MIN, MAXHOURS);

sweepResults = struct('configs', {configs}, 'candidates', {[]}, 'confirmed', {[]});

for cIdx = 1:numel(configs)
    cfg = configs(cIdx);
    folds = nchoosek(cfg.allBPA, cfg.numHold);
    nFolds = size(folds, 1);

    %Bounds: Xi0 in cm; Xi1/Xi2 log10, widened by x12factor; Xi3 raw [xi3min,1]
    x3lb = 0;
    if isfield(cfg, 'xi3min'), x3lb = cfg.xi3min; end
    lb = [-0.02 * 100, log10(g(2)), log10(g(3)), x3lb];
    ub = [0 * 100, log10(cfg.x12factor * g(2)), log10(cfg.x12factor * g(3)), 1];

    fprintf('\n--- Config %d/%d: allBPA=[%s] numHold=%d x12factor=%g | %d folds ---\n', ...
        cIdx, numel(configs), num2str(cfg.allBPA), cfg.numHold, cfg.x12factor, nFolds);

    %Budget guard from a timing probe of the objective
    tProbe = tic;
    extmin1(middleOfBounds(lb, ub), cfg.allBPA, a0);
    perEval = toc(tProbe);
    nEvalPerFold = popEvals(cfg, SOLVER);
    etaFold = perEval * nEvalPerFold / 8;      %parfor across ~8-10 workers
    fprintf('Probe: %.2fs/eval, ~%.0fs/fold, ETA %.2f h for %d folds\n', ...
        perEval, etaFold, etaFold * nFolds / 3600, nFolds);
    folds = folds(1:min(nFolds, max(1, floor(MAXHOURS * 3600 / max(etaFold, 1)))), :);
    if size(folds, 1) < nFolds
        fprintf('WARNING: config truncated to %d/%d folds to fit MAXHOURS\n', size(folds,1), nFolds);
    end

    allCand = [];
    for k = 1:size(folds, 1)
        holdoutIdx = folds(k, :);
        trainIdx = setdiff(cfg.allBPA, holdoutIdx);
        fprintf('  Fold %d/%d: holdout [%s] (train [%s])\n', k, size(folds,1), ...
            num2str(holdoutIdx), num2str(trainIdx));

        opts = optimoptions('gamultiobj', ...
            'UseParallel', true, 'Display', 'off', ...
            'PopulationSize', cfg.pop, 'MaxGenerations', cfg.maxgen, ...
            'MutationFcn', {@mutationadaptfeasible}, ...
            'CrossoverFraction', 0.8, ...
            'CrossoverFcn', {@crossoverscattered}, ...
            'FunctionTolerance', 1e-3);

        tFold = tic;
        switch lower(SOLVER)
            case 'gamultiobj'
                [x, fvals] = gamultiobj(@(X) extmin1(X, trainIdx, a0), 4, [], [], [], [], lb, ub, opts);
            case 'surrogateopt'
                surrOpts = optimoptions('surrogateopt', 'UseParallel', true, ...
                    'Display', 'off', 'MaxFunctionEvaluations', cfg.surrvals);
                sol = surrogateopt(@(X) extmin1scalar(X, trainIdx, a0, W), lb, ub, surrOpts);
                %R2025a: first output IS the solution point (empty if all evals fail)
                assert(~isempty(sol), 'surrogateopt returned empty - all evaluations failed');
                x = reshape(sol, 1, []);
                fvals = extmin1(x, trainIdx, a0);
            otherwise
                error('Unknown FLX2BRK_SOLVER "%s"', SOLVER);
        end

        valF = zeros(size(x,1), 3);
        parfor i = 1:size(x,1)
            valF(i,:) = extmin1(x(i,:), holdoutIdx, a0);
        end
        dist = vecnorm(fvals - valF, 2, 2);
        ind = (1:size(x,1))';
        rows = [ind, repmat(holdoutIdx, size(x,1), numel(holdoutIdx)), x, fvals, valF, dist];
        allCand = [allCand; rows]; %#ok<AGROW>
        fprintf('    fold took %.1f min, %d pareto points\n', toc(tFold)/60, size(x,1));
    end

    %De-normalize: [ind, holdouts, Xi0_cm, log10 Xi1, log10 Xi2, Xi3, train3, val3, dist]
    nH = cfg.numHold;
    iX  = 1 + nH + (1:4);
    iT  = 1 + nH + 4 + (1:3);
    iV  = 1 + nH + 7 + (1:3);
    iD  = 1 + nH + 10 + 1;
    xAct = [allCand(:,iX(1))/100, 10.^allCand(:,iX(2)), 10.^allCand(:,iX(3)), allCand(:,iX(4))];
    cand = [allCand(:,1:1+nH), xAct, allCand(:,iT), allCand(:,iV), allCand(:,iD)];
    %sort: distance, then validation RMSE/FVU/MaxR (last 4 columns are val + dist)
    nCol = size(cand, 2);
    cand = sortrows(cand, [nCol, nCol-3, nCol-2, nCol-1]);

    %Post-filter: beat baseline on the pool, Xi1 >= XI1MIN, Xi3 >= XI3MIN
    keep = false(size(cand,1), 1);
    for ii = 1:size(cand,1)
        fAll = minimizeExtX3(cand(ii,2+nH), cand(ii,3+nH), cand(ii,4+nH), cand(ii,5+nH));
        keep(ii) = all(fAll(POOL,:) <= a0(POOL,:), 'all') ...
            && cand(ii,3+nH) >= XI1MIN && cand(ii,5+nH) >= XI3MIN;
    end
    survivors = cand(keep, :);
    fprintf('Config %d: %d candidates, %d survive (Xi1>=%.2e, Xi3>=%.2f, beats baseline)\n', ...
        cIdx, size(cand,1), size(survivors,1), XI1MIN, XI3MIN);
    sweepResults.candidates{cIdx} = survivors;

    %Confirm top survivors on all 9 pinned + biomimetic 52cm
    %conf columns: [Xi0, Xi1, Xi2, Xi3, trainRMSE, trainFVU, trainMaxR, valDist,
    %               poolRMSE_mean, poolRMSE_min, poolRMSE_max, bio52RMSE]
    nConf = min(5, size(survivors, 1));
    conf = nan(nConf, 12);
    confHoldouts = cell(nConf, 1);
    for ii = 1:nConf
        s = survivors(ii, :);
        confHoldouts{ii} = s(2:1+nH);
        f9 = minimizeExtX3(s(2+nH), s(3+nH), s(4+nH), s(5+nH));
        fb = minimizeExt(s(2+nH), s(3+nH), s(4+nH), s(5+nH), 1);
        conf(ii,:) = [s(2+nH:5+nH), s(6+nH:8+nH), s(end), ...
            mean(f9(POOL,1)), min(f9(POOL,1)), max(f9(POOL,1)), fb(1)];
    end
    sweepResults.confirmed{cIdx} = conf;
    if nConf > 0
        disp(array2table(conf, 'VariableNames', {'Xi0','Xi1','Xi2','Xi3', ...
            'trainRMSE','trainFVU','trainMaxR','valDist','poolRMSE_mean','poolRMSE_min','poolRMSE_max','bio52RMSE'}));
    end
end

stamp = char(string(datetime('now'),'yyyyMMdd'));
outFile = sprintf('sweepExtX3_results_%s.mat', stamp);
save(outFile, 'sweepResults', 'configs', 'XI1MIN', 'XI3MIN', 'POOL', 'SOLVER');
fprintf('\nSweep saved to %s\n', outFile);
end

function ff = extmin1(x, trainIdx, kompare)
    Xi0 = x(1)/100; Xi1 = 10^x(2); Xi2 = 10^x(3); Xi3 = x(4);
    try
        f_all = minimizeExtX3(Xi0, Xi1, Xi2, Xi3, trainIdx);
        fnorm = f_all(trainIdx,:) ./ kompare(trainIdx,:);
        ff = mean(fnorm, 1, 'omitnan');
    catch
        ff = [Inf, Inf, Inf];
    end
end

function fs = extmin1scalar(x, trainIdx, kompare, W)
    f3 = extmin1(x, trainIdx, kompare);
    if any(~isfinite(f3)), fs = Inf; else, fs = W * f3(:); end
end

function x0 = middleOfBounds(lb, ub)
    x0 = lb + 0.5 * (ub - lb);
    x0(2) = lb(2); x0(3) = lb(3);   %start Xi1/Xi2 at the flexor values
end

function n = popEvals(cfg, SOLVER)
    if strcmpi(SOLVER, 'surrogateopt')
        n = cfg.surrvals;
    else
        n = cfg.pop * cfg.maxgen;
    end
end
