function summary = crossPredictFlx(srcMat, pick)
%CROSSPREDICTFLX Cross-prediction check for a flexor pinned-knee solution.
% Given flexor (Xi0, Xi1, Xi2), evaluate:
%   (a) biomimetic flexor   : minimizeFlx   vs its baseline
%   (b) pinned extensor     : minimizeExtX3 (9 tests) - reference (published
%       extensor solution), pure substitution (flexor Xi1/Xi2 + published
%       Xi0/Xi3), and constrained refit (Xi0/Xi3 re-solved with Xi1/Xi2 locked)
%   (c) biomimetic extensor : minimizeExt (52cm) - same three treatments
% Prints and saves tables; returns a summary struct.
%
%   summary = crossPredictFlx()                    % latest 2brk results, else legacy
%   summary = crossPredictFlx(srcMat, pick)        % explicit results .mat + candidate
%
% srcMat: results .mat holding filtered_results + xCols (3-col flexor solution).
% pick:   which filtered candidate to use (default 1).

if nargin < 2, pick = 1; end
if nargin < 1 || isempty(srcMat)
    d = dir([fileparts(mfilename('fullpath')), filesep, 'minimizeFlxPin10_2brk_results_*.mat']);
    if ~isempty(d)
        [~,ix] = max([d.datenum]); srcMat = fullfile(d(ix).folder, d(ix).name);
    else
        srcMat = 'minimizeFlxPin10_results_20260730_2transforms_Z2.mat';
    end
end
fprintf('=== crossPredictFlx | source: %s | pick=%d ===\n', srcMat, pick);
S = load(srcMat, 'filtered_results', 'xCols');
g = S.filtered_results(pick, S.xCols);   % [Xi0 m, Xi1 N/m, Xi2 N/m]
fprintf('Flexor solution: Xi0=%.4f m, Xi1=%.3e N/m, Xi2=%.3e N/m\n', g(1), g(2), g(3));

%% (a) Biomimetic flexor
baseFlx = minimizeFlx(0, Inf, Inf);            %metrics for the 10mm case (h{1})
predFlx = minimizeFlx(g(1), g(2), g(3));
fprintf('\n--- Biomimetic flexor (minimizeFlx, 10mm case) ---\n');
disp(array2table([baseFlx; predFlx], 'VariableNames', {'RMSE','FVU','MaxResidual'}, ...
    'RowNames', {'baseline','flexor sol'}));

%% Published extensor solution (reference)
E = load('minimizeExtPin10_results_20260819_2transforms_Z2.mat', 'filtered_results', 'xCols');
ep = E.filtered_results(pick, E.xCols);  %[Xi0 m, Xi1, Xi2, Xi3]
fprintf('\nPublished extensor solution: Xi0=%.4f, Xi1=%.3e, Xi2=%.3e, Xi3=%.4f\n', ep(1), ep(2), ep(3), ep(4));

%% (b) Pinned extensor, all 9 tests
labels9 = ["40cm","40cm-tendon","42cm","42cm-tendon","43cm","43cm-tendon","46cm","47cm","48cm"];
POOL = [1 2 5 6 7 8];   %tests 3,4,9 excluded per Ben (42cm, 42cm-tendon, 48cm)
a0 = minimizeExtX3(0, Inf, Inf, 0);
refE  = minimizeExtX3(ep(1), ep(2), ep(3), ep(4));
pureE = minimizeExtX3(ep(1), g(2),   g(3),   ep(4));
[rfXi0, rfXi3, refitE] = refitExt(@minimizeExtX3, ep, g, a0, POOL, 4);
fprintf('\nRefit (Xi1/Xi2 locked to flexor): Xi0=%.4f m, Xi3=%.4f\n', rfXi0, rfXi3);
fprintf('\n--- Pinned extensor (minimizeExtX3, 9 tests) ---\n');
disp(array2table([a0, refE, pureE, refitE], 'VariableNames', ...
    {'RMSE_base','FVU_base','MaxR_base','RMSE_pub','FVU_pub','MaxR_pub', ...
     'RMSE_pure','FVU_pure','MaxR_pure','RMSE_refit','FVU_refit','MaxR_refit'}, ...
    'RowNames', cellstr(labels9')));
fprintf('Pool mean RMSE: pub %.3f | pure %.3f | refit %.3f\n', ...
    mean(refE(POOL,1)), mean(pureE(POOL,1)), mean(refitE(POOL,1)));

%% (c) Biomimetic extensor (52cm)
be0   = minimizeExt(0, Inf, Inf, 0, 1);
refB  = minimizeExt(ep(1), ep(2), ep(3), ep(4), 1);
pureB = minimizeExt(ep(1), g(2),   g(3),   ep(4), 1);
[~, ~, refitB] = refitExt(@minimizeExt, ep, g, be0, 1, 4);
fprintf('\n--- Biomimetic extensor (minimizeExt, 52cm) ---\n');
disp(array2table([be0, refB, pureB, refitB], 'VariableNames', ...
    {'RMSE_base','FVU_base','MaxR_base','RMSE_pub','FVU_pub','MaxR_pub', ...
     'RMSE_pure','FVU_pure','MaxR_pure','RMSE_refit','FVU_refit','MaxR_refit'}, ...
    'RowNames', {'52cm'}));

%% Package + save
summary.src = srcMat; summary.pick = pick; summary.flexor = g; summary.extPub = ep;
summary.flxBio = [baseFlx; predFlx]; summary.extPin = {a0, refE, pureE, refitE};
summary.extBio = [be0, refB, pureB, refitB]; summary.pool = POOL;
summary.refit = [rfXi0, rfXi3];
stamp = char(string(datetime('now'),'yyyyMMdd'));
outFile = sprintf('crossPredict_%s.mat', stamp);
save(outFile, '-struct', 'summary');
fprintf('\nSaved %s\n', outFile);
end

function [Xi0out, Xi3out, fitAll] = refitExt(evaluator, ep, g, a0, idxSet, ncol)
%Re-solve (Xi0, Xi3) with Xi1/Xi2 LOCKED to the flexor values, minimizing the
%baseline-normalized mean metric over idxSet. Start from the published solution.
    pool = idxSet(isfinite(idxSet)); %#ok<NASGU>
    fun = @(z) refitObj(z, evaluator, g, a0, idxSet);
    z0 = [ep(1), max(min(ep(4), 1), 0)];
    opts = optimset('Display', 'iter', 'MaxFunEvals', 300, 'MaxIter', 200);
    [z, ~] = fminsearch(fun, z0, opts);
    Xi0out = z(1); Xi3out = z(2);
    fitAll = evalFull(evaluator, Xi0out, g(2), g(3), Xi3out, idxSet, ncol);
end

function f = refitObj(z, evaluator, g, a0, idxSet)
    Xi0 = z(1); Xi3 = z(2);
    if Xi0 < -0.02 || Xi0 > 0 || Xi3 < 0 || Xi3 > 1
        f = 1e6; return
    end
    try
        a = evalFull(evaluator, Xi0, g(2), g(3), Xi3, idxSet, 4);
        f = mean(a(idxSet,:) ./ a0(idxSet,:), 'all');
    catch
        f = 1e6;
    end
end

function a = evalFull(evaluator, Xi0, Xi1, Xi2, Xi3, idxSet, ncol) %#ok<INUSD>
    %evaluate on the union set needed (idxSet, plus all 9 for reporting by caller)
    if numel(idxSet) == 1 && idxSet == 1
        a = evaluator(Xi0, Xi1, Xi2, Xi3, 1);          %minimizeExt (biomimetic)
    else
        a = evaluator(Xi0, Xi1, Xi2, Xi3);             %minimizeExtX3 (all 9 tests)
    end
end
