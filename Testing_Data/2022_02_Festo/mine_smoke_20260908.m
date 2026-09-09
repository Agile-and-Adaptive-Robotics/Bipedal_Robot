%MINE_SMOKE_20260908 read-only inspection of yesterday's 2brk smoke results and
% the 09-07 cross-prediction summary. No optimizer runs, no evaluator calls, no
% parpool - pure loading and printing (mining only, per Ben 2026-09-08).
here = fileparts(mfilename('fullpath'));
cd(here);

S = load('minimizeFlxPin10_results_20260907_2brkt_2trans_smoke.mat');
fprintf('fields in file: %s\n', strjoin(fieldnames(S)', ', '));
fprintf('settings: ALLBPA=[%s] NUMHOLD=%d POP=%d MAXGEN=%d PICK=%d\n', ...
    num2str(S.ALLBPA), S.NUMHOLD, S.POP, S.MAXGEN, S.PICK);
if isfield(S, 'SOLVER'), fprintf('SOLVER=%s\n', S.SOLVER); end
if isfield(S, 'TRANSMODE'), fprintf('TRANSMODE=%s\n', S.TRANSMODE); end
if isfield(S, 'USE_BRACKET2'), fprintf('USE_BRACKET2=%d\n', S.USE_BRACKET2); end
fprintf('candidates: %d total, %d passed the baseline filter\n', ...
    size(S.results_sort_actual, 1), size(S.filtered_results, 1));

n = min(12, size(S.filtered_results, 1));
fprintf('\n top filtered candidates (sorted by validation distance):\n');
fprintf('pick     Xi0(m)        Xi1        Xi2     dist\n');
for p = 1:n
    g = S.filtered_results(p, S.xCols);
    fprintf('%3d  %8.4f  %9.3e  %8.3e  %6.3f\n', p, g(1), g(2), g(3), S.filtered_results(p, end));
end

fprintf('\n pick-1 flexor metrics (f; rows = tests %s):\n', strjoin(string(S.labels(S.ALLBPA)), ', '));
disp(array2table(S.f, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}));
disp(array2table(S.a0, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}));

if isfile('crossPredict_20260907.mat')
    C = load('crossPredict_20260907.mat');
    fprintf('\n=== crossPredict 2026-09-07: src=%s pick=%d ===\n', C.src, C.pick);
    fprintf('flexor  Xi0=%.4f Xi1=%.3e Xi2=%.3e\n', C.flexor(1), C.flexor(2), C.flexor(3));
    fprintf('extPub  Xi0=%.4f Xi1=%.3e Xi2=%.3e Xi3=%.4f\n', C.extPub(1), C.extPub(2), C.extPub(3), C.extPub(4));
    fprintf('refit (Xi0, Xi3) = (%.4f, %.4f)\n', C.refit(1), C.refit(2));
    P = C.pool;
    fprintf('pinned-extensor pool mean RMSE: pub %.2f | pure %.2f | refit %.2f (N*m)\n', ...
        mean(C.extPin{2}(P, 1)), mean(C.extPin{3}(P, 1)), mean(C.extPin{4}(P, 1)));
    fprintf('biomimetic flexor mean RMSE: baseline %.2f -> pick %.2f\n', ...
        mean(C.flxBio(1, 1)), mean(C.flxBio(2, 1)));
    fprintf('biomimetic extensor RMSE: baseline %.2f | pub %.2f | pure %.2f | refit %.2f\n', ...
        C.extBio(1), C.extBio(4), C.extBio(7), C.extBio(10));
end
fprintf('\n(mine_smoke_20260908 done - nothing was executed, only loaded)\n');
