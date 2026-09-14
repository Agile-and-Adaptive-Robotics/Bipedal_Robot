% Dig_FlxPin_mixT_readout_20260913.m
% Readout of the two mixed-convention campaigns (2026-09-13):
%   run 1: minimizeFlxPin10_results_20260913_2brkt_mixT_2folds.mat      (no T5 shift)
%   run 2: minimizeFlxPin10_results_20260912..._mixT_T5y5mm_2folds.mat  (T5 +5mm bracket-y)
% Plus fixed-Xi isolation: evaluate each run's pick under BOTH FLX_T5_YMM settings.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
cd(root); addpath(root); addpath(here);
addpath(genpath(fullfile(fileparts(fileparts(root)), 'Code', 'Matlab')));

mats = {'minimizeFlxPin10_results_20260913_2brkt_mixT_2folds.mat', 'run 1 (no T5 shift)'; ...
        'minimizeFlxPin10_results_20260913_2brkt_mixT_T5y5mm_2folds.mat', 'run 2 (T5 +5mm bracket-y)'};
picks = cell(1,2);
for m = 1:2
    S = load(fullfile(root, mats{m,1}));
    fprintf('\n================ %s =================\n', mats{m,2});
    % folds (only the folds actually run: list = [1 5; 3 4] -> 2)
    hv = [1 5; 3 4];
    for i = 1:2
        d = S.results_cv{i}.distance_all;
        [dmin, ib] = min(d); %#ok<ASGLU>
        x = S.results_cv{i}.optParams_all(ib,:);
        tr = mean(S.results_cv{i}.trainScores_all(ib,:), 1);
        va = mean(S.results_cv{i}.validation_all(ib,:), 1);
        fprintf('fold %d (holdout [%d %d]): best-dist Xi0=%.4f mm  Xi1=%.4g  Xi2=%.4g  ratio=%.3f\n', ...
            i, hv(i,1), hv(i,2), x(1)*10, 10^x(2), 10^x(3), 10^(x(2)-x(3)));
        fprintf('   tRMSE-ratio %.3f  vRMSE-ratio %.3f\n', tr(1), va(1));
    end
    % pooled front
    X = S.filtered_results(:, S.xCols);
    fprintf('xCols = [%s]; col2 range [%.3f %.3f], col3 range [%.3f %.3f]\n', ...
        num2str(S.xCols), min(X(:,2)), max(X(:,2)), min(X(:,3)), max(X(:,3)));
    Xi1 = 10.^X(:,2); Xi2 = 10.^X(:,3);
    rat = Xi1./Xi2;
    fprintf('pooled front: %d rows, median Xi1/Xi2 = %.3f, frac ratio<1.01 = %.2f\n', ...
        size(S.filtered_results,1), median(rat), mean(rat < 1.01));
    % pick
    picks{m} = S.sol_actual;
    fprintf('pick sol_actual: Xi0=%.4f mm  Xi1=%.5g  Xi2=%.5g  ratio=%.4f\n', ...
        S.sol_actual(1)*1000, S.sol_actual(2), S.sol_actual(3), S.sol_actual(2)/S.sol_actual(3));
    fprintf('pick per-test      48cm    46cm    47cm    40cm-t  41cm\n');
    fprintf('pick per-test RMSE: %s\n', num2str(S.f(:,1)', '%.3f  '));
    fprintf('pick per-test FVU : %s\n', num2str(S.f(:,2)', '%.3f  '));
end

%% fixed-Xi shift isolation
fprintf('\n================ T5 shift at FIXED Xi =================\n');
lbl = {'48cm','46cm','47cm','40cm-t','41cm'};
for m = 1:2
    g = picks{m};
    fprintf('%s pick (Xi0=%.4fmm Xi1=%.4g Xi2=%.4g):\n', mats{m,2}, g(1)*10, g(2), g(3));
    for e = {'0','5'}
        setenv('FLX_T5_YMM', e{1});
        [f, ~] = minimizeFlxPin(g(1), g(2), g(3));
        fprintf('   FLX_T5_YMM=%s RMSE: %s\n', e{1}, num2str(f(:,1)', '%6.3f'));
    end
end
setenv('FLX_T5_YMM', '0');
fprintf('\nlabels: %s %s %s %s %s\n', lbl{:});
disp('READOUT DONE');
