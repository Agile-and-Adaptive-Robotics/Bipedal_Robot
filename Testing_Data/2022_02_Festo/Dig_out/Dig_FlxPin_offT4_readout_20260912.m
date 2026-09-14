%% Dig_FlxPin_offT4_readout_20260912.m
% Readout of minimizeFlxPin10_results_20260911_2brkt_2trans_offT4.mat:
% structure check + Ben's preferred split (fold 1: train {3,4,5}, holdout {1,2}).
% NOTE: stored train/val scores are BASELINE-NORMALIZED ratios (min1 output).

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(here)));
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));

S = load('minimizeFlxPin10_results_20260911_2brkt_2trans_offT4.mat');
fprintf('mat vars: %s\n', strjoin(fieldnames(S), ', '));
fprintf('folds in results_cv: %d | ALLBPA = [%s] | numHold = %d\n', ...
    numel(S.results_cv), num2str(S.allBPA), S.numHold);
fprintf('filtered_results: %d rows | pick k1/k2/k3 = %.4g / %.4g / %.4g\n', ...
    size(S.filtered_results,1), S.k1, S.k2, S.k3);
fprintf('pick per-test GoF (RMSE/FVU/MaxRes):\n');
disp(array2table(S.f, 'VariableNames', {'RMSE','FVU','MaxResidual'}, 'RowNames', cellstr(S.labels')));

%% Fold 1 = Ben's preferred split
r1 = S.results_cv{1};
fprintf('\nfold-1 foldIdx (unique rows): [%s]  (train = setdiff([1..5],fold) = [3 4 5])\n', ...
    num2str(unique(r1.foldIdx, 'rows')));
x1 = r1.optParams_all;
xi1 = [x1(:,1)/100, 10.^x1(:,2), 10.^x1(:,3)];   %de-normalize [m, N/m, N/m]
[~, ord] = sort(r1.distance_all);
n1 = min(8, numel(ord));
fprintf('\nfold-1 (train 3,4,5 / holdout 1,2) top-%d by train-vs-validation distance:\n', n1);
fprintf('   Xi0(m)      Xi1       Xi2    | tRMSE tFVU tMax | vRMSE vFVU vMax | dist\n');
for i = ord(1:n1)'   %row vector: a column vector would run the loop ONCE with all indices
    fprintf('%9.5f %8.3e %8.3e | %5.3f %5.3f %5.2f | %5.3f %5.3f %5.2f | %5.3f\n', ...
        xi1(i,:), r1.trainScores_all(i,:), r1.validation_all(i,:), r1.distance_all(i));
end

% Raw per-test GoF at the fold-1 best-distance candidate
best = ord(1);
fprintf('\nRaw per-test GoF at fold-1 best-distance candidate (Xi0 %.5f, Xi1 %.3e, Xi2 %.3e):\n', xi1(best,:));
[fRaw, ~] = minimizeFlxPin(xi1(best,1), xi1(best,2), xi1(best,3));
disp(array2table(fRaw, 'VariableNames', {'RMSE','FVU','MaxResidual'}, 'RowNames', cellstr(S.labels)));
