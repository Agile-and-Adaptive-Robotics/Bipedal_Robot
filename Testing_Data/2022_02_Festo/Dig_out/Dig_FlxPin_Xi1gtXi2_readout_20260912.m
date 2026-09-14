%% Dig_FlxPin_Xi1gtXi2_readout_20260912.m
% Side-by-side readout: constrained (Xi1 > Xi2) vs unconstrained offT4 campaigns.

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(here)));
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));

Su = load('minimizeFlxPin10_results_20260911_2brkt_2trans_offT4.mat');
Sc = load('minimizeFlxPin10_results_20260912_2brkt_2trans_offT4_Xi1gtXi2.mat');

fprintf('unconstrained pick: Xi0 %.2e m | Xi1 %.3e | Xi2 %.3e | %d rows, filter %d/%d\n', ...
    Su.k1, Su.k2, Su.k3, size(Su.filtered_results,1), size(Su.filtered_results,1), numel(Su.all_candidates));
fprintf('constrained   pick: Xi0 %.2e m | Xi1 %.3e | Xi2 %.3e | %d rows, filter %d/%d\n', ...
    Sc.k1, Sc.k2, Sc.k3, size(Sc.filtered_results,1), size(Sc.filtered_results,1), numel(Sc.all_candidates));

fprintf('\nper-test GoF, unconstrained pick vs constrained pick:\n');
fprintf('   test          unc RMSE/FVU/Max        con RMSE/FVU/Max\n');
for j = 1:5
    fprintf('%-11s %6.3f %6.4f %6.3f   |  %6.3f %6.4f %6.3f\n', Su.labels(j), ...
        Su.f(j,:), Sc.f(j,:));
end

%constraint sanity on the constrained pooled front
rc = Sc.filtered_results(:,5) ./ Sc.filtered_results(:,6);
fprintf('\nconstrained front Xi1/Xi2 ratio: min %.3f, median %.3f (must all be >= 1: %s)\n', ...
    min(rc), median(rc), string(all(rc >= 1 - 1e-9)));

%Xi ranges across each pooled front (Xi0 m, Xi1, Xi2 = cols 4,5,6)
fprintf('unc front Xi1 range [%.2e, %.2e], Xi2 range [%.2e, %.2e]\n', ...
    min(Su.filtered_results(:,5)), max(Su.filtered_results(:,5)), min(Su.filtered_results(:,6)), max(Su.filtered_results(:,6)));
fprintf('con front Xi1 range [%.2e, %.2e], Xi2 range [%.2e, %.2e]\n', ...
    min(Sc.filtered_results(:,5)), max(Sc.filtered_results(:,5)), min(Sc.filtered_results(:,6)), max(Sc.filtered_results(:,6)));

%% fold-1 (train 3,4,5 / holdout 1,2) top rows, constrained
r1 = Sc.results_cv{1};
x1 = [r1.optParams_all(:,1)/100, 10.^r1.optParams_all(:,2), 10.^r1.optParams_all(:,3)];
[~, ord] = sort(r1.distance_all);
n1 = min(6, numel(ord));
fprintf('\nconstrained fold-1 top-%d by distance (train {3,4,5} / holdout {1,2}):\n', n1);
fprintf('   Xi0(m)      Xi1       Xi2    Xi1/Xi2 | tRMSE tFVU tMax | vRMSE vFVU vMax | dist\n');
for i = ord(1:n1)'
    fprintf('%9.5f %8.3e %8.3e %7.2f | %5.3f %5.3f %5.2f | %5.3f %5.3f %5.2f | %5.3f\n', ...
        x1(i,:), x1(i,2)/x1(i,3), r1.trainScores_all(i,:), r1.validation_all(i,:), r1.distance_all(i));
end
best = ord(1);
fprintf('\nraw per-test GoF at constrained fold-1 best (Xi0 %.5f, Xi1 %.3e, Xi2 %.3e):\n', x1(best,:));
[fRaw, ~] = minimizeFlxPin(x1(best,1), x1(best,2), x1(best,3));
disp(array2table(fRaw, 'VariableNames', {'RMSE','FVU','MaxResidual'}, 'RowNames', cellstr(Sc.labels)));
