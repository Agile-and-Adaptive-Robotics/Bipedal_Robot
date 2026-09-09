%Mine old "_all" CV mats: per-test held-out predictability patterns.
%RMSE+FVU primary; MaxResidual reported but never ranked (noisy tests).
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');

%% ================= FLEXOR PIN (0817_all, 10 leave-2-out folds) =================
S = load('minimizeFlxPin10_results_20260817_all.mat', 'results_cv');
cv = S.results_cv; nf = numel(cv);
labels5 = ["48cm","46cm","47cm","40cm-tendon","41cm"];
a0f = minimizeFlxPin(0, Inf, Inf);                 %5x3 baseline
normF = zeros(nf, 5);                              %norm RMSE per test, best cand per fold
holdF = false(nf, 5); distF = zeros(nf,1);
for k = 1:nf
    d = cv{k}.distance_all(:); [~,b] = min(d);
    xa = cv{k}.optParams_all(b,:);
    F = minimizeFlxPin(xa(1)/100, 10^xa(2), 10^xa(3));   %all 5 tests
    normF(k,:) = (F(:,1) ./ a0f(:,1))';
    holdF(k,:) = ismember(1:5, cv{k}.foldIdx(b,:));
    distF(k) = d(b);
end
fprintf('===== FLEXOR PIN: held-out normalized RMSE (best candidate per fold) =====\n');
fprintf('%-20s %7s %7s %7s %7s %7s\n', 'holdout', labels5);
for k = 1:nf
    fprintf('%-20s', strjoin(labels5(holdF(k,:)),'+'));
    fprintf(' %7.2f', normF(k,:)); fprintf('   (held-out:'); 
    fprintf(' %4.2f', normF(k, holdF(k,:))); fprintf(')\n');
end
fprintf('%-20s %7s\n', 'HELD-OUT MEAN', ' ');
for j = 1:5
    m = mean(normF(holdF(:,j), j), 'omitnan');
    fprintf('  test %-12s : held-out normRMSE %.2f  (n=%d folds)\n', labels5(j), m, sum(holdF(:,j)));
end

%% ================= EXTENSOR PIN (0809_all, 9 folds, numHold=4) =================
T = load('minimizeExtPin10_results_20260809_all.mat', 'results_sort_actual', 'validLabels');
R = T.results_sort_actual;
lab9 = string(T.validLabels);
a0e = minimizeExtX3(0, Inf, Inf, 0);               %9x3 baseline
combos = unique(R(:, 2:5), 'rows', 'stable');
ne = size(combos,1);
normE = zeros(ne, 9); holdE = false(ne, 9); distE = zeros(ne,1); xiE = zeros(ne,4);
for k = 1:ne
    same = all(R(:,2:5) == combos(k,:), 2);
    rowsE = R(same, :);
    [~,b] = min(rowsE(:,16));
    r = rowsE(b,:);
    Xi = [r(6)/100, 10^r(7), 10^r(8), r(9)];
    xiE(k,:) = Xi;
    F = minimizeExtX3(Xi(1), Xi(2), Xi(3), Xi(4));       %all 9 tests
    normE(k,:) = (F(:,1) ./ a0e(:,1))';
    holdE(k, combos(k, combos(k,:) > 0)) = true;
    distE(k) = r(16);
end
fprintf('\n===== EXTENSOR PIN: held-out normalized RMSE (best candidate per fold) =====\n');
fprintf('%-28s %6s %6s %6s %6s %6s %6s %6s %6s %6s\n', 'holdout(4 of 9)', lab9);
for k = 1:ne
    hs = strjoin(lab9(holdE(k,:)),'+');
    fprintf('%-28s', hs);
    fprintf(' %6.2f', normE(k,:)); fprintf('\n');
end
fprintf('%-28s (cols = tests 1..9)\n', ' ');
fprintf('Per-test held-out predictability (norm RMSE, averaged over folds where test was held out):\n');
for j = 1:9
    if any(holdE(:,j))
        m = mean(normE(holdE(:,j), j), 'omitnan');
        fprintf('  test %-12s : %.2f  (n=%d)\n', lab9(j), m, sum(holdE(:,j)));
    else
        fprintf('  test %-12s : never held out\n', lab9(j));
    end
end
fprintf('\nXi ranges across extensor folds: Xi1 [%.3g, %.3g]  Xi3 [%.2f, %.2f]\n', ...
    min(xiE(:,2)), max(xiE(:,2)), min(xiE(:,4)), max(xiE(:,4)));
