% finish_diag_lock_runs.m — 2026-09-15
% Post-hoc completion of the two DIAG lock runs whose GAs succeeded but whose
% compile step hit the length() bug (fixed in minimizeFlxPin10mmX3_2brkt.m).
% Replicates the driver's compile -> sort -> de-normalize -> filter -> pick.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
setenv('FLXPX3_TANGENCY', 'DIAG');   %match the env the folds ran under

mats = {'minimizeFlxPin10mmX3_2brkt_results_20260915_L107_DIAG.mat', ...
        'minimizeFlxPin10mmX3_2brkt_results_20260915_L1trans_DIAG.mat'};
labels = ["48cm", "46cm", "47cm", "40cm-tendon", "42cm"];

for mi = 1:2
    file = mats{mi};
    fprintf('\n############ %s ############\n', file);
    S = load(file);
    results_cv = S.results_cv; numBPA = S.numBPA; allBPA = S.allBPA;
    numHold = size(results_cv{1}.foldIdx, 2);

    %% compile (fixed ind)
    all_candidates = [];
    for i = 1:numBPA
        fold = results_cv{i}.foldIdx;
        x2 = results_cv{i}.optParams_all;
        train = results_cv{i}.trainScores_all;
        val = results_cv{i}.validation_all;
        dist = results_cv{i}.distance_all;
        ind = (1:size(x2,1)).';
        rows = [ind, fold, x2, train, val, dist];
        all_candidates = [all_candidates; rows]; %#ok<AGROW>
    end
    rankCol = 1; holdCols = 2:(1+numHold);
    xCols = numel(rankCol) + numel(holdCols) + (1:4);
    trainCols = xCols(end) + (1:3); valCols = trainCols(end) + (1:3); distCol = valCols(end) + 1;
    results = all_candidates;
    results_sort = sortrows(results, [distCol valCols(1:2) trainCols(1:2)]);
    x_actual = [results_sort(:,xCols(1))/100, 10.^results_sort(:,xCols(2)), ...
                10.^results_sort(:,xCols(3)), 10.^results_sort(:,xCols(4))];
    results_sort_actual = [results_sort(:,rankCol), results_sort(:,holdCols), x_actual, ...
                           results_sort(:,trainCols), results_sort(:,valCols), results_sort(:,distCol)];
    fprintf('Compiled front: %d rows. Sorted (Xi0 cm, Xi1, Xi2, Xi3, trainRMSEratio, valRMSEratio, dist):\n', ...
        size(results_sort_actual,1));
    for ii = 1:size(results_sort_actual,1)
        fprintf('  %2d | Xi0 %6.4f  Xi1 %9.4g  Xi2 %9.4g  Xi3 %9.4g | train %.3f  val %.3f  dist %.4f\n', ...
            ii, results_sort_actual(ii,4), results_sort_actual(ii,5), results_sort_actual(ii,6), ...
            results_sort_actual(ii,7), results_sort_actual(ii,8), results_sort_actual(ii,11), results_sort_actual(ii,14));
    end

    %% baseline + filter (pass on BPAs 2-5, as in the driver)
    [a0, ~] = minimizeFlxPinX3_2brkt(0, Inf, Inf, Inf);
    baselineScores = a0;
    N = size(results_sort_actual, 1);
    keep = false(N,1);
    for ii = 1:N
        Xi0 = results_sort_actual(ii,xCols(1)); Xi1 = results_sort_actual(ii,xCols(2));
        Xi2 = results_sort_actual(ii,xCols(3)); Xi3 = results_sort_actual(ii,xCols(4));
        f_all = minimizeFlxPinX3_2brkt(Xi0, Xi1, Xi2, Xi3);
        keep(ii) = all(f_all(2,:) <= baselineScores(2,:)) && all(f_all(3,:) <= baselineScores(3,:)) ...
                && all(f_all(4,:) <= baselineScores(4,:)) && all(f_all(5,:) <= baselineScores(5,:));
    end
    filtered_results = results_sort_actual(keep, :);
    fprintf('Filtered %d -> %d candidates.\n', N, sum(keep));

    %% pick (driver uses pick = 2; fall back to 1 if only one row survived)
    if size(filtered_results,1) >= 2, pick = 2; else, pick = 1; end
    sol_actual = filtered_results(pick, xCols);
    k1 = sol_actual(1); k2 = sol_actual(2); k3 = sol_actual(3); k4 = sol_actual(4);
    [f, bpa] = minimizeFlxPinX3_2brkt(k1, k2, k3, k4);
    disp(array2table([k1, k2, k3, k4], 'VariableNames', {'X0_cm','X1','X2','X3'}));
    disp(array2table(f, 'VariableNames', {'RMSE','FVU','MaxResidual'}, 'RowNames', cellstr(labels)));
    fprintf('Mean optimized: RMSE %.4f, FVU %.4f, Max. Residual %.4f (pick %d of %d)\n\n', ...
        mean(f(:,1)), mean(f(:,2)), mean(f(:,3)), pick, size(filtered_results,1));

    %% save back, full workspace
    FINISHED = 'post-hoc compile+filter+pick 2026-09-15 (fixed ind); ENFORCE legs cancelled by Ben';
    save(file);
    fprintf('SAVED FULL WORKSPACE: %s\n', file);
end
fprintf('\nALL DONE\n');
