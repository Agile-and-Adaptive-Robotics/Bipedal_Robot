% finish_diag_lock_runs_pass2.m — 2026-09-15
% Repair the validation columns of the two DIAG lock mats: the driver's valF
% parfor ran on pool workers that never saw the client's FLXPX3_TANGENCY=DIAG
% (workers keep the env from pool spawn), so held-out evals were ENFORCE-rejected
% to Inf. Recompute validation CLIENT-side (DIAG), re-sort, re-filter, re-pick.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
setenv('FLXPX3_TANGENCY', 'DIAG');   %client-side; serial calls see this

mats = {'minimizeFlxPin10mmX3_2brkt_results_20260915_L107_DIAG.mat', ...
        'minimizeFlxPin10mmX3_2brkt_results_20260915_L1trans_DIAG.mat'};
labels = ["48cm", "46cm", "47cm", "40cm-tendon", "42cm"];

for mi = 1:2
    file = mats{mi};
    fprintf('\n############ %s ############\n', file);
    S = load(file);
    results_cv = S.results_cv; numBPA = S.numBPA; allBPA = S.allBPA;
    numHold = size(results_cv{1}.foldIdx, 2);

    % compile (fixed ind)
    all_candidates = [];
    for i = 1:numBPA
        fold = results_cv{i}.foldIdx;
        x2 = results_cv{i}.optParams_all;
        train = results_cv{i}.trainScores_all;
        dist = results_cv{i}.distance_all;   %worker artifact, replaced below
        ind = (1:size(x2,1)).';
        rows = [ind, fold, x2, train, nan(size(train)), dist]; %#ok<AGROW>
        all_candidates = [all_candidates; rows]; %#ok<AGROW>
    end
    holdCols = 2:(1+numHold); xCols = (1+numHold) + (1:4);
    trainCols = xCols(end) + (1:3); valCols = trainCols(end) + (1:3); distCol = valCols(end) + 1;

    % client-side full evaluation per compiled candidate: train already stored;
    % validation = the candidate's own holdout rows of the full 5-test evaluation
    [a0, ~] = minimizeFlxPinX3_2brkt(0, Inf, Inf, Inf);
    N = size(all_candidates, 1);
    fvals = all_candidates(:, trainCols); vals = nan(N, 3); dists = nan(N, 1);
    fAll = cell(N, 1);
    for ii = 1:N
        Xi0 = all_candidates(ii,xCols(1))/100; Xi1 = 10.^all_candidates(ii,xCols(2));
        Xi2 = 10.^all_candidates(ii,xCols(3)); Xi3 = 10.^all_candidates(ii,xCols(4));
        f_all = minimizeFlxPinX3_2brkt(Xi0, Xi1, Xi2, Xi3);
        fAll{ii} = f_all;
        holdIdx = all_candidates(ii, holdCols);
        vals(ii,:) = mean(f_all(holdIdx,:) ./ a0(holdIdx,:), 1, 'omitnan');
        dists(ii) = vecnorm(fvals(ii,:) - vals(ii,:), 2);
    end
    all_candidates(:, valCols) = vals; all_candidates(:, distCol) = dists;

    results_sort = sortrows(all_candidates, [distCol valCols(1) trainCols(1)]);
    x_actual = [results_sort(:,xCols(1))/100, 10.^results_sort(:,xCols(2)), ...
                10.^results_sort(:,xCols(3)), 10.^results_sort(:,xCols(4))];
    results_sort_actual = [results_sort(:,1), results_sort(:,holdCols), x_actual, ...
                           results_sort(:,trainCols), results_sort(:,valCols), results_sort(:,distCol)];
    fprintf('Sorted candidates (Xi0 cm | Xi1 Xi2 Xi3 | trainR valR dist):\n');
    for ii = 1:size(results_sort_actual,1)
        fprintf('  %2d | %6.4f  %9.4g %9.4g %9.4g | %.3f %.3f %.4f\n', ii, ...
            results_sort_actual(ii,4), results_sort_actual(ii,5), results_sort_actual(ii,6), ...
            results_sort_actual(ii,7), results_sort_actual(ii,8), results_sort_actual(ii,11), results_sort_actual(ii,14));
    end

    % filter on BPAs 2-5 vs baseline (driver semantics)
    keep = false(N,1);
    for ii = 1:N
        keep(ii) = all(fAll{ii}(2,:) <= a0(2,:)) && all(fAll{ii}(3,:) <= a0(3,:)) ...
                && all(fAll{ii}(4,:) <= a0(4,:)) && all(fAll{ii}(5,:) <= a0(5,:));
    end
    filtered_results = results_sort_actual(keep, :);
    fprintf('Filtered %d -> %d candidates.\n', N, sum(keep));

    if size(filtered_results,1) >= 2, pick = 2; else, pick = 1; end
    sol_actual = filtered_results(pick, xCols);
    k1 = sol_actual(1); k2 = sol_actual(2); k3 = sol_actual(3); k4 = sol_actual(4);
    [f, bpa] = minimizeFlxPinX3_2brkt(k1, k2, k3, k4);
    disp(array2table([k1, k2, k3, k4], 'VariableNames', {'X0_cm','X1','X2','X3'}));
    disp(array2table(f, 'VariableNames', {'RMSE','FVU','MaxResidual'}, 'RowNames', cellstr(labels)));
    fprintf('Mean optimized: RMSE %.4f, FVU %.4f, Max. Residual %.4f (pick %d of %d)\n', ...
        mean(f(:,1)), mean(f(:,2)), mean(f(:,3)), pick, size(filtered_results,1));

    FINISHED = 'pass2 2026-09-15: validation recomputed client-side (DIAG); ENFORCE legs cancelled by Ben';
    save(file);
    fprintf('SAVED FULL WORKSPACE: %s\n', file);
end
fprintf('\nALL DONE\n');
