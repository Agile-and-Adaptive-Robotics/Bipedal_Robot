%% MinimizeExtPin
%Minimization scheme

clear; clc; close all

%% Initial Baseline Evaluation (for constraint bounds)
[a0, ~] = minimizeExtX3(0, Inf, Inf, 0, 1:4);  % All for baseline
baselineScores = a0;  % RMSE, FVU, Max Residual
fprintf('Baseline: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n', mean(baselineScores(:,1)),mean(baselineScores(:,2)),mean(baselineScores(:,3)));

load minimizeFlxPin10_results.mat sol_actual
sol_actual1 = sol_actual;
[a1, ~] = minimizeExtX3(sol_actual1(1), sol_actual1(2), sol_actual1(3), 0, 1:4);   % Use solution from Flexor bracket, and compare results
clear sol_actual
baselineScores1 = a1;  % RMSE, FVU, Max Residual
fprintf('Baseline using previous opt: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n', mean(baselineScores1(:,1)),mean(baselineScores1(:,2)),mean(baselineScores1(:,3)));

%% Cross-validation setup
numBPA = 4;
allBPA = 1:numBPA;

results_cv = cell(1, numBPA);
scores_cv = zeros(numBPA, 3);  % Will store RMSE, FVU, Max Resid for held-out validation
labels = ["42cm", "42cm-tendon", "46cm", "48cm"];
Color_var = nchoosek(1:4,3);  % Each row = 3 BPAs for training

%% Problem bounds
lb = [-0.01 * 100, 3,3, 0];   % [mm, log10(N/m), log10(N/m)]
ub = [0.03 * 100, 7, 7, 0.02*100];

%% Solver
for holdoutIdx = 1:numBPA
    fprintf('\n---- Cross-validation: Holding out BPA #%d (%s) ----\n', ...
        holdoutIdx, labels(holdoutIdx));

    trainIdx = setdiff(allBPA, holdoutIdx);
    baseline_train = baselineScores(trainIdx,:);

    % GA options
    opts = optimoptions('gamultiobj', ...
        'UseParallel', true, ...
        'Display', 'iter', ...
        'InitialPopulationRange',[lb; ub], ...
        'PopulationSize', 100, ...
        'MaxGenerations', 100, ...
        'MutationFcn', {@mutationadaptfeasible}, ...
        'CrossoverFraction', 0.8, ...
        'CrossoverFcn', {@crossoverscattered}, ...
        'FunctionTolerance', 1e-4);
%     goal = [0 0 0];
%     weight = [1 5 0.75];
    opts.HybridFcn = {@fgoalattain};
%     opts.OutputFcn = {@debugPop};
    % Run optimization
    [x, fvals] = gamultiobj(@(X) min1(X, trainIdx), 4, [], [], [], [], ...
        lb, ub, [], opts);

    % Evaluate each solution on held-out BPA
    valF = zeros(size(x,1), 3);
    for i = 1:size(x,1)
        valF(i,:) = min1(x(i,:), holdoutIdx);
    end

    % Store full set (no bestIdx decision now)
    results_cv{holdoutIdx}.optParams_all = x;        % Nx3
    results_cv{holdoutIdx}.trainScores_all = fvals;  % Nx3
    results_cv{holdoutIdx}.validation_all = valF;    % Nx3
    results_cv{holdoutIdx}.distance_all = vecnorm(fvals - valF, 2, 2);  % Nx1
    results_cv{holdoutIdx}.foldIdx = repmat(holdoutIdx, size(x,1), 1);  % Nx1
end
    
% [sol,fval,Pareto_front, Pareto_Fvals, exitflag,output] = GODLIKE(@min1,lb,ub,[],'NumObjectives',3,...
%                                          'algorithms', {'DE';'GA';'ASA';'PSO'},...
%                                          'display'   , 'plot',...
%                                          'popsize'   , 75 );
% x = Pareto_front;
% fvals = Pareto_Fvals;

%% === Compile All Pareto Candidates from Cross-Validation ===
all_candidates = [];  % Will collect [foldIdx, x(3), fvals(3), valF(3), dist]

for i = 1:numBPA
    fold = results_cv{i}.foldIdx;           % Nx1
    x = results_cv{i}.optParams_all;        % Nx3
    train = results_cv{i}.trainScores_all;  % Nx3
    val = results_cv{i}.validation_all;     % Nx3
    dist = results_cv{i}.distance_all;      % Nx1

    rows = [fold, x, train, val, dist];     % Nx11
    all_candidates = [all_candidates; rows];
end

%% Sort by validation distance first, then validation metrics
results = all_candidates;  % [fold, Xi0, Xi1, Xi2, train (3), val (3), dist]
results_sort = sortrows(results, [12 9 10 11 6 7 8]);  % sort by distance, then val, then train

%% De-normalize to get physical parameters
x_actual = [results_sort(:,2)/100, 10.^results_sort(:,3), 10.^results_sort(:,4), results_sort(:,5)/100];
results_sort_actual = [results_sort(:,1), x_actual, results_sort(:,6:end)];

%% Pick best solution (later, flexible)
pick = 1;
sol_actual = results_sort_actual(pick, 2:5);  % [Xi0, Xi1, Xi2]
[f, bpa] = minimizeExtX3(sol_actual(1), sol_actual(2), sol_actual(3), sol_actual(4), 1:4);  % [f: 4x3], [bpa: full struct]

fprintf('\nPerformance with sol_actual:\n');
disp(array2table(f, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}, ...
                    'RowNames', cellstr(labels')));


%% Plot results

%% --- Define color scheme and labels ---
c = cell(8,1);
c{1} = '#FFD700'; % gold → Hybrid
c{2} = '#FFB14E'; % orange
c{3} = '#FA8775'; % light orange
c{4} = '#EA5F94'; % pink
c{5} = '#CD34B5'; % magenta → Predicted
c{6} = '#9D02D7'; % magenta 2
c{7} = '#0000FF'; % indigo → Measured
c{8} = '#000000'; % black

labels = ["42cm", "42cm-tendon", "46cm", "48cm"];
tileLabels = {'\bf (A)', '\bf (B)', '\bf (C)', '\bf (D)'};
sz = 60;

%% --- Torque Figure with 2x2 tiles ---
figT = figure('Name','Torque - 2x2','Color','w');
figT.Position = [100 100 950 700];
tT = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

for i = 1:4
    ax = nexttile;
    hold on

    % Plot order: hybrid (gold), experimental (indigo), original (gray dash), predicted (magenta)
    scatter(bpa(i).A_h, bpa(i).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{1}, 'DisplayName', 'Hybrid');

    scatter(bpa(i).Aexp, bpa(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');

    plot(bpa(i).Ak, bpa(i).M, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, ...
        'DisplayName', 'Original');

    plot(bpa(i).Ak, bpa(i).M_p(:,3), '-', 'Color', c{5}, 'LineWidth', 2.5, ...
        'DisplayName', 'Predicted');

    % Tile-specific title and annotation label
    title(['\bf ' labels(i)], 'Interpreter', 'tex');
    
    % === Annotation Labels (A–D) ===
    annotation(figT, 'textbox', [0.035  0.89  0.05  0.05], 'String', '\bf (A)', ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    annotation(figT, 'textbox', [0.49  0.89  0.05  0.05], 'String', '\bf (B)', ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    annotation(figT, 'textbox', [0.035  0.41  0.05  0.05], 'String', '\bf (C)', ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    annotation(figT, 'textbox', [0.49  0.41  0.05  0.05], 'String', '\bf (D)', ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');


    % Axis config
    set(gca, ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial', ...
        'LineWidth', 2, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'YLim', [0 15], ...
        'TickLength', [0.025 0.05], ...
        'GridLineStyle','none');

    if ismember(i,[1,3])
        ylabel('\bf Torque, N\cdotm','Interpreter','tex')
    end
    if ismember(i,[3,4])
        xlabel('\bf \theta_{k} , \circ','Interpreter','tex')
    end
end

% Legend in top-right tile only
lg = legend(tT.Children(3));
lg.Location = 'northeast';
lg.FontSize = 8;

%% --- Muscle Length Figure with 2x2 tiles ---
figL = figure('Name','Muscle Length - 2x2','Color','w');
figL.Position = [100 100 950 700];
tL = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

for i = 1:4
    ax = nexttile;
    hold on

    % Predicted
    Lm_p = bpa(i).Lmt_p - 2 * bpa(i).fitn - bpa(i).ten;
    Lm   = bpa(i).Lmt   - 2 * bpa(i).fitn - bpa(i).ten;

    scatter(bpa(i).A_h, bpa(i).Lm_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}); % Hybrid (gold)

%     scatter(bpa(i).Aexp, bpa(i).Lm_exp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
%         'MarkerFaceColor', c{7}); % Measured (indigo)

    plot(bpa(i).Ak, Lm, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);      % Original
    plot(bpa(i).Ak, Lm_p, '-', 'Color', c{5}, 'LineWidth', 2.5);            % Predicted

    title(['\bf ' labels(i)], 'Interpreter', 'tex');

    % Annotation (A–D)
    annLabels = {'\bf (A)', '\bf (B)', '\bf (C)', '\bf (D)'};
    xAnn = [0.035, 0.49, 0.035, 0.49];
    yAnn = [0.89, 0.89, 0.41, 0.41];
    annotation(figL, 'textbox', [xAnn(i) yAnn(i) 0.05 0.05], ...
        'String', annLabels{i}, 'FontSize', 12, ...
        'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'TickLength', [0.025 0.05], 'GridLineStyle','none');

    if ismember(i,[1,3])
        ylabel('\bf Muscle Length, m','Interpreter','tex')
    end
    if ismember(i,[3,4])
        xlabel('\bf \theta_{k} , \circ','Interpreter','tex')
    end
end


%% --- Moment Arm Figure with 2x2 tiles ---
figMA = figure('Name','Moment Arm - 2x2','Color','w');
figMA.Position = [100 100 950 700];
tMA = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

for i = 1:4
    ax = nexttile;
    hold on

    G_p = hypot(bpa(i).mA_p(:,1), bpa(i).mA_p(:,2));

    scatter(bpa(i).A_h, bpa(i).mA_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7});  % Hybrid

%     scatter(bpa(i).Aexp, bpa(i).mAexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
%         'MarkerFaceColor', c{7});  % Measured

    plot(bpa(i).Ak, bpa(i).mA, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);      % Original
    plot(bpa(i).Ak, G_p, '-', 'Color', c{5}, 'LineWidth', 2.5);                   % Predicted

    title(['\bf ' labels(i)], 'Interpreter', 'tex');

    % Annotation (A–D)
    annLabels = {'\bf (A)', '\bf (B)', '\bf (C)', '\bf (D)'};
    xAnn = [0.035, 0.49, 0.035, 0.49];
    yAnn = [0.89, 0.89, 0.41, 0.41];
    annotation(figMA, 'textbox', [xAnn(i) yAnn(i) 0.05 0.05], ...
        'String', annLabels{i}, 'FontSize', 12, ...
        'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'TickLength', [0.025 0.05], 'GridLineStyle','none');

    if ismember(i,[1,3])
        ylabel('\bf Moment Arm, m','Interpreter','tex')
    end
    if ismember(i,[3,4])
        xlabel('\bf \theta_{k} , \circ','Interpreter','tex')
    end
end

%% --- Strain Figure with 2x2 tiles ---
figS = figure('Name','Strain - 2x2','Color','w');
figS.Position = [100 100 950 700];
tS = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

for i = 1:4
    ax = nexttile;
    hold on
    strain_h = (bpa(i).rest - bpa(i).Lm_h)/bpa(i).rest;
    scatter(bpa(i).A_h, strain_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{1});  % Hybrid

%     scatter(bpa(i).Aexp, bpa(i).strain_exp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
%         'MarkerFaceColor', c{7});  % Measured

    plot(bpa(i).Ak, bpa(i).strain, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);     % Original
    plot(bpa(i).Ak, bpa(i).strain_p, '-', 'Color', c{5}, 'LineWidth', 2.5);          % Predicted

    title(['\bf ' labels(i)], 'Interpreter', 'tex');

    % Annotation (A–D)
    annLabels = {'\bf (A)', '\bf (B)', '\bf (C)', '\bf (D)'};
    xAnn = [0.035, 0.49, 0.035, 0.49];
    yAnn = [0.89, 0.89, 0.41, 0.41];
    annotation(figS, 'textbox', [xAnn(i) yAnn(i) 0.05 0.05], ...
        'String', annLabels{i}, 'FontSize', 12, ...
        'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'TickLength', [0.025 0.05], 'GridLineStyle','none');

    if ismember(i,[1,3])
        ylabel('\bf Strain','Interpreter','tex')
    end
    if ismember(i,[3,4])
        xlabel('\bf \theta_{k} , \circ','Interpreter','tex')
    end
end

%% Helper functions
function ff = min1(x, trainIdx)
    if numel(x) == 4 && size(x,1) == 1
        % OK
    else
        error('min1: Input x must be a 1x4 vector');
    end
    Xi0 = x(1) / 100;
    Xi1 = 10^x(2);
    Xi2 = 10^x(3);
    Xi3 = x(4) / 100;
    [f_all, ~] = minimizeExtX3(Xi0, Xi1, Xi2, Xi3, trainIdx); % Nx3 matrix (e.g., 3x3 if 3 training BPAs)
    ff = mean(f_all, 1, 'omitnan'); % Return 1x3: [mean RMSE, mean FVU, mean MaxResidual]
end

%% --- Nonlinear constraint
function [c, ceq] = nonlinc(X, baseline, trainIdx)
    % Inputs:
    %   X         = [Xi0_cm, log10(Xi1), log10(Xi2), Xi3_cm]
    %   baseline  = baselineScores (size: 4x3)
    %   trainIdx  = which BPAs are being optimized
    try
        % Convert to physical parameters
        Xi0 = X(1) / 100;
        Xi1 = 10^X(2);
        Xi2 = 10^X(3);
        Xi3 = X(4) / 100;
        % Call model with current X on the training BPAs
        [f_all, ~] = minimizeExtX3(Xi0, Xi1, Xi2, Xi3, trainIdx);
        % beat the average baseline, not per-BPA
        mean_baseline = mean(baseline(trainIdx,:), 1, 'omitnan');  % 1x3
        mean_model = mean(f_all, 1, 'omitnan');                    % 1x3

        c = mean_model - mean_baseline;  % Element-wise (positive means violation)
        ceq = [];

        if any(c > 0)
            fprintf('[nonlinc] Violated mean at x = %.4f, %.4f, %.4f %.4f\n', X);
        end
    catch
        fprintf('[nonlinc] Error — invalid at x = %.4f, %.4f, %.4f %.4f\n', X);
        c = ones(4, 1) * 1e3;  % Large penalty
        ceq = [];
    end
end

function [state, options, optchanged] = debugPop(options, state, flag)
    optchanged = false;  % required for gamultiobj
    if strcmp(flag, 'iter') || strcmp(flag, 'done')
        if isfield(state, 'Population')
            fprintf('[debugPop] Generation %d\n', state.Generation);
            assignin('base', 'currentPop', state.Population);
            assignin('base', 'currentScores', state.Score);
        end
    end
end

function stop = gaplotpareto3D_simple(options, state, flag)
    stop = false;
    persistent figHandle
    if strcmp(flag, 'init') || isempty(figHandle) || ~isvalid(figHandle)
        figHandle = figure(99); clf; 
    end
    if strcmp(flag, 'iter') || strcmp(flag, 'done')
        scores = state.Score;
        if size(scores,2) == 3
            figure(figHandle); clf;
            scatter3(scores(:,1), scores(:,2), scores(:,3), 40, 'filled');
            xlabel('RMSE'); ylabel('FVU'); zlabel('Max Residual');
            title('Pareto Front (Training Set)');
            grid on; view(135, 30);
        end
    end
end