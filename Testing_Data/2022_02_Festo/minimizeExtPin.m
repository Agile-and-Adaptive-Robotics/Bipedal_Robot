
%Minimization scheme

clear; clc; close all

%% Initial Baseline Evaluation (for constraint bounds)
[a, b, c, d] = minimizeExt10mm(0, Inf, Inf, 1:4);  % All for baseline
baselineScores = [a; b; c; d];  % RMSE, FVU, Max Residual
fprintf('Baseline RMSE: %.4f ± %.4f\n', mean(baselineScores(:,1)), std(baselineScores(:,1)));

load minimizeFlxPin10_results.mat sol_actual
sol_actual1 = sol_actual;
[a1, b1, c1, d1] = minimizeExt10mm(sol_actual1(1), sol_actual1(2), sol_actual1(3), 1:4);   % Use solution from Flexor bracket, and compare results
clear sol_actual
baselineScores1 = [a1, b1, c1, d1];  % RMSE, FVU, Max Residual
fprintf('Baseline RMSE: %.4f ± %.4f\n', mean(baselineScores1(:,1)), std(baselineScores1(:,1)));

%% Cross-validation setup
numBPA = 4;
allBPA = 1:numBPA;

results_cv = cell(1, numBPA);
scores_cv = zeros(numBPA, 3);  % Will store RMSE, FVU, Max Resid for held-out validation

%% Problem bounds
lb = [-0.01 * 1000, 1e3, 1e3];   % [mm, log10(N/m), log10(N/m)]
ub = [0.03 * 1000, 1e8, 1e8];

%% Solver
for holdoutIdx = 1:numBPA
    fprintf('---- Cross-validation: Holding out BPA #%d ----\n', holdoutIdx);
    trainIdx = setdiff(allBPA, holdoutIdx);

    % Solver settings 
    opts = optimoptions('gamultiobj', ...
        'UseParallel', true, ...
        'Display', 'iter', ...
        'PlotFcn', {@gaplotpareto}, ...
        'PopulationSize', 100, ...
        'MaxGenerations', 600, ...
        'FunctionTolerance', 1e-4);
    
    % Optimization settings
    [x, fvals, exitflag, output, population, scores] = gamultiobj( ...
        @(X) min1(X, trainIdx), 3, ...
        [], [], [], [], ...
        lb, ub, ...
        @(X) nonlinc(X), ...
        opts);

    % Evaluate validation on held-out BPA
    valF = zeros(size(xOpt,1), 3);  % Validation RMSE, FVU, Max Resid
    for i = 1:size(xOpt,1)
        valF(i,:) = minimizeExt10mm(xOpt(i,1)/1000, xOpt(i,2), xOpt(i,3), holdoutIdx);
    end
    
    distances = vecnorm(fvals - valF, 2, 2);  % Error between training and validation
    [~, bestIdx] = min(distances);
    
    bestX = xOpt(bestIdx,:);
    scores_cv(holdoutIdx,:) = valF(bestIdx,:);
    
    % Store results
    results_cv{holdoutIdx}.optParams = bestX;
    results_cv{holdoutIdx}.validation = valF(bestIdx,:);
    results_cv{holdoutIdx}.fvals_train = fvals(bestIdx,:);
    results_cv{holdoutIdx}.population = population;
    
    % Plot results
    [~, ~, bpa] = minimizeExt10mm(bestX(1)/1000, bestX(2), bestX(3), 1:4);
    plotBPASet(bpa);
end
    
% [sol,fval,Pareto_front, Pareto_Fvals, exitflag,output] = GODLIKE(@min1,lb,ub,[],'NumObjectives',3,...
%                                          'algorithms', {'DE';'GA';'ASA';'PSO'},...
%                                          'display'   , 'plot',...
%                                          'popsize'   , 75 );
% x = Pareto_front;
% fvals = Pareto_Fvals;
%% Get validation GoF results from pareto front
val_Fvals = zeros(size(fvals,1),3);
for i = 1:length(x)
    val_Fvals(i,:) = min2([x(i,1),x(i,2),x(i,3)]);    %Get validation Fvals for all Pareto_front points
end

%% Sort results
ind = (1:length(x))';  %Index to original Pareto_front and Pareto_Fvals
relate = vecnorm(fvals-val_Fvals,2,2);   %Find the distance between the optimization and validation solutions for the same input
results = [ind, x, fvals, val_Fvals, relate]; 
results_sort = sortrows(results,[11 8 9 10 5 6 7]); %Sort results first on distance between optimization and validation, then on validation columns, then on original Fvals columns.
results_sort_actual = [results_sort(:,1), results_sort(:,2)/100, 10.^results_sort(:,3), 10.^results_sort(:,4), results_sort(:,5:end)];

%% Pick ultimate solution
pick = 1; %Pick the best solution from the sorted results (should be 1)
sol_actual = results_sort_actual(pick, 2:4);  %Best solution                                   
[~, ~, ~, bpa] = minimizeExt10mm(sol_actual(1), sol_actual(2), sol_actual(3), 1:4);           % Now pull bpa structures out       

%% Plot results
labels = ["42cm", "42cm-tendon", "46cm", "48cm"];
cmap = lines(4);
sz = 50;

for i = 1:4
    % --- Torque ---
    figure; hold on; 
    ax = gca; ax.FontWeight = 'bold'; ax.FontSize = 11; ax.LineWidth = 2;
    ax.XMinorTick = 'on'; ax.YMinorTick = 'on';

    plot(bpa(i).Ak, bpa(i).M, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Original');
    scatter(bpa(i).Aexp, bpa(i).Mexp, sz, 'MarkerFaceColor', '#CD34B5', 'DisplayName', 'Experiment');
    scatter(bpa(i).A_h, bpa(i).M_h, sz, 'MarkerFaceColor', '#FFD700', 'DisplayName', 'Hybrid');
    plot(bpa(i).Ak, bpa(i).M_p(:,3), '-', 'Color', cmap(i,:), 'LineWidth', 2.5, 'DisplayName', 'Predicted');

    title(['Torque - ' labels(i)]);
    xlabel('\theta_k, °'); ylabel('Torque (Nm)');
    legend('Location', 'best'); ylim padded;

    % --- Muscle Length ---
    figure; hold on; 
    ax = gca; ax.FontWeight = 'bold'; ax.FontSize = 11; ax.LineWidth = 2;
    ax.XMinorTick = 'on'; ax.YMinorTick = 'on';

    Lm   = bpa(i).Lmt - 2 * bpa(i).fitn - bpa(i).ten;
    Lm_p = bpa(i).Lmt_p - 2 * bpa(i).fitn - bpa(i).ten;

    plot(bpa(i).Ak, Lm_p, '-', 'Color', cmap(i,:), 'LineWidth', 2.5, 'DisplayName', 'Predicted');
    scatter(bpa(i).A_h, bpa(i).Lm_h, sz, 'MarkerFaceColor', '#CD34B5', 'DisplayName', 'Measured');
    plot(bpa(i).Ak, Lm, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Original');

    title(['Muscle Length - ' labels(i)]);
    xlabel('\theta_k, °'); ylabel('Length (m)');
    legend('Location', 'best'); ylim padded;

    % --- Moment Arm ---
    figure; hold on; 
    ax = gca; ax.FontWeight = 'bold'; ax.FontSize = 11; ax.LineWidth = 2;
    ax.XMinorTick = 'on'; ax.YMinorTick = 'on';

    G_p = hypot(bpa(i).mA_p(:,1), bpa(i).mA_p(:,2));
    plot(bpa(i).Ak, G_p, '-', 'Color', cmap(i,:), 'LineWidth', 2.5, 'DisplayName', 'Predicted');
    scatter(bpa(i).A_h, bpa(i).mA_h, sz, 'MarkerFaceColor', '#CD34B5', 'DisplayName', 'Measured');
    plot(bpa(i).Ak, bpa(i).mA, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Original');

    title(['Moment Arm - ' labels(i)]);
    xlabel('\theta_k, °'); ylabel('Moment Arm (m)');
    legend('Location', 'best'); ylim padded;
end

%% Helper functions
function ff = min1(x, trainIdx)
    Xi0 = x(1) / 1000;
    Xi1 = 10^x(2);
    Xi2 = 10^x(3);
    ff = minimizeExt10mm(Xi0, Xi1, Xi2, trainIdx);
end

%% --- Nonlinear constraint
function [c, ceq] = nonlinc(X, baseline, trainIdx)
    f = min1(X, trainIdx);
    c = f - mean(baseline(trainIdx,:),1);
    ceq = [];
end

function plotBPASet(bpa)
    labels = ["42cm", "42cm-tendon", "46cm", "48cm"];
    for i = 1:numel(bpa)
        figure;
        hold on; grid on;
        plot(bpa(i).Ak, bpa(i).M_p(:,3), 'b', 'LineWidth', 2, 'DisplayName', 'Predicted');
        scatter(bpa(i).Aexp, bpa(i).Mexp, 40, 'k', 'filled', 'DisplayName', 'Experiment');
        plot(bpa(i).Ak, bpa(i).M, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName', 'Original');
        title(['Torque Prediction - ', labels(i)]);
        xlabel('\theta_k (°)'); ylabel('Torque (Nm)');
        legend('Location','best'); ylim padded;
        set(gca, 'FontWeight','bold', 'FontSize',12, 'LineWidth',1.5, 'Box','on', 'XMinorTick','on', 'YMinorTick','on');
    end
end