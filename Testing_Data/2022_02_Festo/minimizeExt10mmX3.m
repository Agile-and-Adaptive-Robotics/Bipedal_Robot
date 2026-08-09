%% MinimizeExtPinX3.m
%Minimization scheme

clear; clc; close all

%% Cross-validation setup
allBPA = [1, 2, 3, 4, 5, 6, 7, 8, 9];           %We're going to skip 42cm l_rest with a tendon
labels = ["40cm", "40cm-tendon", "42cm", "42cm-tendon", "43cm", "43cm-tendon", "46cm", "47cm", "48cm"];
validLabels = labels(allBPA);
numBPA = numel(allBPA);

results_cv = cell(1, numBPA);  % Will store RMSE, FVU, Max Resid for BPA(s) optimized
scores_cv = zeros(numBPA, 3);  % Will store RMSE, FVU, Max Resid for BPA(s) held-out for validation

%% Initial Baseline Evaluation (for constraint bounds)
[a0, bpa0] = minimizeExtX3(0, Inf, Inf, 0);  % All for baseline
baselineScores = a0;  % RMSE, FVU, Max Residual
fprintf('\nPerformance with no length offset and infinite stiffness:\n');
disp(array2table(a0, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}, ...
                    'RowNames', cellstr(labels')));
fprintf('Mean Baseline: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n\n', mean(baselineScores(:,1)),mean(baselineScores(:,2)),mean(baselineScores(:,3)));

load minimizeFlxPin10_results_20260729.mat xCols filtered_results
pick = 1; %Pick the best solution from the sorted results (should be 1)
sol_actual = filtered_results(pick, xCols);  %Best solution
g = sol_actual;
[a1, ~] = minimizeExtX3(-g(1), g(2), g(3), 0);   % Use solution from Flexor bracket, and compare results
[a2, bpa2] = minimizeExtX3(-g(1), g(2), g(3), 0.2);   % Use solution from Flexor bracket, reverse length offset, and guess for Xi3
clear sol_actual pick filtered_results xCols
fprintf('\nBaseline using previous opt:\n');
disp(array2table(a1, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}, ...
                    'RowNames', cellstr(labels')));
baselineScores1 = a1./(a0);  % RMSE, FVU, Max Residual, normalized to baselineScores
fprintf('Mean normalized to baseline score: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n\n', mean(baselineScores1(:,1)),mean(baselineScores1(:,2)),mean(baselineScores1(:,3)));
fprintf('\nBaseline using best guess:\n');
disp(array2table(a2, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}, ...
                    'RowNames', cellstr(labels')));
baselineScores2 = a2./(a0); % RMSE, FVU, Max Residual, normalized to baselineScores
fprintf('Mean normalized to baseline score: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n', mean(baselineScores2(:,1)),mean(baselineScores2(:,2)),mean(baselineScores2(:,3)));


%% Problem bounds
lb = [-0.020 * 100, log10(5e3), log10(5e3), 0];   % [cm, log10(N/m), log10(N/m), unitless]
ub = [0.01 * 100, log10(5e7), log10(5e7), 3];

% lb = [-g(1) * 100, log10(g(2)), log10(g(3)), 0];   % [cm, log10(N/m), log10(N/m), unitless]
% ub = [-g(1) * 100, log10(g(2)), log10(g(3)), 3];

% A = [0 -1 1 0; ...              % x2 (bending) is less stiff than x1 (axial), (x2 <= x1)
%      0 0 0 0; ...
%      0 0 0 0];
%  b = [0, 0, 0, 0]';

clear sol_actual
%% Solver
numHold = 4;                        %Number of BPAs held out for validation
list = nchoosek(allBPA,numHold);          %Choose how many BPAs to hold out, the others for training
for k = 1:length(list)
    holdoutIdx = list(k,:);
    for n = 1:size(holdoutIdx,2)
        fprintf('\n---- Cross-validation: Holding out BPA #%d (%s) ----\n', ...
            holdoutIdx(n), labels(holdoutIdx(n)));
    end
    trainIdx = setdiff(allBPA, holdoutIdx);
    baseline_train = baselineScores(trainIdx,:);

    % GA options
    opts = optimoptions('gamultiobj', ...
        'UseParallel', true, ...
        'Display', 'iter', ...
        'PlotFcn', {@gaplotpareto3D_simple}, ...    %'InitialPopulationRange',[-.015*100, log10(8e4), log10(8e3), 0.1; -0.007*100, log10(8e6), log10(8e5), 0.4], ...
        'PopulationSize', 150, ... %was 150
        'MaxGenerations', 750, ... %was 750
        'MutationFcn', {@mutationadaptfeasible}, ...
        'CrossoverFraction', 0.8, ...
        'CrossoverFcn', {@crossoverscattered}, ...
        'FunctionTolerance', 1e-3);
% hybridOpts = optimoptions('fgoalattain', ...
%     'Display','iter');

% opts.HybridFcn = {@fgoalattain, hybridOpts};
%     opts.OutputFcn = {@debugPop};

    % Run optimization
     [x, fvals,exitflag,output,population,scores] = gamultiobj(@(X) min1(X, trainIdx, a0), 4, [], [], [], [], ...
                                                    lb, ub, ... % @(x) nonlcon2(x), ...
                                                    opts);
    
% If I want to use GODLIKE instead:  
%     trainIdx_fixed = trainIdx;  % capture loop var for closure
%     min1_wrapper = @(x) min1(x, trainIdx_fixed);
%     [sol,fval,x, fvals, exitflag,output] = GODLIKE(min1_wrapper,lb,ub,[],'NumObjectives',3,...
%                                              'algorithms', {'PSO';'GA';'DE';'ASA'},...
%                                              'display'   , 'plot', ...
%                                              'UseParallel','true', ...
%                                              'popsize'   , 100 );

    % Evaluate each solution on held-out BPA
    valF = zeros(size(x,1), 3);
    parfor i = 1:size(x,1)
        valF(i,:) = min1(x(i,:), holdoutIdx, a0);
    end

    % Store full set (no bestIdx decision now)
    
    results_cv{k}.optParams_all = x;        % Nx4
    results_cv{k}.trainScores_all = fvals;  % Nx3
    results_cv{k}.validation_all = valF;    % Nx3
    results_cv{k}.distance_all = vecnorm(fvals - valF, 2, 2);  % Nx1
    results_cv{k}.foldIdx = repmat(holdoutIdx, size(x,1), 1);  % Nx1
end

%% === Compile All Pareto Candidates from Cross-Validation ===
all_candidates = [];  % Will collect [foldIdx, x(3), fvals(3), valF(3), dist]

for i = 1:numBPA
    fold = results_cv{i}.foldIdx;           % NxnumHold
    x2 = results_cv{i}.optParams_all;        % Nx4
    train = results_cv{i}.trainScores_all;  % Nx3
    val = results_cv{i}.validation_all;     % Nx3
    dist = results_cv{i}.distance_all;      % Nx1
    ind = 1:length(x2);                         %create an index
    ind = ind';                                 %Make Nx1 column array to show original results order
    rows = [ind, fold, x2, train, val, dist];   % Nx(12+numHold)
    all_candidates = [all_candidates; rows];
end

%% Dynamic column indices before de-normalizing
rankCol   = 1;
holdCols  = numel(rankCol) + (1:numHold);
xCols     = numel(rankCol) + numel(holdCols) + (1:4);
trainCols = numel(rankCol) + numel(holdCols) + numel(xCols) + (1:3);
valCols   = numel(rankCol) + numel(holdCols) + numel(xCols) + numel(trainCols) + (1:3);
distCol   = numel(rankCol) + numel(holdCols) + numel(xCols) + numel(trainCols) + numel(valCols) + 1;

%% Sort by validation distance first, then validation metrics
results = all_candidates;  % [fold, Xi0, Xi1, Xi2, Xi3, train (3), val (3), dist]
results_sort = sortrows(results, [distCol valCols trainCols]);  % sort by distance, then val, then train

%% De-normalize to get physical parameters
x_actual = [results_sort(:,xCols(1))/100, 10.^results_sort(:,xCols(2)), 10.^results_sort(:,xCols(3)), results_sort(:,xCols(4))];
results_sort_actual = [results_sort(:,rankCol), results_sort(:,holdCols), x_actual, results_sort(:,trainCols), results_sort(:,valCols), results_sort(:,distCol)];

%% --- Filter Pareto candidates against baseline on BPAs 1, 3 & 4 ---
N = size(results_sort_actual, 1);
keep = false(N,1);

for ii = 1:N
    % extract decision variables
    Xi0 = results_sort_actual(ii,xCols(1));
    Xi1 = results_sort_actual(ii,xCols(2));
    Xi2 = results_sort_actual(ii,xCols(3));
    Xi3 = results_sort_actual(ii,xCols(4));

    % re-evaluate on all 9 BPAs
    f_all = minimizeExtX3(Xi0, Xi1, Xi2, Xi3);   % returns 4×3 [RMSE, FVU, MaxResidual]

    % compare RMSE & FVU for BPAs 1, 3 & 4 to baselineScores
    pass1 = all( f_all(1,:) <= baselineScores(1,:) );
    pass2 = all( f_all(2,:) <= baselineScores(2,:) );
    pass3 = all( f_all(3,:) <= baselineScores(3,:) );
    pass4 = all( f_all(4,:) <= baselineScores(4,:) );
    pass5 = all( f_all(5,:) <= baselineScores(5,:) );
    pass6 = all( f_all(6,:) <= baselineScores(6,:) );
    pass7 = all( f_all(7,:) <= baselineScores(7,:) );
    pass8 = all( f_all(8,:) <= baselineScores(8,:) );
    pass9 = all( f_all(9,:) <= baselineScores(9,:) );

    keep(ii) = pass1 && pass2 && pass5 && pass6;
end

% keep only the rows that passed all desired checks
filtered_results = results_sort_actual(keep, :);
fprintf('Filtered %d → %d candidates.\n', N, sum(keep));


%% Pick best solution (later, flexible)
 
pick = 1;
sol_actual = filtered_results(pick, xCols);
[f, bpa] = minimizeExtX3(sol_actual(1), sol_actual(2), sol_actual(3), sol_actual(4));  % [f: 4x3], [bpa: full struct]

fprintf('\nPerformance with no length offset and infinite stiffness:\n');
disp(array2table(a0, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}, ...
                    'RowNames', cellstr(labels')));
fprintf('Mean Baseline: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n\n', mean(baselineScores(:,1)),mean(baselineScores(:,2)),mean(baselineScores(:,3)));

fprintf('\nPerformance with sol_actual:\n');
disp(array2table(f, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}, ...
                    'RowNames', cellstr(labels')));
fprintf('Mean GoF: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n\n', mean(f(:,1)),mean(f(:,2)),mean(f(:,3)));

disp(array2table(sol_actual, 'VariableNames', {'X0', 'X1', 'X2','X3'}, ...
'RowNames', {'Solution'}))
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

tileIdxs = [1, 5, 8, 12, 15, 19, 22, 26, 29]; 
% tileIdxs = [1, 5, 8, 12];  % A, B, C, D
% tileIdxs = [1, 5, 10];  % A, B, C
el = numel(tileIdxs);
tileSpans = [1 3];      % Span: [rows cols]
% tileSpans = [1 1];      % Span: [rows cols]
tileOrder = [1, 2, 3, 4, 5, 6, 7, 8, 9];
% tileOrder = [4, 3, 1, 2];
% tileLabels = {'(A)', '(B)', '(C)', '(D)'};
% Annotation positions [x, y] in normalized figure units
% xAnn = [0.035, 0.51, 0.035, 0.51];  % (A), (B), (C), (D)
% yAnn = [0.89, 0.89, 0.41, 0.41];    % (A), (B), (C), (D)
% xAnn = [0.035, 0.51, 0.265];  % (A), (B), (C)
% yAnn = [0.89, 0.89, 0.41];    % (A), (B), (C)
sz = 60;

%for plotting 

%% --- Torque Figure with tiles ---
figT = figure('Name','Torque','Color','w');
figT.Position = [100 100 950 700];
tT = tiledlayout(ceil(el/2),7,'TileSpacing','tight','Padding','tight');

for j = 1:el
    i = tileOrder(j);
    ax = nexttile(tileIdxs(j), tileSpans);
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
    title(['\bf ' labels(i)], 'Interpreter','tex');
    % ylabel('\bf Torque, N\cdotm','Interpreter','tex')
    % xlabel('\bf \theta_{k} , \circ','Interpreter','tex')
    % annotation(figT, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
    %     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

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
end

%shared axes labels
ylabel(tT,'\bf Torque, N\cdotm','Interpreter','tex')
xlabel(tT,'\bf \theta_{k} , \circ','Interpreter','tex')

% Legend in top-right tile only
lg = legend(tT.Children(1));
lg.Location = 'northeast';
lg.FontSize = 8;

%% --- Muscle Length Figure with tiles---
figL = figure('Name','Muscle Length','Color','w');
figL.Position = [100 100 950 700];
tL = tiledlayout(ceil(el/2),7,'TileSpacing','tight','Padding','tight');

for j = 1:el
    i = tileOrder(j);
    ax = nexttile(tileIdxs(j), tileSpans);
    hold on

    % Predicted
    Lm_p = bpa(i).Lmt_p - 2 * bpa(i).fitn - bpa(i).ten;
    Lm   = bpa(i).Lmt   - 2 * bpa(i).fitn - bpa(i).ten;

    scatter(bpa(i).A_h, bpa(i).Lm_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7},'DisplayName', 'Measured'); % Hybrid (gold)
    plot(bpa(i).Ak, Lm, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2,'DisplayName', 'Original');      % Original
    plot(bpa(i).Ak, Lm_p, '-', 'Color', c{5}, 'LineWidth', 2.5,'DisplayName', 'Predicted');            % Predicted
    % ylabel('\bf Length','Interpreter','tex')
    % xlabel('\bf \theta_{k} , \circ','Interpreter','tex')

    % Tile-specific title and annotation label
    title(['\bf ' labels(i)], 'Interpreter','tex');
    % annotation(figL, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
    %     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    % Axis config
    set(gca, ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial', ...
        'LineWidth', 2, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'TickLength', [0.025 0.05], ...
        'GridLineStyle','none');

end

%shared axes labels
ylabel(tL,'\bf Length','Interpreter','tex')
xlabel(tL,'\bf \theta_{k} , \circ','Interpreter','tex')

% Legend in top-right tile only
lg = legend(tL.Children(1));
lg.Location = 'northeast';
lg.FontSize = 8;

%% --- Moment Arm Figure with tiles ---
figMA = figure('Name','Moment Arm','Color','w');
figMA.Position = [100 100 950 700];
tMA = tiledlayout(ceil(el/2),7,'TileSpacing','tight','Padding','tight');

for j = 1:el
    i = tileOrder(j);
    ax = nexttile(tileIdxs(j), tileSpans);
    hold on

    G_p = hypot(bpa(i).mA_p(:,1), bpa(i).mA_p(:,2));

    scatter(bpa(i).A_h, bpa(i).mA_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7},'DisplayName', 'Measured');  % Hybrid
    plot(bpa(i).Ak, bpa(i).mA, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2,'DisplayName', 'Original');      % Original
    plot(bpa(i).Ak, G_p, '-', 'Color', c{5}, 'LineWidth', 2.5,'DisplayName', 'Predicted');                   % Predicted

    % ylabel(ax,'\bf Moment arm, m','Interpreter','tex')
    % xlabel(ax,'\bf \theta_{k} , \circ','Interpreter','tex')

    % Tile-specific title and annotation label
    title(['\bf ' labels(i)], 'Interpreter','tex');
    % % annotation(figMA, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
    % %     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    % Axis config
    set(gca, ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial', ...
        'LineWidth', 2, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'TickLength', [0.025 0.05], ...
        'GridLineStyle','none');
end

%shared axes labels
ylabel(tMA,'\bf Moment arm, m','Interpreter','tex')
xlabel(tMA,'\bf \theta_{k} , \circ','Interpreter','tex')

% Legend in top-right tile only
lg = legend(tMA.Children(1));
lg.Location = 'northeast';
lg.FontSize = 8;

%% --- Strain Figure with tiles ---
figS = figure('Name','Relative Strain','Color','w');
figS.Position = [100 100 950 700];
tS = tiledlayout(ceil(el/2),7,'TileSpacing','tight','Padding','tight');

for j = 1:el
    i = tileOrder(j);
    ax = nexttile(tileIdxs(j), tileSpans);
    hold on
    
    strain_h = (bpa(i).rest - bpa(i).Lm_h)/bpa(i).rest;
    kmax = (bpa(i).rest - bpa(i).Kmax)/bpa(i).rest;
    scatter(bpa(i).A_h, strain_h/kmax, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{7},'DisplayName', 'Measured');
    plot(bpa(i).Ak, bpa(i).strain/kmax, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2,'DisplayName', 'Original');
    plot(bpa(i).Ak, bpa(i).strain_p/kmax, '-', 'Color', '#CD34B5', 'LineWidth', 2.5,'DisplayName', 'Predicted');
    % ylabel(tS,'\bf \epsilon^*','Interpreter','tex')
    % xlabel(tS,'\bf \theta_{k} , \circ','Interpreter','tex')

    % Tile-specific title and annotation label
    title(['\bf ' labels(i)], 'Interpreter','tex');
    % annotation(figS, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
    %     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    % Axis config
    set(gca, ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial', ...
        'LineWidth', 2, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'TickLength', [0.025 0.05], ...
        'GridLineStyle','none');
end

%shared axes labels
ylabel(tS,'\bf \epsilon^*','Interpreter','tex')
xlabel(tS,'\bf \theta_{k} , \circ','Interpreter','tex')

% Legend in top-right tile only
lg = legend(tS.Children(1));
lg.Location = 'northeast';
lg.FontSize = 8;
%% Helper functions
function ff = min1(x, trainIdx, kompare)
    if numel(x) == 4 && size(x,1) == 1
        % OK
    else
        error('min1: Input x must be a 1x4 vector');
    end
    Xi0 = x(1) / 100;
    Xi1 = 10^x(2);
    Xi2 = 10^x(3);
    Xi3 = x(4);
    try
        [f_all, ~] = minimizeExtX3(Xi0, Xi1, Xi2, Xi3, trainIdx); % Nx3 matrix (e.g., 3x3 if 3 training BPAs)
%         ff = mean(f_all, 1, 'omitnan');              % Return 1x3: [mean RMSE, mean FVU, mean MaxResidual]
        fnorm = f_all(trainIdx,:)./kompare(trainIdx,:);     %normalize results before taking the mean
        ff = mean(fnorm, 1, 'omitnan');              % Return 1x3: [mean RMSE, mean FVU, mean MaxResidual]
        if ~isnumeric(ff) || numel(ff) ~= 3
            ff = [Inf, Inf, Inf];  % Defensive return if shape is wrong
        end
    catch
        ff = [Inf, Inf, Inf];      % Defensive return if minimizeExtX3 throws
    end
end


%% --- Nonlinear constraint
function [c, ceq] = nonlinc(X, baseline, trainIdx)
    % Inputs:
    %   X         = [Xi0_cm, log10(Xi1), log10(Xi2), , Xi3_%]
    %   baseline  = baselineScores (size: 4x3)
    %   trainIdx  = which BPAs are being optimized
    try
        % Convert to physical parameters
        Xi0 = X(1) / 100;
        Xi1 = 10^X(2);
        Xi2 = 10^X(3);
        Xi3 = X(4);
        % Call model with current X on the training BPAs
        [f_all, ~] = minimizeExtX3(Xi0, Xi1, Xi2, Xi3, trainIdx);
        % beat the average baseline, not per-BPA
        mean_baseline = mean(baseline(trainIdx,:), 1, 'omitnan');  % 1x3
        mean_model = mean(f_all, 1, 'omitnan');                    % 1x3

        c = mean_model - mean_baseline;  % Element-wise (positive means violation)
        ceq = [];

        if any(c > 0)
%             fprintf('[nonlinc] Violated mean at x = %.4f, %.4f, %.4f %.4f\n', X);
        end
    catch
%         fprintf('[nonlinc] Error — invalid at x = %.4f, %.4f, %.4f %.4f\n', X);
        c = ones(4, 1) * 1e3;  % Large penalty
        ceq = [];
    end
end

function [c, ceq] = nonlcon2(x)
    % Inequality constraints (c <= 0)
    c = x(3) - x(2);  % This ensures x(3) < x(2)

    % No equality constraints
    ceq = [];
end

% function [state, options, optchanged] = debugPop(options, state, flag)
%     optchanged = false;  % required for gamultiobj
%     if strcmp(flag, 'iter') || strcmp(flag, 'done')
%         if isfield(state, 'Population')
%             fprintf('[debugPop] Generation %d\n', state.Generation);
%             assignin('base', 'currentPop', state.Population);
%             assignin('base', 'currentScores', state.Score);
%         end
%     end
% end

function [state, options, optchanged] = gaplotpareto3D_simple(options, state, flag)
    optchanged = false;  % Must be returned even if unchanged
    persistent figHandle
    if strcmp(flag, 'init') || isempty(figHandle) || ~isvalid(figHandle)
        figHandle = figure(99); 
        set(figHandle, 'Name', 'Live Pareto Front', 'NumberTitle', 'off');
    end
    if strcmp(flag, 'iter') || strcmp(flag, 'done')
        scores = state.Score;
        if ~isempty(scores) && isnumeric(scores) && size(scores,2) == 3
            figure(figHandle); 
            scatter3(scores(:,1), scores(:,2), scores(:,3), 50, 'filled');
            xlabel('\bf RMSE', 'FontSize', 12);
            ylabel('\bf FVU', 'FontSize', 12);
            zlabel('\bf Max Residual', 'FontSize', 12);
            title('\bf Pareto Front (Training Set)', 'FontSize', 14);
            grid on;
            view(135, 30);
            drawnow;
        end
    end
end


