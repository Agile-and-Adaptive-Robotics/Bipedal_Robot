%%minimizeFlxPin10mm_2brk.m
%Cross-validation driver for the two-bracket flexor evaluator
%minimizeFlxPin2brk (changes (d) single transform, (e) K=[X1,X2,X1],
%(f) second origin-side bracket). Copy of minimizeFlxPin10mm.m with:
%  - evaluator swapped to minimizeFlxPin2brk (USE_BRACKET2 threaded through)
%  - settings block below (MODE env var: 'smoke' default, 'full' for production)
%  - pick used STRICTLY: evaluated solution is always filtered_results(pick, xCols)
%  - solver switch: FLX2BRK_SOLVER = 'gamultiobj' (default) | 'surrogateopt'
%    (surrogateopt minimizes the baseline-normalized weighted sum)
%  - results saved to minimizeFlxPin10_2brk_results_<yyyymmdd>.mat

clear; clc; close all

%% Settings
MODE = getenv('FLX2BRK_MODE');
if isempty(MODE), MODE = 'smoke'; end
isSmoke = strcmpi(MODE, 'smoke');

SOLVER = getenv('FLX2BRK_SOLVER');
if isempty(SOLVER), SOLVER = 'gamultiobj'; end

USE_BRACKET2 = true;     %false = (d)+(e) ablation: no second bracket
DO_PLOTS = ~batchStartupOptionUsed;   %auto: plots when run interactively
PICK = 1;                %which filtered Pareto candidate to evaluate

if isSmoke
    ALLBPA   = [2, 3, 4, 5];
    NUMHOLD  = 1;
    POP      = 25;
    MAXGEN   = 30;
    SURRVALS = 1500;
else
    ALLBPA   = [1, 2, 3, 4, 5];           % Use if all data are valid
    % ALLBPA = [2, 3, 4, 5];              % Use if old data do not hold up
    NUMHOLD  = 2;
    POP      = 150;
    MAXGEN   = 600;
    SURRVALS = 6000;
end

labels = ["48cm", "46cm", "47cm", "40cm-tendon", "41cm"];
validLabels = labels(ALLBPA);
numBPA = numel(ALLBPA);

fprintf('=== minimizeFlxPin10mm_2brk | MODE=%s | SOLVER=%s | USE_BRACKET2=%d ===\n', MODE, SOLVER, USE_BRACKET2);
fprintf('ALLBPA=[%s] NUMHOLD=%d POP=%d MAXGEN=%d\n', num2str(ALLBPA), NUMHOLD, POP, MAXGEN);

results_cv = cell(1, numBPA);  % Will store RMSE, FVU, Max Resid for BPA(s) optimized
scores_cv = zeros(numBPA, 3);  % Will store RMSE, FVU, Max Resid for BPA(s) held-out for validation

%% Calculate baseline
[a0, bpa0] = minimizeFlxPin2brk(0,Inf,Inf,[],USE_BRACKET2);   %no extra length, infinite bracket stiffness
baselineScores = a0;
fprintf('\nPerformance with no length offset and infinite stiffness:\n');
disp(array2table(a0, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}, ...
                    'RowNames', cellstr(labels')));
fprintf('Mean baseline training: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n\n',mean(baselineScores(:,1)),mean(baselineScores(:,2)),mean(baselineScores(:,3)));


%% Problem bounds
lb = [0*100, log10(5e3), log10(5e3)];
ub = [0.020*100, log10(5e7), log10(5e7)];

%% Solver
list = nchoosek(ALLBPA,NUMHOLD);          %Choose how many BPAs to hold out, the others for training
W = [0.5, 0.3, 0.2];                      %scalarization weights for surrogateopt
for k = 1:length(list)
    holdoutIdx = list(k,:);
    for n = 1:size(holdoutIdx,2)
        fprintf('\n---- Cross-validation: Holding out BPA #%d (%s) ----\n', ...
            holdoutIdx(n), labels(holdoutIdx(n)));
    end
    trainIdx = setdiff(ALLBPA, holdoutIdx);
    baseline_train = baselineScores(trainIdx,:);

    opts = optimoptions('gamultiobj', ...
            'UseParallel', true, ...
            'Display', 'iter', ...
            'PopulationSize', POP, ...
            'MaxGenerations', MAXGEN, ...
            'MutationFcn', {@mutationadaptfeasible}, ...
            'CrossoverFraction', 0.8, ...
            'CrossoverFcn', {@crossoverscattered}, ...
            'FunctionTolerance', 4e-3);
    if DO_PLOTS
        opts.PlotFcn = {@gaplotpareto3D_simple};
    end

    switch lower(SOLVER)
        case 'gamultiobj'
            [x, fvals,exitflag,output,population,scores] = gamultiobj(@(X) min1(X, trainIdx, a0, USE_BRACKET2), 3, [], [], [], [], ...
                                                            lb, ub, opts);
        case 'surrogateopt'
            surrOpts = optimoptions('surrogateopt', ...
                'UseParallel', true, ...
                'Display', 'iter', ...
                'MaxFunctionEvaluations', SURRVALS);
            sol = surrogateopt(@(X) min1scalar(X, trainIdx, a0, W, USE_BRACKET2), lb, ub, surrOpts);
            %R2025a: first output IS the solution point (empty if all evals fail)
            assert(~isempty(sol), 'surrogateopt returned empty - all evaluations failed');
            x = reshape(sol, 1, []);
            fvals = min1(x, trainIdx, a0, USE_BRACKET2);
        otherwise
            error('Unknown FLX2BRK_SOLVER "%s": use gamultiobj or surrogateopt', SOLVER);
    end

    % Evaluate each solution on held-out BPA
    valF = zeros(size(x,1), 3);
    parfor i = 1:size(x,1)
        valF(i,:) = min1(x(i,:), holdoutIdx, a0, USE_BRACKET2);
    end

    % Store full set (no bestIdx decision now)
    results_cv{k}.optParams_all = x;        % Nx3
    results_cv{k}.trainScores_all = fvals;  % Nx3
    results_cv{k}.validation_all = valF;    % Nx3
    results_cv{k}.distance_all = vecnorm(fvals - valF, 2, 2);  % Nx1
    results_cv{k}.foldIdx = repmat(holdoutIdx, size(x,1), 1);  % Nx1
end
%% === Compile All Pareto Candidates from Cross-Validation ===
all_candidates = [];  % Will collect [ind, foldIdx, x(3), fvals(3), valF(3), dist]

for i = 1:numBPA
    fold = results_cv{i}.foldIdx;               % NxnumHold
    x2 = results_cv{i}.optParams_all;           % Nx3
    train = results_cv{i}.trainScores_all;      % Nx3
    val = results_cv{i}.validation_all;         % Nx3
    dist = results_cv{i}.distance_all;          % Nx1
    ind = 1:length(x2);                         %create an index
    ind = ind';                                 %Make Nx1 column array to show original results order
    rows = [ind, fold, x2, train, val, dist];   % Nx(11+numHold)
    all_candidates = [all_candidates; rows];
end

%% Dynamic column indices before de-normalizing
rankCol   = 1;
holdCols  = numel(rankCol) + (1:NUMHOLD);
xCols     = numel(rankCol) + numel(holdCols) + (1:3);
trainCols = numel(rankCol) + numel(holdCols) + numel(xCols) + (1:3);
valCols   = numel(rankCol) + numel(holdCols) + numel(xCols) + numel(trainCols) + (1:3);
distCol   = numel(rankCol) + numel(holdCols) + numel(xCols) + numel(trainCols) + numel(valCols) + 1;

%% Sort by validation distance first, then validation metrics
results = all_candidates;
results_sort = sortrows(results, [distCol valCols(1:2) trainCols(1:2)]);  % sort by distance, then val, then train (excluding Max. Residual)

%% De-normalize to get physical parameters
x_actual = [results_sort(:,xCols(1))/100, 10.^results_sort(:,xCols(2)), 10.^results_sort(:,xCols(3))];
results_sort_actual = [results_sort(:,rankCol), results_sort(:,holdCols), x_actual, results_sort(:,trainCols), results_sort(:,valCols), results_sort(:,distCol)];

%% --- Filter Pareto candidates against baseline on all BPAs ---
N = size(results_sort_actual, 1);
keep = false(N,1);

for ii = 1:N
    % extract decision variables
    Xi0 = results_sort_actual(ii,xCols(1));
    Xi1 = results_sort_actual(ii,xCols(2));
    Xi2 = results_sort_actual(ii,xCols(3));

    % re-evaluate on all BPAs
    f_all = minimizeFlxPin2brk(Xi0, Xi1, Xi2, [], USE_BRACKET2);   % returns nBPA x 3 [RMSE, FVU, MaxResidual]

    pass = true;
    for j = 1:numBPA
        pass = pass && all( f_all(j,1:3) <= baselineScores(j,1:3) );
    end
    keep(ii) = pass;
end

% keep only the rows that passed all checks
filtered_results = results_sort_actual(keep, :);
fprintf('Filtered %d → %d candidates.\n', N, sum(keep));


%% Pick best solution (strictly by PICK; no hardcoded overrides)
if isempty(filtered_results)
    warning('No candidate passed the baseline filter — pick from results_sort_actual instead.');
    sol_actual = results_sort_actual(min(PICK, size(results_sort_actual,1)), xCols);
else
    sol_actual = filtered_results(min(PICK, size(filtered_results,1)), xCols);
end
k1 = sol_actual(1);
k2 = sol_actual(2);
k3 = sol_actual(3);
[f, bpa] = minimizeFlxPin2brk(k1, k2, k3, [], USE_BRACKET2);  % [f: nBPA x 3], [bpa: full struct]

disp(array2table([k1, k2, k3], 'VariableNames', {'X0', 'X1', 'X2'}));

%Show results for pre- and post-optimization
fprintf('\nPerformance with no length offset and infinite stiffness:\n');
disp(array2table(a0, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}, ...
                    'RowNames', cellstr(labels')));

fprintf('\nPerformance with sol_actual (pick=%d):\n', PICK);
disp(array2table(f, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}, ...
                    'RowNames', cellstr(labels')));

fprintf('Mean baseline: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n\n',mean(a0,1,'omitnan'));
fprintf('Mean optimized: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n\n',mean(f,1,'omitnan'));

%% Save results
stamp = char(string(datetime('now'),'yyyyMMdd'));
resultFile = sprintf('minimizeFlxPin10_2brk_results_%s.mat', stamp);
save(resultFile, 'results_cv', 'all_candidates', 'results_sort', 'results_sort_actual', ...
     'filtered_results', 'xCols', 'a0', 'f', 'k1', 'k2', 'k3', 'PICK', ...
     'ALLBPA', 'NUMHOLD', 'POP', 'MAXGEN', 'USE_BRACKET2', 'SOLVER', 'labels', 'W');
fprintf('Results saved to %s\n', resultFile);

%% Plot torque curves, pre- and post-Optimized
if DO_PLOTS

load ForceStrainForFit.mat z
z(isnan(z)) = 0;
X = linspace(0,620,20); %Pressure for interpolation
X = X(2:20);
Y = linspace(0,1,30);   %Relative strain range for interpolation

figTpre = figure('Name','Torque, Pre-Optimized','Color','w');
figTpre.Position = [100 100 950 700];
tTpre = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');

titles = ["\bf 48.5 cm", "\bf 45.7 cm","\bf 47.9 cm", "\bf 40.6 cm", "\bf 41.7 cm"];
subtitles = ["\bf Pre-optimized","\bf Pre-optimized","\bf Pre-optimized","\bf Pre-optimized","\bf Optimized","\bf Optimized","\bf Optimized","\bf Optimized"];

for k = 1:numBPA
    ax = nexttile(k);
    j = ALLBPA(k);
    hold on
    % Pre-optimization: Mold calculation
    Yq = bpa(j).strain./((bpa(j).rest - bpa(j).Kmax)/bpa(j).rest);
    Xq = bpa(j).P.*ones(size(Yq));
    Vq = interp2(X, Y', z, Xq, Yq,'spline');
    Fold = Vq.*bpa(j).unitD;     %Old force calc
    Fq = hypot(Fold(:,1), Fold(:,2)); %Old force magnitude on the XY plane
    Mold = -bpa(j).mA.*Fq;       %Torque with old force calc
    Mold(Yq>1) = 0;
    Mold(bpa(j).strain < 0) = NaN;
    scatter(bpa(j).A_h, bpa(j).M_h, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', '#FFD700', 'DisplayName', 'Hybrid');
    scatter(bpa(j).Aexp, bpa(j).Mexp, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', '#0000FF', 'DisplayName', 'Measured');
    plot(bpa(j).Ak, Mold, '--', 'Color', '#000000', 'LineWidth', 2, 'DisplayName', 'Original');
    plot(bpa(j).Ak, bpa(j).M, '-.', 'Color', '#FA8775', 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');

    clear Vq Fold Fq Mold
    title(titles(k), 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');

    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    xlim([-120 20]); ylim([-25 0]);
end

ylabel(tTpre,'\bf Torque, N\cdotm','Interpreter','tex');
xlabel(tTpre,'\bf \theta_{k} , \circ','Interpreter','tex');


figTpost = figure('Name','Torque, Post-Optimized','Color','w');
figTpost.Position = [100 100 950 700];
tTpost = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');

for k = 1:numBPA
    ax = nexttile(k);
    j = ALLBPA(k);
    hold on

    scatter(bpa(j).Aexp, bpa(j).Mexp, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', '#0000FF', 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).M, '-.', 'Color', '#FA8775', 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, bpa(j).M_p(:,3), '-', 'Color', '#CD34B5', 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');

    title(titles(k), 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');

    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    xlim([-120 20]); ylim([-25 0]);
end

% Shared labels
ylabel(tTpost,'\bf Torque, N\cdotm','Interpreter','tex');
xlabel(tTpost,'\bf \theta_{k} , \circ','Interpreter','tex');

% Legends in 2nd and 4th tile
lg = legend(tTpre.Children(end-1));
lg.Location = 'best';
lg.FontSize = 8;

lg2 = legend(tTpost.Children(end-1));
lg2.Location = 'best';
lg2.FontSize = 8;

%% Plot muscle length, optimization and validation
figL = figure('Name','Muscle Length','Color','w');
figL.Position = [100 100 950 700];
tL = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');

for k = 1:numBPA
    ax = nexttile(k);
    j = ALLBPA(k);
    if bpa(j).ten > 0
        title(sprintf('\\bf l_0 = %0.1f cm, %0.0f mm tendon',bpa(j).rest*100, bpa(j).ten*10^3),'Interpreter','tex')
    else
        title(sprintf('\\bf l_0 = %0.1f cm',bpa(j).rest*100),'Interpreter','tex')
    end
    hold on

    % Calculate predicted
    Lm_p = bpa(j).Lmt_p - 2 * bpa(j).fitn - bpa(j).ten;
    Lm   = bpa(j).Lmt   - 2 * bpa(j).fitn - bpa(j).ten;
    % Measured
    scatter(bpa(j).A_h, bpa(j).Lm_h, 60, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', '#0000FF', 'DisplayName', 'Measured');
    %Old prediction
    plot(bpa(j).Ak, Lm, '-.', 'Color', '#FA8775', 'LineWidth', 2.5, 'DisplayName', 'Original prediction');
    % New prediction
    plot(bpa(j).Ak, Lm_p, '-', 'Color', '#CD34B5', 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');


    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

ylabel(tL,'\bf Muscle length, m','Interpreter','tex');
xlabel(tL,'\bf \theta_{k} , \circ','Interpreter','tex');

lg = legend(tL.Children(end-1));
lg.Location = 'best';
lg.FontSize = 8;

%% Plot moment arm, optimization and validation
figMA = figure('Name','Moment Arm','Color','w');
figMA.Position = [100 100 950 700];
tMA = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');


for k = 1:numBPA
    ax = nexttile(k);
    j = ALLBPA(k);
    if bpa(j).ten > 0
        title(sprintf('\\bf l_0 = %0.1f cm, %0.0f mm tendon',bpa(j).rest*100, bpa(j).ten*10^3),'Interpreter','tex')
    else
        title(sprintf('\\bf l_0 = %0.1f cm',bpa(j).rest*100),'Interpreter','tex')
    end
    hold on

    G_p = hypot(bpa(j).mA_p(:,1), bpa(j).mA_p(:,2));

    scatter(bpa(j).A_h, bpa(j).mA_h, 60, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', '#0000FF', 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).mA, '-.', 'Color', '#FA8775', 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, G_p, '-', 'Color', '#CD34B5', 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');

    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

ylabel(tMA,'\bf Moment arm, m','Interpreter','tex');
xlabel(tMA,'\bf \theta_{k} , \circ','Interpreter','tex');
legend(tMA.Children(end-1),'Location','best','FontSize',8);

%% Plot relative strain, optimization and validation
figS = figure('Name','Relative Strain','Color','w');
figS.Position = [100 100 950 700];
tS = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');

for k = 1:numBPA
    ax = nexttile(k);
    j = ALLBPA(k);
    if bpa(j).ten > 0
        title(sprintf('\\bf l_0 = %0.1f cm, %0.0f mm tendon',bpa(j).rest*100, bpa(j).ten*10^3),'Interpreter','tex')
    else
        title(sprintf('\\bf l_0 = %0.1f cm',bpa(j).rest*100),'Interpreter','tex')
    end
    hold on

    strain_h = (bpa(j).rest - bpa(j).Lm_h)/bpa(j).rest;
    kmax = (bpa(j).rest - bpa(j).Kmax)/bpa(j).rest;
    scatter(bpa(j).A_h, strain_h/kmax, 60, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', '#0000FF', 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).strain/kmax, '-.', 'Color', '#FA8775', 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, bpa(j).strain_p/kmax, '-', 'Color', '#CD34B5', 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');

    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

ylabel(tS,'\bf Relative strain','Interpreter','tex');
xlabel(tS,'\bf \theta_{k} , \circ','Interpreter','tex');
legend(tS.Children(end-1),'Location','best','FontSize',8);

end % DO_PLOTS

%% Helper functions
function ff = min1(x, trainIdx, kompare, useB2)
    if numel(x) == 3 && size(x,1) == 1
        % OK
    else
        error('min1: Input x must be a 1x3 vector');
    end
    Xi0 = x(1) / 100;
    Xi1 = 10^x(2);
    Xi2 = 10^x(3);

    try
        [f_all, ~] = minimizeFlxPin2brk(Xi0, Xi1, Xi2, trainIdx, useB2); % Nx3 matrix for training BPAs
        fnorm = f_all(trainIdx,:)./kompare(trainIdx,:);     %normalize results before taking the mean
        ff = mean(fnorm, 1, 'omitnan');              % Return 1x3: [mean RMSE, mean FVU, mean MaxResidual]
        if ~isnumeric(ff) || numel(ff) ~= 3
            ff = [Inf, Inf, Inf];  % Defensive return if shape is wrong
        end
    catch
        ff = [Inf, Inf, Inf];      % Defensive return if minimizeFlxPin2brk throws
    end
end

function fs = min1scalar(x, trainIdx, kompare, W, useB2)
    f3 = min1(x, trainIdx, kompare, useB2);
    if any(~isfinite(f3))
        fs = Inf;
    else
        fs = W * f3(:);
    end
end


%% --- Nonlinear constraint (unused; kept from original)
function [c, ceq] = nonlinc(X, baseline, trainIdx)
    % Inputs:
    %   X         = [Xi0_cm, log10(Xi1), log10(Xi2)]
    %   baseline  = baselineScores (size: 4x3)
    %   trainIdx  = which BPAs are being optimized
    try
        Xi0 = X(1) / 100;
        Xi1 = 10^X(2);
        Xi2 = 10^X(3);
        [f_all, ~] = minimizeFlxPin2brk(Xi0, Xi1, Xi2, trainIdx);
        mean_baseline = mean(baseline(trainIdx,:), 1, 'omitnan');  % 1x3
        mean_model = mean(f_all, 1, 'omitnan');                    % 1x3

        c = mean_model - mean_baseline;  % Element-wise (positive means violation)
        ceq = [];
    catch
        c = ones(4, 1) * 1e3;  % Large penalty
        ceq = [];
    end
end

% function [state, options, optchanged] = debugPop(options, state, flag)
%     ...
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
