%%minimizeFlxPin10mm.m
%Minimization scheme

clear; clc; close all
% profile on

%% Cross-validation setup
allBPA = [2, 3, 4, 5];           % Use if all data are valid 
% allBPA = [2, 3, 4, 5];              % Use if old data do not hold up
labels = ["48cm", "46cm", "47cm", "40cm-tendon", "42cm"];
validLabels = labels(allBPA);
numBPA = numel(allBPA);

results_cv = cell(1, numBPA);  % Will store RMSE, FVU, Max Resid for BPA(s) optimized
scores_cv = zeros(numBPA, 3);  % Will store RMSE, FVU, Max Resid for BPA(s) held-out for validation

%% Calculate baseline
[a0, bpa0] = minimizeFlxPinX3(0,Inf,Inf,Inf);         %Get current goodness of fit measures with no extra length and infinite bracket stiffness
baselineScores = a0;
fprintf('\nPerformance with no length offset and infinite stiffness:\n');
disp(array2table(a0, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}, ...
                    'RowNames', cellstr(labels')));
fprintf('Mean baseline training: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n\n',mean(baselineScores(:,1)),mean(baselineScores(:,2)),mean(baselineScores(:,3)));


%% Problem bounds
lb = [0*100, log10(5e3), log10(5e3), log10(5e3)];
ub = [0.020*100, log10(5e7), log10(5e7), log10(5e7)];

%% Solver
numHold = 2;                        %Number of BPAs held out for validation
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
            'Display', 'iter', ...        %'InitialPopulationRange',[0.002, 3.9, 3.9; 0.012*100 4.3 4.3], ...
            'PopulationSize', 150, ...  %Originally 75
            'MaxGenerations', 600, ... %Originally 600
            'MutationFcn', {@mutationadaptfeasible}, ...
            'CrossoverFraction', 0.8, ...
            'CrossoverFcn', {@crossoverscattered}, ...
            'PlotFcn', {@gaplotpareto3D_simple}, ...
            'FunctionTolerance', 4e-3);
    % hybridOpts = optimoptions('fgoalattain', ...
    %     'Display','iter');
    % 
    % opts.HybridFcn = {@fgoalattain, hybridOpts};
        
    % Run optimization
     [x, fvals,exitflag,output,population,scores] = gamultiobj(@(X) min1(X, trainIdx, a0), 4, [], [], [], [], ...
                                                    lb, ub, ... % @(x) nonlcon2(x), ...
                                                    opts);

% If I want to use GODLIKE instead:  
%     trainIdx_fixed = trainIdx;  % capture loop var for closure
%     min1_wrapper = @(x) min1(x, trainIdx_fixed, a0);
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
holdCols  = numel(rankCol) + (1:numHold);
xCols     = numel(rankCol) + numel(holdCols) + (1:4);
trainCols = numel(rankCol) + numel(holdCols) + numel(xCols) + (1:3);
valCols   = numel(rankCol) + numel(holdCols) + numel(xCols) + numel(trainCols) + (1:3);
distCol   = numel(rankCol) + numel(holdCols) + numel(xCols) + numel(trainCols) + numel(valCols) + 1;

%% Sort by validation distance first, then validation metrics
results = all_candidates;  % [fold, Xi0, Xi1, Xi2, Xi3, train (3), val (3), dist]
results_sort = sortrows(results, [distCol valCols(1:2) trainCols(1:2)]);  % sort by distance, then val, then train (excluding Max. Residual)

%% De-normalize to get physical parameters
x_actual = [results_sort(:,xCols(1))/100, 10.^results_sort(:,xCols(2)), 10.^results_sort(:,xCols(3)), 10.^results_sort(:,xCols(4))];
results_sort_actual = [results_sort(:,rankCol), results_sort(:,holdCols), x_actual, results_sort(:,trainCols), results_sort(:,valCols), results_sort(:,distCol)];

%% --- Filter Pareto candidates against baseline on BPAs 1, 2, 3 & 4 ---
N = size(results_sort_actual, 1);
keep = false(N,1);

for ii = 1:N
    % extract decision variables
    Xi0 = results_sort_actual(ii,xCols(1));
    Xi1 = results_sort_actual(ii,xCols(2));
    Xi2 = results_sort_actual(ii,xCols(3));
    Xi3 = results_sort_actual(ii,xCols(4));

    % re-evaluate on all 4 BPAs
    f_all = minimizeFlxPinX3(Xi0, Xi1, Xi2, Xi3);   % returns 4×3 [RMSE, FVU, MaxResidual] Assuming number of BPAs = 4.

    % compare RMSE & FVU for BPAs 1, 2, 3 & 4 to baselineScores
    % pass1 = all( f_all(1,1:3) <= baselineScores(1,1:3) );
    pass2 = all( f_all(2,1:3) <= baselineScores(2,1:3) );
    pass3 = all( f_all(3,1:3) <= baselineScores(3,1:3) );
    pass4 = all( f_all(4,1:3) <= baselineScores(4,1:3) );
    pass5 = all( f_all(5,1:3) <= baselineScores(5,1:3) );

    keep(ii) = pass2 && pass3 && pass4 && pass5;
end

% keep only the rows that passed all three checks
filtered_results = results_sort_actual(keep, :);
fprintf('Filtered %d → %d candidates.\n', N, sum(keep));


%% Pick best solution (later, flexible)
 
pick = 2;
sol_actual = filtered_results(pick, xCols);
k1 = sol_actual(1);
k2 = sol_actual(2);
k3 = sol_actual(3);
k4 = sol_actual(4);
% k1 = 0.01;
% k2 = 5e4;
% k3 = 3e3;
% k4 = 5e6;
[f, bpa] = minimizeFlxPinX3(k1, k2, k3, k4);  % [f: 4x3], [bpa: full struct]

disp(array2table([k1, k2, k3, k4], 'VariableNames', {'X0', 'X1', 'X2', 'X3'}));

%Show results for pre- and post-optimization
fprintf('\nPerformance with no length offset and infinite stiffness:\n');
disp(array2table(a0, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}, ...
                    'RowNames', cellstr(labels')));

fprintf('\nPerformance with sol_actual:\n');
disp(array2table(f, 'VariableNames', {'RMSE', 'FVU', 'MaxResidual'}, ...
                    'RowNames', cellstr(labels')));

fprintf('Mean baseline: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n\n',a0(1),a0(2),a0(3));
fprintf('Mean optimized: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n\n',f(1),f(2),f(3));
   

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

tileLabels = {'(A)', '(B)', '(C)', '(D)'};
% Annotation positions [x, y] in normalized figure units
xAnn = [0, 0.48, 0, 0.48];
yAnn = [0.94, 0.94, 0.45, 0.45];
sz = 60;

%% Plot torque curves, pre- and post-Optimized
load ForceStrainForFit.mat z
z(isnan(z)) = 0;
X = linspace(0,620,20); %Pressure for interpolation
X = X(2:20);
Y = linspace(0,1,30);   %Relative strain range for interpolation

figTpre = figure('Name','Torque, Pre-Optimized','Color','w');
figTpre.Position = [100 100 950 700];
tTpre = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');

% titles = ["\bf 48.5 cm", "\bf 45.7 cm","\bf 47.9 cm", "\bf 40.6 cm", "\bf 41.7 cm"];
titles = validLabels;
subtitles = ["\bf Pre-optimized","\bf Pre-optimized","\bf Pre-optimized","\bf Pre-optimized","\bf Optimized","\bf Optimized","\bf Optimized","\bf Optimized"];

for k = 1:numBPA
    ax = nexttile(k);
    j = allBPA(k);
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
    scatter(bpa(j).A_h, bpa(j).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{1}, 'DisplayName', 'Hybrid');
    scatter(bpa(j).Aexp, bpa(j).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, Mold, '--', 'Color', c{8}, 'LineWidth', 2, 'DisplayName', 'Original');
    plot(bpa(j).Ak, bpa(j).M, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');

    clear Vq Fold Fq Mold
    title(titles(k), 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    
    % ylabel('\bf Torque, N \cdot m', 'Interpreter', 'tex', ...
    %         'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    % xlabel('\bf \theta_{k}, \circ', 'Interpreter', 'tex', ...
    %         'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');

    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    % subtitle(subtitles(k), 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold');
    xlim([-120 20]); ylim([-25 0]);
end

ylabel(tTpre,'\bf Torque, N\cdotm','Interpreter','tex');
xlabel(tTpre,'\bf \theta_{k} , \circ','Interpreter','tex');


figTpost = figure('Name','Torque, Post-Optimized','Color','w');
figTpost.Position = [100 100 950 700];
tTpost = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');

for k = 1:numBPA
    ax = nexttile(k);
    j = allBPA(k);
    hold on
    
    scatter(bpa(j).Aexp, bpa(j).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).M, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, bpa(j).M_p(:,3), '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');
    
    title(titles(k), 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    
    % ylabel('\bf Torque, N \cdot m', 'Interpreter', 'tex', ...
    %         'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    % xlabel('\bf \theta_{k}, \circ', 'Interpreter', 'tex', ...
    %         'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    % subtitle(subtitles(k), 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold');
    xlim([-120 20]); ylim([-25 0]);
end

% Shared labels
ylabel(tTpost,'\bf Torque, N\cdotm','Interpreter','tex');
xlabel(tTpost,'\bf \theta_{k} , \circ','Interpreter','tex');

% (A)-(D) annotations
% for j = 1:4
%     annotation(gcf, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
%         'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% end

% Column titles
% annotation(gcf, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% annotation(gcf, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

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
    j = allBPA(k);
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
    scatter(bpa(j).A_h, bpa(j).Lm_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    %Old prediction
    plot(bpa(j).Ak, Lm, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Original prediction');
    % New prediction
    plot(bpa(j).Ak, Lm_p, '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');


    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

ylabel(tL,'\bf Muscle length, m','Interpreter','tex');
xlabel(tL,'\bf \theta_{k} , \circ','Interpreter','tex');


% for j = 1:4
%     annotation(figL, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
%         'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% end
% 
% annotation(figL, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% annotation(figL, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

lg = legend(tL.Children(end-1));
lg.Location = 'best';
lg.FontSize = 8;

% lg2 = legend(tL.Children(end-(numBPA-1)));
% lg2.Location = 'best';
% lg2.FontSize = 8;

%% Plot moment arm, optimization and validation
figMA = figure('Name','Moment Arm','Color','w');
figMA.Position = [100 100 950 700];
tMA = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');


for k = 1:numBPA
    ax = nexttile(k);
    j = allBPA(k);
    if bpa(j).ten > 0
        title(sprintf('\\bf l_0 = %0.1f cm, %0.0f mm tendon',bpa(j).rest*100, bpa(j).ten*10^3),'Interpreter','tex')
    else
        title(sprintf('\\bf l_0 = %0.1f cm',bpa(j).rest*100),'Interpreter','tex')
    end
    hold on
    
    G_p = hypot(bpa(j).mA_p(:,1), bpa(j).mA_p(:,2));
    
    scatter(bpa(j).A_h, bpa(j).mA_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).mA, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, G_p, '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');

    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

ylabel(tMA,'\bf Moment arm, m','Interpreter','tex');
xlabel(tMA,'\bf \theta_{k} , \circ','Interpreter','tex');
% for j = 1:4
%     annotation(figMA, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
%         'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% end
% annotation(figMA, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% annotation(figMA, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
legend(tMA.Children(end-1),'Location','best','FontSize',8);
% legend(tMA.Children(end-(numBPA-1)),'Location','best','FontSize',8);

%% Plot relative strain, optimization and validation
figS = figure('Name','Relative Strain','Color','w');
figS.Position = [100 100 950 700];
tS = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');

for k = 1:numBPA
    ax = nexttile(k);
    j = allBPA(k);
    if bpa(j).ten > 0
        title(sprintf('\\bf l_0 = %0.1f cm, %0.0f mm tendon',bpa(j).rest*100, bpa(j).ten*10^3),'Interpreter','tex')
    else
        title(sprintf('\\bf l_0 = %0.1f cm',bpa(j).rest*100),'Interpreter','tex')
    end
    hold on
    
    strain_h = (bpa(j).rest - bpa(j).Lm_h)/bpa(j).rest;
    kmax = (bpa(j).rest - bpa(j).Kmax)/bpa(j).rest;
    scatter(bpa(j).A_h, strain_h/kmax, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).strain/kmax, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, bpa(j).strain_p/kmax, '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');

    % if j == 1
    %     ylabel('\bf Optimized', 'Interpreter', 'tex', ...
    %         'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    % end
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

ylabel(tS,'\bf Relative strain','Interpreter','tex');
xlabel(tS,'\bf \theta_{k} , \circ','Interpreter','tex');
% for j = 1:4
%     annotation(figS, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
%         'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% end
% annotation(figS, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% annotation(figS, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
legend(tS.Children(end-1),'Location','best','FontSize',8);
% legend(tS.Children(end-(numBPA+1)),'Location','best','FontSize',8);


%% Helper functions
function ff = min1(x, trainIdx, kompare)
    if numel(x) == 4 && size(x,1) == 1
        % OK
    else
        error('min1: Input x must be a 1x3 vector');
    end
    Xi0 = x(1) / 100;
    Xi1 = 10^x(2);
    Xi2 = 10^x(3);
    Xi3 = 10^x(4);

    try
        [f_all, ~] = minimizeFlxPinX3(Xi0, Xi1, Xi2, Xi3, trainIdx); % Nx3 matrix (e.g., 3x3 if 3 training BPAs)
        fnorm = f_all(trainIdx,:)./kompare(trainIdx,:);     %normalize results before taking the mean
        ff = mean(fnorm, 1, 'omitnan');              % Return 1x3: [mean RMSE, mean FVU, mean MaxResidual]
        if ~isnumeric(ff) || numel(ff) ~= 3
            ff = [Inf, Inf, Inf];  % Defensive return if shape is wrong
        end
    catch
        ff = [Inf, Inf, Inf];      % Defensive return if minimizeFlxPin throws
    end
end


%% --- Nonlinear constraint
function [c, ceq] = nonlinc(X, baseline, trainIdx)
    % Inputs:
    %   X         = [Xi0_cm, log10(Xi1), log10(Xi2)]
    %   baseline  = baselineScores (size: 4x3)
    %   trainIdx  = which BPAs are being optimized
    try
        % Convert to physical parameters
        Xi0 = X(1) / 100;
        Xi1 = 10^X(2);
        Xi2 = 10^X(3);
        Xi3 = 10^X(4);
        % Call model with current X on the training BPAs
        [f_all, ~] = minimizeFlxPinX3(Xi0, Xi1, Xi2, Xi3, trainIdx);
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