%%minimizeFlxPin10mm.m
%Minimization scheme

clear; clc; close all
profile on

[a, b, ~] = minimizeFlxPin(0,Inf,Inf);         %Get current goodness of fit measures with no extra length and infinite bracket stiffness
%% Problem setup

lb = [-0.03*100, 3, 3];
ub = [0.03*100, 9, 9];

%% Solve 
opts = optimoptions('gamultiobj', ...
        'UseParallel', true, ...
        'Display', 'iter', ...        %
        'InitialPopulationRange',[-.003*100 log10(700) log10(370); 0.015*100 5 4], ...
        'PopulationSize', 75, ...  %Originally 75
        'MaxGenerations', 600, ... %Originally 600
        'MutationFcn', {@mutationadaptfeasible}, ...
        'CrossoverFraction', 0.8, ...
        'CrossoverFcn', {@crossoverscattered}, ...
        'PlotFcn', {@gaplotpareto}, ...
        'FunctionTolerance', 1e-5);
goal = [0 0 0];
weight = [1/a(1) 1/a(2) 1/a(3)];
opts.HybridFcn = {@fgoalattain, goal, weight};
    
[x, fvals, exitflag, output, population, scores] = gamultiobj( ...
    @min1, 3, ...
    [], [], [], [], ...
    lb, ub, ...    %@(x) nonlinc(x), ...
    opts);

profile viewer
% [sol,fval,Pareto_front, Pareto_Fvals, exitflag,output] = GODLIKE(@min1,lb,ub,[],'NumObjectives',3,...
%                                          'algorithms', {'DE';'GA';'ASA';'PSO'},...
%                                          'display'   , 'plot',...
%                                          'popsize'   , 75 );
% x = Pareto_front;
% fvals = Pareto_Fvals;
%% Get validation GoF results from pareto front
val_Fvals = zeros(size(fvals));
for i = 1:length(x)
    val_Fvals(i,:) = min2([x(i,1),x(i,2),x(i,3)]);    %Get validation Fvals for all Pareto_front points
end

%% Sort results
ind = 1:length(x);  %Index to original Pareto_front and Pareto_Fvals
relate = vecnorm(fvals-val_Fvals,2,2);   %Find the distance between the optimization and validation solutions for the same input
results = [ind', x, fvals, val_Fvals, relate]; 
results_sort = sortrows(results,[11 8 9 10 5 6 7]); %Sort results first on distance between optimization and validation, then on validation columns, then on original Fvals columns.
results_sort_actual = [results_sort(:,1), results_sort(:,2)/100, 10.^results_sort(:,3), 10.^results_sort(:,4), results_sort(:,5:end)];

%% Pick ultimate solution
pick = 1; %Pick the best solution from the sorted results (should be 1)
sol_actual = results_sort_actual(pick, 2:4);  %Best solution                                   
[u,v,bpa] = minimizeFlxPin(sol_actual(1),sol_actual(2),sol_actual(3));           % Now pull bpa structures out       

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

%% Plot torque curves, Optimized and validation 
load ForceStrainForFit.mat z
X = linspace(0,620,20); %Pressure for interpolation
X = X(2:20);
Y = linspace(0,1,30);   %Relative strain range for interpolation

figT = figure('Name','Torque','Color','w');
figT.Position = [100 100 950 700];
tT = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

titles = ["\bf 48.0 cm", "\bf 45.7 cm","\bf 48.0 cm", "\bf 45.7 cm"];
subtitles = ["\bf Pre-optimized","\bf Pre-optimized","\bf Optimized","\bf Validation"];

for j = 1:2
    ax = nexttile(j);
    hold on

    % Pre-optimization: Mold calculation
    Yq = bpa(j).strain./((bpa(j).rest - bpa(j).Kmax)/bpa(j).rest);
    Xq = bpa(j).P;
    Vq = zeros(size(bpa(j).unitD,1),1);
    for k = 1:size(bpa(j).unitD,1)
        if bpa(j).strain(k) >= -0.03 && Yq(k) <= 1
            Vq(k) = interp2(X, Y, z, Xq, Yq(k));
        elseif Yq(k) > 1
            Vq(k) = 0;
        elseif bpa(j).strain(k) < -0.03
            Vq(k) = NaN;
        end
    end
    Fold = Vq.*bpa(j).unitD;
    Fq = sqrt(Fold(:,1).^2 + Fold(:,2).^2);
    Mold = -bpa(j).mA.*Fq;

    scatter(bpa(j).A_h, bpa(j).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{1}, 'DisplayName', 'Hybrid');
    scatter(bpa(j).Aexp, bpa(j).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, Mold, '--', 'Color', c{8}, 'LineWidth', 2, 'DisplayName', 'Original');
    plot(bpa(j).Ak, bpa(j).M, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');

    title(titles(j), 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    
    ylabel('\bf Torque, N \cdot m', 'Interpreter', 'tex', ...
            'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    xlabel('\bf \theta_{k}, \circ', 'Interpreter', 'tex', ...
            'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');

    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    subtitle(subtitles(j), 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold');
    xlim([-120 20]); ylim([-25 0]);
end

for j = 1:2
    ax = nexttile(2 + j);
    hold on
    
    scatter(bpa(j).Aexp, bpa(j).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).M, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, bpa(j).M_p(:,3), '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');
    
    title(titles(j), 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    
    ylabel('\bf Torque, N \cdot m', 'Interpreter', 'tex', ...
            'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    xlabel('\bf \theta_{k}, \circ', 'Interpreter', 'tex', ...
            'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    subtitle(subtitles(j+2), 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold');
    xlim([-120 20]); ylim([-25 0]);
end

% Shared labels
% ylabel(tT,'\bf Torque, N\cdotm','Interpreter','tex');
% xlabel(tT,'\bf \theta_{k} , \circ','Interpreter','tex');

% (A)-(D) annotations
for j = 1:4
    annotation(gcf, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end

% Column titles
% annotation(gcf, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% annotation(gcf, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

% Legends in 2nd and 4th tile
lg = legend(tT.Children(end-1));
lg.Location = 'best';
lg.FontSize = 8;

lg2 = legend(tT.Children(end-3));
lg2.Location = 'best';
lg2.FontSize = 8;

%% Plot muscle length, optimization and validation
figL = figure('Name','Muscle Length','Color','w');
figL.Position = [100 100 950 700];
tL = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

for j = 1:2
    ax = nexttile(j);
    hold on
    
    % Predicted
    Lm_p = bpa(j).Lmt_p - 2 * bpa(j).fitn - bpa(j).ten;
    Lm   = bpa(j).Lmt   - 2 * bpa(j).fitn - bpa(j).ten;
    
    scatter(bpa(j).A_h, bpa(j).Lm_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, Lm, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Original prediction');

    if j == 1
        ylabel('\bf Pre-optimization', 'Interpreter', 'tex', ...
            'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    end
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

for j = 1:2
    ax = nexttile(2 + j);
    hold on
    
    % Predicted
    Lm_p = bpa(j).Lmt_p - 2 * bpa(j).fitn - bpa(j).ten;
    Lm   = bpa(j).Lmt   - 2 * bpa(j).fitn - bpa(j).ten;
    
    scatter(bpa(j).A_h, bpa(j).Lm_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, Lm, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Original prediction');
    plot(bpa(j).Ak, Lm_p, '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');

    if (2 + j) == 3
        ylabel('\bf Optimized', 'Interpreter', 'tex', ...
            'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    end
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

ylabel(tL,'\bf Muscle length, m','Interpreter','tex');
xlabel(tL,'\bf \theta_{k} , \circ','Interpreter','tex');


for j = 1:4
    annotation(figL, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end

annotation(figL, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
annotation(figL, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

lg = legend(tL.Children(end-1));
lg.Location = 'best';
lg.FontSize = 8;

lg2 = legend(tL.Children(end-3));
lg2.Location = 'best';
lg2.FontSize = 8;

%% Plot moment arm, optimization and validation
figMA = figure('Name','Moment Arm','Color','w');
figMA.Position = [100 100 950 700];
tMA = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

for j = 1:2
    ax = nexttile(j);
    hold on
    scatter(bpa(j).A_h, bpa(j).mA_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).mA, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');

    if j == 1
        ylabel('\bf Pre-optimization', 'Interpreter', 'tex', ...
            'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    end
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

for j = 1:2
    ax = nexttile(2 + j);
    hold on
    
    G_p = hypot(bpa(j).mA_p(:,1), bpa(j).mA_p(:,2));
    
    scatter(bpa(j).A_h, bpa(j).mA_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).mA, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, G_p, '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');

    if (2 + j) == 3
        ylabel('\bf Optimized', 'Interpreter', 'tex', ...
            'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    end
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

ylabel(tMA,'\bf Moment arm, m','Interpreter','tex');
xlabel(tMA,'\bf \theta_{k} , \circ','Interpreter','tex');
for j = 1:4
    annotation(figMA, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end
annotation(figMA, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
annotation(figMA, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
legend(tMA.Children(end-1),'Location','best','FontSize',8);
legend(tMA.Children(end-3),'Location','best','FontSize',8);

%% Plot relative strain, optimization and validation
figS = figure('Name','Relative Strain','Color','w');
figS.Position = [100 100 950 700];
tS = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

for j = 1:2
    ax = nexttile(j);
    hold on
    
    strain_h = (bpa(j).rest - bpa(j).Lm_h)/bpa(j).rest;
    kmax = (bpa(j).rest - bpa(j).Kmax)/bpa(j).rest;
    
    scatter(bpa(j).A_h, strain_h/kmax, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).strain/kmax, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Original prediction');

    if j == 1
        ylabel('\bf Pre-optimization', 'Interpreter', 'tex', ...
            'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    end
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

for j = 1:2
    ax = nexttile(2 + j);
    hold on
    
    strain_h = (bpa(j).rest - bpa(j).Lm_h)/bpa(j).rest;
    kmax = (bpa(j).rest - bpa(j).Kmax)/bpa(j).rest;
    scatter(bpa(j).A_h, strain_h/kmax, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).strain/kmax, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, bpa(j).strain_p/kmax, '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');

    if (2 + j) == 3
        ylabel('\bf Optimized', 'Interpreter', 'tex', ...
            'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    end
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

ylabel(tS,'\bf Relative strain','Interpreter','tex');
xlabel(tS,'\bf \theta_{k} , \circ','Interpreter','tex');
for j = 1:4
    annotation(figS, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end
annotation(figS, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
annotation(figS, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
legend(tS.Children(end-1),'Location','best','FontSize',8);
legend(tS.Children(end-3),'Location','best','FontSize',8);


%% Helper functions
function ff = min1(x)
ff = minimizeFlxPin(x(:,1)./100,10.^x(:,2),10.^x(:,3)); %get GOF vector
end


function gg = min2(x)
[~, gg] = minimizeFlxPin(x(:,1)./100,10.^x(:,2),10.^x(:,3)); %get validation vector
end

function [c, ceq] = nonlinc(X)
    % Evaluate the model
    [f, g] = minimizeFlxPin(X(:,1)/100, 10^X(:,2), 10^X(:,3));
    
    % Baseline a, b must be pre-loaded or set globally
    persistent a_local b_local
    if isempty(a_local) || isempty(b_local)
        [a_local, b_local, ~] = minimizeFlxPin(0, Inf, Inf);
    end
    
    % Constraints: f < a and g < b
    c = [f - a_local; g - b_local];
    ceq = [];
end