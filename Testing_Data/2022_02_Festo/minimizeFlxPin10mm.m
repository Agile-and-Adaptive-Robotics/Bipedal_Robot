%%minimizeFlxPin10mm.m
%Minimization scheme

clear; clc; close all

[a, b, ~] = minimizeFlxPin(0,Inf,Inf);         %Get current goodness of fit measures with no extra length and infinite bracket stiffness
%% Problem setup

lb = [-0.03*100, 3, 3];
ub = [0.03*100, 9, 9];

%% Solve 
opts = optimoptions('gamultiobj', ...
        'UseParallel', true, ...
        'Display', 'iter', ...
        'InitialPopulationRange',[-.003*100 log10(700) log10(370); 0.015*100 5 4], ...
        'PopulationSize', 75, ...
        'MaxGenerations', 600, ...
        'MutationFcn', {@mutationadaptfeasible}, ...
        'CrossoverFraction', 0.8, ...
        'CrossoverFcn', {@crossoverscattered}, ...
    'PlotFcn', {@gaplotpareto}, ...
    'FunctionTolerance', 1e-5);
goal = [0 0 0];
weight = [1 3 0.25];
opts.HybridFcn = {@fgoalattain, goal, weight};
    
[x, fvals, exitflag, output, population, scores] = gamultiobj( ...
    @min1, 3, ...
    [], [], [], [], ...
    lb, ub, ...
    @(x) nonlinc(x), ...
    opts);

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
xAnn = [0.002, 0.48, 0.002, 0.48];
yAnn = [0.89, 0.89, 0.41, 0.41];
sz = 60;

%% Plot torque curves, Optimized and validation 
load ForceStrainForFit.mat z
X = linspace(0,620,20); %Pressure for interpolation
X = X(2:20);
Y = linspace(0,1,30);   %Relative strain range for interpolation

figT = figure('Name','Torque','Color','w');
figT.Position = [100 100 950 700];
tT = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

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

    % Add y-label only for tile 1
    if (j) == 1
        ylabel('\bf Pre-optimization', 'Interpreter', 'tex', ...
            'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    end
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    xlim([-120 20]); ylim([-35 0]);
end

for j = 1:2
    ax = nexttile(2 + j);
    hold on
    scatter(bpa(j).A_h, bpa(j).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).M, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, bpa(j).M_p(:,3), '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');
    
    % Add y-label only for tile 3
    if (2 + j) == 3
        ylabel('\bf Optimized', 'Interpreter', 'tex', ...
            'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    end
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    xlim([-120 20]); ylim([-35 0]);
end

% Shared labels
ylabel(tT,'\bf Torque, N\cdotm','Interpreter','tex');
xlabel(tT,'\bf \theta_{k} , \circ','Interpreter','tex');

% (A)-(D) annotations
for j = 1:4
    annotation(gcf, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end

% Column titles
annotation(gcf, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
annotation(gcf, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

% Legends in 2nd and 4th tile
lg = legend(tT.Children(2));
lg.Location = 'northeast';
lg.FontSize = 8;

lg2 = legend(tT.Children(4));
lg2.Location = 'northeast';
lg2.FontSize = 8;

%% Plot muscle length, optimization and validation
for i = 1:2
    figure
    ax = gca;
    Lm = bpa(i).Lmt-2*bpa(i).fitn-bpa(i).ten;      %Original predicted muscle length
    Lm_p = bpa(i).Lmt_p-2*bpa(i).fitn-bpa(i).ten;    %Optimized muscle length (Lmt_p uses sol_actual(1)
    hold on
    scatter(bpa(i).A_h,bpa(i).Lm_h,[],'filled','DisplayName','Measured')
    plot(bpa(i).Ak,Lm_p,'DisplayName','New predict')
    plot(bpa(i).Ak,Lm,'DisplayName','Predict original')
    hold off
    title(sprintf('Muscle length, %s',str(i)))
    xlabel('\theta_{k}, \circ')
    ylabel('Length, m')
    legend
end

%% Plot moment arm, optimization and validation
for i = 1:2
    figure
    ax = gca;
    G_p = (bpa(i).mA_p(:,1).^2+bpa(i).mA_p(:,2).^2).^(1/2);      %z-axis moment arm for optimized
    hold on
    scatter(bpa(i).A_h,bpa(i).mA_h,[],'filled','DisplayName','Measured')
    plot(bpa(i).Ak,G_p,'DisplayName','New predict')
    plot(bpa(i).Ak,bpa(i).mA,'DisplayName','Predict original')
    hold off
    title(sprintf('Moment arm, %s',str(i)))
    xlabel('\theta_{k}, \circ')
    ylabel('Length, m')
    legend
end

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