
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
weight = [1 3 0.5];
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

%% Plot torque curves, Optimized and validation 
load ForceStrainForFit.mat z
X = linspace(0,620,20); %Pressure for interpolation
X = X(2:20);
Y = linspace(0,1,30);   %Relative strain range for interpolation

str = ["Optimization"; "Validation"];
str2 = ["Torque"; "Muscle Length"; "Moment Arm"];
for i = 1:2
    %Find old old torque, then plot everything
    clear Yq Xq Vq Fold Fq Mold
    Yq = bpa(i).strain./((bpa(i).rest-bpa(i).Kmax)/bpa(i).rest);
    Xq = bpa(i).P;
    Vq = zeros(size(bpa(i).unitD,1),1);
    for j = 1:size(bpa(i).unitD, 1)
        if bpa(i).strain(j) >=-0.03 && Yq(j) <=1
            Vq(j) = interp2(X, Y, z, Xq, Yq(j));
%             Vq(j) = f10(Yq(j),Xq);
        elseif Yq(j)>1
            Vq(j) = 0;
        elseif bpa(i).strain(j) < -0.03
            Vq(j) = NaN;
        end
    end   
    Fold = Vq.*bpa(i).unitD;    %Force vector
    Fq = (Fold(:,1).^2+Fold(:,2).^2).^(1/2);
    Mold = -bpa(i).mA.*Fq;
    
    
    figure
    ax = gca;
    hold on
    scatter(bpa(i).Aexp,bpa(i).Mexp,[],'filled','DisplayName','Experiment')
    scatter(bpa(i).A_h,bpa(i).M_h,[],'filled','DisplayName','Hybrid')
    plot(bpa(i).Ak,bpa(i).M_p(:,3),'DisplayName','New predict')
    plot(bpa(i).Ak,bpa(i).M,'DisplayName','Predict original') %"original" is with updated BPA characterization
    plot(bpa(i).Ak,Mold,'DisplayName','Predict old') %"original" is with updated BPA characterization
    hold off
    title(sprintf('Torque, %s',str(i)))
    xlabel('\theta_{k}, \circ')
    ylabel('Torque, N\cdotm')
    legend
end

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