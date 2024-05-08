
%Minimization scheme

clear; clc; close all

[a, b, ~] = minimizeFlx(0,Inf,Inf);         %Get current goodness of fit measures with no extra length and infinite bracket stiffness
%% Problem setup
x1 = optimvar('x1',1,'LowerBound',-0.007,'UpperBound',0.02,'Type','continuous');
kT = optimvar('kT',1,'LowerBound',100,'UpperBound',10^5,'Type','continuous');
kB = optimvar('kB',1,'LowerBound',100,'UpperBound',10^5,'Type','continuous');
lb = [-0.007, 1*10^2, 1*10^2];
ub = [0.02, 10^5, 10^5];
[f, g]   = fcn2optimexpr(@minimizeFlx,x1,kT,kB,'OutputSize',{[1,3],[1,3]},'ReuseEvaluation',true);
fun1 = f(1,1);
fun2 = f(1,2);
fun3 = f(1,3);
prob = optimproblem;
prob.Objective = fun2;
prob.Constraints.better1 = f(1) <= a(1);          %solution must be better than no optimization
prob.Constraints.better1 = f(2) <= a(2);          %solution must be better than no optimization
prob.Constraints.better2 = g(1) <= a(1);          %validation solution must be better than no optimization
prob.Constraints.better2 = g(2) <= a(2);          %validation solution must be better than no optimization
x0.x1 = 0.013;
x0.kT = 1000;
x0.kB = 1000;
% hybridopts = optimoptions('fmincon','OptimalityTolerance',1e-10);
% options = optimoptions(@particleswarm,'Display','iter','UseParallel',true,'PlotFcn','pswplotbestf',...
%                         'OutputFcn',@pswplotranges,'HybridFcn',{'fmincon',hybridopts});
options = 
show(prob)

% fun = @(x,y,z) minimizeFlx;
%% Solve 
[sol,fval,exitflag,output] = solve(prob,x0);
% [sol,fval,Pareto_front, Pareto_Fvals, exitflag,output] = GODLIKE(fun,lb,ub,[],...
%                                          'algorithms', {'DE';'GA';'ASA'},...
%                                          'display'   , 'plot',...
%                                          'popsize'   , 500);

[u,v,bpa] = minimizeFlx(sol.x1,sol.kT,sol.kB);           % Now pull bpa structures out       

%% Plot Optimized fit
figure
hold on
scatter(bpa(1).Aexp,bpa(1).Mexp,[],'filled','DisplayName','Experiment')
plot(bpa(1).Ak,bpa(1).M_p,'DisplayName','New predict')
plot(bpa(1).Ak,bpa(1).M,'DisplayName','Predict original')
hold off
title('Optimization')
xlabel('\theta_{k}, \circ')
ylable('Torque, N\cdotm')
legend

%% Plot validation
figure
hold on
scatter(bpa(2).Aexp,bpa(2).Mexp,[],'filled','DisplayName','Experiment')
plot(bpa(2).Ak,bpa(2).M_p,'DisplayName','New predict')
plot(bpa(2).Ak,bpa(2).M,'DisplayName','Predict original')
hold off
title('Validation')
xlabel('\theta_{k}, \circ')
ylable('Torque, N\cdotm')
legend

%% Helper functions
function ff = min1(x1,kT,kB)
ff = minimizeFlx(x1,kT,kB); %get GOF vector
end


function gg = min2(x1,kT,kB)
[~, gg] = minimizeFlx(x1,kT,kB); %get validation vector
end