
%Minimization scheme

clear; clc; close all

[a, b, ~] = minimizeFlx(0,Inf,Inf);
%% Problem setup
x1 = optimvar('x1',1,'LowerBound',-0.007,'UpperBound',0.02,"Type","continuous");
kT = optimvar('kT',1,'LowerBound',10000,'UpperBound',Inf,"Type","continuous");
kB = optimvar('kB',1,'LowerBound',10000,'UpperBound',Inf,"Type","continuous");
[f, g, bpa]   = fcn2optimexpr(@minimizeFlx,x1,kT,kB,'OutputSize',{[1,3],[1,3],[1,2]});
prob = optimproblem;
prob.Objective = f;

x0.x1 = 0.013;
x0.kT = 300000;
x0.kB = 300000;
show(prob)

[sol,fval,exitflag,output] = solve(prob,x0);


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