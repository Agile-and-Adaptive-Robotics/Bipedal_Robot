
%Minimization scheme

clear; clc; close all

[a, b, ~] = minimizeFlx(0,Inf,Inf);         %Get current goodness of fit measures with no extra length and infinite bracket stiffness
%% Problem setup
% x1 = optimvar('x1',1,'LowerBound',-0.007*10^3,'UpperBound',0.02*10^3,'Type','continuous');
% kT = optimvar('kT',1,'LowerBound',100*10^-3,'UpperBound',10^5*10^-3,'Type','continuous');
% kB = optimvar('kB',1,'LowerBound',100*10^-3,'UpperBound',10^5*10^-3,'Type','continuous');
lb = [-0.01*10^3, 10^3*10^-3, 10^3*10^-3];
ub = [0.03*10^3, 10^8*10^-3, 10^8*10^-3];
% f   = fcn2optimexpr(@min1,x1*10^-3,kT*10^3,kB*10^3,'OutputSize',[1,3],'ReuseEvaluation',true);
% g   = fcn2optimexpr(@min2,x1*10^-3,kT*10^3,kB*10^3,'OutputSize',[1,3],'ReuseEvaluation',true);
% fun1 = f(1,1);
% fun2 = f(1,2);
% fun3 = f(1,3);
% prob = optimproblem;
% prob.Objective = f;
% prob.Constraints.better1 = f <= a;          %solution must be better than no optimization
% prob.Constraints.better2 = g <= b;          %validation solution must be better than no optimization
% x0.x1 = 0.013*10^3;
% x0.kT = 1000*10^-3;
% x0.kB = 1000*10^-3;
% hybridopts = optimoptions('fmincon','OptimalityTolerance',1e-10);
% options = optimoptions(@particleswarm,'Display','iter','UseParallel',true,'PlotFcn','pswplotbestf',...
%                         'OutputFcn',@pswplotranges,'HybridFcn',{'fmincon',hybridopts});
% options = optimoptions(@paretosearch,'Display','iter','PlotFcn',{'psplotparetof','psplotparetox'},'InitialPoints',cell2mat(struct2cell(x0))); 
% show(prob)
% f = @(x,y,z)min1(x,y,z);
% g = @(x,y,z)min2(x,y,z);
% fun = @(x,y,z) minimizeFlx;

%% Solve 
% [sol,fval,exitflag,output] = solve(prob,Solver="paretosearch");
% [sol,fval,exitflag,output] = paretosearch(f,3,[],[],[],[],lb,ub,@nonlinc,options);
[sol,fval,Pareto_front, Pareto_Fvals, exitflag,output] = GODLIKE(@min1,lb,ub,[],'NumObjectives',3,...
                                         'algorithms', {'DE';'GA';'PSO';'ASA'},...
                                         'display'   , 'plot',...
                                         'popsize'   , 300);

sol_actual = [sol(1)*10^-3, sol(2)*10^3, sol(3)*10^3];                                     
[u,v,bpa] = minimizeFlx(sol_actual(1),sol_actual(2),sol_actual(3));           % Now pull bpa structures out       

%% Plot torque curves, Optimized and validation 
load ForceStrainForFit.mat z
X = linspace(0,620,20); %Pressure for interpolation
X = X(2:20);
Y = linspace(0,1,30);   %Relative strain range for interpolation

str = ["Optimization"; "Validation"];
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
    title(str(i))
    xlabel('\theta_{k}, \circ')
    ylabel('Torque, N\cdotm')
    legend
end

%% Plot muscle length, optimization and validation
for i = 1:2
    figure
    ax = gca;
    Lm = bpa(1).Lmt-2*bpa(i).fitn-bpa(i).ten;      %Original predicted muscle length
    Lm_p = bpa(1).Lmt_p-2*bpa(i).fitn-bpa(i).ten;    %Optimized muscle length (Lmt_p uses sol_actual(1)
    hold on
    scatter(bpa(i).A_h,bpa(i).Lm_h,[],'filled','DisplayName','Measured')
    plot(bpa(i).Ak,Lm_p,'DisplayName','New predict')
    plot(bpa(i).Ak,Lm,'DisplayName','Predict original')
    hold off
    title(str(i))
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
    title(str(i))
    xlabel('\theta_{k}, \circ')
    ylabel('Length, m')
    legend
end

%% Helper functions
function ff = min1(x)
ff = minimizeFlx(x(:,1)*10^-3,x(:,2)*10^3,x(:,3)*10^3); %get GOF vector
end


function gg = min2(x)
[~, gg] = minimizeFlx(x(:,1)*10^-3,x(:,2)*10^3,x(:,3)*10^3); %get validation vector
end

function [c, ceq] = nonlinc(f,g)
c = [f-a; g-b];
ceq = [];
end