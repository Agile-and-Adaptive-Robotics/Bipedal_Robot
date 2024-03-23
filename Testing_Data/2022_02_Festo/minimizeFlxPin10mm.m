
%% Optimization for fitting length
% Sections:
%     - Load data from MonoPamDataExplicit and measured data. Make initial fitting length guess. Pass on to minimizer.
%     - Pass on fitting length to another function to validate fit for other muscle length.
%     - Plot results of both muscles using new fitting length.
%     - Nested minimization function to use either Simplx, Pattern search, or Pareto search.
%     - Nested function that calculates SSE and RMSE of new fitting value.
%     - Nested function that calls muscle class and calculates torque.

clear; clc; close all

load KneeFlxPin_10mm_48cm.mat phiD Location Name CrossPoint Dia T rest tendon kmax fitting TorqueR
load Plot_KneeFlxPin_10mm_48cm.mat Angle Torque pres 
pres = mean(pres);                  %Make pressure a scalar value
y = Torque;                        %Make it just the data

x0 = 0.0254;                        %Initial fitting length guess
%Search optimization 48cm, fitting = 32.11 mm for best RMSE and R^2 
%Search optimization 46cm, fitting = 30.62 mm for best RMSE and R^2, 31.3 mm for max error. 

[point,value] = minimizeFlx10mm_nest(x0,y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, pres);

%% Validate with other muscle length in this configuration
A = [point, value];
%sort by RMSE, fvu, and max. residual, respectively
if length(point)>=3
    for i = 1:3
    B(:,:,i) = sortrows(A,i+1);
    end
    else
    B(:,:,1) = (A);
end
%After the fact, we know that first row of B is the best. Use it.
% B = B(1,:,1);
%pick first row of each
for i = 1:size(B,3)
 P(i,:) = B(1,2:4,i);
 par(i,:) = B(1,1,i);
end
% point = mean(point);

for i = 1:size(B,3)
RMSE1(i) = P(i,1);
Rsq1(i) = 1-P(i,2);
maxResid1(i) = P(i,3);
% fprintf('fitting length is %5d with adj. R^2=%5d\n',point,1-value)
fprintf('fitting length is %5d with RMSE=%5d, adj. R^2=%5d, and max residual=%5d\n',par(i),RMSE1(i), Rsq1(i), maxResid1(i))

[RMSE2(i), fvu2(i), maxResid2(i), Ang2, Tmeas2, Tnew2(:,i), Torig2, r2] = validateFlxPin10mm(par(i));
Rsq2(i) = 1 - fvu2(i);
fprintf('Validation returns RMSE=%5d, R^2=%5d, and max residual=%5d\n',RMSE2(i),Rsq2(i),maxResid2(i))
end


%% Plot the results
%Matlab hex color values:
c{1} = '#FFD700'; %gold
c{2} = '#FFB14E'; %orange
c{3} = '#FA8775'; %light orange
c{4} = '#EA5F94'; %pink
c{5} = '#CD34B5'; %magenta
c{6} = '#9D02D7'; %magenta 2
c{7} = '#0000FF'; %indigo
c{8} = '#000000'; %black
sz = 60;        %size of data points

opt = ["RMSE"; "R^2"; "Max Error"];

% %Call class for optimized muscle, use torque
for i = 1:size(B,3)
    Bifemsh_Pam_new(i) = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, par(i), pres);
    Tnew1(:,i) = Bifemsh_Pam_new(i).Torque(:,3);
end
Torig1 = TorqueR(:,3);
% 
%Plot Torque values against theoretical for optimized muscle
figure
ax1 = gca;
scM = scatter(Angle,Torque,sz,'filled','MarkerFaceColor',c{7},'DisplayName','Measured Torque');
hold on
PL1 = plot(phiD', Torig1,'Color',c{2},'Linewidth',2,'DisplayName','Calculated, original');
for i = 1:size(B,3)
    dispO = sprintf('Calculated, optimized');
    PL2 = plot(phiD', Tnew1(:,i),'Color',c{6-i},'Linewidth',2,'DisplayName',dispO);
end
hold off
title(sprintf('Optimized, l_{rest}=%4.3f m',rest),'FontSize',12,'FontWeight','Bold')
xlabel('Knee angle, \circ','FontSize',12,'Interpreter','tex')
ylabel('Torque, N \cdot m','FontSize',12,'Interpreter','tex')
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'southwest';

%Plot validation
figure
ax2 = gca;
hold on
scM2 = scatter(Ang2,Tmeas2,sz,'filled','MarkerFaceColor',c{7},'DisplayName','Measured Torque');
PL3 = plot(phiD', Torig2,'Color',c{2},'Linewidth',2,'DisplayName','Calculated, original');
for i = 1:size(B,3)
    dispV = sprintf('Calculated, optimized');
    PL4 = plot(phiD', Tnew2(:,i),'Color',c{6-i},'Linewidth',2,'DisplayName',dispV);
end
hold off
title(sprintf('Validation, l_{rest}=%4.3f m',r2),'FontSize',12,'FontWeight','Bold')
xlabel('Knee angle, \circ','FontSize',12,'Interpreter','tex')
ylabel('Torque, N \cdot m','FontSize',12,'Interpreter','tex')
ax2.FontSize = 12;
ax2.FontWeight = 'bold';
ax2.FontName = 'Arial';
ax2.YAxis.LineWidth = 2; ax2.YAxis.FontSize = 10;
ax2.XAxis.LineWidth = 2; ax2.XAxis.FontSize = 10;
lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'southwest';


%% Nested function for optimization
   function [X, FVAL] = minimizeFlx10mm_nest(x0,y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, pres)

        fun = @(x)minimize(y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres);
    
%     options = optimset('Display','iter','PlotFcns',@optimplotfval);     %fminsearch options
%     [X, FVAL, flag, put] = fminsearch(fun,x0,options);

%      [X, FVAL] = bnd_con_pen_nelder_mead( fun, x0, [0.0254 0.04], [], 100, [], 1, 1);
    
%     options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestx);        %pattern search options
%     [X, FVAL, flag, put] = patternsearch(fun,x0,[],[],[],[],0.02,.04,[],options);
    
    options = optimoptions('paretosearch','PlotFcn',{@psplotbestf, @psplotfuncount, @psplotparetof},'InitialPoints',x0);
    [X, FVAL, flag, put] = paretosearch(fun,1,[],[],[],[],0.0254,0.04,[],options);

%       options = optimoptions('particleswarm');
%       [X, FVAL, flag, put] = particleswarm(fun,1,0.02,0.04);
    
%     hybridopts = optimoptions('patternsearch','Display','iter');
%     options = optimoptions('gamultiobj','Display','iter','PlotFcn',{@gaplotpareto, @gaplotbestindiv, @gaplotscores},'HybridFcn',{@fgoalattain,hybridopts});
%     [X, FVAL, flag, put] = gamultiobj(fun,1,[],[],[],[],0.02,0.04,[],options);
    
    disp(X)
    disp(FVAL)

%% Nested function that calculates RMSE, FVU, and max error
        function f = minimize(y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres)

        ynew = nestedfun1(phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres);
          
        yresid = y-ynew;                     %residual error
        maxResid = max(abs(yresid));             %largest residual
        SSresid = sum(yresid.^2,'omitnan'); %Sum of squares of the residual
        SStot = sum((y-mean(y)).^2);        %total sum of squares
        n = (sum(~isnan(ynew)));             %number of data points
        RMSE = sqrt(SSresid/n);             % RMSE for function 1
        fvu = SSresid*(n-1)/(SStot*n);      %Fraction of Variance Unexplained (FVU),adjusted, Rsq = 1 - FVU  

        f1 = RMSE;
        f2 = fvu;                       
        f3 = maxResid;
       
%         f = f2;
        f = [f1, f2, f3];

  
            function val = nestedfun1(phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres)                               
           % Find the difference between experimental and calculated results
            fitt = x;
            Bifemsh_Pam_adj = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitt, pres); %Call the class to do the Torque calculation
            Tqz = Bifemsh_Pam_adj.Torque(:,3);    %Calculated Adjusted Torque
            
            %Get model torque values at each angle there is measurement data for
            F = griddedInterpolant(phiD',Tqz);
            val = F(Angle);
            end

    end
end