
%% Optimization for fitting length
% Sections:
%     - Load data from MonoPamDataExplicit_compare and measured data. Make initial fitting length guess. Pass on to minimizer.
%     - Pass on fitting length to another function to validate fit for other muscle length.
%     - Plot results of both muscles using new fitting length.
%     - Nested minimization function to use either Simplx, Pattern search, or Pareto search.
%     - Nested function that calculates SSE and RMSE of new fitting value.
%     - Nested function that calls muscle class and calculates torque.

clear; clc; close all

load KneeFlxPin_10mm_48cm.mat phiD Location Name CrossPoint Dia T rest tendon kmax fitting
load Plot_KneeFlxPin_10mm_48cm.mat Angle Torque pres 
pres = mean(pres);                  %Make pressure a scalar value
y = Torque';                        %Make it just the data

x0 = 0.0351;                        %Initial fitting length guess
%SIMPLX search for SSE, x = 0.0347; for RMSE x = 0.033 
%pattern search for SSE, x = 0.03613; for RMSE x = 0.0358.

[point,value] = minimizeFlx10mm_nest(x0,y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, pres);

%% Validate with other muscle length in this configuration
P = value;
par = point;
point = mean(point);

fprintf('fitting length is %5d with RMSE=%5d, adj. R^2=%5d, and max residual=%5d\n',point,value)

[RMSE, fvu, maxResid, A2, T2, v2, r2] = validateFlxPin10mm(point);
Rsq = 1 - fvu;
fprintf('Validation returns RMSE=%5d, R^2=%5d, and max residual=%5d\n',RMSE,Rsq,maxResid)



%% Plot the results
%Matlab hex color values:
c1 = '#FFD700'; %gold
c2 = '#FFB14E'; %orange
c3 = '#FA8775'; %light orange
c4 = '#EA5F94'; %pink
c5 = '#CD34B5'; %magenta
c6 = '#9D02D7'; %magenta 2
c7 = '#0000FF'; %indigo
c8 = '#000000'; %black
sz = 60;        %size of data points
c = {c1; c2; c3; c4; c5; c6; c7; c8};

%Call class for optimized muscle, use torque
Bifemsh_Pam_new = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, point, pres);
Theoretical = Bifemsh_Pam_new.Torque(:,3);

%Plot Torque values against theoretical for optimized muscle
figure
scM = scatter(Angle,Torque,sz,'filled','MarkerFaceColor',c7,'DisplayName','Measured Torque');
hold on
Disp = sprintf('Theoretical Torque');
PL = plot(phiD, Theoretical,'Color',c5,'Linewidth',2,'DisplayName',Disp);
title(sprintf('Isometric Torque vs Knee Angle, 10mm Flexor, %4.3f m long',rest),'FontSize',12,'FontWeight','Bold')
xlabel('Knee angle, \circ','FontSize',12,'Interpreter','tex')
ylabel('Torque, N \cdot m','FontSize',12,'Interpreter','tex')
ax2 = gca;
ax2.FontSize = 12;
ax2.FontWeight = 'bold';
ax2.FontName = 'Arial';
ax2.YAxis.LineWidth = 2; ax2.YAxis.FontSize = 10;
ax2.XAxis.LineWidth = 2; ax2.XAxis.FontSize = 10;
lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'southwest';
hold off

%Plot validation
figure
scM2 = scatter(A2,T2,sz,'filled','MarkerFaceColor',c7,'DisplayName','Measured Torque');
hold on
Disp2 = sprintf('Theoretical Torque');
PL2 = plot(phiD, v2,'Color',c5,'Linewidth',2,'DisplayName',Disp2);
title(sprintf('Isometric Torque vs Knee Angle, 10mm Flexor, %4.3f m long',r2),'FontSize',12,'FontWeight','Bold')
xlabel('Knee angle, \circ','FontSize',12,'Interpreter','tex')
ylabel('Torque, N \cdot m','FontSize',12,'Interpreter','tex')
ax2 = gca;
ax2.FontSize = 12;
ax2.FontWeight = 'bold';
ax2.FontName = 'Arial';
ax2.YAxis.LineWidth = 2; ax2.YAxis.FontSize = 10;
ax2.XAxis.LineWidth = 2; ax2.XAxis.FontSize = 10;
lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'southwest';
hold off


%% Nested function for optimization
   function [X, FVAL] = minimizeFlx10mm_nest(x0,y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, pres)

    fun = @(x)minimize(y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres);
    
%     options = optimset('Display','iter','PlotFcns',@optimplotfval);     %fminsearch options
%     [X, FVAL, flag, put] = fminsearch(fun,x0,options);

     [X, FVAL] = bnd_con_pen_nelder_mead( fun, x0, [0.03 0.04], [], 100, [], 1, 1);
    
%     options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestx);        %pattern search options
%     [X, FVAL, flag, put] = patternsearch(fun,x0,[],[],[],[],0.02,.04,[],options);
    
%     options = optimoptions('paretosearch','PlotFcn',{@psplotbestf, @psplotfuncount, @psplotparetof},'InitialPoints',x0);
%     [X, FVAL, flag, put] = paretosearch(fun,1,[],[],[],[],0.02,0.04,[],options);
    
%     hybridopts = optimoptions('patternsearch','Display','iter');
%     options = optimoptions('gamultiobj','Display','iter','PlotFcn',{@gaplotpareto, @gaplotdistance, @gaplotscores},'HybridFcn',{@fgoalattain,hybridopts});
%     [X, FVAL, flag, put] = gamultiobj(fun,1,[],[],[],[],0.02,0.04,[],options);
    
    disp(X)
    disp(FVAL)

%% Nested function that calculates RMSE, FVU, and max error
        function f = minimize(y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres)

        ynew = nestedfun1(phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres);
          
        yresid = y-ynew;                     %residual error
        maxResid = max(yresid);             %largest residual
        SSresid = sum(yresid.^2,'omitnan'); %Sum of squares of the residual
        SStot = sum((y-mean(y)).^2);        %total sum of squares
        n = (sum(~isnan(ynew)));             %number of data points
        RMSE = sqrt(SSresid/n);             % RMSE for function 1
        fvu = SSresid*(n-1)/(SStot*n);      %Fraction of Variance Unexplained (FVU),adjusted, Rsq = 1 - FVU  

        f1 = RMSE;
        f2 = fvu;                       
        f3 = maxResid;
       
%         f1 = [f1, f1_other];
%         f = f1_other(4);
          f = [f1, f2, f3];
%         f = mean(f1,'omitnan');     %Combine the error from the different equation fits if not pareto search
  
            function val = nestedfun1(phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres)
                                
           % Find the difference between experimental and calculated results
            fitt = x;
            Bifemsh_Pam_adj = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitt, pres); %Call the class to do the Torque calculation
            Tqz = Bifemsh_Pam_adj.Torque(:,3);    %Calculated Adjusted Torque

            %fit options
            mod_Pam = fittype('cubicinterp');
            Options = fitoptions(mod_Pam);
            Options.Normal = 'on';
            
            %Get model torque values at each angle there is measurement data for
             Options.Exclude = isnan(Tqz);
             mdl_Pam = fit(phiD',Tqz,mod_Pam,Options);
             val = feval(mdl_Pam,Angle');
            end

    end
end