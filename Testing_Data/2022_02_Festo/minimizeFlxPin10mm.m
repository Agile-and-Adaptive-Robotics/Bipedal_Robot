
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
[SSE, RMSE, A2, T2, v2, r2] = validateFlxPin10mm(point);

fprintf('fitting length is %5d with a SSE of %5d\n',point,value)

fprintf('Validation returns, SSE of %5d with an RMSE of %5d\n',SSE(4),RMSE(4))



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
Bifemsh_Pam_new = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, point, pres);
TorqueR = Bifemsh_Pam_new.Torque(:,3,:);
Theoretical = cell(2, size(TorqueR,3));
Theoretical{1,1} = 'Hunt Eq.';
Theoretical{1,2} = 'Exponential Eq.';
Theoretical{1,3} = 'Polynomial Eq.';
Theoretical{1,4} = 'Exponential Eq., Simplified';
Theoretical{1,5} = 'Polynomial Eq., Simplified';
for i = 1:length(Theoretical)
    Theoretical{2,i} = TorqueR(:,:,i)';
end

PL = cell(1, size(Theoretical,2));

%Plot Torque values against theoretical for optimized muscle
figure
scM = scatter(Angle,Torque,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Measured Torque');
hold on
    for i = 4 %1:size(Theoretical,2)
        txt = Theoretical{1,i};
        T1 = 2*i-1;
        Disp{T1} = sprintf('Theoretical Torque');
        PL{i} = plot(phiD, Theoretical{2,i},'Color',c{i},'Linewidth',2,'DisplayName',Disp{T1});
    end
title(sprintf('Isometric Torque vs Knee Angle, 10mm Flexor, %4.3f m long',rest),'FontSize',12,'FontWeight','Bold')
xlabel('Knee angle, \circ','FontSize',12)
ylabel('Torque, N \bullet m','FontSize',12)
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

PL2 = cell(1, size(Theoretical,2));

%Plot validation
figure
scM2 = scatter(A2,T2,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Measured Torque');
hold on
    for i = 4   %1:size(Theoretical,2)
        txt = Theoretical{1,i};
        T1 = 2*i-1;
        Disp2{T1} = sprintf('Theoretical Torque');
        PL2{i} = plot(phiD, v2{i},'Color',c{i},'Linewidth',2,'DisplayName',Disp2{T1});
    end
title(sprintf('Isometric Torque vs Knee Angle, 10mm Flexor, %4.3f m long',r2),'FontSize',12,'FontWeight','Bold')
xlabel('Knee angle, \circ','FontSize',12)
ylabel('Torque, N \bullet m','FontSize',12)
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

%% Nested function that calculates SSE and RMSE
        function f = minimize(y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres)

        ynew = nestedfun1(phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres);
        
        yresid = cell(length(ynew),1);
        SSresid = cell(length(ynew),1);
        fu = cell(length(ynew),1);
        
        for i = 1:length(ynew)
            yresid{i} = y-ynew{i};              %residual error
            SSresid{i} = sum(yresid{i}.^2,'omitnan'); %Sum of squares of the residual
            fu{i} = sqrt(SSresid{i}/length(yresid{i}));        % RMSE for function 1
        end
        
        f1 = zeros(length(SSresid),1);
        f2 = zeros(length(fu),1);
        for i = 1:length(fu)
            f1(i) = SSresid{i};
            f2(i) = fu{i};                       %This works as f for pareto search
        end
       
%        [f1_other, f2_other] =  validateFlxPin10mm(x);       %Combine both muscle lengths
       
%          f1 = [f1, f1_other];
%         f = f1_other(4);
      f = f1(4);
%       f = mean(f1,'omitnan');     %Combine the error from the different equation fits if not pareto search

  
  
        
            function val = nestedfun1(phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres)
                                
            % Find the difference between experimental and calculated results
            % Inputs:
            % fit = x/100; Make x larger than fitting size by 100 so that
            % search functions can work better
            
            fitt = x;
            Bifemsh_Pam_adj = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitt, pres); %Call the class to do the Torque calculation
            
            Tqz = cell(size(Bifemsh_Pam_adj.Torque,3),1);
            for j = 1:size(Bifemsh_Pam_adj.Torque,3)
                Tqz{j} = Bifemsh_Pam_adj.Torque(:,3,j);    %Calculated Adjusted Torque
            end
            

            %fit options
            mod_Pam = fittype('cubicinterp');
            Options = fitoptions(mod_Pam);
            Options.Normal = 'on';
            
            %prepare cells
            mdl_Pam = cell(size(Bifemsh_Pam_adj.Torque,3),1);
            val = cell(length(mdl_Pam),1);
            
            %Get values at each angle there is measurement data for
                for j = 1:length(Tqz)
                 Options.Exclude = isnan(Tqz{j});
                 mdl_Pam{j} = fit(phiD',Tqz{j},mod_Pam,Options);
                 val{j} = feval(mdl_Pam{j},Angle');
                end
            end

    end
end