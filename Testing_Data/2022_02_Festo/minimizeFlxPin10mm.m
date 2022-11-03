
%% Optimization for fitting length

clear; clc; close all

load KneeFlxPin_10mm_48cm.mat phiD Location Name CrossPoint Dia T rest tendon kmax fitting
load Plot_KneeFlxPin_10mm_48cm.mat Angle Torque pres 
pres = mean(pres);                  %Make pressure a scalar value
y = Torque';                        %Make it just the data

x0 = 0.035;
[point,value, exit, out] = minimizeFlx10mm_nest(x0,y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, pres);

%% Validate with other muscle length in this configuration
% P = value;
% par1 = mean(point);
% par2 = mode(point);
[SSE, RMSE] = validateFlxPin10mm(point);

fprintf('fitting length is %5d with a SSE of %5d\n',point,value)

fprintf('Validation returns, SSE of %5d with an RMSE of %5d\n',SSE,RMSE)



%% Plot the results


%% Nested function for optimization
   function [X, FVAL, flag, put] = minimizeFlx10mm_nest(x0,y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, pres)

    fun = @(x)minimize(y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres);
    
%     options = optimset('Display','iter','PlotFcns',@optimplotfval);     %fminsearch options
%     [X, FVAL, flag, put] = fminsearch(fun,x0,options);

    options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);        %pattern search options
    [X, FVAL, flag, put] = patternsearch(fun,x0,[],[],[],[],0.01,.05,[],options);
    
%     options = optimoptions('paretosearch','PlotFcn','psplotparetof','InitialPoints',x0);
%     [X, FVAL, flag, put] = paretosearch(fun,1,[],[],[],[],0.01,0.05,[],options);
    
    disp(X)
    disp(FVAL)

%% Nested function that calculates SSE and RMSE
        function f = minimize(y, phiD, Angle, Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, x, pres)

%        call2 = {X1,Name, Location, CrossPoint, Dia, T, fitting, pres, phiD, y1, y3, c};
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
            f1(i) = SSresid{1};
            f2(i) = fu{i};                       %This works as f for pareto search
        end
        
%        f = [f1; f2];
       f = min(f1);     %Combine the error from the different equation fits if not pareto search

  
  
        
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