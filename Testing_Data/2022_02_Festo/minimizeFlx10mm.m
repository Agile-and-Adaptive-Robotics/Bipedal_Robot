
%Load 

clear; clc; close all

load KneeFlx_10mm_42cm.mat phiD Location Name CrossPoint Dia T_Pam fitting
load Plot_KneeFlx_10mm_42cm.mat X1 Angle Torque mdl1 TorqueMean pres 
pres = mean(pres);                  %Make pressure a scalar value
y1 = Torque';                        %Make it just the data
y = feval(mdl1, X1);               %Make y1 the curve fit
y3 = diff(y);
weight = [100 0];               %Assign relative weight to  functions 1 and 2
weight_norm = norm(weight);         %Find magnitude of vector "weight"
c = weight/weight_norm;             %Weights c1 & c2 for functions 1 & 2, respectively, in unit vector form

rest = 0.420;
tendon = 0.0345;
kmax = 0.350;
x0 = [rest, tendon, kmax];
[point,value] = minimizeFlx10mm_nest(x0,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);


   function [X, FVAL] = minimizeFlx10mm_nest(x0,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c)

    fun = @(x)minimize(x,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);
    
    options = optimset('Display','iter','PlotFcns',@optimplotfval);     %fminsearch options
    [X, FVAL] = fminsearch(fun,x0,options);

%     options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);        %pattern search options
%     [X, FVAL] = patternsearch(fun,x0,[],[],[],[],[0.4 0.01 0.345],[.44 .06 .355],[],options);
    
%     options = optimoptions('paretosearch','PlotFcn','psplotparetof','InitialPoints',x0);
%     [X, FVAL] = paretosearch(fun,3,[1 1 0],0.6,[],[],[0.4 0.001 0.3],[0.55 0.1 0.4],[],options);
    
    disp(X)
    disp(FVAL)

        function f = minimize(x,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c)

%        call2 = {X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c};
        [y2, y4] = nestedfun1(x,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD);
        yresid1 = y1-y2;                     %residual error
        SSresid1 = sum(yresid1.^2);          %Sum of squares of the residual
%         SStotal1 = sum((y1-mean(y1)).^2);    %total sum of squares
%         f1 = SSresid1/SStotal1;              %SSE/SST
        f1 = SSresidl/length(SSresid1);      % RMSE for function 1

        yresid2 = y3-y4;                     %residual error from derivatives
        SSresid2 = sum(yresid2.^2);          %Sum of squares of the residual from the derivatives
%         SStotal2 = sum((y3-mean(y3)).^2);    %total sum of squares
%         f2 = SSresid2/SStotal2;              %SSE/SST
        f2 = SSresid2/length(SSresid2);      % RMSE for function 2

        f = c(1)*f1+c(2)*f2;     %Combine the error from the data and the error of the slopes
%         f = [c(1)*f1; c(2)*f2];  %For pareto search. f1 is data match, f2 is slope match
        
        % Add term Cost = (1/(max_length-length)^2 + 1/(length)   
        
            function [val, der] = nestedfun1(x,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD)
            % Find the difference between experimental and calculated results
            % Inputs:
            %Rest = resting length, meters
            %tendon = tendon length, meters
            %kmax = maximum contracted length, meters (will be converted to percent later)
            Rest = x(1);
            Teender = x(2);
            kmax1 = x(3);
            Bifemsh_Pam_adj = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T_Pam, Rest, kmax1, Teender, fitting, pres); %Call the class to do the Torque calculation
            T2 = Bifemsh_Pam_adj.Torque(:,3);             %Calculated Adjusted Torque
            mod_Pam = fittype('smoothingspline');
            Options = fitoptions(mod_Pam);
            Options.Normal = 'on';
            Options.Exclude = isnan(T2);
            mdl_Pam = fit(phiD',T2,mod_Pam,Options);
            val = feval(mdl_Pam,X1);
            der = diff(val);
                       
            end

    end
end