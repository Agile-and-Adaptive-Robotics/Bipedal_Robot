

rest = 0.415;
tendon = 0.012;
kmax = 0.350;
x0 = [rest, tendon, kmax];
[point,value] = minimizeFlx10mm_nest(x0);


   function [x, fval] = minimizeFlx10mm_nest(x0)

    options = optimset('Display','iter','PlotFcns',@optimplotfval);
%     fun = @(x)minimize;
    [x, fval] = fminsearch(@minimize,x0,options);

        function f = minimize(x)
        load KneeFlx_10mm_42cm.mat phiD Location Name CrossPoint Dia T_Pam fitting
        load Plot_KneeFlx_10mm_42cm.mat X1 Angle Torque mdl1 TorqueMean pres 
        pres = mean(pres);                  %Make pressure a scalar value
        y1 = Torque;                        %Make it just the data
        y = feval(mdl1, X1);               %Make y1 the curve fit
        y3 = diff(y);
        c = [1 1];                          %Weights c1 & c2 for functions 1 & 2, respectively
        
        [y2, y4] = nestedfun1(x);
        yresid1 = y1-y2;                     %residual error
        SSresid1 = sum(yresid1.^2);          %Sum of squares of the residual
        %SStotal = (length(y1)-1)*var(y1);   %total sum of squares

        yresid2 = y3-y4;                     %residual error from derivatives
        SSresid2 = sum(yresid2.^2);          %Sum of squares of the residual from the derivatives

        f = c(1)*SSresid1+c(2)*SSresid2;     %Combine the error from the data and the error of the slopes


            function [val, der] = nestedfun1(x)
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