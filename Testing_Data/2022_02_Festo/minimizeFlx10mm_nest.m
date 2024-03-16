   function f = minimizeFlx10mm_nest(x,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c)

%     options = optimset('Display','iter','PlotFcns',@optimplotfval);
%     [x, fval] = fminsearch(@minimize,x0,options);
% 
%         function f = minimize(x)


        [y2, y4] = nestedfun1(x);
        yresid1 = y1-y2;                     %residual error
        SSresid1 = sum(yresid1.^2);          %Sum of squares of the residual
        SStotal1 = sum((y1-mean(y1)).^2);    %total sum of squares
        f1 = SSresid1/SStotal1;              %SSE/SST, fraction of variance unexplained

        yresid2 = y3-y4;                     %residual error from derivatives
        SSresid2 = sum(yresid2.^2);          %Sum of squares of the residual from the derivatives
        SStotal2 = sum((y3-mean(y3)).^2);    %total sum of squares
        f2 = SSresid2/SStotal2;              %SSE/SST, fraction of variance unexplained


        f = c(1)*f1 + c(2)*f2;     %Combine the error from the data and the error of the slopes
 

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

%     end
end