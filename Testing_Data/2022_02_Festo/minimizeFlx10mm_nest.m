   function f = minimizeFlx10mm_nest(Name, Location, CrossPoint, Dia, T_Pam, Rest, kmax1, Teender, x, pres, phiD, y1, c)

        y2 = nestedfun1(x);
        yresid1 = y1-y2;                     %residual error
        maxError = max(yresid1);
        SSresid1 = sum(yresid1.^2);          %Sum of squares of the residual
        SStotal1 = sum((y1-mean(y1)).^2);    %total sum of squares
        n = (sum(~isnan(val)));             %number of data points
        RMSE = sqrt(SSresid/n);             % RMSE for function 1
        fvu = SSresid1/SStotal1;              %SSE/SST, fraction of variance unexplained
        
        f1 = RMSE;
        f2 = fvu;
        f3 = maxError;
        f = c(1)*f1 + c(2)*f2 + c(3)*f3;     %Combine the error from the data and the weights
 

            function val = nestedfun1(x)
            % Find the difference between experimental and calculated results
            % Inputs:
            %Rest = resting length, meters
            %tendon = tendon length, meters
            %kmax = maximum contracted length, meters (will be converted to percent later)
            fitting = x;
            Bifemsh_Pam_adj = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T_Pam, Rest, kmax1, Teender, fitting, pres); %Call the class to do the Torque calculation
            T2 = Bifemsh_Pam_adj.Torque(:,3);             %Calculated Adjusted Torque
            F = griddedInterpolant(phiD',T2);
            val = F(Angle);
            end

end