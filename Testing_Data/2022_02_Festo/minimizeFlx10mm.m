function [Afinal, M, I] = minimizeFlx10mm
load KneeFlx_10mm_42cm.mat TorqueR phiD Location Name kmax CrossPoint Dia T_Pam fitting
Theoretical = TorqueR(:,3)';   %Load theoretical torque
load Plot_KneeFlx_10mm_42cm.mat X1 Angle Torque modp mdl1 gofp2 TorqueMean pres 
pres = mean(pres);                  %Make pressure a scalar value
%y1 = feval(mdl1, X1);               %Make y1 the curve fit
y1 = Torque;                        %Make it just the data
rest = 0.415;
tendon = 0.012;
kmax = 0.350;
%KMAX = kmax;

deltaTendon = (-5:30)/1000;
deltaRest = (-10:20)/1000;
deltaK = (-10:10)/1000;
A = NaN(length(deltaTendon),length(deltaRest),length(deltaK));
%A = NaN(length(deltaTendon),length(deltaRest));

  for t = 1:length(deltaTendon)         %vary the tendon length
    for u = 1:length(deltaRest)         %vary the resting length
        for v = 1:length(deltaK)       %vary the kmax length
            Tendon(t) = tendon+deltaTendon(t);
            restingL(u) = rest+deltaRest(u);
            KMAX(v) = kmax+deltaK(v);
            y2 = nestedfun1(restingL(u),Tendon(t),KMAX(v));
            yresid = y1-y2;                     %residual error
            SSresid = sum(yresid.^2);            %Sum of squares of the residual
            %SStotal = (length(y1)-1)*var(y1);   %total sum of squares
            A(t,u,v) = SSresid;
        end
    end
  end

Afinal = A
[M,I] = min(A,[],'all','linear');         %M finds the minimum value of A over all dimensions. I is the index location.
M
I
[Irow, Icolumn, I3] = ind2sub([length(deltaTendon),length(deltaRest), length(deltaK)],I)

    function val = nestedfun1(restingL, Tendon, KMAX)
        % Find the difference between experimental and calculated results
        % Inputs:
        %Rest = resting length, meters
        %tendon = tendon length, meters
        %kmax = maximum contracted length, meters (will be converted to percent later)
        Rest = restingL;
        kmax1 = KMAX;
        Teender = Tendon; 
%         pres = mean(pres);         %average pressure, make sure 
        Bifemsh_Pam_adj = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T_Pam, Rest, kmax1, Teender, fitting, pres); %Call the class to do the Torque calculation
        T2 = Bifemsh_Pam_adj.Torque(:,3);             %Calculated Adjusted Torque
        mod_Pam = fittype('smoothingspline');
        Options = fitoptions(mod_Pam);
        Options.Normal = 'on';
        Options.Exclude = isnan(T2);
        mdl_Pam = fit(phiD',T2,mod_Pam,Options);
        val = feval(mdl_Pam,X1);
    end

end