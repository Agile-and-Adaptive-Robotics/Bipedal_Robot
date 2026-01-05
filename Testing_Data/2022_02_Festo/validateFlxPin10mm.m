
function [RMSE, fvu, maxResid, Angle, Torque, Tnew, Torig, rest] = validateFlxPin10mm(fitting)
%Validation of minimizeFlxPin10mm
%   Validate fitting length obtained by the other function.
%   Goodness of fit measures
%   Angle and torque from experiments
%   Calculated torque with new fitting
%   resting length of muscle
load KneeFlxPin_10mm_46cm.mat phiD Location Name CrossPoint Dia T rest tendon kmax TorqueR
load Plot_KneeFlxPin_10mm_46cm.mat Angle Torque pres 
Torig = TorqueR(:,3);
pres = mean(pres);                  %Make pressure a scalar value
y = Torque;                        %Make it just the data

Bifemsh_Pam_adj = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres); %Call the class to do the Torque calculation
Tnew = Bifemsh_Pam_adj.Torque(:,3);
            
%Get expected values at each angle there is experimental data for
F = griddedInterpolant(phiD',Tnew);
val = F(Angle);
        
yresid = y-val;                     %residual error
maxResid = max(abs(yresid));             %largest residual
SSresid = sum(yresid.^2,'omitnan'); %Sum of squares of the residual
SStot = sum((y-mean(y)).^2);        %total sum of squares
n = (sum(~isnan(val)));             %number of data points
RMSE = sqrt(SSresid/n);             % RMSE for function 1
fvu = SSresid*(n-1)/(SStot*n);      %Fraction of Variance Unexplained (FVU),adjusted, Rsq = 1 - FVU       
   
end

