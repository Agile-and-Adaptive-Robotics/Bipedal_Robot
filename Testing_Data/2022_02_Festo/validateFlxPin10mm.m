
function [RMSE, fvu, maxResid, Angle, Torque, Tqz, rest] = validateFlxPin10mm(fitting)
%Validation of minimizeFlxPin10mm
%   Validate fitting length obtained by the other function.
load KneeFlxPin_10mm_46cm.mat phiD Location Name CrossPoint Dia T rest tendon kmax
load Plot_KneeFlxPin_10mm_46cm.mat Angle Torque pres 
% pres = mean(pres);                  %Make pressure a scalar value
y = Torque';                        %Make it just the data

Bifemsh_Pam_adj = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres); %Call the class to do the Torque calculation
  
Tqz = Bifemsh_Pam_adj.Torque(:,3);
%fit options
mod_Pam = fittype('cubicinterp');
Options = fitoptions(mod_Pam);
Options.Normal = 'on';
            
%Get values at each angle there is measurement data for
 Options.Exclude = isnan(Tqz);
 mdl_Pam = fit(phiD',Tqz,mod_Pam,Options);
 val = feval(mdl_Pam,Angle');
        
yresid = y-val;                     %residual error
maxResid = max(yresid);             %largest residual
SSresid = sum(yresid.^2,'omitnan'); %Sum of squares of the residual
SStot = sum((y-mean(y)).^2);        %total sum of squares
n = (sum(~isnan(val)));             %number of data points
RMSE = sqrt(SSresid/n);             % RMSE for function 1
fvu = SSresid*(n-1)/(SStot*n);      %Fraction of Variance Unexplained (FVU),adjusted, Rsq = 1 - FVU       
   
end

