
function [SSE, RMSE] = validateFlxPin10mm(fitting)
%Validation of minimizeFlxPin10mm
%   Validate fitting length obtained by the other function.
load KneeFlxPin_10mm_48cm.mat phiD Location Name CrossPoint Dia T rest tendon kmax
load Plot_KneeFlxPin_10mm_48cm.mat Angle Torque pres 
pres = mean(pres);                  %Make pressure a scalar value
y = Torque';                        %Make it just the data

Bifemsh_Pam_adj = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres); %Call the class to do the Torque calculation
  
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

        yresid = cell(length(val),1);
        SSresid = cell(length(val),1);
        fu = cell(length(val),1);
        
        for i = 1:length(val)
            yresid{i} = y-val{i};              %residual error
            SSresid{i} = sum(yresid{i}.^2,'omitnan'); %Sum of squares of the residual
            fu{i} = sqrt(SSresid{i}/length(yresid{i}));        % RMSE for function 1
        end

        SSE = zeros(length(SSresid),1);
        RMSE = zeros(length(fu),1);
        for i = 1:length(fu)
            SSE(i) = SSresid{i};                   %SSE values
            RMSE(i) = fu{i};                       %RMSE for different equations
        end  
    
end

