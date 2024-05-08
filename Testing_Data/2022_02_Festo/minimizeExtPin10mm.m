%% Optimize predicted torque for extensors.
function f = minimizeExtPin10mm(a)

%% 46cm length
%kf = knee flexor, kf(1) = pinned joint, kf(2) = biomimetic;
%ke = knee extensor, same as above
%kf.L := lengths = [42 46 48] cm
%example: kf(1).L(2).Mz z-axis torque for pinned knee, flexor, 46cm length

load KneeFlxPin_10mm_46cm.mat Bifemsh_Pam phiD
Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
load Plot_KneeFlxPin_10mm_46cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
kf(1).L(2) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
                  'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
                  'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
                  'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
                  'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, ...
                  'Mpre',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
                  'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5));
clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A

%% 48 cm length
load KneeFlxPin_10mm_48cm.mat Bifemsh_Pam phiD
Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
load Plot_KneeFlxPin_10mm_48cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
kf(1).L(3) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
                  'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
                  'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
                  'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
                  'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, ...
                  'Mpre',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
                  'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5));
clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A

%% What to optimize based on selection
Aexp = kf(1).L(a).Mexp;      %Experimental angles
Mexp = kf(1).L(a).Mexp;        %Measured Torque
mA_h = kf(1).L(a).mA_h;        %Measured moment arm
M_h = kf(1).L(a).M_h;          %Torque calculated with the hybrid method
M_p = zeros(length(kf(1).L(a).Mexp)); %Optimized torque prediction

%% asdf


%% Nonlinear least squares fit
x0 = 1;         %Initial guess
fun = @SSE;
% options = optimoptions('fmincon','FunValCheck','off');
x = bnd_con_pen_nelder_mead( fun, x0, [0 pi], [], [], 100, 1, 1);
Tval = nested2;

%% Plot Optimized fit
figure
hold on
scatter(kf(1).L(a).Aexp,kf(1).L(a).Mexp,[],'filled','DisplayName','Experiment')
scatter(kf(1).L(a).Aexp,Tval,[],'filled','DisplayName','New predict')
plot(kf(1).L(a).Ak,kf(1).L(a).Mpre,'DisplayName','Predict original')
hold off
legend

%% Plot validation
figure
hold on
scatter(Aexp,Mexp,[],'filled','DisplayName','Experiment')
scatter(Aexp,Tval,[],'filled','DisplayName','New predict')
plot(,Mpre,'DisplayName','Predict original')
hold off
legend
%% Subfunctions
function t = SSE(k)
     Fnew = F-k*F.*Bifemsh_Pam.UnitDirection(:,1);   
     Tpredict1 = cross(mA,Fnew);
     Tpredict2 = griddedInterpolant(phiD',Tpredict1(:,3));
     M_p = Tpredict2(Ang_meas);
     [RMSE, fvu, maxResid] = Go_OfF(Torque,M_p);
     t = [RMSE, fvu, maxResid];
end

function [Tv, gof] = nested2
    Tv = M_p;
    [RMSE, t, maxResid] = Go_OfF(Torque,M_p);
    gof = [RMSE, t, maxResid];
        
end




end


    