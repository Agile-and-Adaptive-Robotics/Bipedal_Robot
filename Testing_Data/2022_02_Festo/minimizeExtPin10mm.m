%% Optimize predicted torque for extensors.
function x = minimizeExtPin10mm

load KneeFlxPin_10mm_46cm.mat Bifemsh_Pam phiD
theta = acos(Bifemsh_Pam.UnitDirection);
F = Bifemsh_Pam.Force;
mA = Bifemsh_Pam.MomentArm;

load Plot_KneeFlxPin_10mm_46cm.mat Angle Torque ICRtoMuscle TorqueHand
A = sortrows([Angle, Torque],1);
Ang_meas = A(:,1);       %Experimental angles
T_meas = A(:,2);        %Measured Torque
mA_hyb = ICRtoMuscle;
T_hyb = TorqueHand;     %Torque calculated with the hybrid method
Tnew = zeros(length(T_meas),3);

%% Nonlinear least squares fit
x0 = 1;         %Initial guess
fun = @sse;
% options = optimoptions('fmincon','FunValCheck','off');
x = bnd_con_pen_nelder_mead( fun, x0, [0 pi], [], [], 100, 1, 1);
Tval = nested2;

%% Plot
figure
hold on
scatter(Ang_meas,T_meas,[],'filled','DisplayName','Experiment')
scatter(Ang_meas,Tval,[],'filled','DisplayName','New predict')
plot(phiD',Bifemsh_Pam.Torque(:,3),'DisplayName','Predict original')
hold off
legend

%% Nested function
    function t = sse(k)
     Fnew = F-k*F.*Bifemsh_Pam.UnitDirection(:,1);   
     Tpredict1 = cross(mA,Fnew);
     Tpredict2 = griddedInterpolant(phiD',Tpredict1(:,3));
     Tnew = Tpredict2(Ang_meas);
     [RMSE, t, maxResid] = Go_OfF(Torque,Tnew);
%      t = sum((Torque-Tnew).^2,'omitnan');
    end

    function [Tv, gof] = nested2
        Tv = Tnew;
        [RMSE, t, maxResid] = Go_OfF(Torque,Tnew);
        gof = [RMSE, t, maxResid];
        
    end




end


    