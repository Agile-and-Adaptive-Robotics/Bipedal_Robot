
load KneeFlx_10mm_42cm.mat TorqueR phiD Location Name kmax CrossPoint Dia T_Pam fitting
Theoretical = TorqueR(:,3)';   %Load theoretical torque
load Plot_KneeFlx_10mm_42cm.mat X1 Angle Torque modp mdl1 gofp2 TorqueMean pres 
pres = mean(pres);                  %Make pressure a scalar value
y1 = Torque;                        %Make y1 the data points
y = feval(mdl1, X1);                %Make y the curve fit
y3 = diff(y);                       %y3 is the derivative of the curve fit

rest = 0.415;                       %x(1), measured
tendon = 0.012;                     %x(2), guess
kmax = 0.350;                       %x(3), measured
x0 = [rest, tendon, kmax];
c = [1 10];                         %c1 is value of importance for matching data points, c2 is value of importance for matching slope of fitted lines
[X, FVAL] = minimizeFlx10mm_nest(x0,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c)