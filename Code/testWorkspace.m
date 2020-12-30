% clc                 %Command Line Clear
% clear               %Clear the workspace of stored variables
% close all           %Close all open figures
% 

sig = 15000/(1.75*0.25);

Kf = 1 + 0.95*(2.6 - 1);

sigMax = Kf * sig;
sigMax = sigMax*10^(-3);

a = (0.76 * 260)^2/25;

b = -1/3*log10(0.76*260/25);

N = (sigMax/a)^(1/b)