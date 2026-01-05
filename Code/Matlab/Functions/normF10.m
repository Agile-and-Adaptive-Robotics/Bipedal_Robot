%% A function that gives normalized force for a 10mm Festo BPA, given pressure and relative strain
% Inputs:
%   RelStrain = Relative Strain. 
%   Pressure = Pressure (kPa).
% Outputs:
%   - Fn1 = Force, unitless. Exponential function fit.
%   - Fn2 = Force, unitless. Polynomial function fit.
%   - Fn3 = Force, unitless. Simplified exponential function fit.
%   - Fn4 = Force, unitless. Simplified polynomial function fit.

function [Fn1, Fn2, Fn3, Fn4] = normF10(RelStrain, Pressure)
minArgs = 1;
maxArgs =4;
nargoutchk(minArgs,maxArgs)
% disp("You requested " + nargout + " outputs.")

x = RelStrain;
y = Pressure/620;

%Exponential equation
a0 =     -0.6291;
a1 =       3.281;
a2 =      -0.459;
a3 =      0.3834;
a4 =      0.2034;
a5 =     -0.2279;
Fn1 = a0+exp(-a1.*x+a2)+y.*exp(-a3*x.^2+a4.*y)+a5.*y;

%Polynomial function
p01 = 0.7477;
p02 = 0.2525;
p10 = -1.02;
p11 = -1.354;
p12 = 0.5467;
p20 = 1.472;
p30 = -0.645;
Fn2 = p10*x + p01*y + p20*x.^2 + p11*x.*y + p02*y.^2 + p30*x.^3 + p12*x.*y.^2;

%Exponential function, finalized
b0 = -0.6144;
b1 = 3.559;
b2 = 0.6111;
b3 = 0.4963;
Fn3 = b0+b2*exp(-b1.*x)+y.*exp(-b3*(x.^2));

%"Polynomial function, simplified"
d1 = 0.06965;
Fn4 = -(x-y).*(x.^2-x+1)+ d1.*(y.^3-1);

% fprintf('Exponential equation\n Fn1 = %4.2f \n \n', Fn1);
% fprintf('Polynomial function\n Fn2 = %4.2f \n \n',Fn2);
% fprintf('Exponential function, finalized \n Fn3 = %4.2f \n \n',Fn3);
% fprintf('Polynomial function, simplified \n Fn4 = %4.2f \n \n',Fn4);

end
