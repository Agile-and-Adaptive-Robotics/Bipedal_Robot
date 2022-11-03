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
a0 = -0.6464;
a1 = 2.929;
a2 = -0.4217;
a3 = 0.362;
a4 = 0.2403;
a5 = -0.2731;
Fn1 = a0+exp(-a1.*x+a2)+y.*exp(-a3*x.^2+a4.*y)+a5.*y;

%Polynomial function
p01 = 0.7272;
p02 = 0.2743;
p10 = -0.9639;
p11 = -1.633;
p12 = 0.7778;
p20 = 1.559;
p30 = -0.7412;
Fn2 = p10*x + p01*y + p20*x.^2 + p11*x.*y + p02*y.^2 + p30*x.^3 + p12*x.*y.^2;

%Exponential function, simplified
c1 = 1.7;
c2 = 0.2;
Fn3 = -1 + exp(-c1.*x) + y.*exp(-c2*x.^2);

%"Polynomial function, simplified"
d1 = 0.03243;
Fn4 = -(x-y).*(x.^2-x+1)+ d1.*(y.^3-1);

% fprintf('Exponential equation\n Fn1 = %4.2f \n \n', Fn1);
% fprintf('Polynomial function\n Fn2 = %4.2f \n \n',Fn2);
% fprintf('Exponential function, simplified \n Fn3 = %4.2f \n \n',Fn3);
% fprintf('Polynomial function, simplified \n Fn4 = %4.2f \n \n',Fn4);

end
