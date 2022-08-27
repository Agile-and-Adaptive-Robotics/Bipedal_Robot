%% Function f10 solves Dr. Hunt's equation/data for force
% This equation and its parameters was found using the curve fitting app
% and matching the data in the lookup table. Only the first zero value in
% each column was saved.

%Input:
%  x = relative strain
%  y = Pressure
%Output:
%  f10 = Force

function F = f10(x,y)
a0 = -133;
a1 = 7.981;
a2 = 5.139;
a3 = 1.135;
a4 = -0.4361;
a5 = 0.04946;

F = a0 + exp(-a1*x+a2) + y*exp(-(a3*x).^2+a4)+a5*y;