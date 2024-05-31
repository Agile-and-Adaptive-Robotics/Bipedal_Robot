%% Function f10 solves Dr. Hunt's equation/data for force
% This equation and its parameters was found using the curve fitting app
% and matching the data in the lookup table. Only the first zero value in
% each column was saved.

%Input:
%  x = relative strain
%  y = Pressure, kPa
%Output:
%  f10 = Force, N

function F = f10(x,y)
% coefficients w/ 95% confidence interval
a0 = -139.9;        % [-141.9, 137]
a1 = 6.313;         % [6.155, 6.472]
a2 = 5.166;         % [5.154, 5.177]
a3 = 0.9434;        % [0.9111, 0.9757]
a4 = -0.2811;       % [-0.3213, -0.241]
a5 = -0.07965;      % [-0.1109, -0.04839]



F = a0 + exp(-a1.*x+a2) + y.*exp(-(a3*x).^2+a4)+a5.*y;

%% Goodness of fit
% SSE: 1854;
% R-Square: 0.9995; (Same as adjusted)
% RMSE: 2.396