%% Equation for Fmax with Ben and Lawrence's recorded data for 10mm BPA.
%Ben has five 10 mm ID BPAs with force measurements at maximum pressure and
%many relative strains. He also has many BPAs measured at zero strain and
%various pressures. 
%Lawrence has data on five BPAs with high resolution of pressure data but
%only at a few different amounts of relative strain

%Inputs:
%   - Resting length (in meters)
%   - Pressure (in kPa); if no pressure given, use 620;
%Outputs:
%   - Force (Newtons)

function z = Fmax20(restingLength, Pressure)

if nargin ==1 
    Pressure = 620;
end

x = restingLength;
y = Pressure;
%Parameters for Fmax equation
a1 = 1.4877;                    %N/kPa
a2 = 0.0248;                   %1/kPa
% a = [a1 a2];
% a_ub = [1.5011 0.0258];      %Upper bound on coefficient confidence interval
% a_lb = [1.4745 0.0238];      %Lower bound on coefficient confidence interval
                             %Maximum Pressure

z = y*(a1*atan(a2*(x-0.0075)*y));
end