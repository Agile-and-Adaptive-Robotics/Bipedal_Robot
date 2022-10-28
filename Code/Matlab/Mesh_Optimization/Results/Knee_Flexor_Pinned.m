%% Mesh Points Calculation - Bicep Femoris Short Head
% This code will calculate the torque for the pinned knee joint robot leg

%% Freshen up the workspace
clc
clear
close all

%% Add paths to the muscle and pam calculators
current_dir = cd;
all_code = fullfile(current_dir,'../..');
addpath(genpath(all_code));

%% Joint rotation transformation matrices
positions = 100;
fprintf('The algorithm will be calculating Torque at %d different joint positions.\n', positions)

R = zeros(3, 3, positions);
T = zeros(4, 4, positions);

kneeMin = -2.0943951;
kneeMax = 0.17453293;
phi = linspace(kneeMin, kneeMax, positions);
%We want one of our positions to be home position, so let's make the
%smallest value of phi equal to 0
[val, pos] = min(abs(phi));
phi(pos) = 0;

for i = 1:positions
    hipToKnee = [0.0045, -0.3958, 0];
    R(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    
    T(:, :, i) = RpToTrans(R(:, :, i), hipToKnee');
    
end


%% PAM calculation
Name = 'Bicep Femoris (Short Head)';
CrossPoint = 2;
Dia = 10;
Location = zeros(2,3,positions);
%Origin and Insertion from Assem2.75 Solidworks assembly
for i = 1:positions
    Location(:,:,i) = [-0.0795, 0.100, 0.0328;
                    -0.05011, -0.045, 0.0328];
end
rest = 0.457; %resting length, m
kmax = 0.380;     %Length at maximum contraction, m
tendon = 0;       %pinned, no tendon
fitting = 0.020; %fitting length
%Pres = 605.6671;     %average pressure for 49cm
Pres = 604.5573;      %average pressure for 46cm
Bifemsh_Pam = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, Pres);


%% Unstacking the Torques to identify specific rotations

TorqueR = Bifemsh_Pam.Torque;


%% Plotting Torque Results
phiD = phi*180/pi;



%% Plot just robot Z axis Torque
figure
plot(phiD, TorqueR(:, 3))
title(sprintf('BPA Z Torque, Length = %1.3  mm',rest))
xlabel('Knee Extension(+)/Flexion(-), degrees')
ylabel('Torque, Nm')
hold on
for i = 2:size(TorqueR,3)
    plot(phiD, TorqueR(:,3,i))
end
legend('Hunt','Exponential','Polynomial','Exponential, simplified','Polynomial, simplified')
