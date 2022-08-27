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
    
%     hipToKnee_Pam = [0.0045, -0.3958, 0];
%     R_Pam(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;   %Rotation matrix for robot
%                     sin(phi(i)), cos(phi(i)), 0;
%                     0, 0, 1];
%     
%     T_Pam(:, :, i) = RpToTrans(R_Pam(:, :, i), hipToKnee_Pam');     %Transformation matrix for robot
end


%% PAM calculation
Name = 'Bicep Femoris (Short Head)';
CrossPoint = 2;
Dia = 10;
Location = zeros(2,3,positions);
%Origin and Insertion from Assem2.75 Solidworks assembly
for i = 1:positions
    Location(:,:,i) = [-0.075, 0.100, 0.0328;
            -0.05011, -0.045, 0.0328];
end
rest = 0.457; %resting length, m
kmax = 0.380;     %Length at maximum contraction, m
tendon = 0;       %pinned, no tendon
fitting = 0.0254; %fitting length
Pres = 604.8;     %average pressure
Bifemsh_Pam = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, Pres);



% figure
% plot(phi,Bifemsh_Pam.MuscleLength)  %Length including fittings
% 
% max_length = max(Bifemsh_Pam.MuscleLength)
% min_length = min(Bifemsh_Pam.MuscleLength)
% ratio = min_length/max_length


%% Unstacking the Torques to identify specific rotations

TorqueR = Bifemsh_Pam.Torque;


%% Plotting Torque Results
phiD = phi*180/pi;


%% Plotting muscle lengths and moment arms using two different moment arm
%calculations
% ML = Bifemsh.MuscleLength;
% PamL = Bifemsh_Pam.MuscleLength;
% for i = 1:size(Bifemsh.MomentArm,1)
%     MA(i,:) = -norm(Bifemsh.MomentArm(i,1:2));               %Muscle moment arm, Z axis
%     BPAma(i,:) = -norm(Bifemsh_Pam.MomentArm(i,1:2));        %BPA moment arm, Z axis
% end
% dM = diff(Bifemsh.MuscleLength);           %Muscle length difference
% dP = diff(Bifemsh_Pam.MuscleLength);       %PAM length difference
% dO = diff(phiD);                           %Angle difference
% 
% figure
% hold on
% sgtitle('Bicep Femoris Short Head Length and Moment Arm through Knee Flexion and Extension')
% 
% subplot(2, 2, 1)
% plot(phiD, ML, phiD, PamL)
% title('Muscle and PAM Lengths')
% xlabel('Knee Ext(+)/Flx(-), degrees')
% ylabel('Length, m')
% legend('Human', 'PAM')
% 
% subplot(2, 2, 2)
% plot(phiD, MA, phiD, BPAma)
% title('Moment arm, Z axis, vector method')
% xlabel('Knee Ext(+)/Flx(-), degrees')
% ylabel('Length, m')
% legend('Human', 'PAM')
% 
% subplot(2, 2, 3)
% plot(phiD(1:99), -dM./dO', phiD(1:99), -dP./dO')
% title('Moment arm, Z axis, left difference method')
% xlabel('Knee Ext(+)/Flx(-), degrees')
% ylabel('Length, m')
% legend('Human', 'PAM')
% 
% subplot(2, 2, 4)
% plot(phiD(2:100), -dM./dO', phiD(2:100), -dP./dO')
% title('Moment arm, Z axis, right difference method')
% xlabel('Knee Ext(+)/Flx(-), degrees')
% ylabel('Length, m')
% legend('Human', 'PAM')
% 
% hold off


%% Plot just robot Z axis Torque
figure
plot(phiD, TorqueR(:, 3))
title('BPA Z Torque, Length = 457 mm')
xlabel('Knee Extension(+)/Flexion(-), degrees')
ylabel('Torque, Nm')
legend('BPA')
