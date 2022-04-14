%% Mesh Points Calculation - Bicep Femoris Short Head
% This code will calculate the torque difference between all of the points
% from one bone mesh to another to determine the best location for muscle
% placement

%% Freshen up the workspace
clc
clear
close all

%% Add paths to the muscle and pam calculators
current_dir = cd;
all_code = fullfile(current_dir,'../..');
addpath(genpath(all_code));

% addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Human_Data
% addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Robot_Data
% addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Functions
% addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Bone_Mesh_Plots\Open_Sim_Bone_Geometry
% addpath
% C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Mesh_Optimization\results

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

%% Muscle calculation
Name = 'Bicep Femoris (Short Head)';
MIF = 804;
OFL = 0.173; TSL = 0.089; Pennation = 0.40142573;
Location = zeros(3,3,positions);
for i = 1:positions
    Location(:,:,i) = [0.005, -0.211, 0.023;
            -0.03, -0.036, 0.029;
            -0.023, -0.056, 0.034];
end
CrossPoint = 2;
Bifemsh = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

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



figure
plot(phi,Bifemsh_Pam.MuscleLength)  %Length including fittings

max_length = max(Bifemsh_Pam.MuscleLength)
min_length = min(Bifemsh_Pam.MuscleLength)
ratio = min_length/max_length


%% Unstacking the Torques to identify specific rotations
Torque1 = Bifemsh.Torque;
TorqueR = Bifemsh_Pam.Torque;

%% Add Torques from the Muscle Group
TorqueH = Torque1;

%% Plotting Torque Results
phiD = phi*180/pi;

TorqueEx = zeros(size(TorqueH, 1), 1);
TorqueEy = zeros(size(TorqueH, 1), 1);
TorqueEz = zeros(size(TorqueH, 1), 1);

for i = 1:size(TorqueR, 1)
    if TorqueH(i, 1) >= 0
        TorqueEx(i) = TorqueR(i, 1) - TorqueH(i, 1);
    else
        TorqueEx(i) = TorqueH(i, 1) - TorqueR(i, 1);
    end
    
    if TorqueH(i, 2) >= 0
        TorqueEy(i) = TorqueR(i, 2) - TorqueH(i, 2);
    else
        TorqueEy(i) = TorqueH(i, 2) - TorqueR(i, 2);
    end
    
    if TorqueH(i, 3) >= 0
        TorqueEz(i) = TorqueR(i, 3) - TorqueH(i, 3);
    else
        TorqueEz(i) = TorqueH(i, 3) - TorqueR(i, 3);
    end
end

figure
hold on
sgtitle('Bicep Femoris Short Head Torque through Knee Flexion and Extension')

subplot(3, 2, 1)
plot(phiD, TorqueH(:, 3), phiD, TorqueR(:, 3))
title('Muscle and PAM Z Torque')
xlabel('Knee Ext(+)/Flx(-), degrees')
ylabel('Torque, Nm')
legend('Human', 'PAM')

subplot(3, 2, 2)
plot(phiD, TorqueEz)
legend('Optimal PAM Location')
xlabel('Knee Ext(+)/Flx(-), degrees')
ylabel('Torque, Nm')
title('Adjusted Error Z Torque')

subplot(3, 2, 3)
plot(phiD, TorqueH(:, 2), phiD, TorqueR(:, 2))
title('Muscle and PAM Y Torque')
xlabel('Knee Ext(+)/Flx(-), degrees')
ylabel('Torque, Nm')
legend('Human', 'PAM')

subplot(3, 2, 4)
plot(phiD, TorqueEy)
legend('Optimal PAM Location')
xlabel('Knee Ext(+)/Flx(-), degrees')
ylabel('Torque, Nm')
title('Adjusted Error Y Torque')

subplot(3, 2, 5)
plot(phiD, TorqueH(:, 1), phiD, TorqueR(:, 1))
title('Muscle and PAM X Torque')
xlabel('Knee Ext(+)/Flx(-), degrees')
ylabel('Torque, Nm')
legend('Human', 'PAM')

subplot(3, 2, 6)
plot(phiD, TorqueEx)
legend('Optimal PAM Location')
xlabel('Knee Ext(+)/Flx(-), degrees')
ylabel('Torque, Nm')
title('Adjusted Error X Torque')

hold off

%% Plotting muscle lengths and moment arms using two different moment arm
%calculations
ML = Bifemsh.MuscleLength;
PamL = Bifemsh_Pam.MuscleLength;
for i = 1:size(Bifemsh.MomentArm,1)
    MA(i,:) = -norm(Bifemsh.MomentArm(i,1:2));               %Muscle moment arm, Z axis
    BPAma(i,:) = -norm(Bifemsh_Pam.MomentArm(i,1:2));        %BPA moment arm, Z axis
end
dM = diff(Bifemsh.MuscleLength);           %Muscle length difference
dP = diff(Bifemsh_Pam.MuscleLength);       %PAM length difference
dO = diff(phiD);                           %Angle difference

figure
hold on
sgtitle('Bicep Femoris Short Head Length and Moment Arm through Knee Flexion and Extension')

subplot(2, 2, 1)
plot(phiD, ML, phiD, PamL)
title('Muscle and PAM Lengths')
xlabel('Knee Ext(+)/Flx(-), degrees')
ylabel('Length, m')
legend('Human', 'PAM')

subplot(2, 2, 2)
plot(phiD, MA, phiD, BPAma)
title('Moment arm, Z axis, vector method')
xlabel('Knee Ext(+)/Flx(-), degrees')
ylabel('Length, m')
legend('Human', 'PAM')

subplot(2, 2, 3)
plot(phiD(1:99), -dM./dO', phiD(1:99), -dP./dO')
title('Moment arm, Z axis, left difference method')
xlabel('Knee Ext(+)/Flx(-), degrees')
ylabel('Length, m')
legend('Human', 'PAM')

subplot(2, 2, 4)
plot(phiD(2:100), -dM./dO', phiD(2:100), -dP./dO')
title('Moment arm, Z axis, right difference method')
xlabel('Knee Ext(+)/Flx(-), degrees')
ylabel('Length, m')
legend('Human', 'PAM')

hold off
%% Plotting the angle between the vectors

aHR = zeros(size(TorqueH, 1), 1);
aHRH = zeros(size(TorqueH, 1), 1);

for i = 1:size(TorqueH, 1)
    uvecH = TorqueH(i, :)/norm(TorqueH(i, :));
                
    %Sometimes the BPA can't produce any force due to high
    %contraction. We will set it equal to negative the human
    %vector to maximize the penalty. Consider changing later
    if norm(TorqueR(i, :)) == 0
        uvecR = -uvecH;
    else
        uvecR = TorqueR(i, :)/norm(TorqueR(i, :));
    end
    
    aHR(i) = dot(uvecH, uvecR);
end


figure
hold on
title('Angle between the Human Torque Vector and PAM Torque Vectors')
plot(phiD, aHR)
legend('Human and Optimal PAM')
ylabel('Radians')
xlabel('Knee Angle, degree')
hold off

%% Plotting on the Mesh Skeleton

HMuscleLocation = {Bifemsh.Location};
HMuscleCross = {Bifemsh.Cross};

RMuscleLocation = {Bifemsh_Pam.Location};
RMuscleCross = {Bifemsh_Pam.Cross};

Bones = {'Femur', 'Tibia'};

run("MuscleBonePlotting")

%% Plot just robot Z axis Torque
figure
plot(phiD, TorqueR(:, 3))
title('BPA Z Torque, Length = 485 mm')
xlabel('Knee Extension(+)/Flexion(-), degrees')
ylabel('Torque, Nm')
legend('BPA')
