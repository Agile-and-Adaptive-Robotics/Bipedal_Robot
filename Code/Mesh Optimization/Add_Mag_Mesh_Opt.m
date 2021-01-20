%% Mesh Points Calculation
% This code will calculate the torque difference between all of the points
% from one bone mesh to another to determine the best location for muscle
% placement

%% Freshen up the workspace
clc
clear
close all

%% Add paths to the muscle and pam calculators
addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Human_Data
addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Robot_Data
addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Functions
addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Bone_Mesh_Plots\Open_Sim_Bone_Geometry

%% Joint rotation transformation matrices
positions = 25;


Rx = zeros(3, 3, positions);
Ry = zeros(3, 3, positions);
Rz = zeros(3, 3, positions);
T = zeros(4, 4, positions);

% Adduction/Abduction, X rotation
adducMax = 20*pi/180;
adducMin = -45*pi/180;
theta = linspace(adducMin, adducMax, positions);

% Since the adductors primarily rotate in the x direction, we don't need to
% observe all the full rotation of y and z direction

% Internal/External rotation, Y rotation
rotationMax = 40*pi/180;
rotationMin = -45*pi/180;
phi = [rotationMin, 0, rotationMax];

% Flexion/Extension, Z rotation
flexMax = 85*pi/180;
flexMin = -25*pi/180;
gamma = [flexMin, 0, flexMax];

pelvisToHip = [-0.0707, -0.0661, 0.0835];

fprintf('The algorithm will be calculating Torque at %d different joint positions.\n', length(theta)*length(phi)*length(gamma))

j = 1;
for iii = 1:length(gamma)
    for ii = 1:length(phi)
        for i = 1:length(theta)
            Rx(:, :, i) = [1, 0, 0;
                           0, cos(theta(i)), -sin(theta(i));
                           0, sin(theta(i)), cos(theta(i))];
                        
            Ry(:, :, ii) = [cos(phi(ii)), 0, -sin(phi(ii));
                            0, 1, 0;
                            sin(phi(ii)), 0, cos(phi(ii))];
                
            Rz(:, :, iii) = [cos(gamma(iii)), -sin(gamma(iii)), 0;
                            sin(gamma(iii)), cos(gamma(iii)), 0;
                            0, 0, 1];
                        
            R = Rz(:, :, iii)*Ry(:, :, ii)*Rx(:, :, i);
            T(:, :, j) = RpToTrans(R, pelvisToHip');
            j = j + 1;
        end
    end
end

%% Muscle calculation
Name = 'Adductor Magnus 1';
MIF = 381;
OFL = 0.087; TSL = 0.06; Pennation = 0.08726646;
Location = [-0.073, -0.117, 0.025;
            -0.004, -0.121, 0.034];
CrossPoint = 2;
Add_Mag1 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Adductor Magnus 2';
MIF = 343;
OFL = 0.121; TSL = 0.12; Pennation = 0.05235988;
Location2 = [-0.083, -0.119, 0.031;
            0.005, -0.229, 0.023];
CrossPoint = 2;
Add_Mag2 = MonoMuscleData(Name, Location2, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Adductor Magnus 3';
MIF = 488;
OFL = 0.131; TSL = 0.249; Pennation = 0.08726646;
Location = [-0.111, -0.114, 0.049;
            0.007, -0.384, -0.027];
CrossPoint = 2;
Add_Mag3 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

%% PAM calculation
Name = 'Adductor Magnus';
CrossPoint = 2;
Dia = 40;

Add_Mag_PamH = MonoPamData(Name, Location2, CrossPoint, Dia, T);        %Pam placed at adductor Magnus 2 location

Location = [-0.083, -0.119, 0.031;
            0.005, -0.229, 0.023];
Add_Mag_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);

%% Unstacking the Torques to identify specific rotations
Torque1 = zeros(length(theta), 3, length(phi), length(gamma));
Torque2 = zeros(length(theta), 3, length(phi), length(gamma));
Torque3 = zeros(length(theta), 3, length(phi), length(gamma));
TorqueR = zeros(length(theta), 3, length(phi), length(gamma));
TorqueRH = zeros(length(theta), 3, length(phi), length(gamma));

PAMTorque = Add_Mag_Pam.Torque;

PAMTorqueH = Add_Mag_PamH.Torque;

j = 1;
for iii = 1:length(gamma)
    for ii = 1:length(phi)
        for i = 1:length(theta)
            Torque1(i, :, ii, iii) = Add_Mag1.Torque(j, :);
            Torque2(i, :, ii, iii) = Add_Mag2.Torque(j, :);
            Torque3(i, :, ii, iii) = Add_Mag3.Torque(j, :);
            TorqueR(i, :, ii, iii) = PAMTorque(j, :);
            TorqueRH(i, :, ii, iii) = PAMTorqueH(j, :);

            j = j + 1;
        end
    end
end

%% Add Torques from the Muscle Group
TorqueH = Torque1 + Torque2 + Torque3;

%% First cost function calculation
C = costFunction(TorqueH, TorqueR);

%% Generating new points for the PAM based on the bone mesh
Pelvis = xlsread('Pelvis_R_Mesh_Points.xlsx');
Femur = xlsread('Femur_Mesh_Points.xlsx');

fprintf('The algorithm will be calculating Torque between %d different mesh locations.\n', size(Pelvis, 1)*size(Femur, 1))

meshTracker = [0, 0];

iC = 1;                     %Index variable for the cost function
CMaxPrev = 10^5;

for i = 1:size(Femur, 1)
    for ii = 1:size(Pelvis, 1)
        clc
        fprintf('%d \t of %d \n', iC, size(Pelvis, 1)*size(Femur, 1))
        
        Location(1, :) = Pelvis(ii, :);
        Location(2, :) = Femur(i, :);
        
        Add_Mag_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);
        PAMTorque = Add_Mag_Pam.Torque;
        
        j = 1;
        for iG = 1:length(gamma)
            for iP = 1:length(phi)
                for iT = 1:length(theta)
                    TorqueR(iT, :, iP, iG) = PAMTorque(j, :);

                    j = j + 1;
                end
            end
        end

        C = costFunction(TorqueH, TorqueR);
        
        if C < CMaxPrev
            if isequal(Add_Mag_Pam.LengthCheck, 'Usable')
                Tracker = [i, ii];
                CMaxPrev = C;
            end
        end
        
        iC = iC + 1;
    end
end

%% Plotting Torque Results
% Home Position

if exist('Tracker', 'var') == 0
    Location = originalLocation;
else
    Location(1, :) = Pelvis(Tracker(2), :);
    Location(2, :) = Femur(Tracker(1), :);
end

Add_Mag_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);
PAMTorque = Add_Mag_Pam.Torque;

j = 1;
for iG = 1:length(gamma)
    for iP = 1:length(phi)
        for iT = 1:length(theta)
            TorqueR(iT, :, iP, iG) = PAMTorque(j, :);

            j = j + 1;
        end
    end
end

gammaD = gamma*180/pi;
phiD = phi*180/pi;
thetaD = theta*180/pi;

% Adjusted error about pure x rotation

TorqueEx = zeros(size(TorqueH, 1), 1);
TorqueEHx = zeros(size(TorqueH, 1), 1);
TorqueEy = zeros(size(TorqueH, 1), 1);
TorqueEHy = zeros(size(TorqueH, 1), 1);
TorqueEz = zeros(size(TorqueH, 1), 1);
TorqueEHz = zeros(size(TorqueH, 1), 1);

for i = 1:size(TorqueR, 1)
    if TorqueH(i, 1) >= 0
        TorqueEx(i) = TorqueR(i, 1, 2, 2) - TorqueH(i, 1, 2, 2);
        TorqueEHx(i) = TorqueRH(i, 1, 2, 2) - TorqueH(i, 1, 2, 2);
    else
        TorqueEx(i) = TorqueH(i, 1, 2, 2) - TorqueR(i, 1, 2, 2);
        TorqueEHx(i) = TorqueH(i, 1, 2, 2) - TorqueRH(i, 1, 2, 2);
    end
    
    if TorqueH(i, 2) >= 0
        TorqueEy(i) = TorqueR(i, 2, 2, 2) - TorqueH(i, 2, 2, 2);
        TorqueEHy(i) = TorqueRH(i, 2, 2, 2) - TorqueH(i, 2, 1, 1);
    else
        TorqueEy(i) = TorqueH(i, 2, 2, 2) - TorqueR(i, 2, 2, 2);
        TorqueEHy(i) = TorqueH(i, 2, 2, 2) - TorqueRH(i, 2, 2, 2);
    end
    
    if TorqueH(i, 3) >= 0
        TorqueEz(i) = TorqueR(i, 3, 2, 2) - TorqueH(i, 3, 2, 2);
        TorqueEHz(i) = TorqueRH(i, 3, 2, 2) - TorqueH(i, 3, 2, 2);
    else
        TorqueEz(i) = TorqueH(i, 3, 2, 2) - TorqueR(i, 3, 2, 2);
        TorqueEHz(i) = TorqueH(i, 3, 2, 2) - TorqueRH(i, 3, 2, 2);
    end
end

%The Adductor will primarily be rotating the hip about the x axis, so we
%will primarily be looking at the theta rotation. We will later look at
%slices of theta as the z direction rotates which is gamma

figure
hold on
sgtitle('Adductor Magnus Torque through Pure Hip Adduction/Abduction')

subplot(3, 2, 1)
plot(thetaD, TorqueH(:, 1, 2, 2), thetaD, TorqueR(:, 1, 2, 2), thetaD, TorqueRH(:, 1, 2, 2))
title('Muscle and PAM X Torque')
xlabel('Pure Hip Adduction/Abduction, degrees')
ylabel('Torque, Nm')
xlim([min(thetaD), max(thetaD)])
legend('Human', 'PAM', 'PAM at Human Locations')

subplot(3, 2, 2)
plot(thetaD, TorqueEx, thetaD, TorqueEHx)
legend('Optimal PAM Location', 'At Human Locations')
xlabel('Pure Hip Adduction/Abduction, degrees')
ylabel('Torque, Nm')
xlim([min(thetaD), max(thetaD)])
title('Adjusted Error X Torque')

subplot(3, 2, 3)
plot(thetaD, TorqueH(:, 2, 2, 2), thetaD, TorqueR(:, 2, 2, 2), thetaD, TorqueRH(:, 2, 2, 2))
title('Muscle and PAM Y Torque')
xlabel('Pure Hip Adduction/Abduction, degrees')
ylabel('Torque, Nm')
xlim([min(thetaD), max(thetaD)])
legend('Human', 'PAM', 'PAM at Human Locations')

subplot(3, 2, 4)
plot(thetaD, TorqueEy, thetaD, TorqueEHy)
legend('Optimal PAM Location', 'At Human Locations')
xlabel('Pure Hip Adduction/Abduction, degrees')
ylabel('Torque, Nm')
xlim([min(thetaD), max(thetaD)])
title('Adjusted Error Y Torque')

subplot(3, 2, 5)
plot(thetaD, TorqueH(:, 3, 2, 2), thetaD, TorqueR(:, 3, 2, 2), thetaD, TorqueRH(:, 3, 2, 2))
legend('Human Muscle', 'Optimal BPA Location', 'BPA at Human Locations')
title('Muscle and PAM Z Torque')
xlabel('Pure Hip Adduction/Abduction, degrees')
ylabel('Torque, Nm')
xlim([min(thetaD), max(thetaD)])
legend('Human', 'PAM', 'PAM at Human Locations', 'Location', 'best')

subplot(3, 2, 6)
plot(thetaD, TorqueEz, thetaD, TorqueEHz)
legend('Optimal PAM Location', 'At Human Locations')
xlabel('Pure Hip Adduction/Abduction, degrees')
ylabel('Torque, Nm')
xlim([min(thetaD), max(thetaD)])
title('Adjusted Error Z Torque')

hold off


% Plotting at the extreme positions of hip flexion/extension

figure
hold on
sgtitle('Adductor Magnus X Axis Torque through Hip Adduction/Abduction at Extreme Hip Positions')

subplot(2, 2, 1)
plot(thetaD, TorqueH(:, 1, 2, 1), thetaD, TorqueR(:, 1, 2, 1), thetaD, TorqueRH(:, 1, 2, 1))
title('Muscle and PAM X Torque at Maximum Extension')
xlabel('Hip Adduction/Abduction, degrees')
ylabel('Torque, Nm')
xlim([min(thetaD), max(thetaD)])
legend('Human', 'PAM', 'PAM at Human Locations')

subplot(2, 2, 2)
plot(thetaD, TorqueH(:, 1, 2, 3), thetaD, TorqueR(:, 1, 2, 3), thetaD, TorqueRH(:, 1, 2, 3))
title('Muscle and Pam X Torque at Maximum Flexion')
xlabel('Hip Adduction/Abduction, degrees')
ylabel('Torque, Nm')
xlim([min(thetaD), max(thetaD)])
legend('Human', 'Optimal PAM Location', 'At Human Locations')

subplot(2, 2, 3)
plot(thetaD, TorqueH(:, 1, 1, 2), thetaD, TorqueR(:, 1, 1, 2), thetaD, TorqueRH(:, 1, 1, 2))
title('Muscle and PAM X Torque at Maximum External Rotation')
xlabel('Hip Adduction/Abduction, degrees')
ylabel('Torque, Nm')
xlim([min(thetaD), max(thetaD)])
legend('Human', 'PAM', 'PAM at Human Locations')

subplot(2, 2, 4)
plot(thetaD, TorqueH(:, 1, 3, 2), thetaD, TorqueR(:, 1, 3, 2), thetaD, TorqueRH(:, 1, 3, 2))
title('Muscle and PAM X Torque at Maximum Internal Rotation')
xlabel('Hip Adduction/Abduction, degrees')
ylabel('Torque, Nm')
xlim([min(thetaD), max(thetaD)])
legend('Human', 'PAM', 'PAM at Human Locations')

hold off

%% Plotting the Angle between vectors

% Calculating the angle between vectors during only abduction and adduction
aHR = zeros(size(TorqueH, 1), 1);
aHRH = zeros(size(TorqueH, 1), 1);

for i = 1:size(TorqueH, 1)
    uvecH = TorqueH(i, :, 2, 2)/norm(TorqueH(i, :, 2, 2));
                
    %Sometimes the BPA can't produce any force due to high
    %contraction. We will set it equal to negative the human
    %vector to maximize the penalty. Consider changing later
    if norm(TorqueR(i, :)) == 0
        uvecR = -uvecH;
    else
        uvecR = TorqueR(i, :, 2, 2)/norm(TorqueR(i, :, 2, 2));
    end
    
    if norm(TorqueRH(i, :)) == 0
        uvecRH = -uvecH;
    else
        uvecRH = TorqueRH(i, :, 2, 2)/norm(TorqueRH(i, :, 2, 2));
    end
    
    aHR(i) = dot(uvecH, uvecR);
    aHRH(i) = dot(uvecH, uvecRH);
end


figure
hold on
title('Angle between the Human Torque Vector and PAM Torque Vectors during Pure Hip Adduction/Abduction')

plot(thetaD, aHR, thetaD, aHRH)
legend('Human and Optimal PAM', 'Human and PAM at Human Locations')
ylabel('Radians')
xlabel('Hip Adduction/Abduction')
hold off


%% Plotting on the Skeleton

HMuscleLocation = {Add_Mag1.Location, Add_Mag2.Location, Add_Mag3.Location};
HMuscleCross = {Add_Mag1.Cross, Add_Mag2.Cross, Add_Mag3.Cross};

RMuscleLocation = {Add_Mag_Pam.Location};
RMuscleCross = {Add_Mag_Pam.Cross};

Bones = {'Pelvis', 'Femur'};

run("Bone_Mesh_Plots\MuscleBonePlotting")






% %% Plotting Torque Results
% % This looks at the x, y, and z torque when rotating the hip through the x
% % and z axis. 
% 
% %Set up a mesh for x and y coordinates on the plot
% [mTheta, mGamma] = meshgrid(theta, gamma);
% 
% %Create variables for the x, y, and z toque
% xTorqueHxzRotation = zeros(length(gamma), length(theta));
% xTorqueRxzRotation = zeros(length(gamma), length(theta));
% 
% yTorqueHxzRotation = zeros(length(gamma), length(theta));
% yTorqueRxzRotation = zeros(length(gamma), length(theta));
% 
% zTorqueHxzRotation = zeros(length(gamma), length(theta));
% zTorqueRxzRotation = zeros(length(gamma), length(theta));
% 
% Location(1, :) = Pelvis(Tracker(2), :);
% Location(2, :) = Femur(Tracker(1), :);
% 
% Add_Mag_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);
% PAMTorque = Add_Mag_Pam.Torque;
% 
% j = 1;
% for iG = 1:length(gamma)
%     for iP = 1:length(phi)
%         for iT = 1:length(theta)
%             TorqueR(iT, :, iP, iG) = PAMTorque(j, :);
% 
%             j = j + 1;
%         end
%     end
% end
% 
% for iii = 1:length(gamma)
%     for i = 1:length(theta)
%         xTorqueHxzRotation(iii, i) = TorqueH(i, 1, 1, iii);
%         xTorqueRxzRotation(iii, i) = TorqueR(i, 1, 1, iii);
%         
%         yTorqueHxzRotation(iii, i) = TorqueH(i, 2, 1, iii);
%         yTorqueRxzRotation(iii, i) = TorqueR(i, 2, 1, iii);
%         
%         zTorqueHxzRotation(iii, i) = TorqueH(i, 3, 1, iii);
%         zTorqueRxzRotation(iii, i) = TorqueR(i, 3, 1, iii);
%     end
% end
% 
% testAdductorOptimizationPlot
% 

