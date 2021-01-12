%% Mesh Points Calculation
% This code will calculate the torque difference between all of the points
% from one bone mesh to another to determine the best location for muscle
% placement

%% Freshen up the workspace
clc
clear
close all

%% Add paths to the muscle and pam calculators
addpath Human_Data
addpath Robot_Data
addpath Functions
addpath Bone_Mesh_Plots\Open_Sim_Bone_Geometry

%% Joint rotation transformation matrices
positions = 5;
fprintf('The algorithm will be calculating Torque at %d different joint positions.\n', positions*positions*positions)

Rx = zeros(3, 3, positions);
Ry = zeros(3, 3, positions);
Rz = zeros(3, 3, positions);
T = zeros(4, 4, positions);

% Adduction/Abduction, X rotation
adducMax = 20*pi/180;
adducMin = -45*pi/180;
theta = linspace(adducMin, adducMax, positions);

% Internal/External rotation, Y rotation
rotationMax = 40*pi/180;
rotationMin = -45*pi/180;
phi = linspace(rotationMin, rotationMax, positions);

% Flexion/Extension, Z rotation
flexMax = 85*pi/180;
flexMin = -25*pi/180;
gamma = linspace(flexMin, flexMax, positions);

pelvisToHip = [-0.0707, -0.0661, 0.0835];

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
Location = [-0.083, -0.119, 0.031;
            0.005, -0.229, 0.023];
CrossPoint = 2;
Add_Mag2 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Adductor Magnus 3';
MIF = 488;
OFL = 0.131; TSL = 0.249; Pennation = 0.08726646;
Location = [-0.111, -0.114, 0.049;
            0.007, -0.384, -0.027];
CrossPoint = 2;
Add_Mag3 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

%% PAM calculation
Name = 'Adductor Magnus';
Location = [-0.083, -0.119, 0.031;
            0.005, -0.229, 0.023];
CrossPoint = 2;
Dia = 40;
Add_Mag_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);

%% Unstacking the Torques to identify specific rotations
Torque1 = zeros(length(theta), 3, length(phi), length(gamma));
Torque2 = zeros(length(theta), 3, length(phi), length(gamma));
Torque3 = zeros(length(theta), 3, length(phi), length(gamma));
TorqueR = zeros(length(theta), 3, length(phi), length(gamma));

PAMTorque = Add_Mag_Pam.Torque;

j = 1;
for iii = 1:length(gamma)
    for ii = 1:length(phi)
        for i = 1:length(theta)
            Torque1(i, :, ii, iii) = Add_Mag1.Torque(j, :);
            Torque2(i, :, ii, iii) = Add_Mag2.Torque(j, :);
            Torque3(i, :, ii, iii) = Add_Mag3.Torque(j, :);
            TorqueR(i, :, ii, iii) = PAMTorque(j, :);

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
meshTracker = [0, 0];

iC = 2;                     %Index variable for the cost function

for i = 1:size(Femur, 1)
    for ii = 1:size(Pelvis, 1)
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

        C(iC) = costFunction(TorqueH, TorqueR);
        
        if C(iC) == max(C)
            Tracker = [i, ii];
        end

        iC = iC + 1;
    end
end

%% Plotting Torque Results
% This looks at the x, y, and z torque when rotating the hip through the x
% and z axis. 

%Set up a mesh for x and y coordinates on the plot
[mTheta, mGamma] = meshgrid(theta, gamma);

%Create variables for the x, y, and z toque
xTorqueHxzRotation = zeros(length(gamma), length(theta));
xTorqueRxzRotation = zeros(length(gamma), length(theta));

yTorqueHxzRotation = zeros(length(gamma), length(theta));
yTorqueRxzRotation = zeros(length(gamma), length(theta));

zTorqueHxzRotation = zeros(length(gamma), length(theta));
zTorqueRxzRotation = zeros(length(gamma), length(theta));

Location(1, :) = Pelvis(Tracker(2), :);
Location(2, :) = Femur(Tracker(1), :);

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

for iii = 1:length(gamma)
    for i = 1:length(theta)
        xTorqueHxzRotation(iii, i) = TorqueH(i, 1, 1, iii);
        xTorqueRxzRotation(iii, i) = TorqueR(i, 1, 1, iii);
        
        yTorqueHxzRotation(iii, i) = TorqueH(i, 2, 1, iii);
        yTorqueRxzRotation(iii, i) = TorqueR(i, 2, 1, iii);
        
        zTorqueHxzRotation(iii, i) = TorqueH(i, 3, 1, iii);
        zTorqueRxzRotation(iii, i) = TorqueR(i, 3, 1, iii);
    end
end

testAdductorOptimizationPlot

run("Bone_Mesh_Plots\MuscleBonePlotting")
