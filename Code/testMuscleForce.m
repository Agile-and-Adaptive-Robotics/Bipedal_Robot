% This script is for playing with the muscle force that can be generated,
% using Hoy, Thelen, and Millard as references

%% Freshen up the workspace
clc
clear
close all

%% Add paths to the muscle and pam calculators
addpath Human_Data
addpath Robot_Data
addpath Functions


%% Joint rotation transformation matrices
positions = 2;
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

Name = 'Adductor Magnus 2';
MIF = 343; 
TSL = 0.12; Pennation = 0.05235988; OFL = 0.121;
Location = [-0.083, -0.119, 0.031;
            0.005, -0.229, 0.023];
CrossPoint = 2;
Add_Mag2 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Adductor Magnus Bi Test';
MIF = 343; 
TSL = 0.12; Pennation = 0.05235988; OFL = 0.121;
T(:, :, :, 2) = T;
Location = [-0.083, -0.119, 0.031;
            0.005, -0.229, 0.023];
CrossPoint = [2 2];
Add = BiMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

