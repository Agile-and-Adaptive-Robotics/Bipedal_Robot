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
positions = 2;
fprintf('The algorithm will be calculating Torque at %d different joint positions.\n', positions^6)

Rx = zeros(3, 3, positions);
Ry = zeros(3, 3, positions);
Rz = zeros(3, 3, positions);
T = zeros(4, 4, positions);

% Hip rotations
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

%Knee Extension and Flexion
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');

kneeMin = -2.0943951;
kneeMax = 0.17453293;
phi = linspace(kneeMin, kneeMax, positions);

iAngle = 1;             %Tracking location of the knee joint
    
iT = 1;                 %Tracking location to place in the transformation Matrix

for iii = 1:length(gamma)
    for ii = 1:length(theta)
        for i = 1:length(phi)
            hipToKnee = [fcn1(phi(iAngle)), fcn2(phi(iAngle)), 0];
            R(:, :, iT, 2) = [cos(theta(iAngle)), -sin(theta(iAngle)), 0;
                            sin(theta(iAngle)), cos(theta(iAngle)), 0;
                            0, 0, 1];

            T(:, :, iT, 2) = RpToTrans(R(:, :, iT, 2), hipToKnee');
            
            iT = iT + 1;
            
            iAngle = iAngle + 1;
            if iAngle > length(phi)
                iAngle = 1;
            end
        end
    end
end

%% Muscle calculation
Name = 'Bicep Femoris (Long Head)';
MIF = 896;
OFL = 0.109; TSL = 0.326; Pennation = 0.0;
Location = [-0.126, -0.103, 0.069;
            -0.03, -0.036, 0.029;
            -0.023, -0.056, 0.034];
CrossPoint = [2, 2];
Bifemlh = BiMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);


%% PAM calculation
Name = 'Bicep Femoris (Long Head)';
Location = [-0.126, -0.103, 0.069;
            -0.023, -0.056, 0.034];
CrossPoint = [2, 2];
Dia = 20;
Bifemlh_Pam = BiPamData(Name, Location, CrossPoint, Dia, T);

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

iC = 1;                     %Index variable for the cost function
CMaxPrev = -10^5;

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
        
        if C > CMaxPrev
            if isequal(Add_Mag_Pam.LengthCheck, 'Usable')
                Tracker = [i, ii];
                CMaxPrev = C;
            end
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

HMuscleLocation = {Add_Mag1.Location, Add_Mag2.Location, Add_Mag3.Location};
HMuscleCross = {Add_Mag1.Cross, Add_Mag2.Cross, Add_Mag3.Cross};

RMuscleLocation = {Add_Mag_Pam.Location};
RMuscleCross = {Add_Mag_Pam.Cross};

Bones = {'Pelvis', 'Femur'};

run("Bone_Mesh_Plots\MuscleBonePlotting")
