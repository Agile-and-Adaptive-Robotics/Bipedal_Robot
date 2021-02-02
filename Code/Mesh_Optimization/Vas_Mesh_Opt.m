%% Mesh Points Calculation - Bicep Femoris Short Head
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
positions = 100;
fprintf('The algorithm will be calculating Torque at %d different joint positions.\n', positions)

R = zeros(3, 3, positions);
T = zeros(4, 4, positions);

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

for i = 1:positions
    hipToKnee = [fcn1(phi(i)), fcn2(phi(i)), 0];
    R(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    
    T(:, :, i) = RpToTrans(R(:, :, i), hipToKnee');
end

%% Muscle calculation
Name = 'Vastus Intermedius';
MIF = 1365;
OFL = 0.087; TSL = 0.136; Pennation = 0.05235988;
Location = [0.029, -0.192, 0.031;
            0.034, -0.208, 0.029;
            0.034, -0.403, 0.005;
            0.0555, 0.025, 0.0018];
CrossPoint = 4;
Vas_Int = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Vastus Lateralis';
MIF = 1871;
OFL = 0.084; TSL = 0.157; Pennation = 0.08726646;
Location = [0.005, -0.185, 0.035;
            0.027, -0.259, 0.041;
            0.036, -0.403, 0.021;
            0.025, -0.424, 0.018;
            0.06, 0.02, 0.0165];
CrossPoint = 5;
Vas_Lat = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Vastus Medialis';
MIF = 1294;
OFL = 0.089; TSL = 0.126; Pennation = 0.08726646;
Location = [0.014, -0.21, 0.019;
            0.036, -0.277, 0.001;
            0.037, -0.405, -0.013;
            0.027, -0.425, -0.013;
            0.05625, 0.022, -0.0146];
CrossPoint = 5;
Vas_Med = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

%% PAM calculation
Name = 'Vastus Intermedius';
Location = [0.029, -0.192, 0.031;
            0.034, -0.403, 0.005;
            0.0555, 0.025, 0.0018];
CrossPoint = 3;
Dia = 40;
Vas_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);

%% Unstacking the Torques to identify specific rotations
Torque1 = Vas_Int.Torque + Vas_Lat.Torque + Vas_Med.Torque;
TorqueR = Vas_Pam.Torque;

%% Add Torques from the Muscle Group
TorqueH = Torque1;

%% First cost function calculation
C = costFunction(TorqueH, TorqueR);

%% Generating new points for the PAM based on the bone mesh
Femur = xlsread('Femur_Mesh_Points.xlsx');
Tibia = xlsread('Tibia_Mesh_Points.xlsx');

fprintf('The algorithm will be calculating Torque between %d different mesh locations.\n', size(Tibia, 1)*size(Femur, 1))

meshTracker = [0, 0];

originalLocation = Location;

iC = 1;                     %Index variable for the cost function
CMaxPrev = -10^5;

for i = 1:size(Tibia, 1)
    for ii = 1:size(Femur, 1)
        clc
        fprintf('%d \t of %d \n', iC, size(Tibia, 1)*size(Femur, 1))
        Location(1, :) = Femur(ii, :);
        Location(3, :) = Tibia(i, :);
        
        Vas_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);
        TorqueR = Vas_Pam.Torque;
        
        C = costFunction(TorqueH, TorqueR);
        
        if C > CMaxPrev
            if isequal(Vas_Pam.LengthCheck, 'Usable')
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

if exist('Tracker', 'var') == 0
    Location = originalLocation;
else
    Location(1, :) = Femur(Tracker(2), :);
    Location(3, :) = Tibia(Tracker(1), :);
end

%Create variables for the x, y, and z toque
zTorqueHxzRotation = zeros(length(phi));
zTorqueRxzRotation = zeros(length(phi));

Vas_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);
TorqueR = Vas_Pam.Torque;

zTorqueHzRotation = TorqueH(:, 3);
zTorqueRzRotation = TorqueR(:, 3);

yTorqueHzRotation = TorqueH(:, 2);
yTorqueRzRotation = TorqueR(:, 2);

xTorqueHzRotation = TorqueH(:, 1);
xTorqueRzRotation = TorqueR(:, 1);

% oneDofJointTorquePlot

HMuscleLocation = {Vas_Int.Location, Vas_Lat.Location, Vas_Med.Location};
HMuscleCross = {Vas_Int.Cross, Vas_Lat.Cross, Vas_Med.Cross};

RMuscleLocation = {Vas_Pam.Location};
RMuscleCross = {Vas_Pam.Cross};

Bones = {'Femur', 'Tibia'};

run("MuscleBonePlotting")
