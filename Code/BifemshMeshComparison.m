%% Mesh Points Calculation - Bicep Femoris Short Head
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
Name = 'Bicep Femoris (Short Head)';
MIF = 804;
OFL = 0.173; TSL = 0.089; Pennation = 0.40142573;
Location = [0.005, -0.211, 0.023;
            -0.03, -0.036, 0.029;
            -0.023, -0.056, 0.034];
CrossPoint = 2;
Bifemsh = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

%% PAM calculation
Name = 'Bicep Femoris (Short Head)';
Location = [0.005, -0.211, 0.023;
            -0.023, -0.056, 0.034];
CrossPoint = 2;
Dia = 20;
Bifemsh_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);

%% Unstacking the Torques to identify specific rotations
Torque1 = Bifemsh.Torque;
TorqueR = Bifemsh_Pam.Torque;

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
        Location(2, :) = Tibia(i, :);
        
        Bifemsh_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);
        TorqueR = Bifemsh_Pam.Torque;
        
        C = costFunction(TorqueH, TorqueR);
        
        if C > CMaxPrev
            if isequal(Bifemsh_Pam.LengthCheck, 'Usable')
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
    Location(2, :) = Tibia(Tracker(1), :);
end

%Create variables for the x, y, and z toque
zTorqueHxzRotation = zeros(length(phi));
zTorqueRxzRotation = zeros(length(phi));

Bifemsh_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);
TorqueR = Bifemsh_Pam.Torque;

zTorqueHzRotation = TorqueH(:, 3);
zTorqueRzRotation = TorqueR(:, 3);

yTorqueHzRotation = TorqueH(:, 2);
yTorqueRzRotation = TorqueR(:, 2);

xTorqueHzRotation = TorqueH(:, 1);
xTorqueRzRotation = TorqueR(:, 1);

oneDofJointTorquePlot

HMuscleLocation = {Bifemsh.Location};
HMuscleCross = {Bifemsh.Cross};

RMuscleLocation = {Bifemsh_Pam.Location};
RMuscleCross = {Bifemsh_Pam.Cross};

run("Bone_Mesh_Plots\MuscleBonePlotting")
