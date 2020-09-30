% Robot PAM Calculation Optimization
% Author: Connor Morrow
% Date: 09/30/2020
% Description: This script calculate the torque generated across the back joint by
% the PAMs that cross over them. It is intended to be used by the 
% optimization script.

if exist('divisions', 'var') == 0
    divisions = 100;
end

if exist('ChooseJoint', 'var') == 0
    ChooseJoint = 'Uni_Hip';
end

%Flag that tells the code to not run pieces that have been calculated
%previously
if exist('beginOptimization', 'var') == 0
    beginOptimization = 0;
end

if beginOptimization == 0
    %Hip Joint x rotation
    Raxis = [1; 0; 0];
    MaxTheta = 20.6*pi/180;
    MinTheta = -30*pi/180;
    HomeH = [-0.0707; -0.0661; 0.0835];                   %Home position from the Tibia
    Joint1a = JointData('HipX', Raxis, MaxTheta, MinTheta, HomeH, divisions);  

    %Hip Joint z rotation
    Raxis = [0; 0; 1];
    MaxTheta = 110*pi/180;
    MinTheta = -15*pi/180;
    HomeH = [-0.0707; -0.0661; 0.0835];                   %Home position from the Tibia
    Joint1b = JointData('HipZ', Raxis, MaxTheta, MinTheta, HomeH, divisions);

    %Back rotation
    Raxis = [0; 0; 1];
    MaxTheta = 0;
    MinTheta = 0;
    HomeB = [-0.1007; 0.0815; 0];                   %Home position from the Tibia
    Joint2 = JointData('Back', Raxis, MaxTheta, MinTheta, HomeB, divisions);

    R = zeros(3, 3, 1, divisions^2);
    T = zeros(4, 4, 1, divisions^2);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            R(:, :, 1, ii) = Joint1a.RotMat(:, :, j)*Joint1b.RotMat(:, :, i);
            T(:, :, 1, ii) = [R(:, :, 1, ii), HomeH; 0 0 0 1];
            ii = ii+1;
        end
    end
end

%Gluteus Maximus p6 -> f3
if beginOptimization == 0
    Location1 = [-0.162, -0.174, -0.069, -0.046, 0.005;
                  0.04, -0.033, -0.042, -0.116, -0.339;
                  0.023, 0.101, 0.067, 0.06, 0.022];
    GMCrossPoints = 3;
    GMMIF = 5047;  
    Axis1 = [10, 20, 30];
end
Muscle1 = PamData('Gluteus Maximus', Location1, GMCrossPoints, GMMIF, T, Axis1);

%Adductor Magnus p9 -> f4
if beginOptimization == 0
    Location2 = [-0.163, -0.059, -0.059, 0.015;
                  -0.013, -0.108, -0.108, -0.374;
                  0.005, -0.03, -0.03, -0.025];
    AMCrossPoints = 3;
    AMMIF = 2268;  
    Axis1 = [10, 20, 30];
end
Muscle2 = PamData('Adductor Magnus', Location2, AMCrossPoints, AMMIF, T, Axis1);

%Iliacus p1 -> f4
if beginOptimization == 0
    Location3 = [-0.055, -0.024, 0.015;
                0.091, -0.057, -0.374;
                0.085, 0.076, -0.025];         
    ICrossPoints = 3;
    IMIF = 1073;  
    Axis1 = [10, 20, 30];
end
Muscle3 = PamData('Iliacus', Location3, ICrossPoints, IMIF, T, Axis1);

%Include the back to the transformtion matrix for the Psoas
if beginOptimization == 0
    T(:, :, 2, :) = T(:, :, 1, :);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            T(:, :, 1, ii) = Joint2.TransformationMat(:, :, j);
            ii = ii+1;
        end
    end
end

%Psoas b1 -> f4
if beginOptimization == 0
    Location4 = [0.036, -0.024, 0.015;
                  0.007, -0.057, -0.374;
                  0.029, 0.076, -0.025];
    PCrossPoints = [2 3];
    PMIF = 1113;
    Axis1 = [10 , 20, 30;
            10, 20, 30];
end
Muscle4 = PamData('Psoas X', PLocation, PCrossPoints, PMIF, T, Axis1);

HipXTorque = Muscle1.Torque(1, :, 1) + Muscle2.Torque(1, :, 1) + Muscle3.Torque(1, :, 1) + Muscle4.Torque(2, :, 1);
HipYTorque = Muscle1.Torque(1, :, 2) + Muscle2.Torque(1, :, 2) + Muscle3.Torque(1, :, 2) + Muscle4.Torque(2, :, 2);
HipZTorque = Muscle1.Torque(1, :, 3) + Muscle2.Torque(1, :, 3) + Muscle3.Torque(1, :, 3) + Muscle4.Torque(2, :, 3);

%Create the Mesh of Torques to corespond with the joint angles
HipXTorqueMR = zeros(divisions, divisions); HipYTorqueMR = zeros(divisions, divisions); HipZTorqueMR = zeros(divisions, divisions);
for i = 1:divisions
    HipXTorqueMR(:, i) = HipXTorque(((i-1)*divisions)+1:i*divisions);
    HipYTorqueMR(:, i) = HipYTorque(((i-1)*divisions)+1:i*divisions);
    HipZTorqueMR(:, i) = HipZTorque(((i-1)*divisions)+1:i*divisions);
end

%Including Generic variable name for plotting
RobotAxis1 = Joint1b.Theta;
RobotAxis1Label = 'Hip Flexion, Degrees';
RobotAxis2 = Joint1a.Theta;
RobotAxis2Label = 'Adduction/Abduction, Degrees';
RobotTorque1 = HipZTorqueMR;
RobotTorque2 = HipYTorqueMR;
RobotTorque3 = HipXTorqueMR;
RobotTitle1 = 'Robot Hip Torque, uniarticular, Z axis';
RobotTitle2 = 'Robot Hip Torque, uniarticular, Y axis';
RobotTitle3 = 'Robot Hip Torque, uniarticular, X axis';


%Store all necessary variables into a cell array for the optimization
%script
Muscles = {Muscle1, Muscle2, Muscle3, Muscle4};