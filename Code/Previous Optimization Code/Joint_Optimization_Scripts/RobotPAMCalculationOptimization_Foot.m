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
    ChooseJoint = 'Foot';
end

%Flag that tells the code to not run pieces that have been calculated
%previously
if exist('beginOptimization', 'var') == 0
    beginOptimization = 0;
end

if beginOptimization == 0
    Raxis = [-0.10501355; -0.17402245; 0.97912632];
    MaxTheta = 20*pi/180;
    MinTheta = -50*pi/180;
    Home = [0; -0.43; 0];                   %Home position from the Tibia
    Joint1 = JointData('Ankle', Raxis, MaxTheta, MinTheta, Home, divisions);

    Raxis = [0.7871796; 0.60474746; -0.12094949];
    MaxTheta = 35*pi/180;
    MinTheta = -25*pi/180;
    Home = [-0.04877; -0.04195; 0.00792];   %Home position from the Talus
    Joint2 = JointData('Subtalar', Raxis, MaxTheta, MinTheta, Home, divisions);

    T_a = Joint1.TransformationMat(:, :, :);
    T_s = Joint2.TransformationMat(:, :, :);

    %Create a Mesh of the Ankle and MTP. First two dimensions are the
    %values of the transformation matrics. The third dimension determines
    %the joint being observed, and the fourth dimension is its location
    %in the mesh, which will need to be converted from a one dimensional
    %array to two dimensions
    T = zeros(4, 4, 3, divisions^2);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            T(:, :, 1, ii) = T_a(:, :, i);            %Keep the ankle at one point as the MTP changes
            T(:, :, 2, ii) = T_s(:, :, j);            
            ii = ii + 1;
        end
    end
end

%Tibialis Posterior t4 -> a4
if beginOptimization == 0
    Location1 = [  -0.02, -0.014, 0.042, 0.077;
                    -0.085, -0.405, 0.033, 0.016;
                    0.002, -0.023, -0.029, -0.028];
    TPCrossPoints = [3 3];
    TPMIF = 3549;  %max isometric force
    Axis1 = [3; 1];                           %The axis of interest when calculating the moment arm about each joint
end
Muscle1 = PamData('Tibialis Posterior', Location1, TPCrossPoints, TPMIF, T, Axis1);

%Tibialis Anterior t2 -> a3
if beginOptimization == 0
    Location2 = [0.006, 0.033, 0.117;
                -0.062, -0.394, 0.014;
                0.047, 0.007, -0.036];    
    TACrossPoints = [3 3];
    TAMIF = 905;  %max isometric force
    Axis2 = [3; 1];                           %The axis of interest when calculating the moment arm about each joint
end
Muscle2 = PamData('Tibialis Posterior', Location2, TACrossPoints, TAMIF, T, Axis2);

%Peroneus Brevis t2 -> a2
if beginOptimization == 0
    Location3 = [0.006, -0.02, -0.014, 0.047, 0.079;
                -0.062, -0.418, -0.429, 0.027, 0.022;
                0.047, 0.028, 0.029, 0.023, 0.034];
    PBCrossPoints = [4 4];
    PBMIF = 435;  %max isometric force
    Axis3 = [3; 1];                           %The axis of interest when calculating the moment arm about each joint
end
Muscle3 = PamData('Peroneus Brevis', Location3, PBCrossPoints, PBMIF, T, Axis3);

%Peroneus Longus t2 -> a3
if beginOptimization == 0
    Location4 = [0.005, -0.021, -0.016, 0.044, 0.048, 0.085, 0.117;
                -0.062, -0.42, -0.432, 0.023, 0.011, 0.007, 0.014;
                0.047, 0.029, 0.029, 0.022, 0.028, 0.012, -0.036];
    PLCrossPoints = [4 4];
    PLMIF = 943;  %max isometric force
    Axis4 = [3; 1];                           %The axis of interest when calculating the moment arm about each joint
end
Muscle4 = PamData('Peroneus Longus', Location4, PLCrossPoints, PLMIF, T, Axis4);

%Peroneus Tertius t2 -> a2
if beginOptimization == 0
    Location5 = [0.006, 0.023, 0.079;
                  -0.062, -0.407, 0.016;
                  0.079, 0.022, 0.034];        
    PTCrossPoints = [3 3];
    PTMIF = 180;  %max isometric force
    Axis5 = [3; 1];                           %The axis of interest when calculating the moment arm about each joint
end
Muscle5 = PamData('Peroneus Tertius', Location5, PTCrossPoints, PTMIF, T, Axis5);

%Torque Calcs, "R" for robot
Torque1 = Muscle1.Torque(1, :) + Muscle2.Torque(1, :) + Muscle3.Torque(1, :) + Muscle4.Torque(1, :) + Muscle5.Torque(1, :);
Torque2 = Muscle1.Torque(2, :) + Muscle2.Torque(2, :) + Muscle3.Torque(2, :) + Muscle4.Torque(2, :) + Muscle5.Torque(2, :);

Torque1M = zeros(divisions, divisions); Torque2M = zeros(divisions, divisions);
%Create the Mesh of Torques to corespond with the joint angles
for i = 1:divisions
    Torque1M(:, i) = Torque1(((i-1)*divisions)+1:i*divisions);
    Torque2M(:, i) = Torque2(((i-1)*divisions)+1:i*divisions);
end

%Including Generic variable name for plotting
RobotAxis1 = Joint1.Theta;
RobotAxis1Label = 'Ankle Flexion, Degrees';
RobotAxis2 = Joint2.Theta;
RobotAxis2Label = 'MTP Flexion, Degrees';
RobotTorque1 = Torque1M;
RobotTorque2 = Torque2M;
RobotTitle1 = 'Robot Ankle Torque, Z Axis';
RobotTitle2 = 'Robot Subtalar Torque, X Axis';


%Store all necessary variables into a cell array for the optimization
%script
Muscles = {Muscle1, Muscle2, Muscle3, Muscle4, Muscle5};