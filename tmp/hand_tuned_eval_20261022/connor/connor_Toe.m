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
    ChooseJoint = 'Toe';
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
    MaxTheta = 0;
    MinTheta = 0;               %Ben's script lists the min and max, however only 0 is used. Why?
    Home = [-0.04877; -0.04195; 0.00792];   %Home position from the Talus
    Joint2 = JointData('Subtalar', Raxis, MaxTheta, MinTheta, Home, divisions);

    Raxis = [-0.5809544; 0; 0.81393611];
    MaxTheta = 80*pi/180;
    MinTheta = -30*pi/180;
    Home = [0.1788; -0.002; 0.00108];       %Home Position from the Calcn
    Joint3 = JointData('MTP', Raxis, MaxTheta, MinTheta, Home, divisions);

    T_a = Joint1.TransformationMat(:, :, :);
    T_s = Joint2.TransformationMat(:, :, :);
    T_m = Joint3.TransformationMat(:, :, :);

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
            T(:, :, 2, ii) = T_s(:, :, i);            %Subtalar doesn't change
            T(:, :, 3, ii) = T_m(:, :, j);            %Change the MTP while the subtalar stays the same
            ii = ii + 1;
        end
    end
end

% Creating PAM class module
%Looking at the flexor digitorus longus, but now for the bipedal robot
%The robot uses the same physical skeletal geometry, but changes where
%the PAMs are attached when compared to muscles. Because of this,
%transformation matrices from the skeletal structure will be reused
%from above. 

%For the robot, the wrapping locations also include a global origin
%point (in this case, t4) and a global insertion point (in this case,
%a7).


%FDL Insertion points, found from AttachPoints_Robot. The first column is
%the insertion point that shares a point with the tibia (t4). The last
%point is a7, the insertion point for the foot
if beginOptimization == 0
    Location1 = [-0.02, -0.015, 0.044, 0.071, 0.166, -0.002, 0.028, 0.044;
                    -0.085, -0.405, 0.032, 0.018, -0.008, -0.008, -0.007, -0.006;
                    0.002, -0.02, -0.028, -0.026, 0.012, 0.015, 0.022, 0.024];    

    Axis1 = [3; 1; 3];                           %The axis of interest when calculating the moment arm about each joint
    FDLCrossPoints = [3, 3, 6];
    FDLMIF = 310;  %max isometric force
end
Muscle1 = PamData('Flexor Digitorum Longus', Location1, FDLCrossPoints, FDLMIF, T, Axis1);


%FHL t4 to a8
if beginOptimization == 0
    Location2 = [-0.02, -0.019, 0.037, 0.104, 0.173, 0.016, 0.056;
                    -0.085, -.408, 0.028, 0.007, -0.005, -0.006, -0.01;
                    0.002, -0.017, -0.024, -0.026, -0.027, -0.026, -0.018];
    FHLCrossPoints = [3, 3, 6];
    FHLMIF = 322;  %max isometric force
    Axis2 = [3; 1; 3];                           %The axis of interest when calculating the moment arm about each joint
end
Muscle2 = PamData('Flexor Hallucis Longus', Location2, FHLCrossPoints, FHLMIF, T, Axis1);

%EDL, t2 to a5
if beginOptimization == 0
    Location3 = [0.006, 0.029, 0.092, 0.162, 0.005, 0.044;
        -0.062, -0.401, 0.039, 0.006, 0.005, 0.0;
        0.047, 0.007, 0.0, 0.013, 0.015, 0.025];
    EDLCrossPoints = [3, 3, 5];
    EDLMIF = 512;  %max isometric force
end
Muscle3 = PamData('Extensor Digitorum Longus', Location3, EDLCrossPoints, EDLMIF, T, Axis1);


%EHL t2 to a6
if beginOptimization == 0
    Location4 = [0.006, 0.033, 0.097, 0.129, 0.173, 0.03, 0.056;
                -0.062, -0.398, 0.039, 0.031, 0.014, 0.004, 0.003;
                0.047, -0.008, -0.021, -0.026, -0.028, -0.024, -0.019];
    EHLViaPoints = [3, 3, 6];
    EHLMIF = 162;  %max isometric force
end
Muscle4 = PamData('Extensor Hallucis Longus', Location4, EHLViaPoints, EHLMIF, T, Axis1);

%Torque Calcs, "R" for robot
Torque1 = Muscle1.Torque(1, :, :) + Muscle2.Torque(1, :, :) + Muscle3.Torque(1, :, :) + Muscle4.Torque(1, :, :);
Torque2 = Muscle1.Torque(2, :, :) + Muscle2.Torque(2, :, :) + Muscle3.Torque(2, :, :) + Muscle4.Torque(2, :, :);
Torque3 = Muscle1.Torque(3, :, :) + Muscle2.Torque(3, :, :) + Muscle3.Torque(3, :, :) + Muscle4.Torque(3, :, :);

%Create the Mesh of Torques to corespond with the joint angles
Torque1M = zeros(divisions, divisions); Torque2M = zeros(divisions, divisions); Torque3M = zeros(divisions, divisions);
for i = 1:divisions
    Torque1M(:, i) = Torque1(((i-1)*divisions)+1:i*divisions);
    Torque2M(:, i) = Torque2(((i-1)*divisions)+1:i*divisions);
    Torque3M(:, i) = Torque3(((i-1)*divisions)+1:i*divisions);
end

%Including Generic variable name for plotting
RobotAxis1 = Joint1.Theta;
RobotAxis1Label = 'Ankle Flexion, Degrees';
RobotAxis2 = Joint3.Theta;
RobotAxis2Label = 'MTP Flexion, Degrees';
RobotTorque1 = Torque1M;
RobotTorque2 = Torque2M;
RobotTorque3 = Torque3M;
RobotTitle1 = 'Robot Ankle Torque, Z Axis';
RobotTitle2 = 'Robot Subtalar Torque, X Axis';
RobotTitle3 = 'Robot MTP Torque, Z Axis';


%Store all necessary variables into a cell array for the optimization
%script
Muscles = {Muscle1, Muscle2, Muscle3, Muscle4};