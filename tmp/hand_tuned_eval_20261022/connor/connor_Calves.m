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
    ChooseJoint = 'Calves';
end

%Flag that tells the code to not run pieces that have been calculated
%previously
if exist('beginOptimization', 'var') == 0
    beginOptimization = 0;
end

%Create a Joint object that calculates things like transformation
%matrix
%Note: Currently only goes for the matrices at minimum theta. Will
%implement some iteration process to go from minimum to maximum. 
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

    %Knee z rotation
    Raxis = [0; 0; 1];
    MaxTheta = 0.17453293;
    MinTheta = -2.0943951;
    Home = [-0.0707; -0.0661; 0.0835];      %Not actually the home position, but a placeholder
    Joint3 = JointData('Knee', Raxis, MaxTheta, MinTheta, Home, divisions);

    knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
    knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
    KneeXFunc = fit(knee_angle_x,knee_x,'cubicspline');
    knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
    knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
    KneeYFunc = fit(knee_angle_y,knee_y,'cubicspline');

    kneeHome = zeros(3, divisions);
    for j = 1:divisions
        kneeHome(:, j) = [KneeXFunc(Joint3.Theta(j)); KneeYFunc(Joint3.Theta(j)); 0];
    end

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
            T(:, :, 2, ii) = T_s(:, :, i);            %Subtalar doesn't change
            ii = ii + 1;
        end
    end
end

%Soleus a1 -> t1
if beginOptimization == 0
    Location3 = [0.047, 0;
                -0.113, 0.031;
                -0.002, -0.005];
    SCrossPoints = [2 2];
    SMIF = 3549;  %max isometric force
    Axis3 = [3; 1];                           %The axis of interest when calculating the moment arm about each joint        
end


Muscle3 = PamData('Soleus', Location3, SCrossPoints, SMIF, T, Axis3);

T(:, :, 2:3, :) = T(:, :, 1:2, :);
ii = 1;
for i = 1:divisions
    for j = 1:divisions
        T(:, :, 1, ii) = [Joint3.RotMat(:, :, j), kneeHome(:, j); 0, 0, 0, 1]; %Change the knee position
        ii = ii + 1;
    end
end

%Medial Gastroocnemius a1 -> f4      
if beginOptimization == 0
    Location1 = [  0.015, -0.03, 0;
                    -0.374, -0.402, 0.031;
                    -0.025, -0.026, -0.005];
    MGCrossPoints = [3 3 3];
    MGMIF = 1558;  %max isometric force
    Axis1 = [3; 3; 1];                           %The axis of interest when calculating the moment arm about each joint
end  

Muscle1 = PamData('Medial Gastroocnemius', Location1, MGCrossPoints, MGMIF, T, Axis1);

%Lateral Gastrocenemius a1 -> f5
if beginOptimization == 0
    Location2 = [0.009, -0.022, 0;
                    -0.378, -0.395, 0.031;
                    0.027, 0.027, -0.005];

    LGCrossPoints = [3 3 3];
    LGMIF = 683;  %max isometric force
    Axis2 = [3; 3; 1];                           %The axis of interest when calculating the moment arm about each joint
end
Muscle2 = PamData('Lateral Gastrocenemius', Location2, LGCrossPoints, LGMIF, T, Axis2);

%Torque Calcs, "R" for robot
Torque1 = Muscle1.Torque(2, :) + Muscle2.Torque(2, :);   %Torque about the ankle
Torque2 = Muscle1.Torque(1, :) + Muscle2.Torque(1, :);        %Torque about the knee
Torque3 = Muscle3.Torque(1, :);                       %Torque about Ankle due to Soleus

%Create the Mesh of Torques to corespond with the joint angles
Torque1M = zeros(divisions, divisions); Torque3M = zeros(divisions, divisions); Torque2M = zeros(divisions, divisions);
for i = 1:divisions
    Torque1M(:, i) = Torque1(((i-1)*divisions)+1:i*divisions);
    Torque2M(:, i) = Torque2(((i-1)*divisions)+1:i*divisions);
    Torque3M(:, i) = Torque3(((i-1)*divisions)+1:i*divisions);

end

%Including Generic variable name for plotting
RobotAxis1 = Joint1.Theta;
RobotAxis1Label = 'Ankle Flexion, Degrees';
RobotAxis2 = Joint3.Theta;
RobotAxis2Label = 'Knee Flexion, Degrees';
RobotTorque1 = Torque2M;
RobotTorque2 = Torque1M;
RobotTorque3 = Torque3M;
RobotTitle1 = 'Robot Knee Torque, Gastrocnemius, Z axis';
RobotTitle2 = 'Robot Ankle Torque, Gastrocnemius, Z` axis';
RobotTitle3 = 'Robot Ankle Torque, Soleus, Z` axis';


%Store all necessary variables into a cell array for the optimization
%script
Muscles = {Muscle1, Muscle2, Muscle3};