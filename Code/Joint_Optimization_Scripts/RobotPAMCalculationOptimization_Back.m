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
    ChooseJoint = 'Back';
end

%Flag that tells the code to not run pieces that have been calculated
%previously
if exist('beginOptimization', 'var') == 0
    beginOptimization = 0;
end

if beginOptimization == 0
    %Back x axis rotation
    Raxis = [1; 0; 0];
    MaxTheta = 15*pi/180;
    MinTheta = -15*pi/180;
    Home = [-0.1007; 0.0815; 0];                   
    Joint1a = JointData('BackX', Raxis, MaxTheta, MinTheta, Home, divisions);

    %Back y rotation
    Raxis = [0; 1; 0];
    MaxTheta = 0;
    MinTheta =0;               %Ben's script lists the min and max, however only 0 is used. Why?
    Home = [-0.1007; 0.0815; 0];
    Joint1b = JointData('BackY', Raxis, MaxTheta, MinTheta, Home, divisions);

    %Back z rotation
    Raxis = [0; 0; 1];
    MaxTheta = 15*pi/180;
    MinTheta = -45*pi/180;
    Home = [-0.1007; 0.0815; 0];
    Joint1c = JointData('BackZ', Raxis, MaxTheta, MinTheta, Home, divisions);

    R = zeros(3, 3, 1, divisions^2);
    T = zeros(4, 4, 1, divisions^2);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            R(:, :, 1, ii) = Joint1a.RotMat(:, :, i)*Joint1b.RotMat(:, :, i)*Joint1c.RotMat(:, :, j);
            T(:, :, 1, ii) = [R(:, :, 1, ii), Home; 0 0 0 1];
            ii = ii+1;
        end
    end

    T(:, :, 2, :) = T(:, :, 1, :);
    T(:, :, 3, :) = T(:, :, 1, :);
end


if beginOptimization == 0
    %For Erector Spinae
    Location1 = [-0.122, -0.158, -0.072;
                -0.051, 0.054, 0.172;
                0.078, 0.041, 0.031];
    ESCrossPoints = 3;                    %Via points are the points where a transformation matrix is needed. Typically wrap point + 1
    ESMIF = 2500;
    Axis1 = [10, 20, 30];                           %The axis of interest when calculating the moment arm about each joint. The axis is 1, but is listed as 10 so that the cross product doesn't rotate the resulting vector. See PamData > CrossProd

    %For Internal Oblique  p2 -> b3
    Location2 = [-0.04, 0.07;
                 0.07, 0.16;
                 0.116, 0.015];
    IOCrossPoints = 2;                    %Via points are the points where a transformation matrix is needed. Typically wrap point + 1
    IOMIF = 900;
    AxisIO = [10, 20, 30];            

    %For External Oblique p4 -> b4
    Location3 = [-0.03, 0.065;
                 -0.064, 0.11;
                 0.01, 0.11];
    EOCrossPoints = 2;                    %Via points are the points where a transformation matrix is needed. Typically wrap point + 1
    EOMIF = 900;
    AxisEO = [10, 20, 30];        

    Location = {Location1, Location2, Location3};
end

Muscle1 = PamData('Erector Spinae', Location{1}, ESCrossPoints, ESMIF, T, Axis1);
Muscle2 = PamData('Internal Oblique', Location{2}, IOCrossPoints, IOMIF, T, AxisIO);
Muscle3 = PamData('External Oblique', Location{3}, EOCrossPoints, EOMIF, T, AxisEO);

Torque1 = Muscle1.Torque+Muscle2.Torque+Muscle3.Torque;

%Reorganize the torque calculations into a matrix for 3D plotting
for i = 1:divisions
    Torque1M(:, i, :) = Torque1(:, ((i-1)*divisions)+1:i*divisions, :);
end

%Transposing to match Ben's results
for i = 1:size(Torque1M, 3)
    Torque1M(:, :, i) = Torque1M(:, :, i)';
end

%Including a way to make the data generic, so that things can be
%plotted in one location instead of spread out through all of the
%different sections
RobotAxis1 = Joint1c.Theta;
RobotAxis1Label = 'Lumbar Bending, Degrees';
RobotAxis2 = Joint1a.Theta;
RobotAxis2Label = 'Back Flexion, Degrees';
RobotTorque1 = Torque1M(:, :, 1);
RobotTorque2 = Torque1M(:, :, 2);
RobotTorque3 = Torque1M(:, :, 3);
RobotTitle1 = 'Robot Back Torque, Erector Spinae, X Axis';
RobotTitle2 = 'Robot Back Torque, Erector Spinae, Y Axis';
RobotTitle3 = 'Robot Back Torque, Erector Spinae, Z Axis';

%Store all necessary variables into a cell array for the optimization
%script
Muscles = {Muscle1, Muscle2, Muscle3};