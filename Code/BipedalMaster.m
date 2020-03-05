% Bipedal Master
% Author: Connor Morrow
% Date: 1/14/2020
% Description: This script runs calculations to generate the human torque
% plots adn robot torque plots.

if exist('Optimize', 'var') == 0
    clear
    clc
    close all
end

%This script will run both the humanoid model and the robot model scripts
%for generating Torque, allowing us to easily view both and make
%comparisons.

%Choose the Joint that you will to observe
%Options are: Back, Bi_Hip, Calves, Foot, Toe, Uni_Hip
ChooseJoint = 'Back';

% %Choose the number of divisions for the angles of rotation
divisions = 100;

%% ----------------- Setup -------------------------------
%Include relevant folders
addpath('Open_Sim_Bone_Geometry')
addpath('Functions')
addpath('Human_Data')

%% ------------- Humanoid Model --------------
%Runs the humanoid model. Only run if you need to update the data
%run("HumanoidMuscleCalculation.m")

%Loads important data for the human model. Data was previously created and
%then stored, for quicker load time.
%Data saved: Back, Bi_Hip, Calves, Foot, Toe, Uni_Hip
load(strcat('Human_', ChooseJoint, '_Data.mat'));

%% ------------- Robot Model -----------------
%Runs the bipedal model
run("RobotPAMCalculationOptimization.m")    
% run("RobotPAMCalculation.m")    

%% ------------- Error Caclulation ------------
%After running both scripts, this portion checks to see the error between
%the two calculations

Error1 = HumanTorque1 - RobotTorque1;
Error2 = HumanTorque2 - RobotTorque2;

if exist('HumanTorque3', 'var') == 1
    Error3 = HumanTorque3 - RobotTorque3;
    if exist('HumanTorque4', 'var') == 1
        Error4 = HumanTorque4 - RobotTorque4;
    end
end


%% ------------- Plot Creation -----------

%Decide on the colorbar range for all plots
caxisRange = [-40 150];

%% ------------------ Human Plots --------------------
figure
hold on
surf(HumanAxis1*180/pi, HumanAxis2*180/pi, HumanTorque1, 'EdgeColor', 'none')
title(HumanTitle1); xlabel(HumanAxis1Label); ylabel(HumanAxis2Label); zlabel('Torque, N*m')
colorbar; caxis(caxisRange)
hold off

figure
surf(HumanAxis1*180/pi, HumanAxis2*180/pi, HumanTorque2, 'EdgeColor', 'none')
title(HumanTitle2); xlabel(HumanAxis1Label); ylabel(HumanAxis2Label); zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

if exist('HumanTorque3', 'var') == 1
    figure
    surf(HumanAxis1*180/pi, HumanAxis2*180/pi, HumanTorque3, 'EdgeColor', 'none')
    title(HumanTitle3); xlabel(HumanAxis1Label); ylabel(HumanAxis2Label); zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end

if exist('HumanTorque4', 'var') == 1
    figure
    surf(HumanAxis1*180/pi, HumanAxis2*180/pi, HumanTorque4, 'EdgeColor', 'none')
    title(HumanTitle4); xlabel(HumanAxis1Label); ylabel(HumanAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end

%% ---------------- Robot Plots -------------------

figure
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, RobotTorque1, 'EdgeColor', 'none')
title(RobotTitle1); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

figure
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, RobotTorque2, 'EdgeColor', 'none')
title(RobotTitle2); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

if exist('RobotTorque3', 'var') == 1
    figure
    surf(RobotAxis1*180/pi, RobotAxis2*180/pi, RobotTorque3, 'EdgeColor', 'none')
    title(RobotTitle3); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end

if exist('RobotTorque4', 'var') == 1
    figure
    surf(RobotAxis1*180/pi, RobotAxis2*180/pi, RobotTorque4, 'EdgeColor', 'none')
    title(RobotTitle4); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end

%% ---------------- Error Plots ----------------------
figure
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error1, 'EdgeColor', 'none')
title(strcat(RobotTitle1, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

figure
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error2, 'EdgeColor', 'none')
title(strcat(RobotTitle2, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

if exist('HumanTorque3', 'var') == 1
    Error3 = HumanTorque3 - RobotTorque3;
    figure
    surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error3, 'EdgeColor', 'none')
    title(strcat(RobotTitle3, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end

if exist('HumanTorque4', 'var') == 1
    Error4 = HumanTorque4 - RobotTorque4;
    figure
    surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error4, 'EdgeColor', 'none')
    title(strcat(RobotTitle4, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end
    