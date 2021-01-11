% Bone Mesh Code
% Author: Connor Morrow
% Date: 1/14/2020
% Description: This script pulls all of the data captured from OpenSim bone 
% geometries to create plots of the bone points. These points will be used
% in the muscle placement optimization code.

% clc
% clear 
% close all

addpath('Open_Sim_Bone_Geometry')

%Spine
PointsFile = 'Spine_Mesh_Points.xlsx';

Spine = xlsread(PointsFile);

figure
hold on
plot3(Spine(:, 1), Spine(:, 2), Spine(:, 3), '.')
plot3(0, 0, 0, 'o', 'color', 'r')
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
hold off

%Sacrum
PointsFile = 'Sacrum_Mesh_Points.xlsx';

Sacrum = xlsread(PointsFile);

figure
hold on
plot3(Sacrum(:, 1), Sacrum(:, 2), Sacrum(:, 3), '.')
plot3(0, 0, 0, 'o', 'color', 'r')
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
hold off

%Pelvis
PointsFile = 'Pelvis_R_Mesh_Points.xlsx';

Pelvis = xlsread(PointsFile);

figure
hold on
plot3(Pelvis(:, 1), Pelvis(:, 2), Pelvis(:, 3), '.', 'color', 'b')
plot3(0, 0, 0, 'o', 'color', 'r')
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
hold off

%Tibia
PointsFile = 'Tibia_Mesh_Points.xlsx';

Tibia = xlsread(PointsFile);

figure
hold on
plot3(Tibia(:, 1), Tibia(:, 2), Tibia(:, 3), '.')
plot3(0, 0, 0, 'o','color',  'r')
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
hold off

%Femur
PointsFile = 'Femur_Mesh_Points.xlsx';

Femur = xlsread(PointsFile);

figure
hold on
plot3(Femur(:, 1), Femur(:, 2), Femur(:, 3), '.')
plot3(0, 0, 0, 'o', 'color', 'r')
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
hold off

%Talus
PointsFile = 'Talus_Mesh_Points.xlsx';

Talus = xlsread(PointsFile);

figure
hold on
plot3(Talus(:, 1), Talus(:, 2), Talus(:, 3), '.')
plot3(0, 0, 0, 'o', 'color', 'r')
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
hold off

%Calcaneus
PointsFile = 'Calcaneus_Mesh_Points.xlsx';

Calcaneus = xlsread(PointsFile);

figure
hold on
plot3(Calcaneus(:, 1), Calcaneus(:, 2), Calcaneus(:, 3), '.')
plot3(0, 0, 0, 'o', 'color', 'r')
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
hold off

%Toes
PointsFile = 'Toes_Mesh_Points.xlsx';

Toes = xlsread(PointsFile);

figure
hold on
plot3(Toes(:, 1), Toes(:, 2), Toes(:, 3), '.')
plot3(0, 0, 0, 'o', 'color', 'r')
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
hold off

%Bifemoral Head, test to see muscle path
Bifemlh = [-0.126, -0.03, -0.023;
            -0.103, -0.036, -0.056;
            0.069, 0.029, 0.034];
        
Bifemlh = Bifemlh';

%Perform translations on the points, so that they can all be plotted
%together to form the leg
axisLimits = [-1 1 -1 1 -1.25 0.75];

Back = [-0.1007, 0.0815 0];
SpineM = Spine+Back;

Hip = [-0.0707, -0.0661, 0.0835];
FemurM = Femur+Hip;

Knee = [-0.0045 -0.3958 0];
TibiaM = Tibia+Hip+Knee; %error

Ankle = [0 -0.43, 0];   %From the hip?
TalusM = Talus+Hip+Knee+Ankle;

Subtalar = [-0.04877 -0.04195 0.00792];
CalcaneusM = Calcaneus+Hip+Knee+Ankle+Subtalar;

MTP = [0.1788 -0.002 0.00108];
ToesM = Toes+Hip+Knee+Ankle+Subtalar+MTP;

BifemlhM = Bifemlh;
BifemlhM(2:3, :) = Bifemlh(2:3, :)+Hip+Knee;

Connection1 = [Back; Hip];
Connection2 = [Hip; Knee+Hip];
Connection3 = [Knee+Hip; Ankle+Knee+Hip];

figure
hold on
plot3(0, 0, 0, 'o', 'color', 'r')
plot3(SpineM(:, 1), -SpineM(:, 3), SpineM(:, 2), '.', 'color', 'b')
plot3(Sacrum(:, 1), -Sacrum(:, 3), Sacrum(:, 2), '.', 'color', 'b')
plot3(Pelvis(:, 1), -Pelvis(:, 3), Pelvis(:, 2), '.', 'color', 'b')
plot3(FemurM(:, 1), -FemurM(:, 3), FemurM(:, 2), '.', 'color', 'b')
plot3(TibiaM(:, 1), -TibiaM(:, 3), TibiaM(:, 2), '.', 'color', 'b')
plot3(TalusM(:, 1), -TalusM(:, 3), TalusM(:, 2), '.', 'color', 'b')
plot3(CalcaneusM(:, 1), -CalcaneusM(:, 3), CalcaneusM(:, 2), '.', 'color', 'b')
plot3(ToesM(:, 1), -ToesM(:, 3), ToesM(:, 2), '.', 'color', 'b')
plot3(BifemlhM(:, 1), -BifemlhM(:, 3), BifemlhM(:, 2), 'color', 'r', 'LineWidth', 3)
plot3(Connection1(:, 1), -Connection1(:, 3), Connection1(:, 2), 'color', 'g', 'LineWidth', 2)
plot3(Connection2(:, 1), -Connection2(:, 3), Connection2(:, 2), 'color', 'g', 'LineWidth', 2)
plot3(Connection3(:, 1), -Connection3(:, 3), Connection3(:, 2), 'color', 'g', 'LineWidth', 2)
axis(axisLimits)
hold off
