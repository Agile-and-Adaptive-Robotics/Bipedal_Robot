% Bone_Mesh_Plotting
% Author: Connor Morrow
% Date: 10/14/2020
% Description: This script pulls date from bone scans used for OpenSim to
% generate a pot matrix representation of the bone, and then plots the
% muscle pathing on top of it.

addpath('Open_Sim_Bone_Geometry')

%Retrieve the datasets
Spine = xlsread('Spine_Mesh_Points.xlsx');
Sacrum = xlsread('Sacrum_Mesh_Points.xlsx');
Pelvis = xlsread('Pelvis_R_Mesh_Points.xlsx');
Tibia = xlsread('Tibia_Mesh_Points.xlsx');
Femur = xlsread('Femur_Mesh_Points.xlsx');
Talus = xlsread('Talus_Mesh_Points.xlsx');
Calcaneus = xlsread('Calcaneus_Mesh_Points.xlsx');
Toes = xlsread('Toes_Mesh_Points.xlsx');

%Add joint offsets to each of the datasets (Sacrum and Pelvis are already
%in the correct frame)
Back = [-0.1007, 0.0815, 0];
Spine = Spine + Back;

Hip = [-0.0707, -0.0661, 0.0835];
Femur = Femur+Hip;

Knee = [-0.0045, -0.3958, 0];
Tibia = Tibia+Hip+Knee; %error

Ankle = [0, -0.43, 0];   %From the hip?
Talus = Talus+Hip+Knee+Ankle;

Subtalar = [-0.04877, -0.04195, 0.00792];
Calcaneus = Calcaneus+Hip+Knee+Ankle+Subtalar;

MTP = [0.1788, -0.002, 0.00108];
Toes = Toes+Hip+Knee+Ankle+Subtalar+MTP;

%Rotate the points to be in an upward orientation when plotted
RotationM = [1, 0, 0;
            0, 0, 1;
            0, -1, 0];
Spine = Spine*RotationM;
Sacrum = Sacrum*RotationM;
Pelvis = Pelvis*RotationM;
Femur = Femur*RotationM;
Tibia = Tibia*RotationM;
Talus = Talus*RotationM;
Calcaneus = Calcaneus*RotationM;
Toes = Toes*RotationM;

%Take the data from the optimization and put it into the correct frame
%---------- In Progress ------------
Bifemlh = [-0.126, -0.03, -0.023;
            -0.103, -0.036, -0.056;
            0.069, 0.029, 0.034];

Bifemlh = Bifemlh';
Bifemlh(2:3, :) = Bifemlh(2:3, :) + Hip + Knee; %2 is probably a cross point

Bifemlh = Bifemlh*RotationM;

%Plot everything
axisLimits = [-1 1 -1 1 -1.25 0.75];

figure
hold on
plot3(Spine(:, 1), Spine(:, 2), Spine(:, 3), '.', 'color', 'b');
plot3(Sacrum(:, 1), Sacrum(:, 2), Sacrum(:, 3), '.', 'color', 'b');
plot3(Pelvis(:, 1), Pelvis(:, 2), Pelvis(:, 3), '.', 'color', 'b');
plot3(Femur(:, 1), Femur(:, 2), Femur(:, 3), '.', 'color', 'b');
plot3(Tibia(:, 1), Tibia(:, 2), Tibia(:, 3), '.', 'color', 'b');
plot3(Talus(:, 1), Talus(:, 2), Talus(:, 3), '.', 'color', 'b');
plot3(Calcaneus(:, 1), Calcaneus(:, 2), Calcaneus(:, 3), '.', 'color', 'b');
plot3(Toes(:, 1), Toes(:, 2), Toes(:, 3), '.', 'color', 'b');
% plot3(BifemlhM(:, 1), -BifemlhM(:, 3), BifemlhM(:, 2), 'color', 'r', 'LineWidth', 3)

plot3(Bifemlh(:, 1), Bifemlh(:, 2), Bifemlh(:, 3), 'color', 'r', 'LineWidth', 3)
axis(axisLimits)
hold off
