% MuscleBonePlotting
% Author: Connor Morrow
% Date: 11/11/2020
% Description: This script plots the bone skeleton and the muscles that are
% currently being investigated

% clc
% clear all
% close all

addpath('Open_Sim_Bone_Geometry')

%% Retrieve the datasets
Spine = xlsread('Spine_Mesh_Points.xlsx');
Sacrum = xlsread('Sacrum_Mesh_Points.xlsx');
Pelvis = xlsread('Pelvis_R_Mesh_Points.xlsx');
Tibia = xlsread('Tibia_Mesh_Points.xlsx');
Femur = xlsread('Femur_Mesh_Points.xlsx');
Talus = xlsread('Talus_Mesh_Points.xlsx');
Calcaneus = xlsread('Calcaneus_Mesh_Points.xlsx');
Toes = xlsread('Toes_Mesh_Points.xlsx');

%% Add joint offsets to each of the datasets (Sacrum and Pelvis are already
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

%% Calculate Muscle Locations
if isequal(Bones{1}, 'Pelvis')
    offset1 = 0;
end

if isequal(Bones{2}, 'Femur')
    offset2 = Hip;
end


for i = 1:size(HMuscleLocation, 2)
    HMuscleLocation{i}(1:HMuscleCross{i}-1, :) = HMuscleLocation{i}(1:HMuscleCross{i}-1, :) + offset1;
    HMuscleLocation{i}(HMuscleCross{i}:end, :) = HMuscleLocation{i}(HMuscleCross{i}:end, :) + offset1 + offset2;
    HMuscleLocation{i} = HMuscleLocation{i}*RotationM;
end

for i = 1:size(RMuscleLocation, 2)
    RMuscleLocation{i}(1:RMuscleCross{i}-1, :) = RMuscleLocation{i}(1:RMuscleCross{i}-1, :) + offset1;
    RMuscleLocation{i}(RMuscleCross{i}:end, :) = RMuscleLocation{i}(RMuscleCross{i}:end, :) + offset1 + offset2;
    RMuscleLocation{i} = RMuscleLocation{i}*RotationM;
end

% for i = 1:size(HMuscleLocation, 2)
%     HMuscleLocation{i}(1:HMuscleCross{i}-1, :) = HMuscleLocation{i}(1:HMuscleCross{i}-1, :) + Hip;
%     HMuscleLocation{i}(HMuscleCross{i}:end, :) = HMuscleLocation{i}(HMuscleCross{i}:end, :) + Hip + Knee;
%     HMuscleLocation{i} = HMuscleLocation{i}*RotationM;
% end
% 
% for i = 1:size(RMuscleLocation, 2)
%     RMuscleLocation{i}(1:RMuscleCross{i}-1, :) = RMuscleLocation{i}(1:RMuscleCross{i}-1, :) + Hip;
%     RMuscleLocation{i}(RMuscleCross{i}:end, :) = RMuscleLocation{i}(RMuscleCross{i}:end, :) + Hip + Knee;
%     RMuscleLocation{i} = RMuscleLocation{i}*RotationM;
% end

% M1Locations = Add_Mag1.Location;
% M2Locations = Add_Mag2.Location;
% M3Locations = Add_Mag3.Location;
% 
% P1Locations = Add_Mag_Pam.Location;
% 
% M1Locations(2, :) = M1Locations(2, :) + Hip;
% M2Locations(2, :) = M2Locations(2, :) + Hip;
% M3Locations(2, :) = M3Locations(2, :) + Hip;
% P1Locations(2, :) = P1Locations(2, :) + Hip;
% 
% M1Locations = M1Locations*RotationM;
% M2Locations = M2Locations*RotationM;
% M3Locations = M3Locations*RotationM;
% 
% P1Locations = P1Locations*RotationM;

%% Plot everything
axisLimits = [-1 1 -1 1 -1.25 0.75];

figure
hold on
%Bone Plotting
plot3(Spine(:, 1), Spine(:, 2), Spine(:, 3), '.', 'color', 'b');
plot3(Sacrum(:, 1), Sacrum(:, 2), Sacrum(:, 3), '.', 'color', 'b');
plot3(Pelvis(:, 1), Pelvis(:, 2), Pelvis(:, 3), '.', 'color', 'b');
plot3(Femur(:, 1), Femur(:, 2), Femur(:, 3), '.', 'color', 'b');
plot3(Tibia(:, 1), Tibia(:, 2), Tibia(:, 3), '.', 'color', 'b');
plot3(Talus(:, 1), Talus(:, 2), Talus(:, 3), '.', 'color', 'b');
plot3(Calcaneus(:, 1), Calcaneus(:, 2), Calcaneus(:, 3), '.', 'color', 'b');
plot3(Toes(:, 1), Toes(:, 2), Toes(:, 3), '.', 'color', 'b');

%MusclePlotting
for i = 1:size(HMuscleLocation, 2)
    plot3(HMuscleLocation{i}(:, 1), HMuscleLocation{i}(:, 2), HMuscleLocation{i}(:, 3), '.-', 'color', 'r', 'linewidth', 3);
end
for i = 1:size(RMuscleLocation, 2)
        plot3(RMuscleLocation{i}(:, 1), RMuscleLocation{i}(:, 2), RMuscleLocation{i}(:, 3), '.-', 'color', 'g', 'linewidth', 3);
end
% plot3(M1Locations(:, 1), M1Locations(:, 2), M1Locations(:, 3), '.-', 'color', 'r');
% plot3(M2Locations(:, 1), M2Locations(:, 2), M2Locations(:, 3), '.-', 'color', 'r');
% plot3(M3Locations(:, 1), M3Locations(:, 2), M3Locations(:, 3), '.-', 'color', 'r');
% plot3(P1Locations(:, 1), P1Locations(:, 2), P1Locations(:, 3), '.-', 'color', 'g');
axis(axisLimits)
hold off
