% BoneMeshSkeleton
% Author: Connor Morrow
% Date: 11/11/2020
% Description: This Script plots on the skeletal geometry of the human

clc
clear all
close all

current_dir = cd;
all_code = fullfile(current_dir,'../..');
addpath(genpath(all_code));
% addpath('Open_Sim_Bone_Geometry')
% addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Functions
% addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code
% addpath C:\users\Connor\Documents\GitHub\Bipedal_Robot\Code\Human_Data

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
% Femur = Femur;

Knee = [-0.0045, -0.3958, 0];
% Tibia = Tibia; %error

Ankle = [0, -0.43, 0];   %From the hip?
% Talus = Talus+Hip+Knee+Ankle;

Subtalar = [-0.04877, -0.04195, 0.00792];
% Calcaneus = Calcaneus+Hip+Knee+Ankle+Subtalar;

MTP = [0.1788, -0.002, 0.00108];
% Toes = Toes+Hip+Knee+Ankle+Subtalar+MTP;

%Rotate the points to be in an upward orientation when plotted
RotationM = [1, 0, 0;
            0, 0, 1;
            0, -1, 0];
SpineP = Spine*RotationM;
SacrumP = Sacrum*RotationM;
PelvisP = Pelvis*RotationM;
FemurP = Femur*RotationM;
TibiaP = Tibia*RotationM;
TalusP = Talus*RotationM;
CalcaneusP = Calcaneus*RotationM;
ToesP = Toes*RotationM;

%% Plot everything
axisLimits = [-1 1 -1 1 -1.25 0.75];

figure
hold on
plot3(SpineP(:, 1), SpineP(:, 2), SpineP(:, 3), '.', 'color', 'b');
plot3(SacrumP(:, 1), SacrumP(:, 2), SacrumP(:, 3), '.', 'color', 'b');
plot3(PelvisP(:, 1), PelvisP(:, 2), PelvisP(:, 3), '.', 'color', 'b');
plot3(FemurP(:, 1), FemurP(:, 2), FemurP(:, 3), '.', 'color', 'b');
plot3(TibiaP(:, 1), TibiaP(:, 2), TibiaP(:, 3), '.', 'color', 'b');
plot3(TalusP(:, 1), TalusP(:, 2), TalusP(:, 3), '.', 'color', 'b');
plot3(CalcaneusP(:, 1), CalcaneusP(:, 2), CalcaneusP(:, 3), '.', 'color', 'b');
plot3(ToesP(:, 1), ToesP(:, 2), ToesP(:, 3), '.', 'color', 'b');
xlabel('X')
ylabel('Z')
zlabel('Y')
axis(axisLimits)
view([0, -1, 0])
hold off

%% Plot Continuously
iteration = 100;

flexMax = 80*pi/180;
extMax = -20*pi/180;
phi = linspace(extMax, flexMax, iteration);
pelvisToHip = [-0.0707, -0.0661, 0.0835];
for i = 1:iteration
    R(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    
    T(:, :, i) = RpToTrans(R(:, :, i), pelvisToHip');
    T2(:, :, i) = RpToTrans(eye(3), Knee');
end


Name = 'Gluteus Maximus 1';
MIF = 573;
OFL = 0.142; TSL = 0.125; Pennation = 0.08726646;
Location = [-0.119, 0.061, 0.07;
            -0.129, 0.001, 0.089;
            -0.046, -0.025, 0.039;
            -0.028, -0.057, 0.047];
CrossPoint = 3;
Glut_Max1 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T); 

for i = 1:iteration
    Glut_Max1_Location(:, :, i) = Glut_Max1.Location;
end
for ii = 1:iteration
    for i = 1:size(Glut_Max1_Location, 1)
        if i >= CrossPoint
            Glut_Max1_Location(i, :, ii) = RowVecTrans(T(:, :, ii), Glut_Max1_Location(i, :, ii))*RotationM;
        else
            Glut_Max1_Location(i, :, ii) = Glut_Max1_Location(i, :, ii)*RotationM;
        end
    end
end

Name = 'Rectus Femoris (Quadriceps)';   %The Location should be updated to reflect the moving patella points. For now it is static in the knee frame
mif = 1169;
Pennation = 0.08726646;
Location = [-0.029, -0.031, 0.097;
            0.033, -0.403, 0.002;
            0.062, 0.021, 0.0014];
CrossPoint = [2, 3];
for i = 1:iteration
    Rect_Location(:, :, i) = Location;
end
for ii = 1:iteration
    for i = 1:size(Rect_Location, 1)
        if i >= CrossPoint(2)
            Rect_Location(i, :, ii) = RowVecTrans(T(:, :, ii)*T2(:, :, ii), Rect_Location(i, :, ii))*RotationM;
        elseif i>= CrossPoint(1)
            Rect_Location(i, :, ii) = RowVecTrans(T(:, :, ii), Rect_Location(i, :, ii))*RotationM;
        else
            Rect_Location(i, :, ii) = Rect_Location(i, :, ii)*RotationM;
        end
    end
end
    


h = figure
filename = 'testAnimated.gif';
for i = 1:iteration
    cla
    
    %Calculate new skeletal points
    R = [cos(phi(i)), -sin(phi(i)), 0; sin(phi(i)), cos(phi(i)), 0; 0, 0, 1];
    FemurT = RpToTrans(R, Hip');
    TibiaT = RpToTrans(eye(3), Knee');
    TalusT = RpToTrans(eye(3), Ankle');
    CalcaneusT = RpToTrans(eye(3), Subtalar');
    ToesT = RpToTrans(eye(3), MTP');
    
    for ii = 1:size(Femur, 1)
        FemurP(ii, :) = RowVecTrans(FemurT, Femur(ii, :))*RotationM;
    end
    for ii = 1:size(Tibia, 1)    
        TibiaP(ii, :) = RowVecTrans(FemurT*TibiaT, Tibia(ii, :))*RotationM;
    end
    for ii = 1:size(Talus, 1)
        TalusP(ii, :) = RowVecTrans(FemurT*TibiaT*TalusT, Talus(ii, :))*RotationM;
    end
    for ii = 1:size(Calcaneus, 1)
        CalcaneusP(ii, :) = RowVecTrans(FemurT*TibiaT*TalusT*CalcaneusT, Calcaneus(ii, :))*RotationM;
    end
    for ii = 1:size(Toes, 1)
        ToesP(ii, :) = RowVecTrans(FemurT*TibiaT*TalusT*CalcaneusT*ToesT, Toes(ii, :))*RotationM;
    end

    hold on
    plot3(SpineP(:, 1), SpineP(:, 2), SpineP(:, 3), '.', 'color', 'b');
    plot3(SacrumP(:, 1), SacrumP(:, 2), SacrumP(:, 3), '.', 'color', 'b');
    plot3(PelvisP(:, 1), PelvisP(:, 2), PelvisP(:, 3), '.', 'color', 'b');
    plot3(FemurP(:, 1), FemurP(:, 2), FemurP(:, 3), '.', 'color', 'b');
    plot3(TibiaP(:, 1), TibiaP(:, 2), TibiaP(:, 3), '.', 'color', 'b');
    plot3(TalusP(:, 1), TalusP(:, 2), TalusP(:, 3), '.', 'color', 'b');
    plot3(CalcaneusP(:, 1), CalcaneusP(:, 2), CalcaneusP(:, 3), '.', 'color', 'b');
    plot3(ToesP(:, 1), ToesP(:, 2), ToesP(:, 3), '.', 'color', 'b');
    
    plot3(Glut_Max1_Location(:, 1, i), Glut_Max1_Location(:, 2, i), Glut_Max1_Location(:, 3, i), 'color', 'r', 'LineWidth', 4) 
    plot3(Rect_Location(:, 1, i), Rect_Location(:, 2, i), Rect_Location(:, 3, i), 'color', 'r', 'LineWidth', 4)
    
    xlabel('X')
    ylabel('Z')
    zlabel('Y')
    axis(axisLimits)
    view([0, -1, 0])
    pause(0.1);
    hold off
    
    %Create gif
    frame = getframe(h);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf)
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
    end
end



