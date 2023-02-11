%% Mesh Points Calculation - Bicep Femoris Short Head
% This code will calculate the torque difference between all of the points
% from one bone mesh to another to determine the best location for muscle
% placement

%% Freshen up the workspace
clc
clear
close all

%% Add paths to the muscle and pam calculators
current_dir = cd;
all_code = fullfile(current_dir,'../..');
addpath(genpath(all_code));

% addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Human_Data
% addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Robot_Data
% addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Functions
% addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Bone_Mesh_Plots\Open_Sim_Bone_Geometry

%% Joint rotation transformation matrices
positions = 100;
fprintf('The algorithm will be calculating Torque at %d different joint positions.\n', positions)

R = zeros(3, 3, positions);
T = zeros(4, 4, positions);
R_Pam = zeros(3, 3, positions);
T_Pam = zeros(4, 4, positions);
t1toICR = zeros(1,3,positions);
T_t1_ICR = zeros(4, 4, positions);
T_ICR_t1 = zeros(4, 4, positions);

c = pi/180; %Convert from degrees to radians

%Knee Extension and Flexion

%Human knee
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');
%Human Quadriceps muscles
%Rectus Femoris 
%(Insertion point as a function of knee angle due to patella)
rect_fem_x = [0.0156367; 0.0179948; 0.0274274; 0.029683; 0.0306; 0.0366; 0.0422; 0.0451; 0.0484; 0.0533; 0.0617; 0.0634; 0.067; 0.0733];
rect_fem_xD = c*[-120.118; -114.871; -90.068; -83.532; -80; -60; -40; -30; -20; -10; 0; 1.6; 5; 10];
fcn3 = fit(rect_fem_xD,rect_fem_x,'smoothingspline');
rect_fem_y =[0.0234; 0.0238; 0.0251; 0.0253; 0.025284; 0.0249; 0.0243; 0.0239; 0.0234; 0.0228; 0.0210; 0.0206; 0.0192; 0.0160];
rect_fem_yD = c*[-120; -114.6; -90; -83.5; -80.01; -60; -40; -30; -20; -10; 0; 1.6; 5; 10 ];
fcn4 = fit(rect_fem_yD,rect_fem_y,'smoothingspline');
%Vastus Medialis
%(Insertion point as a function of knee angle due to patella)
fcn5 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.009811),'smoothingspline');
fcn6 = fit(rect_fem_yD,rect_fem_y-(0.02346-0.02242),'smoothingspline');
%Vastus Intermedius
%(Insertion point as a function of knee angle due to patella)
fcn7 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.008285),'smoothingspline');
fcn8 = fit(rect_fem_yD,rect_fem_y+(0.0256239-0.02346),'smoothingspline');
%Vastus Lateralis
%(Insertion point as a function of knee angle due to patella)
fcn9 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.0142881),'smoothingspline');
fcn10 = fit(rect_fem_yD,rect_fem_y-(0.02346-0.0215281),'smoothingspline');


%Robot Knee
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam =     ([23.30	22.22	21.55	21.09	19.91	18.70	17.48	13.82	10.44	7.60	5.52	4.35	4.16	5.01	7.04	10.47]')/1000;
fcn11 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam =     ([-416.65	-417.03	-417.19	-417.28	-417.41	-417.41	-417.30	-416.28	-414.36	-411.72	-408.62	-405.32	-402.08	-399.16	-396.85	-395.66]')/1000;
fcn12 = fit(knee_angle,knee_y_Pam,'cubicspline');

%Theta1 to ICR
t1_ICR_x = ([29.66	28.54	27.86	27.40	26.23	25.03	23.81	20.03	16.17	12.34	8.67	5.24	2.04	-1.01	-4.1	-7.58]')/1000;
fcn13 = fit(knee_angle,t1_ICR_x,'cubicspline');
t1_ICR_y = ([25.97	25.74	25.61	25.53	25.35	25.19	25.03	24.57	24.04	23.39	22.66	21.93	21.32	20.99	21.2	22.33]')/1000;
fcn14 = fit(knee_angle,t1_ICR_y,'cubicspline');

kneeMin = -2.0943951;
kneeMax = 0.17453293;
phi = linspace(kneeMin, kneeMax, positions);
%We want one of our positions to be home position, so let's make the
%smallest value of phi equal to 0
[val, pos] = min(abs(phi));
phi(pos) = 0;
phiD = phi*180/pi;

for i = 1:positions
    hipToKnee = [fcn1(phi(i)), fcn2(phi(i)), 0];
    R(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    
    T(:, :, i) = RpToTrans(R(:, :, i), hipToKnee');
    
    hipToKnee_Pam = [fcn11(phi(i)), fcn12(phi(i)), 0];
    R_Pam(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;   %Rotation matrix for robot
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    
    T_Pam(:, :, i) = RpToTrans(R_Pam(:, :, i), hipToKnee_Pam');     %Transformation matrix for robot
    
    %Transformation matrix for ICR to theta1 and inverse
    t1toICR(1,:,i) = [fcn13(phi(i)), fcn14(phi(i)), 0];
    T_t1_ICR(:, :, i) = RpToTrans(eye(3), t1toICR(1,:,i)');    
    T_ICR_t1(:, :, i) = RpToTrans(eye(3), -t1toICR(1,:,i)');
end

%% Muscle calculation
Name = 'Vastus Medialis';
MIF = 1294;
OFL = 0.089; TSL = 0.126; Pennation = 0.08726646;
Location = zeros(5,3,positions);
for i = 1:positions
    if phiD(i) < -101
     Location(:,:,i) = [0.014, -0.21, 0.019;
                 0.036, -0.277, 0.001;
                 0.037, -0.405, -0.013;
                 0.027, -0.425, -0.013;
                 fcn5(phi(i)), fcn6(phi(i)) -0.0146];
    elseif phiD(i) > -101 && phiD(i) < -69
     Location(:,:,i) = [0.014, -0.21, 0.019;
                 0.036, -0.277, 0.001;
                 0.037, -0.405, -0.013;
                 0.037, -0.405, -0.013;
                 fcn5(phi(i)), fcn6(phi(i)) -0.0146];
    else 
     Location(:,:,i) = [0.014, -0.21, 0.019;
                 0.036, -0.277, 0.001;
                 0.036, -0.277, 0.001;
                 0.036, -0.277, 0.001;
                 fcn5(phi(i)), fcn6(phi(i)) -0.0146];
    end
end
CrossPoint = 5;
Vas_Med = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Vastus Intermedius';
MIF = 1365;
OFL = 0.087; TSL = 0.136; Pennation = 0.05235988;
Location = zeros(4,3,positions);
for i = 1:positions
    if phiD(i) < -80
    Location(:,:,i) = [0.029, -0.192, 0.031;
                0.034, -0.208, 0.029;
                0.034, -0.403, 0.005;
                fcn7(phi(i)), fcn8(phi(i)) 0.0018];
    else
    Location(:,:,i) = [0.029, -0.192, 0.031;
                0.034, -0.208, 0.029;
                0.034, -0.208, 0.029;
                fcn7(phi(i)), fcn8(phi(i)) 0.0018];
    end
end
CrossPoint = 4;
Vas_Int = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Vastus Lateralis';
MIF = 1871;
OFL = 0.084; TSL = 0.157; Pennation = 0.08726646;
Location = zeros(5,3,positions);
for i = 1:positions
    if phiD(i) < -110
    Location(:,:,i) = [0.005, -0.185, 0.035;
                0.027, -0.259, 0.041;
                0.036, -0.403, 0.021;
                0.025, -0.424, 0.018;
                fcn9(phi(i)), fcn10(phi(i)) 0.0165];
    elseif phiD(i) > -110 && phiD(i) < -69
    Location(:,:,i) = [0.005, -0.185, 0.035;
                0.027, -0.259, 0.041;
                0.036, -0.403, 0.021;
                0.036, -0.403, 0.021;
                fcn9(phi(i)), fcn10(phi(i)) 0.0165];
    else
    Location(:,:,i) = [0.005, -0.185, 0.035;
                0.027, -0.259, 0.041;
                0.027, -0.259, 0.041;
                0.027, -0.259, 0.041;
                fcn9(phi(i)), fcn10(phi(i)) 0.0165];
    end
end
CrossPoint = 5;
Vas_Lat = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

%% PAM calculation
Name = 'Vastus Intermedius';
Location = zeros(8,3,positions);
%Origin and routing location from Ben
p1 = [0.040, 0.035, 0];                %Origin
p2 = [0.0739, -0.27476, 0];            %BPA contacts mounting base
p3 = [0.05817, -0.37427, 0.000];       %femur channel contact
p4 = [0.04317, -0.41646, 0.000];       %femoral condyle contact
p5 = [0.05871, 0.02916, 0];             %Tibia contact initial
p6 = [0.05871, 0.01646, 0];             %Tibia contact
p7 = [0.05230, -0.00793, 0];            %Tibia tendon contact
p8 = [0.03598, -0.02845, 0];            %patellar ligament ring

%Points for transformation matrices;
v5 = zeros(1,3,positions);
v6 = zeros(1,3,positions);
v7 = zeros(1,3,positions);
v8 = zeros(1,3,positions);

%Set up angle limits (degrees)
AA = -68;      
BB = -32.8;
CC = -8.42;
DD = 6.78;
for i=1:positions
    
    v5(:, :, i) = RowVecTrans(T_ICR_t1(:, :, i),p5);
    v6(:, :, i) = RowVecTrans(T_ICR_t1(:, :, i),p6);
    v7(:, :, i) = RowVecTrans(T_ICR_t1(:, :, i),p7);
    v8(:, :, i) = RowVecTrans(T_ICR_t1(:, :, i),p8);
    
    if phiD(i) > DD 
        Location(:,:,i) = [p1;                           
                          p2;                             
                          p2;                             
                          p2;                            
                          v7(:,:,i);                     
                          v7(:, :, i);                   
                          v7(:, :, i);                   
                          v8(:, :, i)];
                      
    elseif phiD(i) > CC && phiD(i)<= DD
        Location(:,:,i) = [p1;                           
                          p2;                             
                          p2;                             
                          p2;                            
                          v6(:,:,i);                     
                          v6(:, :, i);                   
                          v7(:, :, i);                   
                          v8(:, :, i)]; 
                      
    elseif phiD(i) > BB   && phiD(i) <= CC 
        Location(:,:,i) = [p1;                            
                          p2;                             
                          p2;                             
                          p2;                             
                          v5(:, :, i);                    
                          v6(:, :, i);                   
                          v7(:, :, i);                   
                          v8(:, :, i)];                   
            
    elseif phiD(i) > AA && phiD(i) <= BB
        Location(:,:,i) = [p1;               
                          p2;            
                          p3;       
                          p3;       
                          v5(:, :, i);                             
                          v6(:, :, i);                             %Tibia contact
                          v7(:, :, i);                             %Tibia tendon contact
                          v8(:, :, i)];                            %patellar ligament ring

    elseif phiD(i) <= AA
        Location(:,:,i) = [p1;               %Origin
                          p2;            %BPA contacts mounting base
                          p3;       %femur channel contact
                          p4;       %femoral condyle contact
                          v5(:, :, i);                             %Tibia contact initial
                          v6(:, :, i);                             %Tibia contact
                          v7(:, :, i);                             %Tibia tendon contact
                          v8(:, :, i)];                            %patellar ligament ring
    else
    end
end
CrossPoint = 5;
Dia = 10;
rest = 0.520;
kmax = 0.432;
tendon = 0.030; 
fitting = 0.0254; 
pres = 596.4717;         %average pressure
Vas_Pam = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T_Pam, rest, kmax, tendon, fitting, pres);

%% Unstacking the Torques to identify specific rotations
Force1 = Vas_Int.Force + Vas_Lat.Force + Vas_Med.Force;
Torque1 = Vas_Int.Torque + Vas_Lat.Torque + Vas_Med.Torque;
TorqueR = Vas_Pam.Torque;

%% Add Torques from the Muscle Group
TorqueH = Torque1;


%% Plotting Torque Results

TorqueEx = zeros(size(TorqueH, 1), 1);
TorqueEy = zeros(size(TorqueH, 1), 1);
TorqueEz = zeros(size(TorqueH, 1), 1);

for i = 1:size(TorqueR, 1)
    if TorqueH(i, 1) >= 0
        TorqueEx(i) = TorqueR(i, 1) - TorqueH(i, 1);
    else
        TorqueEx(i) = TorqueH(i, 1) - TorqueR(i, 1);
    end
    
    if TorqueH(i, 2) >= 0
        TorqueEy(i) = TorqueR(i, 2) - TorqueH(i, 2);
    else
        TorqueEy(i) = TorqueH(i, 2) - TorqueR(i, 2);
    end
    
    if TorqueH(i, 3) >= 0
        TorqueEz(i) = TorqueR(i, 3) - TorqueH(i, 3);
    else
        TorqueEz(i) = TorqueH(i, 3) - TorqueR(i, 3);
    end
end

figure
hold on
sgtitle('Bicep Femoris Short Head Torque through Knee Flexion and Extension')

subplot(3, 2, 1)
plot(phiD, TorqueH(:, 3), phiD, TorqueR(:, 3))
legend('Human Muscle', 'Optimal BPA Location')
title('Muscle and PAM Z Torque')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')
legend('Human', 'PAM')

subplot(3, 2, 2)
plot(phiD, TorqueEz)
legend('Optimal PAM Location')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')
title('Adjusted Error Z Torque')

subplot(3, 2, 3)
plot(phiD, TorqueH(:, 2), phiD, TorqueR(:, 2))
title('Muscle and PAM Y Torque')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')
legend('Human', 'PAM')

subplot(3, 2, 4)
plot(phiD, TorqueEy)
legend('Optimal PAM Location')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')
title('Adjusted Error Y Torque')

subplot(3, 2, 5)
plot(phiD, TorqueH(:, 1), phiD, TorqueR(:, 1))
title('Muscle and PAM X Torque')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')
legend('Human', 'PAM')

subplot(3, 2, 6)
plot(phiD, TorqueEx)
legend('Optimal PAM Location')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')
title('Adjusted Error X Torque')

hold off




%% Plotting the angle between the vectors

aHR = zeros(size(TorqueH, 1), 1);
aHRH = zeros(size(TorqueH, 1), 1);

for i = 1:size(TorqueH, 1)
    uvecH = TorqueH(i, :)/norm(TorqueH(i, :));
                
    %Sometimes the BPA can't produce any force due to high
    %contraction. We will set it equal to negative the human
    %vector to maximize the penalty. Consider changing later
    if norm(TorqueR(i, :)) == 0
        uvecR = -uvecH;
    else
        uvecR = TorqueR(i, :)/norm(TorqueR(i, :));
    end
    
    aHR(i) = dot(uvecH, uvecR);
end


figure
hold on
title('Angle between the Human Torque Vector and PAM Torque Vectors')
plot(phiD, aHR)
legend('Human and Optimal PAM')
ylabel('Radians')
xlabel('Knee Angle, degree')
hold off

figure
plot(phiD, TorqueR(:, 3))
legend('BPA')
title('PAM Z Torque')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')

% HMuscleLocation = {Vas_Int.Location, Vas_Lat.Location, Vas_Med.Location};
% HMuscleCross = {Vas_Int.Cross, Vas_Lat.Cross, Vas_Med.Cross};
% 
% RMuscleLocation = {Vas_Pam.Location};
% RMuscleCross = {Vas_Pam.Cross};
% 
% Bones = {'Femur', 'Tibia'};
% 
% run("MuscleBonePlotting")

%% Compare to results
Load = [16.269	11.219	10.593	11.395	10.655	7.7816	3.958	1.909	1.22];     %Load in Newtons
K_ang = [-120	-108.5	-93	-70	-53.5	-49	-21	-9	-2]*c;      %Knee angle
LC_ang = [30.5	30	34	42.5	40	41.5	38	38	35]*c;      %Load Cell angle

d = 312.46/1000;
ang = -91.39;
p_rf = [d*cosd(ang), d*sind(ang), 0]';     %point of reaction force
T_t1_rf = RpToTrans(eye(3),p_rf);   %Tranformation matrix from theta 1 to reaction point
Trk = pagemtimes(TransInv(T_t1_rf),T_t1_ICR);
s1 = Trk(1,4,:);
s1 = squeeze(s1);
fcn15 = fit(phi',s1,'cubicspline');
s2 = Trk(2,4,:);
s2 = squeeze(s2);
fcn16 = fit(phi',s2,'cubicspline');

Trk = zeros(4,4,length(Load));
Fr = zeros(6,1,length(Load));
AdTrk = zeros(6,6,length(Load));
Fk = zeros(6,1,length(Load));

for i=1:length(Load)
    Trk(:,:,i) = RpToTrans(eye(3),[fcn15(K_ang(i)), fcn16(K_ang(i)), 0]');
    Fr(:,:,i) = -[0; 0; 0; Load(i)*cos(pi-LC_ang(i)); Load(i)*sin(pi-LC_ang(i)); 0];
    AdTrk(:,:,i) = Adjoint(Trk(:,:,i));
    Fk(:,:,i) = AdTrk(:,:,i)'*Fr(:,:,i);
    
end

TorqueZ = Fk(3,1,:);
TorqueZ = squeeze(TorqueZ);


figure
plot(phiD, TorqueR(:, 3), K_ang/c, TorqueZ,'o')
legend('Theoretical','Measured')
title('PAM Z Torque')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')