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
knee_x_Pam =     [0.0010	0.0027	0.0038	0.0045	0.0064	0.0084	0.0105	0.0164	0.0213	0.0246	0.0255	0.0239	0.0197	0.0132	0.0052	-0.0036]'; 
fcn11 = fit(knee_angle,knee_x_Pam,'cubicspline');  %Tibia reference frame
knee_y_Pam =     [-0.3982	-0.3969	-0.3962	-0.3958	-0.3950	-0.3944	-0.3942	-0.3951	-0.3984	-0.4035	-0.4099	-0.4167	-0.4228	-0.4274	-0.4298	-0.4292]';
fcn12 = fit(knee_angle,knee_y_Pam,'cubicspline');  %Tibia reference frame
%Patella locations in Femur frame
P_T_x = [0.06055	0.06	0.05956	0.06094	0.05953	0.06072	0.05927	0.05982	0.05955	0.05228	0.04206	0.03099	0.016	-0.00064	-0.01722	-0.03161]'; %Patella Top, x location
fcn13 = fit(knee_angle,P_T_x,'cubicspline');
P_T_y = [-0.37779	-0.38615	-0.38366	-0.3857	-0.38907	-0.39348	-0.39766	-0.41177	-0.42227	-0.43372	-0.44375	-0.45354	-0.45795	-0.45839	-0.45529	-0.44819]'; %Patella Top, y location
fcn14 = fit(knee_angle,P_T_y,'cubicspline');
P_B_x = [0.05898	0.0585	0.05891	0.05828	0.05822	0.05856	0.05842	0.04735	0.038	0.0264	0.01175	-0.00467	-0.02223	-0.03997	-0.05465	-0.06578]'; %Patella Bottom, x location
fcn15 = fit(knee_angle,P_B_x,'cubicspline');
P_B_y = [-0.41712	-0.42073	-0.42306	-0.42447	-0.42847	-0.43298	-0.43779	-0.44568	-0.455	-0.46374	-0.46981	-0.46987	-0.46785	-0.46119	-0.44844	-0.4316]'; %Patella Bottom, y location
fcn16 = fit(knee_angle,P_B_y,'cubicspline');

kneeMin = -2.0943951;
kneeMax = 0.17453293;
phi = linspace(kneeMin, kneeMax, positions);
%We want one of our positions to be home position, so let's make the
%smallest value of phi equal to 0
[val, pos] = min(abs(phi));
phi(pos) = 0;


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
end

%% Muscle calculation
Name = 'Vastus Medialis';
MIF = 1294;
OFL = 0.089; TSL = 0.126; Pennation = 0.08726646;
Location = [0.014, -0.21, 0.019;
            0.036, -0.277, 0.001;
            0.037, -0.405, -0.013;
            0.027, -0.425, -0.013;
            fcn5(phi(i)), fcn6(phi(i)) -0.0146];
CrossPoint = 5;
Vas_Med = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Vastus Intermedius';
MIF = 1365;
OFL = 0.087; TSL = 0.136; Pennation = 0.05235988;
Location = [0.029, -0.192, 0.031;
            0.034, -0.208, 0.029;
            0.034, -0.403, 0.005;
            fcn7(phi(i)), fcn8(phi(i)) 0.0018];
CrossPoint = 4;
Vas_Int = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Vastus Lateralis';
MIF = 1871;
OFL = 0.084; TSL = 0.157; Pennation = 0.08726646;
Location = [0.005, -0.185, 0.035;
            0.027, -0.259, 0.041;
            0.036, -0.403, 0.021;
            0.025, -0.424, 0.018;
            fcn9(phi(i)), fcn10(phi(i)) 0.0165];
CrossPoint = 5;
Vas_Lat = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

%% PAM calculation
Name = 'Vastus Intermedius';

% Origin Location from Ben
Location = [0.030, -0.050, 0;
            0.048, -0.349, 0.000;               %BPA contacts head of socket head cap screw that joins Femur to Femoral end
            fcn13(phi(i)), fcn14(phi(i)), 0;    %Top of patella, as a function of knee angle
            fcn15(phi(i)), fcn16(phi(i)), 0;    %Bottom of patella, as a function of knee angle
            0.04261, 0.07741, 0.000];           %Top of patellar ligament bracket
CrossPoint = 5;
Dia = 40;
Vas_Pam = MonoPamDataPhysicalExtensor(Name, Location, CrossPoint, Dia, T_Pam);

%% Unstacking the Torques to identify specific rotations
Force1 = Vas_Int.Force + Vas_Lat.Force + Vas_Med.Force;
Torque1 = Vas_Int.Torque + Vas_Lat.Torque + Vas_Med.Torque;
TorqueR = Vas_Pam.Torque;

%% Add Torques from the Muscle Group
TorqueH = Torque1;


%% Plotting Torque Results
phiD = phi*180/pi;

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

HMuscleLocation = {Vas_Int.Location, Vas_Lat.Location, Vas_Med.Location};
HMuscleCross = {Vas_Int.Cross, Vas_Lat.Cross, Vas_Med.Cross};

RMuscleLocation = {Vas_Pam.Location};
RMuscleCross = {Vas_Pam.Cross};

Bones = {'Femur', 'Tibia'};

run("MuscleBonePlotting")
