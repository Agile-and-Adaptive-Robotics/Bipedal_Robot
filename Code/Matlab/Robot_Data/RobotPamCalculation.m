%% Robot Pam Data
% This script is a collection of all of the proposed pams for the bipedal
% robot. The script calculates the amount of torque that each on can
% generate about a joint.
clear
clc
close all


addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Functions

%% ------------- Hip Muscles ----------------

%% Abduction & Adduction
iteration = 2;
R = zeros(3, 3, iteration);
T = zeros(4, 4, iteration);
adducMax = 30*pi/180;
adducMin = -55*pi/180;
theta = linspace(adducMin, adducMax, iteration);
pelvisToHip = [-0.0707, -0.0661, 0.0835];

% tibiaToPelvis = [

for i = 1:iteration
    R(:, :, i) = [cos(theta(i)), -sin(theta(i)), 0;
                    sin(theta(i)), cos(theta(i)), 0;
                    0, 0, 1];
    
    T(:, :, i) = RpToTrans(R(:, :, i), pelvisToHip');
end

Name = 'Gluteus Maximus';
Location = [-0.119, 0.061, 0.07;
            -0.046, -0.025, 0.039;
            -0.028, -0.057, 0.047];
CrossPoint = 2;
Dia = 10;
Glut_Max1_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T); 

Name = 'Gluteus Medius';
Location = [-0.041, 0.03, 0.121;
            -0.022, -0.012, 0.056];
CrossPoint = 2;
Dia = 20;
Glut_Med_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);

Name = 'Periformis';
Location = [-0.14, 0, 0.024;
            -0.015, -0.004, 0.44];
CrossPoint = 2;
Dia = 10;
Peri_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);

Name = 'Adductor Brevis';
Location = [-0.059, -0.091, 0.016;
            0.001, -0.12, 0.029];
CrossPoint = 2;
Dia = 10;
Add_Brev_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);

%Knee Extension and Flexion
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');

kneeMin = -2.0943951;
kneeMax = 0.17453293;
phi = linspace(kneeMin, kneeMax, iteration);

for i = 1:iteration
    hipToKnee = [fcn1(phi(i)), fcn2(phi(i)), 0];
    R(:, :, i, 2) = [cos(theta(i)), -sin(theta(i)), 0;
                    sin(theta(i)), cos(theta(i)), 0;
                    0, 0, 1];
    
    T(:, :, i, 2) = RpToTrans(R(:, :, i, 2), hipToKnee');
end

Name = 'Sartorius';
Location = [-0.015, -0.001, 0.124;
            -0.006, -0.042, -0.04;
            0.024, -0.084, -0.0235];
CrossPoint = [2, 2];
Dia = 10;
Sar_Pam = BiPamData(Name, Location, CrossPoint, Dia, T);

%% Hip Flexion and Extension
Name = 'Iliacus';
Location = [-0.067, 0.036, 0.085;
            -0.029, -0.081, 0.082;
            -0.019, -0.062, 0.013];
CrossPoint = 3;
Dia = 40;
Iliacus_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);

Name = 'Tensory Fascia Latae';
Location = [-0.031, 0.021, 0.124;
            0.005, -0.405, 0.036;
            0.006, -0.049, 0.03];
CrossPoint = [2, 3];
Dia = 10;
TFL_Pam = BiPamData(Name, Location, CrossPoint, Dia, T);

Name = 'Rectus Femoris (Quadriceps)';   %The Location should be updated to reflect the moving patella points. For now it is static in the knee frame
Location = [-0.029, -0.031, 0.097;
            0.033, -0.403, 0.002;
            0.062, 0.021, 0.0014];
CrossPoint = [2, 3];
Dia = 20;
Rect_Fem_Pam = BiPamData(Name, Location, CrossPoint, Dia, T);

Name = 'Gracilis';
Location = [-0.074, -0.119, 0.028;
            0.006, -0.084 -0.023];
CrossPoint = [2, 2];
Dia = 10;
Grac_Pam = BiPamData(Name, Location, CrossPoint, Dia, T);

Name = 'Adductor Magnus';
Location = [-0.083, -0.119, 0.031;
            0.005, -0.229, 0.023];
CrossPoint = 2;
Dia = 10;
Add_Mag_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);

Name = 'Bicep Femoris (Long Head)';
Location = [-0.126, -0.103, 0.069;
            -0.023, -0.056, 0.034];
CrossPoint = [2, 2];
Dia = 20;
Bifemlh_Pam = BiPamData(Name, Location, CrossPoint, Dia, T);

Name = 'Semimembranosus';
Location = [-0.119, -0.097, 0.072;
            -0.027, -0.048, -0.02];
CrossPoint = [2, 2];
Dia = 20;
Semimem_Pam = BiPamData(Name, Location, CrossPoint, Dia, T);

%% Hip Internal and External Rotation
Name = 'Quadricep Femoris';             %I don't think thats the actual name for this muscle
Location = [-0.114, -0.115, 0.052;
            -0.038, -0.036, 0.037];
CrossPoint = 2;
Dia = 20;
Quad_Fem_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);


%% Knee Muscles Specifically
Name = 'Bicep Femoris (Short Head)';
Location = [0.005, -0.211, 0.023;
            -0.023, -0.056, 0.034];
CrossPoint = 2;
Dia = 20;
Bifemsh_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);

Name = 'Vastus Intermedius';
Location = [0.029, -0.192, 0.031;
            0.034, -0.403, 0.005;
            0.0555, 0.025, 0.0018];
CrossPoint = 3;
Dia = 40;
Vas_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);

%Biarticular with Femur to Ankle
% Name = 'Lateral Gastrocnemius';
% MIF = 683;
% Location = [-0.022, -0.395, 0.027;
%             -0.03, -0.402, 0.027;
%             0, 0.031, -0.005];
% CrossPoint = [3, 3, 3];
% Lat_Gas = MuscleData(Name, Location, CrossPoint, MIF, T);
% 
% Name = 'Medial Gastrocnemius';
% MIF = 1558;
% Location = [-0.019, -0.393, -0.024;
%             -0.03, -0.402, -0.026;
%             0, 0.031, -0.005];
% CrossPoint = [3, 3, 3];
% Med_Gas = MuscleData(Name, Location, CrossPoint, MIF, T);