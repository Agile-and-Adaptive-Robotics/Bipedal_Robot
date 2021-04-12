%% Human Muscle Data
% This script is a collection of all of the human muscles and ways of
% generating torque information for each of those muscles, based on the
% type of actuation they create
clear
clc
close all


addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Functions

%% ------------- Hip Muscles ----------------

%% Abduction & Adduction
iteration = 100;
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

Name = 'Gluteus Maximus 1';
MIF = 573;
OFL = 0.142; TSL = 0.125; Pennation = 0.08726646;
Location = [-0.119, 0.061, 0.07;
            -0.129, 0.001, 0.089;
            -0.046, -0.025, 0.039;
            -0.028, -0.057, 0.047];
CrossPoint = 3;
Glut_Max1 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T); 

Name = 'Gluteus Maximus 2';
MIF = 819.0;
OFL = 0.147; TSL = 0.127; Pennation = 0.0;
Location = [-0.135, 0.018, 0.056;
            -0.138, -0.052, 0.091;
            -0.043, -0.053, 0.029;
            -0.016, -0.102, 0.042];
CrossPoint = 3;
Glut_Max2 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T); 

Name = 'Gluteus Maximus 3';
MIF = 552.0;
OFL = 0.144; TSL = 0.145; Pennation = 0.08726646;
Location = [-0.156, -0.031, 0.006;
            -0.153, -0.105, 0.04;
            -0.03, -0.104, 0.013;
            -0.006, -0.142, 0.041];
CrossPoint = 3;
Glut_Max3 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T); 

Name = 'Gluteus Medius 1';
MIF = 819.0;
OFL = 0.0535; TSL = 0.078; Pennation = 0.13962634;
Location = [-0.041, 0.03, 0.121;
            -0.022, -0.012, 0.056];
CrossPoint = 2;
Glut_Med1 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Gluteus Medius 2';
MIF = 573;
OFL = 0.0845; TSL = 0.053; Pennation = 0.0;
Location = [-0.086, 0.044, 0.077;
            -0.026, -0.006, 0.053];
CrossPoint = 2;
Glut_Med2 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Gluteus Medius 3';
MIF = 653;
OFL = 0.0646; TSL = 0.053; Pennation = 0.33161256;
Location = [-0.122, 0.011, 0.065;
            -0.031, -0.005, 0.052];
CrossPoint = 2;
Glut_Med3 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Gluteus Minimus 1';
MIF = 270;
OFL = 0.068; TSL = 0.016; Pennation = 0.17453293;
Location = [-0.047, -0.008, 0.106;
            -0.007, -0.01, 0.056];
CrossPoint = 2;
Glut_Min1 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Gluteus Minimus 2';
MIF = 285;
OFL = 0.056; TSL = 0.026; Pennation = 0.0;
Location = [-0.063, -0.006, 0.099;
            -0.01, -0.01, 0.056];
CrossPoint = 2;
Glut_Min2 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Gluteus Minimus 3';
MIF = 323;
OFL = 0.038; TSL = 0.051; Pennation = 0.01745329;
Location = [-0.083, -0.006, 0.086;
            -0.013, -0.008, 0.055];
CrossPoint = 2;
Glut_Min3 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Piriformis';
MIF = 444;
OFL = 0.026; TSL = 0.115; Pennation = 0.17453293;
Location = [-0.14, 0, 0.024;
            -0.119, -0.028, 0.066;
            -0.015, -0.004, 0.44];
CrossPoint = 3;
Peri = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Adductor Brevis';
MIF = 429;
OFL = 0.133; TSL = 0.02; Pennation = 0.0;
Location = [-0.059, -0.091, 0.016;
            0.001, -0.12, 0.029];
CrossPoint = 2;
Add_Brev = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Adductor Longus';
MIF = 627;
OFL = 0.138; TSL = 0.11; Pennation = 0.10471976;
Location = [-0.032, -0.084, 0.017;
            0.005, -0.211, 0.023];
CrossPoint = 2;
Add_Long = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Pectineus';
MIF = 266;
OFL = 0.1; TSL = 0.033; Pennation = 0.0;
Location = [-0.043, -0.077, 0.045;
            -0.012, -0.082, 0.025];
CrossPoint = 2;
Pect = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

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
    R(:, :, i, 2) = [cos(phi(i)), -sin(phi(i)), 0;
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    
    T(:, :, i, 2) = RpToTrans(R(:, :, i, 2), hipToKnee');
end

Name = 'Sartorius';
MIF = 156;
OFL = 0.52; TSL = 0.1; Pennation = 0.0;
Location = [-0.015, -0.001, 0.124;
            -0.003, -0.357, -0.042;
            -0.006, -0.042, -0.04;
            0.006, -0.059, -0.038;
            0.024, -0.084, -0.0235];
CrossPoint = [2, 3];
Sar = BiMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

%% Hip Flexion and Extension
Name = 'Iliacus';
MIF = 1073;
OFL = 0.1; TSL = 0.1; Pennation = 0.12217305;
Location = [-0.067, 0.036, 0.085;
            -0.026, -0.055, 0.081;
            -0.029, -0.081, 0.082;
            0.002, -0.054, 0.006;
            -0.019, -0.062, 0.013];
CrossPoint = 4;
Iliacus = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Psoas';
MIF = 1113;
OFL = 0.1; TSL = 0.16; Pennation = 0.13962634;
Location = [-0.065, 0.089, 0.029;
            -0.024, -0.057, 0.076;
            -0.029, -0.081, 0.082;
            0.002, -0.051, 0.004;
            -0.019, -0.06, 0.01];
CrossPoint = 4;
Psoas = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Tensory Fascia Latae';
MIF = 233;
OFL = 0.095; TSL = 0.425; Pennation = 0.05235988;
Location = [-0.031, 0.021, 0.124;
            0.029, -0.1, 0.06;
            0.005, -0.405, 0.036;
            0.006, -0.049, 0.03];
CrossPoint = [2, 4];
TFL = BiMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Rectus Femoris (Quadriceps)';   %The Location should be updated to reflect the moving patella points. For now it is static in the knee frame
MIF = 1169;
OFL = 0.114; TSL = 0.31; Pennation = 0.08726646;
Location = [-0.029, -0.031, 0.097;
            0.033, -0.403, 0.002;
            0.062, 0.021, 0.0014];
CrossPoint = [2, 3];
Rect_Fem = BiMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Gracilis';
MIF = 162;
OFL = 0.352; TSL = 0.126; Pennation = 0.05235988;
Location = [-0.072, -0.119, 0.028;
            -0.027, -0.032, -0.038;
            -0.019, -0.052, -0.036;
            0.006, -0.084 -0.023];
CrossPoint = [2, 2];
Grac = BiMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Adductor Magnus 1';
MIF = 381;
OFL = 0.087; TSL = 0.06; Pennation = 0.08726646;
Location = [-0.073, -0.117, 0.025;
            -0.004, -0.121, 0.034];
CrossPoint = 2;
Add_Mag1 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Adductor Magnus 2';
MIF = 343;
OFL = 0.121; TSL = 0.12; Pennation = 0.05235988;
Location = [-0.083, -0.119, 0.031;
            0.005, -0.229, 0.023];
CrossPoint = 2;
Add_Mag2 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Adductor Magnus 3';
MIF = 488;
OFL = 0.131; TSL = 0.249; Pennation = 0.08726646;
Location = [-0.111, -0.114, 0.049;
            0.007, -0.384, -0.027];
CrossPoint = 2;
Add_Mag3 = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Bicep Femoris (Long Head)';
MIF = 896;
OFL = 0.109; TSL = 0.326; Pennation = 0.0;
Location = [-0.126, -0.103, 0.069;
            -0.03, -0.036, 0.029;
            -0.023, -0.056, 0.034];
CrossPoint = [2, 2];
Bifemlh = BiMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Semimembranosus';
MIF = 1288;
OFL = 0.08; TSL = 0.359; Pennation = 0.26179939;
Location = [-0.119, -0.097, 0.072;
            -0.035, -0.035, -0.019;
            -0.027, -0.048, -0.02];
CrossPoint = [2, 2];
Semimem = BiMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Semitendinosus';
MIF = 410;
OFL = 0.201; TSL = 0.2555; Pennation = 0.08726646;
Location = [-0.126, -0.11, 0.06;
            -0.042, -0.029, -0.023;
            -0.033, -0.053, -0.023;
            -0.011, -0.075, -0.025;
            0.003, -0.096, -0.019];
CrossPoint = [2, 2];
Semiten = BiMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

%% Hip Internal and External Rotation
Name = 'Gemellus';
MIF = 164;
OFL = 0.024; TSL = 0.039; Pennation = 0.0;
Location = [-0.113, -0.082, 0.071;
            -0.014, -0.003, 0.044];
CrossPoint = 2;
Gem = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Quadricep Femoris';             %I don't think thats the actual name for this muscle
MIF = 381;
OFL = 0.054; TSL = 0.024; Pennation = 0.0;
Location = [-0.114, -0.115, 0.052;
            -0.038, -0.036, 0.037];
CrossPoint = 2;
Quad_Fem = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);


%% Knee Muscles Specifically
Name = 'Bicep Femoris (Short Head)';
MIF = 804;
OFL = 0.173; TSL = 0.089; Pennation = 0.40142573;
Location = [0.005, -0.211, 0.023;
            -0.03, -0.036, 0.029;
            -0.023, -0.056, 0.034];
CrossPoint = 2;
Bifemsh = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Vastus Intermedius';
MIF = 1365;
OFL = 0.087; TSL = 0.136; Pennation = 0.05235988;
Location = [0.029, -0.192, 0.031;
            0.034, -0.208, 0.029;
            0.034, -0.403, 0.005;
            0.0555, 0.025, 0.0018];
CrossPoint = 4;
Vas_Int = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Vastus Lateralis';
MIF = 1871;
OFL = 0.084; TSL = 0.157; Pennation = 0.08726646;
Location = [0.005, -0.185, 0.035;
            0.027, -0.259, 0.041;
            0.036, -0.403, 0.021;
            0.025, -0.424, 0.018;
            0.06, 0.02, 0.0165];
CrossPoint = 5;
Vas_Lat = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

Name = 'Vastus Medialis';
MIF = 1294;
OFL = 0.089; TSL = 0.126; Pennation = 0.08726646;
Location = [0.014, -0.21, 0.019;
            0.036, -0.277, 0.001;
            0.037, -0.405, -0.013;
            0.027, -0.425, -0.013;
            0.05625, 0.022, -0.0146];
CrossPoint = 5;
Vas_Med = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

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




