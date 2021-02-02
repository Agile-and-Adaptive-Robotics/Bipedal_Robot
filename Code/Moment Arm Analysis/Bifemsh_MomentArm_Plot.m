%% Freshen up the workspace
clc
clear
close all

%% Add paths to the muscle and pam calculators
addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Human_Data
addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Functions

%% Joint rotation transformation matrices
positions = 1000;

R = zeros(3, 3, positions);
T = zeros(4, 4, positions);

%Knee Extension and Flexion
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');

% kneeMin = -2.0943951;
% kneeMax = 0.17453293;

kneeMin = -120*pi/180;
kneeMax = 10*pi/180;
phi = linspace(kneeMin, kneeMax, positions);

for i = 1:positions
    hipToKnee = [fcn1(phi(i)), fcn2(phi(i)), 0];
    R(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    
    T(:, :, i) = RpToTrans(R(:, :, i), hipToKnee');
end

%% Muscle calculation
Name = 'Bicep Femoris (Short Head)';
MIF = 804;
OFL = 0.173; TSL = 0.089; Pennation = 0.40142573;
Location = [0.005, -0.211, 0.023;
            -0.03, -0.036, 0.029;
            -0.023, -0.056, 0.034];
CrossPoint = 2;
Bifemsh = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);


%% Creating Moment Arm Variables
BifemshMA = zeros(size(Bifemsh.MomentArm, 1), 1);

for i = 1:size(BifemshMA, 1)
    BifemshMA(i) = norm(Bifemsh.MomentArm(i, :));
    
    maTest(i) = Bifemsh.Torque(i, 3)/norm(Bifemsh.Force(i, :));
end

phiD = phi*180/pi;

figure
plot(phiD, BifemshMA)
title('My Calculated Moment Arm Magnitude')

for i = 1:(size(BifemshMA, 1)-1)
    ma(i) = (Bifemsh.MuscleLength(i)-Bifemsh.MuscleLength(i+1))/(phi(i)-phi(i+1));
end

op = xlsread('Bifemsh_MomentArm_OpenSim.xlsx');
opma = op(:, 3);
opD = op(:, 2);

figure
hold on
plot(phiD(1:end-1), -ma)
plot(opD, opma)
legend('Replicated Moment Arm Calculation', 'OpenSim Moment Arm')
title('OpenSim Moment Arm')
hold off

figure
hold on
plot(phiD, maTest)
plot(opD, opma)
legend('My Torque/Scalar Force', 'OpenSim Moment Arm')
title('OS Moment Arm Compared to Derived Moment Arm')
hold off

figure
plot(phiD, Bifemsh.Torque(:, 3))

TorqueZTest = -ma*norm(Bifemsh.Force(1:end-1, :));

figure
plot(phiD(1:end-1), TorqueZTest)

Force = Bifemsh.Force;

for i = 1:length(Force)
    normForce(i) = norm(Force(i, :));
end

op = xlsread('Bifemsh_Force_OpenSim.xlsx');
opf = op(:, 3);
% opD = op(:, 2);

figure
hold on
plot(phiD, normForce)
plot(opD, opf)
title('Muscle Force')
legend('My Calculated Force', 'OS Force')
hold off

op = xlsread('Bifemsh_Torque_OpenSim.xlsx');
opt = op(:, 3);

figure
plot(opD, opma)
title('OpenSim Moment Arm')

figure
plot(opD, opf)
title('OpenSim Force')

figure
plot(opD, opt)
title('OpenSim Torque')

opTcalc = opma.*opf;

figure
hold on
plot(opD, opt)
plot(opD, opTcalc)
title('OpenSim torque and recreation')
legend('Given T', 'Calced T')
hold off

op = xlsread('Bifemsh_TendonForce_OpenSim.xlsx');
optf = op(:, 3);

opTcalcTendon = opma.*optf;

figure
hold on
plot(opD, opt)
plot(opD, opTcalcTendon)
title('OpenSim Torque and Tendon Recreation')
legend('Given T', 'Calced T from Tendon Force')
hold off

figure
hold on
plot(opD, optf)
plot(phiD, normForce)
hold off

myT = -ma.*normForce(1:end - 1);

theirT = opma.*optf;

figure
hold on
plot(opD, opt)
plot(phiD(1:end-1), myT)
% plot(opD, theirT)
% legend('Given T', 'My Calced T', 'ma*tendon force')
hold off

figure
plot(opD, opma)

for i = 1:size(Bifemsh.Torque, 1)
    Torque(i) = norm(Bifemsh.Torque(i, :));
end

figure
plot(phiD, Torque)