
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

%Knee Extension and Flexion
%Human
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');
%Robot
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam =     [0.0065    0.0083    0.0094    0.0101    0.0120    0.0140    0.0161    0.0220    0.0269    0.0302    0.0311    0.0295    0.0253    0.0189    0.0109    0.0021]';
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam =     [-0.3981   -0.3968   -0.3961   -0.3957   -0.3949   -0.3943   -0.3941   -0.3950   -0.3982   -0.4034   -0.4098   -0.4165   -0.4227   -0.4273   -0.4297   -0.4289]';
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');

%Theta1 to ICR
t1_ICR_x = ([29.66	28.54	27.86	27.40	26.23	25.03	23.81	20.03	16.17	12.34	8.67	5.24	2.04	-1.01	-4.1	-7.58]')/1000;
fcn13 = fit(knee_angle,t1_ICR_x,'cubicspline');
t1_ICR_y = ([25.97	25.74	25.61	25.53	25.35	25.19	25.03	24.57	24.04	23.39	22.66	21.93	21.32	20.99	21.2	22.33]')/1000;
fcn14 = fit(knee_angle,theta1_ICR_y,'cubicspline');

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
    
    hipToKnee_Pam = [fcn3(phi(i)), fcn4(phi(i)), 0];
    R_Pam(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;   %Rotation matrix for robot
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    
    T_Pam(:, :, i) = RpToTrans(R_Pam(:, :, i), hipToKnee_Pam');     %Transformation matrix for robot
    
    t1toICR(1,:,i) = [fcn13(phi(i)), fcn14(phi(i)), 0];
    T_t1_ICR(:, :, i) = RpToTrans(R_Pam(:, :, i), t1toICR(1,:,i)');    
    T_ICR_t1(:, :, i) = TransInv(T_t1_ICR(:, :, i));
end

%% Muscle calculation
Name = 'Bicep Femoris (Short Head)';
MIF = 804;
OFL = 0.173; TSL = 0.089; Pennation = 0.40142573;
Location = zeros(3,3,positions);
for i = 1:positions
    Location(:,:,i) = [0.005, -0.211, 0.023;
            -0.03, -0.036, 0.029;
            -0.023, -0.056, 0.034];
end
CrossPoint = 2;
Bifemsh = MonoMuscleData(Name, Location, CrossPoint, MIF, TSL, Pennation, OFL, T);

%% PAM calculation
Name = 'Bicep Femoris (Short Head)';
CrossPoint = 2;

Location = zeros(2,3,positions);
for i = 1:positions
%Origin and Insertion from Ben
    Location(:,:,i) = [-0.050, 0.035, 0.0328;
            -0.02788, -0.04598, 0.0328];
end

%10 mm Festo
% Dia = 10;
% rest = 0.415;
% kmax = 0.350;
% tendon = 0.011; 
% fitting = 0.0254; 
% pres = 603.5236;         %average pressure
% Bifemsh_Pam = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T_Pam, rest, kmax, tendon, fitting, pres);
% 
% fitting = 0.0352; 
% Bifemsh_Pam_adj = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T_Pam, rest, kmax, tendon, fitting, pres);

%20 mm Festo
Dia = 20;
rest = 0.423;
kmax = 0.322;
tendon = 0.017; 
fitting = 0.0254; 
pres1 = 273.9783;         %average pressure, first test
pres2 = 484.8063;         %average pressure, first test
pres3 = 606.4926;         %average pressure, first test
Bifemsh_Pam1 = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T_Pam, rest, kmax, tendon, fitting, pres1);
Bifemsh_Pam2 = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T_Pam, rest, kmax, tendon, fitting, pres2);
Bifemsh_Pam3 = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T_Pam, rest, kmax, tendon, fitting, pres3);

fitting_adj = 0.0352; 
Bifemsh_Pam_adj1 = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T_Pam, rest, kmax, tendon, fitting_adj, pres1);
Bifemsh_Pam_adj2 = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T_Pam, rest, kmax, tendon, fitting_adj, pres2);
Bifemsh_Pam_adj3 = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T_Pam, rest, kmax, tendon, fitting_adj, pres3);

%% Unstacking the Torques to identify specific rotations
Torque1 = Bifemsh.Torque;
TorqueR = Bifemsh_Pam.Torque(:,:,4);
TorqueR_adj = Bifemsh_Pam_adj.Torque(:,:,4);

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
title('Muscle and PAM Z Torque')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Torque, N \cdot m','Interpreter','tex')
legend('Human', 'PAM')

subplot(3, 2, 2)
plot(phiD, TorqueEz)
legend('Optimal PAM Location')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Torque, N \cdot m','Interpreter','tex')
title('Adjusted Error Z Torque')

subplot(3, 2, 3)
plot(phiD, TorqueH(:, 2), phiD, TorqueR(:, 2))
title('Muscle and PAM Y Torque')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Torque, N \cdot m','Interpreter','tex')
legend('Human', 'PAM')

subplot(3, 2, 4)
plot(phiD, TorqueEy)
legend('Optimal PAM Location')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Torque, N \cdot m','Interpreter','tex')
title('Adjusted Error Y Torque')

subplot(3, 2, 5)
plot(phiD, TorqueH(:, 1), phiD, TorqueR(:, 1))
title('Muscle and PAM X Torque')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Torque, N \cdot m','Interpreter','tex')
legend('Human', 'PAM')

subplot(3, 2, 6)
plot(phiD, TorqueEx)
legend('Optimal PAM Location')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Torque, N \cdot m','Interpreter','tex')
title('Adjusted Error X Torque')

hold off

%% Compare Expected vs Adjusted PAM values
figure
plot(phiD, Bifemsh_Pam_adj.Torque(:, 3), phiD, TorqueR(:, 3))
title('Muscle and PAM Z Torque')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Torque, N \cdot m','Interpreter','tex')
legend('Optimized', 'Original Theoretical')


%% Plotting muscle lengths and moment arms using two different moment arm
%calculations
ML = Bifemsh.MuscleLength;
PamL = Bifemsh_Pam.MuscleLength;
for i = 1:size(Bifemsh.MomentArm,1)
    MA(i,:) = norm(Bifemsh.MomentArm(i,1:2));               %Muscle moment arm, Z axis
    BPAma(i,:) = norm(Bifemsh_Pam.MomentArm(i,1:2));        %BPA moment arm, Z axis
end
dM = diff(Bifemsh.MuscleLength);           %Muscle length difference
dP = diff(Bifemsh_Pam.MuscleLength);       %PAM length difference
dO = diff(phiD);                           %Angle difference

figure
hold on
sgtitle('Bicep Femoris Short Head Length and Moment Arm through Knee Flexion and Extension')

subplot(2, 2, 1)
plot(phiD, ML, phiD, PamL)
title('Muscle and PAM Lengths')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Length, m')
legend('Human', 'PAM')

subplot(2, 2, 2)
plot(phiD, MA, phiD, BPAma)
title('Moment arm, Z axis, vector method')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Length, m')
legend('Human', 'PAM')

subplot(2, 2, 3)
plot(phiD(1:99), -dM./dO', phiD(1:99), -dP./dO')
title('Moment arm, Z axis, left difference method')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Length, m')
legend('Human', 'PAM')

subplot(2, 2, 4)
plot(phiD(2:100), -dM./dO', phiD(2:100), -dP./dO')
title('Moment arm, Z axis, right difference method')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Length, m')
legend('Human', 'PAM')

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

%% Plotting on the Mesh Skeleton

HMuscleLocation = {Bifemsh.Location};
HMuscleCross = {Bifemsh.Cross};

RMuscleLocation = {Bifemsh_Pam.Location};
RMuscleCross = {Bifemsh_Pam.Cross};

Bones = {'Femur', 'Tibia'};

run("MuscleBonePlotting")
