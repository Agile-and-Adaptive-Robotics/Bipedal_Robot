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
knee_x_Pam =     [0.0010	0.0027	0.0038	0.0045	0.0064	0.0084	0.0105	0.0164	0.0213	0.0246	0.0255	0.0239	0.0197	0.0132	0.0052	-0.0036]';
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam =     [-0.3982	-0.3969	-0.3962	-0.3958	-0.3950	-0.3944	-0.3942	-0.3951	-0.3984	-0.4035	-0.4099	-0.4167	-0.4228	-0.4274	-0.4298	-0.4292]';
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');


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
    
    T(:, :, i) = RpToTrans(R(:, :, i), hipToKnee');  %Human Transformation matrix
    
    hipToKnee_Pam = [fcn3(phi(i)), fcn4(phi(i)), 0];
    R_Pam(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;   %Rotation matrix for robot
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    
    T_Pam(:, :, i) = RpToTrans(R_Pam(:, :, i), hipToKnee_Pam');     %Transformation matrix for robot
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

%% PAM calculation
Name = 'Bicep Femoris (Short Head)';
CrossPoint = 2;
Dia = 20;

Bifemsh_PamH = MonoPamData(Name, Location, CrossPoint, Dia, T);

Location = [-0.050, -0.045, 0.0328;
            -0.03239, -0.08217, 0.0328];
Bifemsh_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T_pam);

%% Unstacking the Torques to identify specific rotations
Torque1 = Bifemsh.Torque;
TorqueR = Bifemsh_Pam.Torque;

%% Add Torques from the Muscle Group
TorqueH = Torque1;

%% First cost function calculation
C = costFunctionKnee(TorqueH, TorqueR);

%% Generating new points for the PAM based on the bone mesh
Femur = xlsread('Femur_Mesh_Points.xlsx');
Tibia = xlsread('Tibia_Mesh_Points.xlsx');

fprintf('The algorithm will be calculating Torque between %d different mesh locations.\n', size(Tibia, 1)*size(Femur, 1))

meshTracker = [0, 0];

originalLocation = Location;

iC = 1;                     %Index variable for the cost function
CMaxPrev = 10^5;

for i = 1:size(Tibia, 1)
    for ii = 1:size(Femur, 1)
        clc
        fprintf('%d \t of %d \n', iC, size(Tibia, 1)*size(Femur, 1))
        Location(1, :) = Femur(ii, :);
        Location(2, :) = Tibia(i, :);
        
        Bifemsh_Pam.Location = Location;
%         Bifemsh_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);
        TorqueR = Bifemsh_Pam.Torque;
        
        C = costFunctionKnee(TorqueH, TorqueR);
        
        if C < CMaxPrev
            if isequal(Bifemsh_Pam.LengthCheck, 'Usable')
                Tracker = [i, ii];
                CMaxPrev = C;
            end
        end

        iC = iC + 1;
    end
end

% %% Refining the Mesh Search
if exist('Tracker', 'var') == 0
    Location = originalLocation;
else
    Location(1, :) = Femur(Tracker(2), :);
    Location(2, :) = Tibia(Tracker(1), :);
end

bestMeshLocation = Location;

neg = [-1 1];
epsilon = 0.01;
refine = 0.75;

LocationTracker = Location;

%Change this to a while loop of reducing epsilon after testing 
for i = 1:10
    startingLocation = LocationTracker;
    for k1 = 1:2
        ep1(1) = epsilon*neg(k1);
    for k2 = 1:2
        ep1(2) = epsilon*neg(k2);
    for k3 = 1:2
        ep1(3) = epsilon*neg(k3);
    for k4 = 1:2
        ep2(1) = epsilon*neg(k4);
    for k5 = 1:2
        ep2(2) = epsilon*neg(k5);
    for k6 = 1:2
        ep2(3) = epsilon*neg(k6);
        
        Location(1, :) = startingLocation(1, :) + ep1;
        Location(2, :) = startingLocation(2, :) + ep2;
        
        Bifemsh_Pam.Location = Location;
        
%         Bifemsh_Pam = MonoPamData(Name, Location, CrossPoint, Dia, T);
        TorqueR = Bifemsh_Pam.Torque;
        
        C = costFunctionKnee(TorqueH, TorqueR);
        
        if C < CMaxPrev
            if isequal(Bifemsh_Pam.LengthCheck, 'Usable')
                LocationTracker = Location;
                CMaxPrev = C;
            end
        end
        
    end
    end
    end
    end
    end
    end
    if LocationTracker == startingLocation
        epsilon = epsilon*refine;
    end
end

%% Plotting Torque Results

% Bifemsh_Pam.Location = LocationTracker;
TorqueR = Bifemsh_Pam.Torque;

TorqueRH = Bifemsh_PamH.Torque;

phiD = phi*180/pi;

TorqueEx = zeros(size(TorqueH, 1), 1);
TorqueEHx = zeros(size(TorqueH, 1),1 );
TorqueEy = zeros(size(TorqueH, 1), 1);
TorqueEHy = zeros(size(TorqueH, 1), 1);
TorqueEz = zeros(size(TorqueH, 1), 1);
TorqueEHz = zeros(size(TorqueH, 1), 1);

for i = 1:size(TorqueR, 1)
    if TorqueH(i, 1) >= 0
        TorqueEx(i) = TorqueR(i, 1) - TorqueH(i, 1);
        TorqueEHx(i) = TorqueRH(i, 1) - TorqueH(i, 1);
    else
        TorqueEx(i) = TorqueH(i, 1) - TorqueR(i, 1);
        TorqueEHx(i) = TorqueH(i, 1) - TorqueRH(i, 1);
    end
    
    if TorqueH(i, 2) >= 0
        TorqueEy(i) = TorqueR(i, 2) - TorqueH(i, 2);
        TorqueEHy(i) = TorqueRH(i, 2) - TorqueH(i, 2);
    else
        TorqueEy(i) = TorqueH(i, 2) - TorqueR(i, 2);
        TorqueEHy(i) = TorqueH(i, 2) - TorqueRH(i, 2);
    end
    
    if TorqueH(i, 3) >= 0
        TorqueEz(i) = TorqueR(i, 3) - TorqueH(i, 3);
        TorqueEHz(i) = TorqueRH(i, 3) - TorqueH(i, 3);
    else
        TorqueEz(i) = TorqueH(i, 3) - TorqueR(i, 3);
        TorqueEHz(i) = TorqueH(i, 3) - TorqueRH(i, 3);
    end
end

figure
hold on
sgtitle('Bicep Femoris Short Head Torque through Knee Flexion and Extension')

subplot(3, 2, 1)
plot(phiD, TorqueH(:, 3), phiD, TorqueR(:, 3), phiD, TorqueRH(:, 3))
legend('Human Muscle', 'Optimal BPA Location', 'BPA at Human Locations')
title('Muscle and PAM Z Torque')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')
legend('Human', 'PAM', 'PAM at Human Locations', 'Location', 'best')

subplot(3, 2, 2)
plot(phiD, TorqueEz, phiD, TorqueEHz)
legend('Optimal PAM Location', 'At Human Locations')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')
title('Adjusted Error Z Torque')

subplot(3, 2, 3)
plot(phiD, TorqueH(:, 2), phiD, TorqueR(:, 2), phiD, TorqueRH(:, 2))
title('Muscle and PAM Y Torque')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')
legend('Human', 'PAM', 'PAM at Human Locations')

subplot(3, 2, 4)
plot(phiD, TorqueEy, phiD, TorqueEHy)
legend('Optimal PAM Location', 'At Human Locations')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')
title('Adjusted Error Y Torque')

subplot(3, 2, 5)
plot(phiD, TorqueH(:, 1), phiD, TorqueR(:, 1), phiD, TorqueRH(:, 1))
title('Muscle and PAM X Torque')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')
legend('Human', 'PAM', 'PAM at Human Locations')

subplot(3, 2, 6)
plot(phiD, TorqueEx, phiD, TorqueEHx)
legend('Optimal PAM Location', 'At Human Locations')
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
    
    if norm(TorqueRH(i, :)) == 0
        uvecRH = -uvecH;
    else
        uvecRH = TorqueRH(i, :)/norm(TorqueRH(i, :));
    end
    
    aHR(i) = dot(uvecH, uvecR);
    aHRH(i) = dot(uvecH, uvecRH);
end


figure
hold on
title('Angle between the Human Torque Vector and PAM Torque Vectors')

plot(phiD, aHR, phiD, aHRH)
legend('Human and Optimal PAM', 'Human and PAM at Human Locations')
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
