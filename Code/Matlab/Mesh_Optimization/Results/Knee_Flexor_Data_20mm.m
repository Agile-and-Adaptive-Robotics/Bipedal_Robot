

%% Mesh Points Calculation - Bicep Femoris Short Head
% This code will calculate the torque difference between all of the points
% from one bone mesh to another to determine the best location for muscle
% placement

%% Freshen up the workspace
clc
clear
close all

set(groot, ...
    'defaultAxesFontName','Arial', ...
    'defaultAxesFontSize',10, ...
    'defaultAxesFontWeight','bold', ...
    'defaultAxesLabelFontSizeMultiplier',1, ...
    'defaultAxesTitleFontSizeMultiplier',1.2, ...
    'defaultAxesLineWidth',2, ...
    'defaultAxesBox','off', ...
    'defaultAxesXMinorTick','on', ...
    'defaultAxesYMinorTick','on', ...
    'defaultAxesTickLength',[0.025 0.05], ...
    'defaultAxesXGrid','off', ...
    'defaultAxesYGrid','off', ...
    'defaultLineLineWidth',2, ...
    'defaultLegendFontName','Arial', ...
    'defaultLegendFontSize',8, ...
    'defaultLegendFontWeight','bold', ...
    'defaultLegendBox','off')

%% Add paths to the muscle and pam calculators
% Repo root and path setup (same block as Knee_Extensor_20mm.m): resolves
% MuscleBonePlotting/AnimateKneeBoneMuscle (Bone_Mesh_Plots), the bone
% meshes in Open_Sim_Bone_Geometry, Colors.m, and the
% minimizeFlxPin10_results_20260730 mat + OpenSim_Bifem txt files in
% Testing_Data\2022_02_Festo, regardless of cwd.
scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root)
        error('Could not locate the Bipedal_Robot repo root from %s', scriptDir)
    end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
% Mesh_Optimization must win any shadowing contest against data subfolders.
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
% Append (do not prepend) so Code\Matlab keeps winning name collisions.
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

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
%Human
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');
%Robot
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam =     ([23.30	22.22	21.55	21.09	19.91	18.70	17.48	13.82	10.44	7.60	5.52	4.35	4.16	5.01	7.04	10.47]')/1000;
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam =     ([-416.65	-417.03	-417.19	-417.28	-417.41	-417.41	-417.30	-416.28	-414.36	-411.72	-408.62	-405.32	-402.08	-399.16	-396.85	-395.66]')/1000;
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');

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
    
    t1toICR(1,:,i) = [fcn13(phi(i)), fcn14(phi(i)), 0]; %distance from theta1 to ICR
    T_t1_ICR(:, :, i) = RpToTrans(eye(3), t1toICR(1,:,i)');    %transform from the ICR frame to theta1
    T_ICR_t1(:, :, i) = RpToTrans(eye(3), -t1toICR(1,:,i)');    %transform from t1 frame to ICR
end

phiD = phi*180/pi;

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

%Original origin and insertion from Ben
p10 = [-0.050, 0.035, 0.050];       %Origin
p20 = [-0.01224, -0.00887, 0.02787];  %Insertion distance from theta1

resultData = load('Bifemsh_20mm_Result.mat', ...
    'xBest', 'routeCtx', 'Xi3');

requiredResultFields = {'xBest','routeCtx','Xi3'};
for iField = 1:numel(requiredResultFields)
    if ~isfield(resultData, requiredResultFields{iField})
        error('Knee_Flexor_Data_20mm:MissingResultField', ...
            ['Bifemsh_20mm_Result.mat does not contain %s. ' ...
             'Rerun the updated Opt_run.m.'], requiredResultFields{iField})
    end
end

xBest = resultData.xBest;
routeCtx = resultData.routeCtx;

% Rebuild the current flexor geometry so this script can be used for
% post-optimization hand tuning.  The saved optimizer result supplies the
% starting design and kinematic context, while buildGeoExclusion supplies
% the current hard geometry and current BPA/wrap definitions:
%   bpaRb = nominal BPA centerline distance from the hard surface
%   bpaRs = minimum allowable BPA centerline clearance from the hard surface
%   wRap  = obstacle wrapping radius used for the Xi3 calculation
routeCtx.geo = buildGeoExclusion();

Xi3 = resultData.Xi3;

p1 = xBest(1:3); %Origin
p2 = xBest(4:6); %End/insertion in the t1 frame

% Use the saved kinematic/routing context with the current geometry
% definitions and the Xi3 value saved by Opt_run.
[Location, bendMeasure, routeInfo] = ...
    buildKneeFlexorRoute20mm(p1, p2, xBest(8), routeCtx);

% Original two-point route: no intermediate wrapping point and zero bend.
Location0 = zeros(2,3,positions);
for i = 1:positions
    p20ICR = RowVecTrans(T_ICR_t1(:,:,i), p20);
    Location0(:,:,i) = [p10; p20ICR];
end
bendMeasure0 = zeros(positions,1);

% BPA 2 (Ben, 2026-09-21 asymmetric routing): pEnd{2} keeps the mirrored
% distal attachment and p1{2} shares BPA 1's side of the knee
% (flexorBpa2Endpoints20mm).  Every optimized curve below is the PAIR
% through MonoPam_mult, matching predictKneeFlexor20mm / Opt_run.
[p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2);
[LocationB, bendMeasureB] = ...
    buildKneeFlexorRoute20mm(p1B, p2B, xBest(8), routeCtx);
LocationPair = {Location; LocationB};
bendMeasurePair = {bendMeasure; bendMeasureB};

%20 mm Festo
Dia = 20;
% rest = 0.423; %resting length, m
% kmax = 0.322; %Length at maximum contraction, m
rest0 = 0.415; %resting length, m
kmax0 = (1-.255)*rest0; %Length at maximum contraction, m
tendon0 = 0.015; 
% fitting = 0.021; %Lower profile fittings at this BPA diameter

%from optimization:
rest   =  xBest(7);
tendon = xBest(8);
KMAX = 0.255; %maximum contraction percentage
kmax = (1 - KMAX)*rest;  % (1 - KMAX)*rest; Maximum contraction length.
fitting = 0.021;
%pres1 = 273.9783;         %average pressure, first test
pres1 = 0;
pres2 = 325;         %average pressure, first test
%pres3 = 606.4926;         %average pressure, first test
pres3 = 620;

% Load optimized stiffness parameters.
% Prefer the run's own Xi block (XiUsed = [Xi0 Xi1 Xi2 Xi3], saved by
% Opt_run alongside the displayed design) so the curves match the
% optimization that produced xBest.  Mats without XiUsed fall back to
% the legacy 2026-07-30 pick below.
try
    resultXi = load('Bifemsh_20mm_Result.mat', 'XiUsed');
    XiUsed = resultXi.XiUsed;
    Xi0 = XiUsed(1);
    Xi1 = XiUsed(2);
    Xi2 = XiUsed(3);
    Xi3 = XiUsed(4);
    fprintf(['Stiffness from the result mat XiUsed: ' ...
        'Xi0=%.6g, Xi1=%.6g, Xi2=%.6g, Xi3=%.6g\n'], Xi0, Xi1, Xi2, Xi3)
catch
    load minimizeFlxPin10_results_20260730_2transforms_Z2.mat filtered_results xCols
    pick = 1;
    g = filtered_results(pick,xCols);
    Xi0 = g(1);
    Xi1 = g(2);
    Xi2 = g(3);
    Xi3 = resultData.Xi3;
    fprintf(['Stiffness from legacy 20260730 pick (mat has no XiUsed): ' ...
        'Xi0=%.6g, Xi1=%.6g, Xi2=%.6g, Xi3=%.6g\n'], Xi0, Xi1, Xi2, Xi3)
end

wraps = 6;
BPAcount = 2;   % optimized curves = the BPA PAIR (BPA 2 = derived route)

% Original work: correct two-point Location and exactly zero X3 bend.
% Single-BPA comparator, same as predictOriginalKneeFlexor20mm.
Bifemsh_Pam0 = MonoPamDataExplicit_balanceX3(Name, Location0, CrossPoint, Dia, T_Pam, rest0, kmax0, tendon0, fitting, pres3, Xi0, Xi1, Xi2, Xi3, wraps, phiD, 1, bendMeasure0);

%Optimizer results: BOTH BPAs through MonoPam_mult, as in Opt_run
Bifemsh_Pam1 = MonoPam_mult(Name, LocationPair, CrossPoint, Dia, T_Pam, rest, kmax, tendon, fitting, pres1, Xi0, Xi1, Xi2, Xi3, wraps, phiD, BPAcount, bendMeasurePair);
Bifemsh_Pam2 = MonoPam_mult(Name, LocationPair, CrossPoint, Dia, T_Pam, rest, kmax, tendon, fitting, pres2, Xi0, Xi1, Xi2, Xi3, wraps, phiD, BPAcount, bendMeasurePair);
Bifemsh_Pam3 = MonoPam_mult(Name, LocationPair, CrossPoint, Dia, T_Pam, rest, kmax, tendon, fitting, pres3, Xi0, Xi1, Xi2, Xi3, wraps, phiD, BPAcount, bendMeasurePair);

%% Create strings for later plots
%First pressure
sT1 = sprintf('Theoretical %d kPa',pres1);
sM1 = sprintf('Measured %d kPa',pres1);
%Second pressure
sT2 = sprintf('Theoretical %d kPa',pres2);
sM2 = sprintf('Measured %d kPa',pres2);
%Third pressure
sT0 = sprintf('Original %d kPa',pres3);
sT3 = sprintf('Theoretical %d kPa',pres3);
sM3 = sprintf('Measured %d kPa',pres3);

%% Unstacking the Torques to identify specific rotations
Torque1 = Bifemsh.Torque;
TorqueR = Bifemsh_Pam3.Torque(:,:,1);
TorqueR_adj = Bifemsh_Pam3.Torque_p(:,:,1);

%% Add Torques from the Muscle Group
TorqueH = Torque1;

H = readmatrix('OpenSim_Bifem_Results.txt', ...
    'FileType', 'text', ...
    'NumHeaderLines', 7);

humanAngle = H(:,2);
TorqueHz = H(:,4);

%% Plotting Torque Results
phiD = phi*180/pi;

TorqueEx = zeros(size(TorqueH, 1), 1);
TorqueEy = zeros(size(TorqueH, 1), 1);
TorqueEz = zeros(size(TorqueH, 1), 1);

for i = 1:size(TorqueR, 1)
    if TorqueH(i, 1) >= 0
        TorqueEx(i) = TorqueR_adj(i, 1) - TorqueH(i, 1);
    else
        TorqueEx(i) = TorqueH(i, 1) - TorqueR_adj(i, 1);
    end
    
    if TorqueH(i, 2) >= 0
        TorqueEy(i) = TorqueR_adj(i, 2) - TorqueH(i, 2);
    else
        TorqueEy(i) = TorqueH(i, 2) - TorqueR_adj(i, 2);
    end
    
    if TorqueH(i, 3) >= 0
        TorqueEz(i) = TorqueR_adj(i, 3) - TorqueHz(i);
    else
        TorqueEz(i) = TorqueHz(i) - TorqueR_adj(i, 3);
    end
end

figure
hold on
sgtitle('Bicep Femoris Short Head Torque through Knee Flexion and Extension')

subplot(3, 2, 1)
plot(humanAngle, TorqueHz, phiD, TorqueR_adj(:, 3))
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
plot(phiD, TorqueH(:, 2), phiD, TorqueR_adj(:, 2))
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
plot(phiD, TorqueH(:, 1), phiD, TorqueR_adj(:, 1))
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
plot(phiD, Bifemsh_Pam3.Torque_p(:, 3), phiD, TorqueR(:, 3),phiD, Bifemsh_Pam0.Torque_p(:, 3),phiD, TorqueHz)
title('Muscle and PAM Z Torque')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Torque, N \cdot m','Interpreter','tex')
legend('Optimized stiffness aware', 'Stiffness unaware',"Unoptimized, stiffness aware","Human values")


% %% Plotting muscle lengths and moment arms using two different moment arm
% %calculations
% ML = Bifemsh.MuscleLength;
% PamL = Bifemsh_Pam3.L_p;
% for i = 1:size(Bifemsh.MomentArm,1)
%     MA(i,:) = norm(Bifemsh.MomentArm(i,1:2));               %Muscle moment arm, Z axis
%     BPAma(i,:) = norm(Bifemsh_Pam3.mA_p(i,1:2));        %BPA moment arm, Z axis
% end
% dM = diff(Bifemsh.MuscleLength);           %Muscle length difference
% dP = diff(Bifemsh_Pam3.L_p);       %PAM length difference
% dO = diff(phiD);                           %Angle difference
% 
% figure
% hold on
% sgtitle('Bicep Femoris Short Head Length and Moment Arm through Knee Flexion and Extension')
% 
% subplot(2, 2, 1)
% plot(phiD, ML, phiD, PamL)
% title('Muscle and PAM Lengths')
% xlabel('Knee angle, \circ','Interpreter','tex')
% ylabel('Length, m')
% legend('Human', 'PAM')
% 
% subplot(2, 2, 2)
% plot(phiD, MA, phiD, BPAma)
% title('Moment arm, Z axis, vector method')
% xlabel('Knee angle, \circ','Interpreter','tex')
% ylabel('Length, m')
% legend('Human', 'PAM')
% 
% subplot(2, 2, 3)
% plot(phiD(1:99), -dM./dO', phiD(1:99), -dP./dO')
% title('Moment arm, Z axis, left difference method')
% xlabel('Knee angle, \circ','Interpreter','tex')
% ylabel('Length, m')
% legend('Human', 'PAM')
% 
% subplot(2, 2, 4)
% plot(phiD(2:100), -dM./dO', phiD(2:100), -dP./dO')
% title('Moment arm, Z axis, right difference method')
% xlabel('Knee angle, \circ','Interpreter','tex')
% ylabel('Length, m')
% legend('Human', 'PAM')
% 
% hold off
%% Plotting the angle between the vectors

aHR = zeros(size(TorqueHz, 1), 1);
aHRH = zeros(size(TorqueHz, 1), 1);

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


%% Compare to results
%Longer Tibia
Load = [18.5 18.05 32 40.97 44.4 50.84 62.45 69.4 64.3 70.6 90.15 70];     %Load in Newtons
K_ang = [-125 -114 -98 -83 -75.5 -69 -55.5 -53.001 3 7 -6.5 -32]*c;      %Knee angle
LC_ang = [32.5 30 28 26 24.5 17 20 24 7 5.5 10 13.5]*c;      %Load Cell angle

d = 320/1000;
ang = -82.97;
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
    Fr(:,:,i) = [0; 0; 0; Load(i)*cos(LC_ang(i)+pi); Load(i)*sin(LC_ang(i)+pi); 0];
    AdTrk(:,:,i) = Adjoint(Trk(:,:,i));
    Fk(:,:,i) = AdTrk(:,:,i)'*Fr(:,:,i);
    
end

TorqueZ1 = Fk(3,1,:);
TorqueZ1 = squeeze(TorqueZ1);

%Tibia 2
Load2 = [91.7 114.1 129.8 78 92.78 63.6 97.28 71.92 84.5 93.9];     %Load in Newtons
K_ang2 = [-53.01 -41 -30 -26.01 -26.001 -18.5 -18 -7 -9 0]*c;      %Knee angle
LC_ang2 = [-5 -7 -12 -12 -12.5 -16.5 -15 -19.5 -15 -17]*c;      %Load Cell angle

d2 = 218.29/1000;
ang2 = -79.84;
p_rf2 = [d2*cosd(ang2), d2*sind(ang2), 0]';     %point of reaction force
T_t1_rf2 = RpToTrans(eye(3),p_rf2);   %Tranformation matrix from theta 1 to reaction point
Trk2 = pagemtimes(TransInv(T_t1_rf2),T_t1_ICR);
s3 = Trk2(1,4,:);
s3 = squeeze(s3);
fcn17 = fit(phi',s3,'cubicspline');
s4 = Trk2(2,4,:);
s4 = squeeze(s4);
fcn18 = fit(phi',s4,'cubicspline');

Trk2 = zeros(4,4,length(Load2));
Fr2 = zeros(6,1,length(Load2));
AdTrk2 = zeros(6,6,length(Load2));
Fk2 = zeros(6,1,length(Load2));

for i=1:length(Load2)
    Trk2(:,:,i) = RpToTrans(eye(3),[fcn17(K_ang2(i)), fcn18(K_ang2(i)), 0]');
    Fr2(:,:,i) = [0; 0; 0; Load2(i)*cos(LC_ang2(i)+pi); Load2(i)*sin(LC_ang2(i)+pi); 0];
    AdTrk2(:,:,i) = Adjoint(Trk2(:,:,i));
    Fk2(:,:,i) = AdTrk2(:,:,i)'*Fr2(:,:,i);
    
end

TorqueZ2 = Fk2(3,1,:);
TorqueZ2 = squeeze(TorqueZ2);

TorqueZ = [TorqueZ1(1:8); TorqueZ2; TorqueZ1(9:12)];
K_ang = [K_ang(1:8)'; K_ang2'; K_ang(9:12)'];

Presh = [613	614	614	615	615	615	618	618	620	620	620	385	451	296.6	421	281	324.8	325.58	325.58	325	500	560]; %Measured pressure

%% Create accessible color scheme
c1 = '#FFD700'; %gold
c2 = '#FFB14E'; %orange
c3 = '#FA8775'; %light orange
c4 = '#EA5F94'; %pink
c5 = '#CD34B5'; %magenta
c6 = '#9D02D7'; %magenta 2
c7 = '#0000FF'; %indigo
c8 = '#000000'; %black
sz = 60*Presh/620;        %size of data points
sz2 = sz*0.666; %size of second data points
C = {c1; c2; c3; c4; c5; c6; c7; c8};

%% Plot the results
Cang = K_ang/c;

figure
hold on
plot(phiD, Bifemsh_Pam0.Torque_p(:, 3),'Color',c2)
scatter(K_ang/c, TorqueZ,sz,'filled','MarkerFaceColor',c4)
plot(humanAngle, TorqueHz,'--','Color',c8);
plot(phiD, Bifemsh_Pam3.Torque_p(:,3),'Color',c7)
legend(sT3,sM3,"Human","New theoretical")
title('PAM Z Torque')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Torque, N \cdot m','Interpreter','tex')
ax = gca;
set(ax,'FontWeight','bold','LineWidth',2,'FontSize',10)
kid = ax.Children;
set(kid,'LineWidth',2)
hold off

%% Compare theoretical to OpenSim
TabMA = readmatrix('OpenSim_Bifem_MomentArm.txt');
knee_angle_rMA = TabMA(:,2)';           %Angle values directly from O
Bifemsh_MA = TabMA(:,3)';              %Torque values directly from OpenSim

Tab = readmatrix('OpenSim_Bifem_Results.txt');
knee_angle_rT = Tab(:,2)';           %Angle values directly from O
Bifemsh_T = Tab(:,4)';              %Torque values directly from OpenSim

figure
hold on
%plot(phiD, Bifemsh_Pam3.Torque(:,3),'-b', phiD, Bifemsh_Pam_adj2.Torque(:,3),'--r',phiD, Bifemsh_Pam_adj1.Torque(:,3),'.-g', K_ang/c, TorqueZ,'o', knee_angle_rT, Bifemsh_T,':k','LineWidth',2)
plot(phiD, Bifemsh_Pam3.Torque(:,3),':','Color',c7)
plot(phiD, Bifemsh_Pam3.Torque_p(:,3),'LineStyle','--','Color',c7); 
plot(phiD, Bifemsh_Pam2.Torque_p(:,3),'LineStyle','--','Color',c6);
plot(phiD, Bifemsh_Pam1.Torque_p(:,3),'LineStyle','--','Color',c5);
plot(phiD, Bifemsh_Pam0.Torque_p(:,3),'LineStyle','-.','Color',c3);
plot(knee_angle_rT, Bifemsh_T,'LineStyle',':','Color',c2);
scatter(K_ang/c, TorqueZ,sz,'filled','MarkerFaceColor',c4')
legend('Unoptimized',sT3,sT2,sT1,sT0,'OpenSim Human Torque','Measured','Location','southwest')
title('Knee Torque, Human vs. $\phi$20 mm BPA, $l_{rest}=420$ mm','Interpreter','latex')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('Torque, N \cdot m','Interpreter','tex')
ax = gca;
set(ax,'FontWeight','bold','LineWidth',2,'FontSize',10)
kid = ax.Children;
set(kid,'LineWidth',2)
hold off

%% Plot length and compare
bpa = Bifemsh_Pam3;
bpa0 = Bifemsh_Pam0;

Lm = bpa.RestingL .* (1 - bpa.strain_p{1});
Lm = Lm(:);

Lm0 = bpa0.RestingL .* (1 - bpa0.strain_p);
Lm0 = Lm0(:);

Angle = [-125 -114 -98 -83 -75.5 -69 -55.5 -53.001 -53.01 -41 -30];
InflatedLength = [334 334 338 343 NaN 348 356 356 357 365 368]/1000;

figure
hold on
scatter(Angle, InflatedLength)
plot(phiD(:), Lm)
plot(phiD(:), Lm0)
hold off
legend('Measured','Prediction, adjusted','Prediction, original')
title('BPA length')
xlabel('Knee angle, degrees')
ylabel('Length, m')

%% X3 strain definitions
figure
hold on
plot(phiD(:), Bifemsh_Pam3.strain_f{1}, 'LineWidth', 2)
plot(phiD(:), Bifemsh_Pam3.strain_p{1}, '--', 'LineWidth', 2)
plot(phiD(:), Bifemsh_Pam3.Contraction{1}, '-.', 'LineWidth', 2)
yline(0, ':k', 'Minimum strain')
yline(KMAX, ':', 'KMAX')
hold off
box off
grid off
legend('strain_f, includes Xi3', 'strain_p, excludes Xi3', ...
    'Contraction', 'Location', 'best')
xlabel('Knee angle, degrees')
ylabel('Strain')
title('BPA strain definitions')


%% Muscle Bone Plotting

% Static + animated muscle/bone plots use the SAME route-aware interface
% as Knee_Extensor_20mm.m: one bonePlotArgs cell list, MuscleBonePlotting
% for the static pose, AnimateKneeBoneMuscle for the animation. Location
% rows before CrossPoint are femur-frame; rows CrossPoint:end are already
% in ICR coordinates (exactly what buildKneeFlexorRoute20mm returns and
% what MonoPamDataExplicit_balanceX3 consumes).

% plotIdxUse = pos (phi was forced to 0 there). One pose per route list.
plotIdxUse = pos;

% Robot BPA-1 route for the bone plots is the saved Location itself: p1 is
% in the femur frame; the intermediate wrap and p2 rows are in ICR
% coordinates. Do NOT use Bifemsh_Pam3.Location here, because that stores
% v2 after transforming p2 into the knee/ICR frame for the torque
% calculation.

% Robot BPA-2 route (Ben, 2026-09-21): pEnd{2} keeps the mirrored distal
% attachment, but p1{2} shares BPA 1's side of the knee:
%   p1{2}z = p1{1}z - (pEnd{1}z - pEnd{2}z)
% (flexorBpa2Endpoints20mm). The route is solved, not mirrored, so its
% wrap point follows this asymmetric path.
[p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2);
[LocationB, ~, routeInfoB] = ...
    buildKneeFlexorRoute20mm(p1B, p2B, tendon, routeCtx);

fprintf(['BPA 2 route: p1{2} = [%.6f %.6f %.6f] m (femur), ' ...
    'pEnd{2} = [%.6f %.6f %.6f] m (t1); wrap active at home pose = %d\n'], ...
    p1B, p2B, routeInfoB.active(plotIdxUse))

bonePlotArgs = {T, T_ICR_t1, phi, pos, p1, p2, Bifemsh, ...
    'DisplayRotation',eye(3), ...
    'DisplayAxisMap',eye(3), ...
    'Location',Location,'CrossPoint',CrossPoint,'T_Pam',T_Pam, ...
    'Location2',LocationB, ...
    'HumanLabels',{'Biceps Femoris (Short Head)'}, ...
    'XLim',[], 'YLim',[], 'ZLim',[]};

%Static display
run('MuscleBonePlotting.m')

%Animated continous loop
% AnimateKneeBoneMuscle(bonePlotArgs{:}, ...
%     'PauseTime', 0.18, ...
%     'FrameStep', 1, ...
%     'Loop', true);

%Export as GIF (zero -> full flexion -> full extension -> zero, same
%frame order as Knee_Extensor_20mm.m; empty limits fit both BPA routes)
AnimateKneeBoneMuscle(bonePlotArgs{:}, ...
    'FullSkeleton',true, ... % false = femur and tibia only
    'FrameIndices',[pos:-1:1, 2:positions, positions-1:-1:pos], ...
    'PauseTime',0.02, ...
    'Loop',false, ...
    'ExportGif',true, ...
    'GifFile','Knee_Flexor_20mm.gif', ...
    'FrameRate',20, ...
    'CameraOrbitDeg', -90, ...
    'XLim',[],'YLim',[],'ZLim',[])
clear bonePlotArgs

%Export as MP4
% AnimateKneeBoneMuscle(bonePlotArgs{:}, ...
%     'PauseTime', 0.02, ...
%     'Loop', false, ...
%     'ExportVideo', true, ...
%     'VideoFile', 'Knee_Flexor_20mm.mp4', ...
%     'FrameRate', 20)

%To pause at a specific frame (for example frame 50)
% AnimateKneeBoneMuscle(bonePlotArgs{:}, ...
%     'PauseAtFrames', 50)
