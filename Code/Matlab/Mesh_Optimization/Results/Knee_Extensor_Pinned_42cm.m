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
% R_Pam = zeros(3, 3, positions);
% T_Pam = zeros(4, 4, positions);

c = pi/180; %Convert from degrees to radians

kneeMin = -120*c;
kneeMax = 35*c;
phi = linspace(kneeMin, kneeMax, positions);

%We want one of our positions to be home position, so let's make the
%smallest value of phi equal to 0
[val, pos] = min(abs(phi));
phi(pos) = 0;
phiD = phi*180/pi;       %Knee angle in degrees

for i = 1:positions
    hipToKnee = [0.0045, -0.3958, 0];              %pinned location
    R(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    
    T(:, :, i) = RpToTrans(R(:, :, i), hipToKnee');
    
end

%% PAM calculation
Name = 'Vastus Intermedius Normal';
Location = zeros(5,3,positions);

% Origin Location from Ben
for i = 1:positions
    if phiD(i) >= -74.01 && phiD(i) < -19.6
     Location(:,:,i) = [0.030, -0.050, 0;     %Origin
                0.060, -0.350, 0.000;         %BPA contacts screw that joins femur body with femoral condyles
                0.04128, -0.410,    0;        %Contact point between 19.6 and 74.01 degrees flexion
                0.04128, -0.410,    0;        %Row 4 = Row 3
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    elseif phiD(i) >= -19.6 && phiD(i) < 20
     Location(:,:,i) = [0.030, -0.050, 0;     %Origin
                0.060, -0.350, 0.000;         %BPA contacts screw that joins femur body with femoral condyles
                0.060, -0.350, 0.000;         %Row 3 = Row 2
                0.060, -0.350, 0.000;         %Row 4 = Row 2
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    elseif phiD(i) >= 20
     Location(:,:,i) = [0.030, -0.050, 0;     %Origin
                0.030, -0.050, 0;             %Row 2 = Row 1 (no screw contact)
                0.030, -0.050, 0;             %Row 3 = Row 2
                0.030, -0.050, 0;             %Row 4 = Row 2
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    else
     Location(:,:,i) = [0.030, -0.050, 0;     %Origin
                0.060, -0.350, 0.000;         %BPA contacts screw that joins femur body with femoral condyles
                0.04128, -0.410,    0;        %Contact point between 19.6 and 74.01 degrees flexion
                0.01138, -0.425 0;            %Contact point over 74.01 degrees flexion
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    end
end
CrossPoint = 5;
Dia = 10;
rest = (415)/1000;   %resting length clamp to clamp, minus the barb
kmax = 0.349;          %length at maximum contraction
tendon = 0;             %Tendon length
fitting = 0.0254;           %fitting length
pres = 605;             %Pressure, kPa
Vas_Pam_ideal = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);

Name = 'Vastus Intermedius Realistic, shorter resting length';
Location = zeros(5,3,positions);
for i = 1:positions
    if phiD(i) >= -74.01 && phiD(i) < -19.6
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin
                0.060, -0.350, 0.000;         %BPA contacts screw that joins femur body with femoral condyles
                0.04128, -0.410,    0;        %Contact point between 19.6 and 74.01 degrees flexion
                0.04128, -0.410,    0;        %Row 4 = Row 3
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    elseif phiD(i) >= -19.6 && phiD(i) < 20
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin
                0.060, -0.350, 0.000;        %BPA contacts screw that joins femur body with femoral condyles
                0.060, -0.350, 0.000;        %Row 3 = Row 2
                0.060, -0.350, 0.000;        %Row 4 = Row 2
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    elseif phiD(i) >= 20
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin
                0.030, -0.050, 0;        %Row 2 = Row 1 (no screw contact)
                0.030, -0.050, 0;        %Row 3 = Row 2
                0.030, -0.050, 0;        %Row 4 = Row 2
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    else
     Location(:,:,i) = [0.025, -0.050, 0;             %Origin, origin point goes -5mm at high angles due to load
                0.060, -0.350, 0.000;         %BPA contacts screw that joins femur body with femoral condyles
                0.04128, -0.410,    0;        %Contact point between 19.6 and 74.01 degrees flexion
                0.01568, -0.42093 0;            %Contact point over 74.01 degrees flexion, added compression of BPA
                0.0375, -0.074, 0.000];     %Tibia bracket (insertion)
    end
end      
CrossPoint = 5;
Dia = 10;
rest = (415-8)/1000;   %resting length clamp to clamp, minus the barb
kmax = 0.349;          %length at maximum contraction
tendon = 0;             %Tendon length
fitting = 0.0254;           %fitting length
pres = 605;             %Pressure, kPa
Vas_Pam_real1 = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);

Name = 'Vastus Intermedius Realistic 2, shorter fittings';
rest = 415/1000;
fitting = 0.0254-0.004;           %fitting length
Vas_Pam_real2 = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);

Name = 'Vastus Intermedius Realistic 3, longer resting length';
rest = (415+8)/1000;
fitting = 0.0254;           %fitting length
Vas_Pam_real3 = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);

Name = 'Vastus Intermedius Tendon and slip';
Location = zeros(5,3,positions);
for i = 1:positions
    if phiD(i) >= -74.01 && phiD(i) < -19.6
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin
                0.053, -0.350, 0.020;         %Contact, no screw
                0.04128, -0.410,    0.015;        %Contact point between 19.6 and 74.01 degrees flexion
                0.04128, -0.410,    0.015;        %Row 4 = Row 3
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    elseif phiD(i) >= -19.6 && phiD(i) < 20
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin
                0.064, -0.350, 0;        %Contact, screw
                0.064, -0.350, 0;        %Row 3 = Row 2
                0.064, -0.350, 0;        %Row 4 = Row 2
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    elseif phiD(i) >= 20
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin
                0.030, -0.050, 0;        %Row 2 = Row 1 (no screw contact)
                0.030, -0.050, 0;        %Row 3 = Row 2
                0.030, -0.050, 0;        %Row 4 = Row 2
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    else
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin, origin point goes -5mm at high angles due to load
                0.053, -0.350, 0.020;         %Contact, no screw
                0.04128, -0.410,    0.015;        %Contact point between 19.6 and 74.01 degrees flexion
                0.01568, -0.42093 0.010;            %Contact point over 74.01 degrees flexion, added compression
                0.0375, -0.074, 0.000];     %Tibia bracket (insertion)
    end
end      
CrossPoint = 5;
Dia = 10;
rest = (415)/1000;   %resting length clamp to clamp, minus the barb
kmax = 0.349;          %length at maximum contraction
tendon = 0.02;             %Tendon length
fitting = 0.0254;           %fitting length
pres = 605;             %Pressure, kPa
Vas_Pam_slip = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);

Name = 'Vastus Intermedius Tendon, real, no slip';
Location = zeros(5,3,positions);
for i = 1:positions
    if phiD(i) >= -74.01 && phiD(i) < -19.6
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin
                0.060, -0.350, 0.000;         %BPA contacts screw that joins femur body with femoral condyles
                0.04128, -0.410,    0;        %Contact point between 19.6 and 74.01 degrees flexion
                0.04128, -0.410,    0;        %Row 4 = Row 3
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    elseif phiD(i) >= -19.6 && phiD(i) < 20
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin
                0.060, -0.350, 0.000;        %BPA contacts screw that joins femur body with femoral condyles
                0.060, -0.350, 0.000;        %Row 3 = Row 2
                0.060, -0.350, 0.000;        %Row 4 = Row 2
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    elseif phiD(i) >= 20
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin
                0.030, -0.050, 0;        %Row 2 = Row 1 (no screw contact)
                0.030, -0.050, 0;        %Row 3 = Row 2
                0.030, -0.050, 0;        %Row 4 = Row 2
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    else
     Location(:,:,i) = [0.025, -0.050, 0;             %Origin, origin point goes -5mm at high angles due to load
                0.060, -0.350, 0.000;         %BPA contacts screw that joins femur body with femoral condyles
                0.04128, -0.410,    0;        %Contact point between 19.6 and 74.01 degrees flexion
                0.01568, -0.42093 0;            %Contact point over 74.01 degrees flexion observed compression
                0.0375, -0.074, 0.000];     %Tibia bracket (insertion)
    end
end      
CrossPoint = 5;
Dia = 10;
rest = (415)/1000;   %resting length clamp to clamp, minus the barb
kmax = 0.349;          %length at maximum contraction
tendon = 0.022;             %Tendon length
fitting = 0.0254;           %fitting length
pres = 605;             %Pressure, kPa
Vas_Pam_tendon_real = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);

Name = 'Vastus Intermedius Tendon, ideal';
Location = zeros(5,3,positions);
for i = 1:positions
    if phiD(i) >= -74.01 && phiD(i) < -19.6
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin
                0.060, -0.350, 0.000;         %BPA contacts screw that joins femur body with femoral condyles
                0.04128, -0.410,    0;        %Contact point between 19.6 and 74.01 degrees flexion
                0.04128, -0.410,    0;        %Row 4 = Row 3
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    elseif phiD(i) >= -19.6 && phiD(i) < 20
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin
                0.060, -0.350, 0.000;        %BPA contacts screw that joins femur body with femoral condyles
                0.060, -0.350, 0.000;        %Row 3 = Row 2
                0.060, -0.350, 0.000;        %Row 4 = Row 2
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    elseif phiD(i) >= 20
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin
                0.030, -0.050, 0;        %Row 2 = Row 1 (no screw contact)
                0.030, -0.050, 0;        %Row 3 = Row 2
                0.030, -0.050, 0;        %Row 4 = Row 2
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    else
     Location(:,:,i) = [0.030, -0.050, 0;             %Origin, origin point goes -5mm at high angles due to load
                0.060, -0.350, 0.000;         %BPA contacts screw that joins femur body with femoral condyles
                0.04128, -0.410,    0;        %Contact point between 19.6 and 74.01 degrees flexion
                0.01138, -0.425 0;            %Contact point over 74.01 degrees flexion observed compression
                0.0425, -0.07591, 0.000];     %Tibia bracket (insertion)
    end
end      
CrossPoint = 5;
Dia = 10;
Rest = (415)/1000;   %resting length clamp to clamp, minus the barb
kmax = 0.349;          %length at maximum contraction
tendon = 0.022;             %Tendon length
fitting = 0.0254;           %fitting length
pres = 605;             %Pressure, kPa
Vas_Pam_tendon_ideal = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);

%% Unstacking the Torques to identify specific rotations
% Force1 = Vas_Int.Force + Vas_Lat.Force + Vas_Med.Force;
% Torque1 = Vas_Int.Torque + Vas_Lat.Torque + Vas_Med.Torque;
TorqueR1 = Vas_Pam_ideal.Torque;
TorqueR2 = Vas_Pam_real1.Torque;
TorqueR3 = Vas_Pam_real2.Torque; 
TorqueR4 = Vas_Pam_real3.Torque;
TorqueR5 = Vas_Pam_tendon_ideal.Torque;
TorqueR6 = Vas_Pam_tendon_real.Torque;
TorqueR7 = Vas_Pam_slip.Torque;

%% Add Torques from the Muscle Group
% TorqueH = Torque1;


%% Plotting Torque Results

figure
plot(phiD, TorqueR1(:, 3),'-',phiD, TorqueR2(:, 3),'--^',phiD, TorqueR3(:, 3),'--<',phiD, TorqueR4(:, 3),'--v',phiD, TorqueR5(:, 3),'-',phiD, TorqueR6(:, 3),'-.',phiD, TorqueR7(:,3),':')
title('BPA Z Torque, Length = 415 mm')
xlabel('Knee Extension(+)/Flexion(-), degrees')
ylabel('Torque, Nm')
legend('Ideal','Realistic1','Realistic2','Realistic3','Tendon Ideal','Tendon Real','Slippage+Tendon')

figure
plot(phiD, TorqueR5(:, 1),phiD, TorqueR7(:, 1),'-.')
title('BPA X Torque, Length = 415 mm')
xlabel('Knee Extension(+)/Flexion(-), degrees')
ylabel('Torque, Nm')
legend('Ideal Tendon','Slippage+Tendon')
% TorqueEx = zeros(size(TorqueH, 1), 1);
% TorqueEy = zeros(size(TorqueH, 1), 1);
% TorqueEz = zeros(size(TorqueH, 1), 1);
% 
% for i = 1:size(TorqueR, 1)
%     if TorqueH(i, 1) >= 0
%         TorqueEx(i) = TorqueR(i, 1) - TorqueH(i, 1);
%     else
%         TorqueEx(i) = TorqueH(i, 1) - TorqueR(i, 1);
%     end
%     
%     if TorqueH(i, 2) >= 0
%         TorqueEy(i) = TorqueR(i, 2) - TorqueH(i, 2);
%     else
%         TorqueEy(i) = TorqueH(i, 2) - TorqueR(i, 2);
%     end
%     
%     if TorqueH(i, 3) >= 0
%         TorqueEz(i) = TorqueR(i, 3) - TorqueH(i, 3);
%     else
%         TorqueEz(i) = TorqueH(i, 3) - TorqueR(i, 3);
%     end
% end

% figure
% hold on
% sgtitle('Bicep Femoris Short Head Torque through Knee Flexion and Extension')
% 
% subplot(3, 2, 1)
% plot(phiD, TorqueH(:, 3), phiD, TorqueR(:, 3))
% legend('Human Muscle', 'Optimal BPA Location')
% title('Muscle and PAM Z Torque')
% xlabel('Knee Extension/Rotation, degrees')
% ylabel('Torque, Nm')
% legend('Human', 'PAM')
% 
% subplot(3, 2, 2)
% plot(phiD, TorqueEz)
% legend('Optimal PAM Location')
% xlabel('Knee Extension/Rotation, degrees')
% ylabel('Torque, Nm')
% title('Adjusted Error Z Torque')
% 
% subplot(3, 2, 3)
% plot(phiD, TorqueH(:, 2), phiD, TorqueR(:, 2))
% title('Muscle and PAM Y Torque')
% xlabel('Knee Extension/Rotation, degrees')
% ylabel('Torque, Nm')
% legend('Human', 'PAM')
% 
% subplot(3, 2, 4)
% plot(phiD, TorqueEy)
% legend('Optimal PAM Location')
% xlabel('Knee Extension/Rotation, degrees')
% ylabel('Torque, Nm')
% title('Adjusted Error Y Torque')
% 
% subplot(3, 2, 5)
% plot(phiD, TorqueH(:, 1), phiD, TorqueR(:, 1))
% title('Muscle and PAM X Torque')
% xlabel('Knee Extension/Rotation, degrees')
% ylabel('Torque, Nm')
% legend('Human', 'PAM')
% 
% subplot(3, 2, 6)
% plot(phiD, TorqueEx)
% legend('Optimal PAM Location')
% xlabel('Knee Extension/Rotation, degrees')
% ylabel('Torque, Nm')
% title('Adjusted Error X Torque')
% 
% hold off