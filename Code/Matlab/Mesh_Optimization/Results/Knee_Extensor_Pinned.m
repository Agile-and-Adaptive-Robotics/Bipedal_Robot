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
% R_Pam = zeros(3, 3, positions);
% T_Pam = zeros(4, 4, positions);

c = pi/180; %Convert from degrees to radians

kneeMin = -120*c;
kneeMax = 10*c;
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
Name = 'Vastus Intermedius, Robot';
Location = zeros(5,3,positions);

% Origin Location from Ben
p1 = [0.030, -0.050, 0];             %Extensor Origin
p2 = [0.055, -0.350, 0.000];         %BPA contacts screw that joins femur body with femoral condyles
p2p = [0.053, -0.350, -0.005];         %BPA slips off screw that joins femur body with femoral condyles
p3 = [0.04128, -0.410,    -0.01];        %Contact point between 19.6 and 74.01 degrees flexion
p4 = [0.01138, -0.425, -0.010];            %Contact point over 74.01 degrees flexion
p5 = [0.0425, -0.07591, 0.000];      %%Tibia bracket (insertion)

for i = 1:positions
    if phiD(i) >= -74.01 && phiD(i) < -15.5
     Location(:,:,i) = [p1;...     %Origin
                        p2p;...    %BPA contacts screw that joins femur body with femoral condyles
                        p3;...     %Contact point between 19.6 and 74.01 degrees flexion
                        p3;...     %Row 4 = Row 3
                        p5];       %Tibia bracket (insertion)
    elseif phiD(i) >= -15.5 && phiD(i) < 20
     Location(:,:,i) = [p1;...     %Origin
                        p2;...     %BPA contacts screw that joins femur body with femoral condyles
                        p2;...     %Row 3 = Row 2
                        p2;...     %Row 4 = Row 2
                        p5];       %Tibia bracket (insertion)
    elseif phiD(i) >= 20
     Location(:,:,i) = [p1;...     %Origin
                        p1;...     %Row 2 = Row 1 (no screw contact)
                        p1;...     %Row 3 = Row 2
                        p1;...     %Row 4 = Row 2
                        p5];       %Tibia bracket (insertion)
    else
     Location(:,:,i) = [p1;...     %Origin
                        p2p;...         %BPA slips off screw that joins femur body with femoral condyles
                        p3;...        %Contact point between 19.6 and 74.01 degrees flexion
                        p4;...            %Contact point over 74.01 degrees flexion
                        p5];     %Tibia bracket (insertion)
    end
end
        
CrossPoint = 5;
Dia = 10;
fitting = 0.0254;           %fitting length
% tendon0 = 0;
% tendon22 = 0.022;

%41.5 cm, no tendon
rest = 415/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.350;           %length at maximum contraction
tendon = 0;             %Tendon length
pres = 605.2351;        %Pressure, kPa
Vas_Pam_42cm = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);

%41.5 cm,  tendon
rest = 415/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.350;           %length at maximum contraction
tendon = 0.0384;         %Tendon length
pres = 605.2351;        %Pressure, kPa
Vas_Pam_42cm_tendon = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);

%45.7 cm, no tendon
rest = 457/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.380;           %length at maximum contraction
tendon = 0;             %Tendon length
pres = 602;             %Pressure, kPa
Vas_Pam_46cm = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);

%48.6 cm, no tendon
rest = 485/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.405;           %length at maximum contraction
tendon = 0;             %Tendon length
pres = 605.8523;        %Pressure, kPa
Vas_Pam_48cm = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);

%% Unstacking the Torques to identify specific rotations
Torque_42cm = Vas_Pam_42cm.Torque(:,3,4);
Torque_42cm_ten = Vas_Pam_42cm_tendon.Torque(:,3,4);
Torque_46cm = Vas_Pam_46cm.Torque(:,3,4);
Torque_48cm = Vas_Pam_48cm.Torque(:,3,4);

%% Plotting Torque Results

figure
plot(phiD, Torque_42cm)
title('\bf Iso. Torque vs $\theta_{k}$, $l_{rest}=415\,$ mm, no tendon','Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque, $N \cdot m$','Interpreter','latex')

figure
plot(phiD, Torque_42cm_ten)
title('\bf Iso. Torque vs $\theta_{k}$, $l_{rest}=415\,$ mm, $38\,mm$ tendon','Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque, $N \cdot m$','Interpreter','latex')

figure
plot(phiD, Torque_46cm)
title('\bf Iso. Torque vs $\theta_{k}$, $l_{rest}=457\,$ mm','Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque, $N \cdot m$','Interpreter','latex')

figure
plot(phiD, Torque_48cm)
title('\bf Iso. Torque vs $\theta_{k}$, $l_{rest}=485\,$ mm','Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque, $N \cdot m$','Interpreter','latex')