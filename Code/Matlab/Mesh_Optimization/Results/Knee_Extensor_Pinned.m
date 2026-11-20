%% Mesh Points Calculation - Bicep Femoris Short Head
% This code will calculate the torque difference between all of the points
% from one bone mesh to another to determine the best location for muscle
% placement

%% Freshen up the workspace
clc
clear
close all

%% Add paths to the muscle and pam calculators
% current_dir = cd;
% all_code = fullfile(current_dir,'../..');
% addpath(genpath(all_code));

%% Joint rotation transformation matrices
positions = 100;
fprintf('The algorithm will be calculating Torque at %d different joint positions.\n', positions)

R = zeros(3, 3, positions);
T = zeros(4, 4, positions);
% R_Pam = zeros(3, 3, positions);
% T_Pam = zeros(4, 4, positions);

c = pi/180; %Convert from degrees to radians

kneeMin = -125*c;
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
Name = 'Vastus Intermedius, Robot';
Location = zeros(7,3,positions);
CrossPoint = 7;
Dia = 10;
fitting = 0.0254;           %fitting length

% Locations from Assem2.75
% p1 = [0.030, -0.050, 0];             %Extensor Origin
% p2 = [0.055, -0.350, 0.000];         %BPA contacts screw that joins femur body with femoral condyles
% p2p = [0.053, -0.350, -0.005];         %BPA slips off screw that joins femur body with femoral condyles
% p3 = [0.04128, -0.410,    -0.01];        %Contact point between 19.6 and 74.01 degrees flexion
% p4 = [0.01138, -0.425, -0.010];            %Contact point over 74.01 degrees flexion
% p5 = [0.0425, -0.07591, 0.000];      %%Tibia bracket (insertion)

p1 = [0.030, -0.050, 0];             %Extensor Origin
p2 = [0.065, -0.350, 0];         %BPA contacts screw that joins femur body with femoral condyles
p2p = [0.054, -0.350, -0.02];         %BPA slips off screw that joins femur body with femoral condyles
p3 = [0.04244, -0.40847,  0];        %Contact point 
p4 = [0.02818, -0.42157, 0];            %Contact point
p5 = [0.01052, -0.42519, 0];
p7 = [0.0425, -0.07591, 0.000];      %%Tibia bracket (insertion)

%Angles where via points added
ang1 = 27;
ang2 = -23;
ang3 = -53;
ang4 = -87.5;


for i = 1:positions
    if phiD(i) >= ang1
     Location(:,:,i) = [p1;...     %Origin
                        p1;...     %Row 2 (no screw contact)
                        p1;...     %Row 3
                        p1;...     %Row 4
                        p1;...     %Row 5
                        p1;...      %Row 6
                        p7];       %Tibia bracket (insertion)
    elseif phiD(i) >= ang2 && phiD(i) < ang1
     Location(:,:,i) = [p1;...     %Origin
                        p2;...        %Contact with screw that joins femur body with femoral condyles
                        p2;...        %Contact point
                        p2;...
                        p2;...
                        p2;...
                        p7];     %Tibia bracket (insertion)
    elseif phiD(i) >= ang3 && phiD(i) < ang2
     Location(:,:,i) = [p1;...     %Origin
                        p2;...         %BPA slips off screw that joins femur body with femoral condyles
                        p3;...        %Contact point
                        p3;...
                        p3;...
                        p3;...
                        p7];     %Tibia bracket (insertion)
    elseif phiD(i) >= ang4 && phiD(i) < ang3
     Location(:,:,i) = [p1;...     %Origin
                        p2;...         %BPA slips off screw that joins femur body with femoral condyles
                        p3;...        %Contact point
                        p4;...
                        p4;...
                        p4;...
                        p7];     %Tibia bracket (insertion)
    else
     Location(:,:,i) = [p1;...     %Origin
                        p2;...         %BPA slips off screw that joins femur body with femoral condyles
                        p3;...        %Contact point
                        p4;...
                        p5;...
                        hipToKnee(1)+0.03*cosd(phiD(i)+9.07), hipToKnee(2)+0.03*sind(phiD(i)+9.072), 0;...            %Contact point with -87.5 degrees flexion or more (i.e. more negative), when the angle is 87.5 the angle from horizontal to point 5 is 78.428 deg, so this 9.072 deg is the difference and point 6 follows this minimum circle.
                        p7];    %Tibia bracket (insertion)
                        % p7-[25.4/4*10^-3, 0, 0]];     %Tibia bracket
                        % (insertion) don't really know what this is doing
    end
end
        
%40.5 cm, no tendon
rest40 = 400/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.337;           %length at maximum contraction
tendon = 0;             %Tendon length
pres = 625;        %Pressure, kPa
Vas_Pam_40cm = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest40, kmax, tendon, fitting, pres);

%40.5 cm, 40 mm tendon
rest40 = 400/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.337;           %length at maximum contraction
tendon1 = 0.040;         %Tendon length
pres = 618;        %Pressure, kPa
Vas_Pam_40cm_tendon = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest40, kmax, tendon1, fitting, pres);

%41.5 cm, no tendon
rest42 = 415/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.3445;           %length at maximum contraction
tendon = 0;             %Tendon length
pres = 605.2351;        %Pressure, kPa
Vas_Pam_42cm = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest42, kmax, tendon, fitting, pres);

%41.5 cm,  tendon
rest42 = 415/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.3445;           %length at maximum contraction
tendon2 = 0.040;         %Tendon length
pres = 605.2351;        %Pressure, kPa
Vas_Pam_42cm_tendon = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest42, kmax, tendon2, fitting, pres);

%43.6 cm, no tendon
rest436 = 436/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.371;           %length at maximum contraction
tendon = 0;             %Tendon length
pres = 621.8;        %Pressure, kPa
Vas_Pam_43cm = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest436, kmax, tendon, fitting, pres);

%43.2 cm, 22 mm tendon
rest432 = 432/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.367;           %length at maximum contraction
tendon3 = 0.022;             %Tendon length
pres = 623;        %Pressure, kPa
Vas_Pam_43cm_tendon = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest432, kmax, tendon3, fitting, pres);

%45.7 cm, no tendon
rest46 = 457/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.3815;           %length at maximum contraction
tendon = 0;             %Tendon length
pres = 602;             %Pressure, kPa
Vas_Pam_46cm = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest46, kmax, tendon, fitting, pres);

%46.5 cm, no tendon
rest47 = 465/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.387;           %length at maximum contraction
tendon = 0;             %Tendon length
pres = 622.2;             %Pressure, kPa
Vas_Pam_47cm = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest46, kmax, tendon, fitting, pres);

%48.6 cm, no tendon
rest48 = 480/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.405;           %length at maximum contraction
tendon = 0;             %Tendon length
pres = 605.8523;        %Pressure, kPa
Vas_Pam_48cm = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest48, kmax, tendon, fitting, pres);

%% Unstacking the Torques to identify specific rotations
Torque_40cm = Vas_Pam_40cm.Torque(:,3);
Torque_40cm_ten = Vas_Pam_40cm_tendon.Torque(:,3);
Torque_42cm = Vas_Pam_42cm.Torque(:,3);
Torque_42cm_ten = Vas_Pam_42cm_tendon.Torque(:,3);
Torque_43cm = Vas_Pam_43cm.Torque(:,3);
Torque_43cm_ten = Vas_Pam_43cm_tendon.Torque(:,3);
Torque_46cm = Vas_Pam_46cm.Torque(:,3);
Torque_47cm = Vas_Pam_47cm.Torque(:,3);
Torque_48cm = Vas_Pam_48cm.Torque(:,3);

%% Plotting Torque Results
xLim = [-125 35];
yLim = [0 25];

figure
h = tiledlayout(5,2);

ax1 = nexttile(1);
plot(phiD, Torque_48cm)
title('\bf l_{rest}=48.5 cm','Interpreter','tex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque, N\cdotm','Interpreter','tex')
set(ax1,'XLim',xLim,'YLim',yLim)

ax2 = nexttile(2);
plot(phiD, Torque_47cm)
title('\bf l_{rest}=46.5 cm','Interpreter','tex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque, N\cdotm','Interpreter','tex')
set(ax2,'XLim',xLim,'YLim',yLim)

ax3 = nexttile(3);
plot(phiD, Torque_46cm)
title('\bf l_{rest}=45.7 cm','Interpreter','tex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque, N\cdotm','Interpreter','tex')
set(ax3,'XLim',xLim,'YLim',yLim)

ax5 = nexttile(5);
plot(phiD, Torque_43cm_ten)
title(sprintf('l_{rest} = 43.2 cm, tendon = %d mm',tendon3*10^3),'Interpreter','tex','FontWeight','bold')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque N\cdotm','Interpreter','tex')
set(ax5,'XLim',xLim,'YLim',yLim)

ax6 = nexttile(6);
plot(phiD, Torque_43cm)
title('\bf l_{rest}=43.6 cm, no tendon','Interpreter','tex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque, N\cdotm','Interpreter','tex')
set(ax6,'XLim',xLim,'YLim',yLim)

ax7 = nexttile(7);
plot(phiD, Torque_42cm_ten)
title(sprintf('l_{rest} = 41.5 cm, tendon = %d mm',tendon2*10^3),'Interpreter','tex','FontWeight','bold')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque N\cdotm','Interpreter','tex')
set(ax7,'XLim',xLim,'YLim',yLim)

ax8 = nexttile(8);
plot(phiD, Torque_42cm)
title('\bf l_{rest}=41.5 cm, no tendon','Interpreter','tex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque, N\cdotm','Interpreter','tex')
set(ax8,'XLim',xLim,'YLim',yLim)

ax9 = nexttile(9);
plot(phiD', Torque_40cm_ten)
title(sprintf('l_{rest} = 40.1 cm, tendon = %d mm',tendon1*10^3),'Interpreter','tex','FontWeight','bold')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque N\cdotm','Interpreter','tex')
set(ax9,'XLim',xLim,'YLim',yLim)

ax10 = nexttile(10);
plot(phiD, Torque_40cm)
title('\bf l_{rest}=40.5 cm, no tendon','Interpreter','tex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Torque, N\cdotm','Interpreter','tex')
set(ax10,'XLim',xLim,'YLim',yLim)