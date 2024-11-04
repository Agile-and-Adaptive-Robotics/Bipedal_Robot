%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

% f = fullfile('github/bipedal_robot/code/matlab');
% qt = addpath(genpath(f));

TabMA = readmatrix('OpenSim_Bifem_MomentArm.txt');
knee_angle_rMA = TabMA(:,2)';           %Angle values directly from O
Bifemsh_MA = TabMA(:,3)';              %Torque values directly from OpenSim

Tab = readmatrix('OpenSim_Bifem_Results.txt');
knee_angle_rT = Tab(:,2)';           %Angle values directly from O
Bifemsh_T = Tab(:,4)';              %Torque values directly from OpenSim

load KneeFlx_20mm_42cm.mat
% rest = 0.423; %resting length, m
% kmax = 0.322; %Length at maximum contraction, m
% rest = 0.420; %resting length, m
% kmax = 0.315; %Length at maximum contraction, m
% dia = 20;
% fitting = 0.0254;
KMAX = (rest-kmax)/rest;

Theo1 = Bifemsh_Pam1.Torque(:,3)';          %Original Torque calculations
Theo2 = Bifemsh_Pam2.Torque(:,3)';
Theo3 = Bifemsh_Pam3.Torque(:,3)';
Theo1opt = Bifemsh_Pam_adj1.Torque(:,3)';       %Optimized fitting length from previous study
Theo2opt = Bifemsh_Pam_adj2.Torque(:,3)';
Theo3opt = Bifemsh_Pam_adj3.Torque(:,3)';


%% Test 1 done with CALT load cell
%Test 1 == sheet FlxTest20mm_42cm (3) from Results_table10mm_FullSize
%Test 2 == sheet FlxTest20mm_42cm (2) from Results_table10mm_FullSize
%% Torque calculated from measurements
Angle = [-125	-114	-98	-83	-75.5	-69	-55.5	-53	-53	-41	-30	-26	-26	-18.5	-18	-7	-9	0	3	7.7	-6.5	-32];

Torque = [-8.84322291179594,-7.38771014573315,-6.12295907061857,-5.11416855627929,-4.00577506118469,-3.08271954812317,-3.04251451225018,-2.59384366567520,-1.86178531541841,-0.774849368079027,-4.65440909149842,-3.19703540259264,-2.57966382173969,-0.959098909344151,-0.499572417339300,-1.33299525437727,-0.772724305518311];
%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength = [334	334	338	343	NaN	348	356	356	357	365	368	381	377	387	385	396	392	397	398	401	384	367]/1000;

ICRtoMuscle = [25	25	35	35	NaN	37	40	45	48	45	40	40	40	42	45	46	48	48	42	43	45	38]/1000;

% %load pressure where applicable
% test = 1;
% runsperseries = 17;
% 
% pres = zeros(1,runsperseries);
%  
% 
%      for j = 1:runsperseries
%                 file_name = sprintf('FlxTest%0.0f_%0.0f.mat', test,j);
%                 load(file_name,'Stats')
%                 pres(1,j) = Stats{'Mean',2};
%      end
pres = [613	614	614	615	615	615	618	618	620	620	620	385	451	296.6	421	281	324.8	325.58	325.58	325	500	560];

contract = (rest-InflatedLength)/rest;
realRel = contract./KMAX;

F = festo4(20, realRel, pres);

TorqueHand = -ICRtoMuscle.*F;  %Moment will be negative because it causes flexion

%% Prepare for plots by creating colors and sizes
%Matlab hex color values:
% Create accessible color scheme
c1 = '#FFD700'; %gold
c2 = '#FFB14E'; %orange
c3 = '#FA8775'; %light orange
c4 = '#EA5F94'; %pink
c5 = '#CD34B5'; %magenta
c6 = '#9D02D7'; %magenta 2
c7 = '#0000FF'; %indigo
c8 = '#000000'; %black
sz = 60;        %size of data points
sz2 = sz*pres/620; %size of second data points
c = {c1; c2; c3; c4; c5; c6; c7; c8};

%% Plot expected versus measured moment arm
Ma = Bifemsh_Pam_adj3.MomentArm;                 %Calculated moment arm
G = -(Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque

figure
hold on
pp2 = plot(phiD,G,'DisplayName','\bf Expected $r_{\hat{k}}$');
ss2 = scatter(Angle, -ICRtoMuscle,'DisplayName','\bf Measured $r_{\hat{k}}$');
ppOS = plot(knee_angle_rMA,Bifemsh_MA,'--','DisplayName','\bf OpenSim $r_{\hat{k}}$');
title('\bf Expected vs measured $r_{\hat{k}}$', 'Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Z axis $r_{\hat{k}}$, m','Interpreter','latex')
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 12;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 12;
lgdMa = legend('Interpreter','latex');
lgdMa.FontSize = 12;
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values
strain = Bifemsh_Pam_adj3.Contraction;
relstrain = (strain)./KMAX;

figure
hold on
plot(phiD,relstrain,'DisplayName','\bf Expected \epsilon^*')
scatter(Angle,realRel,'DisplayName','\bf Measured \epsilon^*')
title('Expected vs measured \epsilon^*','Interpreter','tex')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('\epsilon^*','Interpreter','tex')
ax2 = gca;
ax2.FontSize = 12;
ax2.FontWeight = 'bold';
ax2.FontName = 'Arial';
ax2.YAxis.LineWidth = 2; ax2.YAxis.FontSize = 12;
ax2.XAxis.LineWidth = 2; ax2.XAxis.FontSize = 12;
lgdMa = legend('Interpreter','tex');
lgdMa.FontSize = 12;
hold off

%% Plot measured versus expected BPA length
MuscleLength = Bifemsh_Pam_adj3.MuscleLength-2*fitting-tendon_adj;

figure
hold on
plot(phiD,MuscleLength,'DisplayName','\bf Expected $l_{m}$')
scatter(Angle,InflatedLength,'DisplayName','\bf Measured $l_{m}$')
title('\bf Expected vs measured $l_{m}$','Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf $l_{m}$, m','Interpreter','latex')
ax3 = gca;
ax3.FontSize = 12;
ax3.FontWeight = 'bold';
ax3.FontName = 'Arial';
ax3.YAxis.LineWidth = 2; ax3.YAxis.FontSize = 12;
ax3.XAxis.LineWidth = 2; ax3.XAxis.FontSize = 12;
lgdMa = legend('Interpreter','latex');
lgdMa.FontSize = 12;
hold off

%% Plot Torque with interpolated values from MEAN and RMSE for theoretical values
%This way it is just over RoM measured at each pressure

figure
hold on
gca1 = gca;
gcf1 = gcf;


    for i = 1:size(ANG,1)
        T1 = 2*i;
        H1 = 2*i+1;
        Disp1{i} = sprintf('Theoretical, %.0f kPa',PRZ{i});
        Disp2{i} = sprintf('Theoretical optimized, %.0f kPa',PRZ{i});
        Disp3{i} = sprintf('Measured, %.0f kPa',PRZ{i});
        PL{i} = plot(ANG{i}, val{i},'--','Color',c{H1},'Linewidth',2,'DisplayName',Disp1{i});
        PL_opt{i} = plot(ANG{i}, val_opt{i},'Color',c{H1},'Linewidth',2,'DisplayName',Disp2{i});
        sc{i} = scatter(ANG{i},y{i},sz,'d','filled','MarkerFaceColor',c{H1},'DisplayName',Disp3{i});
    end


SC4 = scatter(Angle,TorqueHand,sz,'filled','MarkerFaceColor',c1,'DisplayName','Hybrid calc');
OST = plot(knee_angle_rT,Bifemsh_T,'-.','Color',c8,'LineWidth',1,'DisplayName','Human');

title('Iso. Torque vs {\theta_{k}}, {\phi}20mm Flexor, l_{rest} = 42.3cm','Interpreter','tex')
xlabel('Knee angle, \circ','FontWeight','bold','Interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','Interpreter','tex')
set(gca1,'FontSize', 12, 'LineWidth',2,'FontWeight', 'bold','FontName','Arial','XMinorGrid','off','XMinorTick','off','YMinorGrid','off','YMinorTick','off');
lgd = legend;
hold off

