%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

f = fullfile('github/bipedal_robot/code/matlab');
qt = addpath(genpath(f));

restingLength = 0.557; %resting length, m
kmax = 0.399; %Length at maximum contraction, m

load 'KneeExt40_OldCalcMethod.mat' theta_K vas_int_tqZ vas_lat_tqZ vas_med_tqZ %Human values, Ben's old method
vas_int_H = vas_int_tqZ(:,1); 
vas_lat_H = vas_lat_tqZ(:,1); 
vas_med_H = vas_med_tqZ(:,1);
vas_tqZ = vas_int_H+vas_lat_H+vas_med_H;

load 'KneeExt40_Comparison.mat'
TheoreticalR1 = Vas_Pam_dist.Torque(:,3)';   %distal insertion ring
TheoreticalR2 = Vas_Pam_prox.Torque(:,3)';       %proximal insertion ring
TorqueHnew = Torque1(:,3);                     %Human calculations, class

Tab = readmatrix('OpenSim_Vasti_Results.txt');
knee_angle_r = Tab(:,2)';           %Angle values directly from O
vas_int_r = Tab(:,3)';              %Torque values directly from OpenSim
vas_lat_r = Tab(:,3)';
vas_med_r = Tab(:,4)';
vas_r = vas_int_r+vas_lat_r+vas_med_r; %combine all the OpenSim torque values
%% Test 1 done with CALT load cell. Tests 2 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet ExtTest10mm_2 from Results_table10mm_pinned_LoadCell
%Test 2 == sheet ExtTest10mm_1 from Results_table10mm_pinned_FishScale
%% Torque calculated from measurements

Angle = [-117	-106	-96	-66.5	-74.5	-85.5	-53.5	-53.5	-17	-17	-42.5	-30	-6	6.5 6 -16];
Angle1 = Angle(1:7);            %100 kPa
Angle2 = Angle(8:12);           %230 kPa
Angle3 = Angle(13:14);          %380 kPa
Angle4 = 6;                     %604 kPa
Angle5 = -16;                   %530 kPa

Torque = [6.853017236	3.061644838	1.756743486	0.428753004	0.962918245	1.433251369	0.353097824	6.225298962	1.054104656	1.235368116	6.042042931	3.234337401	10.36051915	6.080949806];
Torque1 = Torque(1:7);
Torque2 = Torque(8:12);
Torque3 = Torque(13:14);
Torque4 = 19.9071246;
Torque5 = 29.4178228;

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

contraction = [0.043087971	0.049371634	0.066427289	0.059245961	0.062836625	0.062836625	0.055655296	0.156193896	0.172351885	0.172351885	0.170556553	0.179533214	0.21005386	0.228007181 0.283662478];
contraction1 = contraction(1:7);
contraction2 = contraction(8:12);
contraction3 = contraction(13:14);
contraction4 = contraction(15);

ICRtoMuscle = [60	65	57	73	62	67	74	65	70	75	62.5	69	60	64 65]/1000;
ICRtoMuscle1 = ICRtoMuscle(1:7);
ICRtoMuscle2 = ICRtoMuscle(8:12);
ICRtoMuscle3 = ICRtoMuscle(13:14);
ICRtoMuscle4 = ICRtoMuscle(15);


F1 = zeros(1,size(ICRtoMuscle1, 2));
F2 = zeros(1,size(ICRtoMuscle2, 2));
F3 = zeros(1,size(ICRtoMuscle3, 2));
F4 = 0;
F = [F1 F2 F3 F4];
F_alt = F;

TorqueHand1 = zeros(1,size(ICRtoMuscle1, 2));
TorqueHand2 = zeros(1,size(ICRtoMuscle2, 2));
TorqueHand3 = zeros(1,size(ICRtoMuscle3, 2));
TorqueHand4 = 0;
TorqueHand = [TorqueHand1, TorqueHand2, TorqueHand3, TorqueHand4];

Hand_alt = zeros(1,size(TorqueHand,2));

%load pressure where applicable
test = 1;
runsperseries = 14;

    pres1 = zeros(1,runsperseries);
    
    for j = 1:runsperseries-1
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test,j);
                load(file_name,'Stats')
                pres1(1,j) = Stats{'Mean',2};
    end

pres2 = 604.65;
pres = [pres1 pres1(13) pres2];


F = festo4(40, contraction, pres);
TorqueHand = ICRtoMuscle.*F;
F_alt = festo4(40, contraction, 600);
Hand_alt = ICRtoMuscle.*F_alt;

TorqueHand1 = TorqueHand(1:7);
TorqueHand2 = TorqueHand(8:12);
TorqueHand3 = TorqueHand(13:14);
TorqueHand4 = TorqueHand(15);
Hand_alt1 = Hand_alt(1:7);
Hand_alt2 = Hand_alt(8:12);
Hand_alt3 = Hand_alt(13:14);
Hand_alt4 = Hand_alt(15);

%% Mean and RMSE
X = linspace(min(Angle),max(Angle),size(Angle,2));      %Range of motion
% mod = 'sin2';
% fitOptions = fitoptions(mod, 'Normalize', 'on');
% [mdl1u, gof1] = fit(Angle',Torque',mod,fitOptions)
% TorqueStdu = gof1.rmse
% TorqueMeanu = feval(mdl1u,X)';
% 
% [mdl2u, gof2] = fit(Angle',TorqueHand',mod,fitOptions);
% HandStdu = gof2.rmse;
% HandMeanu = feval(mdl2u,X)';

mod = 'poly2';
lin = 'poly1';
fitOp1 = fitoptions(mod,'Normalize','on','Robust','on');
fitOp2 = fitoptions(lin,'Normalize','on','Robust','on');

[mdl1, gof1] = fit(Angle1',Torque1',mod,fitOp1); %100 kPa
TorqueStd1 = gof1.rmse;
TorqueMean1 = feval(mdl1,X)';
[mdl2, gof2] = fit(Angle2',Torque2',lin,fitOp2); %230 kPa
TorqueStd2 = gof2.rmse;
TorqueMean2 = feval(mdl2,X)';
[mdl3, gof3] = fit(Angle3',Torque3',lin,fitOp2); %380 kPa
TorqueStd3 = gof3.rmse;
TorqueMean3 = feval(mdl3,X)';

[mdl1H, gof1H] = fit(Angle1',TorqueHand1',mod,fitOp1); %100 kPa
TorqueStd1H = gof1H.rmse;
TorqueMean1H = feval(mdl1H,X)';

[mdl2H, gof2H] = fit(Angle2',TorqueHand2',lin,fitOp2); %230 kPa
TorqueStd2H = gof2H.rmse;
TorqueMean2H = feval(mdl2H,X)';

[mdl3H, gof3H] = fit(Angle3',TorqueHand3',lin,fitOp2); %380 kPa
TorqueStd3H = gof3H.rmse;
TorqueMean3H = feval(mdl3H,X)';

%% Plotting polynomial solver
sz = 50;
c = [0.8 0.2 0.2]; % color

figure
hold on
title('Isometric Torque vs Knee Angle, 40mm Extensor, 55.7cm long','interpreter','latex')
xlabel('degrees Flexion(-),Extension(+)','interpreter','latex')
ylabel('Torque, N*m','interpreter','latex')
gca1 = gca;
gcf1 = gcf;
set(gcf,'Position',[4 50 800 750]);
set(gca,'FontWeight', 'bold','XLim',[-120 20],'XMinorGrid','off','XMinorTick','off','YMinorGrid','off','YMinorTick','off','YLim',[-10 500]);
plot(knee_angle_r, vas_r,'LineStyle','-','Color',[0 0 1],'Linewidth',3,'DisplayName','Direct values from OpenSim')
plot(phiD, TheoreticalR1,'LineStyle','-.','Color',[0 0.2 0.2],'Linewidth',2,'DisplayName','Robot, distal ring')
plot(phiD, TheoreticalR2,'LineStyle',':','Color',[1 0.4470 0.7410],'Linewidth',2,'DisplayName','Robot, proximal ring')
plot(theta_K, vas_tqZ,'LineStyle','--','Color',[0.6 0 0.8],'Linewidth',2,'DisplayName','Human values (Bolen)')
plot(phiD, TorqueHnew,'LineStyle','-.','Color',[0.4 1 0.8],'Linewidth',2,'DisplayName','Human values (Morrow)')
legend

scatter(Angle1,Torque1,sz,'+','MarkerFaceColor','g','LineWidth',3,'DisplayName','100 kPa, Measured');
scatter(Angle2,Torque2,sz,'+','MarkerFaceColor','r','LineWidth',3,'DisplayName','230 kPa, Measured');
scatter(Angle3,Torque3,sz,'+','MarkerFaceColor','b','LineWidth',3,'DisplayName','380 kPa, Measured');
scatter(Angle4,Torque4,sz,'+','CData',c,'LineWidth',3,'DisplayName','604 kPa, Measured');
scatter(Angle1,TorqueHand1,sz,'g','filled','DisplayName','100 kPa, back calc');
scatter(Angle2,TorqueHand2,sz,'r','filled','DisplayName','230 kPa, back calc');
scatter(Angle3,TorqueHand3,sz,'b','filled','DisplayName','380 kPa, back calc');
scatter(Angle4,TorqueHand4,sz,'filled','CData',c,'DisplayName','604 kPA, back calc');

scatter(Angle(1:15),Hand_alt,sz,'Marker','v','CData',[0.6350 0.0780 0.1840],'DisplayName','projected 600 kPa Torque for data');

legend
hold off

