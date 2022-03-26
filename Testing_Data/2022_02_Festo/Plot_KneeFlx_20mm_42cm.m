%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

f = fullfile('github/bipedal_robot/code/matlab');
qt = addpath(genpath(f));

restingLength = 0.423; %resting length, m
kmax = 0.31725; %Length at maximum contraction, m
dia = 20;

load KneeFlx_20mm_42cm.mat
Theoretical = TorqueR(:,3)';
%% Test 1 done with CALT load cell
%Test 1 == sheet FlxTest20mm from Results_table10mm_FullSize

%% Torque calculated from measurements
Angle = [-2	-11.5	-20	-27	-35	-44	-54.5	-50	-63	-64	-73.5	-82	-88	-103.5	-114	-115	-120];
Angle1 = Angle(1:10);
Angle2 = Angle(11:15);
Angle3 = Angle(16:17);

Torque = [-9.952075229	-8.30977619	-6.937208541	-5.83121108	-4.652816908	-3.594317472	-3.575454453	-3.026551477	-2.199767729	-0.922103339	-5.680697771	-4.039259835	-3.25137834	-1.326357902	-0.705989473	-1.826208162	-1.093715313];
Torque1 = Torque(1:10);
Torque2 = Torque(11:15);
Torque3 = Torque(16:17);
%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength = [401	395	389	387	382	380	378	374.5	371.5	369.5	354	352.5	350.5	345	343.5	338.5	339]/1000;
InflatedLength1 = InflatedLength(1:10);
InflatedLength2 = InflatedLength(11:15);
InflatedLength3 = InflatedLength(16:17);

ICRtoMuscle = [366.43	366.43	366.43	366.43	366.43	366.43	366.43	366.43	366.43	366.43	366.43	366.43	366.43	366.43	366.43	366.43	366.43]/1000;
ICRtoMuscle1 = ICRtoMuscle(1:10);
ICRtoMuscle2 = ICRtoMuscle(11:15);
ICRtoMuscle3 = ICRtoMuscle(16:17);

F = zeros(1,size(InflatedLength, 2));
F1 = F(1:10);
F2 = F(11:15);
F3 = F(16:17);

TorqueHand = zeros(1,size(InflatedLength, 2));

%load pressure where applicable
test = 1;
runsperseries = 17;

pres = zeros(1,runsperseries);
 

     for j = 1:runsperseries(i)
                file_name = sprintf('FlxTest%0.0f_%0.0f.mat', test,j);
                load(file_name,'Stats')
                pres(1,j) = Stats{'Mean',2};
     end
pres1 = pres(1:10);
pres2 = pres(11:15);
pres3 = pres(16:17);

for i = 1:size(InflatedLength, 2)
    F(i) = festo4(InflatedLength(i), restingLength, 10, pres(i), kmax);
    F_alt(i) = festo4(InflatedLength(i), restingLength, 10, 600, kmax);
    TorqueHand(i) = -ICRtoMuscle(i)*F(i);  %Moment will be negative because it causes flexion
    TorqueHand_alt(i) = -ICRtoMuscle(i)*F_alt(i);  %Potentially acheivable force
end
TorqueHand1 = TorqueHand(1:10);
TorqueHand2 = TorqueHand(11:15);
TorqueHand3 = TorqueHand(16:17);
%% Mean and RMSE
X = linspace(min(Angle),max(Angle),size(Angle,2));      %Range of motion
% mod = 'gauss1';
% fitOptions = fitoptions(mod, 'Normalize', 'on');
% [mdl1u, gof1] = fit(Angle',Torque',mod,fitOptions)
% TorqueStdu = gof1.rmse
% TorqueMeanu = feval(mdl1u,X)';
% 
% [mdl2u, gof2] = fit(Angle',TorqueHand',mod,fitOptions);
% HandStdu = gof2.rmse;
% HandMeanu = feval(mdl2u,X)';

modp = 'poly3';
fitOp = fitoptions(modp,'Normalize','on','Robust','on');
[mdl1, gofp1] = fit(Angle1',Torque1',modp,fitOp)
TorqueStd = gofp1.rmse
TorqueMean = feval(mdl1,X)'

[mdl2, gofp2] = fit(Angle2',TorqueHand2',modp,fitOp);
HandStd = gofp2.rmse;
HandMean = feval(mdl2,X)'

%% Plotting with polynomial solver
figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Flexor, 48.5cm long')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
gca1 = gca;
gcf1 = gcf;
set(gcf,'Position',[1 384 950 612]);
set(gca,'FontSize', 18, 'FontWeight', 'bold','XMinorGrid','on','XMinorTick','on','YMinorGrid','on','YMinorTick','on');
plot(phiD, Theoretical,'Color',[0 0.4470 0.7410],'Linewidth',2,'DisplayName','Theoretical Calculation')

Xnew=[X,fliplr(X)];
Y1=[TorqueMean+TorqueStd,fliplr(TorqueMean-TorqueStd)];
Y2=[HandMean+HandStd,fliplr(HandMean-HandStd)];
plot(X,TorqueMean,'--k','Linewidth',2,'DisplayName','Torque mean, scale')
fill(Xnew,Y1,[1 0.4 0.8],'DisplayName','Scale torque SD','FaceAlpha',0.25);
plot(X,HandMean,'--r','Linewidth',2,'DisplayName','Torque mean, hand')
fill(Xnew,Y2,[.6 1.0 .6],'DisplayName','Hand torque SD','FaceAlpha',0.25);


sz = 60;
c1 = [0.8500 0.3250 0.0980]; % color, burnt orange
c2 = [0.6350 0.0780 0.1840]; %color, red/violet
c3 = [0 0.4470 0.7410]; %color, navy blue
c4 = [0.4660 0.6740 0.1880]; %color, moss green
c5 = [0.9290 0.6940 0.1250]; %color, dark yellow
sc1 = scatter(Angle,Torque,sz,'d','CData',c1,'DisplayName','BB LC');
sc2 = scatter(Angle,TorqueHand,sz,'filled','CData',c1,'DisplayName','BB hand');

legend
hold off
%% Plotting with nonlinear solution
figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Flexor, 48.5cm long')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
gca2 = gca;
gcf2 = gcf;
set(gcf,'Position',[960 384 950 612]);
set(gca,'FontSize', 18, 'FontWeight', 'bold','XMinorGrid','on','XMinorTick','on','YMinorGrid','on','YMinorTick','on');
plot(phiD, Theoretical,'Color',[0 0.4470 0.7410],'Linewidth',2,'DisplayName','Theoretical Calculation')

Xnew=[X,fliplr(X)];
Y1=[TorqueMeanu+TorqueStdu,fliplr(TorqueMeanu-TorqueStdu)];
Y2=[HandMeanu+HandStdu,fliplr(HandMeanu-HandStdu)];
plot(X,TorqueMeanu,'--k','Linewidth',2,'DisplayName','Torque mean, scale')
fill(Xnew,Y1,[1 0.4 0.8],'DisplayName','Fish scale SD','FaceAlpha',0.25);
plot(X,HandMeanu,'--r','Linewidth',2,'DisplayName','Torque mean, hand')
fill(Xnew,Y2,[.6 1.0 .6],'DisplayName','Hand torque SD','FaceAlpha',0.25);

sc3 = scatter(Angle,Torque,sz,'d','CData',c1,'DisplayName','JM LC');
sc4 = scatter(Angle,TorqueHand,sz,'filled','CData',c1,'DisplayName','JM hand');

legend
hold off
%% Plot error and standard deviation as bar graphs
% xb=categorical({'Theoretical','Scale','Hand'});
% xb = reordercats(xb,{'Theoretical','Scale','Hand'});
% yb = [mean(Theoretical) mean(TorqueMean) mean(TorqueHand)];
% std_dev = [0 mean(TorqueStd) mean(HandStd)];
% figure
% hold on
% b = bar(xb,yb,'FaceColor',[0 0.8 1.0]);
% gca3 = gca;
% gcf3 = gcf;
% ylabel('Torque, N*m')
% title('Mean Torque Values and SD for Each Calculation Method')
% set(gcf,'Position',[0 0 950 612]);
% set(gca,'FontSize', 18, 'FontWeight', 'bold');
% b.CData = [0 0.4470 0.7410; 0  0  0; 1  0  0];
% errb = errorbar(yb,std_dev ,'LineStyle','none','LineWidth',4,'CapSize',20);
% hold off