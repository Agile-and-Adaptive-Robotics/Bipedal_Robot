%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

restingLength = 0.457; %resting length, m
kmax = 0.380; %Length at maximum contraction, m

load KneeFlxPin_10mm_46cm.mat
Theoretical = TorqueR(:,3)';
%% Tests 1 & 2 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet FlxTest10mm_3 from Results_table10mm_pinned_LoadCell
%Test 2 == sheet FlxTest10mm_3 from Results_table10mm_pinned_FishScale
%% Torque calculated from measurements

Angle1 = [-41	-45	-57	-63	-77	-75	-80	-90	-95	-105];
Angle2 = [-30	-32	-36.5	-41	-51	-52.5	-56.5	-63	-66.5	-74.5	-80.5	-80	-82	-86	-90	-85	-81.5	-76	-70	-62	-58	-49	-45	-37	-33.5	-26.5];
Angle = [Angle1, Angle2];

Torque1 = [-17.62398227	-14.31258839	-10.54411571	-7.71470304	-5.263349252	-5.04823533	-2.832979547	-1.807118322	-1.02853932	-0.561617562];
Torque2 = [-20.00095306	-19.62412624	-17.74104856	-14.64175686	-13.34974574	-11.21468685	-8.847415849	-6.225462211	-5.473176171	-4.159138583	-1.821998215	-1.821998215	-1.293882791	-0.526243409	0	-1.293882791	-2.098936493	-3.894225298	-5.473176171	-7.788450595	-11.02208726	-14.67316852	-18.04885624	-21.08709753	-22.77466623	-19.33565341];
Torque = [Torque1, Torque2];

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength1 = [428	418	409.5	405.5	395	396	390	348	381	376]/1000;
InflatedLength2 = [438	438	434	424	419	415	408	401	399	394	388	388	386	381	380	384	386	390	394	400	407	413	422	426	433	438]/1000;
InflatedLength = [InflatedLength1, InflatedLength2];

ICRtoMuscle1 = [80	81	78	76	62	62.5	60	48	34	26]/1000;
ICRtoMuscle2 = [85	85	85	88	90	88	83	80	75	70	60	64	60	60	63	65	65	68	75	79	84	90	92	90	89	85]/1000;
ICRtoMuscle = [ICRtoMuscle1, ICRtoMuscle2];

F1 = zeros(1,size(InflatedLength1, 2));
F2 = zeros(1,size(InflatedLength2, 2));
F = [F1, F2];

TorqueHand1 = zeros(1,size(InflatedLength1, 2));
TorqueHand2 = zeros(1,size(InflatedLength2, 2));
TorqueHand = [TorqueHand1, TorqueHand2];

%load pressure where applicable
test = 3;
runsperseries = 10;

    pres1 = zeros(1,runsperseries);
    
    for j = 1:runsperseries
                file_name = sprintf('FlxTest%0.0f_%0.0f.mat', test,j);
                load(file_name,'Stats')
                pres1(1,j) = Stats{'Mean',2};
    end

pres2 = 606*ones(1,size(InflatedLength, 2));
pres = [pres1 pres2];

for i = 1:size(InflatedLength, 2)
    F(i) = festo3(InflatedLength(i), restingLength, 10, pres(i), kmax);    
    TorqueHand(i) = -ICRtoMuscle(i)*F(i);  %Torque will be negative because it is causing flexion
end
TorqueHand1 = TorqueHand(1:size(TorqueHand1,2));
TorqueHand2 = TorqueHand((size(TorqueHand1,2)+1):size(TorqueHand,2));

%% Mean and RMSE
X = linspace(min(Angle),max(Angle),size(Angle,2));      %Range of motion
mod = 'gauss1';
fitOptions = fitoptions(mod, 'Normalize', 'on');
[mdl1u, gof1] = fit(Angle',Torque',mod,fitOptions);
TorqueStdu = gof1.rmse;
TorqueMeanu = feval(mdl1u,X)';

[mdl2u, gof2] = fit(Angle',TorqueHand',mod,fitOptions);
HandStdu = gof2.rmse;
HandMeanu = feval(mdl2u,X)';

modp = 'poly3';
fitOp = fitoptions(modp,'Normalize','on');
[mdl1, gofp1] = fit(Angle',Torque',modp,fitOptions)
TorqueStd = gofp1.rmse
TorqueMean = feval(mdl1,X)';

[mdl2, gofp2] = fit(Angle',TorqueHand',modp,fitOptions)
HandStd = gofp2.rmse
HandMean = feval(mdl2,X)';

%% Plotting with polynomial solver
figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Flexor, 45.5cm long')
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
fill(Xnew,Y1,[1 0.4 0.8],'DisplayName','Fish scale std','FaceAlpha',0.25);
fill(Xnew,Y2,[.6 1.0 .6],'DisplayName','Hand torque std','FaceAlpha',0.25);
plot(X,TorqueMean,'--k','Linewidth',2,'DisplayName','Torque mean, scale')
plot(X,HandMean,'--r','Linewidth',2,'DisplayName','Torque mean, hand')

c = [0.8500 0.3250 0.0980]; % color
sz = 50;
c1 = [0.8500 0.3250 0.0980]; % color, burnt orange
c2 = [0.6350 0.0780 0.1840]; %color, red/violet
c3 = [0 0.4470 0.7410]; %color, navy blue
c4 = [0.4660 0.6740 0.1880]; %color, moss green
c5 = [0.9290 0.6940 0.1250]; %color, dark yellow
scatter(Angle1,Torque1,sz,'d','CData',c1,'DisplayName','JM LC');
scatter(Angle2,Torque2,sz,'d','CData',c3,'DisplayName','BB FS');
scatter(Angle1,TorqueHand1,sz,'filled','CData',c1,'DisplayName','JM hand');
scatter(Angle2,TorqueHand2,sz,'filled','CData',c3,'DisplayName','BB hand');

legend
hold off

%% Plotting with nonlinear solver
figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Flexor, 45.5cm long')
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
fill(Xnew,Y1,[1 0.4 0.8],'DisplayName','Fish scale std','FaceAlpha',0.25);
fill(Xnew,Y2,[.6 1.0 .6],'DisplayName','Hand torque std','FaceAlpha',0.25);
plot(X,TorqueMeanu,'--k','Linewidth',2,'DisplayName','Torque mean, scale')
plot(X,HandMeanu,'--r','Linewidth',2,'DisplayName','Torque mean, hand')

c = [0.8500 0.3250 0.0980]; % color
sz = 50;
c1 = [0.8500 0.3250 0.0980]; % color, burnt orange
c2 = [0.6350 0.0780 0.1840]; %color, red/violet
c3 = [0 0.4470 0.7410]; %color, navy blue
c4 = [0.4660 0.6740 0.1880]; %color, moss green
c5 = [0.9290 0.6940 0.1250]; %color, dark yellow
scatter(Angle1,Torque1,sz,'d','CData',c1,'DisplayName','JM LC');
scatter(Angle2,Torque2,sz,'d','CData',c3,'DisplayName','BB FS');
scatter(Angle1,TorqueHand1,sz,'filled','CData',c1,'DisplayName','JM hand');
scatter(Angle2,TorqueHand2,sz,'filled','CData',c3,'DisplayName','BB hand');

legend
hold off
%% Plot error and standard deviation as bar graphs
xb=categorical({'Theoretical','Scale','Hand'});
xb = reordercats(xb,{'Theoretical','Scale','Hand'});
yb = [mean(Theoretical) mean(TorqueMean) mean(TorqueHand)];
std_dev = [0 mean(TorqueStd) mean(HandStd)];
figure
hold on
b = bar(xb,yb,'FaceColor',[0 0.8 1.0]);
gca3 = gca;
gcf3 = gcf;
ylabel('Torque, N*m')
title('Mean Torque Values and SD for Each Calculation Method')
set(gcf,'Position',[0 0 950 612]);
set(gca,'FontSize', 18, 'FontWeight', 'bold');
b.CData = [0 0.4470 0.7410; 0  0  0; 1  0  0];
errb = errorbar(yb,std_dev ,'LineStyle','none','LineWidth',4,'CapSize',20);
hold off