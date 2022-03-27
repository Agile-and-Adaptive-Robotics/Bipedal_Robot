%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

restingLength = 0.415; %resting length, m
kmax = 0.349; %Length at maximum contraction, m

load 'KneeExtPin_10mm_42cm_comparisons.mat'
TorqueR1 = Vas_Pam_ideal.Torque;
TorqueR2 = Vas_Pam_real1.Torque;
TorqueR3 = Vas_Pam_real2.Torque; 
TorqueR4 = Vas_Pam_real3.Torque;
TorqueR5 = Vas_Pam_tendon_ideal.Torque;
TorqueR6 = Vas_Pam_tendon_real.Torque;
TorqueR7 = Vas_Pam_slip.Torque;

clear fit
%% Test 1 done with CALT load cell. Tests 2 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet ExtTest10mm_4 from Results_table10mm_pinned_LoadCell
%Test 2 == sheet ExtTest10mm_5 from Results_table10mm_pinned_LoadCell
%Test 3 == sheet ExtTest10mm_6 from Results_table10mm_pinned_LoadCell
%% Torque calculated from measurements

Angle1 = [-76.5	-69	-51.5	-40.5	-39.5	-55.5	-35	-25.5	-20	-12.5	-5	0.5	4	12.5	20	28	35	23	20	11	1	-16	-24	-44	-53.4	-66];
Angle2 = [-120	-109	-101	-97	-89	-87	-81.5	-75	-74	-59.5	-48.5	-3	-30	-18.5];
Angle3 = [1.5	-5.5	-14.5	-20	-28.5	-35.5	-40	-50.5	-60	-66.5	-72.5	-78.5	-88	-94.5	-105	-112.5	-119.5];
Angle = [Angle1, Angle2, Angle3];

Torque1 = [8.50477362	7.539449285	6.619468483	5.835581619	6.120370289	7.005900186	6.178291184	6.391414397	5.96804559	6.099105123	5.93740855	4.880958611	4.375659895	3.388083398	2.795862193	1.758237576	0.306753172	1.725902869	3.204095044	5.145222736	6.304242352	6.817870142	7.511919159	7.911648561	8.110527045	8.862066797];
Torque2 = [4.032777251	3.725074705	3.501289111	3.342545878	3.070623054	2.844765784	2.844621771	2.773344025	2.674409973	2.382410213	1.680520864	0.394827933	1.758713268	1.550440294];
Torque3 = [0.55424608	0.941512328	1.76545047	1.79737823	2.650554833	2.7667878	2.211239585	3.275126705	3.423290801	3.611274245	3.687991538	3.833457994	3.587645288	3.554650905	3.949198619	4.568831269	4.408023267];
Torque = [Torque1, Torque2, Torque3];

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.

InflatedLength1 = [413.5	405.5	402	400	403	412	401	398	392.5	385	385.5	371.5	365	361	356	355	346	348	354	365	376	398	395	400	407	413]/1000;
InflatedLength2 = [408	402	380	394	390.5	387	383.5	383	381	371	363	340	355	353]/1000;
InflatedLength3 = [340	340	350	350.5	354	362	366	363	368	373	376	383	384	394	391	394.5	395]/1000;
InflatedLength = [InflatedLength1, InflatedLength2 InflatedLength3];

ICRtoMuscle1 = [28	30.5	31.5	31	32	31	35	32	41	44.5	45.5	55.5	66.5	62	65	68.5	71	67	60	55	55	42	38	31	29	28]/1000;
ICRtoMuscle2 = [42	42	42	42	41.5	36	37	36	36	30.5	22	42	28	24]/1000;
ICRtoMuscle3 = [55	43.5	42.5	37	35	33	35	34	35	34	36	36	36	36	36	36	36]/1000;
ICRtoMuscle = [ICRtoMuscle1, ICRtoMuscle2, ICRtoMuscle3];

F1 = zeros(1,size(InflatedLength1, 2));
F2 = zeros(1,size(InflatedLength2, 2));
F3 = zeros(1,size(InflatedLength3, 2));
F = [F1, F2 F3];

TorqueHand1 = zeros(1,size(InflatedLength1, 2));
TorqueHand2 = zeros(1,size(InflatedLength2, 2));
TorqueHand3 = zeros(1,size(InflatedLength3, 2));
TorqueHand = [TorqueHand1, TorqueHand2 TorqueHand3];

%load pressure where applicable
test = [4 5 6];
runsperseries = [26 14 17];

pres1 = zeros(1,runsperseries(1));
pres2 = zeros(1,runsperseries(2));
pres3 = zeros(1,runsperseries(3));
    for i = 1:3
        for j = 1:runsperseries(i)
            if j == 1
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres1(1,j) = Stats{'Mean',2};
            elseif j ==2
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres2(1,j) = Stats{'Mean',2};
            else
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres3(1,j) = Stats{'Mean',2};
            end
        end       
    end

pres = [pres1 pres2 pres3];
pressy = nonzeros(pres);
pres = mean(pressy).*ones(1,size(pres,2));

for i = 1:size(InflatedLength, 2)
    F(i) = festo3(InflatedLength(i), restingLength, 10, pres(i), kmax);    
    TorqueHand(i) = ICRtoMuscle(i)*F(i);
end
TorqueHand1 = TorqueHand(1:size(TorqueHand1,2));
TorqueHand2 = TorqueHand(((size(TorqueHand1,2)+1)):(size(TorqueHand1,2)+size(TorqueHand2,2)));
TorqueHand3 = TorqueHand((size(TorqueHand1,2)+size(TorqueHand2,2)+1):size(TorqueHand,2));
%% Mean and RMSE
X1 = linspace(min(Angle1),max(Angle1));      %Range of motion
X2 = linspace(min(Angle2),max(Angle2));      %Range of motion
X3 = linspace(min(Angle3),max(Angle3));      %Range of motion

modp = 'poly3';
fitOp = fitoptions(modp,'Normalize','on','Robust','on');

[mdl1, gof1] = fit(Angle1',Torque1',modp,fitOp)
TorqueStd1 = gof1.rmse
TorqueMean1 = feval(mdl1,X1)';
[mdl1p, gofp1] = fit(Angle1',TorqueHand1',modp,fitOp);
HandStd1 = gofp1.rmse;
HandMean1 = feval(mdl1p,X1)';

[mdl2, gof2] = fit(Angle2',Torque2',modp,fitOp)
TorqueStd2 = gof2.rmse
TorqueMean2 = feval(mdl2,X2)';
[mdl2p, gofp2] = fit(Angle2',TorqueHand2',modp,fitOp);
HandStd2 = gofp2.rmse;
HandMean2 = feval(mdl2p,X2)';

[mdl3, gof3] = fit(Angle3',Torque3',modp,fitOp)
TorqueStd3 = gof3.rmse
TorqueMean3 = feval(mdl3,X3)';
[mdl3p, gofp3] = fit(Angle3',TorqueHand3',modp,fitOp);
HandStd3 = gofp3.rmse;
HandMean3 = feval(mdl3p,X3)';

%% Plotting polynomial fit
figure
hold on
title('Isometric Z axis Torque, 10mm Extensor, 41.5cm long')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
hold on
gca1 = gca;
gcf1 = gcf;
set(gcf,'Position',[1 384 950 612]);
set(gca,'FontSize', 12, 'FontWeight', 'bold','XMinorGrid','on','XMinorTick','on','YMinorGrid','on','YMinorTick','on');
c1 = [0.8500 0.3250 0.0980]; % color, burnt orange
c2 = [0.6350 0.0780 0.1840]; %color, red/violet
c3 = [0 0.4470 0.7410]; %color, navy blue
c4 = [0.4660 0.6740 0.1880]; %color, moss green
c5 = [0.9290 0.6940 0.1250]; %color, dark yellow

plot(phiD, TorqueR1(:, 3),'Color',c1,'Linewidth',2,'DisplayName','Ideal');
plot(phiD,TorqueR2(:, 3),'--^','DisplayName','Realistic 1','Linewidth',2)
plot(phiD, TorqueR3(:, 3),'--<','DisplayName','Realistic 2','Linewidth',2)
plot(phiD, TorqueR4(:, 3),'--v','DisplayName','Realistic 3','Linewidth',2)
plot(phiD, TorqueR5(:, 3),'-','Color',c3,'DisplayName','Tendon Ideal','Linewidth',2)
plot(phiD, TorqueR6(:, 3),'-.','DisplayName','Tendon Real','Linewidth',2)
plot(phiD, TorqueR7(:,3),':','Color',c2,'DisplayName','Slippage+Tendon','Linewidth',2)

X1new=[X1,fliplr(X1)];
Y1=[TorqueMean1+TorqueStd1,fliplr(TorqueMean1-TorqueStd1)];
Y1p=[HandMean1+HandStd1,fliplr(HandMean1-HandStd1)];
plot(X1,TorqueMean1,'--k','Linewidth',2,'DisplayName','Test1 mean, scale')
fill(X1new,Y1,[.4 .4 1],'DisplayName','Test1, scale SD','FaceAlpha',0.25);
plot(X1,HandMean1,'--r','Linewidth',2,'DisplayName','Test1 mean, hand')
fill(X1new,Y1p,[.6 1.0 .6],'DisplayName','Test1 mean, hand SD','FaceAlpha',0.25);

X2new=[X2,fliplr(X2)];
Y2=[TorqueMean2+TorqueStd2,fliplr(TorqueMean2-TorqueStd2)];
Y2p=[HandMean2+HandStd2,fliplr(HandMean2-HandStd2)];
plot(X2,TorqueMean2,'--k','Linewidth',2,'DisplayName','Test2 mean, scale')
fill(X2new,Y2,[0.4 0 0.6],'DisplayName','Test2, scale SD','FaceAlpha',0.25);
plot(X2,HandMean2,'--r','Linewidth',2,'DisplayName','Test2 mean, hand')
fill(X2new,Y2p,[0 0.6 1],'DisplayName','Test2 mean, hand SD','FaceAlpha',0.25);

X3new=[X3,fliplr(X3)];
Y3=[TorqueMean3+TorqueStd3,fliplr(TorqueMean3-TorqueStd3)];
Y3p=[HandMean3+HandStd3,fliplr(HandMean3-HandStd3)];
plot(X3,TorqueMean3,'--k','Linewidth',2,'DisplayName','Test3 mean, scale')
fill(X3new,Y3,[0.6 0 0.6],'DisplayName','Test3, scale SD','FaceAlpha',0.25);
plot(X3,HandMean3,'--r','Linewidth',2,'DisplayName','Test3 mean, hand')
fill(X3new,Y3p,[.2 1.0 1],'DisplayName','Test3 mean, hand SD','FaceAlpha',0.25);

sz = 50;
sz2 = 100;
scatter(Angle1,Torque1,sz2,'+','CData',c1,'DisplayName','No tendon, LC');
scatter(Angle1,TorqueHand1,sz,'filled','CData',c1,'DisplayName','No tendon, hand');
scatter(Angle2,Torque2,sz2,'+','CData',c2,'DisplayName','Tendon+slip, LC');
scatter(Angle2,TorqueHand2,sz,'filled','CData',c2,'DisplayName','Tendon+slip, hand');
scatter(Angle3,Torque3,sz2,'+','CData',c3,'DisplayName','Tendon, LC');
scatter(Angle3,TorqueHand3,sz,'filled','CData',c3,'DisplayName','Tendon, hand');

legend
hold off
%% Plot compare no tendon
figure
hold on
title('Z axis Torque, 10mm Extensor, 41.5cm long, no tendon')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
hold on
gca1 = gca;
gcf1 = gcf;
set(gcf,'Position',[1 384 950 612]);
set(gca,'FontSize', 12, 'FontWeight', 'bold','XMinorGrid','on','XMinorTick','on','YMinorGrid','on','YMinorTick','on');


plot(phiD, TorqueR1(:, 3),'Color',c1,'Linewidth',2,'DisplayName','Ideal');
plot(phiD,TorqueR2(:, 3),'--^','DisplayName','Realistic 1','Linewidth',2)
plot(phiD, TorqueR3(:, 3),'--<','DisplayName','Realistic 2','Linewidth',2)
plot(phiD, TorqueR4(:, 3),'--v','DisplayName','Realistic 3','Linewidth',2)

X1new=[X1,fliplr(X1)];
Y1=[TorqueMean1+TorqueStd1,fliplr(TorqueMean1-TorqueStd1)];
Y1p=[HandMean1+HandStd1,fliplr(HandMean1-HandStd1)];
plot(X1,TorqueMean1,'--k','Linewidth',2,'DisplayName','Test1 mean, scale')
fill(X1new,Y1,[.4 .4 1],'DisplayName','Test1, scale SD','FaceAlpha',0.25);
plot(X1,HandMean1,'--r','Linewidth',2,'DisplayName','Test1 mean, hand')
fill(X1new,Y1p,[.6 1.0 .6],'DisplayName','Test1 mean, hand SD','FaceAlpha',0.25);


scatter(Angle1,Torque1,sz2,'+','CData',c1,'DisplayName','No tendon, LC');
scatter(Angle1,TorqueHand1,sz,'filled','CData',c1,'DisplayName','No tendon, hand');

legend
hold off


%% Plot compare tendon results
figure
hold on
title('Z axis Torque, 10mm Extensor, 41.5cm long, 22 mm tendon')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
hold on
gca1 = gca;
gcf1 = gcf;
set(gcf,'Position',[1 384 950 612]);
set(gca,'FontSize', 12, 'FontWeight', 'bold','XMinorGrid','on','XMinorTick','on','YMinorGrid','on','YMinorTick','on');
c1 = [0.8500 0.3250 0.0980]; % color, burnt orange
c2 = [0.6350 0.0780 0.1840]; %color, red/violet
c3 = [0 0.4470 0.7410]; %color, navy blue
c4 = [0.4660 0.6740 0.1880]; %color, moss green
c5 = [0.9290 0.6940 0.1250]; %color, dark yellow

plot(phiD, TorqueR5(:, 3),'-','Color',c3,'DisplayName','Tendon Ideal','Linewidth',2)
plot(phiD, TorqueR6(:, 3),'-.','DisplayName','Tendon Real','Linewidth',2)
plot(phiD, TorqueR7(:,3),':','Color',c2,'DisplayName','Slippage+Tendon','Linewidth',2)

X2new=[X2,fliplr(X2)];
Y2=[TorqueMean2+TorqueStd2,fliplr(TorqueMean2-TorqueStd2)];
Y2p=[HandMean2+HandStd2,fliplr(HandMean2-HandStd2)];
plot(X2,TorqueMean2,'--k','Linewidth',2,'DisplayName','Test2 mean, scale')
fill(X2new,Y2,[0.4 0 0.6],'DisplayName','Test2, scale SD','FaceAlpha',0.25);
plot(X2,HandMean2,'--r','Linewidth',2,'DisplayName','Test2 mean, hand')
fill(X2new,Y2p,[0 0.6 1],'DisplayName','Test2 mean, hand SD','FaceAlpha',0.25);

X3new=[X3,fliplr(X3)];
Y3=[TorqueMean3+TorqueStd3,fliplr(TorqueMean3-TorqueStd3)];
Y3p=[HandMean3+HandStd3,fliplr(HandMean3-HandStd3)];
plot(X3,TorqueMean3,'--k','Linewidth',2,'DisplayName','Test3 mean, scale')
fill(X3new,Y3,[0.6 0 0.6],'DisplayName','Test3, scale SD','FaceAlpha',0.25);
plot(X3,HandMean3,'--r','Linewidth',2,'DisplayName','Test3 mean, hand')
fill(X3new,Y3p,[.2 1.0 1],'DisplayName','Test3 mean, hand SD','FaceAlpha',0.25);

scatter(Angle2,Torque2,sz2,'+','CData',c2,'DisplayName','Tendon+slip, LC');
scatter(Angle2,TorqueHand2,sz,'filled','CData',c2,'DisplayName','Tendon+slip, hand');
scatter(Angle3,Torque3,sz2,'+','CData',c3,'DisplayName','Tendon, LC');
scatter(Angle3,TorqueHand3,sz,'filled','CData',c3,'DisplayName','Tendon, hand');

legend
hold off


%% Plotting X axis Torque
figure
hold on
plot(phiD, TorqueR5(:, 1),phiD, TorqueR7(:, 1),'-.')
title('BPA X Torque, Length = 415 mm')
xlabel('Knee Extension(+)/Flexion(-), degrees')
ylabel('Torque, Nm')
legend('Ideal Tendon','Slippage+Tendon')
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