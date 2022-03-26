%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

restingLength = 0.415; %resting length, m
kmax = 0.349; %Length at maximum contraction, m

load KneeExtPin_10mm_42cm.mat
Theoretical = TorqueR(:,3)';
%% Test 1 done with CALT load cell. Tests 2 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet ExtTest10mm_4 from Results_table10mm_pinned_LoadCell
%% Torque calculated from measurements

Angle = [-76.5	-69	-51.5	-40.5	-39.5	-55.5	-35	-25.5	-20	-12.5	-5	0.5	4	12.5	20	28	35	23	20	11	1	-16	-24	-44	-53.4	-66];
% Angle2 = [-102	-100	-81.5	-66	-47	-35	-17	-11	-5];
% Angle = [Angle1, Angle2];

Torque = [8.50477362	7.539449285	6.619468483	5.835581619	6.120370289	7.005900186	6.178291184	6.391414397	5.96804559	6.099105123	5.93740855	4.880958611	4.375659895	3.388083398	2.795862193	1.758237576	0.306753172	1.725902869	3.204095044	5.145222736	6.304242352	6.817870142	7.511919159	7.911648561	8.110527045	8.862066797];
% Torque2 = [5.969270431	5.781894229	4.978853364	3.667219951	2.623266826	2.087906249	1.579313702	1.311633413	0];
% Torque = [Torque1, Torque2];

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength = [413.5	405.5	402	400	403	412	401	398	392.5	385	385.5	371.5	365	361	356	355	346	348	354	365	376	398	395	400	407	413]/1000;
% InflatedLength2 = [430	430	420	415	410	401	390	388	385]/1000;
% InflatedLength = [InflatedLength1, InflatedLength2];

ICRtoMuscle = [28	30.5	31.5	31	32	31	35	32	41	44.5	45.5	55.5	66.5	62	65	68.5	71	67	60	55	55	42	38	31	29	28]/1000;
% ICRtoMuscle2 = [30	30	30	30	34	35	45	50	55]/1000;
% ICRtoMuscle = [ICRtoMuscle1, ICRtoMuscle2];

F = zeros(1,size(InflatedLength, 2));
% F2 = zeros(1,size(InflatedLength2, 2));
% F = [F1, F2];

TorqueHand = zeros(1,size(InflatedLength, 2));
% TorqueHand2 = zeros(1,size(InflatedLength2, 2));
% TorqueHand = [TorqueHand1, TorqueHand2];

%load pressure where applicable
test = 4;
runsperseries = 26;

    pres = zeros(1,runsperseries);
    
    for j = 1:runsperseries
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test,j);
                load(file_name,'Stats')
                pres(1,j) = Stats{'Mean',2};
    end

% pres2 = 612*ones(1,size(InflatedLength2, 2));
% pres = [pres1 pres2];

for i = 1:size(InflatedLength, 2)
    F(i) = festo3(InflatedLength(i), restingLength, 10, pres(i), kmax);    
    TorqueHand(i) = ICRtoMuscle(i)*F(i);
end
% TorqueHand1 = TorqueHand(1:size(TorqueHand1,2));
% TorqueHand2 = TorqueHand(((size(TorqueHand1,2)+1)):size(TorqueHand,2));

%% Mean and RMSE
X = linspace(min(Angle),max(Angle),size(Angle,2));      %Range of motion
mod = 'sin2';
fitOptions = fitoptions(mod, 'Normalize', 'on');
[mdl1u, gof1] = fit(Angle',Torque',mod,fitOptions)
TorqueStdu = gof1.rmse
TorqueMeanu = feval(mdl1u,X)';

[mdl2u, gof2] = fit(Angle',TorqueHand',mod,fitOptions);
HandStdu = gof2.rmse;
HandMeanu = feval(mdl2u,X)';

modp = 'poly3';
fitOp = fitoptions(modp,'Normalize','on');
[mdl1, gofp1] = fit(Angle',Torque',modp,fitOp)
TorqueStd = gofp1.rmse
TorqueMean = feval(mdl1,X)';

[mdl2, gofp2] = fit(Angle',TorqueHand',modp,fitOp);
HandStd = gofp2.rmse;
HandMean = feval(mdl2,X)';

%% Plotting polynomial fit
figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Extensor, 41.5cm long')
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
fill(Xnew,Y1,[1 0.4 0.8],'DisplayName','Scale SD','FaceAlpha',0.25);
plot(X,HandMean,'--r','Linewidth',2,'DisplayName','Torque mean, hand')
fill(Xnew,Y2,[.6 1.0 .6],'DisplayName','Hand torque SD','FaceAlpha',0.25);

sz = 50;
c = [0.8500 0.3250 0.0980]; % color
scatter(Angle,Torque,sz,'d','CData',c,'DisplayName','BB&JM LC');
scatter(Angle,TorqueHand,sz,'filled','CData',c,'DisplayName','BB&JM hand');

legend
hold off
%% Plotting nonlinear fit
figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Extensor, 41.5cm long')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
gca2 = gca;
gcf2 = gcf;
set(gcf,'Position',[960 384 950 612]);
set(gca,'FontSize', 18, 'FontWeight', 'bold','XMinorGrid','on','XMinorTick','on','YMinorGrid','on','YMinorTick','on');
plot(phiD, Theoretical,'Color',[0 0.4470 0.7410],'Linewidth',2,'DisplayName','Theoretical Calculation')

Y3=[TorqueMeanu+TorqueStdu,fliplr(TorqueMeanu-TorqueStdu)];
Y4=[HandMeanu+HandStdu,fliplr(HandMeanu-HandStdu)];
plot(X,TorqueMeanu,'--k','Linewidth',2,'DisplayName','Torque mean, scale')
fill(Xnew,Y3,[1 0.4 0.8],'DisplayName','Scale SD','FaceAlpha',0.25);
plot(X,HandMeanu,'--r','Linewidth',2,'DisplayName','Torque mean, hand')
fill(Xnew,Y4,[.6 1.0 .6],'DisplayName','Hand torque SD','FaceAlpha',0.25);

sz = 50;
c = [0.8500 0.3250 0.0980]; % color
scatter(Angle,Torque,sz,'d','CData',c,'DisplayName','BB&JM LC');
scatter(Angle,TorqueHand,sz,'filled','CData',c,'DisplayName','BB&JM hand');

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