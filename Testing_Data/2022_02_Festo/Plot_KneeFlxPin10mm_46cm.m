%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

restingLength = 0.457; %resting length, m
kmax = 0.380; %Length at maximum contraction, m

load KneeFlxPin_10mm_46cm.mat
Theoretical = TorqueR(:,3)';
%% Tests 1 & 2 . 
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
X1 = linspace(min(Angle),max(Angle),size(Angle,2));      %Range of motion
mod = 'gauss1';
fitOptions = fitoptions(mod, 'Normalize', 'on');
[mdl1u, gof1] = fit(Angle',Torque',mod,fitOptions);
TorqueStdu = gof1.rmse;
TorqueMeanu = feval(mdl1u,X1)';

[mdl2u, gof2] = fit(Angle',TorqueHand',mod,fitOptions);
HandStdu = gof2.rmse;
HandMeanu = feval(mdl2u,X1)';

modp = 'poly3';
fitOp = fitoptions(modp,'Normalize','on');
[mdl1, gofp1] = fit(Angle',Torque',modp,fitOptions)
TorqueStd = gofp1.rmse
TorqueMean = feval(mdl1,X1)';

[mdl2, gofp2] = fit(Angle',TorqueHand',modp,fitOptions)
HandStd = gofp2.rmse
HandMean = feval(mdl2,X1)';

%% Plotting with polynomial solver
%Matlab hex color values:
c1 = '#FFD700'; %gold
c2 = '#FFB14E'; %orange
c3 = '#FA8775'; %light orange
c4 = '#EA5F94'; %pink
c5 = '#CD34B5'; %magenta
c6 = '#9D02D7'; %magenta 2
c7 = '#0000FF'; %indigo
sz = 60;        %size of data points


figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Flexor, 45.5cm long')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
gca1 = gca;
gcf1 = gcf;
% set(gcf,'Position',[1 384 950 612]);
% set(gca,'FontSize', 18, 'FontWeight', 'bold','XMinorGrid','on','XMinorTick','on','YMinorGrid','on','YMinorTick','on');
set(gca,'FontSize', 12, 'FontWeight', 'bold')
plot(phiD, Theoretical,'Color',c5,'Linewidth',2,'DisplayName','Expected')
% chr = ['Adjusted fit'];
% plot(phiD, Theo_adj,'Color',c3,'Linewidth',2,'DisplayName',chr)

plot(X1,TorqueMean,'--','Color',c7,'Linewidth',2,'DisplayName','Torque mean, scale')
plot(X1,HandMean,'--','Color',c1,'Linewidth',2,'DisplayName','Torque expected, hand measurements')

sc1 = scatter(Angle,Torque,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Torque data, experimental');
sc2 = scatter(Angle,TorqueHand,sz,'filled','MarkerFaceColor',c1,'DisplayName','Torque data, Hand calc');

legend
hold off


%% Plot error and standard deviation as bar graphs
xb=categorical({'Theoretical','Scale','Hand'});
xb = reordercats(xb,{'Theoretical','Scale','Hand'});
yb = [mean(Theoretical,'omitnan') mean(TorqueMean) mean(TorqueHand)];
std_dev = [0 mean(TorqueStd) mean(HandStd)];
figure
hold on
b = bar(xb,yb,'FaceColor',[0 0.8 1.0]);
gca3 = gca;
gcf3 = gcf;
ylabel('Torque, N*m')
title('Mean Torque Values and SD for Each Calculation Method')
% set(gcf,'Position',[0 0 950 612]);
set(gca,'FontSize', 12, 'FontWeight', 'bold');
b.CData = [0 0.4470 0.7410; 0  0  0; 1  0  0];
errb = errorbar(yb,std_dev ,'LineStyle','none','LineWidth',4,'CapSize',20);
hold off