%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

restingLength = 0.457; %resting length, m
kmax = 0.380; %Length at maximum contraction, m

load KneeFlxPin_10mm_46cm.mat
Theoretical = TorqueR(:,3)';
%% Test 1 done with CALT load cell. Tests 2 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet FlxTest10mm_1 from Results_table10mm_pinned_LoadCell
%Test 2 == sheet FlxTest10mm_2 from Results_table10mm_pinned_LoadCell
%Test 3 == sheet FlxTest10mm_1 from Results_table10mm_pinned_FishScale
%Test 4 == sheet FlxTest10mm_4 from Results_table10mm_pinned_FishScale
%% Torque calculated from measurements

Angle1 = [];
Angle2 = [];
Angle3 = [];
Angle4 = [];
Angle = [Angle1, Angle2, Angle3, Angle4];

Torque1 = [];
Torque2 = [5.969270431	5.781894229	4.978853364	3.667219951	2.623266826	2.087906249	1.579313702	1.311633413	0];
Torque = [Torque1, Torque2];

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength1 = [453.5	452	456	450.5	449.5	442.5	439	435	422.5	420	410.5	398.5	400.5	395	384	381	382	394	408.5	414.5	420.5	425	427	432	438	438	441	445	447]/1000;
InflatedLength2 = [430	430	420	415	410	401	390	388	385]/1000;
InflatedLength = [InflatedLength1, InflatedLength2];

ICRtoMuscle1 = [29.5	30	30	29.5	30.5	30	30	32.5	34	35	37	38.5	44.5	50	54	60	55	48	34.5	33	32.5	31.5	30	30	30	30	30	30	30]/1000;
ICRtoMuscle2 = [30	30	30	30	34	35	45	50	55]/1000;
ICRtoMuscle = [ICRtoMuscle1, ICRtoMuscle2];

F1 = zeros(1,size(InflatedLength1, 2));
F2 = zeros(1,size(InflatedLength2, 2));
F = [F1, F2];

TorqueHand1 = zeros(1,size(InflatedLength1, 2));
TorqueHand2 = zeros(1,size(InflatedLength2, 2));
TorqueHand = [TorqueHand1, TorqueHand2];

%load pressure where applicable
test = 2;
runsperseries = 29;

    pres1 = zeros(1,runsperseries);
    
    for j = 1:runsperseries
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test,j);
                load(file_name,'Stats')
                pres1(1,j) = Stats{'Mean',2};
    end

pres2 = 612*ones(1,size(InflatedLength2, 2));
pres = [pres1 pres2];

for i = 1:size(InflatedLength, 2)
    F(i) = festo3(InflatedLength(i), restingLength, 10, pres(i), 0.380);    
    TorqueHand(i) = ICRtoMuscle(i)*F(i);
end
TorqueHand1 = TorqueHand(1:size(TorqueHand1,2));
TorqueHand2 = TorqueHand(((size(TorqueHand1,2)+1)):size(TorqueHand,2));
%% Mean and RMSE
X1 = linspace(-120,10,size(Angle,2));      %Range of motion
% modelfun = @(b,x)b(1)*cosd(b(2)*x+b(3)) + b(4)*sind(b(5)*x+b(6))+b(7);
% beta0 = [0 1 1 6 1 90 6];
% mdl1 = fitnlm(Angle,Torque,modelfun,beta0)
% TorqueStd = mdl1.RMSE;
% TorqueMean = feval(mdl1,X);

[mdl1,S1,mu1] = polyfit(Angle,Torque,3);
[TorqueMean, TorqueStd] = polyval(mdl1,X1,S1,mu1);

% mdl2 = fitnlm(Angle,TorqueHand,modelfun,beta0)
% HandStd = mdl2.RMSE;
% HandMean = feval(mdl2,X);
X2 = X1;
[mdl2,S2,mu2] = polyfit(Angle,TorqueHand,3);
[HandMean,HandStd] = polyval(mdl2,X2,S2,mu2);

%% Plotting
figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Extensor, 41.5cm long')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
plot(phiD, Theoretical,'Color',[0 0.4470 0.7410],'Linewidth',2,'DisplayName','Theoretical Calculation')

X1new=[X1,fliplr(X1)];
Y1=[TorqueMean+TorqueStd,fliplr(TorqueMean-TorqueStd)];
X2new=[X2,fliplr(X2)];
Y2=[HandMean+HandStd,fliplr(HandMean-HandStd)];
fill(X1new,Y1,[1 0.4 0.8],'DisplayName','Fish scale std','FaceAlpha',0.25);
fill(X2new,Y2,[.6 1.0 .6],'DisplayName','Hand torque std','FaceAlpha',0.25);
plot(X1,TorqueMean,'--k','Linewidth',2,'DisplayName','Torque mean, scale')
plot(X2,HandMean,'--r','Linewidth',2,'DisplayName','Torque mean, hand')

sz = 50;
c = [0.8500 0.3250 0.0980]; % color
scatter(Angle1,Torque1,sz,'d','g','DisplayName','BB&JM LC');
ss = scatter(Angle2,Torque2,sz,'d','CData',c,'DisplayName','BB FS');
scatter(Angle1,TorqueHand1,[],'g','filled','DisplayName','BB&JM hand');
sb = scatter(Angle2,TorqueHand2,sz,'filled','CData',c,'DisplayName','BB hand');


legend
hold off