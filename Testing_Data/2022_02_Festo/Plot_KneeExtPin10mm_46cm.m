%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

restingLength = 0.457; %resting length, m
kmax = 0.380; %Length at maximum contraction, m

load KneeExtPin_10mm_46cm.mat
Theoretical = TorqueR(:,3)';
%% Test 1 done with CALT load cell. Tests 2 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet ExtTest10mm_2 from Results_table10mm_pinned_LoadCell
%Test 2 == sheet ExtTest10mm_1 from Results_table10mm_pinned_FishScale
%% Torque calculated from measurements

Angle1 = [-125	-120.5	-115	-111.5	-108	-101.5	-88	-80	-71	-58	-40	-32	-17.5	-10.5	-4	2	-3.5	-16	-33.5	-48	-58.5	-72	-80	-85	-88.5	-99.5	-111.5	-114.5	-123];
Angle2 = [-102	-100	-81.5	-66	-47	-35	-17	-11	-5];
Angle = [Angle1, Angle2];

Torque1 = [6.168692258	6.072984004	5.91673494	5.983690138	6.220291194	5.580594946	4.981638504	4.559983227	4.059476693	3.15664206	2.56531606	2.127108003	1.785781666	0.774715989	0.309696885	-0.057265164	0.536074229	2.318668963	2.990337108	3.61356103	3.999580201	4.737782449	5.291362695	5.556525636	6.028604398	6.318428803	6.41672146	6.668472107	6.512037525];
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

pres2 = 606*ones(1,size(InflatedLength2, 2));
pres = [pres1 pres2];

for i = 1:size(InflatedLength, 2)
    F(i) = festo3(InflatedLength(i), restingLength, 10, pres(i), 0.380);    
    TorqueHand(i) = ICRtoMuscle(i)*F(i);
end
TorqueHand1 = TorqueHand(1:size(TorqueHand1,2));
TorqueHand2 = TorqueHand(((size(TorqueHand1,2)+1)):size(TorqueHand,2));
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

%% Plotting polynomial solver
figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Extensor, 45.5cm long')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
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
scatter(Angle1,Torque1,sz,'d','g','DisplayName','BB&JM LC');
scatter(Angle2,Torque2,sz,'d','CData',c,'DisplayName','BB FS');
scatter(Angle1,TorqueHand1,sz,'g','filled','DisplayName','BB&JM hand');
scatter(Angle2,TorqueHand2,sz,'filled','CData',c,'DisplayName','BB hand');


legend
hold off

%% Plotting nonlinear solver
figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Extensor, 45.5cm long')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
plot(phiD, Theoretical,'Color',[0 0.4470 0.7410],'Linewidth',2,'DisplayName','Theoretical Calculation')

Xnew=[X,fliplr(X)];
Y3=[TorqueMeanu+TorqueStdu,fliplr(TorqueMeanu-TorqueStdu)];
Y4=[HandMeanu+HandStdu,fliplr(HandMeanu-HandStdu)];
plot(X,TorqueMeanu,'--k','Linewidth',2,'DisplayName','Torque mean, scale')
fill(Xnew,Y3,[1 0.4 0.8],'DisplayName','Scale SD','FaceAlpha',0.25);
plot(X,HandMeanu,'--r','Linewidth',2,'DisplayName','Torque mean, hand')
fill(Xnew,Y4,[.6 1.0 .6],'DisplayName','Hand torque SD','FaceAlpha',0.25);

scatter(Angle1,Torque1,sz,'d','g','DisplayName','BB&JM LC');
scatter(Angle2,Torque2,sz,'d','CData',c,'DisplayName','BB FS');
scatter(Angle1,TorqueHand1,sz,'g','filled','DisplayName','BB&JM hand');
scatter(Angle2,TorqueHand2,sz,'filled','CData',c,'DisplayName','BB hand');


legend
hold off