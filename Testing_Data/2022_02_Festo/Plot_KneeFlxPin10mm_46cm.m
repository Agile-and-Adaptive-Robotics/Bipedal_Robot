%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

restingLength = 0.457; %resting length, m
kmax = 0.380; %Length at maximum contraction, m

load KneeFlxPin_10mm_46cm.mat
Theoretical = TorqueR(:,3)';
%% Tests 1 & 2 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet FlxTest10mm_2 from Results_table10mm_pinned_FishScale
%Test == sheet FlxTest10mm_3 from Results_table10mm_pinned_FishScale
%% Torque calculated from measurements

% Angle1 = [-18.5	-30.5	-39.5	-36.5	-42	-46.5	-52.5	-69.5	-68.5	-69.9	-60.5	-53	-47.5	-37.5	-44.5	-31.5	-26.5	-23];
Angle = [-30	-32	-36.5	-41	-51	-52.5	-56.5	-63	-66.5	-74.5	-80.5	-80	-82	-86	-90	-85	-81.5	-76	-70	-62	-58	-49	-45	-37	-33.5	-26.5];
% Angle = [Angle1, Angle2];

% Torque1 = [-19.06895474	-15.75701078	-16.75788099	-14.15214602	-10.75456087	-8.935879493	-5.252336926	-3.397585151	-1.308275082	-1.5822075	-3.390835418	-6.554396551	-8.661915638	-13.68558121	-16.00098846	-17.05641248	-18.14010372	-18.66468509];
Torque = [-20.00095306	-19.62412624	-17.74104856	-14.64175686	-13.34974574	-11.21468685	-8.847415849	-6.225462211	-5.473176171	-4.159138583	-1.821998215	-1.821998215	-1.293882791	-0.526243409	0	-1.293882791	-2.098936493	-3.894225298	-5.473176171	-7.788450595	-11.02208726	-14.67316852	-18.04885624	-21.08709753	-22.77466623	-19.33565341];
% Torque = [Torque1, Torque2];

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

% InflatedLength1 = [459	448	449	433	423	419.5	410.5	406	396	392.5	401.5	407	411.5	426	432	436	441	442.5]/1000;
InflatedLength = [438	438	434	424	419	415	408	401	399	394	388	388	386	381	380	384	386	390	394	400	407	413	422	426	433	438]/1000;
% InflatedLength = [InflatedLength1, InflatedLength2];

% ICRtoMuscle1 = [76	84	81	85	84	83	81	85	64	66.5	70	81	83.5	85.5	85	86	80.5	79.5]/1000;
ICRtoMuscle = [85	85	85	88	90	88	83	80	75	70	60	64	60	60	63	65	65	68	75	79	84	90	92	90	89	85]/1000;
% ICRtoMuscle = [ICRtoMuscle1, ICRtoMuscle2];

% F1 = zeros(1,size(InflatedLength1, 2));
F = zeros(1,size(InflatedLength, 2));
% F = [F1, F2];

% TorqueHand1 = zeros(1,size(InflatedLength1, 2));
TorqueHand = zeros(1,size(InflatedLength, 2));
% TorqueHand = [TorqueHand1, TorqueHand2];

%load pressure where applicable
% pres1 = 612*ones(1,size(InflatedLength1, 2));
pres = 612*ones(1,size(InflatedLength, 2));
% pres = [pres1 pres2];

for i = 1:size(InflatedLength, 2)
    F(i) = festo3(InflatedLength(i), restingLength, 10, pres(i), kmax);    
    TorqueHand(i) = -ICRtoMuscle(i)*F(i);  %Torque will be negative because it is causing flexion
end
% TorqueHand1 = TorqueHand(1:size(TorqueHand1,2));
% TorqueHand2 = TorqueHand((size(TorqueHand1,2)+1):size(TorqueHand,2));
%% Mean and RMSE
X = linspace(-120,10,size(Angle,2));      %Range of motion
mod = 'gauss1';
fitOptions = fitoptions(mod);
[mdl1u, gof1] = fit(Angle',Torque',mod,fitOptions)
TorqueStdu = gof1.rmse;
TorqueMeanu = feval(mdl1u,X)';

% [mdl1,S1,mu1] = polyfit(Angle,Torque,3);
% [TorqueMean, TorqueStd] = polyval(mdl1,X,S1,mu1);


[mdl2u, gof2] = fit(Angle',TorqueHand',mod,fitOptions)
HandStdu = gof2.rmse
HandMeanu = feval(mdl2u,X)';
% [mdl2,S2,mu2] = polyfit(Angle,TorqueHand,3);
% [HandMean,HandStd] = polyval(mdl2,X,S2,mu2);

%% Plotting with polynomial solver
% figure
% hold on
% title('Isometric Torque vs Knee Angle, 10mm Extensor, 45.5cm long')
% xlabel('degrees Flexion(-),Extension(+)')
% ylabel('Torque, N*m')
% plot(phiD, Theoretical,'Color',[0 0.4470 0.7410],'Linewidth',2,'DisplayName','Theoretical Calculation')
% 
% Xnew=[X,fliplr(X)];
% Y1=[TorqueMean+TorqueStd,fliplr(TorqueMean-TorqueStd)];
% Y2=[HandMean+HandStd,fliplr(HandMean-HandStd)];
% fill(Xnew,Y1,[1 0.4 0.8],'DisplayName','Fish scale std','FaceAlpha',0.25);
% fill(Xnew,Y2,[.6 1.0 .6],'DisplayName','Hand torque std','FaceAlpha',0.25);
% plot(X,TorqueMean,'--k','Linewidth',2,'DisplayName','Torque mean, scale')
% plot(X,HandMean,'--r','Linewidth',2,'DisplayName','Torque mean, hand')
% 
% c = [0.8500 0.3250 0.0980]; % color
% sz = 50;
% c1 = [0.8500 0.3250 0.0980]; % color, burnt orange
% c2 = [0.6350 0.0780 0.1840]; %color, red/violet
% c3 = [0 0.4470 0.7410]; %color, navy blue
% c4 = [0.4660 0.6740 0.1880]; %color, moss green
% c5 = [0.9290 0.6940 0.1250]; %color, dark yellow
% scatter(Angle1,Torque1,sz,'d','CData',c1,'DisplayName','JM FS');
% scatter(Angle2,Torque2,sz,'d','CData',c3,'DisplayName','BB FS');
% scatter(Angle1,TorqueHand1,sz,'filled','CData',c1,'DisplayName','JM hand');
% scatter(Angle2,TorqueHand2,sz,'filled','CData',c3,'DisplayName','BB hand');
% 
% legend
% hold off

%% Plotting with nonlinear solver
figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Extensor, 45.5cm long')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
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
% scatter(Angle1,Torque1,sz,'d','CData',c1,'DisplayName','JM FS');
scatter(Angle,Torque,sz,'d','CData',c3,'DisplayName','BB FS');
% scatter(Angle1,TorqueHand1,sz,'filled','CData',c1,'DisplayName','JM hand');
scatter(Angle,TorqueHand,sz,'filled','CData',c3,'DisplayName','BB hand');

legend
hold off