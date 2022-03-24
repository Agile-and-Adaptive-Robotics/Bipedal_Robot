%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

restingLength = 0.485; %resting length, m
kmax = 0.398; %Length at maximum contraction, m

load KneeFlxPin_10mm_48cm.mat
Theoretical = TorqueR(:,3)';
%% Test 1&2 done with CALT load cell. Tests 3-5 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet FlxTest10mm_1 from Results_table10mm_pinned_LoadCell
%Test 2 == sheet FlxTest10mm_2 from Results_table10mm_pinned_LoadCell
%Test 3 == sheet FlxTest10mm_1 from Results_table10mm_pinned_FishScale
%Test 4 == sheet FlxTest10mm_4 from Results_table10mm_pinned_FishScale
%Test 5 == sheet FlxTest10mm_2 from Results_table10mm_pinned_FishScale [not
%%used due to uncertainty about resting length]
%% Torque calculated from measurements

Angle1 = [-5.5	-10.5	-25.5	-4.5	-14.5	-20	-33	-39.5	-43	-50	-62.5	-60	-55	-53.5	-46	-40	-36.5	-33.5	-29.5	-24.5	-17.5	-14.5	-5.5	-3	-0.5	3.5	11.5	19.5];
Angle2 = [6	5.5	-1	-15	-26.6	-35	-42.5	-56	-56	-63	-64	-55	-48	-42	-35	-26	-19	-11	0.5	11.5];
Angle3 = [-10	-14	-24	-27	-23	-21	-35	-45	-52	-61	-65];
Angle4 = [-22	-25	-32	-35	-37	-42	-44.5	-52	-57	-62	-68	-71	-69.5	-67	-61.5	-48	-33	-34	-25.5	-36	-2	-8	-13.5	-19.5	-24	-24	-20	-19	-15	-13	-2];
%Angle5 = [-18.5	-30.5	-39.5	-36.5	-42	-46.5	-52.5	-69.5	-68.5	-69.9	-60.5	-53	-47.5	-37.5	-44.5	-31.5	-26.5	-23];
Angle = [Angle1, Angle2, Angle3, Angle4];

Torque1 = [-13.55019944	-12.71574633	-12.2274503	-13.57429518	-12.61783318	-11.4271733	-7.828935983	-6.496239587	-4.940487337	-2.826134001	-1.210800042	-0.110215947	-1.513744872	-3.192740824	-4.861595823	-6.86945259	-8.599788984	-10.3582314	-11.8636379	-13.47755504	-13.77813076	-14.88459711	-14.53453947	-14.33814406	-14.02122399	-12.53904665	-11.16785477	-8.28360707];
Torque2 = [-12.02448648	-12.01072043	-12.46707433	-11.42926509	-9.361208447	-7.676179585	-4.957997672	-1.67225365	-0.807464392	-0.003574341	-2.019875	-4.170597119	-6.466679019	-8.166509485	-10.60330251	-12.64002428	-13.93045057	-14.68948292	-14.67482254	-13.15130291];
Torque3 = [-13.93505973	-13.41723648	-11.54052398	-9.443681562	-13.61674065	-13.66679506	-8.935879493	-5.788289674	-3.353532947	-1.833248329	-0.531376327];
Torque4 = [-15.24648357	-13.94170519	-12.36216232	-11.04703867	-10.0579279	-8.249812691	-6.962421539	-4.874859191	-3.096568251	-2.046884098	-1.0207759	0	0	-1.0207759	-2.834147213	-7.995951798	-11.5993906	-13.3993773	-15.9919036	-11.86241533	-16.23568514	-17.29899289	-17.10478208	-16.01872483	-14.22288584	-14.88045526	-15.17480498	-16.23203091	-16.25073881	-17.09230224	-16.15641241];
%Torque5 = [-19.06895474	-15.75701078	-16.75788099	-14.15214602	-10.75456087	-8.935879493	-5.252336926	-3.397585151	-1.308275082	-1.5822075	-3.390835418	-6.554396551	-8.661915638	-13.68558121	-16.00098846	-17.05641248	-18.14010372	-18.66468509];
Torque = [Torque1, Torque2, Torque3, Torque4];

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength1 = [461.5	456	450.5	460.5	452	449	432	428	423	415	409.5	406	410	414.5	420.5	426	430	433.5	441	445	451	454	455.5	458.5	462	461.5	465	467]/1000;
InflatedLength2 = [467	466	463	450	443	431	423	411	408	404	407	413	420	422	433	442	448	455	460	466]/1000;
InflatedLength3 = [460	457	451	445	450	452	438	427	420	415	409]/1000;
InflatedLength4 = [450	445	439	435	431	426	421	410	408	406	402	398	398	400	404	420	429	434	444	431	461	459	453	450	442	442	443	445	449	454	458]/1000;
%InflatedLength5 = [459	448	449	433	423	419.5	410.5	406	396	392.5	401.5	407	411.5	426	432	436	441	442.5]/1000;
InflatedLength = [InflatedLength1, InflatedLength2, InflatedLength3, InflatedLength4];

ICRtoMuscle1 = [62	74	76	64	72	76	79.5	81.5	77	81	74	68	74	77	75	79	80	78	79	77	75.5	72	61.5	62.5	54.5	48	38.5	30]/1000;
ICRtoMuscle2 = [45	50	55	65	70	74	73	72	65	63	66	71	75	75	76	75	69	65	56	44]/1000;
ICRtoMuscle3 = [66	70	75	77	74	75	78	77	75	68	69]/1000;
ICRtoMuscle4 = [75	77	80	83	84	84	84	79	77	75	69	65	65	67	73	84	83	82	79	83	58	64	69	73	78	78	75	75	75	70	62]/1000;
%ICRtoMuscle5 = [76	84	81	85	84	83	81	85	64	66.5	70	81	83.5	85.5	85	86	80.5	79.5]/1000;
ICRtoMuscle = [ICRtoMuscle1, ICRtoMuscle2, ICRtoMuscle3, ICRtoMuscle4];

F1 = zeros(1,size(InflatedLength1, 2));
F2 = zeros(1,size(InflatedLength2, 2));
F3 = zeros(1,size(InflatedLength3, 2));
F4 = zeros(1,size(InflatedLength4, 2));
%F5 = zeros(1,size(InflatedLength5, 2));
F = [F1, F2, F3, F4];

TorqueHand1 = zeros(1,size(InflatedLength1, 2));
TorqueHand2 = zeros(1,size(InflatedLength2, 2));
TorqueHand3 = zeros(1,size(InflatedLength3, 2));
TorqueHand4 = zeros(1,size(InflatedLength4, 2));
%TorqueHand5 = zeros(1,size(InflatedLength5, 2));
TorqueHand = [TorqueHand1, TorqueHand2, TorqueHand3, TorqueHand4];

%load pressure where applicable
test = [1 2];
runsperseries = [28 19];

pres1 = zeros(1,runsperseries(1));
pres2 = zeros(1,runsperseries(2));
 
 for i = 1:size(test, 2)   
   if i == 1 
     for j = 1:runsperseries(i)
                file_name = sprintf('FlexTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres1(1,j) = Stats{'Mean',2};
     end
   else  
    for j = 1:runsperseries(i)
                file_name = sprintf('FlexTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres2(1,j) = Stats{'Mean',2};
    end
   end
 end
presx = mean(pres2);                             %Missing 20th data point for Test 2 
pres3 = 612*ones(1,size(InflatedLength3, 2));
pres4 = 612*ones(1,size(InflatedLength4, 2));
%pres5 = 612*ones(1,size(InflatedLength5, 2));
pres = [pres1 pres2 presx pres3 pres4];

for i = 1:size(InflatedLength, 2)
    F(i) = festo3(InflatedLength(i), restingLength, 10, pres(i), kmax);    
    TorqueHand(i) = -ICRtoMuscle(i)*F(i);  %Moment will be negative because it causes flexion
end
TorqueHand1 = TorqueHand(1:size(TorqueHand1,2));
TorqueHand2 = TorqueHand((size(TorqueHand1,2)+1):(size(TorqueHand1,2)+size(TorqueHand2,2)));
TorqueHand3 = TorqueHand((size(TorqueHand1,2)+size(TorqueHand2,2)+1):(size(TorqueHand1,2)+size(TorqueHand2,2)+size(TorqueHand3,2)));
TorqueHand4 = TorqueHand((size(TorqueHand1,2)+size(TorqueHand2,2)+size(TorqueHand3,2)+1):(size(TorqueHand1,2)+size(TorqueHand2,2)+size(TorqueHand3,2)+size(TorqueHand4,2)));
%TorqueHand5 = TorqueHand((size(TorqueHand1,2)+size(TorqueHand2,2)+size(TorqueHand3,2)+size(TorqueHand4,2)+1):size(TorqueHand,2));
%% Mean and RMSE
X = linspace(min(Angle),max(Angle),size(Angle,2));      %Range of motion
mod = 'gauss1';
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
scatter(Angle1,Torque1,sz,'d','CData',c1,'DisplayName','JM LC');
scatter(Angle2,Torque2,sz,'d','CData',c3,'DisplayName','BB LC');
scatter(Angle3,Torque3,sz,'d','CData',c4,'DisplayName','BB FS');
scatter(Angle4,Torque4,sz,'d','CData',c5,'DisplayName','BB FS');
%scatter(Angle5,Torque5,sz,'d','CData',c2,'DisplayName','JM FS');
scatter(Angle1,TorqueHand1,sz,'filled','CData',c1,'DisplayName','JM hand');
scatter(Angle2,TorqueHand2,sz,'filled','CData',c3,'DisplayName','BB hand');
scatter(Angle3,TorqueHand3,sz,'filled','CData',c4,'DisplayName','BB hand');
scatter(Angle4,TorqueHand4,sz,'filled','CData',c5,'DisplayName','BB hand');
%scatter(Angle5,TorqueHand5,sz,'filled','CData',c2,'DisplayName','JM hand');


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

scatter(Angle1,Torque1,sz,'d','CData',c1,'DisplayName','JM LC');
scatter(Angle2,Torque2,sz,'d','CData',c3,'DisplayName','BB LC');
scatter(Angle3,Torque3,sz,'d','CData',c4,'DisplayName','BB FS');
scatter(Angle4,Torque4,sz,'d','CData',c5,'DisplayName','BB FS');
%scatter(Angle5,Torque5,sz,'d','CData',c2,'DisplayName','JM FS');
scatter(Angle1,TorqueHand1,sz,'filled','CData',c1,'DisplayName','JM hand');
scatter(Angle2,TorqueHand2,sz,'filled','CData',c3,'DisplayName','BB hand');
scatter(Angle3,TorqueHand3,sz,'filled','CData',c4,'DisplayName','BB hand');
scatter(Angle4,TorqueHand4,sz,'filled','CData',c5,'DisplayName','BB hand');
%scatter(Angle5,TorqueHand5,sz,'filled','CData',c2,'DisplayName','JM hand');

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