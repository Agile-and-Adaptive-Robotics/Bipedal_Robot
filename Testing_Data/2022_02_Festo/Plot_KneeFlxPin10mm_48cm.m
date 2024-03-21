%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

load KneeFlxPin_10mm_48cm.mat
Theoretical = TorqueR(:,3)';

restingLength = 0.485; %resting length, m
kmax = 0.398; %Length at maximum contraction, m

%% Test 1&2 done with CALT load cell. Tests 3-5 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet FlxTest10mm_1 from Results_table10mm_pinned_LoadCell
%Test 2 == sheet FlxTest10mm_2 from Results_table10mm_pinned_LoadCell
%Test 3 == sheet FlxTest10mm_1 from Results_table10mm_pinned_FishScale
%Test 4 == sheet FlxTest10mm_4 from Results_table10mm_pinned_FishScale
%Test 5 == sheet FlxTest10mm_2 from Results_table10mm_pinned_FishScale [not
%%used due to uncertainty about resting length]
%% Torque calculated from measurements
Angle{1} = [-5.5	-10.5	-25.5	-4.5	-14.5	-20	-33	-39.5	-43	-50	-62.5	-60	-55	-53.5	-46	-40	-36.5	-33.5	-29.5	-24.5	-17.5	-14.5	-5.5	-3	-0.5	3.5	11.5	19.5]';
L1 = length(Angle{1});
Angle{2} = [6	5.5	-1	-15	-26.6	-35	-42.5	-56	-56	-63	-64	-55	-48	-42	-35	-26	-19	-11	0.5	11.5]';
L2 = length(Angle{2});
Angle{3} = [-10	-14	-24	-27	-23	-21	-35	-45	-52	-61	-65]';
L3 = length(Angle{3});
Angle{4} = [-22	-25	-32	-35	-37	-42	-44.5	-52	-57	-62	-68	-71	-69.5	-67	-61.5	-48	-33	-34	-25.5	-36	-2	-8	-13.5	-19.5	-24	-24	-20	-19	-15	-13	-2]';
L4 = length(Angle{4});
% Angle{5} = [-18.5	-30.5	-39.5	-36.5	-42	-46.5	-52.5	-69.5	-68.5	-69.9	-60.5	-53	-47.5	-37.5	-44.5	-31.5	-26.5	-23]';
% L5 = length(Angle{5});

Torque{1} = [-13.55019944	-12.71574633	-12.2274503	-13.57429518	-12.61783318	-11.4271733	-7.828935983	-6.496239587	-4.940487337	-2.826134001	-1.210800042	-0.110215947	-1.513744872	-3.192740824	-4.861595823	-6.86945259	-8.599788984	-10.3582314	-11.8636379	-13.47755504	-13.77813076	-14.88459711	-14.53453947	-14.33814406	-14.02122399	-12.53904665	-11.16785477	-8.28360707]';
Torque{2} = [-12.02448648	-12.01072043	-12.46707433	-11.42926509	-9.361208447	-7.676179585	-4.957997672	-1.67225365	-0.807464392	-0.003574341	-2.019875	-4.170597119	-6.466679019	-8.166509485	-10.60330251	-12.64002428	-13.93045057	-14.68948292	-14.67482254	-13.15130291]';
Torque{3} = [-13.93505973	-13.41723648	-11.54052398	-9.443681562	-13.61674065	-13.66679506	-8.935879493	-5.788289674	-3.353532947	-1.833248329	-0.531376327]';
Torque{4} = [-15.24648357	-13.94170519	-12.36216232	-11.04703867	-10.0579279	-8.249812691	-6.962421539	-4.874859191	-3.096568251	-2.046884098	-1.0207759	0	0	-1.0207759	-2.834147213	-7.995951798	-11.5993906	-13.3993773	-15.9919036	-11.86241533	-16.23568514	-17.29899289	-17.10478208	-16.01872483	-14.22288584	-14.88045526	-15.17480498	-16.23203091	-16.25073881	-17.09230224	-16.15641241]';
% Torque{5} = [-19.06895474	-15.75701078	-16.75788099	-14.15214602	-10.75456087	-8.935879493	-5.252336926	-3.397585151	-1.308275082	-1.5822075	-3.390835418	-6.554396551	-8.661915638	-13.68558121	-16.00098846	-17.05641248	-18.14010372	-18.66468509]';
%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
InflatedLength{1} = [461.5	456	450.5	460.5	452	449	432	428	423	415	409.5	406	410	414.5	420.5	426	430	433.5	441	445	451	454	455.5	458.5	462	461.5	465	467]'/1000;
InflatedLength{2} = [467	466	463	450	443	431	423	411	408	404	407	413	420	422	433	442	448	455	460	466]'/1000;
InflatedLength{3} = [460	457	451	445	450	452	438	427	420	415	409]'/1000;
InflatedLength{4} = [450	445	439	435	431	426	421	410	408	406	402	398	398	400	404	420	429	434	444	431	461	459	453	450	442	442	443	445	449	454	458]'/1000;
% InflatedLength{5} = [459	448	449	433	423	419.5	410.5	406	396	392.5	401.5	407	411.5	426	432	436	441	442.5]'/1000;

ICRtoMuscle{1} = [62	74	76	64	72	76	79.5	81.5	77	81	74	68	74	77	75	79	80	78	79	77	75.5	72	61.5	62.5	54.5	48	38.5	30]'/1000;
ICRtoMuscle{2} = [45	50	55	65	70	74	73	72	65	63	66	71	75	75	76	75	69	65	56	44]'/1000;
ICRtoMuscle{3} = [66	70	75	77	74	75	78	77	75	68	69]'/1000;
ICRtoMuscle{4} = [75	77	80	83	84	84	84	79	77	75	69	65	65	67	73	84	83	82	79	83	58	64	69	73	78	78	75	75	75	70	62]'/1000;
% ICRtoMuscle{5} = [76	84	81	85	84	83	81	85	64	66.5	70	81	83.5	85.5	85	86	80.5	79.5]'/1000;

%load pressure where applicable
test = [1 2];
runsperseries = [28 19];

pres{1} = zeros(runsperseries(1),1);
pres{2} = zeros(runsperseries(2),1);
 
 for i = 1:size(test, 2)   
   if i == 1 
     for j = 1:runsperseries(i)
                file_name = sprintf('FlexTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres{1}(j) = Stats{'Mean',2};
     end
   else  
    for j = 1:runsperseries(i)
                file_name = sprintf('FlexTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres{2}(j) = Stats{'Mean',2};
    end
   end
 end
pres{2}(length(pres{2})+1,1) = mean(pres{2});                             %Missing 20th data point for Test 2 
pres{3} = 612*ones(length(InflatedLength{3}),1);
pres{4} = 612*ones(length(InflatedLength{4}),1);
% pres{5} = 612*ones(length(InflatedLength{5}),1);

KMAX = (restingLength-kmax)/restingLength;  %Converts to percentage
strainz = cell(length(Angle),1);
rel = cell(length(Angle),1);
F = cell(length(Angle),1);
TorqueHand = cell(length(Angle),1);

for i = 1:length(Angle)
    strainz{i} = ((restingLength-InflatedLength{i})/restingLength);
    rel{i} = strainz{i}/KMAX;
    F{i} = bpaForce10(restingLength,rel{i},pres{i});

    TorqueHand{i} = -ICRtoMuscle{i}.*F{i};  %Torque will be negative because it is causing flexion
end    

%% Plot setup
%Matlab hex color values:
c{1} = '#FFD700'; %gold
c{2} = '#FFB14E'; %orange
c{3} = '#FA8775'; %light orange
c{4} = '#EA5F94'; %pink
c{5} = '#CD34B5'; %magenta
c{6} = '#9D02D7'; %magenta 2
c{7} = '#0000FF'; %indigo
c{8} = '#000000'; %black
sz = 60;        %size of data points


%% Plot the expected value and scatter the data that show which test they come from
Test = ["FlxTest10mm_1...10mm_pinned_LoadCell"; 
        "FlxTest10mm_2...10mm_pinned_LoadCell";
        "FlxTest10mm_1...10mm_pinned_FishScale";
        "FlxTest10mm_4...10mm_pinned_FishScale";
        "FlxTest10mm_2...10mm_pinned_FishScale"];

%% Convert cells to column arrays once bad tests are eliminated
Angle = cell2mat(Angle');
Torque = cell2mat(Torque');
InflatedLength = cell2mat(InflatedLength');
ICRtoMuscle = cell2mat(ICRtoMuscle');
pres = cell2mat(pres');
strainz = cell2mat(strainz);
rel = cell2mat(rel);
F = cell2mat(F);
TorqueHand = cell2mat(TorqueHand);

%% Plot expected versus measured moment arm
Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque

figure
ax1 = gca;
hold on
pp = plot(phiD,G,'Color',c{7},'DisplayName','MA expected');
if ~iscell(Angle)
    ss = scatter(Angle, ICRtoMuscle,sz,'filled','MarkerFaceColor',c{6},'DisplayName','MA measured');
else
    for i = 1:length(Angle)
    ss{i} = scatter(Angle{i}, ICRtoMuscle{i},sz,'filled','MarkerFaceColor',c{7-i},'DisplayName',Test{i});
    end
end
hold off
title('Expected vs measured moment arm')
xlabel('Knee angle, degrees')
ylabel('Moment Arm, z axis (m)')
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 8;
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values
strain = Bifemsh_Pam.Contraction;
relstrain = (strain)./KMAX;

figure
ax1 = gca;
hold on
plot(phiD,relstrain,'DisplayName','Expected Relative Strain')
if ~iscell(Angle)
    scatter(Angle,rel,sz,'filled','MarkerFaceColor',c{6},'DisplayName','Measured Relative Strain')
else
    for i = 1:length(Angle)
        sc_rel{i} = scatter(Angle{i},rel{i},sz,'filled','MarkerFaceColor',c{7-i},'DisplayName',Test(i));
    end
end
hold off
title('Relative strain')
xlabel('Knee angle, \circ')
ylabel('strain/kmax')
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 8;
hold off

%% Plot measured versus expected strain (like above, but not normalized
figure
hold on
plot(phiD,strain,'Color',c{7},'DisplayName','Expected Strain')
if ~iscell(Angle)
    scatter(Angle,strainz,sz,'filled','MarkerFaceColor',c{6},'DisplayName','MeasuredStrain')
else
    for i = 1:length(Angle)
        sc_rel{i} = scatter(Angle{i},strainz{i},sz,'filled','MarkerFaceColor',c{7-i},'DisplayName',Test(i));
    end
end
hold off
title('Absolute strain')
xlabel('Knee angle, \circ')
ylabel('strain')
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 8;
hold off

%% Plot measured versus expected BPA length
MuscleLength = Bifemsh_Pam.MuscleLength;

figure
ax1 = gca;
hold on
plot(phiD,MuscleLength,'Color',c{7},'DisplayName','Expected Muscle Length')
if ~iscell(Angle)
    scatter(Angle,InflatedLength,sz,'filled','MarkerFaceColor',c{6},'DisplayName','Measured Length')
else
    for i = 1:length(Angle)
        sc_mL{i} = scatter(Angle{i},InflatedLength{i},sz,'filled','MarkerFaceColor',c{7-i},'DisplayName','Measured Length');
    end
end
hold off
title('Expected vs measured muscle length')
xlabel('Knee angle, \circ')
ylabel('Length, m')
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 8;
hold off

%% Plotting Torque 
fig2 = figure;
gcf2 = gcf;
ax2 = gca;
hold on
PL = plot(phiD, Theoretical,'Color',c{7},'Linewidth',2,'DisplayName','Expected Torque');
if ~iscell(Angle)
    scH = scatter(Angle,TorqueHand,sz,'MarkerFaceColor',c{2},'DisplayName','Hybrid calc');
    scM = scatter(Angle,Torque,sz,'filled','d','MarkerFaceColor',c{6},'DisplayName','Measured Torque');
else
    for i = 1:length(Angle)
    scM{i} = scatter(Angle{i},Torque{i},sz,'filled','d','MarkerFaceColor',c{7-i},'DisplayName',sprintf('Measured, test%d',i));
    scH{i} = scatter(Angle{i},TorqueHand{i},sz,'MarkerFaceColor',c{7-i},'DisplayName',sprintf('Back calc, test%d',i));
    end
end
hold off
title('Isometric Torque vs Knee Angle, 10mm Flexor, 48.5cm long')
xlabel('Knee angle, \circ')
ylabel('Torque, $N \cdot m$')
% set(gcf2,'Position',[1 384 950 612]);
set(ax2,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on');
lgd = legend;
lgd.FontSize = 8;


