%% Pinned knee, Extensor
%Run and save data from testing results
clear;
clc;
close all;

load KneeExtPin_10mm_all.mat
restingLength = 0.480;      %resting length, m
kmax = 0.3984;      %Length at maximum contraction, m
Vas_Pam_48cm = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);
Theoretical = Vas_Pam_48cm.Torque(:,3);

%% Tests 1 and 4 done with CALT load cell. Tests 2 and 3 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet ExtTest10mm_1 from Results_table10mm_pinned_LoadCell
%Test 2 == sheet ExtTest10mm_2 from Results_table10mm_pinned_FishScale
%Test 3 == sheet ExtTest10mm_3 from Results_table10mm_pinned_FishScale
%Test 4 == sheet ExtTest10mm_4 from Results_table10mm_pinned_LoadCell
%% Torque calculated from measurements

Angle{1} = [-121.5	-108.5	-101.5	-94	-83	-75	-61.5	-65	-32	-18	-32.5	-47.5	-57.5	-69.5	-80	-90.5	-96.5	-104	-113.5	-123]';
Angle{2} = [-110	-82	-71.5	-55	-42	-43	-65.5	-74.5	-72	-82	-88.5	-98	-113.5]';
Angle{3} = [-110	-105	-103.5	-91	-80	-72.5	-62	-48	-40.5	-31	-30	-31.5	-43	-50	-62.5	-69	-78	-88	-93	-100	-106	-107]';
Angle{4} = [-123	-119	-106	-92]';

Torque{1} = [4.62685036	4.503204046	4.29202985	4.002027833	3.417014674	2.826201474	1.538232508	1.214084497	0.439203694	-0.080856656	0.539217698	1.260349872	2.254964352	3.207688111	3.851548806	4.243514283	4.672306017	4.873564883	5.146026298	6.087161098]';
Torque{2} = [5.783034531	4.71210221	3.667943198	1.847358253	1.044159012	1.579625173	2.356051105	3.159250345	3.50730335	4.203409358	4.203409358	4.97983529	4.97983529]';
Torque{3} = [5.783034531	4.71210221	4.471142438	3.935676278	3.667943198	3.159250345	2.623784185	1.847358253	1.419266981	0.53546616	0	0.53546616	1.044159012	2.356051105	2.891517265	3.400210118	3.935676278	4.203409358	4.71210221	4.97983529	6.291727383	6.559460463]';
Torque{4} = [7.013441569	4.72935096	4.346134571	3.454623635]';

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for tests 1 & 2 were done incorrectly and should be disregarded
InflatedLength{1} = [455.5	461	451	451	440	438.5	435	415.5	413	410	409.5	418	426.5	438	442.5	441.5	450	458.5	455.5	463]'/1000;
InflatedLength{2} = [448	445.5	429.5	413.5	401	410	417.5	422.5	430.5	429.5	441.5	439	443.5]'/1000;
InflatedLength{3} = [445	445	433	426	423	421	420	415	405	400	398	400	405	407	416	420	430	433	441	440	450	453]'/1000;
InflatedLength{4} = [459	456	447	438]'/1000;

ICRtoMuscle{1} = [30.5	41.5	39	35.5	37	40.5	37	46.5	38.5	45.5	36.5	38.5	35	38.5	36	42	39	36	38	37]'/1000;
ICRtoMuscle{2} = [46	54.5	53	52	57	49	53	53	54	51	52	53	55.5]'/1000;
ICRtoMuscle{3} = [32	30	30	30	30	30	30	30	35	35	40	35	34	33	30	30	30	30	30	30	30	30]'/1000;
ICRtoMuscle{4} = [30	30	30	30]'/1000;

%load pressure where applicable
test = [1 3];
runsperseries = [20 4];

pres = cell(length(Angle),1);
pres{1} = zeros(runsperseries(1),1);
pres{4} = zeros(runsperseries(2),1);
    for j = 1:runsperseries(1)
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(1),j);
                load(file_name,'Stats')
                pres{1}(j,1) = Stats{'Mean',2};
    end
    for j = 1:runsperseries(2)
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(2),j);
                load(file_name,'Stats')
                pres{4}(j,1) = Stats{'Mean',2};
    end
pres{2} = 612*ones(length(InflatedLength{2}),1);
pres{3} = 612*ones(length(InflatedLength{3}),1);

KMAX = (restingLength-kmax)/restingLength;  %Converts to percentage
strainz = cell(length(Angle),1);
rel = cell(length(Angle),1);
F = cell(length(Angle),1);
TorqueHand = cell(length(Angle),1);

for i = 1:length(Angle)
    strainz{i} = ((restingLength-InflatedLength{i})/restingLength);
    rel{i} = strainz{i}/KMAX;
    F{i} = bpaForce10(restingLength,rel{i},pres{i});
    TorqueHand{i} = ICRtoMuscle{i}.*F{i};  %Torque will be positive because it is causing extension
end 

%% Prepare for plotting

% Since two of the hand measurement recordings was corrupted, use a different 'Angle' for tests 3+4
AngleX = {Angle{3}; Angle{4}};
ICRtoMuscleX = {ICRtoMuscle{3}; ICRtoMuscle{4}};
TorqueHandX = {TorqueHand{3}; TorqueHand{4}};

% Create accessible color scheme
% Matlab hex color values:
c = cell(8,1);
c{1} = '#FFD700'; %gold
c{2} = '#FFB14E'; %orange
c{3} = '#FA8775'; %light orange
c{4} = '#EA5F94'; %pink
c{5} = '#CD34B5'; %magenta
c{6} = '#9D02D7'; %magenta 2
c{7} = '#0000FF'; %indigo
c{8} = '#000000'; %black
sz = 60;        %size of data points

%X Limits
xLim = [-125 40];

%% Plot the expected value and scatter the data that show which test they come from
Test = ["ExtTest10mm-1 pin LoadCell";
        "ExtTest10mm-2 pin FishScale";
        "ExtTest10mm-3 pin FishScale";
        "ExtTest10mm-4 pin LoadCell"];
%% Convert cells to column arrays once bad tests are eliminated
Angle = cell2mat(Angle');
AngleX = cell2mat(AngleX);
Torque = cell2mat(Torque');
InflatedLength = cell2mat(InflatedLength');
ICRtoMuscle = cell2mat(ICRtoMuscle');
ICRtoMuscleX = cell2mat(ICRtoMuscleX);
pres = cell2mat(pres);
strainz = cell2mat(strainz);
rel = cell2mat(rel);
F = cell2mat(F);
TorqueHand = cell2mat(TorqueHand);
TorqueHandX = cell2mat(TorqueHandX);

%% Plot expected versus measured moment arm
Ma = Vas_Pam_48cm.MomentArm;                 %Calculated moment arm
G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque

fig_MA = figure;
ax1 = gca;
hold on
pp = plot(phiD,G,'Color',c{7},'Linewidth',2,'DisplayName','MA expected');
if ~iscell(Angle)
    ss = scatter(Angle, ICRtoMuscle,sz,'filled','MarkerFaceColor',c{5},'DisplayName','MA measured');
else
    for i = 1:length(AngleX)
    ss{i} = scatter(AngleX{i}, ICRtoMuscleX{i},sz,'filled','MarkerFaceColor',c{6-i},'DisplayName',Test{i+2});
    end
end
hold off
title('Expected vs measured moment arm')
xlabel('Knee angle, degrees')
ylabel('Moment Arm, z axis (m)')
set(ax1,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on');
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 8;
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values
strain = Vas_Pam_48cm.Contraction;
relstrain = (strain)./KMAX;

fig_relstrain = figure;
ax2 = gca;
hold on
plot(phiD,relstrain,'Linewidth',2,'DisplayName','Expected Relative Strain')
if ~iscell(Angle)
    sc_rel = scatter(Angle,rel,sz,'filled','MarkerFaceColor',c{5},'DisplayName','Measured Relative Strain');
else
    for i = 1:length(Angle)
        sc_rel{i} = scatter(Angle{i},rel{i},sz,'filled','MarkerFaceColor',c{6-i},'DisplayName',Test(i));
    end
end
hold off
title('Relative strain')
xlabel('Knee angle, \circ')
ylabel('strain/kmax')
set(ax2,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on');
ax2.FontName = 'Arial';
ax2.YAxis.LineWidth = 2; ax2.YAxis.FontSize = 10;
ax2.XAxis.LineWidth = 2; ax2.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 8;
hold off

%% Plot measured versus expected strain (like above, but not normalized
fig_strain = figure;
ax3 = gca;
hold on
pl_strain = plot(phiD,strain,'Color',c{7},'Linewidth',2,'DisplayName','Expected Strain');
if ~iscell(Angle)
    sc_str = scatter(Angle,strainz,sz,'filled','MarkerFaceColor',c{5},'DisplayName','MeasuredStrain');
else
    for i = 1:length(Angle)
        sc_str{i} = scatter(Angle{i},strainz{i},sz,'filled','MarkerFaceColor',c{6-i},'DisplayName',Test(i));
    end
end
hold off
title('Absolute strain')
xlabel('Knee angle, \circ')
ylabel('strain')
set(ax3,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on');
ax3.FontName = 'Arial';
ax3.YAxis.LineWidth = 2; ax3.YAxis.FontSize = 10;
ax3.XAxis.LineWidth = 2; ax3.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 8;
hold off

%% Plot measured versus expected BPA length
MuscleLength = Vas_Pam_48cm.MuscleLength-2*fitting-tendon;

fig_mL = figure;
ax4 = gca;
hold on
plml = plot(phiD,MuscleLength,'Color',c{7},'Linewidth',2,'DisplayName','Expected Muscle Length');
if ~iscell(Angle)
    sc_ML = scatter(Angle,InflatedLength,sz,'filled','MarkerFaceColor',c{5},'DisplayName','Measured Length');
else
    for i = 1:length(Angle)
        sc_mL{i} = scatter(Angle{i},InflatedLength{i},sz,'filled','MarkerFaceColor',c{6-i},'DisplayName',Test(i));
    end
end
hold off
title('Expected vs measured muscle length')
xlabel('Knee angle, \circ')
ylabel('Length, m')
set(ax4,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on');
ax4.FontName = 'Arial';
ax4.YAxis.LineWidth = 2; ax4.YAxis.FontSize = 12;
ax4.XAxis.LineWidth = 2; ax4.XAxis.FontSize = 12;
lgdMa = legend;
lgdMa.FontSize = 8;
hold off

%% Plotting Torque 
fig_T = figure;
gcf_T = gcf;
ax5 = gca;
hold on
PL = plot(phiD, Theoretical,'Color',c{7},'Linewidth',2,'DisplayName','Expected Torque');
if ~iscell(Angle)
    scH = scatter(AngleX,TorqueHandX,sz,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',c{2},'DisplayName','Hybrid calc');
    scM = scatter(Angle,Torque,sz,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',c{5},'DisplayName','Measured Torque');
else
    for i = 1:length(Angle)
    scM{i} = scatter(Angle{i},Torque{i},sz,'filled','d','MarkerFaceColor',c{7-i},'DisplayName',Test(i));
    end
    for j = 1:length(AngleX)
    scH{j} = scatter(AngleX{j},TorqueHandX{j},sz,'MarkerFaceColor',c{5-j},'DisplayName',Test(j+2));
    end
end
hold off
title('l_{rest}=48.0 cm')
xlabel('Knee angle, \circ')
ylabel('Torque, N\cdotm')
% set(gcf2,'Position',[1 384 950 612]);
set(ax5,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on');
ax5.FontName = 'Arial';
ax5.YAxis.LineWidth = 2; ax5.YAxis.FontSize = 12;
ax5.XAxis.LineWidth = 2; ax5.XAxis.FontSize = 12;
lgd = legend;
lgd.FontSize = 8;
