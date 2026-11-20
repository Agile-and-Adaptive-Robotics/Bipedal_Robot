%% Pinned knee, Extensor
%Run and save data from testing results
clear;
clc;
close all;

%Load file with all the extensors using the pinned knee
load KneeExtPin_10mm_all.mat
%Specify values for this BPA
restingLength = 0.465;      %resting length, m
rest = restingLength;
kmax = 0.387;               %Length at maximum contraction, m
tendon = 0;             %Tendon, measured
fitting = 0.0254;
Vas_Pam_47cm = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, restingLength, kmax, tendon, fitting, pres);
Theoretical = Vas_Pam_47cm.Torque(:,3);

%% Test 1 done with CALT load cell. 
%Test 1 == sheet ExtTest10mm_13 from Results_table_10mm_pinned

%% Torque calculated from measurements
Angle = [23.185	17.441	9.184	2.004	-4.817	-10.561	-16.305	-21.331	-26.357	-32.46	-40.717	-50.051	-55.795	-60.821	-67.642	-73.386	-84.156	-90.259	-75.181	-65.847	-51.487	-33.896	-17.023	-4.099	6.312	16.005]';

Torque = [0.177169291	0.235809831	0.318878577	0.723579948	0.892996679	1.941338478	2.420438069	2.732993875	2.89819417	2.962283833	3.090959977	3.555611423	3.850846249	3.629150651	4.182916171	4.464241445	5.351758717	5.524335655	4.643488202	3.933746211	3.236203509	2.428451719	1.795349447	0.662409429	0.277147599	0.238721942]';

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded
InflatedLength = [404	398	393	393	390	394	397	399	400	404	407	410	410	414	415	416	423	424	418	417	410	404	396	387	388	387]'/1000;

ICRtoMuscle = [80	73	65	55	53	50	48	44	40	35	32	30	30	28	28	28	28	27	30	29	33	37	44	50	59	72]'/1000;

%load pressure where applicable
pres = [620.579	623.607	622.093	621.336	621.336	623.607	622.093	622.093	621.336	620.579	620.58	623.607	622.85	623.607	623.607	622.096	620.579	622.85	623.85	622.85	621.336	622.093	623.607	622.85	622.093	621.336]';

KMAX = (restingLength-kmax)/restingLength;  %Converts to percentage
strainz = (restingLength-InflatedLength)./restingLength;
rel = strainz/KMAX;
F = bpaForce10(restingLength,rel,pres);
TorqueHand = ICRtoMuscle.*F;  %Torque will be positive because it is causing extension


%% Plot setup
%Matlab hex color values:
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

%X axis limit
xLim = [-120 35];

%% Plot the expected value and scatter the data that show which test they come from
Test = ["ExtTest10mm-13 10mm pin LoadCell"];

AngleX = Angle;

%% Plot expected versus measured moment arm
Ma = Vas_Pam_47cm.MomentArm;                 %Calculated moment arm
G = hypot(Ma(:,1),Ma(:,2));         %Moment arm for z-axis torque

fig_MA = figure;
ax1 = gca;
hold on
pp = plot(phiD,G,'Color',c{7},'Linewidth',2,'DisplayName','MA expected');
if ~iscell(Angle)
    ss = scatter(AngleX, ICRtoMuscle,sz,'filled','MarkerFaceColor',c{5},'DisplayName','MA measured');
else
    for i = 1:length(Angle)
    ss{i} = scatter(Angle{i}, ICRtoMuscle{i},sz,'filled','MarkerFaceColor',c{7-2*i},'DisplayName',Test{i});
    end
end
hold off
title('Expected vs measured moment arm')
xlabel('Knee angle, degrees')
ylabel('Moment Arm, z axis (m)')
set(ax1,'FontSize', 12, 'XLim',xLim,'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 8;
hold off

%% Plot relative strain, expected and measured
strain = Vas_Pam_47cm.Contraction;
relstrain = (strain)./KMAX;

fig_relstrain = figure;
ax2 = gca;
hold on
plot(phiD,relstrain,'Linewidth',2,'DisplayName','Expected Relative Strain')
if ~iscell(Angle)
    sc_rel = scatter(AngleX,rel,sz,'filled','MarkerFaceColor',c{5},'DisplayName','Measured Relative Strain');
else
    for i = 1:length(Angle)
        sc_rel{i} = scatter(Angle{i},rel{i},sz,'filled','MarkerFaceColor',c{7-2*i},'DisplayName',Test(i));
    end
end
hold off
title('Relative strain')
xlabel('Knee angle, \circ')
ylabel('strain/kmax')
set(ax2,'FontSize', 12, 'XLim',xLim,'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
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
    sc_str = scatter(AngleX,strainz,sz,'filled','MarkerFaceColor',c{5},'DisplayName','MeasuredStrain');
else
    for i = 1:length(Angle)
        sc_str{i} = scatter(Angle{i},strainz{i},sz,'filled','MarkerFaceColor',c{7-2*i},'DisplayName',Test(i));
    end
end
hold off
title('Absolute strain')
xlabel('Knee angle, \circ')
ylabel('strain')
set(ax3,'FontSize', 12, 'XLim',xLim,'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax3.FontName = 'Arial';
ax3.YAxis.LineWidth = 2; ax3.YAxis.FontSize = 10;
ax3.XAxis.LineWidth = 2; ax3.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 8;
hold off

%% Plot measured versus expected BPA length
MuscleLength = Vas_Pam_47cm.MuscleLength-2*fitting-tendon;

fig_mL = figure;
ax4 = gca;
hold on
plml = plot(phiD,MuscleLength,'Color',c{7},'Linewidth',2,'DisplayName','Expected Muscle Length');
if ~iscell(Angle)
    sc_ML = scatter(AngleX,InflatedLength,sz,'filled','MarkerFaceColor',c{5},'DisplayName','Measured Length');
else
    for i = 1:length(Angle)
        sc_mL{i} = scatter(Angle{i},InflatedLength{i},sz,'filled','MarkerFaceColor',c{7-2*i},'DisplayName','Measured Length');
    end
end
hold off
title('Expected vs measured muscle length')
xlabel('Knee angle, \circ')
ylabel('Length, m')
set(ax4,'FontSize', 12, 'XLim',xLim,'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
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
PL = plot(phiD, Theoretical,'Color',c{5},'Linewidth',2,'DisplayName','Expected Torque');
if ~iscell(Angle)
    scH = scatter(AngleX,TorqueHand,sz,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',c{2},'DisplayName','Hybrid calc');
    scM = scatter(AngleX,Torque,sz,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',c{7},'DisplayName','Measured Torque');
else
    for i = 1:length(Angle)
    scM{i} = scatter(Angle{i},Torque{i},sz,'filled','d','MarkerFaceColor',c{7-i},'DisplayName',sprintf('Measured, test%d',i));
    scH{i} = scatter(Angle{i},TorqueHand{i},sz,'MarkerFaceColor',c{7-2*i},'DisplayName',sprintf('Back calc, test%d',i));
    end
end
hold off
title('l_{rest}=46.5cm')
xlabel('Knee angle, \circ')
ylabel('Torque, N\cdotm')
% set(gcf2,'Position',[1 384 950 612]);
set(ax5,'FontSize', 12, 'XLim',xLim,'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax5.FontName = 'Arial';
ax5.YAxis.LineWidth = 2; ax5.YAxis.FontSize = 12;
ax5.XAxis.LineWidth = 2; ax5.XAxis.FontSize = 12;
lgd = legend;
lgd.FontSize = 8;