%% Pinned knee, Extensor
%Run and save data from testing results
clear;
clc;
close all;

%Load file with all the extensors using the pinned knee
load KneeExtPin_10mm_all.mat
%Specify values for this BPA
restingLength = 0.457;      %resting length, m
kmax = 0.380;               %Length at maximum contraction, m
tendon = 0;             %Tendon, measured
fitting = 0.0254;
Vas_Pam_46cm = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, restingLength, kmax, tendon, fitting, pres);
Theoretical = Vas_Pam_46cm.Torque(:,3);

%% Test 1 done with CALT load cell. Tests 2 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet ExtTest10mm_2 from Results_table_10mm_pinned_LoadCell
%Test 2 == sheet ExtTest10mm_1 from Results_table_10mm_pinned_FishScale
%% Torque calculated from measurements
Angle{1} = [-125	-120.5	-115	-111.5	-108	-101.5	-88	-80	-71	-58	-40	-32	-17.5	-10.5	-4	2	-3.5	-16	-33.5	-48	-58.5	-72	-80	-85	-88.5	-99.5	-111.5	-114.5	-123]';
Angle{2} = [-102	-100	-81.5	-66	-47	-35	-17	-11	-5]';

Torque{1,1} = [6.168692258	6.072984004	5.91673494	5.983690138	6.220291194	5.580594946	4.981638504	4.559983227	4.059476693	3.15664206	2.56531606	2.127108003	1.785781666	0.774715989	0.309696885	-0.057265164	0.536074229	2.318668963	2.990337108	3.61356103	3.999580201	4.737782449	5.291362695	5.556525636	6.028604398	6.318428803	6.41672146	6.668472107	6.512037525]';
Torque{2,1} = [5.969270431	5.781894229	4.978853364	3.667219951	2.623266826	2.087906249	1.579313702	1.311633413	0]';

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded
InflatedLength{1} = [453.5	452	456	450.5	449.5	442.5	439	435	422.5	420	410.5	398.5	400.5	395	384	381	382	394	408.5	414.5	420.5	425	427	432	438	438	441	445	447]'/1000;
InflatedLength{2} = [430	430	420	415	410	401	390	388	385]'/1000;

ICRtoMuscle{1} = [29.5	30	30	29.5	30.5	30	30	32.5	34	35	37	38.5	44.5	50	54	60	55	48	34.5	33	32.5	31.5	30	30	30	30	30	30	30]'/1000;
ICRtoMuscle{2} = [30	30	30	30	34	35	45	50	55]'/1000;

%load pressure where applicable
test = 2;
runsperseries = 29;

pres = cell(length(Angle),1);
pres{1} = zeros(runsperseries,1);
    for j = 1:runsperseries
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test,j);
                load(file_name,'Stats')
                pres{1}(j,1) = Stats{'Mean',2};
    end

pres{2} = 606*ones(length(InflatedLength{2}),1);

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
% Test = ["ExtTest10mm-2 10mm pin LoadCell";
%         "ExtTest10mm-1 10mm pin FishScale"];
%% Convert cells to column arrays once bad tests are eliminated
AngleX = Angle{1};
Angle = cell2mat(Angle');
Torque = cell2mat(Torque);
InflatedLength = InflatedLength{1};
ICRtoMuscle = ICRtoMuscle{1}';
pres = pres{1};
strainz = strainz{1};
rel = rel{1};
F = F{1};
TorqueHand = TorqueHand{1};

%% Plot expected versus measured moment arm
Ma = Vas_Pam_46cm.MomentArm;                 %Calculated moment arm
G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque

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
strain = Vas_Pam_46cm.Contraction;
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
MuscleLength = Vas_Pam_46cm.MuscleLength-2*fitting-tendon;

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
    scM = scatter(Angle,Torque,sz,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',c{7},'DisplayName','Measured Torque');
else
    for i = 1:length(Angle)
    scM{i} = scatter(Angle{i},Torque{i},sz,'filled','d','MarkerFaceColor',c{7-i},'DisplayName',sprintf('Measured, test%d',i));
    scH{i} = scatter(Angle{i},TorqueHand{i},sz,'MarkerFaceColor',c{7-2*i},'DisplayName',sprintf('Back calc, test%d',i));
    end
end
hold off
title('l_{rest}=45.7cm')
xlabel('Knee angle, \circ')
ylabel('Torque, N\cdotm')
% set(gcf2,'Position',[1 384 950 612]);
set(ax5,'FontSize', 12, 'XLim',xLim,'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax5.FontName = 'Arial';
ax5.YAxis.LineWidth = 2; ax5.YAxis.FontSize = 12;
ax5.XAxis.LineWidth = 2; ax5.XAxis.FontSize = 12;
lgd = legend;
lgd.FontSize = 8;