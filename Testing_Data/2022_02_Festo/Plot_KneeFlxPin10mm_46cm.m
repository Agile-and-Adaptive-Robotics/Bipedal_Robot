%% Pinned knee, Extensor
%Run and save data from testing results
clear; 
clc; 
close all;

load KneeFlxPin_10mm_46cm.mat
% Theoretical = TorqueR(:,3);

%Put resting length and kmax again in case something else was run for the
%class we just loaded.
rest = 0.457; %resting length, m
kmax = 0.381; %Length at maximum contraction, m
fitting = 0.0254;
tendon = 0;
restingLength = rest;
Bifemsh_Pam = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, Pres);
Theoretical = Bifemsh_Pam.Torque(:,3);
%% Tests 1, 2, and 3 
%Test 1 == sheet FlxTest10mm_3 from Results_table10mm_pinned_LoadCell
%Test 2 == sheet FlxTest10mm_3 from Results_table10mm_pinned_FishScale
%Test 3 == sheet FlxTest10mm_6 from Results_table10mm_pinned in \2026_06_Festo\
%% Torque calculated from measurements
Angle{1} = [-41	-45	-57	-63	-77	-75	-80	-90	-95	-105]';
Angle{2} = [-30	-32	-36.5	-41	-51	-52.5	-56.5	-63	-66.5	-74.5	-80.5	-80	-82	-86	-90	-85	-81.5	-76	-70	-62	-58	-49	-45	-37	-33.5	-26.5]';
Angle{3} = [-86.5	-78.7	-71	-66	-59.7	-53.6	-45.743	-40.7	-33.9	-25.28	-32	-40.5	-44.3	-50.7	-56.5	-66.5	-71.9	-76.258	-80.9	-21.331];

Torque{1} = [-17.62398227	-14.31258839	-10.54411571	-7.71470304	-5.263349252	-5.04823533	-2.832979547	-1.807118322	-1.02853932	-0.561617562]';
Torque{2} = [-20.00095306	-19.62412624	-17.74104856	-14.64175686	-13.34974574	-11.21468685	-8.847415849	-6.225462211	-5.473176171	-4.159138583	-1.821998215	-1.821998215	-1.293882791	-0.526243409	0	-1.293882791	-2.098936493	-3.894225298	-5.473176171	-7.788450595	-11.02208726	-14.67316852	-18.04885624	-21.08709753	-22.77466623	-19.33565341]';
Torque{3} = [-0.762007711	-2.13625802	-4.098758227	-5.991344332	-8.176277427	-11.3506308	-14.41864492	-15.86146464	-17.40131487	-18.39777443	-16.74423281	-14.34520531	-13.05712202	-10.28841506	-7.588760486	-3.971495011	-2.692785979	-1.787030695	-1.101594264	-18.24645003];
%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength{1} = [428	418	409.5	405.5	395	396	390	384	381	376]'/1000;
InflatedLength{2} = [438	438	434	424	419	415	408	401	399	394	388	388	386	381	380	384	386	390	394	400	407	413	422	426	433	438]'/1000;
InflatedLength{3} = [387	391	395	400	405	412	422	425	430	435	431	423	418	411	405	395	390	387	385	435]'/1000;

ICRtoMuscle{1} = [80	81	78	76	62	62.5	60	48	34	26]'/1000;
ICRtoMuscle{2} = [85	85	85	88	90	88	83	80	75	70	60	64	60	60	63	65	65	68	75	79	84	90	92	90	89	85]'/1000;
ICRtoMuscle{3} = [42	53	66	74	80	80	85	85	85	80	82	83	85	85	80	70	65	55	40	75]'/1000;

%load pressure where applicable
test = 3;
runsperseries = 10;

pres{1} = zeros(runsperseries,1);    
    for j = 1:runsperseries
                file_name = sprintf('FlxTest%0.0f_%0.0f.mat', test,j);
                load(file_name,'Stats')
                pres{1}(j,1) = Stats{'Mean',2};
    end
pres{2} = 606*ones(length(InflatedLength{2}),1);
pres{3} = [613.75	613	613	616	616.7	616	614.5	613	613.764	613	613	615	615	615	613.7	613	615.2	614	613	615.2]';

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

%% Plot the expected value and scatter the data that show which test they come from
% Test = ["FlxTest10mm-3 10mm-pinned-LoadCell";
%         "FlxTest10mm-3 10mm-pinned-FishScale";
%         "FlxTest10mm-6 10mm-pinned"];
Test = ["Test 1";
        "Test 2";
        "Test 3"];
%% Convert cells to column arrays once bad tests are eliminated
% Angle = cell2mat(Angle');
% Torque = cell2mat(Torque');
% InflatedLength = cell2mat(InflatedLength');
% ICRtoMuscle = cell2mat(ICRtoMuscle');
% pres = cell2mat(pres');
% strainz = cell2mat(strainz);
% rel = cell2mat(rel);
% F = cell2mat(F);
% TorqueHand = cell2mat(TorqueHand);
Angle = Angle{3}';
Torque = Torque{3}';
InflatedLength = InflatedLength{3}';
ICRtoMuscle = ICRtoMuscle{3}';
pres = pres{3}';
strainz = strainz{3};
rel = rel{3};
F = F{3};
TorqueHand = TorqueHand{3};

%% Plot expected versus measured moment arm
Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
G = hypot(Ma(:,1),Ma(:,2));         %Moment arm for z-axis torque

fig_MA = figure('Name','Moment Arm');
ax1 = gca;
hold on
pp = plot(phiD,G,'Color',c{7},'Linewidth',2,'DisplayName','MA expected');
if ~iscell(Angle)
    ss = scatter(Angle, ICRtoMuscle,sz,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none','MarkerFaceColor',c{6},'DisplayName','MA measured');
else
    for i = 1:length(Angle)
    ss{i} = scatter(Angle{i}, ICRtoMuscle{i},sz,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none','MarkerFaceColor',c{2*i-1},'DisplayName',Test{i});
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
lgdMa = legend('Location','southeast');
lgdMa.FontSize = 8;
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values
strain = Bifemsh_Pam.Contraction;
relstrain = (strain)./KMAX;

fig_relstrain = figure('Name','Relative Strain');
ax2 = gca;
hold on
plot(phiD,relstrain,'Linewidth',2,'DisplayName','Expected Relative Strain')
if ~iscell(Angle)
    sc_rel = scatter(Angle,rel,sz,'filled','MarkerFaceColor',c{6},'MarkerFaceAlpha',0.75,'MarkerEdgeColor','none','DisplayName','Measured Relative Strain');
else
    for i = 1:length(Angle)
        sc_rel{i} = scatter(Angle{i},rel{i},sz,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none','MarkerFaceColor',c{2*i-1},'DisplayName',Test(i));
    end
end
hold off
title('Relative strain')
xlabel('Knee angle, \circ')
ylabel('strain/kmax')
set(ax2,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax2.FontName = 'Arial';
ax2.YAxis.LineWidth = 2; ax2.YAxis.FontSize = 10;
ax2.XAxis.LineWidth = 2; ax2.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 8;


%% Plot measured versus expected strain (like above, but not normalized
fig_strain = figure('Name','Strain');
ax3 = gca;
hold on
pl_strain = plot(phiD,strain,'Color',c{7},'Linewidth',2,'DisplayName','Expected Strain');
if ~iscell(Angle)
    sc_str = scatter(Angle,strainz,sz,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none','MarkerFaceColor',c{6},'DisplayName','MeasuredStrain');
else
    for i = 1:length(Angle)
        sc_str{i} = scatter(Angle{i},strainz{i},sz,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none','MarkerFaceColor',c{2*i-1},'DisplayName',Test(i));
    end
end
hold off
title('Absolute strain')
xlabel('Knee angle, \circ')
ylabel('strain')
set(ax3,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax3.FontName = 'Arial';
ax3.YAxis.LineWidth = 2; ax3.YAxis.FontSize = 10;
ax3.XAxis.LineWidth = 2; ax3.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 8;
hold off

%% Plot measured versus expected BPA length
MuscleLength = Bifemsh_Pam.MuscleLength-2*fitting-tendon;

fig_mL = figure('Name','Lm');
ax4 = gca;
hold on
plml = plot(phiD,MuscleLength,'Color',c{7},'Linewidth',2,'DisplayName','Expected Muscle Length');
if ~iscell(Angle)
    sc_ML = scatter(Angle,InflatedLength,sz,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none','MarkerFaceColor',c{6},'DisplayName','Measured Length');
else
    for i = 1:length(Angle)
        sc_mL{i} = scatter(Angle{i},InflatedLength{i},sz,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none','MarkerFaceColor',c{2*i-1},'DisplayName','Measured Length');
    end
end
hold off
title('Expected vs measured muscle length')
xlabel('Knee angle, \circ')
ylabel('Length, m')
set(ax4,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax4.FontName = 'Arial';
ax4.YAxis.LineWidth = 2; ax4.YAxis.FontSize = 12;
ax4.XAxis.LineWidth = 2; ax4.XAxis.FontSize = 12;
lgdMa = legend('Location','southeast');
lgdMa.FontSize = 8;
hold off

%% Plotting Torque 
fig_T = figure('Name','Torque');
gcf_T = gcf;
ax5 = gca;
hold on
PL = plot(phiD, Theoretical,'Color',c{7},'Linewidth',2,'DisplayName','Expected Torque');
if ~iscell(Angle)
    scH = scatter(Angle,TorqueHand,sz,'MarkerFaceColor',c{2},'MarkerFaceAlpha',0.75,'MarkerEdgeColor','none','DisplayName','Hybrid calc');
    scM = scatter(Angle,Torque,sz,'filled','d','MarkerFaceColor',c{6},'MarkerFaceAlpha',0.75,'MarkerEdgeColor','none','DisplayName','Measured Torque');
else
    for i = 1:length(Angle)
    strH = sprintf('Hybrid, %s',Test(i));
    scH{i} = scatter(Angle{i},TorqueHand{i},sz,'filled','d','MarkerFaceColor',c{2*i-1},'MarkerFaceAlpha',0.75,'MarkerEdgeColor','none','DisplayName',strH);
    end
    for i = 1:length(Angle)
    strM = sprintf('Measured, %s',Test(i));
    scM{i} = scatter(Angle{i},Torque{i},sz,'MarkerFaceColor',c{2*i-1},'MarkerFaceAlpha',0.75,'MarkerEdgeColor','none','DisplayName',strM);
    end
end
hold off
title('Isometric Torque vs Knee Angle, 10mm Flexor, 45.7cm long')
xlabel('Knee angle, \circ')
ylabel('Torque, $N \cdot m$')
% set(gcf2,'Position',[1 384 950 612]);
set(ax5,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax5.FontName = 'Arial';
ax5.YAxis.LineWidth = 2; ax5.YAxis.FontSize = 12;
ax5.XAxis.LineWidth = 2; ax5.XAxis.FontSize = 12;
lgd = legend('Location','southwest');
lgd.FontSize = 8;