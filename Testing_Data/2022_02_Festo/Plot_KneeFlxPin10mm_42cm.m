%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

load KneeFlxPin_10mm_42cm.mat
% Theoretical = TorqueR(:,3);

%Put resting length and kmax again in case something else was run for the
%class we just loaded.
rest = 0.410; %resting length, m
kmax = 0.335; %Length at maximum contraction, m
fitting = 0.0254;
tendon = 0;
restingLength = rest;
Bifemsh_Pam = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, Pres);
Theoretical = Bifemsh_Pam.Torque(:,3);
%% Test 1
%Test 1 == sheet FlxTest10mm_7 from Results_table10mm_pinned in \2026_06_Festo\
%% Torque calculated from measurements
Angle{1} = [-103.9	-94.2	-86.3	-79.1	-73	-68.36	-73.7	-82.72	-96.3	-105.3	-104.9	-100.6	-97	-92.8	-88.1	-83.4	-78.4	-74.8	-70.5	-62.257	-56	-63]';

Torque{1} = [-2.458276599	-4.81455329	-7.862890964	-11.27489776	-14.46271234	-17.04860431	-13.97377874	-8.933820327	-3.755149753	-1.892440543	-2.064795089	-2.949931919	-3.810328728	-5.113215325	-6.586676176	-8.896928268	-11.6336012	-13.7386268	-16.32140252	-20.60550736	-23.04447722	-18.97492061]';
%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength{1} = [351	355	360	365	371	376	370	361	352	350	348	350	352	355	358	360	365	368	374	381	386	379]'/1000;

ICRtoMuscle{1} = [25	35	65	65	78	80	75	67	40	21	25	34	45	52	57	65	70	69	83	83	98	94]'/1000;

pres{1} = [615.2	617.5	616	616	614.5	614.5	613.7	617.5	616	615	616	617	616.8	615.3	615.2	614.5	614.5	613.7	614	613	616.7	617.5]';    


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
% Test = ["FlxTest10mm-7 10mm-pinned"];
Test = "Test 1";
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
title('Isometric Torque vs Knee Angle, 10mm Flexor, 41.7cm long')
xlabel('Knee angle, \circ')
ylabel('Torque, $N \cdot m$')
% set(gcf2,'Position',[1 384 950 612]);
set(ax5,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax5.FontName = 'Arial';
ax5.YAxis.LineWidth = 2; ax5.YAxis.FontSize = 12;
ax5.XAxis.LineWidth = 2; ax5.XAxis.FontSize = 12;
lgd = legend('Location','southwest');
lgd.FontSize = 8;