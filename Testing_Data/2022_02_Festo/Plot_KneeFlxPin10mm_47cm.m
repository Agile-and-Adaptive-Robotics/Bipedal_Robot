%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

load KneeFlxPin_10mm_47cm.mat
% Theoretical = TorqueR(:,3);

%Put resting length and kmax again in case something else was run for the
%class we just loaded.
rest = 0.479; %resting length, m
kmax = 0.403;     %Length at maximum contraction, m  %Originally recorded as 413 mm
tendon = 0;
Pres = 623.17;      %average pressure (kPa) for 46cm
restingLength = rest;
Bifemsh_Pam = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, Pres);
Theoretical = Bifemsh_Pam.Torque(:,3);
clear Angle Torque InflatedLength ICRtoMuscle pres
%% Tests 1
%Test 1 == sheet FlxTest10mm_4 from Results_table_10mm_pinned

%% Torque calculated from measurements
Angle{1} = [-68.719	-61.898	-54.359	-46.461	-38.204	-32.101	-27.8	-24.203	-25.6	-16	-20.6	-70.81	-63	-53.6	-39.281	-29.6	-16.3	-10.2]';

Torque{1} = [-1.988844602	-4.252795635	-7.353988604	-10.52740929	-13.96640219	-15.47303695	-16.44526799	-17.07325998	-14.69719212	-15.35889509	-15.04293788	-0.564334071	-2.444541016	-5.240218474	-10.95843156	-13.90243892	-16.84615502	-9.932079547]';

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.

InflatedLength{1} = [408	410	421	427	437	441	444	447	445	452	449	403	407	415	430	438	450	453]'/1000;

ICRtoMuscle{1} = [62	71	72	82	82	80	79	77	75	74	75	55	65	75	76	75	75	67]'/1000;


%load pressure where applicable
pres{1} = [621.336	623.607	622.85	622.85	621.336	621.336	623.607	622.85	622.8	623.6	622	626	624	622	625.1	625.1	623.6	623]';

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
c{5} = '#9D02D7'; %magenta 2
c{7} = '#0000FF'; %indigo
c{8} = '#000000'; %black
sz = 60;        %size of data points

%% Plot the Predicted value and scatter the data that show which test they come from
Test = ["FlxTest10mm_4--10mm_pinned"];

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

%% Plot Predicted versus measured moment arm
Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque

fig_MA = figure('Name','Moment Arm');
ax1 = gca;
hold on
pp = plot(phiD,G,'Color',c{7},'Linewidth',2,'DisplayName','MA Predicted');
if ~iscell(Angle)
    ss = scatter(Angle, ICRtoMuscle,sz,'filled','MarkerEdgeColor','none','MarkerFaceColor',c{5},'DisplayName','MA measured');
else
    for i = 1:length(Angle)
    ss{i} = scatter(Angle{i}, ICRtoMuscle{i},sz,'filled','MarkerEdgeColor','none','MarkerFaceColor',c{6-i},'DisplayName',Test{i});
    end
end
hold off
title('Predicted vs measured moment arm')
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
strain = Bifemsh_Pam.Contraction;
relstrain = (strain)./KMAX;

fig_relstrain = figure('Name','Relative Strain');
ax2 = gca;
hold on
plot(phiD,relstrain,'Linewidth',2,'DisplayName','Predicted Relative Strain')
if ~iscell(Angle)
    sc_rel = scatter(Angle,rel,sz,'filled','MarkerFaceColor',c{5},'MarkerEdgeColor','none','DisplayName','Measured Relative Strain');
else
    for i = 1:length(Angle)
        sc_rel{i} = scatter(Angle{i},rel{i},sz,'filled','MarkerEdgeColor','none','MarkerFaceColor',c{6-i},'DisplayName',Test(i));
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


%% Plot measured versus Predicted strain (like above, but not normalized
fig_strain = figure('Name','Contraction');
ax3 = gca;
hold on
pl_strain = plot(phiD,strain,'Color',c{7},'Linewidth',2,'DisplayName','Predicted Strain');
if ~iscell(Angle)
    sc_str = scatter(Angle,strainz,sz,'filled','MarkerEdgeColor','none','MarkerFaceColor',c{5},'DisplayName','MeasuredStrain');
else
    for i = 1:length(Angle)
        sc_str{i} = scatter(Angle{i},strainz{i},sz,'filled','MarkerEdgeColor','none','MarkerFaceColor',c{7-2i},'DisplayName',Test(i));
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

%% Plot measured versus Predicted BPA length
MuscleLength = Bifemsh_Pam.MuscleLength-2*fitting-tendon;

fig_mL = figure('Name','Muscle Length');
ax4 = gca;
hold on
plml = plot(phiD,MuscleLength,'Color',c{7},'Linewidth',2,'DisplayName','Predicted Muscle Length');
if ~iscell(Angle)
    sc_ML = scatter(Angle,InflatedLength,sz,'filled','MarkerEdgeColor','none','MarkerFaceColor',c{5},'DisplayName','Measured Length');
else
    for i = 1:length(Angle)
        sc_mL{i} = scatter(Angle{i},InflatedLength{i},sz,'filled','MarkerEdgeColor','none','MarkerFaceColor',c{6-i},'DisplayName','Measured Length');
    end
end
hold off
title('Predicted vs measured muscle length')
xlabel('Knee angle, \circ')
ylabel('Length, m')
set(ax4,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax4.FontName = 'Arial';
ax4.YAxis.LineWidth = 2; ax4.YAxis.FontSize = 12;
ax4.XAxis.LineWidth = 2; ax4.XAxis.FontSize = 12;
lgdMa = legend;
lgdMa.FontSize = 8;
hold off

%% Plotting Torque 
fig_T = figure('Name','Torque');
gcf_T = gcf;
ax5 = gca;
hold on
PL = plot(phiD, Theoretical,'Color',c{7},'Linewidth',2,'DisplayName','Predicted Torque');
if ~iscell(Angle)
    scH = scatter(Angle,TorqueHand,sz,'MarkerFaceColor',c{2},'MarkerEdgeColor','none','DisplayName','Hybrid calc');
    scM = scatter(Angle,Torque,sz,'filled','MarkerFaceColor',c{5},'MarkerEdgeColor','none','DisplayName','Measured Torque');
else
    for i = 1:length(Angle)
    scM{i} = scatter(Angle{i},Torque{i},sz,'filled','d','MarkerFaceColor',c{6-i},'MarkerEdgeColor','none','DisplayName',sprintf('Measured, test%d',i));
    scH{i} = scatter(Angle{i},TorqueHand{i},sz,'MarkerFaceColor',c{1+i},'MarkerEdgeColor','none','DisplayName',sprintf('Back calc, test%d',i));
    end
end
hold off
title('Isometric Torque vs Knee Angle, 10 mm Flexor, 47.9 cm long')
xlabel('Knee angle, \circ')
ylabel('Torque, $N \cdot m$')
% set(gcf2,'Position',[1 384 950 612]);
set(ax5,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax5.FontName = 'Arial';
ax5.YAxis.LineWidth = 2; ax5.YAxis.FontSize = 12;
ax5.XAxis.LineWidth = 2; ax5.XAxis.FontSize = 12;
lgd = legend;
lgd.FontSize = 8;