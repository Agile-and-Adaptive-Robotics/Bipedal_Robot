%% Pinned knee, Extensor
%Run and save data from testing results
clear;
clc;
close all;

load KneeExtPin_10mm_all.mat
rest = cell(2,1);
kmax = cell(2,1);
tendon = cell(2,1);


rest{1} = 0.400;        %resting length clamp to clamp, minus the barb
kmax{1} = 0.341;           %length at maximum contraction
tendon{1} = 0;            %no tendon condition
pres = 624.9;        %Average Pressure, kPa
Vas_Pam_40cm = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest{1}, kmax{1}, tendon{1}, fitting, pres);
Theoretical = Vas_Pam_40cm.Torque(:,3);

rest{2} = 0.400;
kmax{2} = 0.341;
tendon{2} = 0.040;       %tendon
pres = 618.37;         %average pressure, kPa
Vas_Pam_40cm_tendon = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest{2}, kmax{2}, tendon{2}, fitting, pres);
Theoretical_ten = Vas_Pam_40cm_tendon.Torque(:,3);

%% Tests done with CALT load cell.
%Test 1 == sheet ExtTest10mm_9 from Results_table_10mm_pinned_LoadCell   (no tendon)
%Test 2 == sheet ExtTest10mm_10 from Results_table_10mm_pinned_LoadCell   (tendon) 

%% Torque calculated from measurements
Angle{1} = [28.29	18.15	1.29	-14.86	-25.64	-36.76	-42.87	-49.33	-56.87	-68.36	-77.69	-84.16	-90.1	-99.234	-82	-65.48	-47.538	-26.72	-6.253	31.442]';
Angle{2} = [-101	-84.5	-61	-25	-2.5	4	-14.5	-30	-40	-53.5	-63	-85	-96	-106.5]';

Torque{1} = [0.389344556	1.811294821	4.475091751	5.444622252	5.946945384	5.659055943	6.092765212	6.097614209	6.327941585	7.009503763	7.571177124	7.937808716	8.330691929	9.369972199	7.089234174	5.479367043	4.957246849	4.893517243	4.550471361	0.100372374]';
Torque{2} = [4.107694117	3.248828247	2.648288775	1.998570328	0.745228637	0.249401141	2.374751759	3.256929424	2.737149186	3.29731822	3.603151379	3.80467223	4.347910148	5.068763941]';

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
InflatedLength{1} = [367	367	382	398	397	403	407	411	416	422	428	428	434	434	425	420	398	396	386	366]'/1000;
InflatedLength{2} = [405	395	387	375	363	356	370	375	381	388	390	396	401	403]'/1000;

ICRtoMuscle{1} = [78	61	55	48	38	32	30	30	30	29	28	26	23	23	24	29	31	35	50	88]'/1000;
ICRtoMuscle{2} = [30	30	30	43	55	60	45	40	42	30	30	30	30	25]'/1000;

pres = cell(length(Angle),1);
pres{1} = [622.093	620.579	621.336	621.34	621.336	620.58	622.85	622.093	621.336	623.05	623.336	620.579	622.85	621.336	621.336	622.09	622.6	622.85	621.336	620.579]';
pres{2} = [624	622	622	624	623	623	622	622	622	625	624	623	623	622]';

KMAX = cell(length(Angle),1);
strainz = cell(length(Angle),1);
rel = cell(length(Angle),1);
F = cell(length(Angle),1);
TorqueHand = cell(length(Angle),1);

% korr = [0 0 0];           % correction factor
% r = 0.085;                      %radius of curvature
% wR = 15;           % Angle (deg) that wrapping starts to occur
for i = 1:length(Angle)
    KMAX{i} = (rest{i}-kmax{i})./rest{i};
    strainz{i} = ((rest{i}-InflatedLength{i})./rest{i});
    % for j = 1:length(Angle{i})
    %     if Angle{i}(j)<=wR
    %         strainz{i}(j) = ((rest{i}-(InflatedLength{i}(j)-korr(i)*ICRtoMuscle{i}(j)*deg2rad(wR - Angle{i}(j))))./rest{i});
    %     else
    %     end
    % end
    rel{i} = strainz{i}./KMAX{i};
%     for j = 1:length(Angle{i})
    F{i} = bpaForce10(rest{i},rel{i},pres{i});
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
% Test = ["ExtTest10mm-4 10mm pin LoadCell (no tendon)";
%         "ExtTest10mm-5 10mm pin LoadCell (tendon)";
%         "ExtTest10mm-6 10mm pin LoadCell (tendon)"];
Test = ["ExtTest 12";
        "ExtTest 11"];
%% Convert cells to column arrays 

Angle0 = Angle{1};
Angle1 = cell2mat(Angle(2));
Angle = cell2mat(Angle');       %This makes the if statements in the later code work
Torque0 = Torque{1};
Torque1 = cell2mat(Torque(2));
InflatedLength0 = InflatedLength{1};
InflatedLength1 = cell2mat(InflatedLength(2));
ICRtoMuscle0 = ICRtoMuscle{1};
ICRtoMuscle1 = cell2mat(ICRtoMuscle(2));
% pres = cell2mat(pres);
strainz = cell2mat(strainz);
rel0 = rel{1};
rel1 = cell2mat(rel(2));
F = cell2mat(F);
TorqueHand0 = TorqueHand{1};
TorqueHand1 = cell2mat(TorqueHand(2));
KMAX0 = KMAX{1};
KMAX1 = KMAX{2};


%% Plot expected versus measured moment arm
Ma1 = Vas_Pam_40cm.MomentArm;                  %Calculated moment arm
G1 = (Ma1(:,1).^2+Ma1(:,2).^2).^(1/2);         %Moment arm for z-axis torque
Ma2 = Vas_Pam_40cm_tendon.MomentArm;                  %Calculated moment arm
G2 = (Ma2(:,1).^2+Ma2(:,2).^2).^(1/2);         %Moment arm for z-axis torque

fig_MA = figure('Name','Moment Arm');
figMA = tiledlayout(2,1);
ax1_1 = nexttile(1);
hold on
pp1 = plot(phiD,G1,'Color',c{5},'Linewidth',2,'DisplayName','MA expected');
if ~iscell(Angle)
    ss1 = scatter(Angle0, ICRtoMuscle0,sz,'filled','MarkerFaceColor',c{7},'DisplayName','MA measured');
else
    for i = 1
    ss{i} = scatter(Angle{i}, ICRtoMuscle{i},sz,'filled','MarkerFaceColor',c{9-2*i},'DisplayName',Test{i});
    end
end
hold off
title(sprintf(' \bf Torque, l_{rest} = %0.1f cm, %0.1f mm tendon',rest{1}*100, tendon{1}*10^3),'Interpreter','tex')
xlabel(' \bf Knee angle, \circ')
ylabel(' \bf Moment Arm, m')
set(ax1_1,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax1_1.FontName = 'Arial';
ax1_1.YAxis.LineWidth = 2; ax1_1.YAxis.FontSize = 10;
ax1_1.XAxis.LineWidth = 2; ax1_1.XAxis.FontSize = 10;
lgdMa1 = legend('Location','northwest');
lgdMa1.FontSize = 8;

ax1_2 = nexttile(2);
hold on
pp2 = plot(phiD,G2,'Color',c{5},'Linewidth',2,'DisplayName','MA expected');
if ~iscell(Angle)
    ss2 = scatter(Angle1, ICRtoMuscle1,sz,'filled','MarkerFaceColor',c{7},'DisplayName','MA measured');
else
    for i = 2:length(Angle)
    ss{i} = scatter(Angle{i}, ICRtoMuscle{i},sz,'filled','MarkerFaceColor',c{11-2*i},'DisplayName',Test{i});
    end
end
hold off
title(sprintf(' \bf l_{rest} = %0.1f cm, %0.1f mm tendon',rest{2}*100, tendon{2}*10^3),'Interpreter','tex')
xlabel(' \bf Knee angle, \circ')
ylabel(' \bf Moment Arm, m')
set(ax1_2,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax1_2.FontName = 'Arial';
ax1_2.YAxis.LineWidth = 2; ax1_2.YAxis.FontSize = 10;
ax1_2.XAxis.LineWidth = 2; ax1_2.XAxis.FontSize = 10;
lgdMa2 = legend('Location','northwest');
lgdMa2.FontSize = 8;
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values
%Calculated strain
strain = Vas_Pam_40cm.Contraction;          %strain w/o tendon
strain_ten = Vas_Pam_40cm_tendon.Contraction;  %strain w/ tendon
%Calculated relative strain
if ~iscell(Angle)
    relstrain = (strain)./KMAX0;                             %w/o tendon
    relstrain_ten = (strain_ten)./KMAX1;                     %w/ tendon
else
    relstrain = (strain)./KMAX{1};                             %w/o tendon
    relstrain_ten = (strain_ten)./KMAX{2};                     %w/ tendon
end

fig_relstrain = figure('Name','Relative Strain');
figrelstrain = tiledlayout(2,1);
ax2_1 = nexttile(1);
hold on
plot(phiD,relstrain,'Linewidth',2,'DisplayName','Expected Relative Strain')
if ~iscell(Angle)
    sc_rel = scatter(Angle0,rel0,sz,'filled','MarkerFaceColor',c{6},'MarkerFaceAlpha',0.75,'DisplayName','Measured Relative Strain');
else
    for i = 1
        sc_rel{i} = scatter(Angle{i},rel{i},sz,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',c{8-2*i},'DisplayName',Test(i));
    end
end
hold off
title(sprintf(' \bf l_{rest} = %0.1f cm, %0.1f mm tendon',rest{1}*100, tendon{1}*10^3),'Interpreter','tex')
xlabel('Knee angle, \circ')
ylabel('strain/kmax')
set(ax2_1,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05],'YLim',[0 1.5],'XLim',[-125 35]);
ax2_1.FontName = 'Arial';
ax2_1.YAxis.LineWidth = 2; ax2_1.YAxis.FontSize = 10;
ax2_1.XAxis.LineWidth = 2; ax2_1.XAxis.FontSize = 10;
lgd2_1 = legend('Location','northwest');
lgd2_1.FontSize = 8;

ax2_2 = nexttile(2);
hold on
plot(phiD,relstrain_ten,'Linewidth',2,'DisplayName','Expected Relative Strain')
if ~iscell(Angle)
    sc_rel = scatter(Angle1,rel1,sz,'filled','MarkerFaceColor',c{6},'MarkerFaceAlpha',0.75,'DisplayName','Measured Relative Strain');
else
    for i = 2:length(Angle)
        sc_rel{i} = scatter(Angle{i},rel{i},sz,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',c{8-2*i},'DisplayName',Test(i));
    end
end
hold off
title(sprintf(' \bf l_{rest} = %0.1f cm, %0.1f mm tendon',rest{2}*100, tendon{2}*10^3),'Interpreter','tex')
xlabel('Knee angle, \circ')
ylabel('strain/kmax')
set(ax2_2,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05],'YLim',[0 1.5],'XLim',[-125 35]);
ax2_2.FontName = 'Arial';
ax2_2.YAxis.LineWidth = 2; ax2_2.YAxis.FontSize = 10;
ax2_2.XAxis.LineWidth = 2; ax2_2.XAxis.FontSize = 10;
lgd2_2 = legend('Location','northwest');
lgd2_2.FontSize = 8;

%% Plot measured versus expected BPA length
MuscleLength = Vas_Pam_40cm.MuscleLength-2*fitting-tendon{1};
MuscleLength_ten = Vas_Pam_40cm_tendon.MuscleLength-2*fitting-tendon{2};

figure('Name','BPA length')
xx = tiledlayout(2,1);
Lm1 = nexttile(1);
hold on
pLm1_1 = plot(phiD,MuscleLength,'Color',c{5},'DisplayName',' \bf Predicted');
if ~iscell(Angle)
    sLm1_1 = scatter(Angle0,InflatedLength0,sz,'MarkerFaceColor',c{7},'DisplayName',' \bf Measured');
else
    sLm1_1 = scatter(Angle{1},InflatedLength{1},'DisplayName',' \bf Measured');
end
hold off
title(sprintf(' \bf L_{m}, l_{rest} = %0.1f cm, %0.1f mm tendon',rest{1}*100, tendon{1}*10^3),'Interpreter','tex')
xlabel(' \bf Knee angle, \circ','Interpreter','tex')
ylabel(' \bf l_{M}, m','Interpreter','tex')
set(Lm1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial','TickLength',[0.025, 0.05])
set(Lm1,'XLim',xLim,'XMinorTick','on','YMinorTick','on')
lgdLm1 = legend('Interpreter','tex','Location','southwest');
lgdLm1.FontSize = 8;

Lm2 = nexttile(2);
hold on
pLm2_1 = plot(phiD,MuscleLength_ten,'Color',c{5},'DisplayName',' \bf Predicted');
if ~iscell(Angle)
sLm2_1 = scatter(Angle1,InflatedLength1,sz,'MarkerFaceColor',c{7},'DisplayName',' \bf Measured');
else
sLm2_1 = scatter(Angle{2},InflatedLength{2},'DisplayName',' \bf Measured');
end
hold off
title(sprintf(' \bf L_{m}, l_{rest} = %0.1f cm, %0.1f mm tendon',rest{2}*100, tendon{2}*10^3),'Interpreter','tex')
xlabel(' \bf Knee angle, \circ','Interpreter','tex')
ylabel(' \bf l_{m}, m','Interpreter','tex')
set(Lm2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial','TickLength',[0.025, 0.05])
set(Lm2,'XLim',xLim,'XMinorTick','on','YMinorTick','on')
lgdLm2 = legend('Interpreter','tex');
lgdLm2.FontSize = 8;
set(Lm2,'XLim',xLim,'YLim',[0.34 0.46])
hold off

%% Plotting Z axis Torque values
figure
gca1 = gca;
hold on
PL1 = plot(phiD, Theoretical,'Color',c{5},'Linewidth',2,'DisplayName','Theoretical');
if ~iscell(Angle)
    scM = scatter(Angle0,Torque0,sz,'filled','MarkerFaceColor',c{7},'DisplayName','Measured');
    scH = scatter(Angle0,TorqueHand0,sz,'filled','MarkerFaceColor',c{2},'DisplayName','Back calculated');
else
    scM = scatter(Angle{1},Torque{1},sz,'d','filled','MarkerFaceColor',c{7},'DisplayName','Measured');
    scH = scatter(Angle{1},TorqueHand{1},sz,'filled','MarkerFaceColor',c{2},'DisplayName','Back calculated');
end
hold off
title(sprintf(' \bf Torque, l_{rest} = %0.1f cm, %0.1f mm tendon',rest{1}*100, tendon{1}*10^3),'Interpreter','tex')
xlabel('Knee angle, \circ','FontWeight','bold','Interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','Interpreter','tex')
set(gca1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',xLim,'TickLength',[0.025, 0.05])
set(gca1,'XMinorTick','on','YMinorTick','on');
lgd1 = legend;
lgd1.FontSize = 12;
hold off

figure
gca2 = gca;
hold on
PL2 = plot(phiD, Theoretical_ten,'Color',c{5},'Linewidth',2,'DisplayName',' \bf Theoretical, w/ tendon');
if ~iscell(Angle)
 scM1 = scatter(Angle1,Torque1,sz,'filled','MarkerFaceColor',c{7},'DisplayName','Measured');
 scH1 = scatter(Angle1,TorqueHand1,sz,'filled','MarkerFaceColor',c{2},'DisplayName','Hybrid calc'); 
else
 scM1 = scatter(Angle{2},Torque{2},sz,'d','filled','MarkerFaceColor',c{6},'DisplayName',Test(2));
 scH1 = scatter(Angle{2},TorqueHand{2},sz,'filled','MarkerFaceColor',c{4},'DisplayName',Test(2));
end
hold off
title(sprintf(' \bf Torque, l_{rest} = %0.1f cm, %0.1f mm tendon',rest{2}*100, tendon{2}*10^3),'Interpreter','tex')
xlabel('Knee angle, \circ','FontWeight','bold','Interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','Interpreter','tex')
set(gca2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',xLim,'TickLength',[0.025, 0.05])
set(gca2,'XMinorTick','on','YMinorTick','on');
lgd2 = legend('Interpreter','tex');
lgd2.FontSize = 12;
hold off

