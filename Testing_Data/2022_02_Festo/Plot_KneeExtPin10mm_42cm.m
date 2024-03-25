%% Pinned knee, Extensor
%Run and save data from testing results
clear;
clc;
% close all;

load KneeExtPin_10mm_all.mat
rest = cell(3,1);
kmax = cell(3,1);
tendon = cell(3,1);
rest{1} = 0.415;        %resting length clamp to clamp, minus the barb
kmax{1} = 0.340;           %length at maximum contraction
pres = 605.2351;        %Pressure, kPa
tendon{1} = 0;            %no tendon condition
Vas_Pam_42cm = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest{1}, kmax{1}, tendon{1}, fitting, pres);
Theoretical = Vas_Pam_42cm.Torque(:,3);

rest{2} = 0.415;
kmax{2} = 0.340;
tendon{2} = 0.032;       %22 mm tendon
rest{3} = rest{2};         %repeat  
kmax{3} = kmax{2};
tendon{3} = tendon{2}; 
Vas_Pam_42cm_tendon = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest{2}, kmax{2}, tendon{2}, fitting, pres);
Theoretical_ten = Vas_Pam_42cm_tendon.Torque(:,3);

%% Test 1 done with CALT load cell. Tests 2 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet ExtTest10mm_4 from Results_table10mm_pinned_LoadCell   (no tendon)
%Test 2 == sheet ExtTest10mm_5 from Results_table10mm_pinned_LoadCell   (tendon, test w/ wrapping point #4 +20mm Z direction (i.e. BPA slipped from bolt) 
%Test 3 == sheet ExtTest10mm_6 from Results_table10mm_pinned_LoadCell   (tendon)

%% Torque calculated from measurements
Angle{1} = [-76.5	-69	-51.5	-40.5	-39.5	-55.5	-35	-25.5	-20	-12.5	-5	0.5	4	12.5	20	28	35	23	20	11	1	-16	-24	-44	-53.4	-66]';
Angle{2} = [-120	-109	-101	-97	-89	-87	-81.5	-75	-74	-59.5	-48.5	-3	-30	-18.5]';
Angle{3} = [1.5	-5.5	-14.5	-20	-28.5	-35.5	-40	-50.5	-60	-66.5	-72.5	-78.5	-88	-94.5	-105	-112.5	-119.5]';

Torque{1} = [8.50477362	7.539449285	6.619468483	5.835581619	6.120370289	7.005900186	6.178291184	6.391414397	5.96804559	6.099105123	5.93740855	4.880958611	4.375659895	3.388083398	2.795862193	1.758237576	0.306753172	1.725902869	3.204095044	5.145222736	6.304242352	6.817870142	7.511919159	7.911648561	8.110527045	8.862066797]';
Torque{2} = [4.032777251	3.725074705	3.501289111	3.342545878	3.070623054	2.844765784	2.844621771	2.773344025	2.674409973	2.382410213	1.680520864	0.394827933	1.758713268	1.550440294]';
Torque{3} = [0.55424608	0.941512328	1.76545047	1.79737823	2.650554833	2.7667878	2.211239585	3.275126705	3.423290801	3.611274245	3.687991538	3.833457994	3.587645288	3.554650905	3.949198619	4.568831269	4.408023267]';

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
InflatedLength{1} = [413.5	405.5	402	400	403	412	401	398	392.5	385	385.5	371.5	365	361	356	355	346	348	354	365	376	398	395	400	407	413]'/1000;
InflatedLength{2} = [408	402	380	394	390.5	387	383.5	383	381	371	363	340	355	353]'/1000;
InflatedLength{3} = [340	340	350	350.5	354	362	366	363	368	373	376	383	384	394	391	394.5	395]'/1000;

ICRtoMuscle{1} = [28	30.5	31.5	31	32	31	35	32	41	44.5	45.5	55.5	66.5	62	65	68.5	71	67	60	55	55	42	38	31	29	28]'/1000;
ICRtoMuscle{2} = [42	42	42	42	41.5	36	37	36	36	30.5	22	42	28	24]'/1000;
ICRtoMuscle{3} = [55	43.5	42.5	37	35	33	35	34	35	34	36	36	36	36	36	36	36]'/1000;

%load pressure where applicable
test = [4 5 6];
runsperseries = [26 14 17];

pres = cell(length(Angle),1);
pres{1} = zeros(runsperseries(1),1);
pres{2} = zeros(runsperseries(2),1);
pres{3} = zeros(runsperseries(3),1);
    for i = 1:3
        for j = 1:runsperseries(i)
            if i == 1
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres{1}(j,1) = Stats{'Mean',2};
            elseif i ==2
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres{2}(j,1) = Stats{'Mean',2};
            elseif i ==3
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres{3}(j,1) = Stats{'Mean',2};
            else
            end
        end       
    end

KMAX = cell(length(Angle),1);
strainz = cell(length(Angle),1);
rel = cell(length(Angle),1);
F = cell(length(Angle),1);
TorqueHand = cell(length(Angle),1);

for i = 1:length(Angle)
    KMAX{i} = (rest{i}-kmax{i})./rest{i};
    strainz{i} = ((rest{i}-InflatedLength{i})./rest{i});
    rel{i} = strainz{i}./KMAX{i};
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
Test = ["ExtTest 4";
        "ExtTest 5";
        "ExtTest 6"];
%% Convert cells to column arrays once bad tests are eliminated
Angle1 = cell2mat([Angle(2); Angle(3)]);
Angle0 = Angle{1};
Angle = cell2mat(Angle');       %This makes the if statements in the later code work
Torque0 = Torque{1};
Torque1 = cell2mat([Torque(2); Torque(3)]);
InflatedLength0 = InflatedLength{1};
InflatedLength1 = cell2mat([InflatedLength(2); InflatedLength(3)]);
ICRtoMuscle0 = ICRtoMuscle{1};
ICRtoMuscle1 = cell2mat([ICRtoMuscle(2); ICRtoMuscle(3)]);
pres = cell2mat(pres);
strainz = cell2mat(strainz);
rel0 = rel{1};
rel1 = cell2mat([rel(2); rel(3)]);
F = cell2mat(F);
TorqueHand0 = TorqueHand{1};
TorqueHand1 = cell2mat([TorqueHand(2); TorqueHand(3)]);
KMAX0 = KMAX{1};
KMAX1 = KMAX{2};

%% Plot expected versus measured moment arm
Ma1 = Vas_Pam_42cm.MomentArm;                  %Calculated moment arm
G1 = (Ma1(:,1).^2+Ma1(:,2).^2).^(1/2);         %Moment arm for z-axis torque
Ma2 = Vas_Pam_42cm_tendon.MomentArm;                  %Calculated moment arm
G2 = (Ma2(:,1).^2+Ma2(:,2).^2).^(1/2);         %Moment arm for z-axis torque

fig_MA = figure;
ax1_1 = subplot(2,1,1);
hold on
pp = plot(phiD,G1,'Color',c{5},'Linewidth',2,'DisplayName','MA expected');
if ~iscell(Angle)
    ss = scatter(Angle0, ICRtoMuscle0,sz,'filled','MarkerFaceColor',c{7},'DisplayName','MA measured');
else
    for i = 1
    ss{i} = scatter(Angle{i}, ICRtoMuscle{i},sz,'filled','MarkerFaceColor',c{7-2*i},'DisplayName',Test{i});
    end
end
hold off
title('\bf l_{rest} = 41.5cm, no tendon')
xlabel('\bf Knee angle, \circ')
ylabel('\bf Moment Arm, m')
set(ax1_1,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax1_1.FontName = 'Arial';
ax1_1.YAxis.LineWidth = 2; ax1_1.YAxis.FontSize = 10;
ax1_1.XAxis.LineWidth = 2; ax1_1.XAxis.FontSize = 10;
lgdMa1 = legend('Location','northwest');
lgdMa1.FontSize = 8;

ax1_2 = subplot(2,1,2);
hold on
pp = plot(phiD,G1,'Color',c{5},'Linewidth',2,'DisplayName','MA expected');
if ~iscell(Angle)
    ss = scatter(Angle1, ICRtoMuscle1,sz,'filled','MarkerFaceColor',c{7},'DisplayName','MA measured');
else
    for i = 2:length(Angle)
    ss{i} = scatter(Angle{i}, ICRtoMuscle{i},sz,'filled','MarkerFaceColor',c{7-2*i},'DisplayName',Test{i});
    end
end
hold off
title('\bf l_{rest} , no tendon')
xlabel('\bf Knee angle, \circ')
ylabel('\bf Moment Arm, m')
set(ax1_2,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
ax1_2.FontName = 'Arial';
ax1_2.YAxis.LineWidth = 2; ax1_2.YAxis.FontSize = 10;
ax1_2.XAxis.LineWidth = 2; ax1_2.XAxis.FontSize = 10;
lgdMa2 = legend('Location','northwest');
lgdMa2.FontSize = 8;
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values
%Calculated strain
strain = Vas_Pam_42cm.Contraction;          %strain w/o tendon
strain_ten = Vas_Pam_42cm_tendon.Contraction;  %strain w/ tendon
%Calculated relative strain
if ~iscell(Angle)
    relstrain = (strain)./KMAX0;                             %w/o tendon
    relstrain_ten = (strain_ten)./KMAX1;                     %w/ tendon
else
    relstrain = (strain)./KMAX{1};                             %w/o tendon
    relstrain_ten = (strain_ten)./KMAX{2};                     %w/ tendon
end

fig_relstrain = figure;
ax2_1 = subplot(2,1,1);
hold on
plot(phiD,relstrain,'Linewidth',2,'DisplayName','Expected Relative Strain')
if ~iscell(Angle)
    sc_rel = scatter(Angle0,rel0,sz,'filled','MarkerFaceColor',c{6},'MarkerFaceAlpha',0.75,'DisplayName','Measured Relative Strain');
else
    for i = 1
        sc_rel{i} = scatter(Angle{i},rel{i},sz,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',c{7-2*i},'DisplayName',Test(i));
    end
end
hold off
title('Relative strain')
xlabel('Knee angle, \circ')
ylabel('strain/kmax')
set(ax2_1,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05],'YLim',[0 1.5],'XLim',[-125 35]);
ax2_1.FontName = 'Arial';
ax2_1.YAxis.LineWidth = 2; ax2_1.YAxis.FontSize = 10;
ax2_1.XAxis.LineWidth = 2; ax2_1.XAxis.FontSize = 10;
lgd2_1 = legend('Location','northwest');
lgd2_1.FontSize = 8;

ax2_2 = subplot(2,1,2);
hold on
plot(phiD,relstrain_ten,'Linewidth',2,'DisplayName','Expected Relative Strain')
if ~iscell(Angle)
    sc_rel = scatter(Angle1,rel1,sz,'filled','MarkerFaceColor',c{6},'MarkerFaceAlpha',0.75,'DisplayName','Measured Relative Strain');
else
    for i = 2:length(Angle)
        sc_rel{i} = scatter(Angle{i},rel{i},sz,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',c{7-2*i},'DisplayName',Test(i));
    end
end
hold off
title('Relative strain')
xlabel('Knee angle, \circ')
ylabel('strain/kmax')
set(ax2_2,'FontSize', 12, 'FontWeight', 'bold','XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05],'YLim',[0 1.5],'XLim',[-125 35]);
ax2_2.FontName = 'Arial';
ax2_2.YAxis.LineWidth = 2; ax2_2.YAxis.FontSize = 10;
ax2_2.XAxis.LineWidth = 2; ax2_2.XAxis.FontSize = 10;
lgd2_2 = legend('Location','northwest');
lgd2_2.FontSize = 8;

%% Plot measured versus expected BPA length
MuscleLength = Vas_Pam_42cm.MuscleLength-2*fitting-tendon{1};
MuscleLength_ten = Vas_Pam_42cm_tendon.MuscleLength-2*fitting-tendon{2};

figure
Lm1 = subplot(2,1,1);
hold on
pLm1_1 = plot(phiD,MuscleLength,'Color',c{5},'DisplayName','\bf Expected, no tendon');
if ~iscell(Angle)
    sLm1_1 = scatter(Angle0,InflatedLength0,sz,'MarkerFaceColor',c{7},'DisplayName','\bf Measured, no tendon');
else
    sLm1_1 = scatter(Angle{1},InflatedLength{1},'DisplayName','\bf Measured {l_{m}}, no tendon');
end
hold off
title('\bf Expected vs measured l_{m}, no tendon','Interpreter','tex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf l_{M}, m','Interpreter','tex')
set(Lm1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial','TickLength',[0.025, 0.05])
set(Lm1,'XLim',xLim,'XMinorTick','on','YMinorTick','on')
lgdLm1 = legend('Interpreter','tex','Location','southwest');
lgdLm1.FontSize = 8;

Lm2 = subplot(2,1,2);
hold on
pLm2_1 = plot(phiD,MuscleLength_ten,'Color',c{5},'DisplayName','\bf Expected');
if ~iscell(Angle)
sLm2_1 = scatter(Angle1,InflatedLength1,sz,'MarkerFaceColor',c{7},'DisplayName','\bf w/ tendon');
else
sLm2_1 = scatter(Angle{2},InflatedLength{2},'DisplayName','\bf w/ tendon (slip)');
sLm2_2 = scatter(Angle{3},InflatedLength{3},'DisplayName','\bf w/ tendon');
end
hold off
title('\bf Expected vs measured l_{m}, w/ tendon','Interpreter','tex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf l_{m}, m','Interpreter','tex')
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
    scH = scatter(Angle{1},TorqueHand{1},sz,'filled','MarkerFaceColor',c{1},'DisplayName','Back calculated');
end
hold off
title('l_{rest} = 41.5cm, no tendon','Interpreter','tex')
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
PL2 = plot(phiD, Theoretical_ten,'Color',c{5},'Linewidth',2,'DisplayName','\bf Theoretical, w/ tendon');
if ~iscell(Angle)
 scM1 = scatter(Angle1,Torque1,sz,'filled','MarkerFaceColor',c{7},'DisplayName','Measured');
 scH1 = scatter(Angle1,TorqueHand1,sz,'filled','MarkerFaceColor',c{2},'DisplayName','Hybrid calc'); 
else
 scM1 = scatter(Angle{2},Torque{2},sz,'d','filled','MarkerFaceColor',c{6},'DisplayName',Test(2));
 scH1 = scatter(Angle{2},TorqueHand{2},sz,'filled','MarkerFaceColor',c{4},'DisplayName',Test(2));
 scM2 = scatter(Angle{3},Torque{3},sz,'d','filled','MarkerFaceColor',c{7},'DisplayName',Test(3));
 scH2 = scatter(Angle{3},TorqueHand{3},sz,'filled','MarkerFaceColor',c{1},'DisplayName',Test(3));
end
hold off
title(sprintf('Torque, l_{rest} = 41.5cm, tendon = %g mm',tendon{2}*10^3),'Interpreter','tex')
xlabel('Knee angle, \circ','FontWeight','bold','Interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','Interpreter','tex')
set(gca2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',xLim,'TickLength',[0.025, 0.05])
set(gca2,'XMinorTick','on','YMinorTick','on');
lgd2 = legend('Interpreter','tex');
lgd2.FontSize = 12;
hold off

