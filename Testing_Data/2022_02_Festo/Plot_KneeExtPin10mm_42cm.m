%% Pinned knee, Extensor
%Run and save data from testing results
clear;
clc;
close all;

load 'KneeExtPin_10mm_all.mat'
Theoretical = Torque_42cm;

tendon0 = 0;            %no tendon condition
% tendon22 = 0.022;       %22 mm tendon

tendon22 = 0.022;       %Adjust tendon length
rest = 415/1000;        %resting length clamp to clamp, minus the barb
kmax = 0.349;           %length at maximum contraction
pres = 605.2351;        %Pressure, kPa
Vas_Pam_42cm_tendon = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon22, fitting, pres);
Torque_42cm_ten = Vas_Pam_42cm_tendon.Torque(:,3,4);

Theoretical_ten = Torque_42cm_ten;


rest = 0.415; %resting length, m
kmax = 0.349; %Length at maximum contraction, m
%% Test 1 done with CALT load cell. Tests 2 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet ExtTest10mm_4 from Results_table10mm_pinned_LoadCell   (no tendon)
%Test 2 == sheet ExtTest10mm_5 from Results_table10mm_pinned_LoadCell   (tendon, test w/ wrapping point #4 +20mm Z direction (i.e. BPA slipped from bolt) 
%Test 3 == sheet ExtTest10mm_6 from Results_table10mm_pinned_LoadCell   (tendon)

%% Torque calculated from measurements

Angle1 = [-76.5	-69	-51.5	-40.5	-39.5	-55.5	-35	-25.5	-20	-12.5	-5	0.5	4	12.5	20	28	35	23	20	11	1	-16	-24	-44	-53.4	-66];
Angle2 = [-120	-109	-101	-97	-89	-87	-81.5	-75	-74	-59.5	-48.5	-3	-30	-18.5];
Angle3 = [1.5	-5.5	-14.5	-20	-28.5	-35.5	-40	-50.5	-60	-66.5	-72.5	-78.5	-88	-94.5	-105	-112.5	-119.5];
Angle = [Angle1, Angle2, Angle3];

Torque1 = [8.50477362	7.539449285	6.619468483	5.835581619	6.120370289	7.005900186	6.178291184	6.391414397	5.96804559	6.099105123	5.93740855	4.880958611	4.375659895	3.388083398	2.795862193	1.758237576	0.306753172	1.725902869	3.204095044	5.145222736	6.304242352	6.817870142	7.511919159	7.911648561	8.110527045	8.862066797];
Torque2 = [4.032777251	3.725074705	3.501289111	3.342545878	3.070623054	2.844765784	2.844621771	2.773344025	2.674409973	2.382410213	1.680520864	0.394827933	1.758713268	1.550440294];
Torque3 = [0.55424608	0.941512328	1.76545047	1.79737823	2.650554833	2.7667878	2.211239585	3.275126705	3.423290801	3.611274245	3.687991538	3.833457994	3.587645288	3.554650905	3.949198619	4.568831269	4.408023267];
Torque = [Torque1, Torque2, Torque3];

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.

InflatedLength1 = [413.5	405.5	402	400	403	412	401	398	392.5	385	385.5	371.5	365	361	356	355	346	348	354	365	376	398	395	400	407	413]/1000;
InflatedLength2 = [408	402	380	394	390.5	387	383.5	383	381	371	363	340	355	353]/1000;
InflatedLength3 = [340	340	350	350.5	354	362	366	363	368	373	376	383	384	394	391	394.5	395]/1000;
InflatedLength = [InflatedLength1, InflatedLength2 InflatedLength3];

ICRtoMuscle1 = [28	30.5	31.5	31	32	31	35	32	41	44.5	45.5	55.5	66.5	62	65	68.5	71	67	60	55	55	42	38	31	29	28]/1000;
ICRtoMuscle2 = [42	42	42	42	41.5	36	37	36	36	30.5	22	42	28	24]/1000;
ICRtoMuscle3 = [55	43.5	42.5	37	35	33	35	34	35	34	36	36	36	36	36	36	36]/1000;
ICRtoMuscle = [ICRtoMuscle1, ICRtoMuscle2, ICRtoMuscle3];

F1 = zeros(1,size(InflatedLength1, 2));
F2 = zeros(1,size(InflatedLength2, 2));
F3 = zeros(1,size(InflatedLength3, 2));
F = [F1, F2 F3];

TorqueHand1 = zeros(1,size(InflatedLength1, 2));
TorqueHand2 = zeros(1,size(InflatedLength2, 2));
TorqueHand3 = zeros(1,size(InflatedLength3, 2));
TorqueHand = [TorqueHand1, TorqueHand2 TorqueHand3];

%load pressure where applicable
test = [4 5 6];
runsperseries = [26 14 17];

pres1 = zeros(1,runsperseries(1));
pres2 = zeros(1,runsperseries(2));
pres3 = zeros(1,runsperseries(3));
    for i = 1:3
        for j = 1:runsperseries(i)
            if i == 1
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres1(1,j) = Stats{'Mean',2};
            elseif i ==2
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres2(1,j) = Stats{'Mean',2};
            elseif i ==3
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres3(1,j) = Stats{'Mean',2};
            else
            end
        end       
    end

pres = [pres1 pres2 pres3];
% pressy = nonzeros(pres);
% pres = mean(pressy).*ones(1,size(pres,2));

for i = 1:size(InflatedLength, 2)  
    F(1,i,1) = festo3(InflatedLength(i), rest, 10, pres(i), kmax);    
    TorqueHand(i) = ICRtoMuscle(i)*F(i);  %Torque will be positive because it is causing extension
end

KMAX = (rest-kmax)/rest;
rel = ((rest-InflatedLength)/rest)/KMAX;
Fn = bpaForce10(rest,rel,pres);

for i = 2:(size(Fn,3)+1)
    F(:,:,i) = Fn(:,:,i-1);
end

for i = 1:size(F,3)
    TorqueHand(:,:,i) = ICRtoMuscle.*F(1,:,i);  %Torque will be positive because it is causing extension   
    TorqueHand1(:,:,i) = TorqueHand(1,1:size(TorqueHand1,2),i);
    TorqueHand2(:,:,i) = TorqueHand(1,((size(TorqueHand1,2)+1)):(size(TorqueHand1,2)+size(TorqueHand2,2)),i);
    TorqueHand3(:,:,i) = TorqueHand(1,(size(TorqueHand1,2)+size(TorqueHand2,2)+1):size(TorqueHand,2),i);
end

%% Prepare for plotting
% Create accessible color scheme
c1 = '#FFD700'; %gold
c2 = '#FFB14E'; %orange
c3 = '#FA8775'; %light orange
c4 = '#EA5F94'; %pink
c5 = '#CD34B5'; %magenta
c6 = '#9D02D7'; %magenta 2
c7 = '#0000FF'; %indigo
c8 = '#000000'; %black
sz = 60;        %size of data points
sz2 = sz*0.666; %size of second data points
c = {c1; c2; c3; c4; c5; c6; c7; c8};

%% Plot expected versus measured moment arm
Ma1 = Vas_Pam_42cm.MomentArm;                  %Calculated moment arm
G1 = (Ma1(:,1).^2+Ma1(:,2).^2).^(1/2);         %Moment arm for z-axis torque
Ma2 = Vas_Pam_42cm_tendon.MomentArm;                  %Calculated moment arm
G2 = (Ma2(:,1).^2+Ma2(:,2).^2).^(1/2);         %Moment arm for z-axis torque

figure
ax1_1 = subplot(2,1,1);
hold on
title('\bf Expected vs measured $r_{\hat{k}}$, no tendon', 'Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Z axis $r_{\hat{k}}$, m','Interpreter','latex')
pp1_1 = plot(phiD,G1,'Color',c7,'LineWidth',2,'DisplayName','\bf Expected $r_{\hat{k}}$ w/o tendon');
ss1_1 = scatter(Angle1, ICRtoMuscle1,sz,'filled','MarkerFaceColor',c4,'DisplayName','\bf Measured $r_{\hat{k}}$');
set(ax1_1,'FontSize', 12, 'FontWeight', 'bold','LineWidth', 2, 'FontName','Arial')
lgdMa1 = legend('Interpreter','latex','Location','Northwest');
lgdMa1.FontSize = 12;
hold off

ax1_2 = subplot(2,1,2);
hold on
title('\bf Expected vs measured $r_{\hat{k}}$, tendon', 'Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Z axis $r_{\hat{k}}$, m','Interpreter','latex')
pp2_1 = plot(phiD,G2,'Color',c7,'LineWidth',2,'DisplayName','\bf Expected $r_{\hat{k}}$ w/ tendon');
ss2_1 = scatter(Angle2, ICRtoMuscle2,sz,'filled','MarkerFaceColor',c2,'DisplayName','\bf Measured $r_{\hat{k}}$, tendon (slip)');
ss2_2 = scatter(Angle3, ICRtoMuscle3,sz,'filled','MarkerFaceColor',c4,'DisplayName','\bf Measured $r_{\hat{k}}$, tendon');
set(ax1_2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial')
lgdMa12 = legend('Interpreter','latex','Location','Northwest');
lgdMa12.FontSize = 12;
set(ax1_2,'XLim',[-120 40],'YLim',[0.02 0.08])
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values
%Calculated strain
strain = (rest-(Vas_Pam_42cm.MuscleLength-tendon0-2*fitting))/rest;          %strain w/o tendon
strain_ten = (rest-(Vas_Pam_42cm_tendon.MuscleLength-tendon22-2*fitting))/rest;  %strain w/ tendon
%Calculated relative strain
relstrain = (strain)./KMAX;                             %w/o tendon
relstrain_ten = (strain_ten)./KMAX;                     %w/ tendon
%Measured actual relative strain
realRel1 = (rest-InflatedLength1)/rest/KMAX;            %w/o tendon
realRel2 = (rest-InflatedLength2)/rest/KMAX;            %w/ tendon (BPA slipped)
realRel3 = (rest-InflatedLength3)/rest/KMAX;            %w/ tendon

figure %Show same length BPA w/ and w/o tendon               
ax2_1 = subplot(2,1,1);
hold on
title('Expected vs measured \epsilon^*,no tendon','Interpreter','tex')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('\epsilon^*','Interpreter','tex')
pE1_1 = plot(phiD,relstrain,'Color',c7,'LineWidth',2,'DisplayName','\bf Expected \epsilon^*, no tendon');
sE1_1 = scatter(Angle1,realRel1,sz,'filled','MarkerFaceColor',c4,'DisplayName','\bf Measured \epsilon^*, no tendon');
set(ax2_1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdEp1 = legend('Interpreter','tex','Location','Northwest');
lgdEp1.FontSize = 12;
set(ax2_1,'XLim',[-120 40],'YLim',[0 1.2])
hold off

ax2_2 = subplot(2,1,2);
hold on
title('Expected vs measured \epsilon^*, w/tendon','Interpreter','tex')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('\epsilon^*','Interpreter','tex')
pE2_1 = plot(phiD,relstrain_ten,'Color',c7,'LineWidth',2,'DisplayName','\bf Expected \epsilon^* w/ tendon');
sE2_1 = scatter(Angle2,realRel2,sz,'filled','MarkerFaceColor',c2,'DisplayName','\bf Measured \epsilon^* w/ tendon (slip)');
sE2_2 = scatter(Angle3,realRel3,sz,'filled','MarkerFaceColor',c4,'DisplayName','\bf Measured \epsilon^* w/ tendon');
set(ax2_2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdEp2 = legend('Interpreter','tex','Location','Southeast');
lgdEp2.FontSize = 12;
set(ax2_2,'XLim',[-120 40],'YLim',[0 1.2])
hold off

%% Plot measured versus expected BPA length
MuscleLength = Vas_Pam_42cm.MuscleLength-2*fitting-tendon0;
MuscleLength_ten = Vas_Pam_42cm_tendon.MuscleLength-2*fitting-tendon22;

figure
Lm1 = subplot(2,1,1);
hold on
title('\bf Expected vs measured {l_{m}}, no tendon','Interpreter','tex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf {l_{m}}, m','Interpreter','tex')
pLm1_1 = plot(phiD,MuscleLength,'DisplayName','\bf Expected {l_{m}}, no tendon');
sLm1_1 = scatter(Angle1,InflatedLength1,'DisplayName','\bf Measured {l_{m}}, no tendon');
set(Lm1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdLm1 = legend('Interpreter','tex');
lgdLm1.FontSize = 12;
hold off

Lm2 = subplot(2,1,2);
hold on
title('\bf Expected vs measured {l_{m}}, w/ tendon','Interpreter','tex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf {l_{m}}, m','Interpreter','tex')
pLm2_1 = plot(phiD,MuscleLength_ten,'DisplayName','\bf Expected {l_{m}}, w/ tendon');
sLm2_1 = scatter(Angle2,InflatedLength2,'DisplayName','\bf Measured {l_{m}}, w/ tendon (slip)');
sLm2_2 = scatter(Angle3,InflatedLength3,'DisplayName','\bf Measured {l_{m}}, w/ tendon');
set(Lm2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdLm2 = legend('Interpreter','tex');
lgdLm2.FontSize = 12;
set(Lm2,'XLim',[-120 40],'YLim',[0.34 0.46])
hold off

%% Plotting Z axis Torque values
figure
gca1 = subplot(2,1,1);
hold on
title('Iso. Torque vs {\theta_{k}}, Pinned, Extensor, l_{rest} = 45.7cm, no tendon','Interpreter','tex')
xlabel('Knee angle, \circ','FontWeight','bold','Interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','Interpreter','tex')
PL1 = plot(phiD, Theoretical,'Color',c4,'Linewidth',2,'DisplayName','Theoretical Torque');
scM = scatter(Angle1,Torque1,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Measured Torque');
scH = scatter(Angle1,TorqueHand1(:,:,4),sz2,'filled','MarkerFaceColor',c1,'DisplayName','Back calculated Torque');
set(gca1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial')
lgd1 = legend;
lgd1.FontSize = 12;
hold off

gca2 = subplot(2,1,2);
hold on
title('Iso. Torque vs {\theta_{k}}, Pinned, Extensor, l_{rest} = 45.7cm, 22mm tendon','Interpreter','tex')
xlabel('Knee angle, \circ','FontWeight','bold','Interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','Interpreter','tex')
PL2 = plot(phiD, Theoretical_ten,'Color',c7,'Linewidth',2,'DisplayName','\bf Theoretical, 22mm tendon');
scM1 = scatter(Angle2,Torque2,sz,'d','filled','MarkerFaceColor',c6,'DisplayName','\bf Measured (slip)');
scH1 = scatter(Angle2,TorqueHand2(:,:,4),sz2,'filled','MarkerFaceColor',c4,'DisplayName','\bf Back calculated (slip)');
scM2 = scatter(Angle3,Torque3,sz,'d','filled','MarkerFaceColor',c2,'DisplayName','\bf Measured');
scH2 = scatter(Angle3,TorqueHand3(:,:,4),sz2,'filled','MarkerFaceColor',c1,'DisplayName','\bf Back calculated');
set(gca2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial')
lgd2 = legend('Interpreter','tex');
lgd2.FontSize = 12;
set(gca2,'XLim',[-120 40])
hold off

%% Mean and RMSE
Tqz = cell(2,1);
Tqz{1} = Torque_42cm;        %Calculated Torque, no tendon, new simplified exponential equation
Tqz{2} = Torque_42cm_ten;    %Calculated Torque, 22mm tendon, new simplified exponential equation, no model modification for BPA slipping off of screw
Tqz{3} = Torque_42cm_ten;    %Calculated Torque, 22mm tendon, new simplified exponential equation

%fit options
mod_Pam = fittype('cubicinterp');
Options = fitoptions(mod_Pam);
Options.Normal = 'on';

%prepare cells
mdl_Pam = cell(size(Tqz,1));
val = cell(length(mdl_Pam),1);
          
%Get values at each angle there is measurement data for
for j = 1:length(Tqz)
     Options.Exclude = isnan(Tqz{j});
     mdl_Pam{j} = fit(phiD',Tqz{j},mod_Pam,Options);
end

y = cell(3,1);
y{1} = Torque1';        %no tendon
y{2} = Torque2';        %w/ tendon and BPA slip
y{3} = Torque3';        %w/ tendon
ynew{1} = feval(mdl_Pam{1},Angle1');
ynew{2} = feval(mdl_Pam{2},Angle2');
ynew{3} = feval(mdl_Pam{3},Angle3');
        
yresid = cell(length(ynew),1);
SSresid = cell(length(ynew),1);
fu = cell(length(ynew),1);
        
for i = 1:length(ynew)
    yresid{i} = y{i}-ynew{i};              %residual error
    SSresid{i} = sum(yresid{i}.^2,'omitnan'); %Sum of squares of the residual
    fu{i} = sqrt(SSresid{i}/length(yresid{i}));        % RMSE for function 1
end

fprintf('Original torque calculation, no tendon, returns SSE of %5d with an RMSE of %5d\n',SSresid{1},fu{1})
fprintf('Original torque calculation, tendon w/ BPA slip, returns SSE of %5d with an RMSE of %5d\n',SSresid{2},fu{2})
fprintf('Original torque calculation, tendon, returns SSE of %5d with an RMSE of %5d\n',SSresid{3},fu{3})