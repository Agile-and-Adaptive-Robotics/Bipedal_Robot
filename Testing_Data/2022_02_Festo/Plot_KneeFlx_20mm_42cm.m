%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

f = fullfile('github/bipedal_robot/code/matlab');
qt = addpath(genpath(f));

TabMA = readmatrix('OpenSim_Bifem_MomentArm.txt');
knee_angle_rMA = TabMA(:,2)';           %Angle values directly from O
Bifemsh_MA = TabMA(:,3)';              %Torque values directly from OpenSim

Tab = readmatrix('OpenSim_Bifem_Results.txt');
knee_angle_rT = Tab(:,2)';           %Angle values directly from O
Bifemsh_T = Tab(:,4)';              %Torque values directly from OpenSim

load KneeFlx_20mm_42cm.mat
% rest = 0.423; %resting length, m
% kmax = 0.322; %Length at maximum contraction, m
rest = 0.420; %resting length, m
kmax = 0.315; %Length at maximum contraction, m
dia = 20;
fitting = 0.0254;
KMAX = (rest-kmax)/rest;

Theo1 = Bifemsh_Pam1.Torque(:,3)';          %Original Torque calculations
Theo2 = Bifemsh_Pam2.Torque(:,3)';
Theo3 = Bifemsh_Pam3.Torque(:,3)';
Theo1opt = Bifemsh_Pam_adj1.Torque(:,3)';       %Optimized fitting length from previous study
Theo2opt = Bifemsh_Pam_adj2.Torque(:,3)';
Theo3opt = Bifemsh_Pam_adj3.Torque(:,3)';


%% Test 1 done with CALT load cell
%Test 1 == sheet FlxTest20mm from Results_table10mm_FullSize, data 1:10, 274 kPa
%Test 2 == sheet FlxTest20mm from Results_table10mm_FullSize, data 11:15, 484.8 kPa
%Test 3 == sheet FlxTest20mm from Results_table10mm_FullSize, data 16:17, 606.5 kPa
%% Torque calculated from measurements
Angle = [-2	-11.5	-20	-27	-35	-44	-54.5	-50	-63	-64	-73.5	-82	-88	-103.5	-114	-115	-120];
Angle1 = Angle(1:10);
Angle2 = Angle(11:15);
Angle3 = Angle(16:17);

Torque = [-8.84322291179594,-7.38771014573315,-6.12295907061857,-5.11416855627929,-4.00577506118469,-3.08271954812317,-3.04251451225018,-2.59384366567520,-1.86178531541841,-0.774849368079027,-4.65440909149842,-3.19703540259264,-2.57966382173969,-0.959098909344151,-0.499572417339300,-1.33299525437727,-0.772724305518311];
Torque1 = Torque(1:10);
Torque2 = Torque(11:15);
Torque3 = Torque(16:17);
%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength = [401	395	389	387	382	380	378	374.5	371.5	369.5	354	352.5	350.5	345	343.5	338.5	339]/1000;
InflatedLength1 = InflatedLength(1:10);
InflatedLength2 = InflatedLength(11:15);
InflatedLength3 = InflatedLength(16:17);

ICRtoMuscle = [45	47	51	45	40	37	38	39	36.5	35.5	30.5	20	22	16	4.5	 3	10]/1000;
ICRtoMuscle1 = ICRtoMuscle(1:10);
ICRtoMuscle2 = ICRtoMuscle(11:15);
ICRtoMuscle3 = ICRtoMuscle(16:17);

F = zeros(1,size(InflatedLength, 2));
F1 = F(1:10);
F2 = F(11:15);
F3 = F(16:17);

TorqueHand = zeros(1,size(InflatedLength, 2));

%load pressure where applicable
test = 1;
runsperseries = 17;

pres = zeros(1,runsperseries);
 

     for j = 1:runsperseries
                file_name = sprintf('FlxTest%0.0f_%0.0f.mat', test,j);
                load(file_name,'Stats')
                pres(1,j) = Stats{'Mean',2};
     end
pres1 = pres(1:10);
pres2 = pres(11:15);
pres3 = pres(16:17);

contract = (rest-InflatedLength)/rest;
realRel = contract./KMAX;

F = festo4(20, realRel, pres);
F_alt = festo4(20, realRel, 600);
TorqueHand = -ICRtoMuscle.*F;  %Moment will be negative because it causes flexion
TorqueHand_alt = -ICRtoMuscle.*F_alt;  %Potentially acheivable force

TorqueHand1 = TorqueHand(1:10);
TorqueHand2 = TorqueHand(11:15);
TorqueHand3 = TorqueHand(16:17);

%% Mean and RMSE
Tqz = cell(3,1);
Tqz{1} = Theo1';         %Calculated Torque, new simplified exponential equation w/o optimized fitting length
Tqz{2} = Theo2'; 
Tqz{3} = Theo3'; 
Topt = cell(3,1);
Topt{1} = Theo1opt';      %Calculated Torque, adjusted with optimized fitting length
Topt{2} = Theo2opt';
Topt{3} = Theo3opt';
TH = cell(3,1);
TH{1} = TorqueHand1';   %Placeholder in case we want to compare SSE/RMSE of back calculated torque to measured torque
TH{2} = TorqueHand2';
TH{3} = TorqueHand3';

PRZ = cell(3,1);
PRZ{1} = mean(pres1);
PRZ{2} = mean(pres2);
PRZ{3} = mean(pres3);

ANG = cell(3,1);
ANG{1} = Angle1';
ANG{2} = Angle2';
ANG{3} = Angle3';

%fit options
mod_Pam = fittype('cubicinterp');
Options = fitoptions(mod_Pam);
Options.Normal = 'on';

%prepare cells
mdl_Pam = cell(size(Tqz,1));
mdl_Pam_opt = cell(size(Topt,1));
mdl_PamH = cell(size(TH,1));
val = cell(length(mdl_Pam),1);
val_opt = cell(length(mdl_Pam_opt),1);
valH = cell(length(mdl_PamH),1);

%Get values at each angle there is measurement data for
for j = 1:length(Tqz)
     Options.Exclude = isnan(Tqz{j});
     mdl_Pam{j} = fit(phiD',Tqz{j},mod_Pam,Options);
     Options.Exclude = isnan(Topt{j});
     mdl_Pam_opt{j} = fit(phiD',Topt{j},mod_Pam,Options);
     Options.Exclude = isnan(TH{j});
     mdl_PamH{j} = fit(ANG{j},TH{j},mod_Pam,Options);
     val{j} = feval(mdl_Pam{j},ANG{j});
     val_opt{j} = feval(mdl_Pam_opt{j},ANG{j});
     valH{j} = feval(mdl_PamH{j},ANG{j});
end

y = cell(3,1);
y{1} = Torque1';
y{2} = Torque2';
y{3} = Torque3'; 
ynew = val;
ynew_opt = val_opt;
ynewH = valH;
        
yresid = cell(length(ynew),1);
yresid_opt = cell(length(ynew_opt),1);
yresidH = cell(length(ynewH),1);

SSresid = cell(length(ynew),1);
SSresid_opt = cell(length(ynew_opt),1);
SSresidH = cell(length(ynewH),1);

fu = cell(length(ynew),1);
fu_opt = cell(length(ynew_opt),1);
fuH = cell(length(ynewH),1);
        
for i = 1:length(ynew)
    yresid{i} = y{i}-ynew{i};              %residual error
    yresid_opt{i} = y{i}-ynew_opt{i};              %residual error
    yresidH{i} = y{i}-ynewH{i};              %residual error
    SSresid{i} = sum(yresid{i}.^2,'omitnan'); %Sum of squares of the residual
    SSresid_opt{i} = sum(yresid_opt{i}.^2,'omitnan'); %SSE
    SSresidH{i} = sum(yresidH{i}.^2,'omitnan'); %SSE
    fu{i} = sqrt(SSresid{i}/length(yresid{i}));        % RMSE for function 1
    fu_opt{i} = sqrt(SSresid_opt{i}/length(yresid_opt{i}));        % RMSE
    fuH{i} = sqrt(SSresidH{i}/length(yresidH{i}));        % RMSE
end

for i = 1:length(ANG)
    fprintf('Original torque calculation at %.0f kPa returns SSE of %5d with an RMSE of %5d\n',PRZ{i},SSresid{i},fu{i})
    fprintf('Optimized torque calculation at %.0f kPa returns SSE of %5d with an RMSE of %5d\n',PRZ{i},SSresid_opt{i},fu_opt{i})
    fprintf('Back calculated torque at %.0f kPa returns SSE of %5d with an RMSE of %5d\n',PRZ{i},SSresidH{i},fuH{i})
end
%% Prepare for plots by creating colors and sizes
%Matlab hex color values:
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
Ma = Bifemsh_Pam_adj3.MomentArm;                 %Calculated moment arm
G = -(Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque

figure
hold on
pp2 = plot(phiD,G,'DisplayName','\bf Expected $r_{\hat{k}}$');
ss2 = scatter(Angle, -ICRtoMuscle,'DisplayName','\bf Measured $r_{\hat{k}}$');
ppOS = plot(knee_angle_rMA,Bifemsh_MA,'--','DisplayName','\bf OpenSim $r_{\hat{k}}$');
title('\bf Expected vs measured $r_{\hat{k}}$', 'Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Z axis $r_{\hat{k}}$, m','Interpreter','latex')
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 12;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 12;
lgdMa = legend('Interpreter','latex');
lgdMa.FontSize = 12;
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values
strain = Bifemsh_Pam_adj3.Contraction;
relstrain = (strain)./KMAX;

figure
hold on
plot(phiD,relstrain,'DisplayName','\bf Expected \epsilon^*')
scatter(Angle,realRel,'DisplayName','\bf Measured \epsilon^*')
title('Expected vs measured \epsilon^*','Interpreter','tex')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('\epsilon^*','Interpreter','tex')
ax2 = gca;
ax2.FontSize = 12;
ax2.FontWeight = 'bold';
ax2.FontName = 'Arial';
ax2.YAxis.LineWidth = 2; ax2.YAxis.FontSize = 12;
ax2.XAxis.LineWidth = 2; ax2.XAxis.FontSize = 12;
lgdMa = legend('Interpreter','tex');
lgdMa.FontSize = 12;
hold off

%% Plot measured versus expected BPA length
MuscleLength = Bifemsh_Pam_adj3.MuscleLength-2*fitting-tendon;

figure
hold on
plot(phiD,MuscleLength,'DisplayName','\bf Expected $l_{m}$')
scatter(Angle,InflatedLength,'DisplayName','\bf Measured $l_{m}$')
title('\bf Expected vs measured $l_{m}$','Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf $l_{m}$, m','Interpreter','latex')
ax3 = gca;
ax3.FontSize = 12;
ax3.FontWeight = 'bold';
ax3.FontName = 'Arial';
ax3.YAxis.LineWidth = 2; ax3.YAxis.FontSize = 12;
ax3.XAxis.LineWidth = 2; ax3.XAxis.FontSize = 12;
lgdMa = legend('Interpreter','latex');
lgdMa.FontSize = 12;
hold off

%% Plotting Isometric Torque Values
% figure
% hold on
% gca1 = gca;
% gcf1 = gcf;
% 
% pl3 = plot(phiD, Theo3,'Color',c2,'Linewidth',2,'DisplayName','Theoretical, 606 kPa');
% pl2 = plot(phiD, Theo2,'--','Color',c2,'Linewidth',2,'DisplayName','Theoretical, 485 kPa');
% pl1 = plot(phiD, Theo1,'-.','Color',c2,'Linewidth',2,'DisplayName','Theoretical, 274 kPa');
% 
% pl6 = plot(phiD, Theo3opt,'Color',c4,'Linewidth',2,'DisplayName','Theoretical, 606 kPa');
% pl5 = plot(phiD, Theo2opt,'--','Color',c4,'Linewidth',2,'DisplayName','Theoretical, 485 kPa');
% pl4 = plot(phiD, Theo1opt,'-.','Color',c4,'Linewidth',2,'DisplayName','Theoretical, 274 kPa');
% 
% sc1 = scatter(Angle1,Torque1,sz,'d','MarkerFaceColor',c7,'DisplayName','Measured, 274 kPa');
% sc2 = scatter(Angle2,Torque2,sz,'d','MarkerFaceColor',c6,'DisplayName','Measured, 485 kPa');
% sc3 = scatter(Angle3,Torque3,sz,'d','MarkerFaceColor',c5,'DisplayName','Measured, 606 kPa');
% 
% sc4 = scatter(Angle,TorqueHand,sz,'filled','MarkerFaceColor',c1,'DisplayName','Back calculated');
% % sc7 = scatter(Angle,TorqueHand_alt,sz,'filled','MarkerFaceColor',c3,'DisplayName','Back calculated @ 600 kPa');
% 
% title('Iso. Torque vs {\theta_{k}}, {\phi}20mm Extensor, l_{rest} = 42.3cm','Interpreter','tex')
% xlabel('Knee angle, \circ','FontWeight','bold','Interpreter','tex')
% ylabel('Torque, N{\cdot}m','FontWeight','bold','Interpreter','tex')
% set(gca1,'FontSize', 12, 'LineWidth',2,'FontWeight', 'bold','FontName','Arial','XMinorGrid','off','XMinorTick','off','YMinorGrid','off','YMinorTick','off');
% lgd = legend;
% hold off

%% Plot Torque with interpolated values from MEAN and RMSE for theoretical values
%This way it is just over RoM measured at each pressure

Disp1 = cell(3,1);
Disp2 = cell(3,1);
Disp3 = cell(3,1);
PL = cell(3,1);
PL_opt = cell(3,1);
sc = cell(3,1);

figure
hold on
gca1 = gca;
gcf1 = gcf;


    for i = 1:size(ANG,1)
        T1 = 2*i;
        H1 = 2*i+1;
        Disp1{i} = sprintf('Theoretical, %.0f kPa',PRZ{i});
        Disp2{i} = sprintf('Theoretical optimized, %.0f kPa',PRZ{i});
        Disp3{i} = sprintf('Measured, %.0f kPa',PRZ{i});
        PL{i} = plot(ANG{i}, val{i},'--','Color',c{H1},'Linewidth',2,'DisplayName',Disp1{i});
        PL_opt{i} = plot(ANG{i}, val_opt{i},'Color',c{H1},'Linewidth',2,'DisplayName',Disp2{i});
        sc{i} = scatter(ANG{i},y{i},sz,'d','filled','MarkerFaceColor',c{H1},'DisplayName',Disp3{i});
    end


SC4 = scatter(Angle,TorqueHand,sz,'filled','MarkerFaceColor',c1,'DisplayName','Back calculated');
sc7 = scatter(Angle,TorqueHand_alt,sz,'v','MarkerEdgeColor',c8,'DisplayName','Back calculated @ 600 kPa');
OST = plot(knee_angle_rT,Bifemsh_T,'-.','Color',c8,'LineWidth',1,'DisplayName','OpenSim Torque');

title('Iso. Torque vs {\theta_{k}}, {\phi}20mm Flexor, l_{rest} = 42.3cm','Interpreter','tex')
xlabel('Knee angle, \circ','FontWeight','bold','Interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','Interpreter','tex')
set(gca1,'FontSize', 12, 'LineWidth',2,'FontWeight', 'bold','FontName','Arial','XMinorGrid','off','XMinorTick','off','YMinorGrid','off','YMinorTick','off');
lgd = legend;
hold off

